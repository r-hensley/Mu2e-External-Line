"""Check or refresh the current archive inventory using only Python's stdlib.

Paths are relative to the archive root. Inventory control files under manifests/
are excluded. Existing folder SHA256SUMS files are generated subsets of the same
inventory. Refresh requires explicit acceptance of changed listed files.
"""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re
import stat
import sys


FIELDS = ("path", "size_bytes", "modified_time", "mode_octal", "provenance", "sha256")


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def relative_name(name):
    path = Path(name)
    if path.is_absolute() or ".." in path.parts or not path.parts:
        raise ValueError(f"Invalid inventory path: {name!r}")
    return path.as_posix()


def archive_files(root):
    files = {}
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if relative.parts[0] == "manifests":
            continue
        if path.is_symlink():
            raise ValueError(f"Symlink requires an explicit inventory policy: {relative}")
        if path.is_file():
            files[relative.as_posix()] = path
    return files


def load_checksums(path):
    entries = {}
    for line in path.read_text().splitlines():
        escaped = line.startswith("\\")
        if escaped:
            line = line[1:]
        match = re.fullmatch(r"([0-9a-f]{64}) [ *](.+)", line)
        if not match:
            raise ValueError(f"Invalid checksum line in {path}: {line!r}")
        value, name = match.groups()
        if escaped:
            name = re.sub(r"\\([\\n])", lambda m: "\n" if m[1] == "n" else "\\", name)
        name = relative_name(name)
        if name in entries:
            raise ValueError(f"Duplicate checksum path: {name}")
        entries[name] = value
    return entries


def checksum_text(entries):
    lines = []
    for name, value in sorted(entries.items()):
        escaped = "\\" in name or "\n" in name
        name = name.replace("\\", "\\\\").replace("\n", "\\n")
        lines.append(("\\" if escaped else "") + f"{value}  ./{name}\n")
    return "".join(lines)


def load_table(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise ValueError(f"Unexpected TSV columns in {path}")
        rows = {}
        for row in reader:
            name = relative_name(row["path"])
            if name in rows:
                raise ValueError(f"Duplicate TSV path: {name}")
            rows[name] = row
    return rows


def known_sources(root, files):
    origins = {}
    for name, path in files.items():
        if not name.endswith("/provenance/file_manifest.json"):
            continue
        manifest = json.loads(path.read_text())
        bundle = path.parent.parent
        for field, origin in (("supplied_files", "received_supplied_file"),
                              ("reference_tables", "local_reference_copy")):
            for item in manifest[field]:
                source_name = relative_name(item["path"])
                source = bundle / source_name
                full_name = source.relative_to(root).as_posix()
                if source.stat().st_size != item["bytes"] or digest(source) != item["sha256"]:
                    raise ValueError(f"Supplied/reference file differs from its provenance record: {full_name}")
                origins[full_name] = origin
    return origins


def metadata(path, name, value, provenance):
    info = path.stat()
    seconds, nanos = divmod(info.st_mtime_ns, 1_000_000_000)
    date = datetime.fromtimestamp(seconds, timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")
    return {
        "path": name, "size_bytes": str(info.st_size),
        "modified_time": f"{date}.{nanos:09d}+00:00",
        "mode_octal": f"{stat.S_IMODE(info.st_mode):o}",
        "provenance": provenance, "sha256": value,
    }


def new_provenance(name, origins):
    if name in origins:
        return origins[name]
    if name.endswith("/SHA256SUMS"):
        return "generated_checksum_subset"
    if name.endswith((".md", ".txt")):
        return "local_documentation"
    if "/audit/run_" in name:
        return "local_validation_evidence"
    if "/audit/" in name:
        return "local_audit"
    if "/provenance/" in name:
        return "local_provenance"
    return "local_addition_unclassified"


def write_atomic(path, text):
    if path.is_file() and path.read_bytes() == text.encode("utf-8"):
        return
    temporary = path.with_name(path.name + ".tmp")
    try:
        temporary.write_text(text)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def refresh(root, accepted):
    sha_path = root / "manifests/m4_lattice_archive.sha256"
    table_path = root / "manifests/m4_lattice_archive.tsv"
    expected, previous = load_checksums(sha_path), load_table(table_path)
    if set(expected) != set(previous) or any(expected[n] != previous[n]["sha256"] for n in expected):
        raise ValueError("Existing SHA-256 and TSV inventories disagree")
    files = archive_files(root)
    missing = sorted(set(expected) - set(files))
    if missing:
        raise ValueError(f"Listed files are missing; review removals before refreshing: {missing}")
    values = {name: digest(path) for name, path in files.items()}
    changed = sorted(name for name in expected if values[name] != expected[name])
    unaccepted = sorted(set(changed) - accepted)
    if unaccepted:
        raise ValueError(f"Review changed files, then pass --accept-change for each intended change: {unaccepted}")
    origins = known_sources(root, files)
    # Generate portable lists first, deepest directory first, so parent lists
    # and the archive inventory include the final bytes of child lists.
    portable = sorted((p for p in files.values() if p.name == "SHA256SUMS"),
                      key=lambda p: len(p.parts), reverse=True)
    for path in portable:
        entries = {p.relative_to(path.parent).as_posix(): values[name]
                   for name, p in files.items()
                   if p.is_relative_to(path.parent) and p != path}
        write_atomic(path, checksum_text(entries))
        values[path.relative_to(root).as_posix()] = digest(path)
    rows = [metadata(path, name, values[name], previous[name]["provenance"]
                     if name in previous else new_provenance(name, origins))
            for name, path in files.items()]
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    write_atomic(table_path, output.getvalue())
    write_atomic(sha_path, checksum_text(values))
    print(f"Refreshed {len(files)} archive entries and {len(portable)} generated folder subset(s).")


def check(root):
    expected = load_checksums(root / "manifests/m4_lattice_archive.sha256")
    rows = load_table(root / "manifests/m4_lattice_archive.tsv")
    files = archive_files(root)
    missing = sorted(set(expected) - set(files))
    unlisted = sorted(set(files) - set(expected))
    changed = sorted(name for name in set(expected) & set(files)
                     if digest(files[name]) != expected[name])
    inconsistent = sorted(set(expected) ^ set(rows))
    inconsistent += sorted(name for name in set(expected) & set(rows)
                           if expected[name] != rows[name]["sha256"]
                           or name in files and files[name].stat().st_size != int(rows[name]["size_bytes"]))
    problems = {key: value for key, value in (
        ("missing", missing), ("unlisted", unlisted), ("changed", changed),
        ("TSV_inconsistent", inconsistent)) if value}
    if problems:
        print(json.dumps(problems, indent=2))
        return 1
    print(f"Verified {len(expected)} files: contents and complete path coverage match; SHA-256 and TSV inventories agree.")
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("check", "refresh"))
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--accept-change", action="append", default=[], metavar="RELATIVE_PATH",
                        help="Refresh only: explicitly accept an intended change to a listed file")
    args = parser.parse_args()
    if args.action == "check" and args.accept_change:
        parser.error("--accept-change applies only to refresh")
    try:
        root = args.root.resolve()
        if args.action == "refresh":
            refresh(root, {relative_name(name) for name in args.accept_change})
        return check(root)
    except (OSError, ValueError, KeyError) as error:
        print(f"Inventory error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
