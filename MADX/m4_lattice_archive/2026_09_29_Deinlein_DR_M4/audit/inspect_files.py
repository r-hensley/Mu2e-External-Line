"""Recompute the read-only source/table comparisons from the extracted supplied files.

Run from any directory with python3. This does not execute MAD-X. The adjacent
run directories preserve the previously generated validation inputs and outputs
with LF line endings. Original identities are retained in the provenance record.
Only this audit directory's receipt.json and comparison CSV are written.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import re
import shlex

AUDIT = Path(__file__).resolve().parent
BUNDLE = AUDIT.parent
WORKSPACE = BUNDLE.parents[3]
MANIFEST = BUNDLE / "provenance/file_manifest.json"
REFERENCES = BUNDLE / "audit/reference_tables"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def tfs(text):
    header, columns, rows = {}, [], []
    for line in text.splitlines():
        if line.startswith("@"):
            tokens = shlex.split(line)
            header[tokens[1]] = tokens[3]
        elif line.startswith("*"):
            columns = shlex.split(line)[1:]
        elif line.strip() and not line.startswith(("@", "*", "$")):
            tokens = shlex.split(line)
            assert len(tokens) == len(columns)
            rows.append({key: value if key in ("NAME", "KEYWORD") else float(value)
                         for key, value in zip(columns, tokens)})
    return header, columns, rows


def read_tfs(path):
    return tfs(path.read_text())


def active_source(data):
    # All source files here use // and ! line comments; there are no strings
    # containing these comment delimiters in the statements being inspected.
    return "\n".join(line.split("//", 1)[0].split("!", 1)[0]
                     for line in data.decode().splitlines())


manifest = json.loads(MANIFEST.read_text())
expected_files = {item["path"] for item in manifest["supplied_files"]}
assert len(manifest["supplied_files"]) == len(expected_files) == 6
assert expected_files == {f"mu2e-{kind}-{suffix}"
                          for kind in ("dr", "dabs")
                          for suffix in ("lattice-v2026.09.29.madx",
                                         "twiss-v2026.09.29-madx.tfs",
                                         "plots-v2026.09.29-madx.pdf")}
for item in manifest["supplied_files"] + manifest["reference_tables"]:
    data = (BUNDLE / item["path"]).read_bytes()
    assert len(data) == item["bytes"], item["path"]
    assert sha(data) == item["sha256"], item["path"]
receipt = {
    "audit_date": "2026-10-05",
    "recheck_date": "2026-10-06",
    "line_ending_normalization": manifest.get("line_ending_normalization"),
    "run_evidence_note": "MAD-X ran before LF normalization; current file sizes and hashes describe the LF copies.",
    "source_directory": str(BUNDLE),
    "file_manifest": str(MANIFEST), "manifest_sha256": sha(MANIFEST.read_bytes()),
    "provenance_basis": "User-supplied author correspondence summarized in the docs; extracted supplied files and saved validation runs",
    "email_accessed_independently": False,
    "runtime": {
        "executable": str(WORKSPACE / "madx.exe"),
        "sha256": "539815c372bdf3495bb97f4c7e8c75b4e30391649572d75d990da4106d743b6a",
        "version": "MAD-X 5.08.00 (64 bit, Windows)", "release_date": "2022.01.13",
        "temporary_directory": "/tmp/mu2e-lattice-20261005-r8eszn1s",
        "invocation": "Executable with source bytes on standard input; cwd is isolated run directory",
        "initial_sandbox_failure": "WSL ERROR: UtilBindVsockAnyPort:309: socket failed 1",
        "successful_execution": "WSL interoperability allowed outside sandbox",
    },
    "supplied_files": [{"name": item["path"], "bytes": item["bytes"], "sha256": item["sha256"]}
                       for item in manifest["supplied_files"]],
    "all_supplied_file_hashes_match": True, "models": {},
}
supplied_by_model = {}
for kind in ("dabs", "dr"):
    entry = f"mu2e-{kind}-lattice-v2026.09.29.madx"
    raw = (BUNDLE / entry).read_bytes()
    source = active_source(raw)
    head, columns, supplied = read_tfs(BUNDLE / f"mu2e-{kind}-twiss-v2026.09.29-madx.tfs")
    supplied_by_model[kind] = supplied
    generated_path = AUDIT / f"run_{kind}_export/mu2e-{kind}-twiss-madx.tfs"
    generated_head, generated_columns, generated = read_tfs(generated_path)
    assert columns == generated_columns
    assert len(supplied) == len(generated)
    assert [(r["NAME"], r["KEYWORD"]) for r in supplied] == [
        (r["NAME"], r["KEYWORD"]) for r in generated]
    maximum_differences = {column: max(abs(a[column] - b[column])
        for a, b in zip(supplied, generated))
        for column in columns if column not in ("NAME", "KEYWORD")}
    assert max(maximum_differences.values()) == 0
    runs = {}
    for mode in ("as_received", "export"):
        directory = AUDIT / f"run_{kind}_{mode}"
        log = (directory / "run.log").read_text(errors="replace").replace("X:> ", "")
        diagnostics = [line.strip() for line in log.splitlines()
                       if re.search(r"warning|fatal|finished normally", line, re.I)]
        assert "Number of warnings: 0" in diagnostics
        assert any("MAD-X finished normally" in line for line in diagnostics)
        if mode == "export":
            expected = raw.replace(b"!WRITE, TABLE = TWISS", b"WRITE, TABLE = TWISS").replace(
                b"!WRITE, TABLE = SURVEY", b"WRITE, TABLE = SURVEY")
            assert (directory / "validation_export.madx").read_bytes() == expected
        runs[mode] = {"exit_status": 0, "warnings": 0, "fatal": None,
                      "diagnostics": diagnostics,
                      "files": [{"name": p.name, "bytes": p.stat().st_size,
                                 "sha256": sha(p.read_bytes())}
                                for p in sorted(directory.iterdir()) if p.is_file()]}
    survey = read_tfs(AUDIT / f"run_{kind}_export/mu2e-{kind}-survey-madx.tfs")
    receipt["models"][kind] = {
        "entry_point": entry, "external_call_or_read_dependencies": [],
        "active_call_or_read_commands": re.findall(r"\b(?:CALL|READ)\s*,[^;]*;", source, re.I),
        "active_output_commands": re.findall(r"\b(?:WRITE|PLOT)\s*,[^;]*;", source, re.I),
        "active_twiss_commands": re.findall(r"^\s*TWISS\s*,[^;]*;", source, re.I | re.M),
        "supplied_header": head, "supplied_rows": len(supplied),
        "generated_twiss_rows": len(generated), "generated_survey_rows": len(survey[2]),
        "all_row_names_and_keywords_match": True,
        "maximum_absolute_row_differences": maximum_differences,
        "header_differences": {key: {"supplied": value, "generated": generated_head.get(key)}
                               for key, value in head.items() if value != generated_head.get(key)},
        "runs": runs,
    }
    if kind == "dabs":
        new = {row["NAME"]: row for row in supplied}
        receipt["selected_dabs_rows"] = {name: new[name] for name in (
            "FULL$START", "M0_DECMAG", "M1_DECMAG", "Q_DQ901", "Q_DQ909",
            "Q_DQ933", "M1_DQ933", "M0_DHDA01", "D_DHDA01", "Q_DQDA01",
            "Q_DQDA02", "W_DMACE", "M0_M4DABS", "M1_M4DABS", "FULL$END")}

ring_rows = supplied_by_model["dr"]
extraction_rows = supplied_by_model["dabs"]
optics_columns = ("BETX", "BETY", "ALFX", "ALFY", "DX", "DY", "DPX", "DPY")
assert all(ring_rows[0][key] == ring_rows[-1][key] for key in optics_columns)
shared_name_prefix = 0
for ring_row, extraction_row in zip(ring_rows, extraction_rows):
    if ring_row["NAME"] != extraction_row["NAME"]:
        break
    shared_name_prefix += 1
ring_lookup = {row["NAME"]: row for row in ring_rows}
shared_quads = {}
for number in (303, 302, 301, 202):
    name = f"Q_DQ{number}"
    fields = ("S", "L", "TILT", "K1L")
    assert all(ring_lookup[name][key] == new[name][key] for key in fields)
    shared_quads[name] = {key: new[name][key] for key in fields}
receipt["ring_to_extraction_connection"] = {
    "shared_origin": "S=0, one millimetre upstream of Q303 entrance",
    "ring_start_and_end_optics_identical_at_printed_precision": True,
    "ring_start_optics": {key: ring_rows[0][key] for key in optics_columns},
    "extraction_hardcoded_start_optics": {key: extraction_rows[0][key] for key in optics_columns},
    "shared_initial_quads_with_identical_positions_lengths_rolls_and_k1l": shared_quads,
    "matching_initial_row_name_count": shared_name_prefix,
    "first_different_row_name": {"dr": ring_rows[shared_name_prefix]["NAME"],
                                 "dabs": extraction_rows[shared_name_prefix]["NAME"]},
    "first_septum": {kind: next(row for row in rows if row["NAME"] == "D_DESSA")
                     for kind, rows in supplied_by_model.items()},
    "automatic_cross_file_link": False,
    "interpretation": (
        "Separate standalone jobs with a shared Q303 origin and consistent upstream quad settings. "
        "Extraction starts with rounded ring periodic optics supplied as explicit constants. "
        "Septa and combined-function representations follow the extracted reference trajectory. "
        "No automatic Twiss read or particle handoff connects the jobs."
    ),
}

comparison = []
for family, path in (
    ("2018", REFERENCES / "2018_m4_twiss.tfs"),
    ("2024", REFERENCES / "2024_m4_twiss.tfs"),
):
    _, _, old = read_tfs(path)
    old_markers = {row["NAME"].upper(): row for row in old if row["KEYWORD"] == "MARKER"}
    for number in range(901, 934):
        name = f"Q{number}"
        pieces = [row for row in old if row["NAME"].upper() == name
                  and row["KEYWORD"] == "QUADRUPOLE"]
        assert len(pieces) == 2
        before_k1l = sum(row["K1L"] for row in pieces)
        before_length = sum(row["L"] for row in pieces)
        after = new[f"Q_DQ{number}"]
        after_k1l = after["K1L"] * math.cos(2 * after["TILT"])
        comparison.append({
            "family": family, "name": name,
            "old_signed_k1l_per_m": before_k1l,
            "new_signed_k1l_per_m": after_k1l,
            "relative_change_percent": 100 * (after_k1l / before_k1l - 1),
            "old_length_m": before_length, "new_length_m": after["L"],
            "old_center_s_m": old_markers[f"IQ{number}"]["S"],
            "new_center_s_m": after["S"],
            "center_difference_after_cmag_origin_m": (
                after["S"] - new["M0_DECMAG"]["S"] - old_markers[f"IQ{number}"]["S"]),
        })
receipt["comparison_method"] = (
    "Sum the two old half-quadrupole K1L values; compare with new K1L*cos(2*TILT). "
    "Use old IQnnn center markers, avoiding 2018 CENTRE vs 2024 exit convention. "
    "Align origins at new CMAG entrance; new physical lengths and old effective lengths differ."
)
receipt["quadrupole_comparison"] = comparison
receipt["limitations"] = [
    "Author-provided provenance describes DR settings as present and extraction settings as proposed future corrected settings.",
    "September 29 source label and output timestamps do not date every inherited component or prove operational approval.",
    "No full production-target sequence is supplied; tail selection and splice require a separate study.",
    "No particle tracking, space-charge, loss, scattering, or extinction calculation was performed.",
    "No remote machine-setting or approved global-survey comparison was performed.",
    "Header coupling diagnostics differ between MAD-X versions; only row-level numeric equivalence is established.",
]
(AUDIT / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
with (AUDIT / "quadrupole_comparison.csv").open("w", newline="") as file:
    writer = csv.DictWriter(file, fieldnames=list(comparison[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(comparison)
print(json.dumps({"supplied_files_verified": len(receipt["supplied_files"]),
                  "rows": {kind: model["supplied_rows"] for kind, model in receipt["models"].items()},
                  "all_row_numeric_fields_match": True,
                  "quadrupole_comparisons": len(comparison)}, indent=2))
