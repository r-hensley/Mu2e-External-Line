"""Deterministic hashes and runtime metadata for reproducible simulation files.

The helpers in this module deliberately avoid importing the tracking modules.
That keeps provenance collection inexpensive and prevents circular imports when
particle generators and loaders need to fingerprint their realized arrays.
"""

from __future__ import annotations

from hashlib import sha256
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
import platform
import subprocess
import sys
from typing import Iterable

import numpy as np


_DEPENDENCY_DISTRIBUTIONS = (
    "h5py",
    "matplotlib",
    "numpy",
    "xobjects",
    "xpart",
    "xdeps",
    "xsuite",
    "xtrack",
)


def _update_length_prefixed(digest, value: bytes) -> None:
    """Add one unambiguously delimited byte field to a SHA-256 digest."""
    digest.update(len(value).to_bytes(8, byteorder="big", signed=False))
    digest.update(value)


def array_sha256(*arrays: np.ndarray) -> str:
    """Hash NumPy array values together with each array's dtype and shape.

    Numeric arrays are converted to contiguous little-endian storage before
    hashing, so a digest does not depend on host byte order. Object arrays are
    rejected because their raw bytes contain process-specific pointers rather
    than stable values.
    """
    digest = sha256(b"xsuite-recycler-array-sha256-v1\0")
    for index, value in enumerate(arrays):
        array = np.asarray(value)
        if array.dtype.hasobject:
            raise TypeError("Object arrays cannot be hashed deterministically")

        # Use an explicit byte order for multi-byte numeric values. Dtypes such
        # as int8 have byte order '|' and need no conversion.
        canonical_dtype = (
            array.dtype
            if array.dtype.byteorder == "|"
            else array.dtype.newbyteorder("<")
        )
        canonical = np.ascontiguousarray(array.astype(canonical_dtype, copy=False))
        header = json.dumps(
            {
                "index": index,
                "dtype": canonical.dtype.str,
                "shape": list(canonical.shape),
            },
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        _update_length_prefixed(digest, header)
        _update_length_prefixed(digest, canonical.tobytes(order="C"))
    return digest.hexdigest()


def _project_source_files(project_root: Path) -> list[Path]:
    """Return sorted runtime source/configuration files included in the digest."""
    files = [
        path
        for path in (project_root / "main.py", project_root / "pyproject.toml")
        if path.is_file()
    ]
    package_dir = project_root / "xsuite_recycler"
    if package_dir.is_dir():
        files.extend(package_dir.rglob("*.py"))
    return sorted(files, key=lambda path: path.relative_to(project_root).as_posix())


def source_tree_provenance(project_root: str | Path | None = None) -> dict[str, object]:
    """Hash the Python source tree and return its included relative paths.

    Tests, outputs, documentation, and large particle inputs are intentionally
    excluded: they cannot change the executable tracking model. File names and
    contents are length-delimited in the digest to prevent path/content
    boundary ambiguities.
    """
    root = (
        Path(__file__).resolve().parents[1]
        if project_root is None
        else Path(project_root).resolve()
    )
    files = _project_source_files(root)
    digest = sha256(b"xsuite-recycler-source-tree-sha256-v1\0")
    relative_paths: list[str] = []
    for path in files:
        relative = path.relative_to(root).as_posix()
        relative_paths.append(relative)
        _update_length_prefixed(digest, relative.encode("utf-8"))
        _update_length_prefixed(digest, path.read_bytes())
    return {
        "project_root": str(root),
        "source_tree_sha256": digest.hexdigest(),
        "source_files": relative_paths,
    }


def _run_git(project_root: Path, arguments: Iterable[str]) -> subprocess.CompletedProcess[str]:
    """Run one read-only Git query without invoking a shell."""
    return subprocess.run(
        ["git", "-C", str(project_root), *arguments],
        check=False,
        capture_output=True,
        text=True,
    )


def git_provenance(project_root: str | Path | None = None) -> dict[str, object]:
    """Return the enclosing Git revision and project-scoped dirty state.

    A project can still be fingerprinted when Git is absent or the directory is
    not in a repository. In that case the revision and dirty state are ``None``
    and ``available`` is false; the deterministic source-tree digest remains
    available separately.
    """
    root = (
        Path(__file__).resolve().parents[1]
        if project_root is None
        else Path(project_root).resolve()
    )
    try:
        top_level_result = _run_git(root, ("rev-parse", "--show-toplevel"))
    except OSError:
        return {
            "available": False,
            "revision": None,
            "dirty": None,
            "repository_root": None,
            "project_path": None,
            "status_porcelain": None,
        }
    if top_level_result.returncode != 0:
        return {
            "available": False,
            "revision": None,
            "dirty": None,
            "repository_root": None,
            "project_path": None,
            "status_porcelain": None,
        }

    repository_root = Path(top_level_result.stdout.strip()).resolve()
    try:
        project_path = root.relative_to(repository_root).as_posix()
    except ValueError:
        project_path = root.as_posix()
    revision_result = _run_git(root, ("rev-parse", "HEAD"))
    status_result = _run_git(
        repository_root,
        (
            "status",
            "--porcelain=v1",
            "--untracked-files=normal",
            "--",
            project_path,
        ),
    )
    revision = (
        revision_result.stdout.strip() if revision_result.returncode == 0 else None
    )
    status = status_result.stdout.rstrip("\n") if status_result.returncode == 0 else None
    return {
        "available": revision is not None and status is not None,
        "revision": revision,
        "dirty": bool(status) if status is not None else None,
        "repository_root": str(repository_root),
        "project_path": project_path,
        "status_porcelain": status,
    }


def runtime_provenance() -> dict[str, object]:
    """Return Python/platform details and installed simulation dependencies."""
    dependencies: dict[str, str | None] = {}
    for distribution in _DEPENDENCY_DISTRIBUTIONS:
        try:
            dependencies[distribution] = version(distribution)
        except PackageNotFoundError:
            dependencies[distribution] = None
    return {
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "dependencies": dependencies,
    }


def collect_run_provenance(
    project_root: str | Path | None = None,
) -> dict[str, object]:
    """Collect source-tree, Git, and runtime metadata at simulation start."""
    return {
        "source_tree": source_tree_provenance(project_root),
        "git": git_provenance(project_root),
        "runtime": runtime_provenance(),
    }
