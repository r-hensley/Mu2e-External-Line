"""Unit tests for deterministic run-provenance helpers."""

import numpy as np
import pytest

import xsuite_recycler.provenance as provenance
from xsuite_recycler.provenance import array_sha256, source_tree_provenance


def test_array_hash_is_endian_stable_and_schema_sensitive():
    """Equal numeric arrays hash alike across endian representations."""
    native = np.asarray([1.25, -3.5], dtype="<f8")
    big_endian = np.asarray([1.25, -3.5], dtype=">f8")

    assert array_sha256(native) == array_sha256(big_endian)
    assert array_sha256(native) != array_sha256(native.astype("f4"))
    assert array_sha256(native) != array_sha256(native.reshape(1, 2))
    assert array_sha256(native, native) != array_sha256(native)


def test_array_hash_rejects_process_specific_object_storage():
    """Object arrays are not silently fingerprinted by pointer bytes."""
    with pytest.raises(TypeError, match="Object arrays"):
        array_sha256(np.asarray([{"not": "stable"}], dtype=object))


def test_source_tree_hash_covers_runtime_files_but_not_tests(tmp_path):
    """Runtime source edits change the digest while test-only edits do not."""
    package = tmp_path / "xsuite_recycler"
    package.mkdir()
    (tmp_path / "main.py").write_text("print('entry')\n")
    (tmp_path / "pyproject.toml").write_text("[project]\nname='demo'\n")
    module = package / "model.py"
    module.write_text("VALUE = 1\n")

    first = source_tree_provenance(tmp_path)
    assert first["source_files"] == [
        "main.py",
        "pyproject.toml",
        "xsuite_recycler/model.py",
    ]

    tests = tmp_path / "tests"
    tests.mkdir()
    (tests / "test_model.py").write_text("assert True\n")
    assert source_tree_provenance(tmp_path)["source_tree_sha256"] == first[
        "source_tree_sha256"
    ]

    module.write_text("VALUE = 2\n")
    assert source_tree_provenance(tmp_path)["source_tree_sha256"] != first[
        "source_tree_sha256"
    ]


def test_git_provenance_has_explicit_unavailable_state(monkeypatch, tmp_path):
    """A missing Git executable leaves a useful, JSON-safe result."""

    def missing_git(project_root, arguments):
        """Stand in for an unavailable Git executable."""
        raise OSError("git not found")

    monkeypatch.setattr(provenance, "_run_git", missing_git)
    result = provenance.git_provenance(tmp_path)

    assert result == {
        "available": False,
        "revision": None,
        "dirty": None,
        "repository_root": None,
        "project_path": None,
        "status_porcelain": None,
    }


def test_runtime_provenance_reports_tracking_dependencies():
    """The runtime record includes Python and the direct tracking libraries."""
    result = provenance.runtime_provenance()

    assert result["python_version"]
    assert result["python_executable"]
    assert result["dependencies"]["numpy"] == np.__version__
    assert result["dependencies"]["h5py"]
    assert result["dependencies"]["xtrack"]
