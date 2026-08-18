import json

import h5py
import numpy as np
import pytest

from xsuite_recycler.result_comparison import compare_results, load_result_data


def _write_result(
    path,
    *,
    turns=(0, 10, 20),
    edges=(0.0, 1.0, 2.0),
    counts=((1.0, 1.0), (2.0, 0.0), (0.0, 2.0)),
    outside=(0.0, 0.0, 0.0),
    include_final=True,
):
    counts = np.asarray(counts, dtype=float)
    outside = np.asarray(outside, dtype=float)
    with h5py.File(path, "w") as h5:
        h5.attrs["source_rows"] = 2.0
        h5.attrs["sample_particles"] = 2
        machine = h5.create_group("machine")
        machine.attrs["beta0"] = 0.9
        machine.attrs["p0c_ev"] = 8.0e9
        machine.attrs["rf25_frequency_hz"] = 2.5e6
        profiles = h5.create_group("profiles")
        profiles.create_dataset("turn", data=np.asarray(turns, dtype=np.int64))
        profiles.create_dataset("time_edges_ns", data=np.asarray(edges, dtype=float))
        profiles.create_dataset("counts_macro_equivalent", data=counts)
        profiles.create_dataset("outside_macro_equivalent", data=outside)
        if include_final:
            final = h5.create_group("final_particles")
            final.create_dataset("zeta_m", data=np.asarray([0.0, -27.0]))
            final.create_dataset("ptau", data=np.asarray([1e-4, -1e-4]))
            final.create_dataset("state", data=np.asarray([1, 1]))
            final.create_dataset("bunch_index", data=np.asarray([0, 1]))


def _add_strict_metadata(path, *, rf25_initial_v=3000.0):
    with h5py.File(path, "r+") as h5:
        h5.attrs.update(
            {
                "format_version": 1,
                "model": "RF-only test",
                "xtrack_version": "test",
                "source_blond_commit": "abc",
                "source_blond_driver": "driver.py",
                "coordinate_convention": "test-zeta",
                "energy_convention": "test-ptau",
                "n_turns": 20,
                "end_state_time_s": 2e-4,
                "last_kick_time_s": 1.9e-4,
                "first_rf25_kick_turn": 5,
                "record_every": 10,
                "sample_particles": 2,
                "macro_equivalent_weight": 1.0,
                "physical_proton_weight": 1.0,
                "last_kick_rf53_v": 0.0,
                "last_kick_rf25_v": 80_000.0,
            }
        )
        ramp = h5.create_group("rf_program")
        ramp.attrs["rf25_initial_v"] = rf25_initial_v


def test_identical_results_have_zero_distances_and_json_safe_output(tmp_path):
    left = tmp_path / "left.h5"
    right = tmp_path / "right.h5"
    _write_result(left)
    _write_result(right)

    loaded = load_result_data(left)
    assert loaded.profiles.counts_macro_equivalent.shape == (3, 2)

    comparison = compare_results(left, right)
    assert comparison["profiles"]["common_turns"] == [0, 10, 20]
    assert comparison["profiles"]["wasserstein_ns"] == [0.0, 0.0, 0.0]
    assert comparison["profiles"]["jensen_shannon_distance"] == [0.0, 0.0, 0.0]
    assert comparison["conservation"]["reference"]["conserved_to_writer_tolerance"]
    assert comparison["final_particles"]["reference"]["available"]
    json.dumps(comparison, allow_nan=False)


def test_intersects_turns_and_computes_known_disjoint_profile_distance(tmp_path):
    reference = tmp_path / "reference.h5"
    candidate = tmp_path / "candidate.h5"
    _write_result(reference)
    _write_result(
        candidate,
        turns=(0, 20, 30),
        counts=((1.0, 1.0), (2.0, 0.0), (1.0, 1.0)),
        outside=(0.0, 0.0, 0.0),
    )

    comparison = compare_results(reference, candidate)
    profiles = comparison["profiles"]
    assert profiles["common_turns"] == [0, 20]
    assert profiles["wasserstein_ns"] == pytest.approx([0.0, 1.0])
    assert profiles["jensen_shannon_distance"] == pytest.approx([0.0, 1.0])


def test_outside_mass_is_separate_from_in_window_shape_and_checked_for_conservation(tmp_path):
    reference = tmp_path / "reference.h5"
    candidate = tmp_path / "candidate.h5"
    _write_result(
        reference,
        counts=((1.0, 0.0),),
        turns=(0,),
        outside=(1.0,),
    )
    _write_result(
        candidate,
        counts=((1.0, 0.0),),
        turns=(0,),
        outside=(0.0,),
    )

    comparison = compare_results(reference, candidate)
    profiles = comparison["profiles"]
    assert profiles["wasserstein_ns"] == [0.0]
    assert profiles["jensen_shannon_distance"] == [0.0]
    assert profiles["jensen_shannon_with_outside_distance"][0] > 0
    assert profiles["reference_in_window_fraction"] == [0.5]
    assert profiles["candidate_in_window_fraction"] == [1.0]
    assert comparison["conservation"]["reference"]["conserved_to_writer_tolerance"]
    assert not comparison["conservation"]["candidate"]["conserved_to_writer_tolerance"]


def test_missing_final_group_is_reported_without_losing_profile_comparison(tmp_path):
    reference = tmp_path / "reference.h5"
    candidate = tmp_path / "candidate.h5"
    _write_result(reference, include_final=False)
    _write_result(candidate, include_final=False)

    comparison = compare_results(reference, candidate)
    assert comparison["profiles"]["n_common_turns"] == 3
    assert not comparison["final_particles"]["reference"]["available"]
    assert not comparison["final_particles"]["differences"]["available"]


def test_incompatible_profile_bins_fail_clearly(tmp_path):
    reference = tmp_path / "reference.h5"
    candidate = tmp_path / "candidate.h5"
    _write_result(reference)
    _write_result(candidate, edges=(0.0, 1.1, 2.0))

    with pytest.raises(ValueError, match="Incompatible profile time-bin edges"):
        compare_results(reference, candidate)


def test_strict_comparison_rejects_mismatched_rf_program(tmp_path):
    reference = tmp_path / "reference.h5"
    candidate = tmp_path / "candidate.h5"
    _write_result(reference)
    _write_result(candidate)
    _add_strict_metadata(reference, rf25_initial_v=3000.0)
    _add_strict_metadata(candidate, rf25_initial_v=5000.0)

    # Exploratory pairwise comparison remains available when explicitly non-strict.
    compare_results(reference, candidate)
    with pytest.raises(ValueError, match="rf_program attribute rf25_initial_v"):
        compare_results(reference, candidate, require_matching_run=True)
