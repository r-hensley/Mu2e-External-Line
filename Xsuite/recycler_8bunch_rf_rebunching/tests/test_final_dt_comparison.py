import json

import h5py
import numpy as np
import pytest

from xsuite_recycler.config import C_LIGHT_M_PER_S
from xsuite_recycler.final_dt_comparison import (
    compare_final_dt,
    load_final_dt_data,
    write_final_dt_comparison_plot,
    write_final_dt_comparison_report,
)


def _write_result(
    path,
    dt_ns,
    labels,
    *,
    macro_weight=1.0,
    physical_weight=10.0,
    source_rows=None,
    rf25_frequency_hz=2.5e6,
    rf25_initial_v=3000.0,
    gamma_transition=20.0,
    input_sha256="abc123",
):
    dt_ns = np.asarray(dt_ns, dtype=float)
    labels = np.asarray(labels)
    beta0 = 0.9
    if source_rows is None:
        source_rows = int(dt_ns.size * macro_weight)
    with h5py.File(path, "w") as h5:
        h5.attrs.update(
            {
                "format_version": 1,
                "model": "RF-only test",
                "xtrack_version": "test",
                "source_blond_commit": "abc",
                "source_blond_driver": "driver.py",
                "coordinate_convention": "zeta=-beta*c*dt",
                "energy_convention": "ptau=dpop/beta",
                "n_turns": 20,
                "end_state_time_s": 2.0e-4,
                "last_kick_time_s": 1.9e-4,
                "first_rf25_kick_turn": 5,
                "record_every": 10,
                "sample_particles": dt_ns.size,
                "source_rows": source_rows,
                "macro_equivalent_weight": macro_weight,
                "physical_proton_weight": physical_weight,
                "last_kick_rf53_v": 0.0,
                "last_kick_rf25_v": 80_000.0,
                "input_source_kind": "file",
                "input_sha256": input_sha256,
            }
        )
        machine = h5.create_group("machine")
        machine.attrs.update(
            {
                "beta0": beta0,
                "p0c_ev": 8.0e9,
                "revolution_period_s": 1.0e-5,
                "harmonic_25": 28,
                "rf25_frequency_hz": rf25_frequency_hz,
                "gamma_transition": gamma_transition,
                "momentum_compaction_factor": 1.0 / gamma_transition**2,
                "slip_factor": 1.0 / gamma_transition**2 - 0.01,
            }
        )
        ramp = h5.create_group("rf_program")
        ramp.attrs["rf25_initial_v"] = rf25_initial_v
        ramp.attrs["rf25_final_v"] = 80_000.0
        profiles = h5.create_group("profiles")
        profiles.create_dataset("turn", data=np.asarray([0, 10, 20]))
        profiles.create_dataset(
            "time_edges_ns", data=np.asarray([-1000.0, 1400.0, 3800.0])
        )
        profiles.create_dataset(
            "counts_macro_equivalent",
            data=np.full((3, 2), source_rows / 2.0),
        )
        profiles.create_dataset("outside_macro_equivalent", data=np.zeros(3))
        final = h5.create_group("final_particles")
        final.create_dataset(
            "zeta_m", data=-dt_ns * 1e-9 * beta0 * C_LIGHT_M_PER_S
        )
        final.create_dataset("ptau", data=np.zeros(dt_ns.size))
        final.create_dataset("state", data=np.ones(dt_ns.size, dtype=np.int64))
        final.create_dataset("bunch_index", data=labels)


def test_identical_final_dt_is_json_safe_and_has_zero_shape_distances(tmp_path):
    reference = tmp_path / "reference.h5"
    candidate = tmp_path / "candidate.h5"
    dt = [-100.0, -90.0, 300.0, 310.0]
    labels = [0, 0, 1, 1]
    _write_result(reference, dt, labels)
    _write_result(candidate, dt, labels)

    loaded = load_final_dt_data(reference)
    assert loaded.dt_ns == pytest.approx(dt)
    report = compare_final_dt(
        reference,
        candidate,
        absolute_bins=480,
        absolute_range_ns=(-1000.0, 3800.0),
    )

    absolute = report["comparison"]["absolute_histogram"]
    assert absolute["wasserstein_1_ns"] == pytest.approx(0.0)
    assert absolute["jensen_shannon_distance"] == pytest.approx(0.0)
    assert absolute["histogram_ks_distance"] == pytest.approx(0.0)
    assert absolute["total_variation_distance"] == pytest.approx(0.0)
    assert report["comparison"]["common_bunch_labels"] == ["0", "1"]
    json.dumps(report, allow_nan=False)


def test_reference_center_is_shared_and_unwrapped_bucket_migration_is_visible(tmp_path):
    reference = tmp_path / "reference.h5"
    candidate = tmp_path / "candidate.h5"
    # At 2.5 MHz the bucket period is 400 ns.  The candidate has migrated by
    # one full bucket; modulo wrapping would incorrectly erase this difference.
    _write_result(reference, [390.0, 410.0], [0, 0])
    _write_result(candidate, [790.0, 810.0], [0, 0])

    report = compare_final_dt(
        reference,
        candidate,
        absolute_bins=480,
        absolute_range_ns=(-1000.0, 3800.0),
        bucket_relative_bins=240,
        bucket_relative_range_ns=(-600.0, 600.0),
    )
    ref_bunch = report["reference"]["per_bunch"]["0"]
    cand_bunch = report["candidate"]["per_bunch"]["0"]
    differences = report["comparison"]["per_bunch"]["0"]

    assert ref_bunch["reference_h28_bucket_center_ns"] == pytest.approx(400.0)
    assert cand_bunch["reference_h28_bucket_center_ns"] == pytest.approx(400.0)
    assert ref_bunch["bucket_relative_dt"]["mean_ns"] == pytest.approx(0.0)
    assert cand_bunch["bucket_relative_dt"]["mean_ns"] == pytest.approx(400.0)
    assert differences["bucket_relative_summary_differences"][
        "absolute_difference"
    ]["mean_ns"] == pytest.approx(400.0)
    assert ref_bunch["final_h28_bucket_occupancy"]["1"]["raw_count"] == 2
    assert cand_bunch["final_h28_bucket_occupancy"]["2"]["raw_count"] == 2
    transition = differences["final_h28_bucket_occupancy_differences"]
    assert transition["1"]["candidate_minus_reference_raw_count"] == -2
    assert transition["2"]["candidate_minus_reference_raw_count"] == 2
    assert report["comparison"]["absolute_histogram"][
        "wasserstein_1_ns"
    ] == pytest.approx(400.0)
    assert differences["bucket_relative_histogram"][
        "wasserstein_1_ns"
    ] == pytest.approx(400.0)
    assert report["comparison"]["absolute_histogram"][
        "histogram_ks_distance"
    ] == pytest.approx(1.0)
    assert report["comparison"]["absolute_histogram"][
        "total_variation_distance"
    ] == pytest.approx(1.0)
    assert cand_bunch["bucket_relative_tails"]["abs_dt_gt_120_ns"][
        "raw_count"
    ] == 2
    ref_retention = report["reference"][
        "initial_lineage_reference_bucket_retention"
    ]
    cand_retention = report["candidate"][
        "initial_lineage_reference_bucket_retention"
    ]
    retention_difference = report["comparison"][
        "initial_lineage_reference_bucket_retention_differences"
    ]
    assert ref_retention["raw_outside_reference_bucket"] == 0
    assert ref_retention["fraction_outside_reference_bucket"] == 0.0
    assert cand_retention["raw_outside_reference_bucket"] == 2
    assert cand_retention["fraction_outside_reference_bucket"] == 1.0
    assert retention_difference[
        "candidate_minus_reference_fraction_outside_reference_bucket"
    ] == 1.0


def test_unwrapped_125ns_window_reports_weighted_counts_by_initial_lineage(
    tmp_path,
):
    reference = tmp_path / "reference.h5"
    candidate = tmp_path / "candidate.h5"
    # Historical medians snap to the h=28 centers at 400 and 800 ns.  Values
    # exactly at +/-125 ns remain in-time; no modulo operation is allowed.
    _write_result(
        reference,
        [390.0, 400.0, 410.0, 790.0, 800.0, 810.0],
        [0, 0, 0, 1, 1, 1],
        macro_weight=2.0,
        physical_weight=20.0,
    )
    _write_result(
        candidate,
        [274.0, 275.0, 525.0, 674.0, 675.0, 1200.0],
        [0, 0, 0, 1, 1, 1],
        macro_weight=2.0,
        physical_weight=20.0,
    )

    report = compare_final_dt(reference, candidate)
    reference_window = report["reference"][
        "unwrapped_initial_lineage_125ns_window"
    ]
    candidate_window = report["candidate"][
        "unwrapped_initial_lineage_125ns_window"
    ]

    assert candidate_window["window_half_width_ns"] == 125.0
    assert candidate_window["outside_rule"] == (
        "abs(dt_rel_ns) > window_half_width_ns"
    )
    assert candidate_window["modulo_wrapping"] is False
    assert candidate_window["separate_from_reference_bucket_migration"] is True
    assert candidate_window["reference_bucket_migration_metric_path"] == (
        "initial_lineage_reference_bucket_retention"
    )
    assert reference_window["global"]["raw_outside_count"] == 0

    global_counts = candidate_window["global"]
    assert global_counts["raw_total_count"] == 6
    assert global_counts["raw_outside_count"] == 3
    assert global_counts["raw_outside_fraction"] == pytest.approx(0.5)
    assert global_counts["macro_equivalent_total_count"] == pytest.approx(12.0)
    assert global_counts["macro_equivalent_outside_count"] == pytest.approx(6.0)
    assert global_counts["macro_equivalent_outside_fraction"] == pytest.approx(0.5)
    assert global_counts["physical_proton_total_count"] == pytest.approx(120.0)
    assert global_counts["physical_proton_outside_count"] == pytest.approx(60.0)
    assert global_counts["physical_proton_outside_fraction"] == pytest.approx(0.5)

    # Folding to the nearest h=28 center instead treats the 1200 ns migrant as
    # locally centered. Only the two +/-126 ns particles remain out of time;
    # exact +/-125 ns values stay inside.
    local_window = report["candidate"][
        "nearest_h28_bucket_folded_125ns_local_window"
    ]
    assert local_window["modulo_wrapping"] is True
    assert local_window["initial_lineage_used_for_center"] is False
    assert local_window["global"]["raw_outside_count"] == 2
    assert local_window["global"]["raw_outside_fraction"] == pytest.approx(1 / 3)

    lineage_zero = candidate_window["per_initial_lineage"]["0"]
    lineage_one = candidate_window["per_initial_lineage"]["1"]
    assert lineage_zero["reference_h28_bucket_center_ns"] == pytest.approx(400.0)
    assert lineage_zero["raw_total_count"] == 3
    assert lineage_zero["raw_outside_count"] == 1
    assert lineage_zero["raw_outside_fraction"] == pytest.approx(1.0 / 3.0)
    assert lineage_one["reference_h28_bucket_center_ns"] == pytest.approx(800.0)
    assert lineage_one["raw_total_count"] == 3
    assert lineage_one["raw_outside_count"] == 2
    assert lineage_one["raw_outside_fraction"] == pytest.approx(2.0 / 3.0)

    # The unwrapped source window and nearest-bucket migration answer different
    # questions: three particles are outside +/-125 ns, but only the 1200 ns
    # particle has moved out of its lineage's reference h=28 bucket.
    migration = report["candidate"][
        "initial_lineage_reference_bucket_retention"
    ]
    assert migration["raw_outside_reference_bucket"] == 1
    assert migration["fraction_outside_reference_bucket"] == pytest.approx(1.0 / 6.0)

    differences = report["comparison"][
        "unwrapped_initial_lineage_125ns_window_differences"
    ]
    assert differences["common_initial_lineage_labels"] == ["0", "1"]
    assert differences["global"][
        "candidate_minus_reference_raw_outside_count"
    ] == 3
    assert differences["global"][
        "candidate_minus_reference_macro_equivalent_outside_count"
    ] == pytest.approx(6.0)
    assert differences["global"][
        "candidate_minus_reference_physical_proton_outside_count"
    ] == pytest.approx(60.0)
    assert differences["global"][
        "candidate_minus_reference_raw_outside_fraction"
    ] == pytest.approx(0.5)
    assert differences["per_initial_lineage"]["1"][
        "candidate_minus_reference_raw_outside_fraction"
    ] == pytest.approx(2.0 / 3.0)
    json.dumps(report, allow_nan=False)


def test_unequal_samples_and_weights_compare_normalized_shapes(tmp_path):
    reference = tmp_path / "reference.h5"
    candidate = tmp_path / "candidate.h5"
    _write_result(
        reference,
        [0.0, 10.0],
        [0, 0],
        macro_weight=4.0,
        physical_weight=40.0,
        source_rows=8,
    )
    _write_result(
        candidate,
        [0.0, 0.0, 10.0, 10.0],
        [0, 0, 0, 0],
        macro_weight=2.0,
        physical_weight=20.0,
        source_rows=8,
    )

    report = compare_final_dt(
        reference,
        candidate,
        absolute_bins=400,
        absolute_range_ns=(-100.0, 100.0),
    )
    metric = report["comparison"]["absolute_histogram"]
    assert metric["wasserstein_1_ns"] == pytest.approx(0.0)
    assert metric["jensen_shannon_distance"] == pytest.approx(0.0)
    assert report["reference"]["effective_counts"][
        "raw_effective_sample_size"
    ] == 2
    assert report["candidate"]["effective_counts"][
        "raw_effective_sample_size"
    ] == 4
    assert report["reference"]["effective_counts"][
        "surviving_finite_macro_equivalent_count"
    ] == pytest.approx(8.0)
    assert report["candidate"]["effective_counts"][
        "surviving_finite_macro_equivalent_count"
    ] == pytest.approx(8.0)


def test_empirical_wasserstein_retains_sub_bin_displacement(tmp_path):
    reference = tmp_path / "reference.h5"
    candidate = tmp_path / "candidate.h5"
    _write_result(reference, [0.1], [0])
    _write_result(candidate, [0.9], [0])

    report = compare_final_dt(
        reference,
        candidate,
        absolute_bins=1,
        absolute_range_ns=(0.0, 2.0),
    )
    comparison = report["comparison"]
    assert comparison["absolute_histogram"]["wasserstein_1_ns"] == 0.0
    assert comparison["absolute_empirical_wasserstein_1_ns"] == pytest.approx(0.8)
    assert comparison["per_bunch"]["0"][
        "absolute_empirical_wasserstein_1_ns"
    ] == pytest.approx(0.8)


def test_strict_validation_rejects_incompatible_rf_program(tmp_path):
    reference = tmp_path / "reference.h5"
    candidate = tmp_path / "candidate.h5"
    _write_result(reference, [0.0, 1.0], [0, 0])
    _write_result(candidate, [0.0, 1.0], [0, 0], rf25_initial_v=5000.0)

    with pytest.raises(ValueError, match="rf_program attribute rf25_initial_v"):
        compare_final_dt(reference, candidate)


def test_controlled_transition_gamma_comparison_is_narrow_and_explicit(tmp_path):
    reference = tmp_path / "keegan.h5"
    candidate = tmp_path / "werkema.h5"
    _write_result(reference, [390.0, 410.0], [0, 0], gamma_transition=20.0)
    _write_result(candidate, [385.0, 415.0], [0, 0], gamma_transition=21.6)

    with pytest.raises(ValueError, match="machine attribute gamma_transition"):
        compare_final_dt(reference, candidate)

    report = compare_final_dt(
        reference,
        candidate,
        allow_gamma_transition_difference=True,
    )
    validation = report["validation"]
    assert not validation["strict_same_physics_and_endpoint"]
    assert validation["controlled_transition_gamma_comparison"]
    assert validation["same_full_file_backed_input_required"]
    assert not validation["sampling_sizes_and_weights_may_differ"]
    assert validation["allowed_machine_attribute_differences"] == [
        "gamma_transition",
        "momentum_compaction_factor",
        "slip_factor",
    ]
    gamma = report["machine_parameter_comparison"]["gamma_transition"]
    assert gamma["reference"] == pytest.approx(20.0)
    assert gamma["candidate"] == pytest.approx(21.6)
    assert gamma["candidate_minus_reference"] == pytest.approx(1.6)
    assert gamma["candidate_to_reference_ratio"] == pytest.approx(1.08)
    periods = report["derived_final_voltage_synchrotron_period_comparison"]
    assert periods["candidate_ms"] < periods["reference_ms"]

    # The narrow mode must not mask an unrelated machine or input change.
    with h5py.File(candidate, "r+") as h5:
        h5["machine"].attrs["p0c_ev"] = 8.1e9
    with pytest.raises(ValueError, match="machine attribute p0c_ev"):
        compare_final_dt(
            reference,
            candidate,
            allow_gamma_transition_difference=True,
        )


def test_controlled_gamma_pair_reports_all_folded_125ns_transitions(tmp_path):
    reference = tmp_path / "keegan.h5"
    candidate = tmp_path / "werkema.h5"
    # At T=400 ns, +/-125 ns is inside and +/-126 ns is outside. Construct one
    # particle in each paired classification cell.
    _write_result(
        reference,
        [0.0, 0.0, 126.0, 126.0],
        [0, 0, 0, 0],
        gamma_transition=20.0,
    )
    _write_result(
        candidate,
        [0.0, 126.0, 0.0, 130.0],
        [0, 0, 0, 0],
        gamma_transition=21.6,
    )

    report = compare_final_dt(
        reference,
        candidate,
        allow_gamma_transition_difference=True,
    )
    paired = report["comparison"][
        "paired_nearest_h28_bucket_folded_125ns_transitions"
    ]
    assert paired["available"]
    raw = paired["global"]["counts"]["raw"]
    assert raw == {
        "paired_total": 4,
        "inside_to_inside": 1,
        "inside_to_outside": 1,
        "outside_to_inside": 1,
        "outside_to_outside": 1,
        "changed": 2,
        "net_outside": 0,
    }
    assert paired["per_initial_lineage"]["0"]["counts"]["raw"] == raw
    ref_outside = report["reference"][
        "nearest_h28_bucket_folded_125ns_local_window"
    ]["global"]["raw_outside_count"]
    cand_outside = report["candidate"][
        "nearest_h28_bucket_folded_125ns_local_window"
    ]["global"]["raw_outside_count"]
    assert raw["net_outside"] == cand_outside - ref_outside


def test_transition_gamma_comparison_requires_identical_full_input(tmp_path):
    reference = tmp_path / "keegan.h5"
    candidate = tmp_path / "werkema.h5"
    _write_result(reference, [0.0, 1.0], [0, 0], gamma_transition=20.0)
    _write_result(
        candidate,
        [0.0, 1.0],
        [0, 0],
        gamma_transition=21.6,
        input_sha256="different",
    )

    with pytest.raises(ValueError, match="matching root attribute input_sha256"):
        compare_final_dt(
            reference,
            candidate,
            allow_gamma_transition_difference=True,
        )


def test_bunch_labels_are_required(tmp_path):
    reference = tmp_path / "reference.h5"
    candidate = tmp_path / "candidate.h5"
    _write_result(reference, [0.0, 1.0], [0, 0])
    _write_result(candidate, [0.0, 1.0], [0, 0])
    with h5py.File(candidate, "r+") as h5:
        del h5["final_particles/bunch_index"]

    with pytest.raises(ValueError, match="bunch_index is required"):
        compare_final_dt(reference, candidate)


def test_writers_create_valid_report_and_plot(tmp_path):
    reference = tmp_path / "reference.h5"
    candidate = tmp_path / "candidate.h5"
    _write_result(reference, [390.0, 410.0, 790.0, 810.0], [0, 0, 1, 1])
    _write_result(candidate, [385.0, 415.0, 795.0, 805.0], [0, 0, 1, 1])

    report_path = write_final_dt_comparison_report(
        reference,
        candidate,
        tmp_path / "comparison.json",
        absolute_bins=240,
    )
    plot_path = write_final_dt_comparison_plot(
        reference,
        candidate,
        tmp_path / "comparison.png",
        absolute_bins=240,
    )

    with report_path.open(encoding="utf-8") as stream:
        report = json.load(stream)
    assert report["schema_version"] == 1
    assert plot_path.stat().st_size > 10_000
