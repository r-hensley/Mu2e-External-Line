"""Aggregate reduced-particle generator screens against resampling controls.

Generated-input runs are compared with a file-based historical reference.
Equal-size resamples of that same historical input provide a practical baseline
for discrepancies caused by finite tracking samples rather than the generator.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import numpy as np

from .result_comparison import compare_results


def _tail_summary(final: dict[str, Any], threshold_ns: int) -> dict[str, float]:
    """Combine per-bunch tail fractions using their usable particle counts."""

    entries = list(final["per_bunch"].values())
    total = sum(int(entry["n_alive_finite"]) for entry in entries)
    count = sum(
        float(entry[f"fraction_abs_time_gt_{threshold_ns}_ns"])
        * int(entry["n_alive_finite"])
        for entry in entries
    )
    return {
        "fraction": count / total,
        "sample_equivalent_count": count,
        "alive_finite_sample_particles": total,
    }


def _pair_summary(report: dict[str, Any]) -> dict[str, Any]:
    """Reduce a full result comparison to the metrics used by the RF screen."""

    profiles = report["profiles"]
    final = report["final_particles"]
    result: dict[str, Any] = {
        "candidate_path": report["candidate_path"],
        "run_metadata": report["run_metadata"],
        "profile": {
            "n_common_turns": profiles["n_common_turns"],
            "median_wasserstein_ns": profiles["wasserstein_summary_ns"]["median"],
            "mean_wasserstein_ns": profiles["wasserstein_summary_ns"]["mean"],
            "maximum_wasserstein_ns": profiles["wasserstein_summary_ns"]["maximum"],
            "median_jensen_shannon": profiles["jensen_shannon_summary"]["median"],
            "mean_jensen_shannon": profiles["jensen_shannon_summary"]["mean"],
        },
        "conservation": {
            "reference": report["conservation"]["reference"][
                "conserved_to_writer_tolerance"
            ],
            "candidate": report["conservation"]["candidate"][
                "conserved_to_writer_tolerance"
            ],
        },
    }
    if not (
        final["reference"].get("available")
        and final["candidate"].get("available")
        and final["differences"].get("available")
    ):
        result["final_particles"] = {"available": False}
        return result

    # Compare lineage centroids directly, but scale width differences by their
    # historical widths to obtain dimensionless relative errors.
    reference_bunches = final["reference"]["per_bunch"]
    candidate_bunches = final["candidate"]["per_bunch"]
    differences = final["differences"]["per_bunch"]
    labels = sorted(set(reference_bunches) & set(candidate_bunches), key=int)
    centroid_abs = np.asarray(
        [
            abs(differences[label]["candidate_minus_reference_mean_relative_time_ns"])
            for label in labels
        ]
    )
    time_width_relative = np.asarray(
        [
            differences[label]["candidate_minus_reference_rms_relative_time_ns"]
            / reference_bunches[label]["rms_relative_time_ns"]
            for label in labels
        ]
    )
    energy_width_relative = np.asarray(
        [
            differences[label]["candidate_minus_reference_rms_energy_offset_ev"]
            / reference_bunches[label]["rms_energy_offset_ev"]
            for label in labels
        ]
    )
    result["final_particles"] = {
        "available": True,
        "reference_population": {
            "n_total": final["reference"]["n_total"],
            "n_alive": final["reference"]["n_alive"],
            "n_alive_finite": final["reference"]["n_alive_finite"],
        },
        "candidate_population": {
            "n_total": final["candidate"]["n_total"],
            "n_alive": final["candidate"]["n_alive"],
            "n_alive_finite": final["candidate"]["n_alive_finite"],
        },
        "common_bunch_labels": labels,
        "maximum_absolute_centroid_difference_ns": float(np.max(centroid_abs)),
        "mean_absolute_centroid_difference_ns": float(np.mean(centroid_abs)),
        "maximum_absolute_rms_time_relative_difference": float(
            np.max(np.abs(time_width_relative))
        ),
        "mean_absolute_rms_time_relative_difference": float(
            np.mean(np.abs(time_width_relative))
        ),
        "maximum_absolute_rms_energy_relative_difference": float(
            np.max(np.abs(energy_width_relative))
        ),
        "mean_absolute_rms_energy_relative_difference": float(
            np.mean(np.abs(energy_width_relative))
        ),
        "reference_tails": {
            str(threshold): _tail_summary(final["reference"], threshold)
            for threshold in (80, 100, 120)
        },
        "candidate_tails": {
            str(threshold): _tail_summary(final["candidate"], threshold)
            for threshold in (80, 100, 120)
        },
    }
    return result


def _median_pair_metric(pairs: list[dict], section: str, key: str) -> float:
    """Return the median of one nested metric across comparison pairs."""

    values = [float(pair[section][key]) for pair in pairs]
    return float(np.median(values))


def _safe_ratio(value: float, baseline: float) -> float | None:
    """Form a baseline ratio, preserving the meaning of an exact zero baseline."""

    if baseline > 0:
        return value / baseline
    return 0.0 if value == 0 else None


def summarize_rf_screen(
    reference_path: str | Path,
    generated_paths: Sequence[str | Path],
    control_paths: Sequence[str | Path],
) -> dict[str, Any]:
    """Judge generated seeds against equal-size historical resampling controls.

    ``reference_path`` is the common full or high-statistics file-based run.
    Generated candidates are assessed relative to the median discrepancy seen
    in ``control_paths``.  The returned report records each comparison, the
    resampling baseline, explicit pass/fail criteria, and the screen's limits.
    """

    if not generated_paths:
        raise ValueError("At least one generated result is required")
    if not control_paths:
        raise ValueError("At least one historical resampling control is required")

    generated = [
        _pair_summary(
            compare_results(reference_path, path, require_matching_run=True)
        )
        for path in generated_paths
    ]
    controls = [
        _pair_summary(
            compare_results(reference_path, path, require_matching_run=True)
        )
        for path in control_paths
    ]
    if any(
        pair["run_metadata"]["reference_input_source_kind"] != "file"
        or pair["run_metadata"]["candidate_input_source_kind"] != "generated"
        for pair in generated
    ):
        raise ValueError(
            "Generated screen pairs require a file reference and generated candidate"
        )
    if any(
        pair["run_metadata"]["reference_input_source_kind"] != "file"
        or pair["run_metadata"]["candidate_input_source_kind"] != "file"
        for pair in controls
    ):
        raise ValueError(
            "Historical controls require file-based reference and candidate inputs"
        )
    if any(
        not pair["final_particles"].get("available")
        for pair in [*generated, *controls]
    ):
        raise ValueError("RF screening requires final-particle diagnostics in every run")
    reference_particles = int(
        generated[0]["final_particles"]["reference_population"]["n_total"]
    )
    if any(
        int(pair["final_particles"]["reference_population"]["n_total"])
        != reference_particles
        for pair in [*generated, *controls]
    ):
        raise ValueError("Every RF screen comparison must use the same reference size")
    # The median across control resamples estimates the discrepancy expected
    # solely from tracking a reduced random sample of the historical input.
    baseline = {
        "median_profile_wasserstein_ns": _median_pair_metric(
            controls, "profile", "median_wasserstein_ns"
        ),
        "median_profile_jensen_shannon": _median_pair_metric(
            controls, "profile", "median_jensen_shannon"
        ),
        "median_maximum_absolute_centroid_difference_ns": _median_pair_metric(
            controls,
            "final_particles",
            "maximum_absolute_centroid_difference_ns",
        ),
        "median_mean_absolute_centroid_difference_ns": _median_pair_metric(
            controls,
            "final_particles",
            "mean_absolute_centroid_difference_ns",
        ),
    }

    # Ratios near one mean that the generated sample differs from the reference
    # by roughly the same amount as historical resampling noise.
    for pair in generated:
        pair["ratios_to_historical_resampling_baseline"] = {
            "median_profile_wasserstein": _safe_ratio(
                pair["profile"]["median_wasserstein_ns"],
                baseline["median_profile_wasserstein_ns"],
            ),
            "median_profile_jensen_shannon": _safe_ratio(
                pair["profile"]["median_jensen_shannon"],
                baseline["median_profile_jensen_shannon"],
            ),
            "maximum_absolute_centroid_difference": _safe_ratio(
                pair["final_particles"]["maximum_absolute_centroid_difference_ns"],
                baseline["median_maximum_absolute_centroid_difference_ns"],
            ),
            "mean_absolute_centroid_difference": _safe_ratio(
                pair["final_particles"]["mean_absolute_centroid_difference_ns"],
                baseline["median_mean_absolute_centroid_difference_ns"],
            ),
        }

    # Evaluate population integrity, profile shape, final bunch widths and
    # centroids independently so a single aggregate score cannot hide failures.
    all_conserved = all(
        pair["conservation"]["reference"] and pair["conservation"]["candidate"]
        for pair in [*generated, *controls]
    )
    all_alive_and_finite = all(
        population["n_total"]
        == population["n_alive"]
        == population["n_alive_finite"]
        for pair in [*generated, *controls]
        for population in (
            pair["final_particles"]["reference_population"],
            pair["final_particles"]["candidate_population"],
        )
    )
    profile_within_2x = all(
        pair["ratios_to_historical_resampling_baseline"][
            "median_profile_wasserstein"
        ]
        is not None
        and pair["ratios_to_historical_resampling_baseline"][
            "median_profile_wasserstein"
        ]
        <= 2.0
        and pair["ratios_to_historical_resampling_baseline"][
            "median_profile_jensen_shannon"
        ]
        is not None
        and pair["ratios_to_historical_resampling_baseline"][
            "median_profile_jensen_shannon"
        ]
        <= 2.0
        for pair in generated
    )
    widths_within_5_percent = all(
        pair["final_particles"]["maximum_absolute_rms_time_relative_difference"]
        <= 0.05
        and pair["final_particles"][
            "maximum_absolute_rms_energy_relative_difference"
        ]
        <= 0.05
        for pair in generated
    )
    centroids_within_control_envelope = all(
        pair["ratios_to_historical_resampling_baseline"][
            "mean_absolute_centroid_difference"
        ]
        is not None
        and pair["ratios_to_historical_resampling_baseline"][
            "mean_absolute_centroid_difference"
        ]
        <= 2.0
        for pair in generated
    )
    # Only the >80 ns tail has enough reduced-sample particles for pass/fail.
    # The more extreme >100 and >120 ns tails remain in the report as diagnostics.
    tail80_within_20_percent = True
    for pair in generated:
        reference_tail = pair["final_particles"]["reference_tails"]["80"][
            "fraction"
        ]
        candidate_tail = pair["final_particles"]["candidate_tails"]["80"][
            "fraction"
        ]
        if reference_tail == 0:
            tail80_within_20_percent &= candidate_tail == 0
        else:
            tail80_within_20_percent &= (
                abs(candidate_tail / reference_tail - 1.0) <= 0.20
            )
    criteria = {
        "population_conserved_in_all_runs": all_conserved,
        "all_particles_survive_with_finite_final_coordinates": all_alive_and_finite,
        "profile_distances_within_2x_historical_resampling": profile_within_2x,
        "final_rms_time_and_energy_within_5_percent": widths_within_5_percent,
        "mean_absolute_final_centroid_difference_within_2x_resampling": (
            centroids_within_control_envelope
        ),
        "final_abs_time_gt_80ns_fraction_within_20_percent": tail80_within_20_percent,
    }
    return {
        "reference_path": str(reference_path),
        "sample_design": (
            "Equal tracked populations with 168-way historical stratification; "
            "generated samples have equal population in every microbunch."
        ),
        "historical_resampling_baseline": baseline,
        "generated_pairs": generated,
        "historical_control_pairs": controls,
        "assessment": {
            "name": f"{reference_particles:,}-particle RF dynamic screen",
            "passes_all_criteria": all(criteria.values()),
            "criteria": criteria,
            "tail_interpretation": (
                "The >80 ns region has enough sampled particles for a coarse screen. "
                "The >100 ns and >120 ns regions remain count-limited and are reported "
                "but excluded from pass/fail."
            ),
            "scope": (
                "Passing means the compact input is within the declared tolerances "
                "relative to historical subsampling at this reduced RF-only screen. "
                "It is not proof of exact input provenance or publication-level "
                "ghost-bunch agreement."
            ),
        },
    }
