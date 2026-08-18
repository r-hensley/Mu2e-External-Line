"""Plots that summarize quantitative RF-result comparisons."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .plotting import _save, plt
from .result_comparison import compare_results


def plot_result_comparison(
    reference_path: str | Path,
    candidate_path: str | Path,
    output_path: str | Path,
    *,
    require_matching_run: bool = True,
) -> Path:
    """Plot profile-distance histories and final captured-bunch diagnostics.

    Wasserstein-1 reports a physical displacement in nanoseconds, while the
    dimensionless Jensen-Shannon distance compares normalized profile shape.
    The lower panels compare h=28 bucket widths and low-density arrival-time
    tails for the eight original bunch labels.
    """
    report = compare_results(
        reference_path,
        candidate_path,
        require_matching_run=require_matching_run,
    )
    # These metrics compare normalized line profiles, so they measure spatial
    # evolution independently of the tracked sample's absolute population.
    profiles = report["profiles"]
    turns = np.asarray(profiles["common_turns"])
    wasserstein = np.asarray(profiles["wasserstein_ns"], dtype=float)
    js_distance = np.asarray(profiles["jensen_shannon_distance"], dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(13.2, 7.8), constrained_layout=True)
    axes[0, 0].plot(turns, wasserstein, lw=1.2)
    # The profile histogram uses 2 ns bins. This line is a scale reference,
    # not a claim that sub-bin Wasserstein distances cannot be resolved.
    axes[0, 0].axhline(
        2.0,
        color="tab:red",
        ls="--",
        lw=1,
        label="one 2 ns bin (not a sampling floor)",
    )
    axes[0, 0].set(
        xlabel="Turn",
        ylabel="Wasserstein-1 [ns]",
        title="Normalized line-profile displacement",
    )
    axes[0, 0].legend()

    axes[0, 1].plot(turns, js_distance, lw=1.2, color="tab:orange")
    axes[0, 1].set(
        xlabel="Turn",
        ylabel="Jensen-Shannon distance",
        title="Normalized line-profile shape",
        ylim=(0, max(0.05, float(np.nanmax(js_distance)) * 1.08)),
    )

    final = report["final_particles"]
    if final["reference"].get("available") and final["candidate"].get("available"):
        reference_bunches = final["reference"]["per_bunch"]
        candidate_bunches = final["candidate"]["per_bunch"]
        labels = sorted(set(reference_bunches) & set(candidate_bunches), key=int)
        x = np.arange(len(labels))
        ref_rms = [reference_bunches[label]["rms_relative_time_ns"] for label in labels]
        cand_rms = [candidate_bunches[label]["rms_relative_time_ns"] for label in labels]
        width = 0.38
        axes[1, 0].bar(x - width / 2, ref_rms, width, label="Historical")
        axes[1, 0].bar(x + width / 2, cand_rms, width, label="Generated")
        axes[1, 0].set(
            xlabel="Initial bunch label",
            ylabel="Final RMS time [ns]",
            title="Final h=28 bucket widths",
            xticks=x,
            # Simulation labels are zero-based; plots use familiar bunch 1--8.
            xticklabels=[str(int(label) + 1) for label in labels],
        )
        axes[1, 0].legend()

        # Average each bucket-relative tail fraction over the eight initial
        # labels to summarize satellite populations without hiding asymmetry
        # in the per-bunch values retained in the JSON report.
        thresholds = (80, 100, 120)
        ref_tail = [
            np.mean(
                [
                    reference_bunches[label][f"fraction_abs_time_gt_{threshold}_ns"]
                    for label in labels
                ]
            )
            for threshold in thresholds
        ]
        cand_tail = [
            np.mean(
                [
                    candidate_bunches[label][f"fraction_abs_time_gt_{threshold}_ns"]
                    for label in labels
                ]
            )
            for threshold in thresholds
        ]
        x_tail = np.arange(len(thresholds))
        axes[1, 1].bar(x_tail - width / 2, ref_tail, width, label="Historical")
        axes[1, 1].bar(x_tail + width / 2, cand_tail, width, label="Generated")
        axes[1, 1].set(
            xlabel="Bucket-relative tail threshold [ns]",
            ylabel="Mean particle fraction",
            title="Final low-density tails (all 8 labels)",
            xticks=x_tail,
            xticklabels=[str(value) for value in thresholds],
            yscale="log",
        )
        axes[1, 1].legend()
    else:
        for axis in axes[1]:
            axis.text(0.5, 0.5, "Final-particle data unavailable", ha="center", va="center")
            axis.set_axis_off()

    for axis in axes.ravel():
        axis.grid(True, color="0.9", lw=0.6)
    fig.suptitle("Generated-input RF evolution versus historical-input reference")
    return _save(fig, output_path)
