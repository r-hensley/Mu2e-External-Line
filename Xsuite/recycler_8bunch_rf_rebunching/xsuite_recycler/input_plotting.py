"""Diagnostic plots for the historical and generated longitudinal inputs."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .input_analysis import OCCUPIED_GRID_INDICES, characterize_input
from .plotting import _save, plt


def _within_microbunch_residuals(values: np.ndarray) -> np.ndarray:
    """Subtract each ordered microbunch block's own mean from its samples."""
    # Both supported input sources store 168 equal, contiguous microbunch
    # populations, so reshaping preserves the 168 microbunch blocks.
    matrix = np.asarray(values, dtype=float).reshape(168, -1)
    return matrix - np.mean(matrix, axis=1, keepdims=True)


def plot_input_comparison(
    reference_dt_s,
    reference_dpop,
    candidate_dt_s,
    candidate_dpop,
    output_path: str | Path,
) -> Path:
    """Compare historical and generated microbunch geometry and shape.

    Arrival-time arrays are in seconds and momentum arrays are fractional
    ``dp/p``. Each must contain 168 equally populated, ordered microbunches.
    The top panels compare fitted centers and widths; the bottom panels remove
    each microbunch mean and compare the remaining standardized distributions.
    """
    reference = characterize_input(reference_dt_s, reference_dpop)
    candidate = characterize_input(candidate_dt_s, candidate_dpop)
    microbunch = np.arange(168)
    # The historical fit supplies the common RF-grid origin and spacing;
    # OCCUPIED_GRID_INDICES identifies which 53 MHz slots actually contain beam.
    reference_grid_ns = (
        reference.grid_origin_ns
        + reference.grid_spacing_ns * OCCUPIED_GRID_INDICES
    )

    # Normalize both data sets with historical pooled widths. Using one common
    # scale prevents an overly broad candidate from normalizing away its error.
    reference_time_residual = (
        _within_microbunch_residuals(np.asarray(reference_dt_s) * 1e9)
        / reference.pooled_time_sigma_ns
    ).ravel()
    candidate_time_residual = (
        _within_microbunch_residuals(np.asarray(candidate_dt_s) * 1e9)
        / reference.pooled_time_sigma_ns
    ).ravel()
    reference_dpop_residual = (
        _within_microbunch_residuals(np.asarray(reference_dpop))
        / reference.pooled_dpop_sigma
    ).ravel()
    candidate_dpop_residual = (
        _within_microbunch_residuals(np.asarray(candidate_dpop))
        / reference.pooled_dpop_sigma
    ).ravel()

    fig, axes = plt.subplots(2, 2, figsize=(11.4, 7.8), constrained_layout=True)
    axes[0, 0].plot(
        microbunch,
        reference.time_center_ns - reference_grid_ns,
        ".",
        ms=3,
        label="Historical",
    )
    axes[0, 0].plot(
        microbunch,
        candidate.time_center_ns - reference_grid_ns,
        ".",
        ms=3,
        label="Generated",
    )
    axes[0, 0].axhline(0, color="0.5", lw=0.7)
    axes[0, 0].set(
        xlabel="Ordered microbunch index",
        ylabel="Center residual from historical grid [ns]",
        title="Microbunch center geometry",
    )
    axes[0, 0].legend()

    axes[0, 1].plot(
        microbunch,
        100 * (candidate.time_sigma_ns / reference.time_sigma_ns - 1),
        ".",
        ms=3,
        label=r"$\sigma_t$",
    )
    axes[0, 1].plot(
        microbunch,
        100 * (candidate.dpop_sigma / reference.dpop_sigma - 1),
        ".",
        ms=3,
        label=r"$\sigma_{dp/p}$",
    )
    axes[0, 1].axhline(0, color="0.5", lw=0.7)
    axes[0, 1].set(
        xlabel="Ordered microbunch index",
        ylabel="Generated - historical [%]",
        title="Per-microbunch width differences",
    )
    axes[0, 1].legend()

    bins = np.linspace(-4.5, 4.5, 121)
    normal_x = np.linspace(-4.5, 4.5, 500)
    normal_pdf = np.exp(-0.5 * normal_x**2) / np.sqrt(2 * np.pi)
    for axis, historical, generated, label in (
        (
            axes[1, 0],
            reference_time_residual,
            candidate_time_residual,
            r"$(t-\bar{t}_{\rm micro})/\sigma_{t,\rm hist}$",
        ),
        (
            axes[1, 1],
            reference_dpop_residual,
            candidate_dpop_residual,
            r"$(dp/p-\overline{dp/p}_{\rm micro})/\sigma_{dp/p,\rm hist}$",
        ),
    ):
        axis.hist(
            historical,
            bins=bins,
            density=True,
            histtype="step",
            lw=1.5,
            label="Historical",
        )
        axis.hist(
            generated,
            bins=bins,
            density=True,
            histtype="step",
            lw=1.5,
            label="Generated",
        )
        axis.plot(normal_x, normal_pdf, color="0.25", ls="--", lw=1, label="Normal")
        axis.set(xlabel=label, ylabel="Density", yscale="log", ylim=(1e-5, 1))
        axis.grid(True, color="0.9", lw=0.6)
        axis.legend()
    axes[1, 0].set_title("Within-microbunch time shape")
    axes[1, 1].set_title("Within-microbunch momentum shape")
    fig.suptitle("Historical input versus compact generated model")
    return _save(fig, output_path)
