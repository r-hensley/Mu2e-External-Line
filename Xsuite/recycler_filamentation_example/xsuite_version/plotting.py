"""Plots for the Xsuite Recycler filamentation result."""

from __future__ import annotations

import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c as c_light_si
from scipy.constants import e as elementary_charge_c
from xpart.longitudinal.rf_bucket import RFBucket

try:
    from .simulation import (
        RecyclerInputs,
        XsuiteDerived,
        XsuiteSimulationResult,
        bunch_diagnostics,
        rf_stage,
        rf_voltages,
    )
except ImportError:  # Support direct execution through sibling ``run.py``.
    from simulation import (
        RecyclerInputs,
        XsuiteDerived,
        XsuiteSimulationResult,
        bunch_diagnostics,
        rf_stage,
        rf_voltages,
    )


def source_colors(source_id: np.ndarray) -> np.ndarray:
    """Return passive colors tied to each original 53 MHz bunchlet."""

    values = np.asarray(source_id, dtype=float)
    normalized = (values - values.min()) / (values.max() - values.min())
    return mpl.colormaps["turbo"](normalized)


def h28_separatrix(
    inputs: RecyclerInputs,
    derived: XsuiteDerived,
    voltage_eV: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return an h=28 separatrix from Xpart's ``RFBucket``."""

    if voltage_eV <= 0.0:
        empty = np.asarray([], dtype=float)
        return empty, empty, empty

    # Use Xpart's RF Hamiltonian implementation for the plotted bucket rather
    # than reusing the analytical separatrix from the NumPy reference model.
    bucket = RFBucket(
        circumference=inputs.circumference_m,
        gamma=derived.gamma0,
        mass_kg=derived.mass0_eV * elementary_charge_c / c_light_si**2,
        charge_coulomb=elementary_charge_c,
        alpha_array=np.asarray([derived.momentum_compaction_factor]),
        p_increment=0.0,
        harmonic_list=np.asarray([inputs.harmonic_25]),
        voltage_list=np.asarray([float(voltage_eV)]),
        phi_offset_list=np.asarray([0.0]),
        z_offset=0.0,
    )
    zeta = np.linspace(bucket.z_left, bucket.z_right, 700)
    delta = bucket.separatrix(zeta)
    p_plus = derived.p0c_eV * (1.0 + delta)
    p_minus = derived.p0c_eV * (1.0 - delta)
    energy_plus = np.sqrt(p_plus**2 + derived.mass0_eV**2) - derived.energy0_eV
    energy_minus = np.sqrt(p_minus**2 + derived.mass0_eV**2) - derived.energy0_eV

    # RFBucket and Xsuite use zeta, whose sign is opposite to the positive
    # arrival-time ring phase used on the teaching plots.
    theta = -2.0 * np.pi * zeta / inputs.circumference_m
    order = np.argsort(theta)
    return theta[order], energy_plus[order], energy_minus[order]


def _phase_space_axes(ax: plt.Axes, inputs: RecyclerInputs) -> None:
    half_cell_deg = math.degrees(math.pi / inputs.harmonic_25)
    ax.set_xlim(-1.02 * half_cell_deg, 1.02 * half_cell_deg)
    ax.set_ylim(-45.0, 45.0)
    ax.set_xlabel(r"Recycler longitudinal phase $\theta$ [deg]")
    ax.set_ylabel(r"$\Delta E$ [MeV]")
    ax.grid(alpha=0.18)


def _scatter(
    ax: plt.Axes,
    result: XsuiteSimulationResult,
    index: int,
    *,
    size: float = 3.0,
):
    display = result.display_indices
    return ax.scatter(
        np.degrees(result.theta_rad[index, display]),
        result.delta_energy_eV[index, display] / 1e6,
        c=result.source_id[display],
        s=size,
        alpha=0.72,
        cmap="turbo",
        vmin=-10.5,
        vmax=10.5,
        linewidths=0.0,
        rasterized=True,
    )


def _plot_separatrix(
    ax: plt.Axes,
    inputs: RecyclerInputs,
    derived: XsuiteDerived,
    voltage_eV: float,
) -> None:
    theta, upper, lower = h28_separatrix(inputs, derived, voltage_eV)
    if theta.size:
        ax.plot(np.degrees(theta), upper / 1e6, color="#c92f2f", lw=1.6)
        ax.plot(np.degrees(theta), lower / 1e6, color="#c92f2f", lw=1.6)


def plot_rf_program(inputs: RecyclerInputs, derived: XsuiteDerived) -> plt.Figure:
    times = np.linspace(0.0, 105.0e-3, 1200)
    voltage_53, voltage_25 = rf_voltages(times, inputs, derived)
    fig, ax = plt.subplots(figsize=(10.5, 4.2), constrained_layout=True)
    ax.plot(times * 1e3, voltage_53 / 1e3, lw=2.4, label="53 MHz, h=588")
    ax.plot(times * 1e3, voltage_25 / 1e3, lw=2.4, label="2.5 MHz, h=28")
    ax.axvline(derived.time_53_off_s * 1e3, color="0.35", lw=1)
    ax.axvline(derived.time_25_flat_s * 1e3, color="0.35", lw=1)
    ax.set(
        xlim=(0.0, 105.0),
        ylim=(0.0, 90.0),
        xlabel="Time [ms]",
        ylabel="RF voltage [kV]",
        title="RF program used by the Xsuite cavities",
    )
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    return fig


def plot_checkpoints(result: XsuiteSimulationResult) -> plt.Figure:
    requested_ms = [0.0, result.derived.time_53_off_s * 1e3, 25.0, 55.0,
                    result.derived.time_25_flat_s * 1e3]
    fig, axes = plt.subplots(
        2, 3, figsize=(13.0, 8.2), sharex=True, sharey=True,
        constrained_layout=True,
    )
    scatter = None
    for ax, time_ms in zip(axes.ravel(), requested_ms):
        index = result.nearest_snapshot_index(time_ms * 1e-3)
        time_s = float(result.snapshot_times_s[index])
        scatter = _scatter(ax, result, index)
        _, voltage_25 = rf_voltages(time_s, result.inputs, result.derived)
        _plot_separatrix(ax, result.inputs, result.derived, float(voltage_25))
        _phase_space_axes(ax, result.inputs)
        ax.set_title(f"{time_s*1e3:.2f} ms\n{rf_stage(time_s, result.inputs, result.derived)}")
    axes.ravel()[-1].axis("off")
    if scatter is not None:
        cbar = fig.colorbar(scatter, ax=list(axes.ravel()[:-1]), shrink=0.86, pad=0.012)
        cbar.set_label("Original 53 MHz bunch index")
        cbar.set_ticks([-10, -5, 0, 5, 10])
    fig.suptitle("Xsuite: from 53 MHz bunchlets to a filamented 2.5 MHz bunch")
    return fig


def plot_extraction_gallery(result: XsuiteSimulationResult) -> plt.Figure:
    fig, axes = plt.subplots(
        2, 4, figsize=(14.0, 7.0), sharex=True, sharey=True,
        constrained_layout=True,
    )
    scatter = None
    for number, (ax, time_s) in enumerate(
        zip(axes.ravel(), result.derived.extraction_times_s), start=1
    ):
        index = result.nearest_snapshot_index(float(time_s))
        scatter = _scatter(ax, result, index, size=2.2)
        _plot_separatrix(
            ax, result.inputs, result.derived, result.inputs.voltage_25_final_eV
        )
        _phase_space_axes(ax, result.inputs)
        diag = bunch_diagnostics(
            result.theta_rad[index],
            result.delta_energy_eV[index],
            result.inputs,
            result.derived,
        )
        ax.set_title(
            f"Sample {number}: {result.snapshot_times_s[index]*1e3:.1f} ms\n"
            f"σt={diag['sigma_time_s']*1e9:.1f} ns",
            fontsize=9,
        )
    if scatter is not None:
        cbar = fig.colorbar(scatter, ax=list(axes.ravel()), shrink=0.86, pad=0.012)
        cbar.set_label("Original 53 MHz bunch index")
        cbar.set_ticks([-10, -5, 0, 5, 10])
    fig.suptitle("Xsuite states at the eight extraction times")
    return fig


def plot_waterfall(result: XsuiteSimulationResult) -> plt.Figure:
    """Plot normalized projections at the stored 134-turn cadence."""

    mask = result.snapshot_times_s <= result.derived.time_25_flat_s + 0.51 * 134 * result.derived.t_rev0_s
    indices = np.flatnonzero(mask)
    bins = np.linspace(-math.pi / result.inputs.harmonic_25,
                       math.pi / result.inputs.harmonic_25, 121)
    centers = 0.5 * (bins[:-1] + bins[1:])
    fig, ax = plt.subplots(figsize=(11.0, 8.0), constrained_layout=True)
    for index in indices:
        counts, _ = np.histogram(result.theta_rad[index], bins=bins)
        profile = counts / max(counts.max(), 1)
        offset = result.snapshot_times_s[index] * 1e3
        ax.plot(np.degrees(centers), offset + 1.1 * profile, color="#3a9250", lw=0.72)
    ax.axvline(-4.0, color="red", ls="--", lw=1.0)
    ax.axvline(4.0, color="red", ls="--", lw=1.0)
    ax.set(
        xlabel=r"Recycler longitudinal phase $\theta$ [deg]",
        ylabel="Time in RF program [ms]",
        title="Xsuite waterfall of the h=28-cell time projection",
    )
    ax.grid(alpha=0.18)
    return fig


def plot_diagnostic_comparison(
    result: XsuiteSimulationResult,
    reference: np.lib.npyio.NpzFile,
) -> plt.Figure:
    """Overlay bulk diagnostics from Xsuite and the preserved NumPy map."""

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.5), constrained_layout=True)
    x_ms = result.diagnostic_times_s * 1e3
    ref_x_ms = reference["diagnostic_times_s"] * 1e3
    panels = [
        ("sigma_time_s", 1e9, r"$\sigma_t$ [ns]"),
        ("sigma_energy_eV", 1e-6, r"$\sigma_E$ [MeV]"),
        ("covariance_emittance_eVs", 1.0, r"$\epsilon_{cov}$ [eV s]"),
        ("central_95_time_width_s", 1e9, "Central 95% time width [ns]"),
    ]
    for ax, (field, scale, ylabel) in zip(axes.ravel(), panels):
        ax.plot(ref_x_ms, reference[field] * scale, lw=1.5, label="NumPy reference")
        ax.plot(x_ms, getattr(result, field) * scale, lw=1.2, ls="--", label="Xsuite")
        ax.set(xlabel="Time [ms]", ylabel=ylabel)
        ax.grid(alpha=0.22)
    axes[0, 0].legend(frameon=False)
    fig.suptitle("Longitudinal diagnostics: Xsuite versus preserved reference")
    return fig


def plot_phase_space_comparison(
    result: XsuiteSimulationResult,
    reference: np.lib.npyio.NpzFile,
) -> plt.Figure:
    """Place the two 90 ms phase-space endpoints side by side."""

    manual_turns = reference["snapshot_turns"]
    flat_turn = int(round(result.derived.time_25_flat_s / result.derived.t_rev0_s))
    manual_index = int(np.argmin(np.abs(manual_turns - flat_turn)))
    xsuite_index = int(np.argmin(np.abs(result.snapshot_turns - flat_turn)))
    display = result.display_indices
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.3), sharex=True, sharey=True,
                             constrained_layout=True)
    for ax, title, theta, energy in [
        (axes[0], "Preserved NumPy map", reference["theta_rad"][manual_index],
         reference["delta_energy_eV"][manual_index]),
        (axes[1], "Xsuite line", result.theta_rad[xsuite_index],
         result.delta_energy_eV[xsuite_index]),
    ]:
        ax.scatter(
            np.degrees(theta[display]), energy[display] / 1e6,
            c=result.source_id[display], cmap="turbo", vmin=-10.5, vmax=10.5,
            s=2.5, alpha=0.7, linewidths=0, rasterized=True,
        )
        _plot_separatrix(ax, result.inputs, result.derived,
                         result.inputs.voltage_25_final_eV)
        _phase_space_axes(ax, result.inputs)
        ax.set_title(f"{title}\nturn {result.snapshot_turns[xsuite_index]:,}")
    fig.suptitle("90 ms rebunching endpoint comparison")
    return fig
