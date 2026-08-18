"""Standard plots for the RF-only Recycler rebunching simulation."""

from __future__ import annotations

import os
from pathlib import Path
from tempfile import gettempdir

import h5py

# This WSL host's default Matplotlib config directory is read-only. Set a
# package-specific temporary cache before importing Matplotlib; callers can
# still override it explicitly with MPLCONFIGDIR.
_mpl_config_dir = Path(gettempdir()) / "xsuite-recycler-matplotlib"
_mpl_config_dir.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_mpl_config_dir))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from .config import RFProgramConfig, RecyclerConfig
from .distribution import zeta_to_dt
from .rf_program import RFProgram


def _save(fig: plt.Figure, path: str | Path) -> Path:
    """Save and close a noninteractive Matplotlib figure."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_rf_voltage_program(
    output_path: str | Path,
    ramp: RFProgramConfig | None = None,
) -> Path:
    """Plot the h=588 turn-off and h=28 capture ramps in laboratory units.

    The harmonic number is the number of RF cycles per Recycler revolution;
    h=588 is approximately 53 MHz and h=28 is approximately 2.5 MHz.
    """
    ramp = ramp or RFProgramConfig()
    program = RFProgram(ramp)
    time_s = np.linspace(0.0, 0.100, 2001)

    fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    ax.plot(
        time_s * 1e3,
        program.voltage_25(time_s) / 1e3,
        color="red",
        lw=2,
        label="2.5 MHz (h=28)",
    )
    ax.plot(
        time_s * 1e3,
        program.voltage_53(time_s) / 1e3,
        color="royalblue",
        lw=2,
        label="53 MHz (h=588)",
    )
    ax.axvline(ramp.rf53_turnoff_s * 1e3, color="0.35", ls="--", lw=1)
    ax.axvline(ramp.end_s * 1e3, color="0.35", ls="--", lw=1)
    ax.set(
        xlim=(0, 100),
        ylim=(0, 100),
        xlabel="Time [ms]",
        ylabel="Voltage [kV]",
        title="Recycler RF voltage program",
    )
    ax.grid(True, color="0.85")
    ax.legend(loc="upper right")
    return _save(fig, output_path)


def _machine_and_ramp_from_h5(h5: h5py.File) -> tuple[RecyclerConfig, RFProgramConfig]:
    """Reconstruct the machine and RF dataclasses stored in a result file."""
    machine_keys = RecyclerConfig.__dataclass_fields__.keys()
    ramp_keys = RFProgramConfig.__dataclass_fields__.keys()
    machine = RecyclerConfig(**{key: h5["machine"].attrs[key] for key in machine_keys})
    ramp = RFProgramConfig(**{key: h5["rf_program"].attrs[key] for key in ramp_keys})
    return machine, ramp


def plot_final_phase_space(
    results_path: str | Path,
    output_path: str | Path,
    *,
    bunch_index: int = 7,
    max_scatter: int = 100_000,
    seed: int = 202208,
) -> Path:
    """Plot final longitudinal phase space for one initial bunch label.

    ``bunch_index`` is zero-based. Surviving Xsuite particles are folded into
    the nearest h=28 bucket, then plotted as bucket-relative arrival time
    versus energy offset. A large ensemble is sampled only for rendering; the
    saved simulation data are not modified.
    """
    with h5py.File(results_path, "r") as h5:
        machine, _ = _machine_and_ramp_from_h5(h5)
        zeta = h5["final_particles/zeta_m"][:]
        ptau = h5["final_particles/ptau"][:]
        state = h5["final_particles/state"][:]
        labels = h5["final_particles/bunch_index"][:]
        n_turns = int(h5.attrs["n_turns"])

    # Xsuite marks active particles with a positive state; zero or negative
    # states represent particles that have been stopped or lost.
    selected = (labels == bunch_index) & (state > 0)
    if not np.any(selected):
        raise ValueError(f"No surviving particles found for bunch index {bunch_index}")
    # Convert Xsuite's longitudinal distance zeta back to the arrival-time
    # convention used by the historical BLonD input and the paper plots.
    dt_s = zeta_to_dt(zeta[selected], machine.beta0)
    bucket_period_s = 1.0 / machine.rf25_frequency_hz
    bucket_center_s = np.rint(np.median(dt_s) / bucket_period_s) * bucket_period_s
    # Modulo wrapping places every particle in [-T/2, T/2) about the nearest
    # h=28 bucket center, making the captured phase-space shape easy to read.
    dt_relative_s = (
        (dt_s - bucket_center_s + 0.5 * bucket_period_s) % bucket_period_s
        - 0.5 * bucket_period_s
    )
    # Xsuite defines ptau=(E-E0)/p0c, so multiplying by p0c gives eV.
    energy_offset_ev = ptau[selected] * machine.p0c_ev

    if dt_relative_s.size > max_scatter:
        # Downsampling controls image size only and is deterministic by seed.
        rng = np.random.default_rng(seed)
        take = rng.choice(dt_relative_s.size, size=max_scatter, replace=False)
        dt_relative_s = dt_relative_s[take]
        energy_offset_ev = energy_offset_ev[take]

    fig, ax = plt.subplots(figsize=(7.0, 5.2), constrained_layout=True)
    ax.scatter(
        dt_relative_s,
        energy_offset_ev,
        s=0.35,
        alpha=0.42,
        color="#2878b5",
        edgecolors="none",
        rasterized=True,
    )
    ax.set(
        xlim=(-80e-9, 80e-9),
        ylim=(-4e7, 4e7),
        xlabel="Time relative to h=28 bucket center [s]",
        ylabel="Energy offset [eV]",
        title=(
            "Final RF-only longitudinal phase space "
            f"(bunch {bunch_index + 1}, turn {n_turns:,})"
        ),
    )
    ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    ax.grid(True, color="0.9", lw=0.6)
    return _save(fig, output_path)


def _centers_to_edges(centers: np.ndarray, lower: float, upper: float) -> np.ndarray:
    """Convert sampled coordinate centers to edges suitable for pcolormesh."""
    if centers.size == 1:
        return np.asarray([lower, upper], dtype=float)
    edges = np.empty(centers.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = lower
    edges[-1] = upper
    return edges


def plot_rebunching_waterfall(
    results_path: str | Path,
    output_path: str | Path,
) -> Path:
    """Plot longitudinal density throughout the RF-only tracking cycle.

    The HDF5 profiles are histograms in arrival time recorded at selected
    turns. Their values are already scaled from the tracked sample to the
    equivalent historical macroparticle population. Log and linear panels
    expose both low-density satellite structure and the dominant bunches.
    """
    with h5py.File(results_path, "r") as h5:
        turns = h5["profiles/turn"][:]
        time_edges_ns = h5["profiles/time_edges_ns"][:]
        counts = h5["profiles/counts_macro_equivalent"][:].astype(float)
        n_turns = int(h5.attrs["n_turns"])
        sample_particles = int(h5.attrs["sample_particles"])
        source_rows = int(h5.attrs["source_rows"])

    positive = counts[counts > 0]
    if positive.size == 0:
        raise ValueError("Waterfall contains no nonzero bins")
    # A robust upper color limit prevents a few dense bins from hiding the
    # weaker structures that motivate the waterfall diagnostic.
    vmax = float(np.percentile(positive, 99.8))
    vmax = max(vmax, float(positive.max()) * 0.2)
    vmin = max(float(positive.min()), vmax * 1e-4)
    # Profiles are saved at sampled turn centers; pcolormesh needs bin edges.
    turn_edges = _centers_to_edges(turns.astype(float), 0.0, float(n_turns))

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(13.2, 5.4),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    common = dict(shading="auto", cmap="jet", rasterized=True)
    log_mesh = axes[0].pcolormesh(
        time_edges_ns,
        turn_edges,
        np.ma.masked_less_equal(counts, 0),
        norm=LogNorm(vmin=vmin, vmax=vmax),
        **common,
    )
    linear_mesh = axes[1].pcolormesh(
        time_edges_ns,
        turn_edges,
        counts,
        vmin=0.0,
        vmax=vmax,
        **common,
    )
    for ax, title in zip(axes, ("Log density", "Linear density")):
        ax.set(
            xlabel=r"$\Delta t$ [ns]",
            title=title,
            xlim=(time_edges_ns[0], time_edges_ns[-1]),
            ylim=(0, n_turns),
        )
    axes[0].set_ylabel("Turns through Recycler")
    fig.colorbar(log_mesh, ax=axes[0], label="Equivalent macroparticles / 2 ns")
    fig.colorbar(linear_mesh, ax=axes[1], label="Equivalent macroparticles / 2 ns")
    fig.suptitle(
        f"RF-only rebunching: {sample_particles:,} tracked, scaled to {source_rows:,} macroparticles",
        fontsize=12,
    )
    return _save(fig, output_path)


def plot_run_results(results_path: str | Path, output_dir: str | Path) -> list[Path]:
    """Write the standard RF-program, final-phase-space, and waterfall plots.

    The RF configuration is recovered from the HDF5 result so replotting uses
    the same voltage program that produced the tracked particles.
    """
    output_dir = Path(output_dir)
    with h5py.File(results_path, "r") as h5:
        _, ramp = _machine_and_ramp_from_h5(h5)
    return [
        plot_rf_voltage_program(output_dir / "rf_voltage_program.png", ramp),
        plot_final_phase_space(results_path, output_dir / "final_phase_space.png"),
        plot_rebunching_waterfall(
            results_path,
            output_dir / "rebunching_waterfall.png",
        ),
    ]
