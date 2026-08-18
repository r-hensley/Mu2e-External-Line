"""RF-only tracking and compact diagnostics for the Recycler reconstruction."""

from __future__ import annotations

from dataclasses import asdict
import json
from math import ceil
from pathlib import Path
import shlex
from time import monotonic

import h5py
import numpy as np
import xtrack as xt

from .config import PROTON_MASS_EV, RFProgramConfig, RecyclerConfig
from .distribution import InitialDistribution, load_initial_distribution, zeta_to_dt
from .provenance import collect_run_provenance
from .rf_program import RFProgram


def build_line(machine: RecyclerConfig) -> xt.Line:
    """Build the RF-kick-then-slip Xsuite model for one Recycler turn.

    The two :class:`xtrack.Cavity` elements provide the h=588 and h=28 energy
    kicks.  A nonlinear :class:`xtrack.LineSegmentMap` then advances the
    longitudinal position using the Recycler momentum compaction.  The model
    has no transverse dynamics, impedance, space charge, or other collective
    effects.
    """
    line = xt.Line(
        elements=[
            # These voltages start at zero because run_rf_only overwrites them
            # with the time-programmed values before every one-turn track.
            xt.Cavity(harmonic=machine.harmonic_53, voltage=0.0, lag=0.0),
            xt.Cavity(harmonic=machine.harmonic_25, voltage=0.0, lag=0.0),
            xt.LineSegmentMap(
                length=machine.circumference_m,
                longitudinal_mode="nonlinear",
                momentum_compaction_factor=machine.momentum_compaction_factor,
                slippage_length=machine.circumference_m,
                # RF is represented by the explicit cavities above; keeping
                # the map's built-in RF at zero avoids applying a second kick.
                voltage_rf=[0.0],
                frequency_rf=[0.0],
                lag_rf=[0.0],
            ),
        ],
        element_names=["rf53", "rf2p5", "slip"],
    )
    # The reference particle supplies Xsuite with the mass, charge, momentum,
    # and derived beta/gamma used by both the cavities and longitudinal map.
    line.particle_ref = xt.Particles(
        mass0=PROTON_MASS_EV,
        q0=1,
        p0c=machine.p0c_ev,
    )
    line.build_tracker()
    return line


def make_particles(
    distribution: InitialDistribution,
    machine: RecyclerConfig,
) -> xt.Particles:
    """Create Xsuite particles from the prepared longitudinal coordinates.

    Unspecified transverse coordinates remain at their Xsuite defaults.  A
    stable integer ID is assigned to each input row for diagnostics.
    """
    return xt.Particles(
        mass0=PROTON_MASS_EV,
        q0=1,
        p0c=machine.p0c_ev,
        zeta=distribution.zeta_m,
        ptau=distribution.ptau,
        particle_id=np.arange(distribution.n_particles, dtype=np.int64),
    )


def _record_turns(n_turns: int, record_every: int) -> np.ndarray:
    """Return profile turn numbers, always including the initial and final state."""
    if record_every <= 0:
        raise ValueError("record_every must be positive")
    values = [0, *range(record_every, n_turns + 1, record_every)]
    if values[-1] != n_turns:
        values.append(n_turns)
    return np.asarray(values, dtype=np.int64)


def _profile(
    particles: xt.Particles,
    machine: RecyclerConfig,
    time_edges_ns: np.ndarray,
    macro_weight: float,
) -> tuple[np.ndarray, float]:
    """Histogram live particles in historical arrival-time coordinates.

    Bin contents are multiplied by the macro-equivalent sample weight.  The
    separately returned outside count preserves particles beyond the plotted
    time window, allowing population conservation to be checked.
    """
    dt_ns = zeta_to_dt(np.asarray(particles.zeta), machine.beta0) * 1e9
    alive = np.asarray(particles.state) > 0
    counts, _ = np.histogram(dt_ns[alive], bins=time_edges_ns)
    counts = counts * macro_weight
    outside = macro_weight * np.count_nonzero(
        alive & ((dt_ns < time_edges_ns[0]) | (dt_ns >= time_edges_ns[-1]))
    )
    return counts.astype(np.float32), float(outside)


def _write_results(
    output_path: Path,
    machine: RecyclerConfig,
    ramp: RFProgramConfig,
    distribution: InitialDistribution,
    particles: xt.Particles,
    record_turns: np.ndarray,
    time_edges_ns: np.ndarray,
    counts: np.ndarray,
    outside_counts: np.ndarray,
    n_turns: int,
    record_every: int,
    seed: int,
    elapsed_s: float,
    run_provenance: dict[str, object],
    run_configuration: dict[str, object],
    invocation_argv: list[str] | None,
) -> None:
    """Write configuration, streamed profiles, and final particles to HDF5."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_path, "w") as h5:
        h5.attrs["format_version"] = 1
        h5.attrs["model"] = "RF-only: h=588 and h=28; no impedance or space charge"
        h5.attrs["xtrack_version"] = xt.__version__
        h5.attrs["source_blond_commit"] = "3566338e5e1e201028f2ec3ec62a94dcf1fc0031"
        h5.attrs["source_blond_driver"] = "BLonD/multi_mpi.py"
        h5.attrs["coordinate_convention"] = "zeta=-beta0*c*(input_z/3e8)"
        h5.attrs["energy_convention"] = "ptau=input_dpop/0.9944"
        h5.attrs["n_turns"] = n_turns
        h5.attrs["end_state_time_s"] = n_turns * machine.revolution_period_s
        h5.attrs["last_kick_time_s"] = (n_turns - 1) * machine.revolution_period_s
        h5.attrs["first_rf25_kick_turn"] = ceil(
            ramp.rf25_start_s / machine.revolution_period_s
        )
        h5.attrs["record_every"] = record_every
        h5.attrs["random_seed"] = seed
        h5.attrs["elapsed_tracking_s"] = elapsed_s
        h5.attrs["sample_particles"] = distribution.n_particles
        h5.attrs["source_rows"] = distribution.source_rows
        h5.attrs["macro_equivalent_weight"] = distribution.macro_equivalent_weight
        h5.attrs["physical_proton_weight"] = distribution.physical_proton_weight
        h5.attrs["input_path"] = distribution.source_path
        h5.attrs["input_sha256"] = distribution.source_sha256
        h5.attrs["input_source_coordinate_sha256"] = (
            distribution.source_coordinate_sha256
        )
        h5.attrs["input_tracking_coordinate_sha256"] = (
            distribution.tracking_coordinate_sha256
        )
        h5.attrs["input_selected_indices_sha256"] = (
            distribution.selected_indices_sha256
        )
        h5.attrs["input_sampling_method"] = distribution.sampling_method
        # -1 is a typed sentinel only for HDF5 readers; the source metadata JSON
        # preserves the distinction as a native JSON null.
        h5.attrs["input_sampling_seed"] = (
            -1 if distribution.sampling_seed is None else distribution.sampling_seed
        )
        h5.attrs["input_sampling_strata"] = distribution.sampling_strata
        h5.attrs["input_source_kind"] = distribution.source_kind
        h5.attrs["input_source_metadata_json"] = distribution.source_metadata_json

        source_tree = run_provenance["source_tree"]
        git = run_provenance["git"]
        runtime = run_provenance["runtime"]
        h5.attrs["project_source_sha256"] = source_tree["source_tree_sha256"]
        h5.attrs["project_source_files_json"] = json.dumps(
            source_tree["source_files"], separators=(",", ":")
        )
        h5.attrs["project_git_revision"] = git["revision"] or ""
        h5.attrs["project_git_dirty"] = (
            "unknown" if git["dirty"] is None else str(git["dirty"]).lower()
        )
        h5.attrs["project_git_path"] = git["project_path"] or ""
        h5.attrs["project_git_status_porcelain"] = git["status_porcelain"] or ""
        h5.attrs["python_version"] = runtime["python_version"]
        h5.attrs["python_executable"] = runtime["python_executable"]
        h5.attrs["dependency_versions_json"] = json.dumps(
            runtime["dependencies"],
            sort_keys=True,
            separators=(",", ":"),
        )
        h5.attrs["runtime_provenance_json"] = json.dumps(
            runtime,
            sort_keys=True,
            separators=(",", ":"),
        )
        h5.attrs["run_configuration_json"] = json.dumps(
            run_configuration,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        h5.attrs["invocation_argv_json"] = json.dumps(
            invocation_argv,
            separators=(",", ":"),
        )
        h5.attrs["invocation_command"] = (
            shlex.join(invocation_argv) if invocation_argv is not None else ""
        )
        last_voltage_53, last_voltage_25 = RFProgram(ramp).voltages(
            (n_turns - 1) * machine.revolution_period_s
        )
        h5.attrs["last_kick_rf53_v"] = last_voltage_53
        h5.attrs["last_kick_rf25_v"] = last_voltage_25

        # Configuration values live in attributes so analysis code can recover
        # the exact machine and RF program without importing this package.
        machine_group = h5.create_group("machine")
        for key, value in asdict(machine).items():
            machine_group.attrs[key] = value
        for key in (
            "beta0",
            "gamma0",
            "momentum_compaction_factor",
            "slip_factor",
            "revolution_period_s",
            "revolution_frequency_hz",
            "rf53_frequency_hz",
            "rf25_frequency_hz",
        ):
            machine_group.attrs[key] = getattr(machine, key)

        ramp_group = h5.create_group("rf_program")
        for key, value in asdict(ramp).items():
            ramp_group.attrs[key] = value

        # Profiles are the compact turn-by-turn diagnostic; storing them avoids
        # retaining every particle coordinate at every recorded turn.
        profiles = h5.create_group("profiles")
        profiles.create_dataset("turn", data=record_turns)
        profiles.create_dataset("time_edges_ns", data=time_edges_ns)
        profiles.create_dataset(
            "counts_macro_equivalent",
            data=counts,
            compression="gzip",
            shuffle=True,
            chunks=(1, counts.shape[1]),
        )
        profiles.create_dataset("outside_macro_equivalent", data=outside_counts)

        final = h5.create_group("final_particles")
        final.create_dataset("zeta_m", data=np.asarray(particles.zeta), compression="gzip")
        final.create_dataset("ptau", data=np.asarray(particles.ptau), compression="gzip")
        final.create_dataset("state", data=np.asarray(particles.state), compression="gzip")
        final.create_dataset("bunch_index", data=distribution.bunch_index, compression="gzip")


def run_rf_only(
    input_path: str | Path | None,
    output_path: str | Path,
    *,
    machine: RecyclerConfig | None = None,
    ramp: RFProgramConfig | None = None,
    max_particles: int | None = 20_000,
    record_every: int = 10,
    n_turns: int | None = None,
    seed: int = 202208,
    sampling_strata: int | None = None,
    initial_distribution: InitialDistribution | None = None,
    time_min_ns: float = -1000.0,
    time_max_ns: float = 3800.0,
    time_bin_ns: float = 2.0,
    show_progress: bool = True,
    invocation_argv: list[str] | None = None,
) -> Path:
    """Track a file-backed or generated sample and write RF-only diagnostics.

    Supply either ``input_path`` for a historical two-column file or
    ``initial_distribution`` for an already generated ensemble. With an
    explicit distribution, ``max_particles`` must be ``None``, zero, or its
    existing particle count because this function does not resample generated
    coordinates. The RF
    voltages are evaluated at the start of every turn, then Xsuite applies the
    two cavity kicks followed by longitudinal slippage.  Weighted time profiles
    are streamed at ``record_every`` intervals and written with final particle
    coordinates and full provenance to ``output_path``. A CLI caller can pass
    ``invocation_argv`` so the HDF5 file records the exact argument vector in
    addition to the resolved machine, RF, sampling, and tracking configuration.

    Returns
    -------
    pathlib.Path
        The HDF5 output path.
    """
    machine = machine or RecyclerConfig()
    ramp = ramp or RFProgramConfig()
    program = RFProgram(ramp)
    requested_n_turns = n_turns
    n_turns = ramp.turns_to_end(machine) if n_turns is None else int(n_turns)
    if n_turns <= 0:
        raise ValueError("n_turns must be positive")
    if time_max_ns <= time_min_ns or time_bin_ns <= 0:
        raise ValueError("Invalid profile time range or bin width")

    # Exact divisibility prevents a rounded bin count from silently changing
    # the requested diagnostic bin width.
    n_bins_float = (time_max_ns - time_min_ns) / time_bin_ns
    n_bins = int(round(n_bins_float))
    if not np.isclose(n_bins, n_bins_float):
        raise ValueError("Profile range must be an integer multiple of the bin width")
    time_edges_ns = np.linspace(time_min_ns, time_max_ns, n_bins + 1)

    # Capture executable code and environment before tracking starts. This is
    # especially useful when a long run is launched from a dirty worktree.
    run_provenance = collect_run_provenance()
    run_configuration = {
        "schema_version": 1,
        "machine": asdict(machine),
        "rf_program": asdict(ramp),
        "tracking": {
            "requested_n_turns": (
                None if requested_n_turns is None else int(requested_n_turns)
            ),
            "resolved_n_turns": n_turns,
            "record_every": int(record_every),
            "random_seed": int(seed),
            "time_min_ns": float(time_min_ns),
            "time_max_ns": float(time_max_ns),
            "time_bin_ns": float(time_bin_ns),
        },
        "input_request": {
            "input_path": None if input_path is None else str(input_path),
            "max_particles": (
                None if max_particles is None else int(max_particles)
            ),
            "sampling_strata": (
                None if sampling_strata is None else int(sampling_strata)
            ),
            "initial_distribution_supplied": initial_distribution is not None,
        },
    }

    if initial_distribution is None:
        if input_path is None:
            raise ValueError("input_path is required when no generated distribution is supplied")
        distribution = load_initial_distribution(
            input_path,
            machine,
            max_particles=max_particles,
            seed=seed,
            sampling_strata=sampling_strata,
        )
    else:
        if input_path is not None:
            raise ValueError("Supply either input_path or initial_distribution, not both")
        distribution = initial_distribution
        if max_particles not in (None, 0, distribution.n_particles):
            raise ValueError(
                "Generated distributions must be created at the desired particle count; "
                "max_particles cannot subsample them"
            )
    line = build_line(machine)
    particles = make_particles(distribution, machine)
    turns_to_record = _record_turns(n_turns, record_every)
    # Store only weighted profiles during tracking rather than a potentially
    # enormous particle-by-turn history.
    profiles = np.empty((turns_to_record.size, n_bins), dtype=np.float32)
    outside = np.empty(turns_to_record.size, dtype=np.float64)
    profiles[0], outside[0] = _profile(
        particles, machine, time_edges_ns, distribution.macro_equivalent_weight
    )

    if show_progress:
        print(
            f"Tracking {distribution.n_particles:,} particles for {n_turns:,} turns "
            f"({turns_to_record.size:,} recorded profiles)"
        )

    start = monotonic()
    record_index = 1
    next_report = 10
    for turn in range(n_turns):
        # RF voltages are sampled at the start of the turn. The saved state at
        # completed_turn is therefore the result after that kick and one slip.
        time_s = turn * machine.revolution_period_s
        voltage_53, voltage_25 = program.voltages(time_s)
        line["rf53"].voltage = voltage_53
        line["rf2p5"].voltage = voltage_25
        line.track(particles, num_turns=1)

        completed_turn = turn + 1
        if record_index < turns_to_record.size and completed_turn == turns_to_record[record_index]:
            profiles[record_index], outside[record_index] = _profile(
                particles,
                machine,
                time_edges_ns,
                distribution.macro_equivalent_weight,
            )
            record_index += 1

        percent = 100 * completed_turn // n_turns
        if show_progress and percent >= next_report:
            print(f"  {percent:3d}% ({completed_turn:,}/{n_turns:,} turns)")
            next_report += 10

    elapsed_s = monotonic() - start
    # Include particles outside the plotted window so this check tests tracking
    # loss and diagnostic bookkeeping rather than the chosen histogram limits.
    represented = profiles.sum(axis=1, dtype=np.float64) + outside
    if not np.allclose(
        represented,
        distribution.source_rows,
        rtol=1e-6,
        atol=1e-3,
    ):
        raise RuntimeError(
            "Recorded profile population is not conserved; "
            f"range={represented.min()}..{represented.max()}"
        )
    output_path = Path(output_path)
    _write_results(
        output_path,
        machine,
        ramp,
        distribution,
        particles,
        turns_to_record,
        time_edges_ns,
        profiles,
        outside,
        n_turns,
        record_every,
        seed,
        elapsed_s,
        run_provenance,
        run_configuration,
        invocation_argv,
    )
    if show_progress:
        print(f"Wrote {output_path} after {elapsed_s:.1f} s of tracking")
    return output_path
