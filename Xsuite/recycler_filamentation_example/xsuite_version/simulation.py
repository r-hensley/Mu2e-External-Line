"""Xsuite tracking model for the Recycler filamentation demonstration.

The reference implementation in the sibling ``manual_version`` directory
advances ``theta`` and ``Delta E`` with an explicit NumPy map.  This module instead creates Xsuite
``Particles`` and tracks them through a line containing a nonlinear
``LineSegmentMap`` and two harmonic ``Cavity`` elements.  The line ordering is
slip, h=588 cavity, h=28 cavity, matching the reference drift-then-kick
stroboscopic convention.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import math
from time import perf_counter
from typing import Sequence

import numpy as np
import xpart as xp
import xtrack as xt


C_LIGHT_M_PER_S = 299_792_458.0
TWO_PI = 2.0 * math.pi


# Keep supplied machine settings separate from quantities that Xsuite derives.
# This prevents an analytical helper value from silently replacing the value
# associated with the actual Xsuite reference particle or one-turn line.
@dataclass(frozen=True)
class RecyclerInputs:
    """Documented machine and RF inputs supplied to Xsuite."""

    circumference_m: float = 3319.43
    reference_energy_eV: float = 8_884.642e6
    transition_gamma: float = 21.6
    harmonic_53: int = 588
    harmonic_25: int = 28
    voltage_53_eV: float = 80.0e3
    voltage_25_initial_eV: float = 3.0e3
    voltage_25_final_eV: float = 80.0e3
    porch_turns: int = 1
    turnoff_53_s: float = 5.0e-3
    ramp_25_s: float = 85.0e-3
    extraction_interval_s: float = 48.12e-3
    n_extractions: int = 8
    source_bunch_count: int = 21
    in_time_half_window_s: float = 125.0e-9

    @property
    def momentum_compaction_factor(self) -> float:
        return 1.0 / self.transition_gamma**2


@dataclass(frozen=True)
class XsuiteDerived:
    """Reference and lattice quantities derived by Xsuite."""

    mass0_eV: float
    p0c_eV: float
    energy0_eV: float
    beta0: float
    gamma0: float
    line_length_m: float
    t_rev0_s: float
    revolution_frequency_hz: float
    momentum_compaction_factor: float
    slip_factor: float
    rf53_frequency_hz: float
    rf25_frequency_hz: float
    qs_80kv: float
    synchrotron_period_80kv_s: float
    time_53_off_s: float
    time_25_flat_s: float
    extraction_times_s: np.ndarray

    @property
    def angular_revolution_frequency(self) -> float:
        return TWO_PI * self.revolution_frequency_hz

    @property
    def cell_period_rad(self) -> float:
        return TWO_PI / 28.0


@dataclass(frozen=True)
class XsuiteSimulationResult:
    """Selected Xsuite states and scalar diagnostics."""

    inputs: RecyclerInputs
    derived: XsuiteDerived
    source_id: np.ndarray
    snapshot_turns: np.ndarray
    snapshot_times_s: np.ndarray
    theta_rad: np.ndarray
    delta_energy_eV: np.ndarray
    display_indices: np.ndarray
    diagnostic_turns: np.ndarray
    diagnostic_times_s: np.ndarray
    sigma_time_s: np.ndarray
    sigma_energy_eV: np.ndarray
    covariance_emittance_eVs: np.ndarray
    central_95_time_width_s: np.ndarray
    outside_window_fraction: np.ndarray
    centroid_time_s: np.ndarray
    final_zeta_m: np.ndarray
    final_ptau: np.ndarray
    final_state: np.ndarray
    tracking_elapsed_s: float

    @property
    def n_particles(self) -> int:
        return int(self.source_id.size)

    def nearest_snapshot_index(self, time_s: float) -> int:
        return int(np.argmin(np.abs(self.snapshot_times_s - float(time_s))))


def _scalar(value) -> float:
    return float(np.asarray(value).reshape(-1)[0])


def build_line(inputs: RecyclerInputs | None = None) -> xt.Line:
    """Build the Xsuite drift-then-kick one-turn line."""

    p = inputs or RecyclerInputs()

    # Element order is physically significant.  The preserved classroom map
    # first slips at the old energy and then applies both RF kicks, so the
    # Xsuite line must place the longitudinal map before the cavities.
    line = xt.Line(
        elements=[
            xt.LineSegmentMap(
                length=p.circumference_m,
                # Non-zero transverse tunes allow a complete 6D Twiss call.
                # All tracked particles start at zero transverse coordinates,
                # so this bookkeeping choice does not affect longitudinal motion.
                qx=0.31,
                qy=0.32,
                longitudinal_mode="nonlinear",
                momentum_compaction_factor=p.momentum_compaction_factor,
                slippage_length=p.circumference_m,
                # LineSegmentMap requires RF arrays in nonlinear mode.  These
                # zero entries disable its internal RF kick because the real
                # time-dependent kicks are supplied by the Cavity elements.
                voltage_rf=[0.0],
                frequency_rf=[0.0],
                lag_rf=[0.0],
            ),
            xt.Cavity(harmonic=p.harmonic_53, voltage=0.0, lag=0.0),
            xt.Cavity(harmonic=p.harmonic_25, voltage=0.0, lag=0.0),
        ],
        element_names=["recycler_slip", "rf53", "rf2p5"],
    )
    line.particle_ref = xt.Particles(
        mass0=xp.PROTON_MASS_EV,
        q0=1,
        energy0=p.reference_energy_eV,
    )
    line.build_tracker()
    return line


def derive_from_xsuite(
    line: xt.Line,
    inputs: RecyclerInputs | None = None,
) -> XsuiteDerived:
    """Ask Xsuite for the reference and flattop longitudinal quantities."""

    p = inputs or RecyclerInputs()

    # Twiss needs the final stationary bucket to obtain Qs.  Temporarily put
    # the h=28 cavity at flattop, query the Xsuite line, then restore zero volts
    # before the time-dependent tracking run begins.
    line["rf53"].voltage = 0.0
    line["rf2p5"].voltage = p.voltage_25_final_eV
    twiss = line.twiss(method="6d")
    line["rf2p5"].voltage = 0.0

    particle_ref = line.particle_ref
    t_rev0_s = float(twiss.t_rev0)
    time_53_off_s = p.porch_turns * t_rev0_s + p.turnoff_53_s
    time_25_flat_s = time_53_off_s + p.ramp_25_s
    extraction_times_s = time_25_flat_s + p.extraction_interval_s * np.arange(
        p.n_extractions,
        dtype=float,
    )
    qs = float(twiss.qs)
    return XsuiteDerived(
        mass0_eV=_scalar(particle_ref.mass0),
        p0c_eV=_scalar(particle_ref.p0c),
        energy0_eV=_scalar(particle_ref.energy0),
        beta0=_scalar(particle_ref.beta0),
        gamma0=_scalar(particle_ref.gamma0),
        line_length_m=float(line.get_length()),
        t_rev0_s=t_rev0_s,
        revolution_frequency_hz=1.0 / t_rev0_s,
        momentum_compaction_factor=float(twiss.momentum_compaction_factor),
        slip_factor=float(twiss.slip_factor),
        rf53_frequency_hz=p.harmonic_53 / t_rev0_s,
        rf25_frequency_hz=p.harmonic_25 / t_rev0_s,
        qs_80kv=qs,
        synchrotron_period_80kv_s=t_rev0_s / qs,
        time_53_off_s=time_53_off_s,
        time_25_flat_s=time_25_flat_s,
        extraction_times_s=extraction_times_s,
    )


def rf_voltages(
    time_s: float | np.ndarray,
    inputs: RecyclerInputs,
    derived: XsuiteDerived,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the one-turn porch, linear turn-off, and iso-adiabatic ramp."""

    time = np.asarray(time_s, dtype=float)

    # The single porch turn is based on Xsuite's revolution period rather than
    # the independently quoted report frequency.  This keeps every RF event on
    # the same time base as the line that advances the particles.
    porch_duration_s = inputs.porch_turns * derived.t_rev0_s

    v53 = np.zeros_like(time)
    on_porch = time < porch_duration_s
    on_turnoff = (time >= porch_duration_s) & (time < derived.time_53_off_s)
    v53[on_porch] = inputs.voltage_53_eV
    if np.any(on_turnoff):
        fraction = (time[on_turnoff] - porch_duration_s) / inputs.turnoff_53_s
        v53[on_turnoff] = inputs.voltage_53_eV * (1.0 - fraction)

    v25 = np.zeros_like(time)
    on_ramp = (time >= derived.time_53_off_s) & (time < derived.time_25_flat_s)
    on_flat = time >= derived.time_25_flat_s
    if np.any(on_ramp):
        u = (time[on_ramp] - derived.time_53_off_s) / inputs.ramp_25_s

        # Constant-adiabaticity voltage law.  It is intentionally nonlinear in
        # voltage even though the normalized ramp coordinate u is linear.
        root_ratio = math.sqrt(
            inputs.voltage_25_initial_eV / inputs.voltage_25_final_eV
        )
        denominator = 1.0 + u * (root_ratio - 1.0)
        v25[on_ramp] = inputs.voltage_25_initial_eV / denominator**2
    v25[on_flat] = inputs.voltage_25_final_eV
    return v53, v25


def rf_stage(time_s: float, inputs: RecyclerInputs, derived: XsuiteDerived) -> str:
    porch_duration_s = inputs.porch_turns * derived.t_rev0_s
    if time_s < porch_duration_s:
        return "one-turn 53 MHz porch"
    if time_s < derived.time_53_off_s:
        return "53 MHz turn-off"
    if time_s < derived.time_25_flat_s:
        return "2.5 MHz iso-adiabatic ramp"
    return "2.5 MHz flattop"


def wrap_phase(theta_rad: np.ndarray | float, period_rad: float) -> np.ndarray:
    theta = np.asarray(theta_rad)
    return (theta + 0.5 * period_rad) % period_rad - 0.5 * period_rad


def particles_to_local_coordinates(
    particles: xt.Particles,
    inputs: RecyclerInputs,
    derived: XsuiteDerived,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert Xsuite ``zeta, ptau`` to wrapped ring phase and energy offset."""

    # Xsuite uses zeta = s - beta*c*t.  The classroom coordinate uses positive
    # arrival-time phase theta, hence theta = -2*pi*zeta/C.  Wrapping is done
    # only when recording/plotting; Xsuite itself tracks the unwrapped zeta.
    theta = wrap_phase(
        -TWO_PI * np.asarray(particles.zeta) / inputs.circumference_m,
        derived.cell_period_rad,
    )

    # ptau is (E-E0)/p0c, so this conversion is exact for Xsuite particles and
    # does not use a small-delta energy approximation.
    delta_energy_eV = np.asarray(particles.ptau) * derived.p0c_eV
    return theta, delta_energy_eV


def make_particles(
    line: xt.Line,
    theta_rad: np.ndarray,
    delta_energy_eV: np.ndarray,
    inputs: RecyclerInputs,
    derived: XsuiteDerived,
) -> xt.Particles:
    # Invert the coordinate transformation above.  The initial distribution is
    # identical to the preserved NumPy example, isolating the tracking engine
    # as the only intended difference between the two simulations.
    zeta_m = -inputs.circumference_m * np.asarray(theta_rad) / TWO_PI
    ptau = np.asarray(delta_energy_eV) / derived.p0c_eV
    return line.build_particles(
        zeta=zeta_m,
        ptau=ptau,
        particle_id=np.arange(np.asarray(theta_rad).size, dtype=np.int64),
    )


def display_subset_indices(source_id: np.ndarray, per_source: int = 120) -> np.ndarray:
    indices: list[np.ndarray] = []
    for source_value in np.unique(source_id):
        candidates = np.flatnonzero(source_id == source_value)
        if candidates.size <= per_source:
            indices.append(candidates)
        else:
            positions = np.linspace(0, candidates.size - 1, per_source, dtype=int)
            indices.append(candidates[positions])
    return np.concatenate(indices)


def bunch_diagnostics(
    theta_rad: np.ndarray,
    delta_energy_eV: np.ndarray,
    inputs: RecyclerInputs,
    derived: XsuiteDerived,
) -> dict[str, float]:
    theta = wrap_phase(theta_rad, derived.cell_period_rad)
    time_s = theta / derived.angular_revolution_frequency
    energy = np.asarray(delta_energy_eV, dtype=float)
    centered_energy = energy - np.mean(energy)

    # Covariance emittance is a projection diagnostic, not the conserved
    # fine-grained phase-space area.  It is retained to reproduce the teaching
    # comparison made by the original example.
    covariance = np.cov(np.vstack((time_s, centered_energy)), ddof=0)
    determinant = max(float(np.linalg.det(covariance)), 0.0)
    q_low, q_high = np.quantile(time_s, [0.025, 0.975])
    return {
        "sigma_time_s": float(np.std(time_s)),
        "sigma_energy_eV": float(np.std(centered_energy)),
        "covariance_emittance_eVs": math.sqrt(determinant),
        "central_95_time_width_s": float(q_high - q_low),
        "outside_window_fraction": float(
            np.mean(np.abs(time_s) > inputs.in_time_half_window_s)
        ),
        "centroid_time_s": float(np.mean(time_s)),
    }


def default_snapshot_times(
    inputs: RecyclerInputs,
    derived: XsuiteDerived,
) -> np.ndarray:
    waterfall_step_s = 134 * derived.t_rev0_s

    # Combine named RF events, the report-style 134-turn waterfall cadence,
    # several flattop synchrotron periods, and all eight extraction samples.
    times = [
        0.0,
        inputs.porch_turns * derived.t_rev0_s,
        derived.time_53_off_s,
        derived.time_25_flat_s,
    ]
    times.extend(
        np.arange(
            0.0,
            derived.time_25_flat_s + 0.5 * waterfall_step_s,
            waterfall_step_s,
        )
    )
    times.extend(
        np.arange(
            derived.time_25_flat_s + waterfall_step_s,
            derived.time_25_flat_s + 3.0 * derived.synchrotron_period_80kv_s,
            waterfall_step_s,
        )
    )
    times.extend(derived.extraction_times_s)
    return np.unique(np.asarray(times, dtype=float))


def track_rebunching(
    line: xt.Line,
    initial_theta_rad: np.ndarray,
    initial_delta_energy_eV: np.ndarray,
    source_id: np.ndarray,
    inputs: RecyclerInputs,
    derived: XsuiteDerived,
    *,
    diagnostic_stride_turns: int = 20,
    display_per_source: int = 120,
    show_progress: bool = True,
) -> XsuiteSimulationResult:
    """Track the full eight-extraction history with Xsuite."""

    if diagnostic_stride_turns <= 0:
        raise ValueError("diagnostic_stride_turns must be positive")
    stop_time_s = float(derived.extraction_times_s[-1])
    n_turns = int(math.ceil(stop_time_s / derived.t_rev0_s))
    requested_times = default_snapshot_times(inputs, derived)
    requested_turns = np.rint(
        np.clip(requested_times, 0.0, stop_time_s) / derived.t_rev0_s
    ).astype(int)
    snapshot_turns = np.unique(np.concatenate(([0, n_turns], requested_turns)))
    snapshot_lookup = {
        int(turn): index for index, turn in enumerate(snapshot_turns)
    }

    particles = make_particles(
        line,
        initial_theta_rad,
        initial_delta_energy_eV,
        inputs,
        derived,
    )
    n_particles = int(np.asarray(initial_theta_rad).size)

    # Only selected frames are retained.  Float32 mirrors the preserved result
    # and keeps the output compact while Xsuite continues tracking in Float64.
    theta_snapshots = np.empty((snapshot_turns.size, n_particles), dtype=np.float32)
    energy_snapshots = np.empty_like(theta_snapshots)

    theta, energy = particles_to_local_coordinates(particles, inputs, derived)
    theta_snapshots[snapshot_lookup[0]] = theta
    energy_snapshots[snapshot_lookup[0]] = energy

    diagnostic_turns: list[int] = []
    diagnostic_rows: list[dict[str, float]] = []

    def record_diagnostic(turn: int) -> None:
        current_theta, current_energy = particles_to_local_coordinates(
            particles, inputs, derived
        )
        diagnostic_turns.append(turn)
        diagnostic_rows.append(
            bunch_diagnostics(current_theta, current_energy, inputs, derived)
        )

    record_diagnostic(0)
    start = perf_counter()
    next_report = 10
    for turn in range(1, n_turns + 1):
        time_s = turn * derived.t_rev0_s
        voltage_53, voltage_25 = rf_voltages(time_s, inputs, derived)

        # Updating element strengths per turn preserves the exact nonlinear RF
        # ramp while Xsuite performs every particle kick and slip operation.
        line["rf53"].voltage = float(voltage_53)
        line["rf2p5"].voltage = float(voltage_25)
        line.track(particles, num_turns=1)

        if turn in snapshot_lookup:
            theta, energy = particles_to_local_coordinates(particles, inputs, derived)
            index = snapshot_lookup[turn]
            theta_snapshots[index] = theta
            energy_snapshots[index] = energy
        if turn % diagnostic_stride_turns == 0 or turn == n_turns:
            record_diagnostic(turn)

        percent = 100 * turn // n_turns
        if show_progress and percent >= next_report:
            print(f"  Xsuite {percent:3d}% ({turn:,}/{n_turns:,} turns)")
            next_report += 10

    elapsed_s = perf_counter() - start
    rows = {key: np.asarray([row[key] for row in diagnostic_rows]) for key in diagnostic_rows[0]}
    state = np.asarray(particles.state).copy()
    if not np.all(state > 0):
        raise RuntimeError("one or more Xsuite particles were lost")
    if not np.isfinite(theta_snapshots).all() or not np.isfinite(energy_snapshots).all():
        raise RuntimeError("Xsuite tracking produced non-finite coordinates")

    diagnostic_turn_array = np.asarray(diagnostic_turns, dtype=int)
    return XsuiteSimulationResult(
        inputs=inputs,
        derived=derived,
        source_id=np.asarray(source_id).copy(),
        snapshot_turns=snapshot_turns,
        snapshot_times_s=snapshot_turns * derived.t_rev0_s,
        theta_rad=theta_snapshots,
        delta_energy_eV=energy_snapshots,
        display_indices=display_subset_indices(np.asarray(source_id), display_per_source),
        diagnostic_turns=diagnostic_turn_array,
        diagnostic_times_s=diagnostic_turn_array * derived.t_rev0_s,
        sigma_time_s=rows["sigma_time_s"],
        sigma_energy_eV=rows["sigma_energy_eV"],
        covariance_emittance_eVs=rows["covariance_emittance_eVs"],
        central_95_time_width_s=rows["central_95_time_width_s"],
        outside_window_fraction=rows["outside_window_fraction"],
        centroid_time_s=rows["centroid_time_s"],
        final_zeta_m=np.asarray(particles.zeta).copy(),
        final_ptau=np.asarray(particles.ptau).copy(),
        final_state=state,
        tracking_elapsed_s=elapsed_s,
    )


def validate_one_turn(
    line: xt.Line,
    inputs: RecyclerInputs,
    derived: XsuiteDerived,
) -> dict[str, float]:
    """Validate the Xsuite coordinate sign, slip, kick, and map ordering."""

    theta0 = math.radians(1.0)
    energy0 = 2.0e6
    particles = make_particles(
        line,
        np.asarray([theta0]),
        np.asarray([energy0]),
        inputs,
        derived,
    )
    delta0 = _scalar(particles.delta)
    time_s = 25.0e-3
    voltage_53, voltage_25 = rf_voltages(time_s, inputs, derived)
    line["rf53"].voltage = float(voltage_53)
    line["rf2p5"].voltage = float(voltage_25)

    # This independent expectation checks the coordinate sign and confirms
    # that the Xsuite line applies one full slip before the two cavity kicks.
    expected_theta = wrap_phase(
        theta0 + TWO_PI * derived.slip_factor * delta0,
        derived.cell_period_rad,
    )
    expected_energy = (
        energy0
        + float(voltage_53) * math.sin(inputs.harmonic_53 * float(expected_theta))
        + float(voltage_25) * math.sin(inputs.harmonic_25 * float(expected_theta))
    )
    line.track(particles, num_turns=1)
    actual_theta, actual_energy = particles_to_local_coordinates(
        particles, inputs, derived
    )
    line["rf53"].voltage = 0.0
    line["rf2p5"].voltage = 0.0
    return {
        "theta_error_rad": float(actual_theta[0] - expected_theta),
        "energy_error_eV": float(actual_energy[0] - expected_energy),
    }


def result_digest(result: XsuiteSimulationResult) -> str:
    """Hash all trajectory-defining arrays for deterministic rerun checks."""

    digest = sha256()
    arrays = (
        result.source_id,
        result.snapshot_turns,
        result.snapshot_times_s,
        result.theta_rad,
        result.delta_energy_eV,
        result.diagnostic_turns,
        result.sigma_time_s,
        result.sigma_energy_eV,
        result.covariance_emittance_eVs,
        result.central_95_time_width_s,
        result.outside_window_fraction,
        result.centroid_time_s,
        result.final_zeta_m,
        result.final_ptau,
        result.final_state,
    )
    for array in arrays:
        contiguous = np.ascontiguousarray(array)
        digest.update(contiguous.dtype.str.encode("ascii"))
        digest.update(str(contiguous.shape).encode("ascii"))
        digest.update(contiguous.view(np.uint8))
    return digest.hexdigest()
