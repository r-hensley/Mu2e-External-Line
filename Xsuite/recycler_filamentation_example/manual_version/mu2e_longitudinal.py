"""Mu2e Recycler longitudinal rebunching model for the USPAS live demonstration.

This module is a transparent, classroom-scale reimplementation of the
single-particle Recycler model documented in Mu2e-doc-2771.  It keeps the
kick/drift structure used in ``Longitudinal-Dynamics.ipynb`` while extending
the model to a common ring-phase coordinate with simultaneous h=588 and h=28
RF systems and a time-dependent voltage program.

The implementation is intentionally explicit:

* reference parameters come from Werkema, Mu2e-doc-2771, Table 3;
* the RF program follows Appendix B, including the 3-to-80 kV iso-adiabatic
  h=28 ramp;
* one 1/28-ring cell contains 21 independently colored h=588 bunchlets;
* the initial h=588 bunchlets approximate ESME ``KIND=13`` with a matched
  parabolic phase projection and a 0.12 eV-s limiting contour;
* only selected snapshots and scalar diagnostics are retained.

It is not ESME and must not be used as an operational extinction prediction.
The source ESME run used 1,075,200 particles and implementation details that
are not reproduced point-for-point here.  The classroom model also omits
collective effects, as did the Recycler model described in the 2019 report.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import Iterable, Sequence
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ellipk


MODEL_VERSION = "2026.07.23"
TWO_PI = 2.0 * math.pi


@dataclass(frozen=True)
class RecyclerParameters:
    """Published Recycler parameters and RF timing for one Mu2e bunch cell."""

    reference_energy_eV: float = 8_884.642e6
    proton_rest_energy_eV: float = 938.2720813e6
    transition_gamma: float = 21.6
    revolution_frequency_hz: float = 89.8095e3
    circumference_m: float = 3319.43

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

    bunch_emittance_eVs: float = 0.12
    source_bunch_count: int = 21
    in_time_half_window_s: float = 125.0e-9
    random_seed: int = 20190205

    @cached_property
    def gamma(self) -> float:
        return self.reference_energy_eV / self.proton_rest_energy_eV

    @cached_property
    def beta(self) -> float:
        return math.sqrt(1.0 - 1.0 / self.gamma**2)

    @cached_property
    def reference_momentum_eV_c(self) -> float:
        return math.sqrt(
            self.reference_energy_eV**2 - self.proton_rest_energy_eV**2
        )

    @cached_property
    def slip_factor(self) -> float:
        return 1.0 / self.transition_gamma**2 - 1.0 / self.gamma**2

    @cached_property
    def revolution_period_s(self) -> float:
        return 1.0 / self.revolution_frequency_hz

    @cached_property
    def angular_revolution_frequency(self) -> float:
        return TWO_PI * self.revolution_frequency_hz

    @cached_property
    def porch_duration_s(self) -> float:
        return self.porch_turns * self.revolution_period_s

    @cached_property
    def time_53_off_s(self) -> float:
        return self.porch_duration_s + self.turnoff_53_s

    @cached_property
    def time_25_flat_s(self) -> float:
        return self.time_53_off_s + self.ramp_25_s

    @cached_property
    def cell_period_rad(self) -> float:
        return TWO_PI / self.harmonic_25

    @cached_property
    def source_spacing_rad(self) -> float:
        return TWO_PI / self.harmonic_53

    @cached_property
    def extraction_times_s(self) -> np.ndarray:
        return self.time_25_flat_s + self.extraction_interval_s * np.arange(
            self.n_extractions, dtype=float
        )

    @cached_property
    def source_ids(self) -> np.ndarray:
        half = self.source_bunch_count // 2
        return np.arange(-half, half + 1, dtype=int)

    def parameter_table(self) -> pd.DataFrame:
        """Return a compact table with source values and derived checks."""

        rows = [
            ("Total synchronous energy", self.reference_energy_eV / 1e6, "MeV", "Werkema Table 3"),
            ("Reference momentum", self.reference_momentum_eV_c / 1e6, "MeV/c", "derived; Appendix B gives 8834.96"),
            ("Revolution frequency", self.revolution_frequency_hz / 1e3, "kHz", "Werkema Table 3"),
            ("Transition gamma", self.transition_gamma, "", "Werkema Table 3"),
            ("Slip factor", self.slip_factor, "", "derived; report gives -0.00901"),
            ("53 MHz harmonic", self.harmonic_53, "", "Werkema Table 3"),
            (
                "h=588 RF frequency",
                self.harmonic_53 * self.revolution_frequency_hz / 1e6,
                "MHz",
                "derived",
            ),
            ("2.5 MHz harmonic", self.harmonic_25, "", "Werkema Table 3"),
            (
                "h=28 RF frequency",
                self.harmonic_25 * self.revolution_frequency_hz / 1e6,
                "MHz",
                "derived",
            ),
            ("53 MHz voltage", self.voltage_53_eV / 1e3, "kV", "Werkema Table 3"),
            ("2.5 MHz initial voltage", self.voltage_25_initial_eV / 1e3, "kV", "Appendix B"),
            ("2.5 MHz flattop voltage", self.voltage_25_final_eV / 1e3, "kV", "Werkema Table 3"),
            ("53 MHz turn-off", self.turnoff_53_s * 1e3, "ms", "linear"),
            ("2.5 MHz ramp", self.ramp_25_s * 1e3, "ms", "iso-adiabatic"),
            ("53 MHz bunch emittance", self.bunch_emittance_eVs, "eV-s", "ESME SBNCH"),
        ]
        return pd.DataFrame(rows, columns=["quantity", "value", "unit", "provenance"])


@dataclass(frozen=True)
class InitialEnsemble:
    """Initial 21-bunch h=588 ensemble in one h=28 cell."""

    theta_rad: np.ndarray
    delta_energy_eV: np.ndarray
    source_id: np.ndarray
    local_theta_rad: np.ndarray
    contour_phase_amplitude_rad: float
    contour_area_eVs: float
    particles_per_source: int
    seed: int

    @property
    def n_particles(self) -> int:
        return int(self.theta_rad.size)


@dataclass(frozen=True)
class SimulationResult:
    """Streaming tracking output: selected states plus scalar diagnostics."""

    params: RecyclerParameters
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
    linearized_25: bool = False

    @property
    def n_particles(self) -> int:
        return int(self.source_id.size)

    def nearest_snapshot_index(self, time_s: float) -> int:
        return int(np.argmin(np.abs(self.snapshot_times_s - float(time_s))))

    def state_at(self, time_s: float) -> tuple[np.ndarray, np.ndarray, float]:
        index = self.nearest_snapshot_index(time_s)
        return (
            self.theta_rad[index].astype(float, copy=False),
            self.delta_energy_eV[index].astype(float, copy=False),
            float(self.snapshot_times_s[index]),
        )


def wrap_phase(theta_rad: np.ndarray | float, period_rad: float) -> np.ndarray:
    """Wrap phase into ``[-period/2, period/2)``."""

    theta = np.asarray(theta_rad)
    return (theta + 0.5 * period_rad) % period_rad - 0.5 * period_rad


def rf_voltages(
    time_s: np.ndarray | float,
    params: RecyclerParameters | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the document-faithful h=588 and h=28 voltage programs in volts.

    The h=588 source has one porch turn followed by a 5 ms linear ramp to zero.
    The h=28 source appears at 3 kV when the h=588 ramp ends and follows ESME's
    constant-adiabaticity curve to 80 kV in 85 ms.
    """

    p = params or RecyclerParameters()
    time = np.asarray(time_s, dtype=float)

    v53 = np.zeros_like(time)
    on_porch = time < p.porch_duration_s
    on_turnoff = (time >= p.porch_duration_s) & (time < p.time_53_off_s)
    v53[on_porch] = p.voltage_53_eV
    if np.any(on_turnoff):
        fraction = (
            time[on_turnoff] - p.porch_duration_s
        ) / p.turnoff_53_s
        v53[on_turnoff] = p.voltage_53_eV * (1.0 - fraction)

    v25 = np.zeros_like(time)
    on_ramp = (time >= p.time_53_off_s) & (time < p.time_25_flat_s)
    on_flat = time >= p.time_25_flat_s
    if np.any(on_ramp):
        u = (time[on_ramp] - p.time_53_off_s) / p.ramp_25_s
        ratio_root = math.sqrt(
            p.voltage_25_initial_eV / p.voltage_25_final_eV
        )
        denominator = 1.0 + u * (ratio_root - 1.0)
        v25[on_ramp] = p.voltage_25_initial_eV / denominator**2
    v25[on_flat] = p.voltage_25_final_eV
    return v53, v25


def rf_stage(time_s: float, params: RecyclerParameters | None = None) -> str:
    p = params or RecyclerParameters()
    if time_s < p.porch_duration_s:
        return "one-turn 53 MHz porch"
    if time_s < p.time_53_off_s:
        return "53 MHz turn-off"
    if time_s < p.time_25_flat_s:
        return "2.5 MHz iso-adiabatic ramp"
    return "2.5 MHz flattop"


def small_amplitude_tune(
    harmonic: int,
    voltage_eV: float,
    params: RecyclerParameters | None = None,
) -> float:
    """Small-amplitude synchrotron tune for a stationary sinusoidal bucket."""

    p = params or RecyclerParameters()
    argument = (
        -harmonic
        * p.slip_factor
        * float(voltage_eV)
        / (TWO_PI * p.beta**2 * p.reference_energy_eV)
    )
    return math.sqrt(max(argument, 0.0))


def small_amplitude_period_s(
    harmonic: int,
    voltage_eV: float,
    params: RecyclerParameters | None = None,
) -> float:
    p = params or RecyclerParameters()
    tune = small_amplitude_tune(harmonic, voltage_eV, p)
    return math.inf if tune == 0.0 else p.revolution_period_s / tune


def iso_adiabaticity(params: RecyclerParameters | None = None) -> float:
    """Return ``(1/omega_i - 1/omega_f)/T`` for the h=28 ramp."""

    p = params or RecyclerParameters()
    omega_i = TWO_PI * p.revolution_frequency_hz * small_amplitude_tune(
        p.harmonic_25, p.voltage_25_initial_eV, p
    )
    omega_f = TWO_PI * p.revolution_frequency_hz * small_amplitude_tune(
        p.harmonic_25, p.voltage_25_final_eV, p
    )
    return (1.0 / omega_i - 1.0 / omega_f) / p.ramp_25_s


def bucket_half_height_eV(
    harmonic: int,
    voltage_eV: float,
    params: RecyclerParameters | None = None,
) -> float:
    p = params or RecyclerParameters()
    return math.sqrt(
        -2.0
        * p.beta**2
        * p.reference_energy_eV
        * float(voltage_eV)
        / (math.pi * harmonic * p.slip_factor)
    )


def bucket_area_eVs(
    harmonic: int,
    voltage_eV: float,
    params: RecyclerParameters | None = None,
) -> float:
    p = params or RecyclerParameters()
    height = bucket_half_height_eV(harmonic, voltage_eV, p)
    return 8.0 * height / (
        harmonic * TWO_PI * p.revolution_frequency_hz
    )


def separatrix_energy_eV(
    theta_rad: np.ndarray,
    harmonic: int,
    voltage_eV: float,
    params: RecyclerParameters | None = None,
) -> np.ndarray:
    """Positive branch of a stationary-bucket separatrix."""

    p = params or RecyclerParameters()
    theta = np.asarray(theta_rad, dtype=float)
    half_height = bucket_half_height_eV(harmonic, voltage_eV, p)
    phase = harmonic * theta
    branch = half_height * np.cos(0.5 * phase)
    valid = np.abs(phase) <= math.pi
    return np.where(valid, np.maximum(branch, 0.0), np.nan)


def rf_waveform_eV(
    theta_rad: np.ndarray,
    time_s: float,
    params: RecyclerParameters | None = None,
    *,
    linearized_25: bool = False,
) -> np.ndarray:
    p = params or RecyclerParameters()
    theta = np.asarray(theta_rad, dtype=float)
    v53, v25 = rf_voltages(float(time_s), p)
    kick25 = (
        float(v25) * p.harmonic_25 * theta
        if linearized_25
        else float(v25) * np.sin(p.harmonic_25 * theta)
    )
    return float(v53) * np.sin(p.harmonic_53 * theta) + kick25


def _phase_slip_fractional_momentum(
    delta_energy_eV: np.ndarray,
    params: RecyclerParameters,
) -> np.ndarray:
    total_energy = params.reference_energy_eV + delta_energy_eV
    if np.any(total_energy <= params.proton_rest_energy_eV):
        raise ValueError("particle total energy fell below the proton rest energy")
    momentum = np.sqrt(
        total_energy**2 - params.proton_rest_energy_eV**2
    )
    return (
        momentum - params.reference_momentum_eV_c
    ) / params.reference_momentum_eV_c


def one_turn_map(
    theta_rad: np.ndarray,
    delta_energy_eV: np.ndarray,
    time_s: float,
    params: RecyclerParameters | None = None,
    *,
    linearized_25: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply one ESME-style drift-then-kick turn.

    ``LGRTHM=2`` in the published input uses the simplified
    ``2*pi*eta*Delta p/p`` phase slip.  The map is a composition of a drift
    depending only on energy and a kick depending only on phase, so its
    Jacobian determinant is one even when the RF voltage varies with time.
    """

    p = params or RecyclerParameters()
    theta = np.asarray(theta_rad, dtype=float)
    delta_energy = np.asarray(delta_energy_eV, dtype=float)

    fractional_momentum = _phase_slip_fractional_momentum(delta_energy, p)
    drifted_theta = wrap_phase(
        theta + TWO_PI * p.slip_factor * fractional_momentum,
        p.cell_period_rad,
    )
    kick = rf_waveform_eV(
        drifted_theta,
        time_s,
        p,
        linearized_25=linearized_25,
    )
    return drifted_theta, delta_energy + kick


def contour_energy_eV(
    rf_phase_rad: np.ndarray | float,
    phase_amplitude_rad: float,
    harmonic: int,
    voltage_eV: float,
    params: RecyclerParameters | None = None,
) -> np.ndarray:
    """Positive energy branch of a stationary-bucket contour."""

    p = params or RecyclerParameters()
    phase = np.asarray(rf_phase_rad, dtype=float)
    coefficient = TWO_PI * p.slip_factor / (
        p.beta**2 * p.reference_energy_eV
    )
    radicand = (
        2.0
        * float(voltage_eV)
        / (harmonic * abs(coefficient))
        * (np.cos(phase) - math.cos(phase_amplitude_rad))
    )
    return np.sqrt(np.maximum(radicand, 0.0))


def matched_contour_area_eVs(
    phase_amplitude_rad: float,
    harmonic: int,
    voltage_eV: float,
    params: RecyclerParameters | None = None,
) -> float:
    """Area enclosed by one sinusoidal RF contour in eV-s."""

    p = params or RecyclerParameters()

    def integrand(phase: float) -> float:
        return float(
            contour_energy_eV(
                phase,
                phase_amplitude_rad,
                harmonic,
                voltage_eV,
                p,
            )
        )

    integral, _ = quad(
        integrand,
        0.0,
        float(phase_amplitude_rad),
        epsabs=1e-10,
        epsrel=2e-11,
        limit=200,
    )
    return 4.0 * integral / (
        harmonic * p.angular_revolution_frequency
    )


def solve_contour_phase_amplitude(
    target_area_eVs: float,
    harmonic: int,
    voltage_eV: float,
    params: RecyclerParameters | None = None,
) -> float:
    """Find the RF-phase amplitude enclosing ``target_area_eVs``."""

    p = params or RecyclerParameters()
    target = float(target_area_eVs)
    bucket_area = bucket_area_eVs(harmonic, voltage_eV, p)
    if not 0.0 < target < bucket_area:
        raise ValueError(
            f"target area must lie between 0 and the {bucket_area:.6g} eV-s bucket area"
        )

    def residual(amplitude: float) -> float:
        return matched_contour_area_eVs(
            amplitude, harmonic, voltage_eV, p
        ) - target

    return float(brentq(residual, 1e-10, math.pi - 1e-10, xtol=1e-13))


def sample_matched_parabolic_bunch(
    n_particles: int,
    *,
    harmonic: int,
    voltage_eV: float,
    contour_area_eVs: float,
    rng: np.random.Generator,
    params: RecyclerParameters | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Sample an invariant parabolic bunch inside a limiting RF contour.

    The phase-space density is proportional to ``sqrt(1-I/Imax)``.  For a
    small-amplitude elliptical bucket this gives a parabolic phase projection,
    matching the documented intent of ESME ``KIND=13``.
    """

    p = params or RecyclerParameters()
    count = int(n_particles)
    if count <= 0:
        raise ValueError("n_particles must be positive")

    amplitude = solve_contour_phase_amplitude(
        contour_area_eVs, harmonic, voltage_eV, p
    )
    theta_limit = amplitude / harmonic
    energy_limit = float(
        contour_energy_eV(
            0.0, amplitude, harmonic, voltage_eV, p
        )
    )
    invariant_limit = 1.0 - math.cos(amplitude)
    coefficient = (
        harmonic
        * abs(
            TWO_PI
            * p.slip_factor
            / (p.beta**2 * p.reference_energy_eV)
        )
        / (2.0 * float(voltage_eV))
    )

    theta_parts: list[np.ndarray] = []
    energy_parts: list[np.ndarray] = []
    accepted = 0
    while accepted < count:
        remaining = count - accepted
        batch = max(2048, 5 * remaining)
        theta_candidate = rng.uniform(-theta_limit, theta_limit, batch)
        energy_candidate = rng.uniform(-energy_limit, energy_limit, batch)
        phase_candidate = harmonic * theta_candidate
        invariant = (
            coefficient * energy_candidate**2
            + 1.0
            - np.cos(phase_candidate)
        )
        depth = np.clip(
            1.0 - invariant / invariant_limit,
            0.0,
            1.0,
        )
        keep = (
            (invariant <= invariant_limit)
            & (rng.random(batch) < np.sqrt(depth))
        )
        theta_kept = theta_candidate[keep][:remaining]
        energy_kept = energy_candidate[keep][:remaining]
        theta_parts.append(theta_kept)
        energy_parts.append(energy_kept)
        accepted += theta_kept.size

    return (
        np.concatenate(theta_parts)[:count],
        np.concatenate(energy_parts)[:count],
        amplitude,
    )


def make_initial_ensemble(
    particles_per_source: int = 256,
    params: RecyclerParameters | None = None,
    *,
    seed: int | None = None,
) -> InitialEnsemble:
    """Create the deterministically colored 21-bunch initial distribution."""

    p = params or RecyclerParameters()
    if p.source_bunch_count != p.harmonic_53 // p.harmonic_25:
        raise ValueError("source_bunch_count must equal harmonic_53 / harmonic_25")

    used_seed = p.random_seed if seed is None else int(seed)
    rng = np.random.default_rng(used_seed)
    theta_all: list[np.ndarray] = []
    local_all: list[np.ndarray] = []
    energy_all: list[np.ndarray] = []
    source_all: list[np.ndarray] = []
    amplitude = math.nan

    for source in p.source_ids:
        local_theta, energy, amplitude = sample_matched_parabolic_bunch(
            particles_per_source,
            harmonic=p.harmonic_53,
            voltage_eV=p.voltage_53_eV,
            contour_area_eVs=p.bunch_emittance_eVs,
            rng=rng,
            params=p,
        )
        theta_all.append(
            wrap_phase(
                local_theta + source * p.source_spacing_rad,
                p.cell_period_rad,
            )
        )
        local_all.append(local_theta)
        energy_all.append(energy)
        source_all.append(
            np.full(int(particles_per_source), source, dtype=np.int16)
        )

    return InitialEnsemble(
        theta_rad=np.concatenate(theta_all),
        delta_energy_eV=np.concatenate(energy_all),
        source_id=np.concatenate(source_all),
        local_theta_rad=np.concatenate(local_all),
        contour_phase_amplitude_rad=float(amplitude),
        contour_area_eVs=p.bunch_emittance_eVs,
        particles_per_source=int(particles_per_source),
        seed=used_seed,
    )


def centered_bucket_coordinates(
    theta_rad: np.ndarray,
    params: RecyclerParameters | None = None,
) -> tuple[np.ndarray, float]:
    """Return circularly centered ring phase and its centroid."""

    p = params or RecyclerParameters()
    phase = p.harmonic_25 * np.asarray(theta_rad, dtype=float)
    mean_phase = math.atan2(
        float(np.mean(np.sin(phase))),
        float(np.mean(np.cos(phase))),
    )
    centered = np.angle(np.exp(1j * (phase - mean_phase))) / p.harmonic_25
    return centered, mean_phase / p.harmonic_25


def bunch_diagnostics(
    theta_rad: np.ndarray,
    delta_energy_eV: np.ndarray,
    params: RecyclerParameters | None = None,
) -> dict[str, float]:
    """Bunch widths, covariance area, and a fixed-phase ±125 ns fraction."""

    p = params or RecyclerParameters()
    theta = wrap_phase(theta_rad, p.cell_period_rad)
    time = theta / p.angular_revolution_frequency
    energy = np.asarray(delta_energy_eV, dtype=float)
    centered_energy = energy - np.mean(energy)
    covariance = np.cov(
        np.vstack((time, centered_energy)),
        ddof=0,
    )
    determinant = max(float(np.linalg.det(covariance)), 0.0)
    q_low, q_high = np.quantile(time, [0.025, 0.975])
    return {
        "sigma_time_s": float(np.std(time)),
        "sigma_energy_eV": float(np.std(centered_energy)),
        "covariance_emittance_eVs": math.sqrt(determinant),
        "central_95_time_width_s": float(q_high - q_low),
        "outside_window_fraction": float(
            np.mean(np.abs(time) > p.in_time_half_window_s)
        ),
        "centroid_time_s": float(np.mean(time)),
    }


def display_subset_indices(
    source_id: np.ndarray,
    per_source: int = 120,
) -> np.ndarray:
    """Deterministic approximately uniform display sample from every source."""

    source = np.asarray(source_id)
    indices: list[np.ndarray] = []
    for source_value in np.unique(source):
        candidates = np.flatnonzero(source == source_value)
        if candidates.size <= per_source:
            indices.append(candidates)
            continue
        positions = np.linspace(
            0, candidates.size - 1, int(per_source), dtype=int
        )
        indices.append(candidates[positions])
    return np.concatenate(indices)


def default_snapshot_times(
    params: RecyclerParameters | None = None,
    *,
    include_flattop_animation: bool = True,
    include_extractions: bool = True,
) -> np.ndarray:
    """Key-time and waterfall cadence used by the live demonstration."""

    p = params or RecyclerParameters()
    waterfall_step = 134 * p.revolution_period_s
    times = [
        0.0,
        p.porch_duration_s,
        p.time_53_off_s,
        p.time_25_flat_s,
    ]
    times.extend(
        np.arange(
            0.0,
            p.time_25_flat_s + 0.5 * waterfall_step,
            waterfall_step,
        )
    )
    if include_flattop_animation:
        times.extend(
            np.arange(
                p.time_25_flat_s + waterfall_step,
                p.time_25_flat_s + 3.0 * small_amplitude_period_s(
                    p.harmonic_25, p.voltage_25_final_eV, p
                ),
                waterfall_step,
            )
        )
    if include_extractions:
        times.extend(p.extraction_times_s)
    return np.unique(np.asarray(times, dtype=float))


def track_rebunching(
    ensemble: InitialEnsemble,
    *,
    params: RecyclerParameters | None = None,
    stop_time_s: float | None = None,
    snapshot_times_s: Sequence[float] | None = None,
    diagnostic_stride_turns: int = 20,
    display_per_source: int = 120,
    linearized_25: bool = False,
) -> SimulationResult:
    """Track the ensemble without allocating a particle-by-turn history cube."""

    p = params or RecyclerParameters()
    stop_time = (
        float(p.extraction_times_s[-1])
        if stop_time_s is None
        else float(stop_time_s)
    )
    if stop_time <= 0.0:
        raise ValueError("stop_time_s must be positive")
    stride = int(diagnostic_stride_turns)
    if stride <= 0:
        raise ValueError("diagnostic_stride_turns must be positive")

    requested_times = (
        default_snapshot_times(p)
        if snapshot_times_s is None
        else np.asarray(snapshot_times_s, dtype=float)
    )
    n_turns = int(np.ceil(stop_time / p.revolution_period_s))
    requested_turns = np.rint(
        np.clip(requested_times, 0.0, stop_time) / p.revolution_period_s
    ).astype(int)
    snapshot_turns = np.unique(
        np.concatenate(([0, n_turns], requested_turns))
    )
    snapshot_lookup = {
        int(turn): index for index, turn in enumerate(snapshot_turns)
    }

    n_particles = ensemble.n_particles
    theta_snapshots = np.empty(
        (snapshot_turns.size, n_particles), dtype=np.float32
    )
    energy_snapshots = np.empty_like(theta_snapshots)

    theta = ensemble.theta_rad.astype(float, copy=True)
    energy = ensemble.delta_energy_eV.astype(float, copy=True)
    theta_snapshots[snapshot_lookup[0]] = theta
    energy_snapshots[snapshot_lookup[0]] = energy

    diagnostic_turns: list[int] = []
    diagnostic_rows: list[dict[str, float]] = []

    def record_diagnostic(turn: int) -> None:
        diagnostic_turns.append(turn)
        diagnostic_rows.append(bunch_diagnostics(theta, energy, p))

    record_diagnostic(0)

    for turn in range(1, n_turns + 1):
        time_s = turn * p.revolution_period_s
        theta, energy = one_turn_map(
            theta,
            energy,
            time_s,
            p,
            linearized_25=linearized_25,
        )
        if turn in snapshot_lookup:
            index = snapshot_lookup[turn]
            theta_snapshots[index] = theta
            energy_snapshots[index] = energy
        if turn % stride == 0 or turn == n_turns:
            record_diagnostic(turn)

    diagnostic_frame = pd.DataFrame(diagnostic_rows)
    diagnostic_turn_array = np.asarray(diagnostic_turns, dtype=int)
    return SimulationResult(
        params=p,
        source_id=ensemble.source_id.copy(),
        snapshot_turns=snapshot_turns,
        snapshot_times_s=snapshot_turns * p.revolution_period_s,
        theta_rad=theta_snapshots,
        delta_energy_eV=energy_snapshots,
        display_indices=display_subset_indices(
            ensemble.source_id, display_per_source
        ),
        diagnostic_turns=diagnostic_turn_array,
        diagnostic_times_s=diagnostic_turn_array * p.revolution_period_s,
        sigma_time_s=diagnostic_frame["sigma_time_s"].to_numpy(),
        sigma_energy_eV=diagnostic_frame["sigma_energy_eV"].to_numpy(),
        covariance_emittance_eVs=diagnostic_frame[
            "covariance_emittance_eVs"
        ].to_numpy(),
        central_95_time_width_s=diagnostic_frame[
            "central_95_time_width_s"
        ].to_numpy(),
        outside_window_fraction=diagnostic_frame[
            "outside_window_fraction"
        ].to_numpy(),
        centroid_time_s=diagnostic_frame["centroid_time_s"].to_numpy(),
        linearized_25=bool(linearized_25),
    )


def amplitude_dependent_period_s(
    theta_amplitude_rad: np.ndarray | float,
    params: RecyclerParameters | None = None,
    *,
    voltage_eV: float | None = None,
) -> np.ndarray:
    """Exact pendulum period for an h=28 stationary sinusoidal bucket."""

    p = params or RecyclerParameters()
    voltage = (
        p.voltage_25_final_eV if voltage_eV is None else float(voltage_eV)
    )
    amplitude = np.asarray(theta_amplitude_rad, dtype=float)
    rf_amplitude = np.abs(p.harmonic_25 * amplitude)
    parameter = np.sin(0.5 * rf_amplitude) ** 2
    ratio = 2.0 * ellipk(parameter) / math.pi
    return small_amplitude_period_s(
        p.harmonic_25, voltage, p
    ) * ratio


def source_colors(
    source_id: np.ndarray,
    params: RecyclerParameters | None = None,
    *,
    cmap_name: str = "turbo",
) -> np.ndarray:
    p = params or RecyclerParameters()
    normalized = (
        np.asarray(source_id, dtype=float) - p.source_ids.min()
    ) / (p.source_ids.max() - p.source_ids.min())
    return mpl.colormaps[cmap_name](normalized)


def plot_rf_program(
    params: RecyclerParameters | None = None,
    *,
    end_time_ms: float = 105.0,
) -> plt.Figure:
    p = params or RecyclerParameters()
    time = np.linspace(0.0, end_time_ms * 1e-3, 1200)
    v53, v25 = rf_voltages(time, p)
    fig, ax = plt.subplots(figsize=(10.5, 4.2), constrained_layout=True)
    ax.plot(time * 1e3, v53 / 1e3, color="#3156a6", lw=2.4, label="53 MHz, h=588")
    ax.plot(time * 1e3, v25 / 1e3, color="#c43c39", lw=2.4, label="2.5 MHz, h=28")
    ax.axvline(p.time_53_off_s * 1e3, color="0.35", lw=1.0)
    ax.axvline(p.time_25_flat_s * 1e3, color="0.35", lw=1.0)
    ax.text(
        p.time_53_off_s * 1e3,
        84.0,
        "53 MHz off",
        ha="center",
        va="bottom",
        fontsize=9,
    )
    ax.text(
        p.time_25_flat_s * 1e3,
        84.0,
        "2.5 MHz flattop",
        ha="center",
        va="bottom",
        fontsize=9,
    )
    ax.set(
        xlim=(0.0, end_time_ms),
        ylim=(0.0, 90.0),
        xlabel="Time after second-batch injection [ms]",
        ylabel="RF sum voltage [kV]",
        title="Recycler RF voltage program",
    )
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, loc="center right")
    return fig


def _phase_space_axes(ax: plt.Axes, params: RecyclerParameters) -> None:
    half_cell_deg = math.degrees(math.pi / params.harmonic_25)
    ax.set_xlim(-half_cell_deg * 1.02, half_cell_deg * 1.02)
    ax.set_ylim(-45.0, 45.0)
    ax.set_xlabel(r"Recycler longitudinal phase $\theta$ [deg]")
    ax.set_ylabel(r"$\Delta E$ [MeV]")
    ax.grid(alpha=0.18)


def _scatter_state(
    ax: plt.Axes,
    theta_rad: np.ndarray,
    energy_eV: np.ndarray,
    source_id: np.ndarray,
    display_indices: np.ndarray,
    params: RecyclerParameters,
    *,
    size: float = 4.0,
    alpha: float = 0.72,
) -> mpl.collections.PathCollection:
    idx = np.asarray(display_indices, dtype=int)
    return ax.scatter(
        np.degrees(theta_rad[idx]),
        energy_eV[idx] / 1e6,
        c=source_id[idx],
        s=size,
        alpha=alpha,
        cmap="turbo",
        vmin=params.source_ids.min() - 0.5,
        vmax=params.source_ids.max() + 0.5,
        linewidths=0.0,
        rasterized=True,
    )


def _plot_h28_separatrix(
    ax: plt.Axes,
    voltage_eV: float,
    params: RecyclerParameters,
    *,
    label: str | None = None,
) -> None:
    if voltage_eV <= 0.0:
        return
    theta = np.linspace(
        -math.pi / params.harmonic_25,
        math.pi / params.harmonic_25,
        700,
    )
    energy = separatrix_energy_eV(
        theta, params.harmonic_25, voltage_eV, params
    )
    ax.plot(
        np.degrees(theta),
        energy / 1e6,
        color="#c92f2f",
        lw=1.6,
        label=label,
    )
    ax.plot(
        np.degrees(theta),
        -energy / 1e6,
        color="#c92f2f",
        lw=1.6,
    )


def plot_initial_ensemble(
    ensemble: InitialEnsemble,
    params: RecyclerParameters | None = None,
) -> plt.Figure:
    p = params or RecyclerParameters()
    fig = plt.figure(figsize=(12.5, 5.0), constrained_layout=True)
    grid = fig.add_gridspec(1, 3, width_ratios=(2.3, 1.0, 1.0))
    ax_phase = fig.add_subplot(grid[0, 0])
    ax_time = fig.add_subplot(grid[0, 1])
    ax_energy = fig.add_subplot(grid[0, 2])

    display_idx = display_subset_indices(
        ensemble.source_id,
        min(220, ensemble.particles_per_source),
    )
    scatter = _scatter_state(
        ax_phase,
        ensemble.theta_rad,
        ensemble.delta_energy_eV,
        ensemble.source_id,
        display_idx,
        p,
        size=4.0,
        alpha=0.78,
    )
    _phase_space_axes(ax_phase, p)
    ax_phase.set_title("21 matched 53 MHz bunchlets at injection")
    for center in p.source_ids * p.source_spacing_rad:
        ax_phase.axvline(
            math.degrees(center), color="0.75", lw=0.35, zorder=0
        )
    cbar = fig.colorbar(scatter, ax=ax_phase, pad=0.015)
    cbar.set_label("Original 53 MHz bunch index")
    cbar.set_ticks([-10, -5, 0, 5, 10])

    time_ns = (
        ensemble.theta_rad / p.angular_revolution_frequency * 1e9
    )
    ax_time.hist(
        time_ns,
        bins=150,
        histtype="stepfilled",
        color="#3a6ea5",
        alpha=0.65,
    )
    ax_time.set(
        xlabel="Time relative to h=28 bucket center [ns]",
        ylabel="Macroparticles / bin",
        title="Time projection",
    )
    ax_time.grid(alpha=0.18)

    ax_energy.hist(
        ensemble.delta_energy_eV / 1e6,
        bins=100,
        histtype="stepfilled",
        orientation="horizontal",
        color="#7b5aa6",
        alpha=0.65,
    )
    ax_energy.set(
        xlabel="Macroparticles / bin",
        ylabel=r"$\Delta E$ [MeV]",
        title="Energy projection",
    )
    ax_energy.grid(alpha=0.18)
    return fig


def plot_checkpoint_snapshots(
    result: SimulationResult,
    times_ms: Sequence[float] | None = None,
) -> plt.Figure:
    p = result.params
    requested = (
        [
            0.0,
            p.time_53_off_s * 1e3,
            25.0,
            55.0,
            p.time_25_flat_s * 1e3,
        ]
        if times_ms is None
        else [float(value) for value in times_ms]
    )
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(13.0, 8.2),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    axes_flat = axes.ravel()
    scatter = None
    for ax, time_ms in zip(axes_flat, requested):
        index = result.nearest_snapshot_index(time_ms * 1e-3)
        theta = result.theta_rad[index]
        energy = result.delta_energy_eV[index]
        time_s = float(result.snapshot_times_s[index])
        scatter = _scatter_state(
            ax,
            theta,
            energy,
            result.source_id,
            result.display_indices,
            p,
            size=3.0,
            alpha=0.73,
        )
        _, v25 = rf_voltages(time_s, p)
        _plot_h28_separatrix(ax, float(v25), p)
        _phase_space_axes(ax, p)
        ax.set_title(
            f"{time_s*1e3:.2f} ms\n{rf_stage(time_s, p)}",
            fontsize=10,
        )
    for ax in axes_flat[len(requested):]:
        ax.axis("off")
    if scatter is not None:
        cbar = fig.colorbar(
            scatter,
            ax=list(axes_flat[: len(requested)]),
            shrink=0.86,
            pad=0.012,
        )
        cbar.set_label("Original 53 MHz bunch index")
        cbar.set_ticks([-10, -5, 0, 5, 10])
    fig.suptitle(
        "From 53 MHz bunchlets to a filamented 2.5 MHz bunch",
        fontsize=14,
    )
    return fig


def plot_snapshot_dashboard(
    result: SimulationResult,
    snapshot_index: int,
) -> plt.Figure:
    """Linked state view used by the notebook's time slider."""

    p = result.params
    index = int(np.clip(snapshot_index, 0, result.snapshot_times_s.size - 1))
    time_s = float(result.snapshot_times_s[index])
    theta = result.theta_rad[index].astype(float, copy=False)
    energy = result.delta_energy_eV[index].astype(float, copy=False)
    _, v25 = rf_voltages(time_s, p)

    fig = plt.figure(figsize=(12.5, 8.3), constrained_layout=True)
    grid = fig.add_gridspec(
        3,
        3,
        height_ratios=(0.55, 2.0, 1.1),
        width_ratios=(1.45, 1.45, 1.0),
    )
    ax_program = fig.add_subplot(grid[0, :])
    ax_phase = fig.add_subplot(grid[1:, :2])
    ax_time = fig.add_subplot(grid[1, 2])
    ax_rf = fig.add_subplot(grid[2, 2])

    program_time = np.linspace(
        0.0,
        max(p.time_25_flat_s * 1.08, time_s * 1.03),
        900,
    )
    v53_program, v25_program = rf_voltages(program_time, p)
    ax_program.plot(
        program_time * 1e3,
        v53_program / 1e3,
        color="#3156a6",
        lw=1.8,
        label="53 MHz",
    )
    ax_program.plot(
        program_time * 1e3,
        v25_program / 1e3,
        color="#c43c39",
        lw=1.8,
        label="2.5 MHz",
    )
    ax_program.axvline(time_s * 1e3, color="0.1", lw=1.4)
    ax_program.set(
        ylabel="RF [kV]",
        xlim=(0.0, program_time[-1] * 1e3),
        ylim=(0.0, 88.0),
    )
    ax_program.grid(alpha=0.2)
    ax_program.legend(frameon=False, ncol=2, loc="upper center")

    scatter = _scatter_state(
        ax_phase,
        theta,
        energy,
        result.source_id,
        result.display_indices,
        p,
        size=4.0,
        alpha=0.76,
    )
    _plot_h28_separatrix(
        ax_phase,
        float(v25),
        p,
        label="instantaneous h=28 separatrix" if float(v25) > 0 else None,
    )
    _phase_space_axes(ax_phase, p)
    if float(v25) > 0:
        ax_phase.legend(frameon=False, loc="upper right", fontsize=8)
    cbar = fig.colorbar(scatter, ax=ax_phase, pad=0.012, shrink=0.82)
    cbar.set_label("Original 53 MHz bunch index")
    cbar.set_ticks([-10, -5, 0, 5, 10])

    time_ns = (
        wrap_phase(theta, p.cell_period_rad)
        / p.angular_revolution_frequency
        * 1e9
    )
    ax_time.hist(
        time_ns,
        bins=100,
        color="#3a6ea5",
        alpha=0.68,
        histtype="stepfilled",
    )
    ax_time.axvline(-125.0, color="#c92f2f", lw=1.1, ls="--")
    ax_time.axvline(125.0, color="#c92f2f", lw=1.1, ls="--")
    ax_time.set(
        xlim=(-205.0, 205.0),
        xlabel="Time from synchronous phase [ns]",
        ylabel="Macroparticles / bin",
        title="Time projection",
    )
    ax_time.grid(alpha=0.18)

    theta_grid = np.linspace(
        -math.pi / p.harmonic_25,
        math.pi / p.harmonic_25,
        2000,
    )
    waveform = rf_waveform_eV(
        theta_grid,
        time_s,
        p,
        linearized_25=result.linearized_25,
    )
    ax_rf.plot(
        np.degrees(theta_grid),
        waveform / 1e3,
        color="#3b8f6b",
        lw=1.2,
    )
    ax_rf.axhline(0.0, color="0.45", lw=0.7)
    ax_rf.set(
        xlim=(
            -math.degrees(math.pi / p.harmonic_25),
            math.degrees(math.pi / p.harmonic_25),
        ),
        xlabel=r"$\theta$ [deg]",
        ylabel="Net RF kick [keV]",
        title="Instantaneous RF waveform",
    )
    ax_rf.grid(alpha=0.18)

    diagnostics = bunch_diagnostics(theta, energy, p)
    fig.suptitle(
        f"{time_s*1e3:.2f} ms — {rf_stage(time_s, p)}"
        f" | σt={diagnostics['sigma_time_s']*1e9:.1f} ns"
        f" | σE={diagnostics['sigma_energy_eV']/1e6:.1f} MeV",
        fontsize=13,
    )
    return fig


def plot_preloaded_snapshot_dashboard(
    result: SimulationResult,
    *,
    display_per_source: int = 120,
    phase_bins: int = 100,
    frame_duration_ms: int = 140,
) -> go.Figure:
    """Build a browser-side Plotly animation containing every saved snapshot.

    Tracking, histogramming, RF waveforms, and separatrices are all computed
    before the figure is returned. Plotly slider changes then select serialized
    frames without invoking Python or rebuilding a Matplotlib figure.
    """

    p = result.params
    if result.snapshot_times_s.size == 0:
        raise ValueError("result must contain at least one saved snapshot")
    if int(display_per_source) <= 0:
        raise ValueError("display_per_source must be positive")
    if int(phase_bins) < 20:
        raise ValueError("phase_bins must be at least 20")
    if int(frame_duration_ms) < 0:
        raise ValueError("frame_duration_ms must be nonnegative")

    display_indices = display_subset_indices(
        result.source_id,
        per_source=int(display_per_source),
    )
    source_id = result.source_id[display_indices]
    phase_deg = np.degrees(
        result.theta_rad[:, display_indices]
    ).astype(np.float32)
    energy_mev = (
        result.delta_energy_eV[:, display_indices] / 1e6
    ).astype(np.float32)

    time_edges_ns = np.linspace(-205.0, 205.0, int(phase_bins) + 1)
    time_centers_ns = (
        0.5 * (time_edges_ns[:-1] + time_edges_ns[1:])
    ).astype(np.float32)
    time_profiles = np.empty(
        (result.snapshot_times_s.size, int(phase_bins)),
        dtype=np.float32,
    )
    for index, theta in enumerate(result.theta_rad):
        time_ns = (
            wrap_phase(theta, p.cell_period_rad)
            / p.angular_revolution_frequency
            * 1e9
        )
        time_profiles[index] = np.histogram(
            time_ns,
            bins=time_edges_ns,
        )[0]

    theta_grid = np.linspace(
        -math.pi / p.harmonic_25,
        math.pi / p.harmonic_25,
        420,
    )
    theta_grid_deg = np.degrees(theta_grid).astype(np.float32)
    waveform_kev = np.empty(
        (result.snapshot_times_s.size, theta_grid.size),
        dtype=np.float32,
    )
    separatrix_mev = np.full_like(waveform_kev, np.nan)
    _, voltage_25 = rf_voltages(result.snapshot_times_s, p)
    for index, time_s in enumerate(result.snapshot_times_s):
        waveform_kev[index] = (
            rf_waveform_eV(
                theta_grid,
                float(time_s),
                p,
                linearized_25=result.linearized_25,
            )
            / 1e3
        )
        if float(voltage_25[index]) > 0.0 and not result.linearized_25:
            separatrix_mev[index] = (
                separatrix_energy_eV(
                    theta_grid,
                    p.harmonic_25,
                    float(voltage_25[index]),
                    p,
                )
                / 1e6
            )

    figure = make_subplots(
        rows=2,
        cols=2,
        specs=[[{"rowspan": 2}, {}], [None, {}]],
        row_heights=[0.58, 0.42],
        column_widths=[0.68, 0.32],
        horizontal_spacing=0.12,
        vertical_spacing=0.15,
        subplot_titles=(
            "Longitudinal phase space",
            "Time projection",
            "Instantaneous RF waveform",
        ),
    )

    phase_trace_index = len(figure.data)
    figure.add_trace(
        # Plotly's redraw=False path repaints SVG Scatter traces but not
        # Scattergl, so this type is required for responsive preloaded frames.
        go.Scatter(
            x=phase_deg[0],
            y=energy_mev[0],
            mode="markers",
            marker={
                "size": 3.5,
                "opacity": 0.76,
                "color": source_id,
                "colorscale": "Turbo",
                "cmin": p.source_ids.min() - 0.5,
                "cmax": p.source_ids.max() + 0.5,
                "colorbar": {
                    "title": "Original 53 MHz<br>bunch index",
                    "x": 0.655,
                    "y": 0.50,
                    "len": 0.72,
                    "thickness": 13,
                    "tickvals": [-10, -5, 0, 5, 10],
                },
            },
            customdata=source_id,
            name="tracked macroparticles",
            showlegend=False,
            hovertemplate=(
                "source bunch = %{customdata}<br>"
                "θ = %{x:.4f} deg<br>"
                "ΔE = %{y:.3f} MeV<extra></extra>"
            ),
        ),
        row=1,
        col=1,
    )

    separatrix_upper_index = len(figure.data)
    figure.add_trace(
        go.Scatter(
            x=theta_grid_deg,
            y=separatrix_mev[0],
            mode="lines",
            line={"color": "#c92f2f", "width": 1.8},
            name="instantaneous h=28 separatrix",
            hoverinfo="skip",
        ),
        row=1,
        col=1,
    )
    separatrix_lower_index = len(figure.data)
    figure.add_trace(
        go.Scatter(
            x=theta_grid_deg,
            y=-separatrix_mev[0],
            mode="lines",
            line={"color": "#c92f2f", "width": 1.8},
            name="instantaneous h=28 separatrix",
            showlegend=False,
            hoverinfo="skip",
        ),
        row=1,
        col=1,
    )

    profile_trace_index = len(figure.data)
    figure.add_trace(
        go.Scatter(
            x=time_centers_ns,
            y=time_profiles[0],
            mode="lines",
            fill="tozeroy",
            line={"color": "#3a6ea5", "width": 1.4, "shape": "hv"},
            fillcolor="rgba(58, 110, 165, 0.28)",
            name="time projection",
            showlegend=False,
            hovertemplate=(
                "time = %{x:.1f} ns<br>"
                "macroparticles / bin = %{y:.0f}<extra></extra>"
            ),
        ),
        row=1,
        col=2,
    )

    waveform_trace_index = len(figure.data)
    figure.add_trace(
        go.Scatter(
            x=theta_grid_deg,
            y=waveform_kev[0],
            mode="lines",
            line={"color": "#3b8f6b", "width": 1.5},
            name="net RF kick",
            showlegend=False,
            hovertemplate=(
                "θ = %{x:.3f} deg<br>"
                "net kick = %{y:.2f} keV<extra></extra>"
            ),
        ),
        row=2,
        col=2,
    )

    figure.add_vline(
        x=-p.in_time_half_window_s * 1e9,
        line={"color": "#c92f2f", "dash": "dash", "width": 1.1},
        row=1,
        col=2,
    )
    figure.add_vline(
        x=p.in_time_half_window_s * 1e9,
        line={"color": "#c92f2f", "dash": "dash", "width": 1.1},
        row=1,
        col=2,
    )
    figure.add_hline(
        y=0.0,
        line={"color": "rgba(80, 80, 80, 0.55)", "width": 0.8},
        row=2,
        col=2,
    )

    dynamic_trace_indices = [
        phase_trace_index,
        separatrix_upper_index,
        separatrix_lower_index,
        profile_trace_index,
        waveform_trace_index,
    ]
    frames: list[go.Frame] = []
    for index, time_s in enumerate(result.snapshot_times_s):
        frames.append(
            go.Frame(
                name=str(index),
                data=[
                    go.Scatter(
                        x=phase_deg[index],
                        y=energy_mev[index],
                    ),
                    go.Scatter(y=separatrix_mev[index]),
                    go.Scatter(y=-separatrix_mev[index]),
                    go.Scatter(y=time_profiles[index]),
                    go.Scatter(y=waveform_kev[index]),
                ],
                traces=dynamic_trace_indices,
            )
        )
    figure.frames = frames

    slider_steps = [
        {
            "method": "animate",
            "args": [
                [frame.name],
                {
                    "mode": "immediate",
                    "frame": {"duration": 0, "redraw": False},
                    "transition": {"duration": 0},
                },
            ],
            "label": (
                f"time = {time_s*1e3:.3f} ms"
                f" — {rf_stage(float(time_s), p)}"
            ),
        }
        for frame, time_s in zip(
            frames,
            result.snapshot_times_s,
            strict=True,
        )
    ]

    phase_domain = tuple(figure.layout.xaxis.domain)
    figure.update_layout(
        title={
            "text": "Preloaded Mu2e longitudinal phase-space evolution",
            "x": 0.5,
            "xanchor": "center",
        },
        template="plotly_white",
        autosize=True,
        width=None,
        height=690,
        margin={"t": 105, "r": 95, "b": 105, "l": 70},
        legend={
            "orientation": "h",
            "x": phase_domain[0],
            "y": 1.01,
            "xanchor": "left",
            "yanchor": "bottom",
        },
        sliders=[
            {
                "active": 0,
                "x": phase_domain[0],
                "len": phase_domain[1] - phase_domain[0],
                "y": -0.05,
                "pad": {"t": 0, "b": 0},
                "font": {"size": 1, "color": "rgba(0,0,0,0)"},
                "ticklen": 0,
                "minorticklen": 0,
                "tickwidth": 0,
                "tickcolor": "rgba(0,0,0,0)",
                "currentvalue": {
                    "prefix": "",
                    "suffix": "",
                    "font": {"size": 13, "color": "rgba(45,55,70,1)"},
                },
                "steps": slider_steps,
            }
        ],
        updatemenus=[
            {
                "type": "buttons",
                "direction": "left",
                "x": phase_domain[0],
                "y": 1.08,
                "showactive": False,
                "buttons": [
                    {
                        "label": "Play",
                        "method": "animate",
                        "args": [
                            None,
                            {
                                "fromcurrent": True,
                                "frame": {
                                    "duration": int(frame_duration_ms),
                                    "redraw": False,
                                },
                                "transition": {"duration": 0},
                            },
                        ],
                    },
                    {
                        "label": "Pause",
                        "method": "animate",
                        "args": [
                            [None],
                            {
                                "mode": "immediate",
                                "frame": {"duration": 0, "redraw": False},
                                "transition": {"duration": 0},
                            },
                        ],
                    },
                ],
            }
        ],
    )
    half_cell_deg = math.degrees(math.pi / p.harmonic_25)
    figure.update_xaxes(
        title_text="Recycler longitudinal phase θ [deg]",
        range=[-1.02 * half_cell_deg, 1.02 * half_cell_deg],
        row=1,
        col=1,
    )
    figure.update_yaxes(
        title_text="ΔE [MeV]",
        range=[-45.0, 45.0],
        row=1,
        col=1,
    )
    figure.update_xaxes(
        title_text="Time from synchronous phase [ns]",
        range=[-205.0, 205.0],
        row=1,
        col=2,
    )
    figure.update_yaxes(
        title_text="Macroparticles / bin",
        range=[0.0, 1.08 * float(np.max(time_profiles))],
        row=1,
        col=2,
    )
    figure.update_xaxes(
        title_text="Recycler longitudinal phase θ [deg]",
        range=[-half_cell_deg, half_cell_deg],
        row=2,
        col=2,
    )
    figure.update_yaxes(
        title_text="Net RF kick [keV]",
        range=[-88.0, 88.0],
        row=2,
        col=2,
    )
    return figure


def plot_waterfall(
    result: SimulationResult,
    *,
    start_time_s: float | None = None,
    end_time_s: float | None = None,
    cadence_turns: int = 134,
    bins: int = 150,
) -> plt.Figure:
    """Recreate the report's stacked time-profile waterfall."""

    p = result.params
    start = 0.0 if start_time_s is None else float(start_time_s)
    end = p.time_25_flat_s if end_time_s is None else float(end_time_s)
    cadence_s = int(cadence_turns) * p.revolution_period_s
    requested = np.unique(
        np.concatenate(
            (
                np.arange(start, end + 0.5 * cadence_s, cadence_s),
                np.array([start, end]),
            )
        )
    )
    indices = np.unique(
        [result.nearest_snapshot_index(time) for time in requested]
    )
    edges = np.linspace(
        -math.degrees(math.pi / p.harmonic_25),
        math.degrees(math.pi / p.harmonic_25),
        int(bins) + 1,
    )
    centers = 0.5 * (edges[:-1] + edges[1:])

    fig, ax = plt.subplots(figsize=(11.0, 8.0), constrained_layout=True)
    if len(indices) > 1:
        baseline_step = np.median(
            np.diff(result.snapshot_times_s[indices])
        ) * 1e3
    else:
        baseline_step = cadence_s * 1e3
    scale = 0.82 * baseline_step
    for index in indices:
        theta = result.theta_rad[index].astype(float, copy=False)
        counts, _ = np.histogram(np.degrees(theta), bins=edges)
        profile = counts / max(float(counts.max()), 1.0)
        baseline = result.snapshot_times_s[index] * 1e3
        ax.plot(
            centers,
            baseline + scale * profile,
            color="#4f9961",
            lw=0.75,
            alpha=0.82,
        )

    window_theta_deg = math.degrees(
        p.angular_revolution_frequency * p.in_time_half_window_s
    )
    ax.axvline(-window_theta_deg, color="#c92f2f", ls="--", lw=1.2)
    ax.axvline(window_theta_deg, color="#c92f2f", ls="--", lw=1.2)
    ax.axhline(
        p.time_53_off_s * 1e3,
        color="0.35",
        ls=":",
        lw=1.0,
    )
    ax.axhline(
        p.time_25_flat_s * 1e3,
        color="0.35",
        ls=":",
        lw=1.0,
    )
    ax.text(
        centers[-1],
        p.time_53_off_s * 1e3,
        " 53 MHz off",
        va="bottom",
        ha="right",
        fontsize=9,
    )
    ax.text(
        centers[-1],
        p.time_25_flat_s * 1e3,
        " 2.5 MHz flattop",
        va="bottom",
        ha="right",
        fontsize=9,
    )
    ax.set(
        xlim=(edges[0], edges[-1]),
        ylim=(start * 1e3 - baseline_step, end * 1e3 + 1.2 * scale),
        xlabel=r"Recycler longitudinal phase $\theta$ [deg]",
        ylabel="Time in RF program [ms]",
        title=(
            "Waterfall of the h=28-cell time projection\n"
            f"regular cadence: {cadence_turns} turns = "
            f"{cadence_s*1e3:.3f} ms; flattop endpoint added"
        ),
    )
    ax.grid(axis="x", alpha=0.15)
    return fig


def plot_tune_spread(
    params: RecyclerParameters | None = None,
) -> plt.Figure:
    p = params or RecyclerParameters()
    bucket_edge = math.pi / p.harmonic_25
    amplitude = np.linspace(0.0, 0.995 * bucket_edge, 700)
    period = amplitude_dependent_period_s(amplitude, p)
    tune_ratio = period[0] / period

    fig, axes = plt.subplots(
        1, 2, figsize=(11.5, 4.2), constrained_layout=True
    )
    axes[0].plot(
        np.degrees(amplitude),
        tune_ratio,
        color="#3156a6",
        lw=2.2,
    )
    axes[0].set(
        xlim=(0.0, math.degrees(bucket_edge)),
        ylim=(0.0, 1.03),
        xlabel=r"Oscillation amplitude $|\theta|_{\max}$ [deg]",
        ylabel=r"$f_s/f_{s0}$",
        title="Amplitude-dependent synchrotron frequency",
    )
    axes[0].grid(alpha=0.22)

    axes[1].plot(
        np.degrees(amplitude),
        period * 1e3,
        color="#8b4c9e",
        lw=2.2,
    )
    axes[1].axhline(
        small_amplitude_period_s(
            p.harmonic_25, p.voltage_25_final_eV, p
        )
        * 1e3,
        color="0.35",
        ls="--",
        lw=1.0,
        label="small-amplitude period",
    )
    axes[1].set(
        xlim=(0.0, math.degrees(bucket_edge)),
        ylim=(0.0, min(70.0, np.nanmax(period * 1e3))),
        xlabel=r"Oscillation amplitude $|\theta|_{\max}$ [deg]",
        ylabel="Synchrotron period [ms]",
        title="Outer particles fall behind",
    )
    axes[1].grid(alpha=0.22)
    axes[1].legend(frameon=False)
    return fig


def plot_nonlinear_control_comparison(
    nonlinear: SimulationResult,
    linearized: SimulationResult,
    *,
    time_s: float,
) -> plt.Figure:
    """Compare the physical sinusoidal RF with an amplitude-linear control."""

    if nonlinear.params != linearized.params:
        raise ValueError("comparison results must use the same parameters")
    if not linearized.linearized_25:
        raise ValueError("the control result must have linearized_25=True")

    p = nonlinear.params
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(11.5, 5.0),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    entries = [
        (
            axes[0],
            nonlinear,
            "Physical sinusoidal h=28 RF\namplitude-dependent tune",
        ),
        (
            axes[1],
            linearized,
            "Causal control: linear h=28 kick\none tune at every amplitude",
        ),
    ]
    scatter = None
    actual_times: list[float] = []
    for ax, result, title in entries:
        index = result.nearest_snapshot_index(time_s)
        actual_time = float(result.snapshot_times_s[index])
        actual_times.append(actual_time)
        scatter = _scatter_state(
            ax,
            result.theta_rad[index],
            result.delta_energy_eV[index],
            result.source_id,
            result.display_indices,
            p,
            size=3.2,
            alpha=0.72,
        )
        _phase_space_axes(ax, p)
        ax.set_title(title, fontsize=11)

    _plot_h28_separatrix(
        axes[0],
        p.voltage_25_final_eV,
        p,
        label="physical h=28 separatrix",
    )
    axes[0].legend(frameon=False, loc="upper right", fontsize=8)
    if scatter is not None:
        cbar = fig.colorbar(
            scatter,
            ax=list(axes),
            shrink=0.86,
            pad=0.018,
        )
        cbar.set_label("Original 53 MHz bunch index")
        cbar.set_ticks([-10, -5, 0, 5, 10])
    fig.suptitle(
        "Remove the tune spread and the progressive color shear disappears"
        f"\ncomparison at {np.mean(actual_times)*1e3:.2f} ms",
        fontsize=13,
    )
    return fig


def plot_diagnostics(
    result: SimulationResult,
    *,
    end_time_ms: float | None = None,
) -> plt.Figure:
    p = result.params
    time_ms = result.diagnostic_times_s * 1e3
    mask = (
        np.ones_like(time_ms, dtype=bool)
        if end_time_ms is None
        else time_ms <= float(end_time_ms)
    )
    time_ms = time_ms[mask]

    fig, axes = plt.subplots(
        4,
        1,
        figsize=(11.5, 9.2),
        sharex=True,
        constrained_layout=True,
    )
    axes[0].plot(
        time_ms,
        result.sigma_time_s[mask] * 1e9,
        color="#3156a6",
        lw=1.2,
    )
    axes[0].set_ylabel(r"$\sigma_t$ [ns]")

    axes[1].plot(
        time_ms,
        result.sigma_energy_eV[mask] / 1e6,
        color="#8b4c9e",
        lw=1.2,
    )
    axes[1].set_ylabel(r"$\sigma_E$ [MeV]")

    emittance0 = max(result.covariance_emittance_eVs[0], 1e-30)
    axes[2].plot(
        time_ms,
        result.covariance_emittance_eVs[mask] / emittance0,
        color="#3b8f6b",
        lw=1.2,
    )
    axes[2].set_ylabel(
        r"Centered $\epsilon_{\rm cov}/\epsilon_{\rm cov,0}$"
    )

    n_particles = result.n_particles
    fraction = result.outside_window_fraction[mask]
    floor = 1.0 / n_particles
    axes[3].plot(
        time_ms,
        fraction,
        color="#c43c39",
        lw=1.1,
    )
    axes[3].axhline(
        floor,
        color="0.4",
        ls="--",
        lw=0.9,
        label=f"one-macroparticle floor = {floor:.2g}",
    )
    axes[3].set(
        xlabel="Time after second-batch injection [ms]",
        ylabel="Fraction outside ±125 ns",
        ylim=(0.0, max(0.5, 1.15 * np.max(fraction))),
    )
    axes[3].legend(frameon=False, loc="upper right", fontsize=8)

    max_time = time_ms[-1] if time_ms.size else 0.0
    for ax in axes:
        ax.grid(alpha=0.2)
        ax.axvline(p.time_53_off_s * 1e3, color="0.45", ls=":", lw=0.8)
        ax.axvline(p.time_25_flat_s * 1e3, color="0.45", ls=":", lw=0.8)
        for extraction_time in p.extraction_times_s * 1e3:
            if extraction_time <= max_time + 1e-9:
                ax.axvline(
                    extraction_time,
                    color="#9a7b32",
                    alpha=0.28,
                    lw=0.7,
                )
    axes[0].set_title(
        "Projection breathing, covariance evolution, and finite-statistics tails"
    )
    return fig


def plot_extraction_gallery(
    result: SimulationResult,
) -> plt.Figure:
    p = result.params
    fig, axes = plt.subplots(
        2,
        4,
        figsize=(14.0, 7.0),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    scatter = None
    for number, (ax, time_s) in enumerate(
        zip(axes.ravel(), p.extraction_times_s), start=1
    ):
        index = result.nearest_snapshot_index(float(time_s))
        scatter = _scatter_state(
            ax,
            result.theta_rad[index],
            result.delta_energy_eV[index],
            result.source_id,
            result.display_indices,
            p,
            size=2.2,
            alpha=0.65,
        )
        _plot_h28_separatrix(ax, p.voltage_25_final_eV, p)
        _phase_space_axes(ax, p)
        diag = bunch_diagnostics(
            result.theta_rad[index],
            result.delta_energy_eV[index],
            p,
        )
        ax.set_title(
            f"Bunch {number}: {result.snapshot_times_s[index]*1e3:.1f} ms\n"
            f"σt={diag['sigma_time_s']*1e9:.1f} ns",
            fontsize=9,
        )
    if scatter is not None:
        cbar = fig.colorbar(
            scatter,
            ax=list(axes.ravel()),
            shrink=0.86,
            pad=0.012,
        )
        cbar.set_label("Original 53 MHz bunch index")
        cbar.set_ticks([-10, -5, 0, 5, 10])
    fig.suptitle(
        "Eight extraction times sample different phases of the residual motion",
        fontsize=14,
    )
    return fig


def validation_table(
    params: RecyclerParameters | None = None,
) -> pd.DataFrame:
    """Analytical/source anchors used as notebook acceptance checks."""

    p = params or RecyclerParameters()
    values = [
        (
            "h53 / h2.5",
            p.harmonic_53 / p.harmonic_25,
            21.0,
            "",
        ),
        (
            "1 degree of ring phase",
            p.revolution_period_s / 360.0 * 1e9,
            30.92,
            "ns",
        ),
        (
            "2.5 MHz flattop bucket half-height",
            bucket_half_height_eV(
                p.harmonic_25, p.voltage_25_final_eV, p
            )
            / 1e6,
            42.1,
            "MeV",
        ),
        (
            "2.5 MHz flattop small-amplitude period",
            small_amplitude_period_s(
                p.harmonic_25, p.voltage_25_final_eV, p
            )
            * 1e3,
            18.4,
            "ms",
        ),
        (
            "Iso-adiabatic ramp parameter",
            iso_adiabaticity(p),
            0.144,
            "",
        ),
        (
            "53 MHz full bucket area",
            bucket_area_eVs(
                p.harmonic_53, p.voltage_53_eV, p
            ),
            0.22159,
            "eV-s",
        ),
    ]
    frame = pd.DataFrame(
        values,
        columns=["check", "classroom model", "report/reference", "unit"],
    )
    frame["relative difference"] = (
        frame["classroom model"] - frame["report/reference"]
    ) / frame["report/reference"]
    return frame


def numerical_map_jacobian(
    theta_rad: float,
    delta_energy_eV: float,
    time_s: float,
    params: RecyclerParameters | None = None,
) -> np.ndarray:
    """Finite-difference Jacobian for one-turn symplecticity validation."""

    p = params or RecyclerParameters()
    dtheta = 2e-8
    denergy = 200.0

    def mapped(theta: float, energy: float) -> np.ndarray:
        theta_out, energy_out = one_turn_map(
            np.array([theta]),
            np.array([energy]),
            time_s,
            p,
        )
        return np.array([theta_out[0], energy_out[0]])

    plus_theta = mapped(theta_rad + dtheta, delta_energy_eV)
    minus_theta = mapped(theta_rad - dtheta, delta_energy_eV)
    plus_energy = mapped(theta_rad, delta_energy_eV + denergy)
    minus_energy = mapped(theta_rad, delta_energy_eV - denergy)
    return np.column_stack(
        (
            (plus_theta - minus_theta) / (2.0 * dtheta),
            (plus_energy - minus_energy) / (2.0 * denergy),
        )
    )


def measured_small_amplitude_period_s(
    params: RecyclerParameters | None = None,
    *,
    amplitude_deg: float = 0.05,
    periods_to_measure: int = 6,
) -> float:
    """Measure the flattop period from repeated positive-to-negative crossings."""

    p = params or RecyclerParameters()
    theta = np.array([math.radians(float(amplitude_deg))])
    energy = np.array([0.0])
    predicted_turns = int(
        math.ceil(
            small_amplitude_period_s(
                p.harmonic_25, p.voltage_25_final_eV, p
            )
            / p.revolution_period_s
        )
    )
    crossings: list[int] = []
    previous = float(theta[0])
    total_turns = predicted_turns * (periods_to_measure + 2)
    static_time = p.time_25_flat_s + p.revolution_period_s
    for turn in range(1, total_turns + 1):
        theta, energy = one_turn_map(theta, energy, static_time, p)
        current = float(theta[0])
        if previous > 0.0 and current <= 0.0:
            crossings.append(turn)
        previous = current
        if len(crossings) >= periods_to_measure + 1:
            break
    if len(crossings) < 2:
        raise RuntimeError("insufficient zero crossings to measure period")
    return float(np.mean(np.diff(crossings)) * p.revolution_period_s)


__all__ = [
    "MODEL_VERSION",
    "InitialEnsemble",
    "RecyclerParameters",
    "SimulationResult",
    "amplitude_dependent_period_s",
    "bucket_area_eVs",
    "bucket_half_height_eV",
    "bunch_diagnostics",
    "centered_bucket_coordinates",
    "contour_energy_eV",
    "default_snapshot_times",
    "display_subset_indices",
    "iso_adiabaticity",
    "make_initial_ensemble",
    "matched_contour_area_eVs",
    "measured_small_amplitude_period_s",
    "numerical_map_jacobian",
    "one_turn_map",
    "plot_checkpoint_snapshots",
    "plot_diagnostics",
    "plot_extraction_gallery",
    "plot_initial_ensemble",
    "plot_nonlinear_control_comparison",
    "plot_preloaded_snapshot_dashboard",
    "plot_rf_program",
    "plot_snapshot_dashboard",
    "plot_tune_spread",
    "plot_waterfall",
    "rf_stage",
    "rf_voltages",
    "rf_waveform_eV",
    "sample_matched_parabolic_bunch",
    "separatrix_energy_eV",
    "small_amplitude_period_s",
    "small_amplitude_tune",
    "solve_contour_phase_amplitude",
    "source_colors",
    "track_rebunching",
    "validation_table",
    "wrap_phase",
]
