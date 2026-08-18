"""Machine and RF-program constants with their provenance made explicit."""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil, sqrt


C_LIGHT_M_PER_S = 299_792_458.0
# BLonD constructs the proton mass from scipy.constants.m_p*c**2/e. This is
# the value from the SciPy environment used for the source audit; it differs
# from Xpart's built-in proton mass by 1.27 eV.
PROTON_MASS_EV = 938_272_089.4282572


@dataclass(frozen=True)
class RecyclerConfig:
    """Recycler reference values used by the historical BLonD simulation.

    Distances are in metres, momenta and energies are in eV-based Xsuite
    units, and frequencies are derived from the reference proton.  The
    ``input_time_denominator_m_per_s`` and ``blond_energy_beta`` fields retain
    the approximations used when the historical input file was produced; they
    are intentionally distinct from the more precise ``C_LIGHT_M_PER_S`` and
    computed :attr:`beta0` used by Xsuite.
    """

    circumference_m: float = 3319.41882346
    p0c_ev: float = 8.83532e9
    gamma_transition: float = 20.257643837730637
    harmonic_53: int = 588
    harmonic_25: int = 28
    intensity_protons: float = 1.05e12
    n_bunches: int = 8
    input_time_denominator_m_per_s: float = 3.0e8
    blond_energy_beta: float = 0.9944

    @property
    def energy0_ev(self) -> float:
        """Return the reference proton's total energy in eV."""
        return sqrt(self.p0c_ev**2 + PROTON_MASS_EV**2)

    @property
    def gamma0(self) -> float:
        """Return the reference proton's relativistic Lorentz factor."""
        return self.energy0_ev / PROTON_MASS_EV

    @property
    def beta0(self) -> float:
        """Return the reference proton's speed as a fraction of light speed."""
        return self.p0c_ev / self.energy0_ev

    @property
    def momentum_compaction_factor(self) -> float:
        """Return the first-order momentum compaction ``alpha_c=1/gamma_t^2``."""
        return 1.0 / self.gamma_transition**2

    @property
    def slip_factor(self) -> float:
        """Return ``eta=alpha_c-1/gamma0^2`` in the convention used here."""
        return self.momentum_compaction_factor - 1.0 / self.gamma0**2

    @property
    def revolution_period_s(self) -> float:
        """Return the reference proton's one-turn travel time in seconds."""
        return self.circumference_m / (self.beta0 * C_LIGHT_M_PER_S)

    @property
    def revolution_frequency_hz(self) -> float:
        """Return the reciprocal of the reference revolution period in hertz."""
        return 1.0 / self.revolution_period_s

    @property
    def rf53_frequency_hz(self) -> float:
        """Return the h=588 RF frequency, approximately 53 MHz."""
        return self.harmonic_53 * self.revolution_frequency_hz

    @property
    def rf25_frequency_hz(self) -> float:
        """Return the h=28 RF frequency, approximately 2.5 MHz."""
        return self.harmonic_25 * self.revolution_frequency_hz


@dataclass(frozen=True)
class RFProgramConfig:
    """Sequential h=588 turn-off and h=28 capture-voltage program.

    The default h=28 starting voltage is 3 kV, matching the publication and
    ESME RF program.  Keegan's committed BLonD simulation code used 5 kV;
    that legacy-code sensitivity remains available by constructing this class
    with ``rf25_initial_v=5_000.0`` or through the command line's
    ``--rf25-start-kv 5`` option.
    """

    rf53_initial_v: float = 80_000.0
    rf53_turnoff_s: float = 0.005
    rf25_start_s: float = 0.005
    rf25_ramp_duration_s: float = 0.085
    rf25_initial_v: float = 3_000.0
    rf25_final_v: float = 80_000.0

    @property
    def end_s(self) -> float:
        """Return the time when the h=28 capture ramp ends, in seconds."""
        return self.rf25_start_s + self.rf25_ramp_duration_s

    def turns_to_end(self, machine: RecyclerConfig) -> int:
        """Return turns needed for an end-of-turn state to reach ``end_s``.

        The ceiling is intentional: the final saved state is the first turn
        boundary at or beyond the requested RF-program duration.
        """
        return ceil(self.end_s / machine.revolution_period_s)
