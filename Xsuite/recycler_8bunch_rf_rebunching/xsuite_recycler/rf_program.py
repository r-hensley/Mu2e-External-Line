"""Time-dependent voltage functions for the two Recycler RF systems."""

from __future__ import annotations

import numpy as np

from .config import RFProgramConfig


class RFProgram:
    """Evaluate the sequential h=588 turn-off and h=28 capture ramp."""

    def __init__(self, config: RFProgramConfig | None = None):
        """Create a voltage evaluator, using the standard program by default."""
        self.config = config or RFProgramConfig()

    def voltage_53(self, time_s):
        """Return the h=588 voltage in volts at scalar or array times in seconds.

        The voltage falls linearly from its configured initial value to zero
        over the turn-off interval and remains zero afterward.
        """
        t = np.asarray(time_s, dtype=float)
        fraction = np.clip(1.0 - t / self.config.rf53_turnoff_s, 0.0, 1.0)
        result = self.config.rf53_initial_v * fraction
        return float(result) if result.ndim == 0 else result

    def voltage_25(self, time_s):
        """Return the h=28 voltage in volts at scalar or array times in seconds.

        This reproduces the BLonD inverse-square adiabatic ramp, shifted to the
        configured capture start.  The h=28 system is off before that time.
        """
        t = np.asarray(time_s, dtype=float)
        u = np.clip(
            (t - self.config.rf25_start_s)
            / self.config.rf25_ramp_duration_s,
            0.0,
            1.0,
        )
        v1 = self.config.rf25_initial_v
        v2 = self.config.rf25_final_v
        # This inverse-square interpolation is exactly v1 at u=0 and v2 at
        # u=1. Clipping u also holds the final voltage after the ramp ends.
        denominator = 1.0 - u * (np.sqrt(v2) - np.sqrt(v1)) / np.sqrt(v2)
        ramp = v1 / denominator**2
        # The clipped ramp would otherwise equal v1 before capture begins.
        result = np.where(t < self.config.rf25_start_s, 0.0, ramp)
        return float(result) if result.ndim == 0 else result

    def voltages(self, time_s):
        """Return ``(h=588, h=28)`` voltages for the supplied time or times."""
        return self.voltage_53(time_s), self.voltage_25(time_s)
