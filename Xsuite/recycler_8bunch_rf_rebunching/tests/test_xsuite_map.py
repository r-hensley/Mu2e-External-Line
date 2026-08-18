import numpy as np

from xsuite_recycler.config import RecyclerConfig
from xsuite_recycler.distribution import dt_to_zeta
from xsuite_recycler.simulation import build_line


def test_rf_kick_phase_and_sign():
    machine = RecyclerConfig()
    line = build_line(machine)
    dt_s = 100e-9
    particles = line.build_particles(
        zeta=dt_to_zeta(dt_s, machine.beta0),
        ptau=0.0,
    )
    line["rf53"].voltage = 80_000.0
    line["rf2p5"].voltage = 0.0

    line.track(particles, num_turns=1)

    actual_kick_ev = float(particles.ptau[0]) * machine.p0c_ev
    expected_kick_ev = 80_000.0 * np.sin(2 * np.pi * machine.rf53_frequency_hz * dt_s)
    assert np.isclose(actual_kick_ev, expected_kick_ev, rtol=1e-11, atol=1e-7)


def test_zero_voltage_map_has_one_full_nonlinear_slip():
    machine = RecyclerConfig()
    line = build_line(machine)
    particles = line.build_particles(zeta=0.25, ptau=1e-3)
    initial_zeta = float(particles.zeta[0])
    initial_delta = float(particles.delta[0])

    line.track(particles, num_turns=1)

    expected_zeta = (
        initial_zeta
        - machine.slip_factor * machine.circumference_m * initial_delta
    )
    assert np.isclose(float(particles.zeta[0]), expected_zeta, atol=1e-14)
