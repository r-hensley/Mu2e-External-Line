import numpy as np

from xsuite_recycler.config import RFProgramConfig, RecyclerConfig
from xsuite_recycler.rf_program import RFProgram


def test_publication_default_rf_program_endpoints():
    config = RFProgramConfig()
    program = RFProgram(config)

    assert program.voltage_53(0.0) == 80_000.0
    assert program.voltage_53(0.0025) == 40_000.0
    assert program.voltage_53(0.005) == 0.0
    assert program.voltage_53(0.100) == 0.0

    assert program.voltage_25(np.nextafter(0.005, 0.0)) == 0.0
    assert np.isclose(program.voltage_25(0.005), 3_000.0)
    assert np.isclose(program.voltage_25(0.090), 80_000.0)
    assert np.isclose(program.voltage_25(0.100), 80_000.0)


def test_legacy_keegan_code_starting_voltage_remains_explicitly_available():
    """The legacy Keegan-code 5 kV sensitivity remains selectable."""
    program = RFProgram(RFProgramConfig(rf25_initial_v=5_000.0))

    assert np.isclose(program.voltage_25(0.005), 5_000.0)
    assert np.isclose(program.voltage_25(0.090), 80_000.0)


def test_recycler_derived_frequencies_and_turn_count():
    machine = RecyclerConfig()
    ramp = RFProgramConfig()

    assert np.isclose(machine.revolution_frequency_hz, 89_809.74785452797)
    assert np.isclose(machine.rf53_frequency_hz, 52_808_131.73846245)
    assert np.isclose(machine.rf25_frequency_hz, 2_514_672.939926783)
    assert ramp.turns_to_end(machine) == 8083
