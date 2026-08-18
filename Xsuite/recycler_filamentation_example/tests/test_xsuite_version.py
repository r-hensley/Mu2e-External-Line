"""Focused regression tests for the Xsuite Recycler teaching implementation."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from xsuite_version.run import build_parser
from xsuite_version.simulation import (
    C_LIGHT_M_PER_S,
    RecyclerInputs,
    build_line,
    default_snapshot_times,
    derive_from_xsuite,
    make_particles,
    particles_to_local_coordinates,
    rf_stage,
    rf_voltages,
    track_rebunching,
    validate_one_turn,
)


@pytest.fixture(scope="module")
def xsuite_model():
    """Build the shared line once because Xsuite tracker construction is costly."""

    inputs = RecyclerInputs()
    line = build_line(inputs)
    derived = derive_from_xsuite(line, inputs)
    return inputs, line, derived


@pytest.fixture()
def clean_xsuite_model(xsuite_model):
    """Reset time-dependent cavity strengths before every line-using test."""

    inputs, line, derived = xsuite_model
    line["rf53"].voltage = 0.0
    line["rf2p5"].voltage = 0.0
    return inputs, line, derived


def test_default_reference_tracks_current_directory_layout() -> None:
    """The default comparison file should live under ``manual_version``."""

    args = build_parser().parse_args([])
    project_root = Path(__file__).resolve().parents[1]
    expected = project_root / "manual_version" / "outputs" / "nominal_result.npz"
    assert args.reference_result == expected
    assert args.reference_result.is_file()


def test_line_order_and_xsuite_derived_quantities(xsuite_model) -> None:
    """Check the intended drift-then-kick line and Xsuite reference results."""

    inputs, line, derived = xsuite_model
    assert tuple(line.element_names) == ("recycler_slip", "rf53", "rf2p5")
    assert derived.line_length_m == pytest.approx(inputs.circumference_m)
    assert derived.energy0_eV == pytest.approx(inputs.reference_energy_eV)
    assert derived.t_rev0_s == pytest.approx(
        inputs.circumference_m / (derived.beta0 * C_LIGHT_M_PER_S)
    )
    assert derived.rf53_frequency_hz == pytest.approx(
        inputs.harmonic_53 / derived.t_rev0_s
    )
    assert derived.rf25_frequency_hz == pytest.approx(
        inputs.harmonic_25 / derived.t_rev0_s
    )
    assert derived.slip_factor < 0.0
    assert derived.synchrotron_period_80kv_s == pytest.approx(18.4e-3, abs=0.1e-3)


def test_rf_program_boundaries_and_monotonic_ramp(xsuite_model) -> None:
    """Exercise the porch, turn-off, nonlinear ramp, and flattop boundaries."""

    inputs, _, derived = xsuite_model
    porch_end = inputs.porch_turns * derived.t_rev0_s
    midpoint = 0.5 * (derived.time_53_off_s + derived.time_25_flat_s)
    times = np.asarray(
        [0.0, porch_end, derived.time_53_off_s, midpoint, derived.time_25_flat_s]
    )
    voltage_53, voltage_25 = rf_voltages(times, inputs, derived)

    np.testing.assert_allclose(
        voltage_53,
        [inputs.voltage_53_eV, inputs.voltage_53_eV, 0.0, 0.0, 0.0],
    )
    assert voltage_25[0] == 0.0
    assert voltage_25[1] == 0.0
    assert voltage_25[2] == pytest.approx(inputs.voltage_25_initial_eV)
    assert inputs.voltage_25_initial_eV < voltage_25[3] < inputs.voltage_25_final_eV
    assert voltage_25[4] == pytest.approx(inputs.voltage_25_final_eV)

    ramp_times = np.linspace(derived.time_53_off_s, derived.time_25_flat_s, 101)
    _, ramp_voltage = rf_voltages(ramp_times, inputs, derived)
    assert np.all(np.diff(ramp_voltage) > 0.0)
    assert rf_stage(0.0, inputs, derived) == "one-turn 53 MHz porch"
    assert rf_stage(porch_end, inputs, derived) == "53 MHz turn-off"
    assert rf_stage(derived.time_53_off_s, inputs, derived) == "2.5 MHz iso-adiabatic ramp"
    assert rf_stage(derived.time_25_flat_s, inputs, derived) == "2.5 MHz flattop"


def test_coordinate_round_trip(clean_xsuite_model) -> None:
    """Manual phase/energy coordinates should round-trip through Xsuite exactly."""

    inputs, line, derived = clean_xsuite_model
    theta = np.asarray([-0.08, 0.0, 0.07])
    energy = np.asarray([-2.5e6, 0.0, 3.0e6])
    particles = make_particles(line, theta, energy, inputs, derived)
    recovered_theta, recovered_energy = particles_to_local_coordinates(
        particles, inputs, derived
    )

    np.testing.assert_allclose(recovered_theta, theta, rtol=0.0, atol=2e-15)
    np.testing.assert_allclose(recovered_energy, energy, rtol=2e-14, atol=1e-7)


def test_one_turn_map_matches_independent_expectation(clean_xsuite_model) -> None:
    """Validate coordinate signs, slip placement, and the two cavity kicks."""

    inputs, line, derived = clean_xsuite_model
    error = validate_one_turn(line, inputs, derived)
    assert abs(error["theta_error_rad"]) < 2e-12
    assert abs(error["energy_error_eV"]) < 1e-4


def test_snapshot_schedule_contains_all_named_events(xsuite_model) -> None:
    """The stored schedule must include RF transitions and every observation."""

    inputs, _, derived = xsuite_model
    times = default_snapshot_times(inputs, derived)
    assert times[0] == 0.0
    assert np.all(np.diff(times) > 0.0)
    assert np.any(np.isclose(times, derived.time_53_off_s, rtol=0.0, atol=1e-15))
    assert np.any(np.isclose(times, derived.time_25_flat_s, rtol=0.0, atol=1e-15))
    for extraction_time in derived.extraction_times_s:
        assert np.any(np.isclose(times, extraction_time, rtol=0.0, atol=1e-15))


def test_short_multiturn_tracking_smoke(clean_xsuite_model) -> None:
    """Track a few turns to exercise snapshots, diagnostics, and survival checks."""

    inputs, line, derived = clean_xsuite_model
    # Replacing only the final observation time shortens this test without
    # changing the production simulation's default eight-observation history.
    short_derived = replace(
        derived,
        extraction_times_s=np.asarray([5.0 * derived.t_rev0_s]),
    )
    theta = np.asarray([-0.02, 0.0, 0.02])
    energy = np.asarray([-1.0e6, 0.0, 1.0e6])
    source_id = np.asarray([-1, 0, 1])

    result = track_rebunching(
        line,
        theta,
        energy,
        source_id,
        inputs,
        short_derived,
        diagnostic_stride_turns=2,
        display_per_source=1,
        show_progress=False,
    )

    assert result.snapshot_turns[0] == 0
    assert result.snapshot_turns[-1] == 5
    np.testing.assert_array_equal(result.diagnostic_turns, [0, 2, 4, 5])
    np.testing.assert_allclose(result.theta_rad[0], theta, rtol=0.0, atol=2e-9)
    np.testing.assert_allclose(result.delta_energy_eV[0], energy, rtol=0.0, atol=0.2)
    assert np.all(result.final_state > 0)
    assert np.isfinite(result.theta_rad).all()
    assert np.isfinite(result.delta_energy_eV).all()
