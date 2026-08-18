import json

import h5py
import numpy as np

from xsuite_recycler.config import C_LIGHT_M_PER_S
from xsuite_recycler.screening import summarize_rf_screen


def _write_identical_screen_result(path, source_kind):
    beta0 = 0.9
    labels = np.repeat(np.arange(8), 2)
    time_s = np.tile(np.asarray([-10e-9, 10e-9]), 8)
    with h5py.File(path, "w") as h5:
        h5.attrs.update(
            {
                "format_version": 1,
                "model": "RF-only test",
                "xtrack_version": "test",
                "source_blond_commit": "abc",
                "source_blond_driver": "driver.py",
                "coordinate_convention": "test-zeta",
                "energy_convention": "test-ptau",
                "n_turns": 10,
                "end_state_time_s": 1e-4,
                "last_kick_time_s": 9e-5,
                "first_rf25_kick_turn": 5,
                "record_every": 10,
                "sample_particles": 16,
                "source_rows": 16,
                "macro_equivalent_weight": 1.0,
                "physical_proton_weight": 1.0,
                "last_kick_rf53_v": 0.0,
                "last_kick_rf25_v": 80_000.0,
                "input_source_kind": source_kind,
            }
        )
        machine = h5.create_group("machine")
        machine.attrs["beta0"] = beta0
        machine.attrs["p0c_ev"] = 8e9
        machine.attrs["rf25_frequency_hz"] = 2.5e6
        ramp = h5.create_group("rf_program")
        ramp.attrs["rf25_initial_v"] = 3000.0
        profiles = h5.create_group("profiles")
        profiles.create_dataset("turn", data=[0, 10])
        profiles.create_dataset("time_edges_ns", data=[-100.0, 0.0, 100.0])
        profiles.create_dataset("counts_macro_equivalent", data=[[8, 8], [8, 8]])
        profiles.create_dataset("outside_macro_equivalent", data=[0.0, 0.0])
        final = h5.create_group("final_particles")
        final.create_dataset("zeta_m", data=-beta0 * C_LIGHT_M_PER_S * time_s)
        final.create_dataset("ptau", data=np.tile([-1e-3, 1e-3], 8))
        final.create_dataset("state", data=np.ones(16, dtype=int))
        final.create_dataset("bunch_index", data=labels)


def test_identical_zero_noise_screen_is_json_safe_and_passes(tmp_path):
    reference = tmp_path / "reference.h5"
    generated = tmp_path / "generated.h5"
    control = tmp_path / "control.h5"
    _write_identical_screen_result(reference, "file")
    _write_identical_screen_result(generated, "generated")
    _write_identical_screen_result(control, "file")

    report = summarize_rf_screen(reference, [generated], [control])

    assert report["assessment"]["passes_all_criteria"]
    assert report["assessment"]["name"] == "16-particle RF dynamic screen"
    assert (
        report["generated_pairs"][0]["ratios_to_historical_resampling_baseline"][
            "median_profile_wasserstein"
        ]
        == 0.0
    )
    json.dumps(report, allow_nan=False)
