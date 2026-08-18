from pathlib import Path

import h5py

from xsuite_recycler.config import RFProgramConfig, RecyclerConfig
import xsuite_recycler.plotting as plotting


def test_standard_plot_filenames_are_domain_named(monkeypatch, tmp_path):
    results = tmp_path / "results.h5"
    with h5py.File(results, "w"):
        pass

    monkeypatch.setattr(
        plotting,
        "_machine_and_ramp_from_h5",
        lambda h5: (RecyclerConfig(), RFProgramConfig()),
    )

    written = []

    def fake_rf(output_path, ramp):
        written.append(Path(output_path))
        return Path(output_path)

    def fake_result_plot(results_path, output_path):
        assert Path(results_path) == results
        written.append(Path(output_path))
        return Path(output_path)

    monkeypatch.setattr(plotting, "plot_rf_voltage_program", fake_rf)
    monkeypatch.setattr(plotting, "plot_final_phase_space", fake_result_plot)
    monkeypatch.setattr(plotting, "plot_rebunching_waterfall", fake_result_plot)

    output_dir = tmp_path / "plots"
    paths = plotting.plot_run_results(results, output_dir)

    assert paths == written == [
        output_dir / "rf_voltage_program.png",
        output_dir / "final_phase_space.png",
        output_dir / "rebunching_waterfall.png",
    ]


def test_phase_space_title_uses_selected_structural_bunch(monkeypatch, tmp_path):
    results = tmp_path / "results.h5"
    with h5py.File(results, "w") as h5:
        h5.attrs["n_turns"] = 10
        final = h5.create_group("final_particles")
        final.create_dataset("zeta_m", data=[0.0])
        final.create_dataset("ptau", data=[0.0])
        final.create_dataset("state", data=[1])
        final.create_dataset("bunch_index", data=[3])

    monkeypatch.setattr(
        plotting,
        "_machine_and_ramp_from_h5",
        lambda h5: (RecyclerConfig(), RFProgramConfig()),
    )
    captured = {}

    def capture_title(fig, output_path):
        captured["title"] = fig.axes[0].get_title()
        plotting.plt.close(fig)
        return Path(output_path)

    monkeypatch.setattr(plotting, "_save", capture_title)
    plotting.plot_final_phase_space(
        results,
        tmp_path / "phase_space.png",
        bunch_index=3,
    )

    assert "bunch 4" in captured["title"]
