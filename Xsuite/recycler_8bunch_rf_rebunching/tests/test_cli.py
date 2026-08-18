from pathlib import Path
from types import SimpleNamespace

import pytest

import xsuite_recycler.cli as cli
from xsuite_recycler.cli import _static_assessment


def _otherwise_matching_report(global_dpop_mean_difference):
    distances = {
        "time_residual_wasserstein_standardized": 1.0,
        "dpop_residual_wasserstein_standardized": 1.0,
        "sliced_wasserstein_2d_standardized_mean": 1.0,
    }
    return {
        "distribution_distances": distances,
        "reference_split_half_noise": {
            "distribution_distances": dict(distances),
        },
        "center_geometry": {
            "origin_difference_ns": 0.0,
            "spacing_difference_ns": 0.0,
            "center_rmse_ns": 0.0,
        },
        "pooled_shape": {
            "time_sigma_relative_difference": 0.0,
            "dpop_sigma_relative_difference": 0.0,
            "candidate_correlation": 0.0,
            "reference_correlation": 0.0,
            "global_dpop_mean_difference": global_dpop_mean_difference,
        },
    }


def test_static_screen_rejects_large_global_momentum_shift():
    shifted = _static_assessment(_otherwise_matching_report(0.01))
    matched = _static_assessment(_otherwise_matching_report(1e-7))

    assert not shifted["passes_all_criteria"]
    assert not shifted["criteria"][
        "absolute_global_dpop_mean_difference_le_5e-7"
    ]
    assert matched["passes_all_criteria"]


def _patch_run_pipeline(monkeypatch):
    captured = {}
    generated = object()
    distribution = SimpleNamespace(n_particles=33_600)

    def fake_generate(config):
        captured["generator_config"] = config
        return generated

    def fake_convert(value, machine):
        assert value is generated
        captured["machine"] = machine
        return distribution

    def fake_run(*args, **kwargs):
        captured["run_args"] = args
        captured["run_kwargs"] = kwargs
        return Path(args[1])

    monkeypatch.setattr(cli, "generate_longitudinal_distribution", fake_generate)
    monkeypatch.setattr(cli, "initial_distribution_from_generated", fake_convert)
    monkeypatch.setattr(cli, "run_rf_only", fake_run)
    return captured, distribution


def test_top_level_command_list_is_small_and_domain_named():
    parser = cli.build_parser()
    subparsers = next(
        action
        for action in parser._actions
        if isinstance(action, cli.argparse._SubParsersAction)
    )

    assert list(subparsers.choices) == [
        "run",
        "plot",
        "compare-inputs",
        "compare-runs",
        "compare-final-time",
        "summarize-screen",
    ]


def test_compare_final_time_forwards_controlled_gamma_option(
    monkeypatch, tmp_path, capsys
):
    """The CLI exposes the narrow gamma_t compatibility mode to both outputs."""
    captured = {}
    report = {
        "comparison": {
            "absolute_empirical_wasserstein_1_ns": 1.25,
            "absolute_histogram": {"jensen_shannon_distance": 0.02},
        },
        "reference": {
            "nearest_h28_bucket_folded_125ns_local_window": {
                "global": {"raw_outside_count": 1, "raw_outside_fraction": 0.01}
            },
            "unwrapped_initial_lineage_125ns_window": {
                "global": {"raw_outside_count": 10, "raw_outside_fraction": 0.1}
            }
        },
        "candidate": {
            "nearest_h28_bucket_folded_125ns_local_window": {
                "global": {"raw_outside_count": 2, "raw_outside_fraction": 0.02}
            },
            "unwrapped_initial_lineage_125ns_window": {
                "global": {"raw_outside_count": 12, "raw_outside_fraction": 0.12}
            }
        },
    }

    def fake_compare(reference, candidate, **options):
        captured["compare"] = (reference, candidate, options)
        return report

    def fake_plot(reference, candidate, output, **options):
        captured["plot"] = (reference, candidate, output, options)
        return output

    monkeypatch.setattr(cli, "compare_final_dt", fake_compare)
    monkeypatch.setattr(cli, "write_final_dt_comparison_plot", fake_plot)
    output_dir = tmp_path / "comparison"
    cli.main(
        [
            "compare-final-time",
            "keegan.h5",
            "werkema.h5",
            "--allow-gamma-transition-difference",
            "--output-dir",
            str(output_dir),
        ]
    )

    assert captured["compare"][2]["allow_gamma_transition_difference"] is True
    plot_options = captured["plot"][3]
    assert plot_options["allow_gamma_transition_difference"] is True
    assert plot_options["reference_label"] == "Keegan γₜ"
    assert plot_options["candidate_label"] == "Werkema γₜ"
    output = capsys.readouterr().out
    assert "nearest-bucket folded |dt_local| > 125 ns" in output
    assert "unwrapped lineage |dt_rel| > 125 ns" in output


def test_bare_run_generates_default_distribution(monkeypatch, tmp_path, capsys):
    captured, distribution = _patch_run_pipeline(monkeypatch)
    output_dir = tmp_path / "generated"

    cli.main(["run", "--no-figures", "--output-dir", str(output_dir)])

    config = captured["generator_config"]
    assert config.particles_per_microbunch == 200
    assert config.seed == 202208
    assert captured["run_args"] == (
        None,
        output_dir / "rf_only_results.h5",
    )
    assert captured["run_kwargs"]["initial_distribution"] is distribution
    assert captured["run_kwargs"]["max_particles"] is None
    assert captured["run_kwargs"]["sampling_strata"] is None
    assert captured["run_kwargs"]["invocation_argv"] == [
        "main.py",
        "run",
        "--no-figures",
        "--output-dir",
        str(output_dir),
    ]
    assert captured["run_kwargs"]["ramp"].rf25_initial_v == 3_000.0
    assert captured["machine"].gamma_transition == pytest.approx(
        20.257643837730637
    )
    assert "33,600 particles" in capsys.readouterr().out


def test_run_accepts_legacy_keegan_code_5kv_sensitivity(monkeypatch, tmp_path):
    """An explicit CLI option selects the legacy Keegan-code 5 kV start."""
    captured, _ = _patch_run_pipeline(monkeypatch)

    cli.main(
        [
            "run",
            "--rf25-start-kv",
            "5",
            "--no-figures",
            "--output-dir",
            str(tmp_path / "keegan_5kv_sensitivity"),
        ]
    )

    assert captured["run_kwargs"]["ramp"].rf25_initial_v == 5_000.0


def test_run_accepts_werkema_transition_gamma(monkeypatch, tmp_path):
    """The gamma_t study changes slip factor without changing the RF ramp."""
    captured, _ = _patch_run_pipeline(monkeypatch)

    cli.main(
        [
            "run",
            "--gamma-transition",
            "21.6",
            "--no-figures",
            "--output-dir",
            str(tmp_path / "werkema_gamma"),
        ]
    )

    machine = captured["machine"]
    assert machine.gamma_transition == pytest.approx(21.6)
    assert machine.slip_factor == pytest.approx(-0.00900839413100925)
    assert captured["run_kwargs"]["ramp"].rf25_initial_v == 3_000.0


def test_run_with_input_uses_file_defaults(monkeypatch, tmp_path, capsys):
    captured, _ = _patch_run_pipeline(monkeypatch)
    output_dir = tmp_path / "from_file"
    input_path = tmp_path / "input.txt"

    cli.main(
        [
            "run",
            "--input",
            str(input_path),
            "--no-figures",
            "--output-dir",
            str(output_dir),
        ]
    )

    assert "generator_config" not in captured
    assert captured["run_args"] == (
        input_path,
        output_dir / "rf_only_results.h5",
    )
    assert captured["run_kwargs"]["initial_distribution"] is None
    assert captured["run_kwargs"]["max_particles"] == 20_000
    assert captured["run_kwargs"]["sampling_strata"] == 8
    assert f"Input source: file {input_path}" in capsys.readouterr().out


def test_run_forwards_source_specific_overrides(monkeypatch, tmp_path):
    captured, _ = _patch_run_pipeline(monkeypatch)
    cli.main(
        [
            "run",
            "--particles-per-microbunch",
            "17",
            "--spacing-ns",
            "20.5",
            "--seed",
            "9",
            "--no-figures",
            "--output-dir",
            str(tmp_path / "generated"),
        ]
    )
    assert captured["generator_config"].particles_per_microbunch == 17
    assert captured["generator_config"].spacing_ns == 20.5
    assert captured["generator_config"].seed == 9

    captured, _ = _patch_run_pipeline(monkeypatch)
    cli.main(
        [
            "run",
            "--input",
            str(tmp_path / "input.txt"),
            "--max-particles",
            "0",
            "--sampling-strata",
            "168",
            "--seed",
            "9",
            "--no-figures",
            "--output-dir",
            str(tmp_path / "file"),
        ]
    )
    assert captured["run_kwargs"]["max_particles"] == 0
    assert captured["run_kwargs"]["sampling_strata"] == 168
    assert captured["run_kwargs"]["seed"] == 9


def test_run_uses_source_specific_output_defaults(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    captured, _ = _patch_run_pipeline(monkeypatch)
    cli.main(["run", "--no-figures"])
    assert captured["run_args"][1] == Path("outputs/generated/rf_only_results.h5")

    captured, _ = _patch_run_pipeline(monkeypatch)
    cli.main(["run", "--input", "input.txt", "--no-figures"])
    assert captured["run_args"][1] == Path("outputs/from_file/rf_only_results.h5")


@pytest.mark.parametrize(
    ("option", "value"),
    [("--max-particles", "100"), ("--sampling-strata", "8")],
)
def test_generated_run_rejects_file_only_options(option, value, tmp_path):
    output_dir = tmp_path / "should_not_exist"
    with pytest.raises(SystemExit, match="2"):
        cli.main(["run", option, value, "--output-dir", str(output_dir)])
    assert not output_dir.exists()


@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("--particles-per-microbunch", "10"),
        ("--time-origin-ns", "0"),
        ("--spacing-ns", "20"),
        ("--sigma-time-ns", "5"),
        ("--mean-dpop", "0"),
        ("--sigma-dpop", "0.001"),
        ("--correlation", "0.1"),
    ],
)
def test_file_run_rejects_generator_only_options(option, value, tmp_path):
    output_dir = tmp_path / "should_not_exist"
    with pytest.raises(SystemExit, match="2"):
        cli.main(
            [
                "run",
                "--input",
                str(tmp_path / "input.txt"),
                option,
                value,
                "--output-dir",
                str(output_dir),
            ]
        )
    assert not output_dir.exists()


@pytest.mark.parametrize(
    "arguments",
    [
        ["--particles-per-microbunch", "0"],
        ["--seed", "-1"],
        ["--input", "input.txt", "--seed", "-1"],
        ["--record-every", "0"],
        ["--turns", "0"],
        ["--gamma-transition", "1"],
        ["--gamma-transition", "nan"],
        ["--input", "input.txt", "--max-particles", "-1"],
        ["--input", "input.txt", "--sampling-strata", "7"],
        ["--input", "input.txt", "--max-particles", "1"],
        [
            "--input",
            "input.txt",
            "--max-particles",
            "20",
            "--sampling-strata",
            "168",
        ],
    ],
)
def test_run_rejects_invalid_source_parameters(arguments, tmp_path):
    output_dir = tmp_path / "should_not_exist"
    with pytest.raises(SystemExit, match="2"):
        cli.main(["run", *arguments, "--output-dir", str(output_dir)])
    assert not output_dir.exists()


def test_compare_inputs_rejects_invalid_generator_before_output(tmp_path):
    output_dir = tmp_path / "should_not_exist"
    with pytest.raises(SystemExit, match="2"):
        cli.main(
            [
                "compare-inputs",
                "--input",
                str(tmp_path / "input.txt"),
                "--particles-per-microbunch",
                "0",
                "--output-dir",
                str(output_dir),
            ]
        )
    assert not output_dir.exists()


@pytest.mark.parametrize(
    "command",
    ["figure5", "run-generated", "analyze-input", "compare-results", "compare-final-dt"],
)
def test_old_commands_are_removed(command):
    with pytest.raises(SystemExit, match="2"):
        cli.build_parser().parse_args([command])
