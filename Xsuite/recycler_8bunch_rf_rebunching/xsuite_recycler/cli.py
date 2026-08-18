"""Command-line interface for the Recycler rebunching simulation."""

from __future__ import annotations

import argparse
import json
from math import isfinite
from pathlib import Path
import sys

from .config import RFProgramConfig, RecyclerConfig
from .comparison_plotting import plot_result_comparison
from .distribution import (
    initial_distribution_from_generated,
    load_initial_distribution,
    zeta_to_dt,
)
from .final_dt_comparison import compare_final_dt, write_final_dt_comparison_plot
from .generator import (
    MicrobunchGeneratorConfig,
    generate_longitudinal_distribution,
)
from .input_analysis import compare_inputs
from .input_plotting import plot_input_comparison
from .plotting import plot_run_results
from .result_comparison import compare_results
from .screening import summarize_rf_screen
from .simulation import run_rf_only


def _ramp_from_args(args: argparse.Namespace) -> RFProgramConfig:
    """Build an RF program, converting the CLI's kilovolts to SI volts."""
    return RFProgramConfig(rf25_initial_v=args.rf25_start_kv * 1e3)


def _add_ramp_argument(parser: argparse.ArgumentParser) -> None:
    """Add the shared h=28 starting-voltage option to a parser."""
    parser.add_argument(
        "--rf25-start-kv",
        type=float,
        default=3.0,
        help=(
            "h=28 voltage at 5 ms in kV "
            "(default: 3, from the publication/ESME program; "
            "use 5 for the legacy Keegan-code sensitivity)"
        ),
    )


def _add_tracking_arguments(parser: argparse.ArgumentParser) -> None:
    """Add tracking duration, recording cadence, and plotting options."""
    parser.add_argument("--record-every", type=int, default=10)
    parser.add_argument(
        "--turns",
        type=int,
        default=None,
        help="override the 90 ms/8083-turn endpoint",
    )
    parser.add_argument(
        "--gamma-transition",
        type=float,
        default=RecyclerConfig().gamma_transition,
        help=(
            "Recycler transition gamma "
            f"(default: {RecyclerConfig().gamma_transition:.15g}, "
            "from Keegan's BLonD model; use 21.6 for Werkema's model)"
        ),
    )
    parser.add_argument("--no-figures", action="store_true")
    _add_ramp_argument(parser)


def _add_generator_arguments(
    parser: argparse.ArgumentParser,
    *,
    particles_per_microbunch: int | None,
    include_seed: bool = True,
) -> None:
    """Add options for the compact 168-microbunch input generator.

    ``particles_per_microbunch=None`` is used by the dual-mode ``run`` command
    so it can distinguish an omitted generator option from one explicitly
    supplied alongside a file input.
    """
    defaults = MicrobunchGeneratorConfig()
    generated = parser.add_argument_group("generated-input options")
    generated.add_argument(
        "--particles-per-microbunch",
        type=int,
        default=particles_per_microbunch,
        help=(
            "population in each of 168 generated microbunches "
            + (
                "(run default: 200)"
                if particles_per_microbunch is None
                else f"(default: {particles_per_microbunch})"
            )
        ),
    )
    if include_seed:
        generated.add_argument("--seed", type=int, default=defaults.seed)

    # None is an intentional sentinel in `run`; real generator defaults are
    # applied only after the command knows that generated input was selected.
    use_generator_defaults = particles_per_microbunch is not None

    def default(value: float) -> float | None:
        return value if use_generator_defaults else None

    generated.add_argument(
        "--time-origin-ns",
        type=float,
        default=default(defaults.time_origin_ns),
    )
    generated.add_argument(
        "--spacing-ns", type=float, default=default(defaults.spacing_ns)
    )
    generated.add_argument(
        "--sigma-time-ns", type=float, default=default(defaults.sigma_time_ns)
    )
    generated.add_argument(
        "--mean-dpop", type=float, default=default(defaults.mean_dpop)
    )
    generated.add_argument(
        "--sigma-dpop", type=float, default=default(defaults.sigma_dpop)
    )
    generated.add_argument(
        "--correlation", type=float, default=default(defaults.correlation)
    )


def _generator_config_from_args(
    args: argparse.Namespace,
    *,
    default_particles_per_microbunch: int,
) -> MicrobunchGeneratorConfig:
    """Merge CLI values with generator defaults for the selected workflow."""
    defaults = MicrobunchGeneratorConfig()

    def value(name: str):
        supplied = getattr(args, name)
        return getattr(defaults, name) if supplied is None else supplied

    return MicrobunchGeneratorConfig(
        particles_per_microbunch=(
            default_particles_per_microbunch
            if args.particles_per_microbunch is None
            else args.particles_per_microbunch
        ),
        time_origin_ns=value("time_origin_ns"),
        spacing_ns=value("spacing_ns"),
        sigma_time_ns=value("sigma_time_ns"),
        mean_dpop=value("mean_dpop"),
        sigma_dpop=value("sigma_dpop"),
        correlation=value("correlation"),
        seed=args.seed,
    )


_GENERATOR_OPTION_DESTS = (
    "particles_per_microbunch",
    "time_origin_ns",
    "spacing_ns",
    "sigma_time_ns",
    "mean_dpop",
    "sigma_dpop",
    "correlation",
)


def _supplied_generator_options(args: argparse.Namespace) -> list[str]:
    """Return CLI spellings of generator-only options explicitly supplied."""
    return [
        "--" + name.replace("_", "-")
        for name in _GENERATOR_OPTION_DESTS
        if getattr(args, name) is not None
    ]


def _write_json(path: Path, payload: dict) -> Path:
    """Write deterministic, standards-compliant JSON and return its path."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return path


def _static_assessment(report: dict) -> dict:
    """Apply a transparent first-screen envelope to the fitted generator."""
    distances = report["distribution_distances"]
    noise_distances = report["reference_split_half_noise"]["distribution_distances"]
    distance_keys = (
        "time_residual_wasserstein_standardized",
        "dpop_residual_wasserstein_standardized",
        "sliced_wasserstein_2d_standardized_mean",
    )
    ratios = {
        key: distances[key] / noise_distances[key]
        if noise_distances[key] > 0
        else None
        for key in distance_keys
    }
    geometry = report["center_geometry"]
    shape = report["pooled_shape"]
    criteria = {
        "absolute_origin_difference_le_0p2_ns": abs(geometry["origin_difference_ns"])
        <= 0.2,
        "absolute_spacing_difference_le_0p002_ns": abs(
            geometry["spacing_difference_ns"]
        )
        <= 0.002,
        "center_rmse_le_0p2_ns": geometry["center_rmse_ns"] <= 0.2,
        "absolute_time_sigma_difference_le_2_percent": abs(
            shape["time_sigma_relative_difference"]
        )
        <= 0.02,
        "absolute_dpop_sigma_difference_le_2_percent": abs(
            shape["dpop_sigma_relative_difference"]
        )
        <= 0.02,
        "absolute_correlation_difference_le_0p01": abs(
            shape["candidate_correlation"] - shape["reference_correlation"]
        )
        <= 0.01,
        "absolute_global_dpop_mean_difference_le_5e-7": abs(
            shape["global_dpop_mean_difference"]
        )
        <= 5e-7,
        "all_selected_shape_distances_le_2x_reference_split_control": all(
            value is not None and value <= 2.0 for value in ratios.values()
        ),
    }
    return {
        "name": "first static calibration screen",
        "passes_all_criteria": all(criteria.values()),
        "criteria": criteria,
        "candidate_to_reference_split_distance_ratios": ratios,
        "warning": (
            "The default parameters were inferred from this historical file. "
            "Passing is a calibration/adequacy check, not independent validation "
            "and not particle-by-particle reproduction."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    """Construct the command-line parser and all supported subcommands."""
    parser = argparse.ArgumentParser(
        description="RF-only Xsuite simulation of Recycler rebunching",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser(
        "run",
        help="track generated particles by default, or load particles with --input",
    )
    run.add_argument(
        "--input",
        type=Path,
        default=None,
        help="two-column z [m], dp/p file; omit to generate particles in memory",
    )
    run.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="destination (default: outputs/generated or outputs/from_file)",
    )
    run.add_argument(
        "--seed",
        type=int,
        default=MicrobunchGeneratorConfig().seed,
        help="generator RNG seed or file-sampling seed",
    )
    file_options = run.add_argument_group("file-input options")
    file_options.add_argument(
        "--max-particles",
        type=int,
        default=None,
        help="deterministic stratified sample (file default: 20000; 0 means all)",
    )
    file_options.add_argument(
        "--sampling-strata",
        type=int,
        default=None,
        help="ordered sampling strata (file default: 8; use 168 for microbunch balance)",
    )
    _add_generator_arguments(
        run,
        particles_per_microbunch=None,
        include_seed=False,
    )
    _add_tracking_arguments(run)

    plot = sub.add_parser("plot", help="regenerate standard plots from an HDF5 result")
    plot.add_argument("results", type=Path)
    plot.add_argument("--output-dir", type=Path, default=None)

    compare_inputs_parser = sub.add_parser(
        "compare-inputs",
        help="compare generated particles statistically with a reference input file",
    )
    compare_inputs_parser.add_argument("--input", type=Path, required=True)
    compare_inputs_parser.add_argument(
        "--output-dir", type=Path, default=Path("outputs/input_model")
    )
    _add_generator_arguments(compare_inputs_parser, particles_per_microbunch=6_400)

    compare_runs = sub.add_parser(
        "compare-runs",
        help="compare turn profiles and final bunch properties from two runs",
    )
    compare_runs.add_argument("reference", type=Path)
    compare_runs.add_argument("candidate", type=Path)
    compare_runs.add_argument(
        "--output-dir", type=Path, default=Path("outputs/result_comparison")
    )
    compare_runs.add_argument("--no-plot", action="store_true")
    compare_runs.add_argument(
        "--allow-partial-turn-overlap",
        action="store_true",
        help="permit exploratory comparisons with non-identical run/profile metadata",
    )

    compare_final_time = sub.add_parser(
        "compare-final-time",
        help="compare final arrival-time distributions and h=28 bucket occupancy",
    )
    compare_final_time.add_argument("reference", type=Path)
    compare_final_time.add_argument("candidate", type=Path)
    compare_final_time.add_argument(
        "--output-dir", type=Path, default=Path("outputs/final_dt_comparison")
    )
    compare_final_time.add_argument("--absolute-bins", type=int, default=2400)
    compare_final_time.add_argument("--absolute-min-ns", type=float, default=-1000.0)
    compare_final_time.add_argument("--absolute-max-ns", type=float, default=3800.0)
    compare_final_time.add_argument("--bucket-relative-bins", type=int, default=400)
    compare_final_time.add_argument(
        "--allow-gamma-transition-difference",
        action="store_true",
        help=(
            "controlled full-input study allowing only gamma_transition, "
            "momentum_compaction_factor, and slip_factor to differ"
        ),
    )

    screen = sub.add_parser(
        "summarize-screen",
        help="aggregate generated RF runs against historical resampling controls",
    )
    screen.add_argument("--reference", type=Path, required=True)
    screen.add_argument(
        "--generated", type=Path, action="append", required=True, dest="generated"
    )
    screen.add_argument(
        "--control", type=Path, action="append", required=True, dest="controls"
    )
    screen.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/rf_screen_33k6/rf_screen_summary.json"),
    )

    return parser


def main(argv: list[str] | None = None) -> None:
    """Parse command-line arguments and execute the requested workflow."""
    parser = build_parser()
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(raw_argv)
    # A programmatic caller does not have a meaningful sys.argv[0], so use the
    # repository's documented entry point in that case. The Python executable
    # itself is recorded separately with the runtime provenance.
    invocation_argv = [sys.argv[0] if argv is None else "main.py", *raw_argv]

    if args.command == "run":
        if args.seed < 0 or args.seed >= 2**64:
            parser.error("--seed must be in the range [0, 2**64)")
        if args.record_every <= 0:
            parser.error("--record-every must be positive")
        if args.turns is not None and args.turns <= 0:
            parser.error("--turns must be positive")
        if not isfinite(args.gamma_transition) or args.gamma_transition <= 1.0:
            parser.error("--gamma-transition must be finite and greater than 1")

        # Changing gamma_t changes alpha_c and the slip factor while leaving
        # momentum, circumference, RF frequencies, and input coordinates fixed.
        # This keeps the Werkema comparison to one intentional machine change.
        machine = RecyclerConfig(gamma_transition=args.gamma_transition)
        # `run` has two mutually exclusive input paths. Keeping their unused
        # options at None lets us report accidental cross-mode combinations.
        if args.input is None:
            incompatible = [
                option
                for option, value in (
                    ("--max-particles", args.max_particles),
                    ("--sampling-strata", args.sampling_strata),
                )
                if value is not None
            ]
            if incompatible:
                parser.error(
                    f"{', '.join(incompatible)} require --input; "
                    "omit them when generating particles"
                )
            try:
                generator_config = _generator_config_from_args(
                    args,
                    default_particles_per_microbunch=200,
                )
            except (TypeError, ValueError) as exc:
                parser.error(str(exc))
            generated = generate_longitudinal_distribution(generator_config)
            # The generator works in arrival time and dp/p. Convert those
            # familiar beam coordinates to Xsuite's zeta and ptau here.
            distribution = initial_distribution_from_generated(generated, machine)
            output_dir = args.output_dir or Path("outputs/generated")
            input_path = None
            max_particles = None
            sampling_strata = None
            print(
                "Input source: generated "
                f"(168 x {generator_config.particles_per_microbunch:,} = "
                f"{distribution.n_particles:,} particles)"
            )
        else:
            incompatible = _supplied_generator_options(args)
            if incompatible:
                parser.error(
                    f"{', '.join(incompatible)} cannot be used with --input"
                )
            max_particles = 20_000 if args.max_particles is None else args.max_particles
            sampling_strata = (
                8 if args.sampling_strata is None else args.sampling_strata
            )
            if max_particles < 0:
                parser.error("--max-particles must be nonnegative; use 0 for all rows")
            if sampling_strata <= 0 or sampling_strata % machine.n_bunches:
                parser.error("--sampling-strata must be a positive multiple of 8")
            if 0 < max_particles < sampling_strata:
                parser.error(
                    "--max-particles must be 0 or at least --sampling-strata"
                )
            distribution = None
            output_dir = args.output_dir or Path("outputs/from_file")
            input_path = args.input
            print(f"Input source: file {args.input}")

        output_dir.mkdir(parents=True, exist_ok=True)
        results = run_rf_only(
            input_path,
            output_dir / "rf_only_results.h5",
            machine=machine,
            ramp=_ramp_from_args(args),
            max_particles=max_particles,
            record_every=args.record_every,
            n_turns=args.turns,
            seed=args.seed,
            sampling_strata=sampling_strata,
            initial_distribution=distribution,
            invocation_argv=invocation_argv,
        )
        if not args.no_figures:
            for path in plot_run_results(results, output_dir):
                print(path)
        return

    if args.command == "plot":
        output_dir = args.output_dir or args.results.parent
        for path in plot_run_results(args.results, output_dir):
            print(path)
        return

    if args.command == "compare-inputs":
        try:
            generator_config = _generator_config_from_args(
                args,
                default_particles_per_microbunch=6_400,
            )
        except (TypeError, ValueError) as exc:
            parser.error(str(exc))
        args.output_dir.mkdir(parents=True, exist_ok=True)
        machine = RecyclerConfig()
        historical = load_initial_distribution(args.input, machine, max_particles=0)
        reference_dt_s = zeta_to_dt(historical.zeta_m, machine.beta0)
        # Undo the code-faithful dp/p -> ptau conversion used by the loader so
        # both inputs are compared in the generator's original coordinates.
        reference_dpop = historical.ptau * machine.blond_energy_beta
        generated = generate_longitudinal_distribution(generator_config)
        report = compare_inputs(
            reference_dt_s,
            reference_dpop,
            generated.dt_s,
            generated.dpop,
        )
        report["historical_source"] = {
            "path": historical.source_path,
            "sha256": historical.source_sha256,
        }
        report["generator"] = generated.metadata
        report["assessment"] = _static_assessment(report)
        report_path = _write_json(
            args.output_dir / "input_model_comparison.json", report
        )
        plot_path = plot_input_comparison(
            reference_dt_s,
            reference_dpop,
            generated.dt_s,
            generated.dpop,
            args.output_dir / "input_model_diagnostics.png",
        )
        print(report_path)
        print(plot_path)
        print(
            "static screen:",
            "PASS" if report["assessment"]["passes_all_criteria"] else "CHECK",
        )
        return

    if args.command == "compare-runs":
        args.output_dir.mkdir(parents=True, exist_ok=True)
        require_matching_run = not args.allow_partial_turn_overlap
        report = compare_results(
            args.reference,
            args.candidate,
            require_matching_run=require_matching_run,
        )
        report_path = _write_json(
            args.output_dir / "rf_result_comparison.json", report
        )
        print(report_path)
        if not args.no_plot:
            print(
                plot_result_comparison(
                    args.reference,
                    args.candidate,
                    args.output_dir / "rf_result_comparison.png",
                    require_matching_run=require_matching_run,
                )
            )
        summary = report["profiles"]
        print(
            "median profile W1:",
            f"{summary['wasserstein_summary_ns']['median']:.4g} ns;",
            "median JS:",
            f"{summary['jensen_shannon_summary']['median']:.4g}",
        )
        return

    if args.command == "compare-final-time":
        args.output_dir.mkdir(parents=True, exist_ok=True)
        options = {
            "absolute_bins": args.absolute_bins,
            "absolute_range_ns": (
                args.absolute_min_ns,
                args.absolute_max_ns,
            ),
            "bucket_relative_bins": args.bucket_relative_bins,
            "allow_gamma_transition_difference": (
                args.allow_gamma_transition_difference
            ),
        }
        report = compare_final_dt(args.reference, args.candidate, **options)
        report_path = _write_json(
            args.output_dir / "final_dt_comparison.json", report
        )
        plot_path = write_final_dt_comparison_plot(
            args.reference,
            args.candidate,
            args.output_dir / "final_dt_comparison.png",
            reference_label=(
                "Keegan γₜ"
                if args.allow_gamma_transition_difference
                else "Historical input"
            ),
            candidate_label=(
                "Werkema γₜ"
                if args.allow_gamma_transition_difference
                else "Generated input"
            ),
            **options,
        )
        print(report_path)
        print(plot_path)
        comparison = report["comparison"]
        print(
            "final absolute-dt empirical W1:",
            f"{comparison['absolute_empirical_wasserstein_1_ns']:.6g} ns;",
            "absolute-dt histogram JS:",
            f"{comparison['absolute_histogram']['jensen_shannon_distance']:.6g}",
        )
        reference_local = report["reference"][
            "nearest_h28_bucket_folded_125ns_local_window"
        ]["global"]
        candidate_local = report["candidate"][
            "nearest_h28_bucket_folded_125ns_local_window"
        ]["global"]
        print(
            "nearest-bucket folded |dt_local| > 125 ns:",
            f"reference {reference_local['raw_outside_count']:,} "
            f"({reference_local['raw_outside_fraction']:.6%});",
            f"candidate {candidate_local['raw_outside_count']:,} "
            f"({candidate_local['raw_outside_fraction']:.6%})",
        )
        reference_window = report["reference"][
            "unwrapped_initial_lineage_125ns_window"
        ]["global"]
        candidate_window = report["candidate"][
            "unwrapped_initial_lineage_125ns_window"
        ]["global"]
        print(
            "unwrapped lineage |dt_rel| > 125 ns (tail or migration):",
            f"reference {reference_window['raw_outside_count']:,} "
            f"({reference_window['raw_outside_fraction']:.6%});",
            f"candidate {candidate_window['raw_outside_count']:,} "
            f"({candidate_window['raw_outside_fraction']:.6%})",
        )
        return

    if args.command == "summarize-screen":
        report = summarize_rf_screen(
            args.reference,
            args.generated,
            args.controls,
        )
        print(_write_json(args.output, report))
        print(
            "dynamic screen:",
            "PASS" if report["assessment"]["passes_all_criteria"] else "CHECK",
        )
        return

    parser.error(f"Unknown command: {args.command}")
