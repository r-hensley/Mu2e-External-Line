#!/usr/bin/env python3
"""Collect the completed 5 kV convergence study into one report and figure.

Run this after ``analyze_convergence_5kv.sh``.  The script deliberately knows
the expected study matrix: missing HDF5 runs or analysis reports are treated as
errors instead of silently producing a partial convergence summary.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


PROJECT_DIR = Path(__file__).resolve().parent
OUTPUT_ROOT = PROJECT_DIR / "outputs" / "convergence_5kv"
ANALYSIS_ROOT = OUTPUT_ROOT / "analysis"
PARTICLE_COUNTS = (33_600, 100_800, 336_000, 1_075_200)
REDUCED_COUNTS = PARTICLE_COUNTS[:-1]
SEEDS = (202_208, 202_209, 202_210, 202_211, 202_212)
SOURCE_LABELS = ("historical", "generated")
HISTORICAL_INPUT_SHA256 = (
    "3b591b8381cb7cf57a88d7db08d535f2de182720d9fcb3c335800d7b3283e7e1"
)


def _relative(path: Path) -> str:
    """Return a portable project-relative path for the JSON report."""

    return path.relative_to(PROJECT_DIR).as_posix()


def _load_json(path: Path) -> dict[str, Any]:
    """Load one required JSON object and report a clear path-level error."""

    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty analysis report: {_relative(path)}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in {_relative(path)}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {_relative(path)}")
    return value


def _required(mapping: dict[str, Any], *keys: str) -> Any:
    """Read a nested report value while naming any missing key chain."""

    value: Any = mapping
    traversed: list[str] = []
    for key in keys:
        traversed.append(key)
        if not isinstance(value, dict) or key not in value:
            raise KeyError(f"Missing required report field: {'.'.join(traversed)}")
        value = value[key]
    return value


def _result_path(source_label: str, particle_count: int, seed: int) -> Path:
    """Resolve one HDF5 result in the fixed convergence-matrix layout."""

    directory = f"{source_label}_n{particle_count}_seed{seed}"
    return OUTPUT_ROOT / directory / "rf_only_results.h5"


def _expected_runs() -> list[tuple[str, int, int, Path]]:
    """List the eight count-ladder and eight additional seed-control runs."""

    runs = [
        (source, particle_count, 202_208, _result_path(source, particle_count, 202_208))
        for particle_count in PARTICLE_COUNTS
        for source in SOURCE_LABELS
    ]
    runs.extend(
        (source, 100_800, seed, _result_path(source, 100_800, seed))
        for seed in SEEDS[1:]
        for source in SOURCE_LABELS
    )
    if len(runs) != 16:
        raise RuntimeError(f"Internal error: expected 16 runs, constructed {len(runs)}")
    return runs


def _read_run_metadata(
    source_label: str,
    particle_count: int,
    seed: int,
    path: Path,
) -> dict[str, Any]:
    """Validate one HDF5 run and retain its runtime and provenance metadata."""

    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty HDF5 result: {_relative(path)}")

    expected_source_kind = "file" if source_label == "historical" else "generated"
    try:
        with h5py.File(path, "r") as h5:
            for dataset in (
                "profiles/counts_macro_equivalent",
                "profiles/outside_macro_equivalent",
                "final_particles/state",
                "final_particles/zeta_m",
                "final_particles/ptau",
            ):
                if dataset not in h5:
                    raise KeyError(f"missing dataset {dataset}")

            actual_particles = int(h5.attrs["sample_particles"])
            actual_seed = int(h5.attrs["random_seed"])
            actual_source_kind = str(h5.attrs["input_source_kind"])
            if actual_particles != particle_count:
                raise ValueError(
                    f"sample_particles={actual_particles}, expected {particle_count}"
                )
            if actual_seed != seed:
                raise ValueError(f"random_seed={actual_seed}, expected {seed}")
            if actual_source_kind != expected_source_kind:
                raise ValueError(
                    f"input_source_kind={actual_source_kind!r}, "
                    f"expected {expected_source_kind!r}"
                )

            # The writer records profiles in macro-particle-equivalent units.
            # Their in-window plus outside sum must remain the source population.
            counts = np.asarray(h5["profiles/counts_macro_equivalent"], dtype=float)
            outside = np.asarray(h5["profiles/outside_macro_equivalent"], dtype=float)
            represented = counts.sum(axis=1, dtype=np.float64) + outside
            target_population = float(h5.attrs["source_rows"])
            population_conserved = bool(
                np.allclose(represented, target_population, rtol=1e-6, atol=1e-3)
            )

            state = np.asarray(h5["final_particles/state"])
            zeta = np.asarray(h5["final_particles/zeta_m"])
            ptau = np.asarray(h5["final_particles/ptau"])
            alive_finite = int(
                np.count_nonzero((state > 0) & np.isfinite(zeta) & np.isfinite(ptau))
            )
            all_particles_survive = alive_finite == particle_count

            project_source_sha256 = str(h5.attrs["project_source_sha256"])
            if not re.fullmatch(r"[0-9a-f]{64}", project_source_sha256):
                raise ValueError("project_source_sha256 is not a lowercase SHA-256")
            input_sha256 = str(h5.attrs["input_sha256"])
            if source_label == "historical" and input_sha256 != HISTORICAL_INPUT_SHA256:
                raise ValueError("historical input SHA-256 does not match the reference file")

            return {
                "source": source_label,
                "input_source_kind": actual_source_kind,
                "particle_count": particle_count,
                "seed": seed,
                "path": _relative(path),
                "elapsed_tracking_s": float(h5.attrs["elapsed_tracking_s"]),
                "surviving_finite_particles": alive_finite,
                "all_particles_survive_with_finite_final_coordinates": (
                    all_particles_survive
                ),
                "population_conserved": population_conserved,
                "maximum_absolute_population_error": float(
                    np.max(np.abs(represented - target_population))
                ),
                "project_source_sha256": project_source_sha256,
                "input_sha256": input_sha256,
                "input_source_coordinate_sha256": str(
                    h5.attrs["input_source_coordinate_sha256"]
                ),
                "input_tracking_coordinate_sha256": str(
                    h5.attrs["input_tracking_coordinate_sha256"]
                ),
                "invocation_command": str(h5.attrs["invocation_command"]),
            }
    except (KeyError, OSError, TypeError, ValueError) as exc:
        raise ValueError(f"Invalid HDF5 result {_relative(path)}: {exc}") from exc


def _matched_count_summary(particle_count: int) -> dict[str, Any]:
    """Summarize generated-versus-historical agreement at one sample count."""

    report_dir = ANALYSIS_ROOT / f"matched_n{particle_count}_seed202208"
    final_path = report_dir / "final_dt_comparison.json"
    profile_path = report_dir / "rf_result_comparison.json"
    final = _load_json(final_path)
    profile = _load_json(profile_path)

    reference_outside = float(
        _required(
            final,
            "reference",
            "initial_lineage_reference_bucket_retention",
            "fraction_outside_reference_bucket",
        )
    )
    candidate_outside = float(
        _required(
            final,
            "candidate",
            "initial_lineage_reference_bucket_retention",
            "fraction_outside_reference_bucket",
        )
    )
    outside_difference = candidate_outside - reference_outside

    return {
        "particle_count": particle_count,
        "seed": 202_208,
        "final_report": _relative(final_path),
        "profile_report": _relative(profile_path),
        "final": {
            "absolute_empirical_wasserstein_1_ns": float(
                _required(final, "comparison", "absolute_empirical_wasserstein_1_ns")
            ),
            "absolute_histogram_jensen_shannon_distance": float(
                _required(
                    final,
                    "comparison",
                    "absolute_histogram",
                    "jensen_shannon_distance",
                )
            ),
            "absolute_histogram_total_variation_distance": float(
                _required(
                    final,
                    "comparison",
                    "absolute_histogram",
                    "total_variation_distance",
                )
            ),
            "bucket_relative_empirical_wasserstein_1_ns": float(
                _required(
                    final,
                    "comparison",
                    "bucket_relative_empirical_wasserstein_1_ns",
                )
            ),
        },
        "profiles": {
            "median_wasserstein_ns": float(
                _required(profile, "profiles", "wasserstein_summary_ns", "median")
            ),
            "median_jensen_shannon_distance": float(
                _required(profile, "profiles", "jensen_shannon_summary", "median")
            ),
            "common_turns": int(_required(profile, "profiles", "n_common_turns")),
        },
        "outside_initial_reference_bucket": {
            "historical_fraction": reference_outside,
            "generated_fraction": candidate_outside,
            "generated_minus_historical_fraction": outside_difference,
            "absolute_fraction_difference": abs(outside_difference),
        },
    }


def _count_convergence_summary(source_label: str, particle_count: int) -> dict[str, Any]:
    """Summarize a reduced sample's endpoint distance from its full run."""

    report_path = (
        ANALYSIS_ROOT
        / "convergence_to_full"
        / f"{source_label}_n{particle_count}_to_n1075200"
        / "final_dt_comparison.json"
    )
    report = _load_json(report_path)
    return {
        "particle_count": particle_count,
        "full_particle_count": 1_075_200,
        "seed": 202_208,
        "report": _relative(report_path),
        "absolute_empirical_wasserstein_1_ns": float(
            _required(report, "comparison", "absolute_empirical_wasserstein_1_ns")
        ),
        "bucket_relative_empirical_wasserstein_1_ns": float(
            _required(
                report,
                "comparison",
                "bucket_relative_empirical_wasserstein_1_ns",
            )
        ),
    }


def _seed_from_candidate_path(path: str) -> int:
    """Extract the seed encoded in a screen candidate directory name."""

    match = re.search(r"_seed([0-9]+)/rf_only_results\.h5$", path)
    if match is None:
        raise ValueError(f"Cannot extract seed from screen candidate path: {path}")
    return int(match.group(1))


def _seed_screen_summary() -> dict[str, Any]:
    """Reduce the five-seed screen to its baseline, ratios, and assessment."""

    report_path = (
        ANALYSIS_ROOT
        / "screen_n100800_seeds202208_to_202212"
        / "rf_screen_summary.json"
    )
    report = _load_json(report_path)
    generated_pairs = _required(report, "generated_pairs")
    control_pairs = _required(report, "historical_control_pairs")
    if not isinstance(generated_pairs, list) or len(generated_pairs) != 5:
        raise ValueError("The 100,800-particle screen must contain five generated pairs")
    if not isinstance(control_pairs, list) or len(control_pairs) != 4:
        raise ValueError("The 100,800-particle screen must contain four controls")

    generated_ratios: list[dict[str, Any]] = []
    for pair in generated_pairs:
        candidate_path = str(_required(pair, "candidate_path"))
        ratios = _required(pair, "ratios_to_historical_resampling_baseline")
        generated_ratios.append(
            {
                "seed": _seed_from_candidate_path(candidate_path),
                "candidate_path": candidate_path,
                "ratios": {key: float(value) for key, value in ratios.items()},
            }
        )
    generated_ratios.sort(key=lambda item: item["seed"])
    if [item["seed"] for item in generated_ratios] != list(SEEDS):
        raise ValueError("Generated screen seeds do not match 202208 through 202212")

    ratio_keys = tuple(generated_ratios[0]["ratios"])
    ratio_ranges = {
        key: {
            "minimum": min(item["ratios"][key] for item in generated_ratios),
            "maximum": max(item["ratios"][key] for item in generated_ratios),
        }
        for key in ratio_keys
    }
    control_seeds = sorted(
        _seed_from_candidate_path(str(_required(pair, "candidate_path")))
        for pair in control_pairs
    )

    return {
        "particle_count": 100_800,
        "historical_reference_seed": 202_208,
        "generated_seeds": list(SEEDS),
        "historical_control_seeds": control_seeds,
        "report": _relative(report_path),
        "historical_resampling_baseline": _required(
            report, "historical_resampling_baseline"
        ),
        "generated_ratios_to_historical_resampling_baseline": generated_ratios,
        "generated_ratio_ranges": ratio_ranges,
        "assessment": _required(report, "assessment"),
    }


def _runtime_summary(runs: list[dict[str, Any]]) -> dict[str, Any]:
    """Aggregate tracking time while retaining each run's exact duration."""

    elapsed = np.asarray([run["elapsed_tracking_s"] for run in runs], dtype=float)
    by_source: dict[str, dict[str, float]] = {}
    for source in SOURCE_LABELS:
        values = np.asarray(
            [
                run["elapsed_tracking_s"]
                for run in runs
                if run["source"] == source
            ],
            dtype=float,
        )
        by_source[source] = {
            "run_count": int(values.size),
            "total_tracking_s": float(values.sum()),
            "median_tracking_s": float(np.median(values)),
            "maximum_tracking_s": float(values.max()),
        }
    return {
        "run_count": len(runs),
        "total_tracking_s": float(elapsed.sum()),
        "minimum_tracking_s": float(elapsed.min()),
        "median_tracking_s": float(np.median(elapsed)),
        "maximum_tracking_s": float(elapsed.max()),
        "by_source": by_source,
    }


def build_summary() -> dict[str, Any]:
    """Validate all required artifacts and assemble the machine-readable report."""

    runs = [
        _read_run_metadata(source, particle_count, seed, path)
        for source, particle_count, seed, path in _expected_runs()
    ]
    project_source_hashes = sorted({run["project_source_sha256"] for run in runs})
    if len(project_source_hashes) != 1:
        raise ValueError(
            "The matrix mixes executable source hashes: "
            + ", ".join(project_source_hashes)
        )

    matched_counts = [_matched_count_summary(count) for count in PARTICLE_COUNTS]
    convergence_to_full = {
        source: [_count_convergence_summary(source, count) for count in REDUCED_COUNTS]
        for source in SOURCE_LABELS
    }
    seed_screen = _seed_screen_summary()
    all_survive = all(
        run["all_particles_survive_with_finite_final_coordinates"] for run in runs
    )
    all_conserved = all(run["population_conserved"] for run in runs)

    return {
        "schema_version": 1,
        "study": {
            "name": "5 kV RF-only particle-count and seed convergence matrix",
            "rf25_initial_v": 5000.0,
            "particle_counts": list(PARTICLE_COUNTS),
            "seed_202208_matched_count_pairs": len(PARTICLE_COUNTS),
            "seed_screen_particle_count": 100_800,
            "seed_screen_seeds": list(SEEDS),
            "expected_hdf5_run_count": 16,
        },
        "matched_generated_to_historical": matched_counts,
        "particle_count_convergence_to_full": convergence_to_full,
        "seed_screen": seed_screen,
        "integrity": {
            "all_16_runs_survive_with_finite_final_coordinates": all_survive,
            "all_16_runs_conserve_macro_equivalent_population": all_conserved,
            "all_16_runs_share_project_source_sha256": True,
            "project_source_sha256": project_source_hashes[0],
            "historical_input_sha256": HISTORICAL_INPUT_SHA256,
        },
        "runtime": _runtime_summary(runs),
        "runs": runs,
    }


def _plot_dual_axis(
    axis: plt.Axes,
    counts: np.ndarray,
    left_values: np.ndarray,
    right_values: np.ndarray,
    left_label: str,
    right_label: str,
    title: str,
) -> None:
    """Draw two count-dependent metrics with separate, color-matched y axes."""

    right_axis = axis.twinx()
    left_line = axis.plot(
        counts, left_values, "o-", color="tab:blue", label=left_label
    )[0]
    right_line = right_axis.plot(
        counts, right_values, "s-", color="tab:orange", label=right_label
    )[0]
    axis.set_xscale("log")
    axis.set_xlabel("Tracked particles")
    axis.set_ylabel(left_label, color="tab:blue")
    right_axis.set_ylabel(right_label, color="tab:orange")
    axis.tick_params(axis="y", labelcolor="tab:blue")
    right_axis.tick_params(axis="y", labelcolor="tab:orange")
    axis.grid(True, which="both", alpha=0.25)
    axis.set_title(title)
    axis.legend(handles=[left_line, right_line], loc="best", fontsize=8)


def make_figure(summary: dict[str, Any], output_path: Path) -> None:
    """Plot count convergence and the five generated-to-control seed ratios."""

    matched = summary["matched_generated_to_historical"]
    counts = np.asarray([row["particle_count"] for row in matched], dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    final_w1 = np.asarray(
        [row["final"]["absolute_empirical_wasserstein_1_ns"] for row in matched]
    )
    final_js = np.asarray(
        [
            row["final"]["absolute_histogram_jensen_shannon_distance"]
            for row in matched
        ]
    )
    _plot_dual_axis(
        axes[0, 0],
        counts,
        final_w1,
        final_js,
        "Final empirical W1 [ns]",
        "Final histogram JS distance",
        "Generated vs historical final distribution",
    )

    profile_w1 = np.asarray(
        [row["profiles"]["median_wasserstein_ns"] for row in matched]
    )
    profile_js = np.asarray(
        [row["profiles"]["median_jensen_shannon_distance"] for row in matched]
    )
    _plot_dual_axis(
        axes[0, 1],
        counts,
        profile_w1,
        profile_js,
        "Median profile W1 [ns]",
        "Median profile JS distance",
        "Generated vs historical time-resolved profiles",
    )

    count_axis = axes[1, 0]
    convergence = summary["particle_count_convergence_to_full"]
    for source, marker, color in (
        ("historical", "o", "tab:blue"),
        ("generated", "s", "tab:orange"),
    ):
        rows = convergence[source]
        source_counts = [row["particle_count"] for row in rows] + [1_075_200]
        values = [row["absolute_empirical_wasserstein_1_ns"] for row in rows] + [0.0]
        count_axis.plot(
            source_counts,
            values,
            marker=marker,
            color=color,
            label=f"{source.capitalize()} reduced vs full",
        )
    count_axis.set_xscale("log")
    count_axis.set_xlabel("Tracked particles")
    count_axis.set_ylabel("Final empirical W1 to 1,075,200 run [ns]")
    count_axis.set_title("Particle-count convergence within each input source")
    count_axis.grid(True, which="both", alpha=0.25)
    count_axis.legend(fontsize=8)

    ratio_axis = axes[1, 1]
    seed_screen = summary["seed_screen"]
    ratio_labels = {
        "median_profile_wasserstein": "Profile W1",
        "median_profile_jensen_shannon": "Profile JS",
        "mean_absolute_centroid_difference": "Mean centroid",
        "maximum_absolute_centroid_difference": "Max centroid",
    }
    keys = list(ratio_labels)
    x_positions = np.arange(len(keys), dtype=float)
    seed_rows = seed_screen["generated_ratios_to_historical_resampling_baseline"]
    for row in seed_rows:
        values = [row["ratios"][key] for key in keys]
        ratio_axis.plot(
            x_positions,
            values,
            "o-",
            alpha=0.55,
            linewidth=1.0,
            label=str(row["seed"]),
        )
    ranges = seed_screen["generated_ratio_ranges"]
    for x_position, key in zip(x_positions, keys, strict=True):
        ratio_axis.vlines(
            x_position,
            ranges[key]["minimum"],
            ranges[key]["maximum"],
            color="black",
            linewidth=4,
            alpha=0.35,
        )
    ratio_axis.axhline(1.0, color="black", linestyle=":", label="Control median")
    ratio_axis.axhline(2.0, color="tab:red", linestyle="--", label="2x screen limit")
    ratio_axis.set_xticks(x_positions, [ratio_labels[key] for key in keys])
    ratio_axis.set_ylabel("Generated metric / historical-resampling median")
    ratio_axis.set_title("100,800-particle generated-seed stability")
    ratio_axis.grid(True, axis="y", alpha=0.25)
    ratio_axis.legend(ncol=2, fontsize=8)

    integrity = summary["integrity"]
    runtime_hours = summary["runtime"]["total_tracking_s"] / 3600.0
    pass_text = "PASS" if seed_screen["assessment"]["passes_all_criteria"] else "FAIL"
    fig.suptitle(
        "Recycler RF-only 5 kV convergence summary\n"
        f"Seed screen: {pass_text}; all runs survive/conserve: "
        f"{integrity['all_16_runs_survive_with_finite_final_coordinates']}/"
        f"{integrity['all_16_runs_conserve_macro_equivalent_population']}; "
        f"summed tracking time: {runtime_hours:.2f} h; "
        f"source {integrity['project_source_sha256'][:12]}…",
        fontsize=13,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=170)
    plt.close(fig)


def main() -> None:
    """Write the combined JSON summary and its convergence overview plot."""

    summary = build_summary()
    ANALYSIS_ROOT.mkdir(parents=True, exist_ok=True)
    json_path = ANALYSIS_ROOT / "convergence_summary.json"
    figure_path = ANALYSIS_ROOT / "convergence_summary.png"
    json_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    make_figure(summary, figure_path)
    print(f"Wrote {_relative(json_path)}")
    print(f"Wrote {_relative(figure_path)}")


if __name__ == "__main__":
    main()
