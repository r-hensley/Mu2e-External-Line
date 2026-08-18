"""Run the notebook-independent Recycler filamentation demonstration.

This is a direct Python entry point for the numerical model and plots formerly
orchestrated by ``Mu2e_Recycler_Rebunching_Live_Demo.ipynb``.  The underlying
model remains the explicit NumPy drift--kick map; it has not yet been converted
to Xsuite.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import time
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

try:
    from . import mu2e_longitudinal as m2
except ImportError:  # Support direct execution with ``python run.py``.
    import mu2e_longitudinal as m2


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return parsed


def _json_value(value: Any) -> Any:
    """Convert NumPy/Pandas scalar containers to strict JSON values."""

    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return value


def _save_figure(figure: plt.Figure, path: Path) -> None:
    figure.savefig(path, bbox_inches="tight")
    plt.close(figure)
    print(f"Wrote {path}")


def _portable_path(path: Path) -> str:
    """Prefer a project-relative path in metadata written inside the project."""

    resolved = path.resolve()
    project_root = Path(__file__).resolve().parents[1]
    try:
        return str(resolved.relative_to(project_root))
    except ValueError:
        return str(resolved)


def _array_digest(*arrays: np.ndarray) -> str:
    digest = sha256()
    for array in arrays:
        contiguous = np.ascontiguousarray(array)
        digest.update(contiguous.dtype.str.encode("ascii"))
        digest.update(str(contiguous.shape).encode("ascii"))
        digest.update(contiguous.view(np.uint8))
    return digest.hexdigest()


def _save_nominal_result(path: Path, result: m2.SimulationResult) -> None:
    np.savez_compressed(
        path,
        source_id=result.source_id,
        snapshot_turns=result.snapshot_turns,
        snapshot_times_s=result.snapshot_times_s,
        theta_rad=result.theta_rad,
        delta_energy_eV=result.delta_energy_eV,
        display_indices=result.display_indices,
        diagnostic_turns=result.diagnostic_turns,
        diagnostic_times_s=result.diagnostic_times_s,
        sigma_time_s=result.sigma_time_s,
        sigma_energy_eV=result.sigma_energy_eV,
        covariance_emittance_eVs=result.covariance_emittance_eVs,
        central_95_time_width_s=result.central_95_time_width_s,
        outside_window_fraction=result.outside_window_fraction,
        centroid_time_s=result.centroid_time_s,
    )
    print(f"Wrote {path}")


def _acceptance_checks(
    initial: m2.InitialEnsemble,
    nominal: m2.SimulationResult,
    params: m2.RecyclerParameters,
) -> tuple[list[dict[str, Any]], float, float, float]:
    jacobian = m2.numerical_map_jacobian(
        theta_rad=math.radians(1.0),
        delta_energy_eV=2.0e6,
        time_s=params.time_25_flat_s + params.revolution_period_s,
        params=params,
    )
    determinant = float(np.linalg.det(jacobian))
    measured_period_ms = m2.measured_small_amplitude_period_s(params) * 1e3
    analytical_period_ms = m2.small_amplitude_period_s(
        params.harmonic_25,
        params.voltage_25_final_eV,
        params,
    ) * 1e3
    solved_area = m2.matched_contour_area_eVs(
        initial.contour_phase_amplitude_rad,
        params.harmonic_53,
        params.voltage_53_eV,
        params,
    )

    checks = [
        {
            "check": "source bunch count",
            "value": initial.n_particles // initial.particles_per_source,
            "target": params.source_bunch_count,
        },
        {
            "check": "initial contour area [eV-s]",
            "value": solved_area,
            "target": params.bunch_emittance_eVs,
        },
        {
            "check": "one-turn det(J)",
            "value": determinant,
            "target": 1.0,
        },
        {
            "check": "tracked small-amplitude period [ms]",
            "value": measured_period_ms,
            "target": analytical_period_ms,
        },
        {
            "check": "all nominal coordinates finite",
            "value": float(
                np.isfinite(nominal.theta_rad).all()
                and np.isfinite(nominal.delta_energy_eV).all()
            ),
            "target": 1.0,
        },
    ]
    for check in checks:
        check["difference"] = float(check["value"] - check["target"])

    assert initial.n_particles == 21 * initial.particles_per_source
    assert abs(solved_area - params.bunch_emittance_eVs) < 2e-10
    assert abs(determinant - 1.0) < 2e-7
    assert abs(measured_period_ms - analytical_period_ms) < 0.03
    assert np.isfinite(nominal.theta_rad).all()
    assert np.isfinite(nominal.delta_energy_eV).all()
    assert np.max(np.abs(nominal.theta_rad)) <= (
        math.pi / params.harmonic_25 + 1e-7
    )
    return checks, determinant, measured_period_ms, analytical_period_ms


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run the standalone NumPy Recycler rebunching and filamentation "
            "example and save its figures and numerical results."
        )
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "outputs",
        help="directory for PNG, HTML, NPZ, and JSON products",
    )
    parser.add_argument(
        "--particles-per-source",
        type=_positive_int,
        default=256,
        help="macroparticles in each of the 21 source bunchlets (default: 256)",
    )
    parser.add_argument(
        "--diagnostic-stride-turns",
        type=_positive_int,
        default=20,
        help="turn spacing for scalar diagnostics (default: 20)",
    )
    parser.add_argument(
        "--display-per-source",
        type=_positive_int,
        default=120,
        help="display sample retained per source bunchlet (default: 120)",
    )
    parser.add_argument(
        "--no-dashboard",
        action="store_true",
        help="skip the self-contained Plotly dashboard HTML",
    )
    return parser


def run(args: argparse.Namespace) -> dict[str, Any]:
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.dpi": 150,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
        }
    )

    params = m2.RecyclerParameters()
    harmonic_ratio = params.harmonic_53 / params.harmonic_25
    bucket_height_mev = m2.bucket_half_height_eV(
        params.harmonic_25,
        params.voltage_25_final_eV,
        params,
    ) / 1e6
    analytical_period_ms = m2.small_amplitude_period_s(
        params.harmonic_25,
        params.voltage_25_final_eV,
        params,
    ) * 1e3
    assert harmonic_ratio == 21
    assert abs(bucket_height_mev - 42.1) < 0.1
    assert abs(analytical_period_ms - 18.4) < 0.1

    print(f"Model version: {m2.MODEL_VERSION}")
    print(f"Helper source: {m2.__file__}")
    print(f"NumPy: {np.__version__}")
    print("Analytical acceptance checks: PASS")

    initial = m2.make_initial_ensemble(
        particles_per_source=args.particles_per_source,
        params=params,
    )

    _save_figure(m2.plot_rf_program(params), output_dir / "01_rf_program.png")
    _save_figure(
        m2.plot_initial_ensemble(initial, params),
        output_dir / "02_initial_ensemble.png",
    )

    tracking_start = time.perf_counter()
    nominal = m2.track_rebunching(
        initial,
        params=params,
        diagnostic_stride_turns=args.diagnostic_stride_turns,
        display_per_source=args.display_per_source,
    )
    tracking_elapsed_s = time.perf_counter() - tracking_start
    print(
        f"Tracked {nominal.n_particles:,} particles for "
        f"{nominal.snapshot_turns[-1]:,} turns in {tracking_elapsed_s:.3f} s"
    )

    _save_figure(
        m2.plot_checkpoint_snapshots(nominal),
        output_dir / "03_checkpoint_snapshots.png",
    )
    if not args.no_dashboard:
        dashboard = m2.plot_preloaded_snapshot_dashboard(
            nominal,
            display_per_source=args.display_per_source,
            phase_bins=100,
            frame_duration_ms=140,
        )
        dashboard_path = output_dir / "04_preloaded_dashboard.html"
        dashboard.write_html(
            dashboard_path,
            include_plotlyjs=True,
            full_html=True,
            auto_play=False,
        )
        print(f"Wrote {dashboard_path}")

    _save_figure(
        m2.plot_tune_spread(params),
        output_dir / "05_tune_spread.png",
    )

    control_time_s = float(params.extraction_times_s[1])
    linear_control = m2.track_rebunching(
        initial,
        params=params,
        stop_time_s=control_time_s,
        snapshot_times_s=[control_time_s],
        diagnostic_stride_turns=100,
        display_per_source=args.display_per_source,
        linearized_25=True,
    )
    _save_figure(
        m2.plot_nonlinear_control_comparison(
            nominal,
            linear_control,
            time_s=control_time_s,
        ),
        output_dir / "06_nonlinear_control.png",
    )
    _save_figure(
        m2.plot_waterfall(nominal),
        output_dir / "07_waterfall.png",
    )
    _save_figure(
        m2.plot_diagnostics(nominal),
        output_dir / "08_diagnostics.png",
    )
    _save_figure(
        m2.plot_extraction_gallery(nominal),
        output_dir / "09_extraction_gallery.png",
    )

    checks, determinant, measured_period_ms, analytical_period_ms = (
        _acceptance_checks(initial, nominal, params)
    )
    print("All numerical acceptance checks: PASS")

    result_path = output_dir / "nominal_result.npz"
    _save_nominal_result(result_path, nominal)

    helper_path = Path(m2.__file__).resolve()
    source_digest = sha256(helper_path.read_bytes()).hexdigest()
    result_digest = _array_digest(
        nominal.source_id,
        nominal.snapshot_turns,
        nominal.snapshot_times_s,
        nominal.theta_rad,
        nominal.delta_energy_eV,
        nominal.diagnostic_turns,
        nominal.sigma_time_s,
        nominal.sigma_energy_eV,
        nominal.covariance_emittance_eVs,
        nominal.central_95_time_width_s,
        nominal.outside_window_fraction,
        nominal.centroid_time_s,
    )
    final_diagnostics = m2.bunch_diagnostics(
        nominal.theta_rad[-1].astype(float),
        nominal.delta_energy_eV[-1].astype(float),
        params,
    )
    summary: dict[str, Any] = {
        "model": "standalone NumPy drift-kick Recycler filamentation example",
        "model_version": m2.MODEL_VERSION,
        "helper_source": _portable_path(helper_path),
        "helper_sha256": source_digest,
        "result_sha256": result_digest,
        "numpy_version": np.__version__,
        "particles_per_source": initial.particles_per_source,
        "source_bunch_count": params.source_bunch_count,
        "n_particles": nominal.n_particles,
        "n_turns": int(nominal.snapshot_turns[-1]),
        "end_time_ms": float(nominal.snapshot_times_s[-1] * 1e3),
        "tracking_elapsed_s": tracking_elapsed_s,
        "harmonic_ratio": harmonic_ratio,
        "bucket_half_height_MeV": bucket_height_mev,
        "analytical_period_ms": analytical_period_ms,
        "measured_period_ms": measured_period_ms,
        "one_turn_jacobian_determinant": determinant,
        "initial_contour_area_eVs": initial.contour_area_eVs,
        "final_diagnostics": final_diagnostics,
        "acceptance_checks": checks,
        "all_acceptance_checks_passed": True,
    }
    summary_path = output_dir / "run_summary.json"
    summary_path.write_text(
        json.dumps(_json_value(summary), indent=2, sort_keys=True, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {summary_path}")
    print(f"Result SHA-256: {result_digest}")
    return summary


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
