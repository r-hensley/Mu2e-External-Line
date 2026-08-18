"""Run and validate the Xsuite Recycler filamentation simulation."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import xpart as xp
import xtrack as xt
import xsuite

try:
    # The preserved module is read-only here.  It supplies the exact same
    # deterministic starting ensemble so the tracking engines can be compared
    # without conflating an input-distribution change.
    from ..manual_version import mu2e_longitudinal as reference_model
except ImportError:  # Support direct execution with ``python run.py``.
    manual_directory = Path(__file__).resolve().parents[1] / "manual_version"
    if str(manual_directory) not in sys.path:
        sys.path.insert(0, str(manual_directory))
    import mu2e_longitudinal as reference_model

try:
    from .plotting import (
        plot_checkpoints,
        plot_diagnostic_comparison,
        plot_extraction_gallery,
        plot_phase_space_comparison,
        plot_rf_program,
        plot_waterfall,
    )
    from .simulation import (
        RecyclerInputs,
        build_line,
        bunch_diagnostics,
        derive_from_xsuite,
        result_digest,
        track_rebunching,
        validate_one_turn,
    )
except ImportError:  # Support direct execution with ``python run.py``.
    from plotting import (
        plot_checkpoints,
        plot_diagnostic_comparison,
        plot_extraction_gallery,
        plot_phase_space_comparison,
        plot_rf_program,
        plot_waterfall,
    )
    from simulation import (
        RecyclerInputs,
        build_line,
        bunch_diagnostics,
        derive_from_xsuite,
        result_digest,
        track_rebunching,
        validate_one_turn,
    )


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return parsed


def _json_value(value: Any) -> Any:
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
        # Preserve an absolute path when the user deliberately supplies a
        # reference file from outside this example directory.
        return str(resolved)


def _save_result(path: Path, result) -> None:
    """Persist the Xsuite trajectory using the reference result field names."""

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
        final_zeta_m=result.final_zeta_m,
        final_ptau=result.final_ptau,
        final_state=result.final_state,
    )
    print(f"Wrote {path}")


def compare_to_reference(result, reference: np.lib.npyio.NpzFile) -> dict[str, Any]:
    """Compare particle trajectories and distribution moments on common turns."""

    if not np.array_equal(result.snapshot_turns, reference["snapshot_turns"]):
        raise RuntimeError("Xsuite and reference snapshot turn grids differ")
    if not np.array_equal(result.diagnostic_turns, reference["diagnostic_turns"]):
        raise RuntimeError("Xsuite and reference diagnostic turn grids differ")

    h = result.inputs.harmonic_25
    # Use a circular phase difference so equivalent positions on opposite
    # sides of the wrapped h=28 cell are not reported as a full-cell error.
    phase_difference = np.angle(
        np.exp(1j * h * (result.theta_rad.astype(float) - reference["theta_rad"].astype(float)))
    ) / h
    time_difference_ns = (
        phase_difference / result.derived.angular_revolution_frequency * 1e9
    )
    energy_difference_eV = (
        result.delta_energy_eV.astype(float)
        - reference["delta_energy_eV"].astype(float)
    )

    flat_turn = int(round(result.derived.time_25_flat_s / result.derived.t_rev0_s))
    flat_index = int(np.argmin(np.abs(result.snapshot_turns - flat_turn)))
    final_index = -1

    def endpoint(index: int) -> dict[str, float]:
        xs_diag = bunch_diagnostics(
            result.theta_rad[index],
            result.delta_energy_eV[index],
            result.inputs,
            result.derived,
        )
        ref_diag = bunch_diagnostics(
            reference["theta_rad"][index],
            reference["delta_energy_eV"][index],
            result.inputs,
            result.derived,
        )
        return {
            "turn": int(result.snapshot_turns[index]),
            "xsuite_sigma_time_ns": xs_diag["sigma_time_s"] * 1e9,
            "reference_sigma_time_ns": ref_diag["sigma_time_s"] * 1e9,
            "sigma_time_relative_difference": (
                xs_diag["sigma_time_s"] / ref_diag["sigma_time_s"] - 1.0
            ),
            "xsuite_sigma_energy_MeV": xs_diag["sigma_energy_eV"] / 1e6,
            "reference_sigma_energy_MeV": ref_diag["sigma_energy_eV"] / 1e6,
            "sigma_energy_relative_difference": (
                xs_diag["sigma_energy_eV"] / ref_diag["sigma_energy_eV"] - 1.0
            ),
            "particlewise_time_rms_difference_ns": float(
                np.sqrt(np.mean(time_difference_ns[index] ** 2))
            ),
            "particlewise_time_max_abs_difference_ns": float(
                np.max(np.abs(time_difference_ns[index]))
            ),
            "particlewise_energy_rms_difference_eV": float(
                np.sqrt(np.mean(energy_difference_eV[index] ** 2))
            ),
            "particlewise_energy_max_abs_difference_eV": float(
                np.max(np.abs(energy_difference_eV[index]))
            ),
        }

    comparison = {
        "snapshot_turn_grids_identical": True,
        "diagnostic_turn_grids_identical": True,
        "flattop": endpoint(flat_index),
        "last_extraction": endpoint(final_index),
        "all_snapshots_particlewise_time_rms_difference_ns": float(
            np.sqrt(np.mean(time_difference_ns**2))
        ),
        "all_snapshots_particlewise_energy_rms_difference_eV": float(
            np.sqrt(np.mean(energy_difference_eV**2))
        ),
    }
    return comparison


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the Xsuite version of the Recycler filamentation example."
    )
    folder = Path(__file__).resolve().parent
    parser.add_argument("--output-dir", type=Path, default=folder / "outputs")
    parser.add_argument("--particles-per-source", type=_positive_int, default=256)
    parser.add_argument("--diagnostic-stride-turns", type=_positive_int, default=20)
    parser.add_argument(
        "--reference-result",
        type=Path,
        default=folder.parent / "manual_version" / "outputs" / "nominal_result.npz",
        help="manual NumPy result used for the trajectory comparison",
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

    inputs = RecyclerInputs()

    # Xsuite owns the reference-particle calculations, longitudinal map, RF
    # kicks, and Twiss-derived ring quantities from this point onward.
    line = build_line(inputs)
    derived = derive_from_xsuite(line, inputs)
    one_turn = validate_one_turn(line, inputs, derived)
    assert abs(one_turn["theta_error_rad"]) < 2e-12
    assert abs(one_turn["energy_error_eV"]) < 1e-4
    assert inputs.harmonic_53 // inputs.harmonic_25 == 21
    assert abs(derived.synchrotron_period_80kv_s * 1e3 - 18.4) < 0.1

    print(
        f"Xsuite {xsuite.__version__}; "
        f"Xtrack {xt.__version__}; Xpart {xp.__version__}"
    )
    print(
        "Xsuite-derived reference: "
        f"p0c={derived.p0c_eV/1e9:.9f} GeV/c, "
        f"beta0={derived.beta0:.12f}, gamma0={derived.gamma0:.12f}"
    )
    print(
        "Xsuite-derived ring: "
        f"Trev={derived.t_rev0_s*1e6:.9f} us, "
        f"eta={derived.slip_factor:.12g}, "
        f"Qs(80 kV)={derived.qs_80kv:.12g}"
    )
    print("One-turn Xsuite convention check: PASS")

    reference_parameters = reference_model.RecyclerParameters()
    initial = reference_model.make_initial_ensemble(
        particles_per_source=args.particles_per_source,
        params=reference_parameters,
    )
    result = track_rebunching(
        line,
        initial.theta_rad,
        initial.delta_energy_eV,
        initial.source_id,
        inputs,
        derived,
        diagnostic_stride_turns=args.diagnostic_stride_turns,
        display_per_source=120,
        show_progress=True,
    )
    print(
        f"Tracked {result.n_particles:,} particles for "
        f"{result.snapshot_turns[-1]:,} turns in {result.tracking_elapsed_s:.3f} s"
    )

    reference_path = args.reference_result.resolve()
    if not reference_path.is_file():
        raise FileNotFoundError(f"preserved reference result not found: {reference_path}")
    reference = np.load(reference_path)
    comparison = compare_to_reference(result, reference)

    # Filamentation makes individual trajectories sensitive to tiny numerical
    # differences over many synchrotron periods.  Acceptance is therefore set
    # on the physically relevant distribution widths, while particlewise
    # differences are still recorded for transparency.
    for key in ("flattop", "last_extraction"):
        assert abs(comparison[key]["sigma_time_relative_difference"]) < 0.02
        assert abs(comparison[key]["sigma_energy_relative_difference"]) < 0.02
    print("Distribution-level comparison with preserved NumPy result: PASS")

    _save_figure(plot_rf_program(inputs, derived), output_dir / "01_rf_program.png")
    _save_figure(plot_checkpoints(result), output_dir / "02_xsuite_checkpoints.png")
    _save_figure(plot_waterfall(result), output_dir / "03_xsuite_waterfall.png")
    _save_figure(
        plot_extraction_gallery(result), output_dir / "04_xsuite_extraction_gallery.png"
    )
    _save_figure(
        plot_phase_space_comparison(result, reference),
        output_dir / "05_flattop_phase_space_comparison.png",
    )
    _save_figure(
        plot_diagnostic_comparison(result, reference),
        output_dir / "06_diagnostic_comparison.png",
    )

    result_path = output_dir / "xsuite_result.npz"
    _save_result(result_path, result)
    comparison_path = output_dir / "comparison_with_numpy.json"
    comparison_path.write_text(
        json.dumps(_json_value(comparison), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {comparison_path}")

    digest = result_digest(result)
    summary = {
        "model": "Xsuite LineSegmentMap plus harmonic Cavity elements",
        "xsuite_version": xsuite.__version__,
        "xtrack_version": xt.__version__,
        "xpart_version": xp.__version__,
        "inputs": asdict(inputs),
        "xsuite_derived": asdict(derived),
        "line_element_order": ["recycler_slip", "rf53", "rf2p5"],
        "initial_distribution_source": _portable_path(Path(reference_model.__file__)),
        "reference_result": _portable_path(reference_path),
        "one_turn_validation": one_turn,
        "n_particles": result.n_particles,
        "n_turns": int(result.snapshot_turns[-1]),
        "tracking_elapsed_s": result.tracking_elapsed_s,
        "all_particles_survived": bool(np.all(result.final_state > 0)),
        "all_coordinates_finite": bool(
            np.isfinite(result.theta_rad).all()
            and np.isfinite(result.delta_energy_eV).all()
        ),
        "comparison": comparison,
        "result_sha256": digest,
        "all_acceptance_checks_passed": True,
    }
    summary_path = output_dir / "run_summary.json"
    summary_path.write_text(
        json.dumps(_json_value(summary), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {summary_path}")
    print(f"Xsuite result SHA-256: {digest}")
    reference.close()
    return summary


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
