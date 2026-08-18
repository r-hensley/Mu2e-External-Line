"""Quantitative comparisons for two RF-only reconstruction result files.

The profile comparison is intentionally independent of tracked sample size:
histograms are normalized row by row after the recorded-turn grids are
intersected.  Population conservation and the fraction outside the profile
window are reported separately so normalization cannot hide lost or clipped
mass.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from .config import C_LIGHT_M_PER_S


@dataclass(frozen=True)
class ProfileData:
    """Turn-indexed longitudinal histograms and their out-of-window population."""

    turns: np.ndarray
    time_edges_ns: np.ndarray
    counts_macro_equivalent: np.ndarray
    outside_macro_equivalent: np.ndarray


@dataclass(frozen=True)
class FinalParticleData:
    """Final Xsuite coordinates, survival state, and optional input lineage.

    ``zeta_m`` is Xsuite's longitudinal position, ``ptau`` is the normalized
    energy coordinate, and ``bunch_index`` identifies the initial structural
    bunch (0--7). The 0--167 microbunch index is not stored in result files.
    """

    zeta_m: np.ndarray
    ptau: np.ndarray
    state: np.ndarray
    bunch_index: np.ndarray | None


@dataclass(frozen=True)
class ResultData:
    """In-memory subset of the HDF5 schema needed for comparisons."""

    path: Path
    root_attrs: dict[str, Any]
    machine_attrs: dict[str, Any]
    rf_program_attrs: dict[str, Any]
    profiles: ProfileData
    final_particles: FinalParticleData | None


def _plain_scalar(value: Any) -> Any:
    """Convert NumPy or byte-valued HDF5 attributes to ordinary Python values."""

    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return value


def _attrs(group: h5py.Group | h5py.File) -> dict[str, Any]:
    """Read all attributes from an HDF5 group into a JSON-friendly dictionary."""

    return {key: _plain_scalar(value) for key, value in group.attrs.items()}


def _require_dataset(h5: h5py.File, name: str) -> np.ndarray:
    """Read one required HDF5 dataset or raise a schema-oriented error."""

    if name not in h5:
        raise ValueError(f"Missing required HDF5 dataset: {name}")
    return np.asarray(h5[name][:])


def _validate_profiles(path: Path, profiles: ProfileData) -> None:
    """Validate profile dimensions, monotonic grids, and population values."""

    turns = profiles.turns
    edges = profiles.time_edges_ns
    counts = profiles.counts_macro_equivalent
    outside = profiles.outside_macro_equivalent

    if turns.ndim != 1 or turns.size == 0:
        raise ValueError(f"{path}: profiles/turn must be a nonempty 1D array")
    if np.any(np.diff(turns) <= 0):
        raise ValueError(f"{path}: profiles/turn must be strictly increasing")
    if edges.ndim != 1 or edges.size < 2 or not np.all(np.isfinite(edges)):
        raise ValueError(f"{path}: profiles/time_edges_ns must contain finite bin edges")
    if np.any(np.diff(edges) <= 0):
        raise ValueError(f"{path}: profile time-bin edges must be strictly increasing")
    if counts.shape != (turns.size, edges.size - 1):
        raise ValueError(
            f"{path}: profile counts shape {counts.shape} does not match "
            f"{turns.size} turns and {edges.size - 1} bins"
        )
    if outside.shape != (turns.size,):
        raise ValueError(f"{path}: profiles/outside_macro_equivalent has the wrong shape")
    if not np.all(np.isfinite(counts)) or not np.all(np.isfinite(outside)):
        raise ValueError(f"{path}: profile populations must be finite")
    if np.any(counts < 0) or np.any(outside < 0):
        raise ValueError(f"{path}: profile populations must be nonnegative")


def _load_final_particles(h5: h5py.File, path: Path) -> FinalParticleData | None:
    """Load and validate optional final-particle coordinates from a result file."""

    if "final_particles" not in h5:
        return None
    required = ["zeta_m", "ptau", "state"]
    missing = [name for name in required if f"final_particles/{name}" not in h5]
    if missing:
        raise ValueError(f"{path}: incomplete final_particles group; missing {missing}")

    zeta = np.asarray(h5["final_particles/zeta_m"][:], dtype=float)
    ptau = np.asarray(h5["final_particles/ptau"][:], dtype=float)
    state = np.asarray(h5["final_particles/state"][:])
    labels = None
    if "final_particles/bunch_index" in h5:
        labels = np.asarray(h5["final_particles/bunch_index"][:])

    if zeta.ndim != 1 or ptau.ndim != 1 or state.ndim != 1:
        raise ValueError(f"{path}: final-particle datasets must be 1D")
    if not (zeta.size == ptau.size == state.size):
        raise ValueError(f"{path}: final-particle datasets have inconsistent lengths")
    if labels is not None and (labels.ndim != 1 or labels.size != zeta.size):
        raise ValueError(f"{path}: final bunch labels have an inconsistent length")
    return FinalParticleData(zeta_m=zeta, ptau=ptau, state=state, bunch_index=labels)


def load_result_data(path: str | Path) -> ResultData:
    """Load and validate the comparison-relevant subset of a result file.

    Profiles are stored as one histogram row per recorded turn.  Particle
    coordinates are optional because older or profile-only files may omit them.
    """

    path = Path(path)
    with h5py.File(path, "r") as h5:
        # Keep only the common, compact schema needed by all comparison tools.
        profiles = ProfileData(
            turns=_require_dataset(h5, "profiles/turn").astype(np.int64, copy=False),
            time_edges_ns=_require_dataset(h5, "profiles/time_edges_ns").astype(
                float, copy=False
            ),
            counts_macro_equivalent=_require_dataset(
                h5, "profiles/counts_macro_equivalent"
            ).astype(float, copy=False),
            outside_macro_equivalent=_require_dataset(
                h5, "profiles/outside_macro_equivalent"
            ).astype(float, copy=False),
        )
        _validate_profiles(path, profiles)
        machine_attrs = _attrs(h5["machine"]) if "machine" in h5 else {}
        return ResultData(
            path=path,
            root_attrs=_attrs(h5),
            machine_attrs=machine_attrs,
            rf_program_attrs=_attrs(h5["rf_program"])
            if "rf_program" in h5
            else {},
            profiles=profiles,
            final_particles=_load_final_particles(h5, path),
        )


def _normalized(values: np.ndarray) -> np.ndarray | None:
    """Normalize nonzero histogram mass to one, or return ``None`` if empty."""

    total = float(np.sum(values, dtype=np.float64))
    if total <= 0:
        return None
    return np.asarray(values, dtype=float) / total


def _histogram_wasserstein_ns(
    reference: np.ndarray,
    candidate: np.ndarray,
    time_edges_ns: np.ndarray,
) -> float | None:
    """Wasserstein-1 distance for masses placed at histogram-bin centers."""
    p = _normalized(reference)
    q = _normalized(candidate)
    if p is None or q is None:
        return None
    if p.size == 1:
        return 0.0
    centers = 0.5 * (time_edges_ns[:-1] + time_edges_ns[1:])
    return float(np.sum(np.abs(np.cumsum(p - q)[:-1]) * np.diff(centers)))


def _jensen_shannon_distance(
    reference: np.ndarray,
    candidate: np.ndarray,
) -> float | None:
    """Base-2 Jensen-Shannon distance, bounded between zero and one."""
    p = _normalized(reference)
    q = _normalized(candidate)
    if p is None or q is None:
        return None
    midpoint = 0.5 * (p + q)
    p_mask = p > 0
    q_mask = q > 0
    divergence = 0.5 * np.sum(p[p_mask] * np.log2(p[p_mask] / midpoint[p_mask]))
    divergence += 0.5 * np.sum(q[q_mask] * np.log2(q[q_mask] / midpoint[q_mask]))
    return float(np.sqrt(max(float(divergence), 0.0)))


def _finite_summary(values: list[float | None]) -> dict[str, float | int | None]:
    """Summarize available scalar metrics while ignoring unavailable entries."""

    finite = np.asarray([value for value in values if value is not None], dtype=float)
    if finite.size == 0:
        return {"n_valid": 0, "mean": None, "median": None, "maximum": None}
    return {
        "n_valid": int(finite.size),
        "mean": float(np.mean(finite)),
        "median": float(np.median(finite)),
        "maximum": float(np.max(finite)),
    }


def _conservation(result: ResultData) -> dict[str, Any]:
    """Check represented profile population against the source distribution."""

    profiles = result.profiles
    # The writer partitions every recorded turn into histogrammed mass and mass
    # outside the plotting window.  Adding both prevents clipping from looking
    # like particle loss when compared with the source-row count.
    represented = np.sum(profiles.counts_macro_equivalent, axis=1, dtype=np.float64)
    represented += profiles.outside_macro_equivalent
    target_raw = result.root_attrs.get("source_rows")
    target = float(target_raw) if target_raw is not None else None
    answer: dict[str, Any] = {
        "represented_min": float(np.min(represented)),
        "represented_max": float(np.max(represented)),
        "max_turn_to_turn_range": float(np.max(represented) - np.min(represented)),
        "target_macro_equivalent": target,
        "max_absolute_error": None,
        "max_relative_error": None,
        "conserved_to_writer_tolerance": None,
    }
    if target is not None:
        errors = np.abs(represented - target)
        answer["max_absolute_error"] = float(np.max(errors))
        answer["max_relative_error"] = (
            float(np.max(errors) / abs(target)) if target != 0 else None
        )
        answer["conserved_to_writer_tolerance"] = bool(
            np.allclose(represented, target, rtol=1e-6, atol=1e-3)
        )
    return answer


def _require_compatible_machine(reference: ResultData, candidate: ResultData) -> None:
    """Require the minimal machine quantities needed for a loose comparison."""

    for key in ("beta0", "p0c_ev", "rf25_frequency_hz"):
        left = reference.machine_attrs.get(key)
        right = candidate.machine_attrs.get(key)
        if left is None or right is None:
            continue
        if not np.isclose(float(left), float(right), rtol=1e-12, atol=0.0):
            raise ValueError(
                f"Incompatible machine attribute {key}: reference={left}, candidate={right}"
            )


def _matching_value(left: Any, right: Any) -> bool:
    """Compare metadata scalars with tight tolerance for numeric roundoff."""

    if isinstance(left, (int, float, np.number)) and isinstance(
        right, (int, float, np.number)
    ):
        return bool(np.isclose(float(left), float(right), rtol=1e-12, atol=1e-15))
    return left == right


def _require_matching_attrs(
    name: str,
    reference: dict[str, Any],
    candidate: dict[str, Any],
    *,
    keys: tuple[str, ...] | None = None,
) -> None:
    """Require selected metadata dictionaries to contain matching values."""

    selected = sorted(set(reference) | set(candidate)) if keys is None else keys
    for key in selected:
        if key not in reference or key not in candidate:
            raise ValueError(f"Strict comparison requires matching {name} attribute {key}")
        if not _matching_value(reference[key], candidate[key]):
            raise ValueError(
                f"Strict comparison mismatch in {name} attribute {key}: "
                f"reference={reference[key]!r}, candidate={candidate[key]!r}"
            )


def require_strictly_compatible_results(
    reference: ResultData,
    candidate: ResultData,
) -> None:
    """Require the same map, turn/profile grid, and population scaling."""
    _require_matching_attrs("machine", reference.machine_attrs, candidate.machine_attrs)
    _require_matching_attrs(
        "rf_program", reference.rf_program_attrs, candidate.rf_program_attrs
    )
    _require_matching_attrs(
        "root",
        reference.root_attrs,
        candidate.root_attrs,
        keys=(
            "format_version",
            "model",
            "xtrack_version",
            "source_blond_commit",
            "source_blond_driver",
            "coordinate_convention",
            "energy_convention",
            "n_turns",
            "end_state_time_s",
            "last_kick_time_s",
            "first_rf25_kick_turn",
            "record_every",
            "sample_particles",
            "source_rows",
            "macro_equivalent_weight",
            "physical_proton_weight",
            "last_kick_rf53_v",
            "last_kick_rf25_v",
        ),
    )
    if not np.array_equal(reference.profiles.turns, candidate.profiles.turns):
        raise ValueError("Strict comparison requires identical recorded-turn grids")
    if not np.array_equal(
        reference.profiles.time_edges_ns, candidate.profiles.time_edges_ns
    ):
        raise ValueError("Strict comparison requires identical profile time-bin edges")


def _profile_comparison(reference: ResultData, candidate: ResultData) -> dict[str, Any]:
    """Compare normalized longitudinal profiles on their common recorded turns.

    Wasserstein distance retains the time-axis geometry inside the profile
    window.  Jensen-Shannon distance is reported both inside the window and
    with all outside-window mass represented as one additional category.
    """

    ref = reference.profiles
    cand = candidate.profiles
    if ref.time_edges_ns.shape != cand.time_edges_ns.shape or not np.allclose(
        ref.time_edges_ns, cand.time_edges_ns, rtol=0.0, atol=1e-12
    ):
        raise ValueError("Incompatible profile time-bin edges")

    # Runs may have different recording intervals, so compare only turns that
    # are present in both files rather than assuming row indices align.
    common, ref_index, cand_index = np.intersect1d(
        ref.turns, cand.turns, assume_unique=True, return_indices=True
    )
    if common.size == 0:
        raise ValueError("Result files have no common recorded turns")

    wasserstein: list[float | None] = []
    js_in_window: list[float | None] = []
    js_with_outside: list[float | None] = []
    reference_in_window: list[float | None] = []
    candidate_in_window: list[float | None] = []
    for left_i, right_i in zip(ref_index, cand_index):
        # Each distance normalizes its two rows independently.  Conservation and
        # in-window fractions below keep sample scaling and clipping visible.
        left = ref.counts_macro_equivalent[left_i]
        right = cand.counts_macro_equivalent[right_i]
        wasserstein.append(_histogram_wasserstein_ns(left, right, ref.time_edges_ns))
        js_in_window.append(_jensen_shannon_distance(left, right))

        # W1 has a physical time axis, so undefined outside positions cannot be
        # included.  For JS, outside mass can be one explicit categorical bin.
        left_with_outside = np.append(left, ref.outside_macro_equivalent[left_i])
        right_with_outside = np.append(right, cand.outside_macro_equivalent[right_i])
        js_with_outside.append(
            _jensen_shannon_distance(left_with_outside, right_with_outside)
        )
        left_total = float(np.sum(left_with_outside))
        right_total = float(np.sum(right_with_outside))
        reference_in_window.append(
            float(np.sum(left) / left_total) if left_total > 0 else None
        )
        candidate_in_window.append(
            float(np.sum(right) / right_total) if right_total > 0 else None
        )

    return {
        "common_turns": [int(turn) for turn in common],
        "n_common_turns": int(common.size),
        "wasserstein_ns": wasserstein,
        "jensen_shannon_distance": js_in_window,
        "jensen_shannon_with_outside_distance": js_with_outside,
        "reference_in_window_fraction": reference_in_window,
        "candidate_in_window_fraction": candidate_in_window,
        "wasserstein_summary_ns": _finite_summary(wasserstein),
        "jensen_shannon_summary": _finite_summary(js_in_window),
        "jensen_shannon_with_outside_summary": _finite_summary(js_with_outside),
        "normalization_note": (
            "Wasserstein and jensen_shannon_distance normalize each in-window "
            "profile row; the with_outside metric adds outside mass as one category."
        ),
    }


def _moments(values: np.ndarray) -> dict[str, float | None]:
    """Return the mean and population-standard-deviation beam width."""

    if values.size == 0:
        return {"mean": None, "rms": None}
    return {"mean": float(np.mean(values)), "rms": float(np.std(values, ddof=0))}


def _bucket_summary(
    final: FinalParticleData,
    machine_attrs: dict[str, Any],
) -> dict[str, Any]:
    """Summarize final survival and h=28 bucket-relative moments by lineage.

    Xsuite ``zeta`` is converted to the project's BLonD-compatible arrival-time
    convention. Each initial structural-bunch label is summarized around the
    h=28 bucket center nearest its own median and wrapped into one RF period.
    """

    missing = [
        key
        for key in ("beta0", "p0c_ev", "rf25_frequency_hz")
        if key not in machine_attrs
    ]
    if missing:
        return {
            "available": False,
            "reason": f"missing machine attributes: {', '.join(missing)}",
        }

    beta0 = float(machine_attrs["beta0"])
    p0c_ev = float(machine_attrs["p0c_ev"])
    rf25_frequency_hz = float(machine_attrs["rf25_frequency_hz"])
    if beta0 <= 0 or p0c_ev <= 0 or rf25_frequency_hz <= 0:
        return {"available": False, "reason": "invalid machine attributes"}

    # Convert Xsuite coordinates to the time and energy conventions used by the
    # Recycler analysis.  Positive ``state`` values identify surviving particles.
    time_s = -final.zeta_m / (beta0 * C_LIGHT_M_PER_S)
    energy_ev = final.ptau * p0c_ev
    finite = np.isfinite(time_s) & np.isfinite(energy_ev)
    alive = final.state > 0
    # Files without lineage labels can still produce a global summary; treating
    # them as label zero makes the existing per-bunch structure reusable.
    labels = (
        np.zeros(final.zeta_m.size, dtype=np.int64)
        if final.bunch_index is None
        else final.bunch_index
    )
    label_values = np.unique(labels)
    bucket_period_s = 1.0 / rf25_frequency_hz

    per_bunch: dict[str, Any] = {}
    for label in label_values:
        member = labels == label
        usable = member & alive & finite
        n_total = int(np.count_nonzero(member))
        n_alive = int(np.count_nonzero(member & alive))
        n_usable = int(np.count_nonzero(usable))
        entry: dict[str, Any] = {
            "n_total": n_total,
            "n_alive": n_alive,
            "n_alive_finite": n_usable,
            "survival_fraction": float(n_alive / n_total) if n_total else None,
            "finite_fraction_of_alive": float(n_usable / n_alive) if n_alive else None,
        }
        if n_usable:
            # Snap the median to the nearest h=28 bucket, then wrap offsets into
            # +/- half a period.  This describes bunch width, not bucket migration.
            bucket_center_s = (
                np.rint(np.median(time_s[usable]) / bucket_period_s) * bucket_period_s
            )
            relative_s = (
                (time_s[usable] - bucket_center_s + 0.5 * bucket_period_s)
                % bucket_period_s
                - 0.5 * bucket_period_s
            )
            time_moments = _moments(relative_s * 1e9)
            energy_moments = _moments(energy_ev[usable])
            entry.update(
                {
                    "bucket_center_ns": float(bucket_center_s * 1e9),
                    "mean_relative_time_ns": time_moments["mean"],
                    "rms_relative_time_ns": time_moments["rms"],
                    "mean_energy_offset_ev": energy_moments["mean"],
                    "rms_energy_offset_ev": energy_moments["rms"],
                }
            )
            for threshold_ns in (80, 100, 120):
                entry[f"fraction_abs_time_gt_{threshold_ns}_ns"] = float(
                    np.mean(np.abs(relative_s) > threshold_ns * 1e-9)
                )
        per_bunch[str(_plain_scalar(label))] = entry

    n_total = int(final.zeta_m.size)
    n_alive = int(np.count_nonzero(alive))
    return {
        "available": True,
        "n_total": n_total,
        "n_alive": n_alive,
        "n_alive_finite": int(np.count_nonzero(alive & finite)),
        "survival_fraction": float(n_alive / n_total) if n_total else None,
        "bucket_period_ns": float(bucket_period_s * 1e9),
        "per_bunch": per_bunch,
        "region_note": (
            "Bucket-relative thresholds use time wrapped around each initial "
            "bunch label's nearest h=28 bucket center."
        ),
    }


def _final_summary(result: ResultData) -> dict[str, Any]:
    """Return a final-particle summary or an explicit unavailable marker."""

    if result.final_particles is None:
        return {"available": False, "reason": "final_particles group is absent"}
    return _bucket_summary(result.final_particles, result.machine_attrs)


def _final_differences(reference: dict[str, Any], candidate: dict[str, Any]) -> dict[str, Any]:
    """Subtract matching per-lineage final metrics as candidate minus reference."""

    if not reference.get("available") or not candidate.get("available"):
        return {"available": False, "reason": "one or both final summaries are unavailable"}
    common_labels = sorted(set(reference["per_bunch"]) & set(candidate["per_bunch"]))
    fields = (
        "survival_fraction",
        "mean_relative_time_ns",
        "rms_relative_time_ns",
        "mean_energy_offset_ev",
        "rms_energy_offset_ev",
        "fraction_abs_time_gt_80_ns",
        "fraction_abs_time_gt_100_ns",
        "fraction_abs_time_gt_120_ns",
    )
    per_bunch: dict[str, Any] = {}
    for label in common_labels:
        left = reference["per_bunch"][label]
        right = candidate["per_bunch"][label]
        deltas: dict[str, float | None] = {}
        for field in fields:
            left_value = left.get(field)
            right_value = right.get(field)
            deltas[f"candidate_minus_reference_{field}"] = (
                float(right_value - left_value)
                if left_value is not None and right_value is not None
                else None
            )
        per_bunch[label] = deltas
    return {
        "available": True,
        "common_bunch_labels": common_labels,
        "per_bunch": per_bunch,
    }


def compare_results(
    reference_path: str | Path,
    candidate_path: str | Path,
    *,
    require_matching_run: bool = False,
) -> dict[str, Any]:
    """Compare profiles, population conservation, and final bunch properties.

    By default only the machine quantities required to interpret the coordinates
    must match.  ``require_matching_run=True`` additionally enforces identical
    maps, RF programs, endpoints, profile grids, and population scaling.
    Returned values are safe to serialize directly as JSON.
    """

    reference = load_result_data(reference_path)
    candidate = load_result_data(candidate_path)
    _require_compatible_machine(reference, candidate)
    if require_matching_run:
        require_strictly_compatible_results(reference, candidate)

    reference_final = _final_summary(reference)
    candidate_final = _final_summary(candidate)
    return {
        "reference_path": str(reference.path),
        "candidate_path": str(candidate.path),
        "run_metadata": {
            "reference_input_source_kind": reference.root_attrs.get(
                "input_source_kind"
            ),
            "candidate_input_source_kind": candidate.root_attrs.get(
                "input_source_kind"
            ),
            "strict_matching_run_required": require_matching_run,
        },
        "profiles": _profile_comparison(reference, candidate),
        "conservation": {
            "reference": _conservation(reference),
            "candidate": _conservation(candidate),
        },
        "final_particles": {
            "reference": reference_final,
            "candidate": candidate_final,
            "differences": _final_differences(reference_final, candidate_final),
        },
    }
