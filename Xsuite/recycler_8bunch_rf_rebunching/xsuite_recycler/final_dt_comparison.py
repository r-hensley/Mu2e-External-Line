r"""Strict final-:math:`\Delta t` comparisons for RF-only result files.

The comparison operates on surviving particles with finite final ``zeta`` and
uses the BLonD-compatible convention already recorded by this project,

``dt = -zeta / (beta0 * c)``.

Histogram distances are computed after independently normalizing the two
samples.  Raw, macro-equivalent, physical-proton, and Kish-effective counts are
reported separately so unequal tracking populations cannot be mistaken for a
shape discrepancy.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, replace
from pathlib import Path
from tempfile import gettempdir
from typing import Any

import numpy as np

from .config import C_LIGHT_M_PER_S
from .result_comparison import (
    ResultData,
    load_result_data,
    require_strictly_compatible_results,
)


_TAIL_THRESHOLDS_NS = (80, 100, 120)
_LINEAGE_WINDOW_HALF_WIDTH_NS = 125.0
_SAMPLING_ROOT_ATTRS = (
    "sample_particles",
    "macro_equivalent_weight",
    "physical_proton_weight",
)
_GAMMA_TRANSITION_MACHINE_ATTRS = (
    "gamma_transition",
    "momentum_compaction_factor",
    "slip_factor",
)
_SUMMARY_FIELDS = (
    "mean_ns",
    "rms_ns",
    "minimum_ns",
    "q01_ns",
    "q05_ns",
    "q10_ns",
    "q25_ns",
    "median_ns",
    "q75_ns",
    "q90_ns",
    "q95_ns",
    "q99_ns",
    "maximum_ns",
)


@dataclass(frozen=True)
class FinalDtData:
    """Validated final times, lineage labels, and constant particle weights.

    Only particles with positive Xsuite state and finite ``zeta`` contribute to
    ``dt_ns``.  Stored and surviving counts remain available for diagnostics.
    """

    result: ResultData
    dt_ns: np.ndarray
    bunch_index: np.ndarray
    macro_equivalent_weight: float
    physical_proton_weight: float
    n_stored: int
    n_alive: int
    n_alive_finite_dt: int


def _positive_finite_attr(result: ResultData, key: str) -> float:
    """Read one required positive finite numeric root attribute."""

    if key not in result.root_attrs:
        raise ValueError(f"{result.path}: missing required root attribute {key}")
    value = float(result.root_attrs[key])
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{result.path}: {key} must be positive and finite")
    return value


def _validate_sampling_metadata(result: ResultData) -> tuple[float, float]:
    """Validate sample length and return macro and physical particle weights."""

    final = result.final_particles
    if final is None:
        raise ValueError(f"{result.path}: final_particles group is required")

    sample_particles = _positive_finite_attr(result, "sample_particles")
    if not sample_particles.is_integer():
        raise ValueError(f"{result.path}: sample_particles must be an integer")
    if int(sample_particles) != final.zeta_m.size:
        raise ValueError(
            f"{result.path}: sample_particles does not match final-particle length"
        )

    macro_weight = _positive_finite_attr(result, "macro_equivalent_weight")
    physical_weight = _positive_finite_attr(result, "physical_proton_weight")
    return macro_weight, physical_weight


def _strict_physics_compatibility(
    reference: ResultData,
    candidate: ResultData,
    *,
    allow_gamma_transition_difference: bool = False,
) -> None:
    """Validate an ordinary comparison or a controlled gamma-transition study.

    ``require_strictly_compatible_results`` deliberately treats sampling and
    per-particle weights as part of a matching run.  Final-distribution shape
    comparisons also need to support, for example, a 100k reference against a
    33.6k candidate.  The three sampling-only attributes are therefore
    validated independently and normalized in a temporary immutable view.  A
    controlled gamma-transition study instead requires the same full file input
    and permits only gamma transition plus its derived alpha_c and eta to differ.
    """

    _validate_sampling_metadata(reference)
    _validate_sampling_metadata(candidate)

    if allow_gamma_transition_difference:
        for key in _GAMMA_TRANSITION_MACHINE_ATTRS:
            if key not in reference.machine_attrs or key not in candidate.machine_attrs:
                raise ValueError(
                    "Transition-gamma comparison requires machine attribute "
                    f"{key} in both files"
                )
        for key in ("input_source_kind", "input_sha256"):
            left = reference.root_attrs.get(key)
            right = candidate.root_attrs.get(key)
            if left is None or right is None or left != right:
                raise ValueError(
                    "Transition-gamma comparison requires matching root attribute "
                    f"{key}: reference={left!r}, candidate={right!r}"
                )
        if reference.root_attrs["input_source_kind"] != "file":
            raise ValueError(
                "Transition-gamma comparison requires the same file-backed input"
            )
        for result in (reference, candidate):
            sample = int(result.root_attrs["sample_particles"])
            source = int(result.root_attrs["source_rows"])
            if sample != source:
                raise ValueError(
                    f"{result.path}: transition-gamma comparison requires all "
                    f"source rows, but sample_particles={sample} and source_rows={source}"
                )
            final = result.final_particles
            if final is None:  # Already checked by sampling validation.
                raise ValueError(f"{result.path}: final_particles group is required")
            usable = (final.state > 0) & np.isfinite(final.zeta_m)
            if not np.all(usable):
                raise ValueError(
                    f"{result.path}: transition-gamma pairing requires every "
                    "particle to survive with finite final zeta"
                )
        reference_labels = reference.final_particles.bunch_index
        candidate_labels = candidate.final_particles.bunch_index
        if reference_labels is None or candidate_labels is None:
            raise ValueError(
                "Transition-gamma comparison requires final bunch_index arrays"
            )
        if not np.array_equal(reference_labels, candidate_labels):
            raise ValueError(
                "Transition-gamma comparison requires rowwise matching bunch_index"
            )

        # Normalize only the three declared study parameters in a temporary
        # view. The strict validator then catches every unintended difference.
        normalized_machine = dict(candidate.machine_attrs)
        for key in _GAMMA_TRANSITION_MACHINE_ATTRS:
            normalized_machine[key] = reference.machine_attrs[key]
        require_strictly_compatible_results(
            reference,
            replace(candidate, machine_attrs=normalized_machine),
        )
        return

    normalized_attrs = dict(candidate.root_attrs)
    for key in _SAMPLING_ROOT_ATTRS:
        normalized_attrs[key] = reference.root_attrs[key]
    normalized_candidate = replace(candidate, root_attrs=normalized_attrs)
    require_strictly_compatible_results(reference, normalized_candidate)


def _integer_labels(path: Path, labels: np.ndarray | None) -> np.ndarray:
    """Validate required initial structural-bunch labels as integers."""

    if labels is None:
        raise ValueError(f"{path}: final_particles/bunch_index is required")
    if labels.size == 0:
        return np.asarray(labels, dtype=np.int64)
    if not np.issubdtype(labels.dtype, np.number):
        raise ValueError(f"{path}: bunch labels must be numeric integers")
    numeric = np.asarray(labels, dtype=float)
    if not np.all(np.isfinite(numeric)) or not np.all(numeric == np.rint(numeric)):
        raise ValueError(f"{path}: bunch labels must be finite integers")
    if np.any(numeric < 0):
        raise ValueError(f"{path}: bunch labels must be nonnegative")
    return numeric.astype(np.int64)


def _extract_final_dt(result: ResultData) -> FinalDtData:
    """Select usable survivors and convert Xsuite ``zeta`` to final arrival time.

    The conversion follows the BLonD-compatible project convention
    ``dt = -zeta / (beta0 * c)``.  Non-survivors and non-finite longitudinal
    coordinates are counted for diagnostics but excluded from distributions.
    """

    final = result.final_particles
    if final is None:  # Kept explicit for callers of this private helper.
        raise ValueError(f"{result.path}: final_particles group is required")
    if "beta0" not in result.machine_attrs:
        raise ValueError(f"{result.path}: missing machine attribute beta0")
    beta0 = float(result.machine_attrs["beta0"])
    if not np.isfinite(beta0) or beta0 <= 0:
        raise ValueError(f"{result.path}: beta0 must be positive and finite")

    macro_weight, physical_weight = _validate_sampling_metadata(result)
    labels = _integer_labels(result.path, final.bunch_index)
    if not np.all(np.isfinite(final.state)):
        raise ValueError(f"{result.path}: final particle state must be finite")
    # In Xsuite, a positive state denotes an active particle.  Keep the finite
    # ``zeta`` requirement separate so excluded survivors remain countable.
    alive = final.state > 0
    finite_dt = np.isfinite(final.zeta_m)
    usable = alive & finite_dt
    if not np.any(usable):
        raise ValueError(f"{result.path}: no surviving particles have finite final zeta")

    # Xsuite zeta is a longitudinal distance.  Dividing by beta*c converts it
    # to time, and the minus sign matches the historical BLonD convention.
    dt_ns = -np.asarray(final.zeta_m[usable], dtype=float) / (
        beta0 * C_LIGHT_M_PER_S
    ) * 1e9
    return FinalDtData(
        result=result,
        dt_ns=dt_ns,
        bunch_index=labels[usable],
        macro_equivalent_weight=macro_weight,
        physical_proton_weight=physical_weight,
        n_stored=int(final.zeta_m.size),
        n_alive=int(np.count_nonzero(alive)),
        n_alive_finite_dt=int(np.count_nonzero(usable)),
    )


def load_final_dt_data(path: str | Path) -> FinalDtData:
    """Load one result through the project schema validator and extract dt."""

    return _extract_final_dt(load_result_data(path))


def _validate_histogram_spec(
    bins: int,
    value_range: tuple[float, float],
    name: str,
) -> np.ndarray:
    """Validate a histogram bin count and range, then return uniform edges."""

    if isinstance(bins, bool) or int(bins) != bins or bins < 1:
        raise ValueError(f"{name} bins must be a positive integer")
    lower, upper = (float(value_range[0]), float(value_range[1]))
    if not np.isfinite(lower) or not np.isfinite(upper) or lower >= upper:
        raise ValueError(f"{name} range must contain finite increasing limits")
    return np.linspace(lower, upper, int(bins) + 1)


def _histogram(values: np.ndarray, weight: float, edges: np.ndarray) -> np.ndarray:
    """Histogram samples using one constant macro-equivalent particle weight."""

    return np.histogram(
        values,
        bins=edges,
        weights=np.full(values.size, weight, dtype=float),
    )[0].astype(float, copy=False)


def _normalize(masses: np.ndarray) -> np.ndarray | None:
    """Normalize nonzero histogram masses to one, or return ``None`` if empty."""

    total = float(np.sum(masses, dtype=np.float64))
    return None if total <= 0 else masses / total


def _histogram_distances(
    reference_values: np.ndarray,
    candidate_values: np.ndarray,
    reference_weight: float,
    candidate_weight: float,
    edges: np.ndarray,
) -> dict[str, float | int | None]:
    """Compare two weighted histograms using four normalized shape distances.

    Wasserstein-1, Jensen-Shannon, histogram KS, and total variation are based
    on independently normalized in-range mass.  The returned in-range fractions
    separately expose particles clipped by the requested plotting range.
    """

    # Constant per-particle weights recover macro-equivalent counts.  They
    # cancel from each sample's normalized shape but remain relevant to the
    # independently reported fraction captured by these edges.
    left_mass = _histogram(reference_values, reference_weight, edges)
    right_mass = _histogram(candidate_values, candidate_weight, edges)
    left = _normalize(left_mass)
    right = _normalize(right_mass)
    left_in = float(np.sum(left_mass, dtype=np.float64))
    right_in = float(np.sum(right_mass, dtype=np.float64))
    left_total = float(reference_values.size * reference_weight)
    right_total = float(candidate_values.size * candidate_weight)
    answer: dict[str, float | int | None] = {
        "n_bins": int(edges.size - 1),
        "range_lower_ns": float(edges[0]),
        "range_upper_ns": float(edges[-1]),
        "reference_in_range_fraction": left_in / left_total if left_total else None,
        "candidate_in_range_fraction": right_in / right_total if right_total else None,
        "wasserstein_1_ns": None,
        "jensen_shannon_distance": None,
        "histogram_ks_distance": None,
        "total_variation_distance": None,
    }
    if left is None or right is None:
        return answer

    # These are discrete approximations on common bin centers.  The empirical
    # W1 reported elsewhere avoids this histogram discretization.
    centers = 0.5 * (edges[:-1] + edges[1:])
    if centers.size == 1:
        wasserstein = 0.0
    else:
        wasserstein = float(
            np.sum(np.abs(np.cumsum(left - right)[:-1]) * np.diff(centers))
        )
    midpoint = 0.5 * (left + right)
    left_mask = left > 0
    right_mask = right > 0
    divergence = 0.5 * np.sum(
        left[left_mask] * np.log2(left[left_mask] / midpoint[left_mask])
    )
    divergence += 0.5 * np.sum(
        right[right_mask] * np.log2(right[right_mask] / midpoint[right_mask])
    )
    answer["wasserstein_1_ns"] = wasserstein
    answer["jensen_shannon_distance"] = float(
        np.sqrt(max(float(divergence), 0.0))
    )
    answer["histogram_ks_distance"] = float(
        np.max(np.abs(np.cumsum(left - right)))
    )
    answer["total_variation_distance"] = float(0.5 * np.sum(np.abs(left - right)))
    return answer


def _empirical_wasserstein_1_ns(
    reference_values: np.ndarray,
    candidate_values: np.ndarray,
) -> float | None:
    """Exact 1D W1 for independently normalized empirical samples.

    Each result has one constant per-particle weight, which cancels when its
    empirical distribution is normalized.  Unequal sample sizes are retained
    through the distinct ``1/n`` CDF steps.
    """

    if reference_values.size == 0 or candidate_values.size == 0:
        return None
    left = np.sort(np.asarray(reference_values, dtype=float))
    right = np.sort(np.asarray(candidate_values, dtype=float))
    support = np.sort(np.concatenate((left, right)))
    widths = np.diff(support)
    if widths.size == 0:
        return 0.0
    left_cdf = np.searchsorted(left, support[:-1], side="right") / left.size
    right_cdf = np.searchsorted(right, support[:-1], side="right") / right.size
    return float(np.sum(np.abs(left_cdf - right_cdf) * widths))


def _summary(values_ns: np.ndarray) -> dict[str, float | int | None]:
    """Return count, quantiles, and population-standard-deviation beam width."""

    if values_ns.size == 0:
        return {"n": 0, **{field: None for field in _SUMMARY_FIELDS}}
    quantiles = np.quantile(
        values_ns,
        (0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99),
    )
    return {
        "n": int(values_ns.size),
        "mean_ns": float(np.mean(values_ns)),
        "rms_ns": float(np.std(values_ns, ddof=0)),
        "minimum_ns": float(np.min(values_ns)),
        "q01_ns": float(quantiles[0]),
        "q05_ns": float(quantiles[1]),
        "q10_ns": float(quantiles[2]),
        "q25_ns": float(quantiles[3]),
        "median_ns": float(quantiles[4]),
        "q75_ns": float(quantiles[5]),
        "q90_ns": float(quantiles[6]),
        "q95_ns": float(quantiles[7]),
        "q99_ns": float(quantiles[8]),
        "maximum_ns": float(np.max(values_ns)),
    }


def _summary_differences(
    reference: dict[str, float | int | None],
    candidate: dict[str, float | int | None],
) -> dict[str, dict[str, float | None]]:
    """Return signed and absolute candidate-minus-reference summary changes."""

    signed: dict[str, float | None] = {}
    absolute: dict[str, float | None] = {}
    for field in _SUMMARY_FIELDS:
        left = reference.get(field)
        right = candidate.get(field)
        difference = (
            float(right - left) if left is not None and right is not None else None
        )
        signed[field] = difference
        absolute[field] = abs(difference) if difference is not None else None
    return {
        "candidate_minus_reference": signed,
        "absolute_difference": absolute,
    }


def _nearest_bucket_center_ns(values_ns: np.ndarray, period_ns: float) -> float | None:
    """Snap a sample median to the nearest h=28 RF-bucket center."""

    if values_ns.size == 0:
        return None
    return float(np.rint(np.median(values_ns) / period_ns) * period_ns)


def _reference_bucket_centers(
    reference: FinalDtData,
    period_ns: float,
) -> dict[int, float]:
    """Define one historical-reference h=28 center for each initial lineage.

    These same centers are later applied to both samples, so a candidate timing
    shift cannot disappear through independent recentering.
    """

    centers: dict[int, float] = {}
    for label in np.unique(reference.bunch_index):
        values = reference.dt_ns[reference.bunch_index == label]
        center = _nearest_bucket_center_ns(values, period_ns)
        if center is not None:
            centers[int(label)] = center
    return centers


def _relative_to_reference_centers(
    data: FinalDtData,
    centers_ns: dict[int, float],
) -> np.ndarray:
    """Subtract shared reference centers without wrapping into one RF period."""

    relative = np.empty_like(data.dt_ns)
    for label in np.unique(data.bunch_index):
        key = int(label)
        if key not in centers_ns:
            raise ValueError(
                f"{data.result.path}: surviving bunch label {key} has no reference "
                "h=28 bucket center"
            )
        selected = data.bunch_index == label
        # Deliberately do not wrap modulo the h=28 period: adjacent-bucket and
        # satellite migration must remain visible in the comparison and tails.
        relative[selected] = data.dt_ns[selected] - centers_ns[key]
    return relative


def _bucket_occupancy(
    values_ns: np.ndarray,
    period_ns: float,
    macro_weight: float,
    physical_weight: float,
) -> dict[str, dict[str, float | int]]:
    """Report weighted occupancy of the nearest final h=28 buckets."""

    bucket_ids = np.rint(values_ns / period_ns).astype(np.int64)
    answer: dict[str, dict[str, float | int]] = {}
    for bucket_id in np.unique(bucket_ids):
        count = int(np.count_nonzero(bucket_ids == bucket_id))
        answer[str(int(bucket_id))] = {
            "raw_count": count,
            "fraction": float(count / values_ns.size) if values_ns.size else 0.0,
            "macro_equivalent_count": float(count * macro_weight),
            "physical_proton_count": float(count * physical_weight),
        }
    return answer


def _tail_summary(
    relative_ns: np.ndarray,
    macro_weight: float,
    physical_weight: float,
) -> dict[str, dict[str, float | int]]:
    """Count particles outside each configured reference-relative time region."""

    answer: dict[str, dict[str, float | int]] = {}
    for threshold in _TAIL_THRESHOLDS_NS:
        raw = int(np.count_nonzero(np.abs(relative_ns) > threshold))
        answer[f"abs_dt_gt_{threshold}_ns"] = {
            "raw_count": raw,
            "fraction": float(raw / relative_ns.size) if relative_ns.size else 0.0,
            "macro_equivalent_count": float(raw * macro_weight),
            "physical_proton_count": float(raw * physical_weight),
        }
    return answer


def _outside_lineage_window_counts(
    relative_ns: np.ndarray,
    macro_weight: float,
    physical_weight: float,
) -> dict[str, float | int]:
    """Count particles strictly outside a supplied +/-125 ns coordinate.

    ``relative_ns`` may be either an unwrapped initial-lineage coordinate or a
    nearest-bucket folded local coordinate; the calling diagnostic records that
    convention. Equality at either boundary is in-time because the rule is
    ``abs(dt) > 125 ns``. The three fractions are numerically equal for a result
    with constant particle weights, but reporting them explicitly keeps every
    count's denominator unambiguous for downstream consumers.
    """

    raw_total = int(relative_ns.size)
    raw_outside = int(
        np.count_nonzero(np.abs(relative_ns) > _LINEAGE_WINDOW_HALF_WIDTH_NS)
    )
    fraction = float(raw_outside / raw_total) if raw_total else 0.0
    return {
        "raw_total_count": raw_total,
        "raw_outside_count": raw_outside,
        "raw_outside_fraction": fraction,
        "macro_equivalent_total_count": float(raw_total * macro_weight),
        "macro_equivalent_outside_count": float(raw_outside * macro_weight),
        "macro_equivalent_outside_fraction": fraction,
        "physical_proton_total_count": float(raw_total * physical_weight),
        "physical_proton_outside_count": float(raw_outside * physical_weight),
        "physical_proton_outside_fraction": fraction,
    }


def _unwrapped_lineage_window_diagnostic(
    data: FinalDtData,
    reference_centers_ns: dict[int, float],
) -> dict[str, Any]:
    """Build an unwrapped lineage-window plus migration diagnostic.

    A particle is outside when it lies more than 125 ns from the historical
    h=28 center assigned to its initial lineage.  Therefore an otherwise-local
    particle that migrates to an adjacent RF bucket is also outside.  This is a
    useful conservative combined diagnostic, not Werkema's local pulse window
    and not an extinction prediction.
    """

    relative = _relative_to_reference_centers(data, reference_centers_ns)
    per_lineage: dict[str, Any] = {}
    for label in np.unique(data.bunch_index):
        selected = data.bunch_index == label
        per_lineage[str(int(label))] = {
            "reference_h28_bucket_center_ns": reference_centers_ns[int(label)],
            **_outside_lineage_window_counts(
                relative[selected],
                data.macro_equivalent_weight,
                data.physical_proton_weight,
            ),
        }
    return {
        "window_half_width_ns": _LINEAGE_WINDOW_HALF_WIDTH_NS,
        "outside_rule": "abs(dt_rel_ns) > window_half_width_ns",
        "dt_rel_definition": (
            "dt_ns - historical_reference_h28_center_ns_for_initial_lineage"
        ),
        "survivor_population": "state>0 and finite(zeta_m)",
        "modulo_wrapping": False,
        "separate_from_reference_bucket_migration": True,
        "interpretation": (
            "combined local-tail-or-reference-bucket-migration diagnostic"
        ),
        "model_scope_limitation": (
            "Recycler endpoint only; not a Delivery Ring, extraction-gate, "
            "or proton-target extinction prediction"
        ),
        "reference_bucket_migration_metric_path": (
            "initial_lineage_reference_bucket_retention"
        ),
        "global": _outside_lineage_window_counts(
            relative,
            data.macro_equivalent_weight,
            data.physical_proton_weight,
        ),
        "per_initial_lineage": per_lineage,
    }


def _fold_to_nearest_h28_center(values_ns: np.ndarray, period_ns: float) -> np.ndarray:
    """Fold absolute arrival times onto the shared nearest h=28 center grid."""

    return (values_ns + 0.5 * period_ns) % period_ns - 0.5 * period_ns


def _folded_local_window_diagnostic(
    data: FinalDtData,
    period_ns: float,
) -> dict[str, Any]:
    """Build the nearest-bucket folded +/-125 ns local-pulse diagnostic.

    This is the Recycler-endpoint quantity aligned with Werkema's local RF-pulse
    window: each particle is first referenced to its nearest h=28 center on the
    common zero-phase grid. Initial lineage labels only stratify the result and
    do not choose the center. The metric is not a downstream extinction model.
    """

    local = _fold_to_nearest_h28_center(data.dt_ns, period_ns)
    per_lineage: dict[str, Any] = {}
    for label in np.unique(data.bunch_index):
        selected = data.bunch_index == label
        per_lineage[str(int(label))] = _outside_lineage_window_counts(
            local[selected],
            data.macro_equivalent_weight,
            data.physical_proton_weight,
        )
    return {
        "window_half_width_ns": _LINEAGE_WINDOW_HALF_WIDTH_NS,
        "outside_rule": "abs(dt_local_ns) > window_half_width_ns",
        "dt_local_definition": (
            "((dt_ns + 0.5*h28_bucket_period_ns) modulo "
            "h28_bucket_period_ns) - 0.5*h28_bucket_period_ns"
        ),
        "h28_bucket_period_ns": float(period_ns),
        "bucket_center_grid_origin_ns": 0.0,
        "nearest_bucket_center_rule": (
            "half-open modulo grid [-h28_bucket_period_ns/2, "
            "h28_bucket_period_ns/2)"
        ),
        "modulo_wrapping": True,
        "survivor_population": "state>0 and finite(zeta_m)",
        "initial_lineage_used_for_center": False,
        "per_initial_lineage_is_stratification_only": True,
        "werkema_alignment": (
            "Recycler-endpoint proxy for the local +/-125 ns RF-pulse window"
        ),
        "model_scope_limitation": (
            "not a Delivery Ring, extraction-gate, or proton-target extinction "
            "prediction"
        ),
        "separate_from_unwrapped_initial_lineage_window": True,
        "separate_from_reference_bucket_migration": True,
        "global": _outside_lineage_window_counts(
            local,
            data.macro_equivalent_weight,
            data.physical_proton_weight,
        ),
        "per_initial_lineage": per_lineage,
    }


def _effective_counts(data: FinalDtData) -> dict[str, float | int]:
    """Report raw, weighted, and Kish-effective survivor populations.

    Every particle in one result has the same weight, so the Kish effective
    sample size equals the raw usable sample size.
    """

    usable = data.n_alive_finite_dt
    return {
        "stored_particles": data.n_stored,
        "surviving_particles": data.n_alive,
        "surviving_finite_dt_particles": usable,
        "excluded_alive_nonfinite_dt_particles": data.n_alive - usable,
        "raw_effective_sample_size": usable,
        "kish_effective_sample_size": float(usable),
        "macro_equivalent_weight_per_particle": data.macro_equivalent_weight,
        "physical_proton_weight_per_particle": data.physical_proton_weight,
        "surviving_finite_macro_equivalent_count": float(
            usable * data.macro_equivalent_weight
        ),
        "surviving_finite_physical_proton_count": float(
            usable * data.physical_proton_weight
        ),
    }


def _reference_bucket_retention(
    per_bunch: dict[str, Any],
    macro_weight: float,
    physical_weight: float,
) -> dict[str, float | int]:
    """Count each initial lineage inside versus outside its reference bucket."""
    total = sum(
        int(entry["effective_counts"]["raw_effective_sample_size"])
        for entry in per_bunch.values()
    )
    retained = 0
    for entry in per_bunch.values():
        reference_bucket = str(entry["reference_h28_bucket_id"])
        retained += int(
            entry["final_h28_bucket_occupancy"]
            .get(reference_bucket, {})
            .get("raw_count", 0)
        )
    outside = total - retained
    return {
        "raw_total": total,
        "raw_in_reference_bucket": retained,
        "raw_outside_reference_bucket": outside,
        "fraction_in_reference_bucket": float(retained / total) if total else 0.0,
        "fraction_outside_reference_bucket": float(outside / total) if total else 0.0,
        "outside_reference_bucket_macro_equivalent_count": float(
            outside * macro_weight
        ),
        "outside_reference_bucket_physical_proton_count": float(
            outside * physical_weight
        ),
    }


def _distribution_summary(
    data: FinalDtData,
    period_ns: float,
    reference_centers_ns: dict[int, float],
) -> dict[str, Any]:
    """Build aggregate and per-lineage summaries in absolute and shared frames."""

    relative = _relative_to_reference_centers(data, reference_centers_ns)
    per_bunch: dict[str, Any] = {}
    for label in np.unique(data.bunch_index):
        # ``bunch_index`` follows a particle's initial structural bunch; final
        # bucket occupancy may therefore reveal migration or satellite structure.
        selected = data.bunch_index == label
        absolute = data.dt_ns[selected]
        reference_center = reference_centers_ns[int(label)]
        bucket_relative = absolute - reference_center
        per_bunch[str(int(label))] = {
            "effective_counts": {
                "raw_effective_sample_size": int(absolute.size),
                "kish_effective_sample_size": float(absolute.size),
                "macro_equivalent_count": float(
                    absolute.size * data.macro_equivalent_weight
                ),
                "physical_proton_count": float(
                    absolute.size * data.physical_proton_weight
                ),
            },
            "absolute_dt": _summary(absolute),
            "reference_h28_bucket_center_ns": reference_center,
            "reference_h28_bucket_id": int(np.rint(reference_center / period_ns)),
            "sample_median_nearest_h28_bucket_center_ns": (
                _nearest_bucket_center_ns(absolute, period_ns)
            ),
            "bucket_relative_dt": _summary(bucket_relative),
            "bucket_relative_tails": _tail_summary(
                bucket_relative,
                data.macro_equivalent_weight,
                data.physical_proton_weight,
            ),
            "final_h28_bucket_occupancy": _bucket_occupancy(
                absolute,
                period_ns,
                data.macro_equivalent_weight,
                data.physical_proton_weight,
            ),
        }
    return {
        "effective_counts": _effective_counts(data),
        "absolute_dt": _summary(data.dt_ns),
        "bucket_relative_dt": _summary(relative),
        "bucket_relative_tails": _tail_summary(
            relative, data.macro_equivalent_weight, data.physical_proton_weight
        ),
        "final_h28_bucket_occupancy": _bucket_occupancy(
            data.dt_ns,
            period_ns,
            data.macro_equivalent_weight,
            data.physical_proton_weight,
        ),
        "initial_lineage_reference_bucket_retention": _reference_bucket_retention(
            per_bunch,
            data.macro_equivalent_weight,
            data.physical_proton_weight,
        ),
        "unwrapped_initial_lineage_125ns_window": (
            _unwrapped_lineage_window_diagnostic(data, reference_centers_ns)
        ),
        "nearest_h28_bucket_folded_125ns_local_window": (
            _folded_local_window_diagnostic(data, period_ns)
        ),
        "per_bunch": per_bunch,
    }


def _occupancy_differences(
    reference: dict[str, dict[str, float | int]],
    candidate: dict[str, dict[str, float | int]],
) -> dict[str, dict[str, float | int]]:
    """Compare final RF-bucket occupancy over the union of occupied buckets."""

    answer: dict[str, dict[str, float | int]] = {}
    for bucket_id in sorted(set(reference) | set(candidate), key=int):
        left = reference.get(bucket_id)
        right = candidate.get(bucket_id)
        left_raw = int(left["raw_count"]) if left is not None else 0
        right_raw = int(right["raw_count"]) if right is not None else 0
        left_fraction = float(left["fraction"]) if left is not None else 0.0
        right_fraction = float(right["fraction"]) if right is not None else 0.0
        left_macro = (
            float(left["macro_equivalent_count"]) if left is not None else 0.0
        )
        right_macro = (
            float(right["macro_equivalent_count"]) if right is not None else 0.0
        )
        answer[bucket_id] = {
            "reference_raw_count": left_raw,
            "candidate_raw_count": right_raw,
            "candidate_minus_reference_raw_count": right_raw - left_raw,
            "candidate_minus_reference_fraction": right_fraction - left_fraction,
            "candidate_minus_reference_macro_equivalent_count": (
                right_macro - left_macro
            ),
        }
    return answer


def _tail_differences(
    reference: dict[str, dict[str, float | int]],
    candidate: dict[str, dict[str, float | int]],
) -> dict[str, dict[str, float]]:
    """Compare matching bucket-relative tail regions between two summaries."""

    answer: dict[str, dict[str, float]] = {}
    for region in sorted(set(reference) & set(candidate)):
        answer[region] = {
            "candidate_minus_reference_raw_count": int(
                candidate[region]["raw_count"] - reference[region]["raw_count"]
            ),
            "candidate_minus_reference_fraction": float(
                candidate[region]["fraction"] - reference[region]["fraction"]
            ),
            "candidate_minus_reference_macro_equivalent_count": float(
                candidate[region]["macro_equivalent_count"]
                - reference[region]["macro_equivalent_count"]
            ),
            "absolute_fraction_difference": float(
                abs(candidate[region]["fraction"] - reference[region]["fraction"])
            ),
        }
    return answer


def _lineage_window_count_differences(
    reference: dict[str, float | int],
    candidate: dict[str, float | int],
) -> dict[str, float | int]:
    """Compare one pair of unwrapped 125 ns count-and-fraction summaries."""

    answer: dict[str, float | int] = {}
    for count_kind in ("raw", "macro_equivalent", "physical_proton"):
        total_key = f"{count_kind}_total_count"
        outside_key = f"{count_kind}_outside_count"
        fraction_key = f"{count_kind}_outside_fraction"
        total_difference = candidate[total_key] - reference[total_key]
        outside_difference = candidate[outside_key] - reference[outside_key]
        fraction_difference = candidate[fraction_key] - reference[fraction_key]
        # Preserve integer semantics for raw particle counts; weighted counts
        # and all fractions remain floating point values.
        answer[f"candidate_minus_reference_{total_key}"] = (
            int(total_difference)
            if count_kind == "raw"
            else float(total_difference)
        )
        answer[f"candidate_minus_reference_{outside_key}"] = (
            int(outside_difference)
            if count_kind == "raw"
            else float(outside_difference)
        )
        answer[f"candidate_minus_reference_{fraction_key}"] = float(
            fraction_difference
        )
        answer[f"absolute_{count_kind}_outside_fraction_difference"] = float(
            abs(fraction_difference)
        )
    return answer


def _lineage_window_differences(
    reference: dict[str, Any],
    candidate: dict[str, Any],
) -> dict[str, Any]:
    """Compare global and per-lineage 125 ns window diagnostics.

    Missing lineages are listed rather than silently compared with zero.  This
    mirrors the report's existing bunch-label handling and keeps the physical
    meaning of each initial lineage explicit.
    """

    reference_per_lineage = reference["per_initial_lineage"]
    candidate_per_lineage = candidate["per_initial_lineage"]
    common = sorted(
        set(reference_per_lineage) & set(candidate_per_lineage), key=int
    )
    return {
        "global": _lineage_window_count_differences(
            reference["global"], candidate["global"]
        ),
        "common_initial_lineage_labels": common,
        "reference_only_initial_lineage_labels": sorted(
            set(reference_per_lineage) - set(candidate_per_lineage), key=int
        ),
        "candidate_only_initial_lineage_labels": sorted(
            set(candidate_per_lineage) - set(reference_per_lineage), key=int
        ),
        "per_initial_lineage": {
            label: _lineage_window_count_differences(
                reference_per_lineage[label], candidate_per_lineage[label]
            )
            for label in common
        },
    }


def _paired_window_transition_block(
    reference_outside: np.ndarray,
    candidate_outside: np.ndarray,
    macro_weight: float,
    physical_weight: float,
) -> dict[str, Any]:
    """Count paired inside/outside transitions for one particle selection."""

    paired_total = int(reference_outside.size)
    raw = {
        "paired_total": paired_total,
        "inside_to_inside": int(
            np.count_nonzero(~reference_outside & ~candidate_outside)
        ),
        "inside_to_outside": int(
            np.count_nonzero(~reference_outside & candidate_outside)
        ),
        "outside_to_inside": int(
            np.count_nonzero(reference_outside & ~candidate_outside)
        ),
        "outside_to_outside": int(
            np.count_nonzero(reference_outside & candidate_outside)
        ),
    }
    raw["changed"] = raw["inside_to_outside"] + raw["outside_to_inside"]
    raw["net_outside"] = raw["inside_to_outside"] - raw["outside_to_inside"]

    def scaled_counts(weight: float) -> dict[str, float]:
        """Scale raw paired counts by one constant particle weight."""

        return {key: float(value * weight) for key, value in raw.items()}

    fractions = {
        key: float(value / paired_total) if paired_total else 0.0
        for key, value in raw.items()
    }
    return {
        "counts": {
            "raw": raw,
            "macro_equivalent": scaled_counts(macro_weight),
            "physical_proton": scaled_counts(physical_weight),
        },
        "fractions_of_paired": {
            "raw": fractions,
            "macro_equivalent": dict(fractions),
            "physical_proton": dict(fractions),
        },
    }


def _paired_folded_window_transitions(
    reference: FinalDtData,
    candidate: FinalDtData,
    period_ns: float,
) -> dict[str, Any]:
    """Pair local folded-window classifications by full-input row position.

    The controlled validator guarantees equal full-file populations, rowwise
    lineage labels, and all-survivor finite arrays. Older baseline files do not
    store particle IDs, so this is an explicit deterministic-loader contract,
    not a cryptographically verified particle-ID join.
    """

    if reference.dt_ns.size != candidate.dt_ns.size or not np.array_equal(
        reference.bunch_index, candidate.bunch_index
    ):
        raise ValueError("Controlled transition-gamma particle rows do not align")
    ref_outside = (
        np.abs(_fold_to_nearest_h28_center(reference.dt_ns, period_ns))
        > _LINEAGE_WINDOW_HALF_WIDTH_NS
    )
    cand_outside = (
        np.abs(_fold_to_nearest_h28_center(candidate.dt_ns, period_ns))
        > _LINEAGE_WINDOW_HALF_WIDTH_NS
    )
    per_lineage: dict[str, Any] = {}
    for label in np.unique(reference.bunch_index):
        selected = reference.bunch_index == label
        per_lineage[str(int(label))] = _paired_window_transition_block(
            ref_outside[selected],
            cand_outside[selected],
            reference.macro_equivalent_weight,
            reference.physical_proton_weight,
        )
    return {
        "available": True,
        "direction": "reference_to_candidate",
        "metric": "nearest_h28_bucket_folded_125ns_local_window",
        "pairing": (
            "row-index pairing under same-full-file deterministic-loader contract; "
            "no stored particle_id in legacy baseline"
        ),
        "global": _paired_window_transition_block(
            ref_outside,
            cand_outside,
            reference.macro_equivalent_weight,
            reference.physical_proton_weight,
        ),
        "per_initial_lineage": per_lineage,
    }


def _retention_differences(
    reference: dict[str, float | int],
    candidate: dict[str, float | int],
) -> dict[str, float | int]:
    """Compare aggregate retention in each lineage's reference h=28 bucket."""

    return {
        "candidate_minus_reference_raw_in_reference_bucket": int(
            candidate["raw_in_reference_bucket"]
            - reference["raw_in_reference_bucket"]
        ),
        "candidate_minus_reference_raw_outside_reference_bucket": int(
            candidate["raw_outside_reference_bucket"]
            - reference["raw_outside_reference_bucket"]
        ),
        "candidate_minus_reference_fraction_in_reference_bucket": float(
            candidate["fraction_in_reference_bucket"]
            - reference["fraction_in_reference_bucket"]
        ),
        "candidate_minus_reference_fraction_outside_reference_bucket": float(
            candidate["fraction_outside_reference_bucket"]
            - reference["fraction_outside_reference_bucket"]
        ),
        "absolute_fraction_outside_reference_bucket_difference": float(
            abs(
                candidate["fraction_outside_reference_bucket"]
                - reference["fraction_outside_reference_bucket"]
            )
        ),
        "candidate_minus_reference_outside_reference_bucket_macro_equivalent_count": float(
            candidate["outside_reference_bucket_macro_equivalent_count"]
            - reference["outside_reference_bucket_macro_equivalent_count"]
        ),
    }


def _prepare_pair(
    reference_path: str | Path,
    candidate_path: str | Path,
    *,
    allow_gamma_transition_difference: bool = False,
) -> tuple[FinalDtData, FinalDtData, float]:
    """Load a compatible result pair and return final data plus h=28 period."""

    reference_result = load_result_data(reference_path)
    candidate_result = load_result_data(candidate_path)
    _strict_physics_compatibility(
        reference_result,
        candidate_result,
        allow_gamma_transition_difference=allow_gamma_transition_difference,
    )
    if "rf25_frequency_hz" not in reference_result.machine_attrs:
        raise ValueError(
            f"{reference_result.path}: missing machine attribute rf25_frequency_hz"
        )
    frequency_hz = float(reference_result.machine_attrs["rf25_frequency_hz"])
    if not np.isfinite(frequency_hz) or frequency_hz <= 0:
        raise ValueError("rf25_frequency_hz must be positive and finite")
    return (
        _extract_final_dt(reference_result),
        _extract_final_dt(candidate_result),
        1e9 / frequency_hz,
    )


def _small_amplitude_synchrotron_period_ms(data: FinalDtData) -> float | None:
    """Estimate the final h=28 small-amplitude period when metadata permits."""

    machine = data.result.machine_attrs
    ramp = data.result.rf_program_attrs
    required_machine = (
        "harmonic_25",
        "slip_factor",
        "beta0",
        "p0c_ev",
        "revolution_period_s",
    )
    missing = [key for key in required_machine if key not in machine]
    if missing or "rf25_final_v" not in ramp:
        return None
    tune = np.sqrt(
        float(machine["harmonic_25"])
        * abs(float(machine["slip_factor"]))
        * float(ramp["rf25_final_v"])
        / (
            2.0
            * np.pi
            * float(machine["beta0"]) ** 2
            * (float(machine["p0c_ev"]) / float(machine["beta0"]))
        )
    )
    return float(float(machine["revolution_period_s"]) / tune * 1e3)


def compare_final_dt(
    reference_path: str | Path,
    candidate_path: str | Path,
    *,
    absolute_bins: int = 2400,
    absolute_range_ns: tuple[float, float] = (-1000.0, 3800.0),
    bucket_relative_bins: int = 400,
    bucket_relative_range_ns: tuple[float, float] | None = None,
    allow_gamma_transition_difference: bool = False,
) -> dict[str, Any]:
    """Return a JSON-safe strict comparison of final BLonD-convention dt.

    The map, RF program, endpoint, source scaling, and recorded profile grid
    must match.  Tracked sample sizes and their associated constant weights may
    differ and are recorded explicitly.  Aggregate and per-lineage comparisons
    include absolute time, a common historical bucket-relative frame, tails,
    and final h=28 bucket occupancy. The nearest-bucket folded +/-125 ns metric
    is the Werkema-aligned Recycler-endpoint local-window proxy. A separate
    unwrapped lineage metric combines local tails with whole-bucket migration.
    Neither is a downstream extinction prediction. Set
    ``allow_gamma_transition_difference`` only for the controlled full-input
    study; all differences except gamma transition, alpha_c, and eta still fail.
    """

    reference, candidate, period_ns = _prepare_pair(
        reference_path,
        candidate_path,
        allow_gamma_transition_difference=allow_gamma_transition_difference,
    )
    # Establish centers from the historical sample only, then reuse them for
    # both samples so centroid motion and adjacent-bucket migration remain visible.
    reference_centers = _reference_bucket_centers(reference, period_ns)
    if bucket_relative_range_ns is None:
        bucket_relative_range_ns = (-1.5 * period_ns, 1.5 * period_ns)
    absolute_edges = _validate_histogram_spec(
        absolute_bins, absolute_range_ns, "absolute-dt histogram"
    )
    relative_edges = _validate_histogram_spec(
        bucket_relative_bins,
        bucket_relative_range_ns,
        "bucket-relative histogram",
    )

    reference_summary = _distribution_summary(
        reference, period_ns, reference_centers
    )
    candidate_summary = _distribution_summary(
        candidate, period_ns, reference_centers
    )
    ref_relative = _relative_to_reference_centers(reference, reference_centers)
    cand_relative = _relative_to_reference_centers(candidate, reference_centers)
    common_labels = sorted(
        set(reference_summary["per_bunch"]) & set(candidate_summary["per_bunch"]),
        key=int,
    )
    per_bunch: dict[str, Any] = {}
    for label in common_labels:
        label_value = int(label)
        ref_absolute = reference.dt_ns[reference.bunch_index == label_value]
        cand_absolute = candidate.dt_ns[candidate.bunch_index == label_value]
        ref_entry = reference_summary["per_bunch"][label]
        cand_entry = candidate_summary["per_bunch"][label]
        reference_center = reference_centers[label_value]
        ref_bucket = ref_absolute - reference_center
        cand_bucket = cand_absolute - reference_center
        per_bunch[label] = {
            "absolute_empirical_wasserstein_1_ns": _empirical_wasserstein_1_ns(
                ref_absolute, cand_absolute
            ),
            "absolute_histogram": _histogram_distances(
                ref_absolute,
                cand_absolute,
                reference.macro_equivalent_weight,
                candidate.macro_equivalent_weight,
                absolute_edges,
            ),
            "absolute_summary_differences": _summary_differences(
                ref_entry["absolute_dt"], cand_entry["absolute_dt"]
            ),
            "bucket_relative_histogram": _histogram_distances(
                ref_bucket,
                cand_bucket,
                reference.macro_equivalent_weight,
                candidate.macro_equivalent_weight,
                relative_edges,
            ),
            "bucket_relative_empirical_wasserstein_1_ns": (
                _empirical_wasserstein_1_ns(ref_bucket, cand_bucket)
            ),
            "bucket_relative_summary_differences": _summary_differences(
                ref_entry["bucket_relative_dt"], cand_entry["bucket_relative_dt"]
            ),
            "sample_median_bucket_center_candidate_minus_reference_ns": float(
                cand_entry["sample_median_nearest_h28_bucket_center_ns"]
                - ref_entry["sample_median_nearest_h28_bucket_center_ns"]
            ),
            "bucket_relative_tail_differences": _tail_differences(
                ref_entry["bucket_relative_tails"],
                cand_entry["bucket_relative_tails"],
            ),
            "final_h28_bucket_occupancy_differences": _occupancy_differences(
                ref_entry["final_h28_bucket_occupancy"],
                cand_entry["final_h28_bucket_occupancy"],
            ),
        }

    reference_period_ms = _small_amplitude_synchrotron_period_ms(reference)
    candidate_period_ms = _small_amplitude_synchrotron_period_ms(candidate)
    period_difference_ms = (
        float(candidate_period_ms - reference_period_ms)
        if reference_period_ms is not None and candidate_period_ms is not None
        else None
    )

    return {
        "schema_version": 1,
        "reference_path": str(reference.result.path),
        "candidate_path": str(candidate.result.path),
        "validation": {
            "strict_same_physics_and_endpoint": (
                not allow_gamma_transition_difference
            ),
            "controlled_transition_gamma_comparison": (
                allow_gamma_transition_difference
            ),
            "independent_machine_parameter_changed": (
                "gamma_transition" if allow_gamma_transition_difference else None
            ),
            "derived_machine_parameter_changes": (
                ["momentum_compaction_factor", "slip_factor"]
                if allow_gamma_transition_difference
                else []
            ),
            "allowed_machine_attribute_differences": (
                list(_GAMMA_TRANSITION_MACHINE_ATTRS)
                if allow_gamma_transition_difference
                else []
            ),
            "same_full_file_backed_input_required": (
                allow_gamma_transition_difference
            ),
            "sampling_sizes_and_weights_may_differ": (
                not allow_gamma_transition_difference
            ),
            "coordinate_convention": "dt_ns=-zeta_m/(beta0*c)*1e9",
            "survivor_selection": "state>0 and finite(zeta_m)",
        },
        "machine_parameter_comparison": {
            key: {
                "reference": float(reference.result.machine_attrs[key]),
                "candidate": float(candidate.result.machine_attrs[key]),
                "candidate_minus_reference": float(
                    candidate.result.machine_attrs[key]
                    - reference.result.machine_attrs[key]
                ),
                "candidate_to_reference_ratio": float(
                    candidate.result.machine_attrs[key]
                    / reference.result.machine_attrs[key]
                ),
                "candidate_minus_reference_percent_of_reference": float(
                    100.0
                    * (
                        candidate.result.machine_attrs[key]
                        / reference.result.machine_attrs[key]
                        - 1.0
                    )
                ),
            }
            for key in _GAMMA_TRANSITION_MACHINE_ATTRS
            if key in reference.result.machine_attrs
            and key in candidate.result.machine_attrs
        },
        "derived_final_voltage_synchrotron_period_comparison": {
            "approximation": "small-amplitude h=28 period at rf25_final_v",
            "available": period_difference_ms is not None,
            "reference_ms": reference_period_ms,
            "candidate_ms": candidate_period_ms,
            "candidate_minus_reference_ms": period_difference_ms,
        },
        "configuration": {
            "absolute_histogram_bins": int(absolute_bins),
            "absolute_histogram_range_ns": [
                float(absolute_edges[0]),
                float(absolute_edges[-1]),
            ],
            "bucket_relative_histogram_bins": int(bucket_relative_bins),
            "bucket_relative_histogram_range_ns": [
                float(relative_edges[0]),
                float(relative_edges[-1]),
            ],
            "h28_bucket_period_ns": float(period_ns),
            "histogram_normalization": (
                "Each sample is normalized to unit in-range mass before W1 and JS; "
                "in-range fractions and effective counts are reported separately."
            ),
            "bucket_relative_convention": (
                "For each initial bunch label, both samples subtract the historical "
                "sample's nearest h=28 center without modulo wrapping."
            ),
        },
        "reference": reference_summary,
        "candidate": candidate_summary,
        "comparison": {
            "absolute_empirical_wasserstein_1_ns": _empirical_wasserstein_1_ns(
                reference.dt_ns, candidate.dt_ns
            ),
            "absolute_histogram": _histogram_distances(
                reference.dt_ns,
                candidate.dt_ns,
                reference.macro_equivalent_weight,
                candidate.macro_equivalent_weight,
                absolute_edges,
            ),
            "absolute_summary_differences": _summary_differences(
                reference_summary["absolute_dt"], candidate_summary["absolute_dt"]
            ),
            "bucket_relative_histogram": _histogram_distances(
                ref_relative,
                cand_relative,
                reference.macro_equivalent_weight,
                candidate.macro_equivalent_weight,
                relative_edges,
            ),
            "bucket_relative_empirical_wasserstein_1_ns": (
                _empirical_wasserstein_1_ns(ref_relative, cand_relative)
            ),
            "bucket_relative_summary_differences": _summary_differences(
                reference_summary["bucket_relative_dt"],
                candidate_summary["bucket_relative_dt"],
            ),
            "bucket_relative_tail_differences": _tail_differences(
                reference_summary["bucket_relative_tails"],
                candidate_summary["bucket_relative_tails"],
            ),
            "initial_lineage_reference_bucket_retention_differences": (
                _retention_differences(
                    reference_summary[
                        "initial_lineage_reference_bucket_retention"
                    ],
                    candidate_summary[
                        "initial_lineage_reference_bucket_retention"
                    ],
                )
            ),
            "unwrapped_initial_lineage_125ns_window_differences": (
                _lineage_window_differences(
                    reference_summary["unwrapped_initial_lineage_125ns_window"],
                    candidate_summary["unwrapped_initial_lineage_125ns_window"],
                )
            ),
            "nearest_h28_bucket_folded_125ns_local_window_differences": (
                _lineage_window_differences(
                    reference_summary[
                        "nearest_h28_bucket_folded_125ns_local_window"
                    ],
                    candidate_summary[
                        "nearest_h28_bucket_folded_125ns_local_window"
                    ],
                )
            ),
            "paired_nearest_h28_bucket_folded_125ns_transitions": (
                _paired_folded_window_transitions(reference, candidate, period_ns)
                if allow_gamma_transition_difference
                else {
                    "available": False,
                    "reason": (
                        "paired row-index transitions are only validated for the "
                        "controlled full-input transition-gamma comparison"
                    ),
                }
            ),
            "final_h28_bucket_occupancy_differences": _occupancy_differences(
                reference_summary["final_h28_bucket_occupancy"],
                candidate_summary["final_h28_bucket_occupancy"],
            ),
            "common_bunch_labels": common_labels,
            "reference_only_bunch_labels": sorted(
                set(reference_summary["per_bunch"])
                - set(candidate_summary["per_bunch"]),
                key=int,
            ),
            "candidate_only_bunch_labels": sorted(
                set(candidate_summary["per_bunch"])
                - set(reference_summary["per_bunch"]),
                key=int,
            ),
            "per_bunch": per_bunch,
        },
    }


def write_final_dt_comparison_report(
    reference_path: str | Path,
    candidate_path: str | Path,
    output_path: str | Path,
    **comparison_options: Any,
) -> Path:
    """Compute :func:`compare_final_dt` and write standards-compliant JSON."""

    report = compare_final_dt(reference_path, candidate_path, **comparison_options)
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return path


def write_final_dt_comparison_plot(
    reference_path: str | Path,
    candidate_path: str | Path,
    output_path: str | Path,
    *,
    absolute_bins: int = 2400,
    absolute_range_ns: tuple[float, float] = (-1000.0, 3800.0),
    bucket_relative_bins: int = 400,
    bucket_relative_range_ns: tuple[float, float] | None = None,
    reference_label: str = "Historical input",
    candidate_label: str = "Generated input",
    allow_gamma_transition_difference: bool = False,
) -> Path:
    """Plot endpoint density, residual, tail survival, and local-window rates.

    Absolute densities are independently normalized. The lower panels separate
    the nearest-bucket folded local pulse from the unwrapped lineage coordinate
    that retains bucket migration, and show the exact folded +/-125 ns rate by
    initial lineage. The optional transition-gamma flag uses the same narrow
    validation as the JSON comparison.
    """

    reference, candidate, period_ns = _prepare_pair(
        reference_path,
        candidate_path,
        allow_gamma_transition_difference=allow_gamma_transition_difference,
    )
    reference_centers = _reference_bucket_centers(reference, period_ns)
    if bucket_relative_range_ns is None:
        bucket_relative_range_ns = (-1.5 * period_ns, 1.5 * period_ns)
    absolute_edges = _validate_histogram_spec(
        absolute_bins, absolute_range_ns, "absolute-dt histogram"
    )
    relative_edges = _validate_histogram_spec(
        bucket_relative_bins,
        bucket_relative_range_ns,
        "bucket-relative histogram",
    )

    # Configure the same safe noninteractive Matplotlib behavior as plotting.py.
    mpl_dir = Path(gettempdir()) / "xsuite-recycler-matplotlib"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_dir))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Normalize each sample independently: this figure compares longitudinal
    # shape rather than the number of macroparticles selected for tracking.
    ref_abs = _normalize(
        _histogram(reference.dt_ns, reference.macro_equivalent_weight, absolute_edges)
    )
    cand_abs = _normalize(
        _histogram(candidate.dt_ns, candidate.macro_equivalent_weight, absolute_edges)
    )
    ref_relative_values = _relative_to_reference_centers(
        reference, reference_centers
    )
    cand_relative_values = _relative_to_reference_centers(
        candidate, reference_centers
    )
    ref_local_values = _fold_to_nearest_h28_center(reference.dt_ns, period_ns)
    cand_local_values = _fold_to_nearest_h28_center(candidate.dt_ns, period_ns)
    if ref_abs is None or cand_abs is None:
        raise ValueError("Plot ranges must contain particles from both result files")

    abs_centers = 0.5 * (absolute_edges[:-1] + absolute_edges[1:])
    abs_width = np.diff(absolute_edges)
    ref_density = ref_abs / abs_width
    cand_density = cand_abs / abs_width
    tail_max_ns = max(
        float(np.max(np.abs(ref_relative_values))),
        float(np.max(np.abs(cand_relative_values))),
        abs(float(relative_edges[0])),
        abs(float(relative_edges[-1])),
    )
    tail_thresholds_ns = np.linspace(
        0.0, tail_max_ns, max(int(bucket_relative_bins), 2)
    )

    def survival_fraction(values_ns: np.ndarray) -> np.ndarray:
        """Return P(|dt| > threshold) on the shared plotting grid."""

        ordered = np.sort(np.abs(values_ns))
        return (
            ordered.size
            - np.searchsorted(ordered, tail_thresholds_ns, side="right")
        ) / ordered.size

    ref_local_survival = survival_fraction(ref_local_values)
    cand_local_survival = survival_fraction(cand_local_values)
    ref_unwrapped_survival = survival_fraction(ref_relative_values)
    cand_unwrapped_survival = survival_fraction(cand_relative_values)
    distances = _histogram_distances(
        reference.dt_ns,
        candidate.dt_ns,
        reference.macro_equivalent_weight,
        candidate.macro_equivalent_weight,
        absolute_edges,
    )
    empirical_w1_ns = _empirical_wasserstein_1_ns(
        reference.dt_ns, candidate.dt_ns
    )

    colors = ("#1f5a94", "#d45d2c")
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(13.4, 8.8),
        constrained_layout=True,
    )
    # Panel 1 overlays absolute final-time probability densities.
    axes[0, 0].step(
        abs_centers,
        ref_density,
        where="mid",
        lw=1.35,
        color=colors[0],
        label=f"{reference_label} (n={reference.n_alive_finite_dt:,})",
    )
    axes[0, 0].step(
        abs_centers,
        cand_density,
        where="mid",
        lw=1.2,
        color=colors[1],
        label=f"{candidate_label} (n={candidate.n_alive_finite_dt:,})",
    )
    axes[0, 0].set(
        xlim=(absolute_edges[0], absolute_edges[-1]),
        ylabel="Normalized density [1/ns]",
        title="Final absolute longitudinal time distribution",
    )
    axes[0, 0].legend(loc="upper right", frameon=False, ncol=1)
    axes[0, 0].text(
        0.015,
        0.95,
        f"Empirical W1 = {empirical_w1_ns:.3f} ns\n"
        f"Histogram W1 = {distances['wasserstein_1_ns']:.3f} ns\n"
        f"JS distance = {distances['jensen_shannon_distance']:.4f}\n"
        "Plotted in-range fractions = "
        f"{distances['reference_in_range_fraction']:.6f}, "
        f"{distances['candidate_in_range_fraction']:.6f}",
        transform=axes[0, 0].transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "0.8"},
    )

    # Panel 2 makes small local excesses and deficits visible directly.
    difference = cand_density - ref_density
    axes[0, 1].axhline(0.0, color="0.25", lw=0.8)
    axes[0, 1].step(
        abs_centers, difference, where="mid", color="#6a3d9a", lw=1.0
    )
    axes[0, 1].fill_between(
        abs_centers,
        0.0,
        difference,
        step="mid",
        color="#6a3d9a",
        alpha=0.24,
        linewidth=0,
    )
    axes[0, 1].set(
        xlim=(absolute_edges[0], absolute_edges[-1]),
        xlabel=r"Final $\Delta t$ [ns] (BLonD convention)",
        ylabel=f"{candidate_label} - {reference_label}\n[1/ns]",
        title="Candidate-minus-reference density residual",
    )

    # Panel 3 separates the local folded pulse window (solid) from the combined
    # unwrapped local-tail-or-bucket-migration diagnostic (dashed).
    positive = tail_thresholds_ns > 0
    for values, color, label, linestyle in (
        (ref_local_survival, colors[0], f"{reference_label}: local folded", "-"),
        (cand_local_survival, colors[1], f"{candidate_label}: local folded", "-"),
        (
            ref_unwrapped_survival,
            colors[0],
            f"{reference_label}: unwrapped lineage",
            "--",
        ),
        (
            cand_unwrapped_survival,
            colors[1],
            f"{candidate_label}: unwrapped lineage",
            "--",
        ),
    ):
        visible = positive & (values > 0)
        axes[1, 0].plot(
            tail_thresholds_ns[visible],
            values[visible],
            color=color,
            ls=linestyle,
            lw=1.35,
            label=label,
        )
    axes[1, 0].axvline(
        _LINEAGE_WINDOW_HALF_WIDTH_NS,
        color="black",
        ls=":",
        lw=1.2,
        label="125 ns window",
    )
    axes[1, 0].set(
        xlim=(0.0, tail_max_ns),
        yscale="log",
        xlabel=r"Threshold $|\Delta t|$ [ns]",
        ylabel="Fraction beyond threshold",
        title="Folded local tail versus unwrapped combined diagnostic",
    )
    axes[1, 0].legend(loc="best", frameon=False, fontsize=8)

    # Panel 4 reports the Werkema-aligned local folded +/-125 ns proxy for each
    # original structural bunch. Labels stratify particles but do not recenter.
    ref_local_summary = _folded_local_window_diagnostic(reference, period_ns)
    cand_local_summary = _folded_local_window_diagnostic(candidate, period_ns)
    labels = sorted(
        set(ref_local_summary["per_initial_lineage"])
        & set(cand_local_summary["per_initial_lineage"]),
        key=int,
    )
    positions = np.arange(len(labels), dtype=float)
    width = 0.38
    ref_fractions = [
        ref_local_summary["per_initial_lineage"][label]["raw_outside_fraction"]
        for label in labels
    ]
    cand_fractions = [
        cand_local_summary["per_initial_lineage"][label]["raw_outside_fraction"]
        for label in labels
    ]
    axes[1, 1].bar(
        positions - width / 2,
        ref_fractions,
        width,
        color=colors[0],
        label=reference_label,
    )
    axes[1, 1].bar(
        positions + width / 2,
        cand_fractions,
        width,
        color=colors[1],
        label=candidate_label,
    )
    axes[1, 1].set(
        xticks=positions,
        xticklabels=[str(int(label) + 1) for label in labels],
        xlabel="Initial structural bunch",
        ylabel="Fraction outside local +/-125 ns",
        title="Nearest-h=28-bucket local-window proxy",
    )
    axes[1, 1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axes[1, 1].legend(loc="best", frameon=False)
    for ax in axes.ravel():
        ax.grid(True, color="0.9", lw=0.55)

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path
