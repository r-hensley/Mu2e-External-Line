"""Statistical characterization of ordered Recycler input distributions.

The historical input is ordered as 168 consecutive microbunches with an equal
number of macroparticles in each microbunch.  This module deliberately works
with arrival time and ``dp/p`` arrays rather than Xsuite particle objects so it
can be used both for fitting a generator and for validating its output.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


N_MICROBUNCHES = 168
# The source distribution contains two trains of 84 occupied 53 MHz buckets,
# separated by 42 empty bucket slots.  These indices retain that empty gap when
# the measured bunch centers are fitted to a regular time grid.
OCCUPIED_GRID_INDICES = np.concatenate(
    (np.arange(84, dtype=np.float64), np.arange(126, 210, dtype=np.float64))
)


def _as_ordered_matrix(values, name: str, n_microbunches: int) -> np.ndarray:
    """Validate a flat ordered array and reshape it to one row per microbunch."""

    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional array")
    if array.size == 0 or array.size % n_microbunches:
        raise ValueError(
            f"{name} length must be a nonzero multiple of {n_microbunches}"
        )
    if not np.isfinite(array).all():
        raise ValueError(f"{name} contains non-finite values")
    return array.reshape(n_microbunches, -1)


def _validate_pair(dt_s, dpop, n_microbunches: int) -> tuple[np.ndarray, np.ndarray]:
    """Return compatible time and momentum matrices with usable row samples."""

    time = _as_ordered_matrix(dt_s, "dt_s", n_microbunches)
    momentum = _as_ordered_matrix(dpop, "dpop", n_microbunches)
    if time.shape != momentum.shape:
        raise ValueError("dt_s and dpop must have the same length")
    if time.shape[1] < 2:
        raise ValueError("At least two particles per microbunch are required")
    return time, momentum


def _grid_indices(n_microbunches: int) -> np.ndarray:
    """Return occupied RF-grid indices for the fixed two-train input layout."""

    if n_microbunches != N_MICROBUNCHES:
        raise ValueError(
            "The two-train grid fit is defined only for 168 ordered microbunches"
        )
    return OCCUPIED_GRID_INDICES


@dataclass(frozen=True)
class InputCharacterization:
    """Moments and fixed-grid fit for one ordered input distribution."""

    n_particles: int
    n_microbunches: int
    particles_per_microbunch: int
    time_center_ns: np.ndarray
    time_sigma_ns: np.ndarray
    dpop_mean: np.ndarray
    dpop_sigma: np.ndarray
    correlation: np.ndarray
    grid_origin_ns: float
    grid_spacing_ns: float
    grid_center_residual_rms_ns: float
    pooled_time_sigma_ns: float
    pooled_dpop_sigma: float
    pooled_correlation: float
    global_dpop_mean: float

    def to_summary(self, *, include_per_microbunch: bool = False) -> dict:
        """Return JSON-safe aggregates and, optionally, all microbunch moments."""

        result = {
            "n_particles": self.n_particles,
            "n_microbunches": self.n_microbunches,
            "particles_per_microbunch": self.particles_per_microbunch,
            "grid": {
                "origin_ns": self.grid_origin_ns,
                "spacing_ns": self.grid_spacing_ns,
                "center_residual_rms_ns": self.grid_center_residual_rms_ns,
            },
            "pooled_within_microbunch": {
                "time_sigma_ns": self.pooled_time_sigma_ns,
                "dpop_sigma": self.pooled_dpop_sigma,
                "correlation": self.pooled_correlation,
            },
            "global_dpop_mean": self.global_dpop_mean,
            "per_microbunch_aggregates": {
                "time_sigma_mean_ns": float(np.mean(self.time_sigma_ns)),
                "time_sigma_std_ns": float(np.std(self.time_sigma_ns, ddof=1)),
                "dpop_sigma_mean": float(np.mean(self.dpop_sigma)),
                "dpop_sigma_std": float(np.std(self.dpop_sigma, ddof=1)),
                "dpop_mean_std": float(np.std(self.dpop_mean, ddof=1)),
                "correlation_rms": float(np.sqrt(np.mean(self.correlation**2))),
            },
        }
        if include_per_microbunch:
            result["per_microbunch"] = {
                "time_center_ns": self.time_center_ns.tolist(),
                "time_sigma_ns": self.time_sigma_ns.tolist(),
                "dpop_mean": self.dpop_mean.tolist(),
                "dpop_sigma": self.dpop_sigma.tolist(),
                "correlation": self.correlation.tolist(),
            }
        return result


def characterize_input(
    dt_s,
    dpop,
    *,
    n_microbunches: int = N_MICROBUNCHES,
) -> InputCharacterization:
    """Measure geometry and within-microbunch moments of an ordered distribution.

    Consecutive equal-size blocks are interpreted as individual microbunches.
    The returned widths use sample standard deviations, while pooled values
    combine only within-microbunch residuals and exclude train geometry.
    """

    time, momentum = _validate_pair(dt_s, dpop, n_microbunches)
    per_microbunch = time.shape[1]

    # Remove each microbunch centroid before measuring its intrinsic time-energy
    # shape. A global width would be dominated by microbunch separation.
    time_center_s = np.mean(time, axis=1)
    dpop_mean = np.mean(momentum, axis=1)
    time_residual_s = time - time_center_s[:, None]
    dpop_residual = momentum - dpop_mean[:, None]

    time_ss = np.sum(time_residual_s**2, axis=1)
    dpop_ss = np.sum(dpop_residual**2, axis=1)
    cross = np.sum(time_residual_s * dpop_residual, axis=1)
    time_sigma_ns = np.sqrt(time_ss / (per_microbunch - 1)) * 1e9
    dpop_sigma = np.sqrt(dpop_ss / (per_microbunch - 1))
    correlation = np.divide(
        cross,
        np.sqrt(time_ss * dpop_ss),
        out=np.zeros_like(cross),
        where=(time_ss > 0.0) & (dpop_ss > 0.0),
    )

    time_center_ns = time_center_s * 1e9
    grid = _grid_indices(n_microbunches)
    # Fit center = origin + spacing * occupied_grid_index.  The skipped indices
    # in ``grid`` ensure the empty interval between the two trains is retained.
    design = np.column_stack((np.ones(n_microbunches), grid))
    grid_origin_ns, grid_spacing_ns = np.linalg.lstsq(
        design, time_center_ns, rcond=None
    )[0]
    grid_residual_ns = time_center_ns - (
        grid_origin_ns + grid_spacing_ns * grid
    )

    # Each fitted microbunch mean consumes one degree of freedom. Pooling the
    # residual sums of squares gives one representative within-microbunch width.
    pooled_dof = time.size - n_microbunches
    pooled_time_sigma_ns = float(
        np.sqrt(np.sum(time_residual_s**2) / pooled_dof) * 1e9
    )
    pooled_dpop_sigma = float(np.sqrt(np.sum(dpop_residual**2) / pooled_dof))
    pooled_correlation = float(
        np.sum(time_residual_s * dpop_residual)
        / np.sqrt(np.sum(time_residual_s**2) * np.sum(dpop_residual**2))
    )

    return InputCharacterization(
        n_particles=int(time.size),
        n_microbunches=n_microbunches,
        particles_per_microbunch=per_microbunch,
        time_center_ns=time_center_ns,
        time_sigma_ns=time_sigma_ns,
        dpop_mean=dpop_mean,
        dpop_sigma=dpop_sigma,
        correlation=correlation,
        grid_origin_ns=float(grid_origin_ns),
        grid_spacing_ns=float(grid_spacing_ns),
        grid_center_residual_rms_ns=float(
            np.sqrt(np.mean(grid_residual_ns**2))
        ),
        pooled_time_sigma_ns=pooled_time_sigma_ns,
        pooled_dpop_sigma=pooled_dpop_sigma,
        pooled_correlation=pooled_correlation,
        global_dpop_mean=float(np.mean(momentum)),
    )


def _wasserstein_1d(first, second) -> float:
    """Return exact empirical 1D Wasserstein distance for any sample sizes.

    Equal-size samples use the sorted-pair shortcut.  Unequal-size samples use
    the equivalent integral of the absolute difference between empirical CDFs.
    """

    first = np.sort(np.asarray(first, dtype=np.float64).ravel())
    second = np.sort(np.asarray(second, dtype=np.float64).ravel())
    if first.size == 0 or second.size == 0:
        raise ValueError("Wasserstein inputs must be nonempty")
    if first.size == second.size:
        return float(np.mean(np.abs(first - second)))

    support = np.sort(np.concatenate((first, second)))
    widths = np.diff(support)
    if widths.size == 0:
        return 0.0
    first_cdf = np.searchsorted(first, support[:-1], side="right") / first.size
    second_cdf = np.searchsorted(second, support[:-1], side="right") / second.size
    return float(np.sum(np.abs(first_cdf - second_cdf) * widths))


def _deterministic_subsample(points: np.ndarray, maximum: int, seed: int) -> np.ndarray:
    """Cap a point cloud reproducibly without changing already-small inputs."""

    if points.shape[0] <= maximum:
        return points
    rng = np.random.default_rng(seed)
    indices = np.sort(rng.choice(points.shape[0], size=maximum, replace=False))
    return points[indices]


def _sliced_wasserstein(
    reference: np.ndarray,
    candidate: np.ndarray,
    *,
    n_angles: int,
    max_samples: int,
) -> tuple[float, float]:
    """Return mean and maximum sliced-Wasserstein distances over 2D projections.

    The projection axes evenly cover the unique half-circle of directions.
    Inputs may be deterministically capped to bound comparison cost.
    """

    if n_angles < 1:
        raise ValueError("n_angles must be positive")
    if max_samples < 2:
        raise ValueError("max_samples must be at least two")
    reference = _deterministic_subsample(reference, max_samples, seed=173)
    candidate = _deterministic_subsample(candidate, max_samples, seed=947)
    angles = (np.arange(n_angles, dtype=np.float64) + 0.5) * np.pi / n_angles
    distances = np.empty(n_angles, dtype=np.float64)
    for ii, angle in enumerate(angles):
        direction = np.asarray((np.cos(angle), np.sin(angle)))
        distances[ii] = _wasserstein_1d(
            reference @ direction, candidate @ direction
        )
    return float(np.mean(distances)), float(np.max(distances))


def _rms(values) -> float:
    """Return root mean square about zero for an array of discrepancies."""

    values = np.asarray(values, dtype=np.float64)
    return float(np.sqrt(np.mean(values**2)))


def _comparison_without_noise(
    reference_dt_s,
    reference_dpop,
    candidate_dt_s,
    candidate_dpop,
    *,
    n_microbunches: int,
    sliced_angles: int,
    sliced_max_samples: int,
) -> dict:
    """Build comparison metrics without estimating reference sampling noise.

    Train placement, per-microbunch moments, and residual-distribution shape are
    reported separately.  Both residual point clouds are scaled by reference
    pooled widths, so the distance units describe discrepancies relative to the
    historical input rather than allowing the candidate to set its own scale.
    """

    ref_time, ref_momentum = _validate_pair(
        reference_dt_s, reference_dpop, n_microbunches
    )
    cand_time, cand_momentum = _validate_pair(
        candidate_dt_s, candidate_dpop, n_microbunches
    )
    reference = characterize_input(
        ref_time.ravel(), ref_momentum.ravel(), n_microbunches=n_microbunches
    )
    candidate = characterize_input(
        cand_time.ravel(), cand_momentum.ravel(), n_microbunches=n_microbunches
    )

    ref_time_residual_ns = (
        ref_time * 1e9 - reference.time_center_ns[:, None]
    )
    cand_time_residual_ns = (
        cand_time * 1e9 - candidate.time_center_ns[:, None]
    )
    ref_dpop_residual = ref_momentum - reference.dpop_mean[:, None]
    cand_dpop_residual = cand_momentum - candidate.dpop_mean[:, None]

    # Center within each microbunch and use the reference widths for *both*
    # samples.  Separate geometry metrics below retain any centroid mismatch.
    ref_standardized = np.column_stack(
        (
            (ref_time_residual_ns / reference.pooled_time_sigma_ns).ravel(),
            (ref_dpop_residual / reference.pooled_dpop_sigma).ravel(),
        )
    )
    cand_standardized = np.column_stack(
        (
            (cand_time_residual_ns / reference.pooled_time_sigma_ns).ravel(),
            (cand_dpop_residual / reference.pooled_dpop_sigma).ravel(),
        )
    )
    sliced_mean, sliced_max = _sliced_wasserstein(
        ref_standardized,
        cand_standardized,
        n_angles=sliced_angles,
        max_samples=sliced_max_samples,
    )

    center_difference_ns = candidate.time_center_ns - reference.time_center_ns
    time_sigma_relative = (
        candidate.time_sigma_ns - reference.time_sigma_ns
    ) / reference.time_sigma_ns
    dpop_sigma_relative = (
        candidate.dpop_sigma - reference.dpop_sigma
    ) / reference.dpop_sigma
    dpop_mean_difference = candidate.dpop_mean - reference.dpop_mean
    correlation_difference = candidate.correlation - reference.correlation

    ref_fraction = np.full(n_microbunches, 1.0 / n_microbunches)
    cand_fraction = np.full(n_microbunches, 1.0 / n_microbunches)
    return {
        "reference": reference.to_summary(),
        "candidate": candidate.to_summary(),
        "center_geometry": {
            "origin_difference_ns": (
                candidate.grid_origin_ns - reference.grid_origin_ns
            ),
            "spacing_difference_ns": (
                candidate.grid_spacing_ns - reference.grid_spacing_ns
            ),
            "center_bias_ns": float(np.mean(center_difference_ns)),
            "center_rmse_ns": _rms(center_difference_ns),
            "center_max_abs_ns": float(np.max(np.abs(center_difference_ns))),
            "reference_grid_residual_rms_ns": (
                reference.grid_center_residual_rms_ns
            ),
            "candidate_grid_residual_rms_ns": (
                candidate.grid_center_residual_rms_ns
            ),
        },
        "pooled_shape": {
            "reference_time_sigma_ns": reference.pooled_time_sigma_ns,
            "candidate_time_sigma_ns": candidate.pooled_time_sigma_ns,
            "time_sigma_relative_difference": (
                candidate.pooled_time_sigma_ns / reference.pooled_time_sigma_ns
                - 1.0
            ),
            "reference_dpop_sigma": reference.pooled_dpop_sigma,
            "candidate_dpop_sigma": candidate.pooled_dpop_sigma,
            "dpop_sigma_relative_difference": (
                candidate.pooled_dpop_sigma / reference.pooled_dpop_sigma - 1.0
            ),
            "reference_correlation": reference.pooled_correlation,
            "candidate_correlation": candidate.pooled_correlation,
            "global_dpop_mean_difference": (
                candidate.global_dpop_mean - reference.global_dpop_mean
            ),
        },
        "distribution_distances": {
            "time_residual_wasserstein_standardized": _wasserstein_1d(
                ref_standardized[:, 0], cand_standardized[:, 0]
            ),
            "dpop_residual_wasserstein_standardized": _wasserstein_1d(
                ref_standardized[:, 1], cand_standardized[:, 1]
            ),
            "sliced_wasserstein_2d_standardized_mean": sliced_mean,
            "sliced_wasserstein_2d_standardized_max": sliced_max,
            "sliced_wasserstein_angles": sliced_angles,
            "sliced_wasserstein_max_samples_per_input": sliced_max_samples,
        },
        "per_microbunch_discrepancies": {
            "time_sigma_relative_rmse": _rms(time_sigma_relative),
            "time_sigma_relative_max_abs": float(
                np.max(np.abs(time_sigma_relative))
            ),
            "dpop_mean_rmse": _rms(dpop_mean_difference),
            "dpop_mean_rmse_in_reference_sigma": (
                _rms(dpop_mean_difference) / reference.pooled_dpop_sigma
            ),
            "dpop_sigma_relative_rmse": _rms(dpop_sigma_relative),
            "dpop_sigma_relative_max_abs": float(
                np.max(np.abs(dpop_sigma_relative))
            ),
            "correlation_rmse": _rms(correlation_difference),
        },
        "population": {
            "reference_particles_per_microbunch": (
                reference.particles_per_microbunch
            ),
            "candidate_particles_per_microbunch": (
                candidate.particles_per_microbunch
            ),
            "candidate_to_reference_particles_per_microbunch": (
                candidate.particles_per_microbunch
                / reference.particles_per_microbunch
            ),
            "normalized_fraction_l1": float(
                np.sum(np.abs(cand_fraction - ref_fraction))
            ),
            "normalized_fraction_max_abs": float(
                np.max(np.abs(cand_fraction - ref_fraction))
            ),
        },
    }


def split_half_noise_metrics(
    dt_s,
    dpop,
    *,
    n_microbunches: int = N_MICROBUNCHES,
    seed: int = 202208,
    sliced_angles: int = 16,
    sliced_max_samples: int = 50_000,
) -> dict:
    """Estimate a conservative half-sample control from two balanced halves.

    Every microbunch contributes the same number of particles to each half, so
    the result measures finite-sample variation without changing train weights.
    """

    time, momentum = _validate_pair(dt_s, dpop, n_microbunches)
    half_size = time.shape[1] // 2
    if half_size < 2:
        raise ValueError("At least four particles per microbunch are required")

    rng = np.random.default_rng(seed)
    first_time = np.empty((n_microbunches, half_size), dtype=np.float64)
    second_time = np.empty_like(first_time)
    first_dpop = np.empty_like(first_time)
    second_dpop = np.empty_like(first_time)
    for microbunch in range(n_microbunches):
        # Split a fresh within-microbunch permutation into disjoint balanced halves.
        # If the row length is odd, its one unpaired particle is intentionally
        # left out rather than making the control populations unequal.
        order = rng.permutation(time.shape[1])
        first = order[:half_size]
        second = order[half_size : 2 * half_size]
        first_time[microbunch] = time[microbunch, first]
        first_dpop[microbunch] = momentum[microbunch, first]
        second_time[microbunch] = time[microbunch, second]
        second_dpop[microbunch] = momentum[microbunch, second]

    comparison = _comparison_without_noise(
        first_time.ravel(),
        first_dpop.ravel(),
        second_time.ravel(),
        second_dpop.ravel(),
        n_microbunches=n_microbunches,
        sliced_angles=sliced_angles,
        sliced_max_samples=sliced_max_samples,
    )
    return {
        "particles_per_microbunch_per_half": half_size,
        "center_geometry": comparison["center_geometry"],
        "pooled_shape": comparison["pooled_shape"],
        "distribution_distances": comparison["distribution_distances"],
        "per_microbunch_discrepancies": comparison[
            "per_microbunch_discrepancies"
        ],
    }


def compare_inputs(
    reference_dt_s,
    reference_dpop,
    candidate_dt_s,
    candidate_dpop,
    *,
    n_microbunches: int = N_MICROBUNCHES,
    sliced_angles: int = 32,
    sliced_max_samples: int = 100_000,
    include_reference_split_half: bool = True,
    split_seed: int = 202208,
) -> dict:
    """Compare a candidate input against a historical reference distribution.

    Shape distances use within-microbunch-centered coordinates scaled by the
    reference pooled widths.  Grid and per-microbunch moment differences are
    reported separately, so a good shape score cannot hide incorrect train
    placement or widths.
    """
    result = _comparison_without_noise(
        reference_dt_s,
        reference_dpop,
        candidate_dt_s,
        candidate_dpop,
        n_microbunches=n_microbunches,
        sliced_angles=sliced_angles,
        sliced_max_samples=sliced_max_samples,
    )
    if include_reference_split_half:
        result["reference_split_half_noise"] = split_half_noise_metrics(
            reference_dt_s,
            reference_dpop,
            n_microbunches=n_microbunches,
            seed=split_seed,
            sliced_angles=sliced_angles,
            sliced_max_samples=sliced_max_samples,
        )
        result["reference_split_half_note"] = (
            "This is a conservative half-sample control, not a matched-size "
            "noise floor for the full reference and candidate inputs."
        )
    return result
