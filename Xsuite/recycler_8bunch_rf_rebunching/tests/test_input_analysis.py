import json

import numpy as np

from xsuite_recycler.input_analysis import characterize_input, compare_inputs


def _synthetic_input(
    seed,
    *,
    per_microbunch=240,
    origin_ns=-593.8,
    spacing_ns=19.1547,
    sigma_time_ns=5.025,
    sigma_dpop=1.7594e-4,
    correlation=0.0,
):
    rng = np.random.default_rng(seed)
    occupied = np.concatenate((np.arange(84), np.arange(126, 210)))
    centers_ns = origin_ns + spacing_ns * occupied
    first = rng.standard_normal((168, per_microbunch))
    second = rng.standard_normal((168, per_microbunch))
    dt_s = (
        centers_ns[:, None] + sigma_time_ns * first
    ).ravel() * 1e-9
    dpop = (
        sigma_dpop
        * (correlation * first + np.sqrt(1.0 - correlation**2) * second)
    ).ravel()
    return dt_s, dpop


def test_characterization_recovers_grid_and_common_gaussian_widths():
    dt_s, dpop = _synthetic_input(7, per_microbunch=500, correlation=0.12)
    summary = characterize_input(dt_s, dpop)

    assert summary.n_particles == 168 * 500
    assert abs(summary.grid_origin_ns - (-593.8)) < 0.12
    assert abs(summary.grid_spacing_ns - 19.1547) < 0.001
    assert abs(summary.pooled_time_sigma_ns / 5.025 - 1.0) < 0.01
    assert abs(summary.pooled_dpop_sigma / 1.7594e-4 - 1.0) < 0.01
    assert abs(summary.pooled_correlation - 0.12) < 0.015
    assert summary.grid_center_residual_rms_ns < 0.3
    json.dumps(summary.to_summary(include_per_microbunch=True))


def test_comparison_distinguishes_width_spacing_and_correlation_mismatch():
    reference = _synthetic_input(10)
    matched = _synthetic_input(11)
    mismatched = _synthetic_input(
        12,
        spacing_ns=19.30,
        sigma_time_ns=6.8,
        sigma_dpop=2.3e-4,
        correlation=0.35,
    )

    close = compare_inputs(
        *reference,
        *matched,
        sliced_angles=12,
        sliced_max_samples=12_000,
    )
    far = compare_inputs(
        *reference,
        *mismatched,
        sliced_angles=12,
        sliced_max_samples=12_000,
    )

    assert abs(close["center_geometry"]["spacing_difference_ns"]) < 0.003
    assert far["center_geometry"]["spacing_difference_ns"] > 0.14
    assert abs(close["pooled_shape"]["time_sigma_relative_difference"]) < 0.03
    assert far["pooled_shape"]["time_sigma_relative_difference"] > 0.30
    assert far["pooled_shape"]["dpop_sigma_relative_difference"] > 0.25
    assert far["pooled_shape"]["candidate_correlation"] > 0.30
    assert (
        far["distribution_distances"]["sliced_wasserstein_2d_standardized_mean"]
        > 3
        * close["distribution_distances"][
            "sliced_wasserstein_2d_standardized_mean"
        ]
    )
    assert close["reference_split_half_noise"][
        "particles_per_microbunch_per_half"
    ] == 120
    json.dumps(far)


def test_rejects_nonordered_or_too_small_inputs():
    with np.testing.assert_raises_regex(ValueError, "multiple of 168"):
        characterize_input(np.zeros(169), np.zeros(169))
    with np.testing.assert_raises_regex(ValueError, "At least two"):
        characterize_input(np.zeros(168), np.zeros(168))
