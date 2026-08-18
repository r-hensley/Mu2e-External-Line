import numpy as np
import pytest

from xsuite_recycler.generator import (
    MICROBUNCHES_PER_BUNCH,
    N_BUNCHES,
    N_MICROBUNCHES,
    OCCUPIED_GRID_INDICES,
    MicrobunchGeneratorConfig,
    generate_longitudinal_distribution,
)


def test_generated_distribution_is_reproducible_and_fingerprinted():
    config = MicrobunchGeneratorConfig(particles_per_microbunch=7, seed=12345)

    first = generate_longitudinal_distribution(config)
    second = generate_longitudinal_distribution(config)

    assert np.array_equal(first.dt_s, second.dt_s)
    assert np.array_equal(first.dpop, second.dpop)
    assert first.config_fingerprint == second.config_fingerprint
    assert len(first.config_fingerprint) == 64
    assert first.metadata["config_fingerprint"] == first.config_fingerprint
    assert first.metadata["n_particles"] == N_MICROBUNCHES * 7
    assert first.realized_coordinate_sha256 == second.realized_coordinate_sha256
    assert len(first.realized_coordinate_sha256) == 64
    assert (
        first.metadata["realized_coordinate_sha256"]
        == first.realized_coordinate_sha256
    )
    assert first.metadata["numpy_version"] == np.__version__
    assert first.metadata["rng_bit_generator"] == "PCG64"

    changed_seed = MicrobunchGeneratorConfig(particles_per_microbunch=7, seed=12346)
    changed = generate_longitudinal_distribution(changed_seed)
    assert changed.config_fingerprint != first.config_fingerprint
    assert changed.realized_coordinate_sha256 != first.realized_coordinate_sha256
    assert not np.array_equal(changed.dt_s, first.dt_s)


def test_exact_microbunch_grid_populations_and_bunch_labels():
    particles_per_microbunch = 5
    config = MicrobunchGeneratorConfig(
        particles_per_microbunch=particles_per_microbunch, seed=8
    )
    generated = generate_longitudinal_distribution(config)

    expected_grid = np.concatenate((np.arange(84), np.arange(126, 210)))
    expected_centers_ns = config.time_origin_ns + expected_grid * config.spacing_ns

    assert generated.n_particles == N_MICROBUNCHES * particles_per_microbunch
    assert np.array_equal(generated.occupied_grid_indices, expected_grid)
    assert np.array_equal(generated.occupied_grid_indices, OCCUPIED_GRID_INDICES)
    assert np.array_equal(
        np.bincount(generated.microbunch_index),
        np.full(N_MICROBUNCHES, particles_per_microbunch),
    )
    assert np.array_equal(
        np.bincount(generated.bunch_index),
        np.full(N_BUNCHES, MICROBUNCHES_PER_BUNCH * particles_per_microbunch),
    )
    assert np.allclose(generated.microbunch_center_times_s * 1e9, expected_centers_ns)
    assert np.array_equal(
        generated.center_time_s,
        np.repeat(generated.microbunch_center_times_s, particles_per_microbunch),
    )


def test_generated_moments_follow_configured_bivariate_gaussian():
    config = MicrobunchGeneratorConfig(
        particles_per_microbunch=1_000,
        sigma_time_ns=4.7,
        mean_dpop=2.5e-6,
        sigma_dpop=2.1e-4,
        correlation=-0.31,
        seed=991,
    )
    generated = generate_longitudinal_distribution(config)

    standardized_time = (
        (generated.dt_s - generated.center_time_s) * 1e9 / config.sigma_time_ns
    )
    standardized_dpop = (generated.dpop - config.mean_dpop) / config.sigma_dpop

    assert abs(np.mean(standardized_time)) < 0.01
    assert abs(np.std(standardized_time) - 1.0) < 0.01
    assert abs(np.mean(standardized_dpop)) < 0.01
    assert abs(np.std(standardized_dpop) - 1.0) < 0.01
    assert abs(np.corrcoef(standardized_time, standardized_dpop)[0, 1] - (-0.31)) < 0.01


@pytest.mark.parametrize(
    "kwargs,exception",
    [
        ({"particles_per_microbunch": 0}, ValueError),
        ({"particles_per_microbunch": 2.5}, TypeError),
        ({"seed": -1}, ValueError),
        ({"seed": True}, TypeError),
        ({"time_origin_ns": np.nan}, ValueError),
        ({"spacing_ns": 0.0}, ValueError),
        ({"sigma_time_ns": -1.0}, ValueError),
        ({"sigma_dpop": 0.0}, ValueError),
        ({"correlation": 1.01}, ValueError),
    ],
)
def test_generator_config_rejects_invalid_values(kwargs, exception):
    with pytest.raises(exception):
        MicrobunchGeneratorConfig(**kwargs)


def test_generator_rejects_wrong_config_type():
    with pytest.raises(TypeError):
        generate_longitudinal_distribution({"seed": 3})
