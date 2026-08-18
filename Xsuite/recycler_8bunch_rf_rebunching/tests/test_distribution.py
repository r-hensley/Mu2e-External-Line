import numpy as np

from xsuite_recycler.config import RecyclerConfig
from xsuite_recycler.distribution import (
    dt_to_zeta,
    initial_distribution_from_generated,
    load_initial_distribution,
    zeta_to_dt,
)
from xsuite_recycler.generator import (
    MicrobunchGeneratorConfig,
    generate_longitudinal_distribution,
)


def test_time_coordinate_round_trip():
    machine = RecyclerConfig()
    dt = np.asarray([-1e-6, 0.0, 3.8e-6])
    assert np.allclose(zeta_to_dt(dt_to_zeta(dt, machine.beta0), machine.beta0), dt)


def test_stratified_load_and_blond_energy_mapping(tmp_path):
    machine = RecyclerConfig()
    raw = np.column_stack((np.arange(80, dtype=float), np.linspace(-1e-3, 1e-3, 80)))
    path = tmp_path / "initial.txt"
    np.savetxt(path, raw)

    loaded = load_initial_distribution(path, machine, max_particles=16, seed=4)

    assert loaded.n_particles == 16
    assert np.array_equal(np.bincount(loaded.bunch_index), np.full(8, 2))
    assert loaded.macro_equivalent_weight == 5.0
    assert len(loaded.source_sha256) == 64
    assert len(loaded.source_coordinate_sha256) == 64
    assert len(loaded.tracking_coordinate_sha256) == 64
    assert len(loaded.selected_indices_sha256) == 64
    assert loaded.sampling_method == (
        "contiguous-stratified-random-without-replacement"
    )
    assert loaded.sampling_seed == 4
    assert loaded.sampling_strata == 8
    assert loaded.source_path == str(path.resolve())
    assert np.all(np.isfinite(loaded.zeta_m))
    assert np.all(np.isfinite(loaded.ptau))


def test_168_way_sampling_preserves_every_microbunch(tmp_path):
    machine = RecyclerConfig()
    raw = np.column_stack(
        (
            np.arange(168 * 4, dtype=float),
            np.zeros(168 * 4),
        )
    )
    path = tmp_path / "initial.txt"
    np.savetxt(path, raw)

    loaded = load_initial_distribution(
        path,
        machine,
        max_particles=168,
        seed=11,
        sampling_strata=168,
    )
    selected_z = (
        zeta_to_dt(loaded.zeta_m, machine.beta0)
        * machine.input_time_denominator_m_per_s
    )

    selected_rows = np.rint(selected_z).astype(int)
    assert np.array_equal(np.unique(selected_rows // 4), np.arange(168))
    assert np.array_equal(np.bincount(loaded.bunch_index), np.full(8, 21))


def test_generated_distribution_maps_and_scales_like_historical_reference():
    machine = RecyclerConfig()
    generated = generate_longitudinal_distribution(
        MicrobunchGeneratorConfig(particles_per_microbunch=2, seed=7)
    )
    loaded = initial_distribution_from_generated(generated, machine)

    assert loaded.n_particles == 336
    assert loaded.source_rows == 1_075_200
    assert loaded.macro_equivalent_weight == 3_200.0
    assert loaded.source_kind == "generated"
    assert loaded.source_sha256 == generated.config_fingerprint
    assert loaded.source_coordinate_sha256 == generated.realized_coordinate_sha256
    assert len(loaded.tracking_coordinate_sha256) == 64
    assert loaded.selected_indices_sha256 == ""
    assert loaded.sampling_method == "generated-model"
    assert loaded.sampling_seed == 7
    assert loaded.sampling_strata == 0
    assert np.allclose(
        zeta_to_dt(loaded.zeta_m, machine.beta0),
        generated.dt_s,
    )
    assert np.allclose(loaded.ptau, generated.dpop / machine.blond_energy_beta)
