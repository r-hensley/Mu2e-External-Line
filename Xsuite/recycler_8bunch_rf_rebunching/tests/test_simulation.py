import json

import h5py
import numpy as np

from xsuite_recycler.config import RecyclerConfig
from xsuite_recycler.distribution import initial_distribution_from_generated
from xsuite_recycler.generator import (
    MicrobunchGeneratorConfig,
    generate_longitudinal_distribution,
)
from xsuite_recycler.simulation import run_rf_only


def test_short_rf_only_run_writes_conserved_profiles(tmp_path):
    raw = np.column_stack(
        (
            np.linspace(-180.0, 1020.0, 80),
            np.linspace(-2e-4, 2e-4, 80),
        )
    )
    input_path = tmp_path / "initial.txt"
    output_path = tmp_path / "result.h5"
    np.savetxt(input_path, raw)

    run_rf_only(
        input_path,
        output_path,
        max_particles=16,
        record_every=2,
        n_turns=4,
        show_progress=False,
        invocation_argv=["main.py", "run", "--turns", "4", "--no-figures"],
    )

    with h5py.File(output_path, "r") as h5:
        counts = h5["profiles/counts_macro_equivalent"][:]
        outside = h5["profiles/outside_macro_equivalent"][:]
        assert np.array_equal(h5["profiles/turn"][:], [0, 2, 4])
        assert np.allclose(counts.sum(axis=1) + outside, 80)
        assert np.isfinite(h5["final_particles/zeta_m"][:]).all()
        assert np.isfinite(h5["final_particles/ptau"][:]).all()
        assert h5.attrs["sample_particles"] == 16
        assert h5.attrs["model"].startswith("RF-only")
        assert h5.attrs["input_sampling_method"] == (
            "contiguous-stratified-random-without-replacement"
        )
        assert h5.attrs["input_sampling_seed"] == 202208
        assert h5.attrs["input_sampling_strata"] == 8
        assert len(h5.attrs["input_selected_indices_sha256"]) == 64
        assert len(h5.attrs["input_source_coordinate_sha256"]) == 64
        assert len(h5.attrs["input_tracking_coordinate_sha256"]) == 64
        assert len(h5.attrs["project_source_sha256"]) == 64
        assert h5.attrs["project_git_dirty"] in {"true", "false", "unknown"}
        dependencies = json.loads(h5.attrs["dependency_versions_json"])
        assert dependencies["numpy"] == np.__version__
        assert dependencies["xtrack"]
        assert json.loads(h5.attrs["invocation_argv_json"]) == [
            "main.py",
            "run",
            "--turns",
            "4",
            "--no-figures",
        ]
        assert h5.attrs["invocation_command"] == (
            "main.py run --turns 4 --no-figures"
        )
        run_configuration = json.loads(h5.attrs["run_configuration_json"])
        assert run_configuration["tracking"]["resolved_n_turns"] == 4
        assert run_configuration["tracking"]["record_every"] == 2
        assert run_configuration["input_request"]["max_particles"] == 16


def test_short_generated_run_records_generator_provenance(tmp_path):
    machine = RecyclerConfig()
    generated = generate_longitudinal_distribution(
        MicrobunchGeneratorConfig(particles_per_microbunch=1, seed=99)
    )
    distribution = initial_distribution_from_generated(generated, machine)
    output_path = tmp_path / "generated.h5"

    run_rf_only(
        None,
        output_path,
        machine=machine,
        max_particles=distribution.n_particles,
        record_every=1,
        n_turns=1,
        initial_distribution=distribution,
        show_progress=False,
    )

    with h5py.File(output_path, "r") as h5:
        assert h5.attrs["input_source_kind"] == "generated"
        assert h5.attrs["input_sha256"] == generated.config_fingerprint
        assert h5.attrs["input_path"].startswith("generated://")
        assert h5.attrs["input_source_coordinate_sha256"] == (
            generated.realized_coordinate_sha256
        )
        assert h5.attrs["input_sampling_method"] == "generated-model"
        assert h5.attrs["input_sampling_seed"] == 99
        assert h5.attrs["input_sampling_strata"] == 0
        metadata = json.loads(h5.attrs["input_source_metadata_json"])
        assert metadata["numpy_version"] == np.__version__
        assert metadata["rng_bit_generator"] == "PCG64"
        assert metadata["realized_coordinate_sha256"] == (
            generated.realized_coordinate_sha256
        )
        assert h5.attrs["source_rows"] == 1_075_200
        assert np.allclose(
            h5["profiles/counts_macro_equivalent"][:].sum(axis=1)
            + h5["profiles/outside_macro_equivalent"][:],
            1_075_200,
        )
