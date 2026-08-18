# Curated reference results

This directory keeps small, portable checkpoints that can be versioned without
tracking the much larger run HDF5 files under `outputs/`. The copied files are
unchanged from their source output paths. The archived 3 kV files are described
by [`manifest.json`](manifest.json), while the legacy Keegan-code 5 kV
convergence bundle has its own [README](keegan_5kv/README.md),
[manifest](keegan_5kv/manifest.json), and checksum list.

The RF-dependent artifacts under `publication_3kv/` use the **3 kV h=28
publication/ESME default**. The directory name is retained as provenance from
when 3 kV was an explicit variant. Reproduce the artifacts numerically with the
explicit `--rf25-start-kv 3` arguments in the manifest. Exact JSON/HDF5 bytes
can change when provenance fields, absolute paths, or dependency versions
change; the archived files and their hashes remain the fixed reference.

`input_model/` is independent of the RF program. It compares all 1,075,200
historical rows with a seed-202208 generated realization of the same size.

## Contents

- [`input_model/`](input_model/): full-row static input comparison report and
  diagnostics.
- [`publication_3kv/standard_historical_100k/`](publication_3kv/standard_historical_100k/):
  representative voltage, final-phase-space, and waterfall plots from a
  100,000-particle historical sample.
- [`publication_3kv/dynamic_33k6/`](publication_3kv/dynamic_33k6/): matched
  33,600-particle generated/historical comparison and the three-seed aggregate
  screening report.
- [`publication_3kv/final_time_full/`](publication_3kv/final_time_full/):
  final-time comparison of the full 1,075,200-particle historical and generated
  runs, including folded local and unwrapped lineage 125 ns diagnostics.
- [`publication_3kv/gamma_transition_full/`](publication_3kv/gamma_transition_full/):
  controlled full-input comparison of Keegan's transition gamma with Werkema's
  21.6 value, with its own [guide](publication_3kv/gamma_transition_full/README.md).
- [`keegan_5kv/`](keegan_5kv/): legacy-code 5 kV sensitivity convergence study,
  including full-statistics plots/comparisons, the five-seed reduced screen,
  and the 3 kV-versus-5 kV endpoint report.

No HDF5 file is copied here. Regeneration writes those intermediate files to
the ignored `outputs/` tree. It requires the local `1075200_init_dist.txt`,
whose expected SHA-256 is
`3b591b8381cb7cf57a88d7db08d535f2de182720d9fcb3c335800d7b3283e7e1`.

From the project root, verify every listed result artifact with:

```bash
sha256sum --check reference_results/SHA256SUMS
sha256sum --check reference_results/keegan_5kv/SHA256SUMS
```

The command recipes in `manifest.json` are fully explicit workflow recipes
reconstructed from the stored run metadata and project records. They reproduce
the model settings and numerical analysis, not necessarily byte-identical files
after the code or environment changes. The older HDF5 files did not store the
literal shell command, so the recipes are also not transcripts of the original
shell history.
