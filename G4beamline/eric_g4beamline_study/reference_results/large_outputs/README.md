# Large reference outputs (not committed)

These phase-space PDFs are retained locally but ignored by Git. Both are inside
Eric's study export tar file, available [here](https://prebys.physics.ucdavis.edu/misc/AAAreadme/g4beamline_study/g4beamline_study_export.tgz).

| File | Bytes | SHA-256 |
|---|---:|---|
| `us_phase_space_col.pdf` | 29,988,789 | `27d011577131a77d5e5e9066259773b2fbacbc4e8e9045f29c51ac026a28f2fa` |
| `us_phase_space_nocol.pdf` | 30,027,796 | `96112a39ea8db2427cef5e894cf3a88e38eed66210659dd1103a90c4779b64e9` |

The expected digests are also recorded in
[`../../provenance/import_manifest.tsv`](../../provenance/import_manifest.tsv).
They can be restored together with the ignored particle inputs by running
[`../../data/download_large_files.sh`](../../data/download_large_files.sh).
