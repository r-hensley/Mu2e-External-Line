# Large particle inputs (not committed)

The files in this directory are ignored by Git. Eric's server still exposes the
standalone study data directory [here](https://prebys.physics.ucdavis.edu/misc/AAAreadme/g4beamline_study/).
The later export tar file, which contains the two `.dat` files, can be downloaded
[here](https://prebys.physics.ucdavis.edu/misc/AAAreadme/g4beamline_study/g4beamline_study_export.tgz).

| File | Bytes | SHA-256 | Server source |
|---|---:|---|---|
| `ExtrBeamDSC-magL12noDiff.dat` | 19,262,735 | `0136b4eb07135ca7c0147fa2dc2d26b5c7728e79229327652e375f1e50f3d94b` | [export tar](https://prebys.physics.ucdavis.edu/misc/AAAreadme/g4beamline_study/g4beamline_study_export.tgz) |
| `ExtrBeamTransToCMAG-D.dat` | 49,465,532 | `cf7c836f923dd814f87ca7a6645edff3a02c0c9c05bcf021186002ed9e6b8149` | [download](https://prebys.physics.ucdavis.edu/misc/AAAreadme/g4beamline_study/ExtrBeamTransToCMAG-D.dat) |
| `exttracks_rot_1800000.root` | 45,192,824 | `7df4c14fcb651c9259293b33ca9e5ae645d9ad8a6b61a13095af64c36f055edb` | [download](https://prebys.physics.ucdavis.edu/misc/AAAreadme/g4beamline_study/exttracks_rot_1800000.root) |
| `exttracks_rot_e30_e40_1800000.root` | 13,577,449 | `583286959bb1cdd94df638a62c52c10412e8e3f6a852d0cc9bb4a2fbf130cf74` | [download](https://prebys.physics.ucdavis.edu/misc/AAAreadme/g4beamline_study/exttracks_rot_e30_e40_1800000.root) |
| `mu2e_downstream.root` | 17,985,327 | `7e5cefcb2fad3f2eb615dbfa4e9d47a804a717ee541f6d2eedef96cd748ffc0a` | [download](https://prebys.physics.ucdavis.edu/misc/AAAreadme/g4beamline_study/mu2e_downstream.root) |
| `mu2e_upstream.root` | 18,010,516 | `3853c89068ffd225a358b52309d87f812cf8d8034ecddd8fade34ff75f848968` | [download](https://prebys.physics.ucdavis.edu/misc/AAAreadme/g4beamline_study/mu2e_upstream.root) |
| `mu2e_upstream_gaussian.root` | 16,404,038 | `50404b83ff69620be6bfdb06ab0ba8f4d4e403bddb38714461a1efd5f73c4196` | [download](https://prebys.physics.ucdavis.edu/misc/AAAreadme/g4beamline_study/mu2e_upstream_gaussian.root) |

After downloading a file, compare it with the digest above or the complete
[`../../provenance/import_manifest.tsv`](../../provenance/import_manifest.tsv).

For an automated, verified restoration of these inputs and the two ignored
reference PDFs, run [`../download_large_files.sh`](../download_large_files.sh).
It extracts only the required large members from the export archive and
downloads the five ROOT inputs directly because they are not inside either
server tarball.
