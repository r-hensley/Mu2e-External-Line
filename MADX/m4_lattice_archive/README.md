# Mu2e M4 MAD-X lattice archive

This archive contains Mu2e M4 MAD-X lattice packages and related Delivery Ring
files. Most M4 packages follow the route to the production target. Directory
dates identify the approximate lattice or study generation; they do not imply
that a package is current or operationally approved.

## Directory guide

| Directory | Contents |
|---|---|
| `2015_May/` | Early complete production-target lattice and a MAD-X optics plot |
| `2015_Oct/` | Later 2015 production-target iteration, including its original ZIP |
| `2016_Mar_11/` | Production sequence plus a diagnostic sequence; the production driver is missing `quads_final.dat` |
| `2016_Apr/` | Runnable production-target survey/Twiss-fitted iteration |
| `2016_May-Nov_Storage/` | Mixed-date production working directory with historical diagnostic-line source files, a related recovered Prebys G4Beamline package, and a documented runnable repair; see its `README.md` and `!NOTE.md` |
| `2017_Jan-May/` | Surviving mixed working set using 2017-era files and two shared August 2018 definitions; see `!WARNING.txt` |
| `2018_Aug/` | Complete August 2018 production-target reference package, including its original ZIP |
| `2022_Nov_Diagnostic_Absorber/` | Received November 2022 diagnostic Twiss table, a local source reconstruction, regenerated Twiss output, and two explicitly estimated surveys; see `!NOTE.md` |
| `2024_Oct/` | October 2024 TIEIN/steering study, a DUSAF target-marker variant, and a November 2024 mechanical survey file |
| `2025_Aug_Mad4Tail/` | Q937U-to-target slice of the 2024 lattice with matching and survey experiments; see `!NOTE.txt` |
| [2026_09_29_Deinlein_DR_M4/](2026_09_29_Deinlein_DR_M4/README.md) | September 2026 Delivery Ring and Q303-to-diagnostic-absorber source models, each with supplied Twiss tables and plots; includes the common M4 trunk through Q933 but no production-target continuation. Overview in its `README.md`; full evidence in its `audit/AUDIT.md`. |
| `Undated_mACE_Diagnostic_Absorber/` | Undated mACE tight-focus optics table and a local runnable first-order reconstruction; the best working date estimate is circa 2025; see `!NOTE.md` |
| `manifests/` | Historical merge snapshots with SHA-256 checksums and tab-separated inventories; later additions have separate checksum inventories, as explained in `manifests/README.txt` |

## Diagnostic-absorber models and outputs

| File | Contents |
|---|---|
| `2022_Nov_Diagnostic_Absorber/mu2e-m4da-twiss-2022-11-22.tfs` | Received MAD-X 5.06.00 Twiss output for the M4 diagnostic-absorber route. It follows the August 2018 M4 trunk through Q933, then uses HDA1, QDA01, QDA02, and the absorber. It is retained in the dated reconstruction folder; the original generating source files remain unavailable. |
| `Undated_mACE_Diagnostic_Absorber/mu2e-dabs-v2-functions.csv` | Undated optics-function export for a modified version of the same diagnostic layout, intended to reduce the beam size near the absorber for mACE. It directly contains optics and geometry rather than magnet strengths. The accompanying local reconstruction converts all 692 consecutive-row maps into a runnable MAD-X first-order transport model and reproduces every one of the 693 samples at the CSV's six-decimal precision. This does not recover a unique physical magnet lattice. |
| [2026_09_29_Deinlein_DR_M4/mu2e-dabs-lattice-v2026.09.29.madx](2026_09_29_Deinlein_DR_M4/mu2e-dabs-lattice-v2026.09.29.madx) | Supplied physical-element source model from Q303 through the extraction region and M4 to the diagnostic absorber; 226.137 m. Its revised settings differ from the older diagnostic outputs. The folder also contains a separate complete Delivery Ring source and supplied TFS/PDF files for both models; see the [collection overview](2026_09_29_Deinlein_DR_M4/README.md) and [audit](2026_09_29_Deinlein_DR_M4/audit/AUDIT.md). |

## Important relationships

- The eight files in `2018_Aug/` are byte-identical to the August 2018 MAD-X
  inputs in the Mu2e-External-Line repository.
- `2024_Oct/for_Igor/` is a complete runnable package that was independently
  received from two sources with identical files.
- `2024_Oct/Eliana_DUSAF_Target_Variant/` preserves the same October 2024
  magnet lattice through V943 but changes the terminal target-marker
  bookkeeping. Its `tien_shift` reference appears to be a typo; see its
  `!WARNING.txt`.
- `2025_Aug_Mad4Tail/m4tail.seq` is an extracted tail of the 2024 Igor lattice,
  not a newer complete M4 lattice.
- The received 2022 Twiss table and undated mACE CSV are related geometrically
  but are not the same optics solution. The CSV is the tighter-focus variant. Its sampled
  first-order optics can now be reproduced exactly, although its original
  source and unique hardware-level magnet representation remain unknown.
- `2022_Nov_Diagnostic_Absorber/` reconstructs the missing source from the
  exact August 2018 upstream inputs and the received 2022 table. A MAD-X 5.08.00
  rerun reproduces all 331 element names, types, coordinates, lengths, angles,
  and strengths at printed precision; the largest beta difference is 8e-8 m.
  This is strong numerical reconstruction evidence, not a recovered original.
- The September 2026 DR and extraction jobs share a Q303 starting point and
  upstream quadrupole settings, but run independently. The extraction job uses
  explicit starting optics approximating the ring periodic solution; neither
  job simulates resonant extraction or transfers particles between files.
  The author described the extraction settings as proposed future corrections.
  All 33 common M4 quadrupole strengths differ from the preserved 2018/2024
  values, while several bend values agree with the 2024 study. This establishes
  a revised model, with operational adoption and exact ancestry unresolved.
  A production-target model requires a checked join before the diagnostic bend.
- Files named `!WARNING.txt`, `!NOTE.txt`, or `!NOTE.md`, this guide, and the manifests are
  explanatory documentation rather than received lattice inputs.
- `2016_May-Nov_Storage/diag_line_repaired_2026-08-15.madx` is a local repaired
  copy, not a received lattice. It adds the historically supported 0.76831 m
  `ABSDRIFT` and 1.829 m absorber pair and runs normally with MAD-X 5.08.00;
  the received `diag_line.madx` remains unchanged.
- The recovered Prebys `MAD_Eliana` G4Beamline working package belongs to the
  same May-November 2016 `storage` family. Its `m4.def` is identical here, and
  its `values.dat` is identical to `values_May_2016.dat`; its sequence and
  driver are related but distinct variants.

The current archive-wide inventories include every file outside `manifests/`,
including the September 2026 collection and this guide. The pre-merge manifests
remain historical provenance records. From this directory, verify the archive:

```bash
sha256sum -c manifests/m4_lattice_archive.sha256
```

The September collection also has a portable `SHA256SUMS` list. See
[manifests/README.txt](manifests/README.txt) for the inventory fields.

Archive text files use Linux (LF) line endings. The original and normalized
hashes for converted files are recorded in
[manifests/line_ending_normalization.tsv](manifests/line_ending_normalization.tsv).
The conversion changes line terminators only; lattice statements and numeric
values are unchanged. Binary files retain their received bytes.

The broader chronology and model relationships are documented in
[madx_file_provenance.md](../../madx_file_provenance.md). Full evidence for the
new models is in their [audit document](2026_09_29_Deinlein_DR_M4/audit/AUDIT.md).
