# Mu2e External M4 Beam-Line Models

This repository archives reference lattice and particle-tracking models of the
Mu2e M4 external beam line along the route to the production target. It includes
implementations for MAD-X, G4Beamline, and BDSIM.

> [!IMPORTANT]
> The files in this repository are research reference models. They should not
> be assumed to represent the currently approved, as-built, or as-operated beam
> line without confirmation from the responsible lattice and accelerator teams.

The tracked MAD-X package represents the August 2018 production-target lattice.
The files were added to Git in 2026; those commit dates record repository
publication, not a 2026 lattice update.

## Published Model Snapshots

| Model | Reference entry point | Basis and role |
|---|---|---|
| MAD-X | [`MADX/job_aug2018.madx`](MADX/job_aug2018.madx) | Complete August 2018 optics and survey model for M4 to the production target |
| G4Beamline | [`G4beamline/G4_M4_Mu2e_03.g4bl`](G4beamline/G4_M4_Mu2e_03.g4bl) | Tracking model based on the August 2018 lattice, with model-specific revisions through version 03 |
| BDSIM | [`BDSIM/M4_Mu2e_BDSIM_v1.gmad`](BDSIM/M4_Mu2e_BDSIM_v1.gmad) | BDSIM conversion derived from the G4Beamline version 03 deck |

The documentation labels `madx-v2018.08`, `g4bl-v0.3`, and `bdsim-v1.0`
identify these snapshots. They are model labels, not claims about the current
machine configuration. Git release tags have not yet been published for them.

All three tracked models follow the production-target route:

```text
MSTART ... Q933 -> zero-angle HDA1 -> Q934 ... Q943 -> TIEIN -> production target
```

## August 2018 MAD-X Reference

The MAD-X sequence is 244.900 m long and is associated with Dr. Eliana
Gianfelice-Wendt. Comparison with an independently preserved August 2018 source
package found all eight tracked MAD-X inputs to be byte-for-byte identical.
This anchors the repository model to the August 2018 lattice rather than to its
later Git publication date.

| File | Purpose |
|---|---|
| [`job_aug2018.madx`](MADX/job_aug2018.madx) | Main MAD-X driver; selects the sequence and writes survey and Twiss tables |
| [`m4_aug2018a.seq`](MADX/m4_aug2018a.seq) | Longitudinal element placements for the production-target sequence |
| [`m4_new.def`](MADX/m4_new.def) | Element and magnet definitions |
| [`magnet_eff_lengths.def`](MADX/magnet_eff_lengths.def) | Effective magnetic lengths used for optics |
| [`magnet_mech_lengths.def`](MADX/magnet_mech_lengths.def) | Mechanical lengths used for survey comparisons |
| [`bends_values_aug2018.dat`](MADX/bends_values_aug2018.dat) | August 2018 bend settings |
| [`quads_aug2018.dat`](MADX/quads_aug2018.dat) | Main-line quadrupole settings |
| [`FF_quads_aug2018.dat`](MADX/FF_quads_aug2018.dat) | Final-focus quadrupole settings |

Reference Twiss and survey products are stored under
[`MADX/MAD Output/`](MADX/MAD%20Output/).

With MAD-X installed, run the reference model from the `MADX` directory so its
relative input paths resolve correctly:

```bash
cd MADX
madx < job_aug2018.madx
```

The job writes `m4_survey.tfs` and `m4_twiss.tfs`. The archived August 2018
outputs identify their generator as MAD-X 5.02.00 for 64-bit Darwin. The
supplied source package also runs successfully with MAD-X 5.08.00.

## Relationship Between the Three Implementations

The three models share a common production-lattice basis, but they are not
interchangeable file-for-file translations:

- The G4Beamline file records that version 01 was based on the August 2018
  MAD-X lattice.
- G4Beamline version 02 added code-specific AC-dipole, timing, ordering, and
  tracking changes in 2019.
- G4Beamline version 03 added the first four trims in February 2022.
- The BDSIM file identifies itself as a conversion of a shortened G4Beamline
  version 03 deck and uses BDSIM-specific bend, roll, aperture, and unit
  conventions.

### Earlier Prebys G4Beamline studies

Eric Prebys's [DocDB 4054 v4](https://mu2e-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=4054),
dated 2014-05-07, documents an earlier G4Beamline transmission/extinction
study covering only the AC-dipole-to-production-target section. It describes a
then-current design with a collimator at 928 and records that its MAD lattice
conversion used a Perl script written and run by Jean-Francois Ostiguy. This is
historical study evidence, not the generating chain for the repository's
August 2018-based G4Beamline model.

Prebys's later [G4Beamline study notes](https://prebys.physics.ucdavis.edu/misc/AAAreadme/AAA_g4beamline_study.html)
document a different workflow. His Python `mad2g4bl.py` translated Eliana
Gianfelice-Wendt's MAD-X optics table into G4Beamline fragments. A separate
ROOT macro, `DRGen.C`, converted Vladimir Nagaslaev's particle distribution at
the downstream end of the C-magnet and optionally rematched only its Y phase
space to Eliana's optics. The Y-only operation applies to the input-particle
distribution, not to the lattice geometry. Prebys has confirmed the notes as
his documentation and authorized incorporation of relevant information into
this repository.

A curated copy is now preserved under
[`G4beamline/eric_g4beamline_study/`](G4beamline/eric_g4beamline_study/README.md).
It includes the original documentation, converter, G4Beamline scan decks,
analysis macros, translated lattice fragments, small historical results, and
SHA-256 provenance maps. Large particle inputs and two large result PDFs are
documented but ignored by Git; their READMEs link to Prebys's server and give
the expected sizes and hashes. The accompanying
[`download_large_files.sh`](G4beamline/eric_g4beamline_study/data/download_large_files.sh)
restores only those nine ignored dependencies: four selected export-archive
members and five standalone ROOT inputs. It verifies the server archives and
every installed file by SHA-256 and refuses to overwrite conflicting data.

The recovered `MAD_Eliana` package belongs to the May-November 2016 `storage`
production-target family, not the August 2018 or October 2024 lattice. Its
`m4.def` is byte-identical to the corresponding
[`2016_May-Nov_Storage`](MADX/m4_lattice_archive/2016_May-Nov_Storage/README.md)
archive file, and its `values.dat` is byte-identical to that archive's
`values_May_2016.dat`. The sequence and driver are related working variants,
not exact duplicates. Later rerun/export timestamps do not redefine the
lattice generation.

## Production-Lattice Chronology

Historical source packages were reviewed to establish the lineage below. The
historical MAD-X archive and recovered Prebys study are now included as
comparison material; they are not the repository's published August 2018
reference snapshot.

| Period | Production-target model | Relationship to the published reference |
|---|---|---|
| May 2015 | Complete early M4 model | Includes earlier quadrupole, corrector, and collimator choices |
| October 2015 | Complete optics and placement iteration | Revises several placements and the left-bend-region tune |
| March 2016 | Incomplete packaged production job | Preserves a production sequence but is missing its selected strength file |
| April 2016 | Complete survey/Twiss-fitted model | Introduces the initial Twiss conditions retained by later reference jobs |
| Late 2016 | Complete working production model | Includes additional TIEIN and target-path development |
| 2017–early 2018 | Rematching and alternative setting files | Survives as a mixed working collection rather than a clean release |
| August 2018 | Complete production reference | Exact source package published in this repository |
| October 2024 | TIEIN and trajectory study | Retains the 2018 quadrupoles while revising selected steering settings |
| August 2025 | Q937U-to-target study | Exact extraction of the 2024 production tail, not a newer complete M4 lattice |

Directory names, file timestamps, and Git commits are useful provenance clues,
but they do not by themselves establish that a lattice was approved or used in
operations.

## Known Later Production Studies

The separately preserved October 2024 study package is a deliberate
post-August-2018 configuration, not merely a repackaging of the repository
files. The comparison found:

- Two independently preserved copies of the five-file Igor study package are
  byte-for-byte identical, confirming a common artifact in circulation.
- All 43 common quadrupole strengths remain at their August 2018 values.
- Most principal magnet centers are unchanged.
- The H910/H911/H917/H918 and H912/H916 left-bend settings are revised.
- H936A, V936, and V943 are revised for the downstream trajectory.
- In the Igor package, the later TIEIN and target markers are 10.497 mm
  downstream of the older markers.
- A separate standalone sequence keeps the same magnet lattice through V943
  but replaces the two target markers with one DUSAF-study marker 5.579521 mm
  upstream of the Igor target. Its apparent downstream-shift variable is
  misspelled and inactive, so a reproduced run changes only target-marker
  bookkeeping, not the magnet optics.

These changes were developed for trajectory and TIEIN studies. The available
evidence does not establish which settings, if any, were ultimately adopted for
operation. The August 2025 tail package reproduces the Q937U-to-target portion
of this 2024 study exactly; it is not an independent full-line revision. Its
first direct attempt to constrain the geometric `Z` coordinate uses a keyword
that standard MAD-X does not support. A later same-day driver instead runs
`SURVEY` inside `MATCH, USE_MACRO`, reads the TIEIN `X`, `Y`, and `Z` table
values, and constrains those expressions. That macro-based driver is the
compatible form; the package contains no supplied output identifying which
MAD-X build its author used. These matching-driver changes do not alter the
tail sequence or any magnet setting.

Neither the G4Beamline nor BDSIM model currently carries these identified 2024
steering changes.

## Status and Contributions

The repository currently publishes one reference snapshot for each simulation
code. A future model update should document:

- The represented lattice/configuration date.
- Whether it follows the production-target route or another branch.
- The source and rationale for geometry or magnet-setting changes.
- Code-specific approximations, added physics, apertures, and scoring elements.
- A reproducible entry point and representative validation outputs.
- The review or approval basis for any claim that a model is current or
  operational.

Model updates should keep the MAD-X, G4Beamline, and BDSIM version histories
independent. A code-specific change does not automatically create a new lattice
version in the other implementations.

## Credits

- August 2018 MAD-X reference package: Dr. Eliana Gianfelice-Wendt
- G4Beamline model: Diktys Stratakis
- BDSIM model: Diktys Stratakis

Additional development history is preserved in source comments and the Git
commit history.

## Diagnostic-Absorber Material

The published models in this repository do not contain a complete runnable
diagnostic-absorber branch. The August 2018 MAD-X driver contains only a
commented call to a missing `diag.seq`, together with supporting element
definitions.

Separate provenance work identified the following comparison material, which
is not currently part of the published model package:

| Material | Finding |
|---|---|
| Historical 2016 `diag.seq` files | Preserve complete source-level diagnostic placements, but use older upstream settings and do not reproduce the later output |
| November 2022 diagnostic Twiss table | Matches the August 2018 M4 trunk through Q933, then follows HDA1, QDA01, QDA02, and the absorber. Its generating source job remains missing, but a separately archived local reconstruction using the 2018 inputs reproduces all 331 rows' identities, coordinates, lengths, angles, and strengths at printed precision, with optics differences below `1e-7 m`; that reconstruction is not a recovered or authoritative source file. |
| Modified mACE optics table | Closely matches the 2022 diagnostic geometry after its coordinate offset, but represents substantially different tight-focus optics. It contains no explicit magnet strengths or source date; the best working date estimate is circa 2025 because its inferred H910/H912 bend ratio closely matches the October 2024 study rather than August 2018. A separately archived 692-map MAD-X reconstruction reproduces all 693 table samples at six-decimal precision, but is a first-order transport equivalent rather than a recovered or unique physical magnet lattice. |

The October 2024 production study contains no QDA01/QDA02 or absorber branch.
The virtual scoring planes named `VD_Diagnostic_*` in G4Beamline and BDSIM are
not implementations of this diagnostic-absorber line.
