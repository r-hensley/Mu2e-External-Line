# Eric Prebys G4Beamline study archive

This directory is a content-deduplicated, documented copy of Eric Prebys's
Mu2e G4Beamline studies. It preserves the recoverable source decks, lattice
translations, ROOT analysis macros, small reference results, original notes,
and local copies of the large data files. The material spans several working
periods rather than one frozen software release.

Start with Eric's [original study notes](historical_docs/originals/AAA_g4beamline_study.html)
or the shorter [documentation guide](historical_docs/README.md). The associated
2014 transmission study is preserved as [DocDB 4054 v4](historical_docs/originals/prebys_docdb-4054-v4%20Beam%20Transmission%20and%20Extinction%20in%20the%20Mu2e%20Beam%20Line.pdf).

## What is here

| Directory | Contents |
|---|---|
| [`source/`](source/README.md) | G4Beamline decks, translated lattice fragments, ROOT macros, and scan generators |
| [`conversion/`](conversion/README.md) | Eric's `mad2g4bl.py` converter and a concise guide |
| [`data/`](data/README.md) | Large particle inputs retained locally and ignored by Git |
| [`reference_results/`](reference_results/README.md) | Tables, plots, and selected historical outputs |
| [`historical_docs/`](historical_docs/README.md) | Eric's original HTML notes, the DocDB paper, and email provenance |
| [`historical_snapshots/`](historical_snapshots/README.md) | Unique files from the January 2017 archive that differ from the later export |
| [`provenance/`](provenance/README.md) | Complete source-to-destination, SHA-256, timestamp, and deduplication records |

## Lattice generation

The [`MAD_Eliana`](source/lattice/MAD_Eliana/README.md) files belong to the
May-November 2016 `storage` production-target lattice family. They are not the
August 2018 reference lattice, the October 2024 study, or a 2026 lattice. Some
tables and drivers have 2025-2026 rerun/export timestamps, which record later
processing rather than the represented lattice generation.

## Reproducibility status

The archive contains the important ingredients for reconstructing Eric's later
2015-2017 workflow: source decks, generated lattice segments, analysis macros,
scan decks, and particle inputs. It does not freeze the historical executable
stack or grid environment, so an exact bit-for-bit reproduction of old outputs
is not established. A current local installation provides ROOT 6.32.06 and
G4Beamline 3.08 (Geant4 11.0.3), but no production rerun has been performed as
part of this import.

Large files are intentionally ignored by the folder-local [`.gitignore`](.gitignore).
Their READMEs link back to Eric's server, and the provenance manifest records
the expected byte count and SHA-256 digest for verification after download.
Run [`data/download_large_files.sh`](data/download_large_files.sh) to restore
all nine ignored files without unpacking any other archive content.

## Import summary

The import examined 308 source-path occurrences from the later export, the
January 2017 archive, the server's standalone data directory, and the selected
documentation. It retained 285 unique SHA-256 objects and omitted 23 exact
duplicate occurrences. Originals in the recovery workspace were not moved or
modified; copied historical files retained their source modification times.
