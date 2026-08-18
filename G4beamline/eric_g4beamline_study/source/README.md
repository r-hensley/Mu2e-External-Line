# Simulation sources

This directory separates the historical study's executable inputs by role:

- [`lattice/`](lattice/README.md) contains magnet definitions, the translated
  Eliana MAD-X lattice fragments, and the small embedded MAD-X archive.
- [`scan_decks/`](scan_decks/README.md) contains parameterized G4Beamline decks
  for the standard and component-removal/absorber scans.
- [`analysis/`](analysis/README.md) contains ROOT macros and the deck-generation
  helper used to prepare distributions and analyze transmission.

These are preserved working files, not a claim that the historical study runs
unchanged with current ROOT, G4Beamline, or grid software.
