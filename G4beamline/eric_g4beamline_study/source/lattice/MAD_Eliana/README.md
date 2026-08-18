# `MAD_Eliana` lattice generation

This is a mixed working model from the **May-November 2016 `storage`
production-target lattice family**. Its starting survey/Twiss fit traces to
March 21, 2016, while `line_shifted_new.seq` and the setting files reflect the
later 2016 evolution through Q943 and the production target.

It is not the August 2018 reference lattice, the October 2024 study, or a 2026
lattice. The 2025-2026 dates on `m4.madx`, `m4.log`, `survey_new.tfs`, and
`twiss_new.tfs` are later reruns/exports and do not change the represented
lattice generation.

The related MAD-X archive is documented at
[`2016_May-Nov_Storage`](../../../../../MADX/m4_lattice_archive/2016_May-Nov_Storage/README.md).
The relationship is supported directly by content:

- `m4.def` is byte-for-byte identical in both locations.
- This directory's `values.dat` is byte-for-byte identical to the archive's
  `values_May_2016.dat`.
- `line_shifted_new.seq` and `m4.madx` are related working variants, not exact
  duplicates.

The later export also contains generated G4Beamline fragments for the upstream
and downstream parts of the line. Four older MAD-X table/driver variants are
kept separately in the [January 2017 snapshot](../../../historical_snapshots/2017_01_03/README.md)
because their contents differ from this export.
