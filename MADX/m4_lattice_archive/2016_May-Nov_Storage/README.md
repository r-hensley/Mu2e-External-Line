# May-November 2016 storage-family lattice

This is a mixed-date working directory for the 2016 `storage` production-target
lattice family, not a clean frozen release. It combines the March 21, 2016
survey/Twiss-fit lineage with later 2016 sequence and strength evolution through
Q943 and the production target. See [`!NOTE.md`](!NOTE.md) for the separately
documented diagnostic-line repair.

Eric Prebys's recovered G4Beamline study contains a related `MAD_Eliana`
working package and generated G4Beamline fragments:
[`G4beamline/eric_g4beamline_study/source/lattice/MAD_Eliana`](../../../G4beamline/eric_g4beamline_study/source/lattice/MAD_Eliana/README.md).

The packages are demonstrably related but are not globally byte-identical:

- `m4.def` is byte-for-byte identical.
- The G4Beamline package's `values.dat` is byte-for-byte identical to this
  directory's `values_May_2016.dat`.
- Their `line_shifted_new.seq` and `m4.madx` files are distinct working variants.

Later file, rerun, export, and Git dates do not make either package a newer
lattice generation. Neither directory should be described as current,
official, as-built, or operationally approved without independent evidence.
