# MAD-X to G4Beamline conversion

[`mad2g4bl.py`](mad2g4bl.py) is Eric Prebys's converter for producing
G4Beamline fragments from a MAD/MAD-X optics table. The preserved March 2020
revision supports Python 3 and `TWISS,CENTRE` output.

Read Eric's [original converter notes](../historical_docs/originals/AAAREADME_mad2g4bl.html)
for the detailed input format and transformation conventions. In outline, the
script reads tabular element positions and orientations, maps recognized MAD-X
element types to G4Beamline geometry, handles coordinate and bend-center
conventions, and divides output into upstream/downstream segment files.

The generated files preserved under [`../source/lattice/MAD_Eliana/`](../source/lattice/MAD_Eliana/README.md)
are historical products. Regenerating them requires confirming the expected
input-table columns and local configuration in the script; do not assume that
a modern MAD-X table can be substituted without checking conventions.
