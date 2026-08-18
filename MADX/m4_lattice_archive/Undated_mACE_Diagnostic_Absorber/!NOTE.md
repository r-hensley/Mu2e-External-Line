# Undated mACE diagnostic-absorber optics

Archive organization and reconstruction work performed: **2026-08-15**

## Date estimate

The received CSV contains no embedded creation date, author, MAD-X version, or
source-file reference. Its original generation year is therefore **not
confirmed**.

The best working estimate is **circa 2025**, now with moderate rather than low
confidence. A reasonable inferred window is **late 2024 through June 2026**.
The accompanying correspondence establishes only that Ryan Hensley had
received the CSV from George by June 2026, but the recoverable first-order bend
maps provide an additional chronology clue: the inferred H910-family/H912-
family bend-angle ratio is `0.784544`, extremely close to the `0.784244` ratio
in Eliana's October 2024 study and unlike the August 2018 ratio `0.805195`.
This strongly suggests that the export incorporates the left-bend correction
identified in 2024, although it is still an inference rather than an embedded
file date. The CSV's filesystem modification time is a later attachment/local-
copy timestamp and does not date the lattice.

## Received file

`mu2e-dabs-v2-functions.csv` is preserved byte-for-byte with its pre-move
modification timestamp:

```text
size:    86704 bytes
mtime:   2026-08-13 12:04:02.951287300 -0700
SHA-256: bcd4283d2800673eb85c2ebe9c35c69655783761a0652538fd7e31d32ddec2aa
```

The file is a 693-row optics-function export for the route from DQ303 through
the M4 diagnostic absorber. It gives beginning, midpoint, and endpoint values
for each element but does not directly provide magnet strengths, bend angles,
beam energy, or a generating program/version.

## Reconstruction method

The additional MAD-X files in this directory are a local first-order
reconstruction, not recovered historical source. The CSV's beta, alpha, phase,
and dispersion values uniquely determine an uncoupled first-order transfer map
between each consecutive pair of rows. The generator converts all 692
row-to-row intervals into MAD-X `MATRIX` elements. This includes zero-length
maps between one element's `END` row and the next element's `BEG` row, which
are necessary to retain the table's bend-edge convention.

The driver uses a proton total energy of `8.93828 GeV`, borrowed from the
surviving M4 MAD-X jobs because the CSV does not state a beam energy. All
geometry, initial Twiss/dispersion values, and transfer maps come from the CSV.

The same samples allow approximate hard-edge parameters to be inferred. The
report contains estimates for all 43 quadrupoles and 23 dipoles. Examples are:

- QDA01: `K1 = -0.0665381337 m^-2` over `3.048 m`;
- QDA02: `K1 = +0.1130198610 m^-2` over `3.048 m`;
- HDA1: inferred angle `-0.0872992590 rad`; and
- H910- and H912-family averages `-0.1085321674 rad` and
  `-0.1383379339 rad`, respectively.

These estimates describe first-order maps derived from rounded values. They
are not recovered power-supply settings, field maps, or unique physical
magnet definitions.

## MAD-X reproduction result

The reconstruction was run in an isolated temporary directory with MAD-X
`5.08.00` for 64-bit Windows. It exited normally with status 0 and **zero
warnings**. The output has 1,387 rows: 693 named CSV sample markers, 692
intervening `MATRIX` elements, and the MAD-X sequence endpoints.

All 693 source sample names were found. Every compared `S`, beta, alpha,
dispersion, dispersion-derivative, and phase value matches at the CSV's
six-decimal precision. Before rounding, the largest observed differences are:

| Quantity | Maximum absolute difference |
|---|---:|
| `BETX`, `BETY`, `S`, `DPX`, `MUX`, `MUY` | `0` at stored output precision |
| `ALFX` | `1.11e-14` |
| `ALFY` | `1.00e-14` |
| `DX` | `1.00e-14 m` |
| `DY` | `1.00e-15 m` |
| `DPY` | `3.00e-17` |

## Scope and limitations

This is a complete runnable **first-order transport reconstruction** from DQ303
through the diagnostic absorber. It is useful for reproducing the supplied
Twiss functions and for linear transport studies. It is not a unique physical
MAD-X magnet lattice: the CSV does not preserve a survey frame, apertures,
materials, nonlinear fields, exact magnet/edge definitions, power-supply
values, or the generating program and version. A faithful 3D geometry or
hardware-level model still requires the missing original source.

## Files

- `mu2e-dabs-v2-functions.csv`: byte-preserved received source table.
- `reconstruct_from_csv.mjs`: validates the table and regenerates the MAD-X
  definitions, sequence, driver, and parameter report.
- `mace_transport_maps.def`: 692 inferred first-order transfer maps.
- `mace_transport_reconstructed.seq`: the complete DQ303-to-absorber sequence.
- `job_mace_transport_reconstructed.madx`: runnable MAD-X entry point.
- `inferred_first_order_parameters.txt`: approximate quadrupole K1 and dipole
  angle/tilt estimates.
- `twiss_mace_transport_reconstructed.tfs`: MAD-X 5.08.00 regenerated output.
- `compare_csv_to_tfs.mjs`: reusable named-marker comparison script.
- `reconstruction_comparison.txt`: saved numerical comparison summary.
- `reconstruction_run_5.08.log`: complete successful MAD-X run log.

From this directory, rerun with:

```bash
../../madx.exe < job_mace_transport_reconstructed.madx
node compare_csv_to_tfs.mjs twiss_mace_transport_reconstructed.tfs
```

Rerunning in place overwrites only the locally generated TFS output, not the
received CSV.
