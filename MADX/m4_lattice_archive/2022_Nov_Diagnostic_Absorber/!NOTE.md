# November 2022 diagnostic-absorber reconstruction

Reconstruction made: **2026-08-15**

## Status

This directory contains the received November 2022 Twiss table and a local
reconstruction of the missing MAD-X source that generated it. The
reconstruction is not a recovered original file and must not be represented as
the authoritative November 2022 source.

The received table is `mu2e-m4da-twiss-2022-11-22.tfs`. It is retained
unchanged in this dated folder, with SHA-256:

```text
0783163ce3e1d80ca3ecbc7dd5dfad88019a3346a905fcf7a7244942981ece80
```

## What is confirmed

The received table directly provides the beam energy, initial Twiss and
dispersion values, 197.645 m sequence length, complete element order,
longitudinal coordinates, hard-edge lengths, bend angles, and integrated
quadrupole strengths at its printed precision.

The complete route through `Q933D` matches the August 2018 source lattice. The
four copied input files below are byte-identical to their `2018_Aug/`
counterparts and retain their 2018 modification times:

- `magnet_eff_lengths.def`
- `m4_new.def`
- `bends_values_aug2018.dat`
- `quads_aug2018.dat`

The diagnostic branch observed in the table uses:

- `HDA1` at -0.0873 rad;
- two 1.524 m halves of `QDA01`, each with `K1L=-0.06653187`;
- two 1.524 m halves of `QDA02`, each with `K1L=+0.08649142`;
- a 1.829 m absorber; and
- no named `MWDA002`, `ABSDRIFT`, or `ABSIN` element.

## What is reconstructed or estimated

`diag_2022_reconstructed.seq` retains the exact August 2018 upstream sequence
through `HDA1`, comments out its production-target tail, and reconstructs the
diagnostic tail from the received table. Thick-element center placements were
derived from the reported exit coordinate and length. The QDA `K1` values in
`job_2022_reconstructed.madx` were derived from the reported `K1L` values and
the 1.524 m half lengths; their magnitudes also match the older stored QDA
definitions, with both signs reversed.

The received Twiss table contains no global survey frame. Consequently, the
two generated survey files are explicitly estimates:

- `survey_estimate_2018_frame.tfs` uses the frame from the August 2018 driver.
- `survey_estimate_2016_diag_frame.tfs` uses the frame from the received 2016
  diagnostic driver.

Neither survey estimate can be validated as the frame used in November 2022
without a received 2022 survey or its missing driver.

## Reproduction result

The reconstruction was run with MAD-X 5.08.00 (64-bit Windows). It exited with
status 0 and reported that MAD-X finished normally. It produced three
warnings: one existing August 2018 zero-angle `V906` length-expression warning
and two identical coupling warnings.

The locally generated `twiss_reconstructed_5.08.tfs` reproduces the received
MAD-X 5.06.00 table extremely closely:

| Comparison | Result |
|---|---|
| Data rows | 331 in each file |
| Element names and order | Identical |
| Element keywords | Identical |
| `S`, `L`, `ANGLE`, and `K1L` | Identical at the printed eight-decimal precision |
| Maximum absolute `BETX` difference | 0.00000008 m |
| Maximum absolute `BETY` difference | 0.00000003 m |
| Maximum absolute `ALFX` difference | 0.00000001 |
| Maximum absolute `DX` difference | 0.00000001 m |
| `DPX`, `DY`, `DPY`, `MUX`, and `MUY` | Identical at the printed precision |

The remaining differences are the original versus local run date/time,
MAD-X 5.06.00 versus 5.08.00 origin string, four extra global header values
written by 5.08.00, and the sub-`1e-7` numerical differences above. This is
strong evidence that the reconstructed physics inputs and sequence topology
match the missing source, while not proving identical source text or digits
beyond those preserved in the received output.

## Files

- `job_2022_reconstructed.madx`: local reconstructed entry point.
- `diag_2022_reconstructed.seq`: local reconstructed sequence.
- `twiss_reconstructed_5.08.tfs`: locally regenerated Twiss table.
- `survey_estimate_2018_frame.tfs`: estimated survey using the 2018 frame.
- `survey_estimate_2016_diag_frame.tfs`: estimated survey using the 2016 frame.
- `reconstruction_run_5.08.log`: complete local MAD-X run log.
- `mu2e-m4da-twiss-2022-11-22.tfs`: byte-identical copy of the received
  reference output.
- The four `.def`/`.dat` files: exact August 2018 dependencies listed above.

From this directory, the reconstruction can be rerun with:

```bash
../../madx.exe < job_2022_reconstructed.madx
```

Rerunning in place will overwrite the locally generated reconstruction
outputs. It will not modify the received reference Twiss table.
