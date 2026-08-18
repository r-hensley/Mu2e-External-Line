# Model provenance

This note records which parts of the reconstruction come from Keegan Harrig's
committed BLonD project, which come from the publication, and which are modeling
choices made in this Xsuite port. It is intentionally short enough to ship with
the code.

## Historical source material

- The audited BLonD project is identified by commit
  `3566338e5e1e201028f2ec3ec62a94dcf1fc0031` (`reorg file structure`,
  2023-01-13). Its custom simulation drivers were substantively written in
  2019--2020 and moved without content changes in 2023.
- `BLonD/multi_mpi.py` is the committed RF-only production-shaped driver. It
  uses Recycler harmonic 28, zero RF phase, a 5 kV to 80 kV inverse-square
  voltage ramp lasting 85 ms, and then holds 80 kV. It does not contain the
  preceding harmonic-588 turn-off stage.
- The 5 ms harmonic-588 turn-off and the sequential timing of the two systems
  come from the paper/dissertation and were explicitly confirmed for this
  reconstruction. The publication plots label the harmonic-28 starting voltage
  as 3 kV. Werkema's documented ESME input deck likewise specifies a 0.003 MV
  initial harmonic-28 voltage at its precise 5.012 ms program boundary.
- The Xsuite standard is therefore locked to the sequential 5 ms h=588
  turn-off followed by the 85 ms h=28 inverse-square ramp from 3 kV to 80 kV.
  Zero phase still follows the committed driver. Its 5 kV start remains
  available only as the explicitly labeled legacy-code sensitivity
  `--rf25-start-kv 5`; the durations and 80 kV endpoint are not scan parameters
  in the standard workflow.

## Transition-gamma source comparison

Keegan's committed BLonD machine definition uses
`gamma_t = 20.257643837730637`; this remains the Xsuite default. Werkema's
documented model uses `gamma_t = 21.6`. The controlled CLI override
`--gamma-transition 21.6` changes only transition gamma as an independent
input. Momentum compaction and slip factor are then derived rather than tuned
independently:

| Source value | `gamma_t` | `alpha_c` | `eta` | Small-amplitude h=28 period at 80 kV |
|---|---:|---:|---:|---:|
| Keegan | 20.257643837730637 | 0.0024368126329700006 | -0.008714928548793705 | 18.72423519130096 ms |
| Werkema | 21.6 | 0.002143347050754458 | -0.009008394131009248 | 18.41672119458596 ms |

The period values are derived small-amplitude estimates at the fixed final
80 kV, not results of a transition-gamma scan. The corresponding comparison
mode requires the same complete file-backed input and matching RF program,
endpoint, profile grid, and population scaling; it permits only `gamma_t` and
its derived `alpha_c` and `eta` metadata to differ.

## Particle input

The historical input is [`1075200_init_dist.txt`](1075200_init_dist.txt), a
1,075,200-row, two-column text file. Its SHA-256 is
`3b591b8381cb7cf57a88d7db08d535f2de182720d9fcb3c335800d7b3283e7e1`.
The committed MPI drivers interpret the columns as longitudinal position in
metres and fractional momentum offset `dp/p`; this port preserves their
coordinate and energy conversions.

The in-memory generator was fitted to the geometry and within-microbunch moments
of that file. It creates 168 equally populated Gaussian microbunches in two
trains and assigns each consecutive group of 21 to one of eight structural
bunches. It does not insert ghost bunches, satellite particles, special tails,
or per-microbunch offsets. Agreement with the historical file is therefore a
calibration result, not proof of Keegan's unavailable original generator or
random seed.

## Scope of the port

The current model contains two time-programmed RF cavities followed by the
Recycler nonlinear longitudinal slip map. It contains no impedance, space
charge, transverse dynamics, feedback, aperture, or loss model. Keegan's
committed impedance driver appears to apply two numerical representations of
the same resonator kick and contains unresolved 4x/6x/9x normalization evidence,
so collective effects are deliberately deferred to a separately labeled
sensitivity study.

The transition-gamma endpoint report retains two distinct +/-125 ns
definitions. The primary Werkema-aligned local-pulse proxy folds each particle
around its nearest h=28 bucket. The unwrapped diagnostic instead retains the
particle's initial-lineage reference center and therefore combines local tails
with migration to other buckets. Neither is an extinction prediction: the
model has no extraction sequence, Delivery Ring, extraction gate, M4 transport,
or production-target model, and no multi-value gamma scan or tail ensemble was
performed.

The simulation's physical weighting represents Keegan's `1.05e12` total
intensity. Werkema documents approximately `1e12` per h=28 bunch, or `8e12`
total, so the comparison uses endpoint shapes, fractions, and raw simulation
counts rather than interpreting stored physical-proton counts as Werkema
absolute predictions.

The full Keegan baseline was reused from an older HDF5 file that predates the
current executable-source hash and exact invocation fields. Its input SHA-256,
full row count, machine/RF metadata, endpoint, profile grid, and coordinate
conventions are checked against the new Werkema run, but byte-identical
executable source cannot be established from that legacy file alone.

The curated, small outputs that support the checked results are in
[`reference_results/`](reference_results/). Large HDF5 run files remain
regenerable local artifacts under the ignored `outputs/` directory. The
portable transition-gamma [report](reference_results/publication_3kv/gamma_transition_full/final_dt_comparison.json)
and [plot](reference_results/publication_3kv/gamma_transition_full/final_dt_comparison.png)
record the full 3 kV endpoint comparison.
