# Keegan-code 5 kV reference results

These are the portable outputs from the default RF-only convergence study. The
h=28 system starts at 5 kV, matching Keegan's committed BLonD driver, while the
sequential 5 ms h=588 turn-off plus 85 ms h=28 ramp follows the reconstructed
timing described in [`../../provenance.md`](../../provenance.md).

The study contains 16 complete 8,083-turn simulations:

- matched historical/generated pairs at 33,600, 100,800, 336,000, and
  1,075,200 particles for seed 202208;
- four additional matched pairs at 100,800 particles for seeds 202209--202212;
- 810 longitudinal profiles per run, recorded every 10 turns.

All 16 runs conserve the represented 1,075,200-macroparticle population, retain
every tracked particle with finite final coordinates, and share tracking-source
SHA-256 `c3ad1372c76659bf6dc637c3150736883b378e65d0bff5991b1bef365b107d5f`.
The five-seed 100,800-particle screen passes all declared coarse criteria.

## Contents

- [`convergence_summary.png`](convergence_summary.png) and
  [`convergence_summary.json`](convergence_summary.json): the particle-count
  ladder, reduced-to-full checks, five-seed ratios, run integrity, provenance,
  and timing.
- [`full_historical/`](full_historical/): standard voltage, phase-space, and
  waterfall plots from all 1,075,200 historical particles.
- [`full_comparison/`](full_comparison/): time-resolved and final-arrival-time
  comparison plots/reports for the full historical and generated runs.
- [`screen_n100800/rf_screen_summary.json`](screen_n100800/rf_screen_summary.json):
  five generated seeds against four historical-resampling controls.
- [`rf_start_voltage_sensitivity.json`](rf_start_voltage_sensitivity.json): a
  same-particle 3 kV versus 5 kV endpoint comparison.

The 16 HDF5 files total about 99 MB and are intentionally omitted. Their hashes
remain in [`manifest.json`](manifest.json), while each HDF5 internally records
its exact coordinates, sampling, dependencies, resolved configuration, Git
state, executable-source hash, and invocation. Recreate and analyze the matrix
from the project root with:

```bash
./run_convergence_5kv.sh
./analyze_convergence_5kv.sh
```

Verify the listed portable result artifacts with:

```bash
sha256sum --check reference_results/keegan_5kv/SHA256SUMS
```

The current source-tree hash differs from the tracking hash only because the
screen's hard-coded old `33,600-particle` label was made data-driven after the
runs. That analysis-only correction did not alter the stored trajectories.
