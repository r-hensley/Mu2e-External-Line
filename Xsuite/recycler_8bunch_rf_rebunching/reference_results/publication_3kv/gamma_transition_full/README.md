# Full 3 kV transition-gamma comparison

This portable checkpoint compares the same 1,075,200 historical input
particles tracked for all 8,083 turns with two Recycler transition-gamma
values. The reference uses Keegan's BLonD value
`gamma_transition=20.257643837730637`; the candidate uses Werkema's
`gamma_transition=21.6`. Both use the default 3 kV h=28 starting voltage.

The controlled comparison permits only transition gamma and its two derived
machine quantities, momentum compaction and slip factor, to differ. The input
file, selected rows, RF program, endpoint, profile grid, particle weights, and
all other tracked machine quantities must match.

## Portable files

- [`final_dt_comparison.png`](final_dt_comparison.png) overlays the endpoints,
  shows their density residual, separates local-pulse tails from lineage
  migration, and gives the folded local-window fraction by initial bunch.
- [`final_dt_comparison.json`](final_dt_comparison.json) records the strict
  validation, machine-parameter changes, distribution distances, endpoint
  summaries, bucket migration, and both 125 ns diagnostics.

The exact final absolute-time Wasserstein-1 distance is 1.1791 ns and the
2 ns-bin Jensen--Shannon distance is 0.08262. The final-voltage small-amplitude
synchrotron-period approximation changes from 18.7242 ms to 18.4167 ms.

The Werkema-aligned Recycler-endpoint proxy first folds every particle about
its nearest h=28 bucket and then applies the strict rule
`abs(dt_local) > 125 ns`. It finds 36 reference and 39 candidate particles
outside, or 0.003348% and 0.003627%. Because the same input rows are paired,
the report also shows 39 inside-to-outside and 36 outside-to-inside changes,
for a net increase of three particles.

The separate unwrapped lineage diagnostic uses the reference endpoint's h=28
center for each initial lineage without folding. It finds 156,202 reference
and 158,302 candidate particles outside 125 ns. That much larger number mixes
local tails with whole-bucket migration. Neither diagnostic is a Delivery Ring,
extraction-gate, or proton-target extinction prediction; this model remains
Recycler RF-only, without impedance or space charge.

## Omitted HDF5 intermediates

| Role | Local path | Bytes | SHA-256 |
|---|---|---:|---|
| Keegan-gamma reference | `outputs/final_dt_full_historical/rf_only_results.h5` | 16,108,475 | `9601d0643d172ce718a4bc9b27b4ed83d454d0700d8bfd1d0e4ed915e2295052` |
| Werkema-gamma candidate | `outputs/gamma_transition_3kv/werkema_gt21p6_historical_full/rf_only_results.h5` | 16,106,793 | `b97750b61ef57c6d40fc2b303730d32a7764c3f6c460ec8a6bb77890bf3ed5aa` |

The reference HDF5 is a legacy run made before the strengthened invocation,
source-tree, and Git provenance fields were added. Its exact HDF5 hash and
resolved physics metadata are retained, but the reference command below is a
reconstruction rather than a stored shell transcript. The candidate records
the exact invocation, project source SHA-256
`f3703067ecef427d848e133b72c67923bccd1d8e5ed09cd2f782bbd4b44fae2d`,
and Git revision `741d0e6e522525048b6e535f99a2df703fe53ce0` with a dirty, untracked project
directory.

No HDF5 file is copied into `reference_results/`. Recreate the intermediates
and comparison from the project root with:

```bash
python3 main.py run \
  --input 1075200_init_dist.txt \
  --max-particles 0 \
  --sampling-strata 168 \
  --seed 202208 \
  --turns 8083 \
  --record-every 8083 \
  --rf25-start-kv 3 \
  --gamma-transition 20.257643837730637 \
  --no-figures \
  --output-dir outputs/final_dt_full_historical

python3 main.py run \
  --input 1075200_init_dist.txt \
  --max-particles 0 \
  --sampling-strata 168 \
  --seed 202208 \
  --turns 8083 \
  --record-every 8083 \
  --rf25-start-kv 3 \
  --gamma-transition 21.6 \
  --no-figures \
  --output-dir outputs/gamma_transition_3kv/werkema_gt21p6_historical_full

python3 main.py compare-final-time \
  outputs/final_dt_full_historical/rf_only_results.h5 \
  outputs/gamma_transition_3kv/werkema_gt21p6_historical_full/rf_only_results.h5 \
  --allow-gamma-transition-difference \
  --absolute-bins 2400 \
  --absolute-min-ns -1000 \
  --absolute-max-ns 3800 \
  --bucket-relative-bins 400 \
  --output-dir outputs/gamma_transition_3kv/full_comparison
```

The input file must have SHA-256
`3b591b8381cb7cf57a88d7db08d535f2de182720d9fcb3c335800d7b3283e7e1`.
The root [`manifest.json`](../../manifest.json) records the same commands,
omitted-file hashes, and portable-artifact hashes.
