# Xsuite Recycler RF rebunching

This folder contains an RF-only Xsuite model of Recycler bunch formation. It
tracks the programmed 53 MHz and 2.5 MHz RF systems through the 90 ms
rebunching cycle. Longitudinal space charge and machine impedance are not yet
included. The standard reconstruction is locked to the publication/ESME
3 kV h=28 start: h=588 turns off over the first 5 ms, then h=28 ramps from
3 kV to 80 kV over the next 85 ms.

See [results.md](results.md) for the plot guide, publication cross-references,
and validated numerical results. [provenance.md](provenance.md) separates the
historical-code inputs, publication inputs, and reconstruction choices. Small,
versionable result checkpoints live in
[`reference_results/`](reference_results/README.md).

## Setup

Python 3.10 or newer and the dependencies in `pyproject.toml` are required:
NumPy, h5py, Matplotlib, and Xsuite. From this directory, use the top-level
`main.py` directly:

```bash
python3 main.py --help
```

No editable install or `PYTHONPATH` setting is needed. Run the tests with:

```bash
python3 -m pytest -q
```

## Quick start

With no input file, `run` generates 200 particles in each of 168 microbunches
(33,600 particles total) and tracks the complete RF cycle:

```bash
python3 main.py run
```

The default output directory is `outputs/generated/`. A short mechanics check
can be run with:

```bash
python3 main.py run \
  --turns 100 \
  --record-every 10 \
  --output-dir outputs/smoke_generated
```

To use the historical particle list instead, provide `--input`:

```bash
python3 main.py run \
  --input 1075200_init_dist.txt \
  --max-particles 20000 \
  --record-every 10 \
  --output-dir outputs/historical_20k
```

File mode samples 20,000 particles by default and otherwise writes to
`outputs/from_file/`. Use `--max-particles 0` to track every input row.

## Commands

| Command | Purpose |
|---|---|
| `run` | Generate particles by default, or load them from `--input`; track and save results. |
| `plot` | Regenerate the standard plots from an existing run HDF5 file. |
| `compare-inputs` | Compare generated particles statistically with a reference input file without tracking. |
| `compare-runs` | Compare turn profiles and final per-bunch properties from two compatible runs. |
| `compare-final-time` | Compare final arrival-time distributions and h=28 bucket occupancy. |
| `summarize-screen` | Aggregate several generated runs and historical-resampling controls. |

Use `python3 main.py COMMAND --help` for the current argument list.

## Code map

The implementation is kept in a small set of single-purpose modules:

| Module | Role |
|---|---|
| `config.py` | Recycler constants and the time-dependent RF-program settings. |
| `generator.py` | Statistical construction of the 168 initial microbunches. |
| `distribution.py` | Historical-file loading, sampling, weights, and conversion to Xsuite `zeta`/`ptau`. |
| `rf_program.py` | Voltage functions evaluated at the start of each turn. |
| `simulation.py` | Xsuite one-turn line, particle tracking, profiles, and HDF5 output. |
| `provenance.py` | Stable coordinate/source hashes plus Git, runtime, and dependency metadata. |
| `plotting.py` | Standard RF-program, final-phase-space, and waterfall plots. |
| `input_plotting.py` / `comparison_plotting.py` | Visual diagnostics for input and tracking comparisons. |
| `input_analysis.py` / `result_comparison.py` / `final_dt_comparison.py` | Statistical comparisons before and after tracking. |
| `screening.py` | Multi-seed dynamic screening against historical resampling controls. |
| `cli.py` | The command parsers and dispatch from `main.py`. |
| `make_convergence_summary.py` | Condense the completed archival 5 kV run matrix into one JSON report and four-panel plot. |

Package-level functions and class methods have docstrings explaining their
roles. Comments in the tracking and analysis paths call out Xsuite conventions
and physics choices that are not obvious from the NumPy operations alone.

## `run` inputs and parameters

### Common options

| Option | Default | Meaning |
|---|---:|---|
| `--input PATH` | omitted | Two-column particle file. When omitted, particles are generated in memory. |
| `--output-dir PATH` | `outputs/generated` or `outputs/from_file` | Destination for the HDF5 file and plots. |
| `--seed N` | `202208` | Generator RNG seed or deterministic file-sampling seed; valid range is `0 <= N < 2**64`. |
| `--record-every N` | `10` | Save one 2 ns longitudinal profile every N turns; the initial and final profiles are always included. |
| `--turns N` | `8083` when omitted | Override the complete 90 ms endpoint. |
| `--no-figures` | off | Write only the HDF5 result. Plots can be made later with `plot`. |
| `--rf25-start-kv V` | `3` | h=28 voltage at 5 ms. The standard default is locked to the publication/ESME value; use `5` only for the explicitly labeled legacy Keegan-code sensitivity. |
| `--gamma-transition VALUE` | `20.257643837730637` | Recycler transition gamma from Keegan's BLonD model. Use `21.6` only for the controlled Werkema comparison; this also derives a different `alpha_c` and slip factor. |

The stored full-run state is at 90.00137 ms. The final applied kick is sampled
at 89.99023 ms because 90 ms is not an integer number of Recycler revolutions.

### Generated-input options

These options apply only when `--input` is omitted.

| Option | `run` default | Meaning |
|---|---:|---|
| `--particles-per-microbunch` | `200` | Population in each of 168 microbunches; total count is 168 times this value. |
| `--time-origin-ns` | `-593.7974058` | Center time for occupied grid index zero. |
| `--spacing-ns` | `19.15473645` | Spacing between microbunch grid indices. |
| `--sigma-time-ns` | `5.0249773` | Gaussian time width within each microbunch. |
| `--mean-dpop` | `8.2046563e-8` | Mean relative momentum offset. |
| `--sigma-dpop` | `1.7593985e-4` | Gaussian `dp/p` width. |
| `--correlation` | `0` | Gaussian time-momentum correlation, constrained to `[-1, 1]`. |

The occupied grid indices are fixed to `0..83` and `126..209`: two trains of
84 microbunches separated by 42 empty grid slots. Consecutive groups of 21 are
assigned the eight structural bunch labels. The generator adds no explicit
ghost population, unequal population, per-microbunch jitter, or special tail
model.

### File-input options

These options require `--input`.

| Option | Default | Meaning |
|---|---:|---|
| `--max-particles N` | `20000` | Deterministic stratified sample size. Use `0` for every row. |
| `--sampling-strata N` | `8` | Ordered sampling strata; must be a positive multiple of eight. Use `168` to balance every historical microbunch in comparisons. |

Supplying a generated-only option with `--input`, or a file-only option without
`--input`, is an error. A nonzero sample size must be at least the number of
strata. `--seed` is intentionally valid in both modes.

## Runtime input files

| Input | Used by | Contents |
|---|---|---|
| `1075200_init_dist.txt` | `run --input`, `compare-inputs` | Historical 1,075,200-row distribution. Column 1 is longitudinal `z` in metres; column 2 is `dp/p`. |
| Generated distribution in memory | `run` without `--input`, generated side of `compare-inputs` | 168 equal-population Gaussian microbunches arranged as eight groups of 21. |
| `rf_only_results.h5` | `plot`, `compare-runs`, `compare-final-time`, `summarize-screen` | Run metadata, profiles, final particles, input provenance, and weights. |

The local historical input is byte-identical to
`../../../blond/1075200_init_dist.txt`; both have SHA-256
`3b591b8381cb7cf57a88d7db08d535f2de182720d9fcb3c335800d7b3283e7e1`.
The BLonD `table.dat` impedance table is not read by this RF-only model.

The input file is about 28.5 MB and is intentionally not ignored. It is below
GitHub's 100 MB per-file limit, so ordinary Git is the simplest option for this
single reference input; Git LFS is optional rather than required. Generated
HDF5 runs stay under the ignored `outputs/` directory. The selected plots and
JSON reports in `reference_results/` are small enough for ordinary Git.
Verify the input after a clone or transfer with:

```bash
sha256sum --check 1075200_init_dist.txt.sha256
```

### Historical coordinate conversion

The two input columns use the checked-in BLonD conversion:

```text
dt_blond    = input_z / 3e8
zeta_xsuite = -beta0 * c * dt_blond
dE          = input_dpop * p0c / 0.9944
ptau_xsuite = input_dpop / 0.9944
```

The minus sign converts BLonD arrival time to Xsuite's
`zeta = s - beta*c*t` convention.

## Standard machine and fixed RF model

The standard values are defined in `xsuite_recycler/config.py`. Transition
gamma has one controlled CLI override for source-comparison studies; the RF
durations, endpoint voltage, circumference, and momentum are not scan
parameters in the standard workflow.

| Quantity | Value |
|---|---:|
| Recycler circumference | 3319.41882346 m |
| Reference momentum | 8.83532 GeV/c |
| Transition gamma | 20.257643837730637 (Keegan default) |
| Represented intensity | `1.05e12` protons |
| Structural bunches | 8 |
| RF harmonics | 588 and 28 |
| Derived RF frequencies | 52.8081 MHz and 2.51467 MHz |
| h=588 voltage | 80 kV to 0 linearly over 0--5 ms |
| h=28 voltage | 0 before 5 ms; fixed inverse-square ramp from 3 kV to 80 kV over 5--90 ms (85 ms) |
| RF lag | 0 degrees for both systems |
| Profile window | -1000 to 3800 ns |
| Profile bin width | 2 ns (2,400 bins) |

Each turn applies the two RF kicks followed by a nonlinear
momentum-compaction slip map. There is no wake, impedance, space charge,
transverse motion, aperture, feedback, noise, or loss element.

For the two documented transition-gamma values, the derived quantities are:

| Source value | `gamma_t` | `alpha_c = 1/gamma_t^2` | Slip factor `eta` | Small-amplitude h=28 period at 80 kV |
|---|---:|---:|---:|---:|
| Keegan BLonD model (default) | 20.257643837730637 | 0.0024368126329700006 | -0.008714928548793705 | 18.72423519130096 ms |
| Werkema model (controlled override) | 21.6 | 0.002143347050754458 | -0.009008394131009248 | 18.41672119458596 ms |

Only `gamma_t` is independently changed in this comparison; `alpha_c`, `eta`,
and the synchrotron-period estimate follow from it. The period is a derived
small-amplitude approximation at the final 80 kV, not a tracked-particle scan.

## Outputs and plotting

A normal run writes:

| Output | Description |
|---|---|
| `rf_only_results.h5` | Machine/RF/source metadata, recorded profiles, final `zeta`, `ptau`, particle state, and initial bunch labels. |
| `rf_voltage_program.png` | 53 MHz turn-off and 2.5 MHz capture-voltage program. |
| `final_phase_space.png` | Final phase space of structural bunch 8, relative to its nearest h=28 center. |
| `rebunching_waterfall.png` | Turn-by-turn longitudinal density, logarithmic on the left and linear on the right. |

Reduced runs are weighted to represent the full 1,075,200-macroparticle
source. The HDF5 file records both macroparticle-equivalent and physical-proton
weights. It also records the source-file and selected-row identities (or the
realized generated-coordinate identity), RNG implementation, package versions,
command line, and project source/Git state needed to audit the run.

Regenerate all three plots without tracking again:

```bash
python3 main.py plot \
  outputs/historical_20k/rf_only_results.h5 \
  --output-dir outputs/historical_20k
```

## Comparisons

Compare the generated model with every historical input row. This command uses
6,400 generated particles per microbunch by default so both samples contain
1,075,200 particles:

```bash
python3 main.py compare-inputs \
  --input 1075200_init_dist.txt \
  --output-dir outputs/input_model_default
```

For a matched 33,600-particle tracking comparison:

```bash
python3 main.py run \
  --input 1075200_init_dist.txt \
  --max-particles 33600 \
  --sampling-strata 168 \
  --output-dir outputs/compare_historical_33k6

python3 main.py run \
  --particles-per-microbunch 200 \
  --output-dir outputs/compare_generated_33k6

python3 main.py compare-runs \
  outputs/compare_historical_33k6/rf_only_results.h5 \
  outputs/compare_generated_33k6/rf_only_results.h5 \
  --output-dir outputs/run_comparison_33k6
```

`compare-runs` requires matching machine, RF, endpoint, profile grid, and
population scaling by default. `--allow-partial-turn-overlap` is available only
for intentionally exploratory comparisons.

For the full final-arrival-time comparison:

```bash
python3 main.py run \
  --input 1075200_init_dist.txt \
  --max-particles 0 \
  --sampling-strata 168 \
  --record-every 8083 \
  --no-figures \
  --output-dir outputs/final_time_historical

python3 main.py run \
  --particles-per-microbunch 6400 \
  --record-every 8083 \
  --no-figures \
  --output-dir outputs/final_time_generated

python3 main.py compare-final-time \
  outputs/final_time_historical/rf_only_results.h5 \
  outputs/final_time_generated/rf_only_results.h5 \
  --output-dir outputs/final_time_comparison
```

For its unwrapped branch, the final-time comparison uses each historical
initial label's final h=28 bucket center for both samples, without independently
recentering the candidate. It separately reports a nearest-bucket folded local
window so local pulse shape is not conflated with whole-bucket migration.

For the controlled Keegan-versus-Werkema transition-gamma comparison, track the
same complete historical input and change only `--gamma-transition`:

```bash
python3 main.py run \
  --input 1075200_init_dist.txt \
  --max-particles 0 \
  --sampling-strata 168 \
  --record-every 8083 \
  --rf25-start-kv 3 \
  --gamma-transition 21.6 \
  --no-figures \
  --output-dir outputs/gamma_transition_3kv/werkema_gt21p6_historical_full

python3 main.py compare-final-time \
  outputs/final_dt_full_historical/rf_only_results.h5 \
  outputs/gamma_transition_3kv/werkema_gt21p6_historical_full/rf_only_results.h5 \
  --allow-gamma-transition-difference \
  --output-dir outputs/gamma_transition_3kv/full_comparison
```

The comparison flag is deliberately narrow: it requires the same full
file-backed input and permits only `gamma_transition` plus its derived
`momentum_compaction_factor` and `slip_factor` metadata to differ. The primary
125 ns pulse-shape proxy folds every survivor around its nearest h=28 bucket;
the separate unwrapped initial-lineage diagnostic includes both local tails and
migration away from that lineage's reference bucket. Neither quantity is a
Delivery Ring, extraction-gate, or proton-target extinction prediction.

### Reproduce the archival 5 kV convergence study

[`run_convergence_5kv.sh`](run_convergence_5kv.sh) runs a matched
historical/generated particle-count ladder at 33,600, 100,800, 336,000, and
1,075,200 particles. It also runs four additional historical-resampling and
generated seeds at 100,800 particles, giving five seeds total at that size.
Every run tracks all 8,083 turns, records profiles every 10 turns, and states
the non-default 5 kV RF setting explicitly. This matrix remains a useful
legacy-code sensitivity and convergence result; it does not define the current
3 kV standard configuration.

```bash
./run_convergence_5kv.sh
```

Each matched pair runs concurrently, but pairs are processed one at a time to
keep memory use modest. On a resumed run, an existing HDF5 file is skipped only
after its particle count, endpoint, RF setting, source identity, profiles, and
provenance are validated; an incomplete or incompatible file is refused rather
than overwritten. This is a long local calculation, and the full pair dominates
the runtime. Outputs and individual logs are written below
`outputs/convergence_5kv/`.

After all 16 HDF5 files exist, validate and analyze the whole matrix with:

```bash
./analyze_convergence_5kv.sh
```

That script validates all 16 runs, makes the four matched-size
evolution/endpoint comparisons, six reduced-to-full endpoint comparisons, the
five-seed screen with four historical resampling controls, standard plots from
the full historical run, and `convergence_summary.{json,png}`. It writes only
below `outputs/convergence_5kv/analysis/`.

The completed portable summary and selected full-statistics plots/reports are
in [`reference_results/keegan_5kv/`](reference_results/keegan_5kv/README.md).
The 16 local HDF5 intermediates total about 99 MB and remain ignored.

## Provenance and scope

The numerical machine specification and historical-coordinate conversion come
from the BLonD repository at commit
`3566338e5e1e201028f2ec3ec62a94dcf1fc0031`. The paper and dissertation supply
the h=588 turn-off stage and fixed sequential 5 ms plus 85 ms timing to 80 kV.
The locked standard h=28 starting voltage follows the publication/ESME
program's 3 kV value, while zero phase follows Keegan's committed driver. The
driver's 5 kV start remains available as the explicit legacy-code sensitivity
`--rf25-start-kv 5`; `--gamma-transition 21.6` similarly enables a controlled
Werkema machine-parameter comparison. The original production driver and exact
input-generation script are unavailable.

See the in-project [provenance note](provenance.md) for the durable source
audit boundary. Curated results, their SHA-256 hashes, and exact regeneration
recipes are in [`reference_results/`](reference_results/README.md).
This should be described as an RF-only, dynamically equivalent model rather
than a publication-identical reproduction.
