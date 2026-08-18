# Xsuite implementation

This directory reimplements the manual Recycler filamentation example with
Xsuite while retaining the same initial ensemble, RF program, turn ordering,
diagnostics, and observation schedule. The manual result is used only as a
validation reference; Xsuite performs the particle tracking.

## Files

- `simulation.py` defines the supplied inputs, constructs the Xsuite line,
  derives line/reference quantities, converts coordinates, tracks particles,
  and calculates diagnostics.
- `plotting.py` generates the Xsuite phase-space, waterfall, RF-program, and
  comparison figures.
- `run.py` provides the command-line workflow, reads the manual reference,
  enforces acceptance tolerances, and saves results and metadata.
- `__main__.py` provides the package-style entry point.

## Xsuite line

The one-turn line contains, in order:

1. a nonlinear `LineSegmentMap` representing the full-ring longitudinal slip;
2. an h=588 `Cavity`;
3. an h=28 `Cavity`.

That slip-then-kick order deliberately matches the manual teaching map. Both
cavity voltages are updated before every turn from the shared time-dependent
RF program. The `LineSegmentMap` internal RF arrays are set to zero so the
explicit cavities are the only sources of energy kicks.

The map is supplied with small nonzero transverse tunes so Xsuite can perform
a complete 6D Twiss calculation. This is bookkeeping only: the example does
not model Recycler transverse optics and all tracked transverse coordinates
start at zero.

## Supplied and derived quantities

The line is supplied with the published/classroom circumference, reference
energy, transition gamma, h=588 and h=28 harmonics, RF voltages, and RF timing.
The momentum compaction passed to the lumped map is `1/gamma_transition**2`.

Xsuite then supplies the reference-particle momentum, beta and gamma, and its
Twiss calculation supplies the revolution period, slip factor, RF frequencies,
and flattop synchrotron tune. Thus the tracking and derived quantities are
self-consistent with the Xsuite line, although they are not derived from a
magnet-by-magnet Recycler lattice.

## Coordinate convention

The manual distribution is converted to Xsuite with:

```text
zeta = -C*theta/(2*pi)
ptau = Delta E/p0c
```

Xsuite tracks unwrapped `zeta`. Saved snapshots convert back with
`theta = -2*pi*zeta/C` and wrap only the reported phase into one h=28 cell.
The energy conversion `Delta E = ptau*p0c` uses Xsuite's exact `ptau`
definition rather than a small-momentum approximation.

## Running

Generate the manual reference first. Then, from the project root, run:

```bash
python3 xsuite_version/run.py
```

The default reference is
`manual_version/outputs/nominal_result.npz`. Override it with
`--reference-result` when intentionally comparing with another compatible
manual run. Other useful options are `--particles-per-source`,
`--diagnostic-stride-turns`, and `--output-dir`.

Run the focused test suite from the project root:

```bash
python3 -m pytest -q
```

## Validation

Before the full run, `validate_one_turn` independently calculates the expected
slip and both RF kicks for a probe particle. This checks the phase sign,
`zeta`/`ptau` conversion, element ordering, and energy-kick convention.

After tracking, the runner requires all particles to survive, all saved
coordinates to remain finite, and the manual and Xsuite distribution widths to
agree at the 90 ms flattop and eighth observation. Particlewise differences
are recorded for transparency but are not used as the long-time acceptance
criterion because filamentation amplifies tiny floating-point phase shifts.

The main products are `outputs/xsuite_result.npz`,
`outputs/comparison_with_numpy.json`, `outputs/run_summary.json`, and six PNG
figures.
