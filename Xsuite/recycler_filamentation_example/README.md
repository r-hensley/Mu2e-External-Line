# Recycler filamentation example

This project contains two implementations of the same longitudinal Recycler
teaching simulation:

- [`manual_version/`](manual_version/) evaluates an explicit NumPy
  drift--kick map derived from the USPAS longitudinal-dynamics demonstration.
- [`xsuite_version/`](xsuite_version/) represents the same one-turn map with
  Xsuite particles, a nonlinear `LineSegmentMap`, and two `Cavity` elements.

The purpose of keeping both versions is to test the implementation change in
isolation. They use the same deterministic initial ensemble, RF voltage
program, snapshot schedule, and longitudinal diagnostics. This makes the
manual trajectory a direct numerical reference for the Xsuite trajectory.

This is a classroom-scale single-particle model. It is not an operational
Recycler or extinction simulation.

## Simulated cycle

One h=28 ring cell initially contains 21 h=588 bunchlets. The model applies:

1. one turn at 80 kV on the h=588 system;
2. a 5 ms linear h=588 turn-off;
3. an 85 ms iso-adiabatic h=28 ramp from 3 kV to 80 kV;
4. continued tracking on the h=28 flattop through eight observation times,
   separated by 48.12 ms.

The eight observation times are called extraction times because they follow
the timing used by the source demonstration. No extraction kicker, aperture,
particle removal, transverse motion, or beam-line transport is modeled.

## Running the simulations

Run commands from this directory. Generate the manual reference first, then
run Xsuite against it:

```bash
python3 manual_version/run.py
python3 xsuite_version/run.py
python3 -m pytest -q
```

The manual dashboard is the largest generated product. It can be skipped for
a faster noninteractive run:

```bash
python3 manual_version/run.py --no-dashboard
```

Use `python3 manual_version/run.py --help` or
`python3 xsuite_version/run.py --help` for the available output and sampling
options.

The manual code requires NumPy, SciPy, pandas, Matplotlib, and Plotly. The
Xsuite code additionally requires Xsuite/Xtrack/Xpart, and the test suite uses
pytest. The current verified run used NumPy 2.2.6, Xsuite 0.50.1, Xtrack
0.103.0, and Xpart 0.23.10.

## Implementation differences

Both implementations use the same drift-then-kick stroboscopic convention:
longitudinal slip is applied first, followed by the h=588 and h=28 RF kicks.

The manual version stores each particle as ring phase `theta` and energy
offset `Delta E`, evaluates the map explicitly with NumPy, and wraps phase into
one h=28 cell.

The Xsuite version stores particles in `zeta` and `ptau`. It tracks unwrapped
`zeta` through an Xsuite line and converts to the classroom coordinates only
when saving diagnostics:

```text
theta   = -2*pi*zeta/C
Delta E = ptau*p0c
```

The supplied Recycler circumference, reference energy, transition gamma,
harmonics, and RF program remain model inputs. Xsuite constructs the reference
particle and derives `p0c`, beta, gamma, revolution period, slip factor, RF
frequencies, and flattop synchrotron tune consistently from its line. Because
the line is a lumped `LineSegmentMap`, these are not independently calculated
from a magnet-by-magnet Recycler lattice.

## Verified comparison

The nominal comparison tracks 5,376 particles for 38,336 turns. All particles
survive and all recorded coordinates remain finite.

| Quantity | Manual | Xsuite | Relative difference |
|---|---:|---:|---:|
| Flattop sigma time | 35.904460 ns | 35.903945 ns | -0.001435% |
| Flattop sigma energy | 10.154413 MeV | 10.154583 MeV | +0.001678% |
| Eighth-observation sigma time | 34.402049 ns | 34.403783 ns | +0.005040% |
| Eighth-observation sigma energy | 10.666993 MeV | 10.666619 MeV | -0.003504% |

At the 90 ms flattop, the particlewise RMS differences are 0.0089 ns and
3.46 keV. At the eighth observation they are 0.0312 ns and 16.2 keV. The
independent one-turn convention check agrees to approximately
`6.9e-16 rad` and `0 eV`.

Particlewise differences grow slowly during filamentation because tiny
floating-point differences shift fine phase-space filaments. The distribution
widths are therefore the primary long-time acceptance quantities.

## Outputs

Each implementation writes into its own `outputs/` directory. The manual
version produces the NumPy reference trajectory, figures, diagnostics, and an
optional self-contained dashboard. The Xsuite version produces its trajectory,
figures, a JSON comparison with the manual reference, and a run summary.

The summaries store project-relative source paths so moving the whole example
does not immediately invalidate their provenance fields.

For implementation-specific details, see
[`manual_version/README.md`](manual_version/README.md) and
[`xsuite_version/README.md`](xsuite_version/README.md).

## Model boundary

The model omits collective effects, impedance, longitudinal space charge,
feedback, RF noise, transverse optics, apertures, losses, and extraction
hardware. The nonzero transverse tunes supplied to the Xsuite
`LineSegmentMap` only make a complete 6D Twiss calculation possible; tracked
particles have zero transverse coordinates.
