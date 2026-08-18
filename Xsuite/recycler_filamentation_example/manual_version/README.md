# Manual NumPy implementation

This directory contains the notebook-independent reference implementation of
the Recycler longitudinal filamentation demonstration. It evaluates the
accelerator equations explicitly with NumPy and supplies the deterministic
reference trajectory used to validate the sibling Xsuite implementation.

## Files

- `mu2e_longitudinal.py` contains the parameters, initial-distribution
  generator, one-turn map, diagnostics, validation helpers, and plotting
  functions ported from the USPAS live demonstration.
- `run.py` executes the nominal simulation, performs its numerical acceptance
  checks, writes the trajectory, and generates the figures.
- `__main__.py` provides the package-style entry point.

The longitudinal helper is intentionally preserved as the manual reference.
Changes to it alter both the reference physics and the initial ensemble used
by the Xsuite comparison.

## Coordinates and map

Particles are represented by:

- `theta_rad`: arrival-time phase around the Recycler ring, wrapped into one
  h=28 cell;
- `delta_energy_eV`: total-energy offset from the reference particle.

One turn first advances phase using the nonlinear momentum-dependent slip and
then applies the sinusoidal h=588 and h=28 energy kicks. The RF voltages are
evaluated at each turn from the one-turn porch, linear h=588 turn-off, and
iso-adiabatic h=28 ramp.

The initial state consists of 21 independently labeled h=588 bunchlets inside
one h=28 cell. Each bunchlet is sampled deterministically from the matched
parabolic classroom distribution, using 256 particles per source by default.

## Running

From the project root:

```bash
python3 manual_version/run.py
```

To omit the self-contained Plotly dashboard:

```bash
python3 manual_version/run.py --no-dashboard
```

Useful options include `--particles-per-source`,
`--diagnostic-stride-turns`, `--display-per-source`, and `--output-dir`.

## Validation and outputs

The runner checks the initial contour area, the numerical one-turn Jacobian,
the tracked small-amplitude synchrotron period, finite coordinates, and phase
wrapping. Its primary numerical product is `outputs/nominal_result.npz`, which
the Xsuite runner reads by default.

The saved trajectory contains selected phase-space snapshots rather than every
particle on every turn. Scalar diagnostics are recorded more frequently and
include the time and energy widths, covariance emittance, central 95% time
width, out-of-window fraction, and centroid.

This code is the transparent teaching calculation, not an Xsuite model and not
an operational Recycler prediction.
