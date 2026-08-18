# ROOT analysis and input-generation code

The macros here generate G4Beamline-format particle distributions, convert
delivery-ring/MARS text input, inspect survival and profiles, and calculate
transmission. In particular:

- `EmitGen.C` generates Gaussian or uniform-emittance input distributions.
- `DRGen.C` converts Vladimir Nagaslaev's downstream-C-magnet distribution,
  changes units, optionally rematches Y phase space, and can smear events.
- `makedtscripts.py` expands the `xxxx` angle-scan templates.
- The remaining ROOT macros implement backtracking, cuts, profiles, survival,
  track generation, and transmission calculations.

Several constants, paths, tree names, and optics values are hardcoded. Read
Eric's [original study notes](../../historical_docs/originals/AAA_g4beamline_study.html)
before attempting a rerun.
