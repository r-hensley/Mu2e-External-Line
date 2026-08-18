# G4Beamline scan decks

These directories preserve Eric's generated angle-scan decks. Templates use
`xxxx` in a filename and in `deltaTheta`; `makedtscripts.py` generated concrete
values from 2.0 through 3.0 in 0.1 steps and then 4.0 through 10.0.

- `mu2e-1-1-STANDARD`: complete standard configuration.
- `mu2e-1-1-NOCLEAN`: cleaning-related element(s) removed.
- `mu2e-1-1-NOTAILCOL`: tail collimator removed.
- `mu2e-1-1-NOUSCOL`: upstream collimator removed.
- `mu2e-0-1-W`: downstream-only tungsten scan.
- `mu2e-0-1-SS`: downstream-only stainless-steel scan.

The exact meaning of run switches and input files is documented in Eric's
[original study notes](../../historical_docs/originals/AAA_g4beamline_study.html).
The small plots and logs produced by some scans are under
[`../../reference_results/scan_outputs/`](../../reference_results/scan_outputs/).
