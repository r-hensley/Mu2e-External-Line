# Mu2e External Beam Line Simulation Models

This repository contains accelerator lattice and particle-tracking models of the
Mu2e External Beam Line (M4 line).

The beam line is currently modeled using three simulation codes:

- MAD-X
- G4Beamline
- BDSIM

The purpose of this repository is to maintain validated versions of the lattice
models, provide version control as the models evolve, and make the models
available for use within the Mu2e collaboration.

## Repository Structure

### MAD-X

Contains the MAD-X lattice model and associated input and output files.

The original MAD-X model was created by **Dr. Eliana Gianfelice-Wendt in 2018**.

The MAD-X directory includes the lattice definitions, magnet parameters, survey
information, and Twiss output used to describe the M4 beam line.

Current validated version:

`madx-v2018.08`

### G4Beamline

Contains the G4Beamline model of the M4 beam line together with the input
particle distribution used for tracking studies.

The G4Beamline model was created by **Diktys Stratakis**.

Current validated version:

`g4bl-v0.3`

Main lattice file:

`G4_M4_Mu2e_03.g4bl`

### BDSIM

Contains the BDSIM model of the M4 beam line together with the corresponding
input particle distribution.

The BDSIM model was created by **Diktys Stratakis**.

Current validated version:

`bdsim-v1.0`

Main lattice file:

`M4_Mu2e_BDSIM_v1.gmad`

## Version Control

The MAD-X, G4Beamline, and BDSIM models are maintained independently. A new
version of one model does not require a new version of the other models.

Validated versions are identified using Git tags. The files on the `main`
branch represent the current validated versions of the simulation models.

Current validated tags are:

- `madx-v2018.08` — MAD-X reference model from August 2018
- `g4bl-v0.3` — G4Beamline version 0.3
- `bdsim-v1.0` — BDSIM version 1.0

Previous validated versions are preserved through Git tags and can be retrieved
from the repository history.

When a new model version is validated, the corresponding files and version
information in this README should be updated in the same Git commit. A new Git
tag should then be created for that validated version.

Development versions and archived local working files are not included in the
repository.

## Model Authors

- **MAD-X:** Dr. Eliana Gianfelice-Wendt
- **G4Beamline:** Diktys Stratakis
- **BDSIM:** Diktys Stratakis

## Notes

These models are intended for accelerator and beam-transport studies of the
Mu2e External Beam Line.

The models may continue to evolve as the lattice, magnet parameters, apertures,
and simulation requirements are updated.

Contributions and updates to the models should be documented through Git commits
and version tags to preserve the development history.
