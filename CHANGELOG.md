\# Mu2e External Beam Line Model Change Log



This file summarizes significant changes between validated versions of the

Mu2e External Beam Line simulation models.



Only validated versions are documented here. Development versions are maintained

locally and are not included in the repository.



\---



\## G4Beamline



\### v0.3 — 2026-08-10



Initial validated G4Beamline model included in this repository.



\- Added the G4Beamline lattice model.

\- Added the input particle distribution used for tracking studies.

\- Established this version as the initial reference model for future development.



\*\*Main lattice file:\*\* `G4\_M4\_Mu2e\_03.g4bl`



\*\*Git tag:\*\* `g4bl-v0.3`



\---



\## BDSIM



\### v1.0 — 2026-08-10



Initial validated BDSIM model included in this repository.



\- Added the BDSIM lattice model.

\- Added the corresponding input particle distribution.

\- Established this version as the initial reference model for future development.



\*\*Main lattice file:\*\* `M4\_Mu2e\_BDSIM\_v1.gmad`



\*\*Git tag:\*\* `bdsim-v1.0`



\---



\## MAD-X



\### August 2018 Reference Model



Original MAD-X model created by \*\*Dr. Eliana Gianfelice-Wendt in 2018\*\*.



\- Added the original MAD-X lattice and magnet definition files.

\- Added the associated magnet parameter files.

\- Added the reference survey and Twiss output files.

\- Established the August 2018 model as the current MAD-X reference model.



\*\*Git tag:\*\* `madx-v2018.08`



\---



\## Future Updates



For each new validated model version, a new entry should be added describing

the significant changes relative to the previous validated version.



Each validated update should include:



\- The new version number and validation date.

\- A concise description of the changes.

\- The main files that were modified or replaced.

\- Any significant changes to the lattice, magnet parameters, apertures,

&#x20; geometry, beam distributions, or simulation configuration.

\- The corresponding Git tag.



Detailed file-level changes remain available through the Git commit history.

