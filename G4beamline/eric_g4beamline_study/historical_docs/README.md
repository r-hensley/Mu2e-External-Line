# Documentation guide

The files under [`originals/`](originals/) are preserved historical documents,
not rewritten pages. This guide provides a shorter route through them.

- [G4Beamline study notes](originals/AAA_g4beamline_study.html): the main
  workflow—emittance generation, conversion of Vladimir Nagaslaev's particle
  data, Y-plane phase-space rematching, standard run parameters, and scan sets.
- [General Mu2e G4Beamline notes](originals/AAAREADME_g4beamline.html): geometry
  conventions, translated MAD-X fragments, bend placement, field-versus-iron
  lengths, and step-size cautions.
- [`mad2g4bl.py` notes](originals/AAAREADME_mad2g4bl.html): input table
  requirements, coordinate transformations, element handling, and generated
  segment files. The converter itself is in [`../conversion/`](../conversion/README.md).
- [Grid notes](originals/AAAREADME_g4beamline_grid.html): the historical Mu2e
  grid procedure; paths and commands should be treated as obsolete unless
  independently checked.
- [DocDB 4054 v4](originals/prebys_docdb-4054-v4%20Beam%20Transmission%20and%20Extinction%20in%20the%20Mu2e%20Beam%20Line.pdf): Eric's May 2014 transmission
  and extinction study, covering the then-current downstream design.
- [Saved server-directory index](originals/site_directory_index.html): the
  large-file listing captured during recovery.

The later study notes and DocDB paper are related but not the same generation
of work. DocDB 4054 describes a 2014 downstream-only study and a Perl-based MAD
conversion performed by Jean-Francois Ostiguy. The later notes document Eric's
Python converter and Eliana Gianfelice-Wendt's MAD-X optics.

Eric's 2026 clarification and permission to preserve this material are recorded
in [`EMAIL_PROVENANCE.md`](EMAIL_PROVENANCE.md).
