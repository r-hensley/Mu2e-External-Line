# September 2026 Deinlein Delivery Ring and M4 models

This folder collects the `v2026.09.29` Delivery Ring and extraction models
provided by George Deinlein. The extraction model runs from Q303 to the
diagnostic absorber; it does not include a production-target extension.

The text files use Linux (LF) line endings; PDF bytes are unchanged. Current
and original sizes/hashes are recorded in
[provenance/file_manifest.json](provenance/file_manifest.json). The full
[line-ending record](../manifests/line_ending_normalization.tsv) covers the archive.

## Main points

- **Settings and date:** The author described the Delivery Ring (DR) model as
  using the existing quadrupole settings and the extraction model as using
  proposed future corrections. The September 29, 2026 labels and output dates
  establish a dated calculation; operational adoption of the extraction
  corrections remains unconfirmed.
- **How DR and M4 connect:** Both jobs start near Q303 inside the ring and share
  the upstream quadrupole settings. The extraction file contains explicit,
  rounded starting optics consistent with the ring's periodic solution, then
  follows the extraction septa and CMAG into M4. Each file runs independently.
  For example, changing and rerunning the DR file will not automatically update
  the extraction file's starting values.
- **What is simulated:** These are static optics and survey calculations,
  describing the reference trajectory, beam envelope, and dispersion. They do
  not track a bunch, simulate time-dependent resonant extraction, or transfer
  particles between files. There is therefore no simulated instantaneous bunch
  transfer down M4. Extraction efficiency, losses, scattering, and extinction
  require additional simulation inputs and calculations.
- **Relationship to earlier lattices:** All 33 common Q901–Q933 integrated
  quadrupole strengths differ from the preserved 2018 and 2024 references after
  accounting for roll and sign conventions. Several bend values agree with the
  2024 study, and Q909 is about 30.5 cm farther downstream in this model. These
  files contain revised settings and geometry, but their exact ancestry and
  whether they supersede intermediate studies remain unresolved.
- **Production-target coverage:** The shared M4 trunk reaches Q933 about 170 m
  from the CMAG entrance. The extraction file's full 226.137 m also includes
  the upstream ring extraction region and the diagnostic branch. The author
  suggested an Eliana model for the remaining approximately 70 m to the target,
  without specifying a tail version. A combined production model needs a
  checked join after Q933 and before the diagnostic bend; the target tail
  cannot simply be appended after the absorber.

The saved audit establishes that both sources run under MAD-X 5.08.00 with
zero warnings and reproduce every printed numeric field in the supplied Twiss
rows. Some header diagnostics differ between versions. Full comparisons,
validation records, and integration questions are in [audit/AUDIT.md](audit/AUDIT.md).

## Supplied files

| File | Contents |
|---|---|
| [mu2e-dr-lattice-v2026.09.29.madx](mu2e-dr-lattice-v2026.09.29.madx) | Complete Delivery Ring optics/survey job; 505.2822 m |
| [mu2e-dr-twiss-v2026.09.29-madx.tfs](mu2e-dr-twiss-v2026.09.29-madx.tfs) | Supplied ring Twiss table; 1,997 rows |
| [mu2e-dr-plots-v2026.09.29-madx.pdf](mu2e-dr-plots-v2026.09.29-madx.pdf) | Four supplied ring survey/optics plots |
| [mu2e-dabs-lattice-v2026.09.29.madx](mu2e-dabs-lattice-v2026.09.29.madx) | Q303 through extraction and M4 to the diagnostic absorber; 226.137 m |
| [mu2e-dabs-twiss-v2026.09.29-madx.tfs](mu2e-dabs-twiss-v2026.09.29-madx.tfs) | Supplied extraction Twiss table; 491 rows |
| [mu2e-dabs-plots-v2026.09.29-madx.pdf](mu2e-dabs-plots-v2026.09.29-madx.pdf) | Four supplied extraction survey/optics plots |

## Documentation and evidence

| Location | Contents |
|---|---|
| [audit/AUDIT.md](audit/AUDIT.md) | Detailed provenance, DR–M4 connection, simulation scope, comparisons, validation, and unresolved questions |
| `audit/` | Receipts, comparison tables, run logs, generated outputs, and the audit helper |
| [provenance/file_manifest.json](provenance/file_manifest.json) | Supplied-file sizes, hashes, and provenance record |
| [SHA256SUMS](SHA256SUMS) | Checksums of the collected files |

The broader model history is in
[madx_file_provenance.md](../../../madx_file_provenance.md).
