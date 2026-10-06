# September 2026 Delivery Ring and M4 audit

Audit completed on **October 5, 2026**. Provenance statements below are
paraphrased from author correspondence; no message bodies are reproduced.

## Findings

The contents support a September 2026 updated calculation, with actual changes
from the preserved August 2018 and October 2024 models. The extraction file is
a self-contained Q303-to-diagnostic-absorber model, not a complete M4 model to
the production target. Its reproducibility is established; adoption of the
proposed corrected extraction settings is not established.

## Author-provided provenance

On September 29, 2026, George Deinlein identified the ring model with the
existing quadrupole configuration. He described the extraction settings as
planned corrections for a future configuration, rather than settings already
confirmed to be in operation.

On October 5, he clarified that the extraction model begins near Q303 inside
the ring and follows the extracted route to the diagnostic absorber. He
described its shared M4 coverage as approximately 170 m. His extension to the
production target was unfinished, so he suggested obtaining the remaining
approximately 70 m from an Eliana model. No particular tail file or approved
joining plane was specified.

These are paraphrased provenance statements, not an independently retrieved
mailbox record. Route coverage was checked against the files; operational
adoption and the appropriate production-tail version remain unresolved.

## Supplied files and dates

The six supplied files are kept together in
`MADX/m4_lattice_archive/2026_09_29_Deinlein_DR_M4/`. Each model has
a source `.madx`, a Twiss `.tfs`, and a four-page plot PDF:

- `mu2e-dr-*-v2026.09.29*`: Delivery Ring, 505.2822 m, 1,997 Twiss rows.
- `mu2e-dabs-*-v2026.09.29*`: Q303 to diagnostic absorber, 226.137 m,
  491 Twiss rows.

The supplied text copies use Linux (LF) line endings at Ryan's request on
October 6, 2026. Only CRLF terminators were replaced; lattice statements, table
values, and PDF bytes are unchanged. `provenance/file_manifest.json` records
current sizes/hashes and the original identities as `original_bytes` and
`original_sha256`. Current checks read the LF copies directly. The archive-wide
conversion record is [line_ending_normalization.tsv](../../manifests/line_ending_normalization.tsv).

Both supplied TFS files report ORIGIN="5.09.01 Linux 64" and DATE="29/09/26".
Ring TIME="12.42.05"; extraction TIME="13.18.05". TFS timestamps give no zone.
The PDFs also carry September 29 creation metadata. Their encoded offsets do
not warrant assuming the TFS time zone; their shared minute/second values are
additional consistency evidence rather than an independent source date.
Source headers contain no separate explicit creation-date history. The
September 29 labels and output dates support an updated calculation at that
date, not a claim that every element or calibration was developed in the
preceding weeks.

## MAD-X validation

Runtime: existing MAD-X 5.08.00 (64 bit, Windows), release 2022.01.13.
Executable SHA-256:
539815c372bdf3495bb97f4c7e8c75b4e30391649572d75d990da4106d743b6a
Entry points: each received .madx, submitted on standard input.
Dependency set: the selected .madx and MAD-X executable only; no external
CALL/READ inputs. Both sources run as received with exit status 0 and exactly
"Number of warnings: 0" / "MAD-X finished normally". No fatal messages.

The supplied WRITE/PLOT commands are commented out, so running the originals
computes Twiss/survey tables internally without producing TFS/PDF files.
Separate temporary validation copies enable only the two existing WRITE
commands. Both export runs also exit 0 with zero warnings. The validation
copies retain the original bytes apart from deleting the two leading "!"
characters on those WRITE commands; no lattice setting is changed.

The following sizes were recorded when MAD-X ran, before LF normalization:

- Generated extraction Twiss: 491 rows, 210,219 bytes.
- Generated extraction survey: 491 rows, 114,636 bytes.
- Generated ring Twiss: 1,997 rows, 845,725 bytes.
- Generated ring survey: 1,997 rows, 464,002 bytes.

All supplied row names, keywords, and all 20 numeric fields match the generated
Twiss tables exactly at their printed precision for both models. This includes
positions, lengths, strengths, optics, dispersion, phase, and beam-size columns.
Header DQMIN/DQMIN_PHASE values differ between versions; tiny ring orbit-header
values also differ. Thus the comparison establishes row-level reproducibility,
not equivalence of every MAD-X header coupling diagnostic. The receipt records
all header differences explicitly. No survey table was supplied for comparison.

Experiments were isolated under /tmp/mu2e-lattice-20261005-r8eszn1s. Copies of
the resulting logs, generated tables, and validation drivers are preserved in
this audit directory. The first sandboxed attempts failed before MAD-X with
WSL "UtilBindVsockAnyPort:309: socket failed 1"; allowing Windows
interoperability outside the sandbox resolved the launch issue. This was not a
lattice failure or a missing dependency.

## Route and coordinate interpretation

```text
  Q303 -> extraction septa -> CMAG -> Q901 ... Q933 -> HDA1 -> QDA01/QDA02
                                                        -> mACE -> absorber
                                    production continuation -> Q934 ... Q943
                                                           -> TIEIN -> target
```

The extraction model has no Q934-through-Q943, TIEIN, or production-target
continuation. Its 226.137 m starts 1 mm upstream of Q303, not at the M4 origin
used by Eliana's model. Useful received-table coordinates, in metres:
```text
  CMAG entrance M0_DECMAG:       27.8736
  CMAG exit M1_DECMAG:           30.0436
  Q933 center Q_DQ933:          197.6835
  Q933 exit M1_DQ933:           197.9864
  HDA1 entrance M0_DHDA01:      200.7170
  mACE center W_DMACE:          223.5833
  Absorber entrance M0_M4DABS:  224.6120
  Absorber exit M1_M4DABS:      226.1360
  Sequence end:                226.1370
```

Relative to the CMAG entrance, Q933 center is 169.8099 m and HDA1 entrance
is 172.8434 m. This explains the email's approximately 170 m common trunk:
the file also contains the separate diagnostic branch beyond that trunk.
Its full length is not a claim that the production target is 226 m from Q303.

Combining with Eliana requires selecting an identified production-tail version,
joining at a common plane after Q933 and before the diagnostic bend, and
checking coordinates, frames, optics, dispersion, lengths, and magnet settings.
Eliana's production route uses an unbent HDA1 placeholder instead of the
87.3 mrad diagnostic bend. Appending a production tail after the absorber
would follow the wrong route. No merged model or rematch was attempted here.

## Comparison with earlier references

The new source specifies power-supply currents, polynomial current-to-strength
conversions, family/individual correction factors, physical magnet lengths,
and explicit rolls. The source calls several factors "fudge factors"; their
calibration authority is not independently established by this inspection.

All 33 Q901-Q933 integrated quadrupole strengths differ at printed precision
from both the preserved 2018 and 2024 tables; 18 differ by more than 1%.
The 2024 common quad settings themselves equal the 2018 settings. Differences
range from -16.3241% (Q902 magnitude) to +14.0532% (Q909 magnitude); selected
examples are Q901 -3.4843%, Q903 -14.6094%, and Q925 +13.3110%.
The comparison sums the old half-magnets and evaluates new K1L*cos(2*TILT),
since the new file expresses many defocusing quadrupoles using positive K1
plus a pi/2 roll. A naive sign comparison would be misleading.

New modeled quad lengths are about 34.24-34.30 mm shorter than the old effective
lengths. This is a modeling-convention difference, not evidence that hardware
became shorter. For example Q901 is 0.4229 m versus 0.4572 m previously.

After aligning at CMAG entrance, 32 of the 33 common quadrupole centers agree
with both references within 0.150 mm. Q909 is the exception: its center is
0.30473061 m farther downstream, about one foot. That is a confirmed model
placement difference; whether it is a deliberate geometry correction remains
unresolved. The CSV records every comparison. Old IQnnn markers establish
magnet centers, avoiding the 2018 CENTRE versus 2024 exit table conventions.

The H910/H911/H917/H918 and H912/H916 bend magnitudes equal the 2024 study
values to roughly 1e-8 rad and differ from the 2018 values. New positive bends
with pi roll encode the same horizontal orientation as the old negative bends.
This confirms use of those bend values in this study, not their operational
adoption. V901 is also changed (0.11171438 versus 0.11341001 rad), and the
source uses its own upstream extraction optics rather than Eliana's original
CMAG-entrance initial conditions.

This source is not the missing historical generator of the 2022 TFS or the
undated mACE CSV. Those older artifacts have different settings/initial
conditions and diagnostic placements. The new file instead supplies a dated,
independently reproducible physical-element model of that route.

The 2022 diagnostic output, undated mACE optics, and 2025 production-tail
files cover related variants or partial studies. In particular, the 2025
Q937U-to-target tail inherits the 2024 production configuration. This new
package does not supersede all of those studies, and matching later bend
values does not establish its exact source-file ancestry.

## DR–M4 connection and simulation scope

Ring-to-extraction connection checked directly on 2026-10-05:
Both jobs use S=0 one millimetre upstream of Q303 entrance. Their Q303, Q302,
Q301, and Q202 table entries have exactly the same positions, lengths, rolls,
and integrated quadrupole strengths. The upstream quad supply settings also
agree. The ring's starting and ending beta/alpha/dispersion values are identical
at printed precision, consistent with its periodic Twiss calculation.
The extraction job instead supplies explicit initial values approximating that
ring solution: for example BETX=14.1246 versus the ring's 14.12461348 m,
BETY=5.1552 versus 5.15518930 m, and DX=-0.0235 versus -0.02352730 m.
Small vertical dispersion values are set to zero in the extraction seed.
These values are hard-coded, not imported from a ring output table.

The first 96 supplied table row names agree; identical names do not imply
identical transport throughout that overlap. DESSA/DESSB angles are zero on
the circulating-ring reference path and nonzero on the extraction path.
Q203/Q204 are combined-function bends in the extraction job, and Q205 is
represented by two bending halves, before the path proceeds through the
Lambertson and CMAG into M4. Thus the extraction job already includes the DR
extraction region; it is not simply an M4 file starting at CMAG. The two jobs
represent connected routes and consistent initial optics, but remain separate
calculations with no automatic parameter synchronization, Twiss transfer, or
particle handoff. Rerunning a modified ring does not update the extraction
file's starting constants. The connection receipt preserves the exact values.

The ring periodic starting optics closely match the rounded explicit seeds in
the extraction file, providing a documented initial-optics connection. This
does not supply a realistic extracted six-dimensional particle distribution.
The extraction and ring decks are optics/survey calculations; the extraction
deck has no active tracking or matching commands. AC929 bends are zero, trim
deviations are zero, and no AC waveform is configured. Named collimators lack
explicit aperture values, and the absorber is an INSTRUMENT element with a
length rather than a material-interaction model. NPART=1e12 does not add
collective effects to the shown Twiss calculation. These files alone do not
predict extraction efficiency, collimator losses, scattering, or extinction.

## Questions for future integration

- Which exact Eliana production-tail package and joining plane should be used?
- Is Q909's approximately 30.5 cm downstream displacement intentional?
- When will the proposed extraction quad corrections be implemented, and which
  calibration and length convention should tracking use?

## Preservation and evidence files

Binary files retain their original bytes. Text inputs, saved run evidence,
and comparison tables preserve all content except CRLF-to-LF terminators.
Original and normalized hashes are retained in the provenance and conversion
records. `receipt_from_initial_review.json` remains the historical pre-conversion
receipt; the current `receipt.json` describes the normalized saved files.
The initial source-validation audit used local files; subsequent documentation
and archive changes are recorded in repository history. The repository README
is unchanged. The local context final inventory and
public-facing provenance companion are updated for these new findings.

| File or directory | Contents |
|---|---|
| [receipt.json](receipt.json) | Supplied-file hashes, run receipts, row comparisons, header differences, and the DR–extraction interface |
| [quadrupole_comparison.csv](quadrupole_comparison.csv) | All 33 common magnets versus each reference |
| [inspect_files.py](inspect_files.py) | Repeats the saved comparisons without executing MAD-X |
| `run_*` | Actual run logs, generated validation tables, and export-only drivers |
| `reference_tables/` | Unchanged copies of the older 2018 and 2024 comparison tables |
| [receipt_from_initial_review.json](receipt_from_initial_review.json) | Historical initial-review snapshot; current checks use receipt.json |

## Recheck the saved evidence

From the parent collection folder, run:

```bash
sha256sum -c SHA256SUMS
python3 audit/inspect_files.py
```

The script verifies the supplied-file and comparison-table hashes against
`provenance/file_manifest.json`, then compares the saved validation outputs
with the supplied tables and older references. It reads the extracted files
directly and does not execute MAD-X. Rechecking writes this audit directory's
`receipt.json` and `quadrupole_comparison.csv`. The initial-review snapshot
is retained as a historical record and is not an input to these checks.
