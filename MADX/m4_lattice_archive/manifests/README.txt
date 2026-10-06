CURRENT ARCHIVE INVENTORY

m4_lattice_archive.sha256
m4_lattice_archive.tsv
    Maintained inventories of every regular file outside manifests/, including
    all dated folders, supplied files, generated audit outputs, and guides.
    The manifests/ control directory is excluded to avoid self-reference.
    SHA-256 paths are relative to the archive root. TSV columns record path,
    size in bytes, modification time in UTC with nanosecond precision, file
    mode, provenance, and SHA-256. Modification times are file metadata, not
    lattice configuration dates. Existing provenance labels are retained.

From the archive root, check listed file contents:

    sha256sum -c manifests/m4_lattice_archive.sha256

For content checks plus complete path coverage and SHA-256/TSV agreement:

    python3 manifests/archive_inventory.py check

The second command fails for changed, missing, or unlisted files and for
inconsistent checksum/TSV records. Recorded timestamps and modes are metadata;
the check does not require checkout timestamps or permissions to be identical.

REFRESHING AFTER INTENDED ADDITIONS OR EDITS

    python3 manifests/archive_inventory.py refresh

Refresh validates supplied/reference-file identities from file_manifest.json,
preserves existing provenance labels, regenerates any existing folder-level
SHA256SUMS subsets, and writes both current archive inventories. It refuses
missing listed files and changes to listed contents unless explicitly accepted.
For an intentional archive-guide edit, for example:

    python3 manifests/archive_inventory.py refresh --accept-change README.md

Repeat --accept-change for each reviewed change, using archive-relative paths.
Changes to supplied/reference files still fail their provenance identity check.
Review removals explicitly rather than silently dropping missing listed files.
The helper uses only Python's standard library and does not execute MAD-X.

OPTIONAL PORTABLE SUBSETS

2026_09_29_Deinlein_DR_M4/SHA256SUMS
    Generated subset for verifying a standalone copy of this folder. The
    archive-wide inventory includes this subset file itself. Regenerate it
    with the archive refresh command, never as a separately maintained list.
    The folder's provenance/file_manifest.json records the fixed identities
    of its supplied files and comparison references, along with provenance.

Overview: 2026_09_29_Deinlein_DR_M4/README.md
Full audit: 2026_09_29_Deinlein_DR_M4/audit/AUDIT.md

HISTORICAL PRE-MERGE RECORDS

dean_files_premerge.sha256
dean_files_premerge.tsv
    Preserved snapshot of the former dean_files/ tree before the merge.

eliana_files_premerge.sha256
eliana_files_premerge.tsv
    Preserved snapshot of the former eliana_files/ tree before the merge.
    The TSV maps each original path to its retained merged-archive path.

These four records are historical provenance evidence and are not refreshed
by the current inventory helper. Their SHA paths are relative to the former
source roots; the Eliana TSV supplies the merged-path mapping.

All checksum lists use the standard "hash  path" format accepted by sha256sum.

LINE ENDINGS

Archive text files use Linux (LF) line endings. line_ending_normalization.tsv
records the original and normalized sizes/hashes of the CRLF-to-LF conversions.
Current inventories describe the normalized copies. Historical pre-merge records
still describe the original deliveries and must not be treated as current hashes
for converted text files. Binary files are unchanged.
