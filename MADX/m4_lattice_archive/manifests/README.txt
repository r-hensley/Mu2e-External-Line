ARCHIVE INVENTORIES

m4_lattice_archive.sha256
m4_lattice_archive.tsv
    Current inventories of every regular file outside manifests/, including
    the dated lattice folders, September 2026 collection, saved audit evidence,
    and guides. The manifests/ directory is excluded to avoid self-reference.
    SHA-256 paths are relative to the archive root. TSV columns record path,
    size in bytes, modification time in UTC with nanosecond precision, file
    mode, provenance, and SHA-256. File times are not lattice configuration dates.

From the archive root, verify the listed contents:

    sha256sum -c manifests/m4_lattice_archive.sha256

2026_09_29_Deinlein_DR_M4/SHA256SUMS
    Portable checksum list for this collection, excluding the list itself.
    Its provenance/file_manifest.json records the fixed supplied-file and
    comparison-table identities. The archive-wide list includes this subset.

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

These four historical records are not refreshed with the current inventories.
Their SHA paths are relative to the former source roots; the Eliana TSV provides
its merged-path mapping. All checksum lists use sha256sum's "hash  path" format.

LINE ENDINGS

Archive text files use Linux (LF) line endings. line_ending_normalization.tsv
records the original and normalized sizes/hashes of the CRLF-to-LF conversions.
Current inventories describe the normalized copies. Historical pre-merge records
still describe the original deliveries and must not be treated as current hashes
for converted text files. Binary files are unchanged.
