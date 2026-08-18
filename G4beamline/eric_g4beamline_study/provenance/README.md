# Import provenance

The import used SHA-256 content identity and `copy2` timestamp preservation.

- `import_manifest.tsv` records each unique copied historical file's curated
  destination, byte count, source modification time in nanoseconds, mode,
  SHA-256 digest, and canonical source path.
- `source_path_map.tsv` records all 308 examined source-path occurrences,
  including exact duplicates that were not copied again.
- `deduplication.tsv` maps each of the 23 omitted duplicate occurrences to its
  retained canonical file.
- `DOWNLOAD_MANIFEST_2026-08-15.md` is the preserved local record of the
  original website download and archive extraction.

Canonical precedence was: later export, standalone server-directory files,
January 2017 archive, selected website documentation/tool, local download
record, and the supplied DocDB PDF. The recovery sources remain untouched in
their original workspace locations.
