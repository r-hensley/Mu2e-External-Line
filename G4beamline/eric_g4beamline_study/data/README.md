# Particle input data

[`large_inputs/`](large_inputs/README.md) contains the recovered text and ROOT
particle distributions used by the study. They are retained locally but ignored
by Git because they total about 172 MiB. Download links, sizes, and verification
instructions are in that directory's README.

Restore every ignored large input and reference output with:

```bash
./download_large_files.sh
```

The script downloads and verifies both remote `.tgz` archives when archive
members are missing, selectively extracts only the two required `.dat` inputs
and two large reference PDFs, and directly downloads the five ROOT inputs that
are not present in either archive. Existing matching files are retained;
existing mismatched files cause a safe failure instead of being overwritten.
