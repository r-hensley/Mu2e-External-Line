MANIFEST CONTENTS

dean_files_premerge.sha256
dean_files_premerge.tsv
    Snapshot of the former dean_files/ tree immediately before the merge.

eliana_files_premerge.sha256
eliana_files_premerge.tsv
    Snapshot of the former eliana_files/ tree immediately before the merge.
    The TSV maps each original path to its retained merged-archive path.

m4_lattice_archive.sha256
m4_lattice_archive.tsv
    Snapshot of the final merged archive, excluding the manifests/ directory
    itself so the SHA-256 list remains non-recursive and checkable.

The .sha256 files use the standard "hash  path" format. Pre-merge SHA paths
are relative to their former source roots; the Eliana TSV supplies the current
merged-path mapping for its deduplicated files. The TSV files add sizes, exact
original modification timestamps, file modes, provenance, disposition, and
SHA-256 values.
