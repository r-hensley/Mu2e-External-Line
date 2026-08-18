#!/usr/bin/env bash
# Analyze the completed 5 kV convergence matrix from run_convergence_5kv.sh.
#
# This script deliberately uses the existing CLI instead of adding another
# Python workflow. It first validates every expected HDF5 result, so a partial
# or inconsistent simulation matrix cannot produce a misleading report tree.

set -euo pipefail

project_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$project_dir"

output_root=outputs/convergence_5kv
analysis_root="$output_root/analysis"

particle_counts=(33600 100800 336000 1075200)
extra_seeds=(202209 202210 202211 202212)
expected_results=()

# Seed 202208 is the four-point particle-count ladder for both input sources.
for particle_count in "${particle_counts[@]}"; do
    expected_results+=(
        "$output_root/historical_n${particle_count}_seed202208/rf_only_results.h5"
        "$output_root/generated_n${particle_count}_seed202208/rf_only_results.h5"
    )
done

# Four additional 100,800-particle pairs supply the seed-stability controls.
for seed in "${extra_seeds[@]}"; do
    expected_results+=(
        "$output_root/historical_n100800_seed${seed}/rf_only_results.h5"
        "$output_root/generated_n100800_seed${seed}/rf_only_results.h5"
    )
done

if [[ ${#expected_results[@]} -ne 16 ]]; then
    echo "Internal error: expected-results list has ${#expected_results[@]} entries, not 16" >&2
    exit 2
fi

missing_results=()
for result_path in "${expected_results[@]}"; do
    if [[ ! -s "$result_path" ]]; then
        missing_results+=("$result_path")
    fi
done

if [[ ${#missing_results[@]} -ne 0 ]]; then
    echo "Cannot analyze an incomplete matrix; missing or empty HDF5 results:" >&2
    printf '  %s\n' "${missing_results[@]}" >&2
    echo "Run ./run_convergence_5kv.sh, then retry this script." >&2
    exit 1
fi

# Opening every file catches truncated HDF5 output before any analysis starts.
# The recorded voltage and particle count also guard against mixing an older
# 3 kV run or an incorrectly sampled run into this named 5 kV study.
python3 - "${expected_results[@]}" <<'PY'
from pathlib import Path
import re
import sys

import h5py
import numpy as np

problems = []
project_source_hashes = set()
for argument in sys.argv[1:]:
    path = Path(argument)
    match = re.search(r"^(historical|generated)_n([0-9]+)_seed([0-9]+)$", path.parent.name)
    if match is None:
        problems.append(f"{path}: cannot infer source, count, and seed from directory name")
        continue

    expected_source = "file" if match.group(1) == "historical" else "generated"
    expected_particles = int(match.group(2))
    expected_seed = int(match.group(3))
    try:
        with h5py.File(path, "r") as h5:
            required = (
                "profiles/turn",
                "profiles/counts_macro_equivalent",
                "final_particles/zeta_m",
                "rf_program",
            )
            absent = [name for name in required if name not in h5]
            if absent:
                problems.append(f"{path}: missing {', '.join(absent)}")
                continue

            actual_particles = int(h5.attrs["sample_particles"])
            rf25_initial_v = float(h5["rf_program"].attrs["rf25_initial_v"])
            source_kind = str(h5.attrs["input_source_kind"])
            project_source_hash = str(h5.attrs["project_source_sha256"])
            project_source_hashes.add(project_source_hash)
            turns = np.asarray(h5["profiles/turn"])
            counts = np.asarray(h5["profiles/counts_macro_equivalent"])
            outside = np.asarray(h5["profiles/outside_macro_equivalent"])
            represented = counts.sum(axis=1, dtype=np.float64) + outside
            state = np.asarray(h5["final_particles/state"])
            zeta = np.asarray(h5["final_particles/zeta_m"])
            ptau = np.asarray(h5["final_particles/ptau"])
            alive_finite = int(
                np.count_nonzero(
                    (state > 0) & np.isfinite(zeta) & np.isfinite(ptau)
                )
            )
            if actual_particles != expected_particles:
                problems.append(
                    f"{path}: sample_particles={actual_particles}, "
                    f"expected {expected_particles}"
                )
            if rf25_initial_v != 5000.0:
                problems.append(
                    f"{path}: rf25_initial_v={rf25_initial_v:g} V, expected 5000 V"
                )
            if int(h5.attrs["n_turns"]) != 8083 or int(h5.attrs["record_every"]) != 10:
                problems.append(f"{path}: expected 8083 turns and record_every=10")
            if turns.size != 810 or int(turns[-1]) != 8083:
                problems.append(f"{path}: expected 810 profiles ending at turn 8083")
            if int(h5.attrs["random_seed"]) != expected_seed:
                problems.append(f"{path}: stored input seed does not match directory")
            if source_kind != expected_source:
                problems.append(
                    f"{path}: input_source_kind={source_kind}, expected {expected_source}"
                )
            if len(project_source_hash) != 64:
                problems.append(f"{path}: invalid project source hash")
            if not np.allclose(represented, 1_075_200, rtol=1e-6, atol=1e-3):
                problems.append(f"{path}: represented population is not conserved")
            if alive_finite != expected_particles:
                problems.append(
                    f"{path}: only {alive_finite}/{expected_particles} particles survive finite"
                )
            if expected_source == "file" and str(h5.attrs["input_sha256"]) != (
                "3b591b8381cb7cf57a88d7db08d535f2de182720d9fcb3c335800d7b3283e7e1"
            ):
                problems.append(f"{path}: historical input SHA-256 mismatch")
    except (KeyError, OSError, TypeError, ValueError) as exc:
        problems.append(f"{path}: {exc}")

if len(project_source_hashes) != 1:
    problems.append(
        "matrix mixes executable source hashes: " + ", ".join(sorted(project_source_hashes))
    )

if problems:
    print("HDF5 validation failed:", file=sys.stderr)
    for problem in problems:
        print(f"  {problem}", file=sys.stderr)
    raise SystemExit(1)

print(
    f"Validated {len(sys.argv) - 1} complete 5 kV HDF5 results; "
    f"tracking source SHA-256 {next(iter(project_source_hashes))}."
)
PY

mkdir -p "$analysis_root"

# Compare the historical and generated inputs at each point in the particle
# ladder. The turn-profile comparison and endpoint-only comparison answer
# complementary questions, so both are retained in the same result directory.
for particle_count in "${particle_counts[@]}"; do
    historical="$output_root/historical_n${particle_count}_seed202208/rf_only_results.h5"
    generated="$output_root/generated_n${particle_count}_seed202208/rf_only_results.h5"
    comparison_dir="$analysis_root/matched_n${particle_count}_seed202208"

    python3 main.py compare-runs \
        "$historical" \
        "$generated" \
        --output-dir "$comparison_dir"
    python3 main.py compare-final-time \
        "$historical" \
        "$generated" \
        --output-dir "$comparison_dir"
done

# Measure particle-count convergence separately for each input source. The
# full 1,075,200-particle endpoint is always the reference distribution.
for source_kind in historical generated; do
    full_result="$output_root/${source_kind}_n1075200_seed202208/rf_only_results.h5"
    for particle_count in 33600 100800 336000; do
        reduced_result="$output_root/${source_kind}_n${particle_count}_seed202208/rf_only_results.h5"
        comparison_dir="$analysis_root/convergence_to_full/${source_kind}_n${particle_count}_to_n1075200"

        python3 main.py compare-final-time \
            "$full_result" \
            "$reduced_result" \
            --output-dir "$comparison_dir"
    done
done

# Seed 202208 is the historical reference. All five generated seeds are the
# candidates, while the other four historical resamples establish the
# finite-sampling control envelope at the same 100,800-particle count.
screen_args=(
    --reference "$output_root/historical_n100800_seed202208/rf_only_results.h5"
)
for seed in 202208 "${extra_seeds[@]}"; do
    screen_args+=(
        --generated "$output_root/generated_n100800_seed${seed}/rf_only_results.h5"
    )
done
for seed in "${extra_seeds[@]}"; do
    screen_args+=(
        --control "$output_root/historical_n100800_seed${seed}/rf_only_results.h5"
    )
done
python3 main.py summarize-screen \
    "${screen_args[@]}" \
    --output "$analysis_root/screen_n100800_seeds202208_to_202212/rf_screen_summary.json"

# The full historical run is the standard visual reference for this matrix.
python3 main.py plot \
    "$output_root/historical_n1075200_seed202208/rf_only_results.h5" \
    --output-dir "$analysis_root/full_historical_n1075200_seed202208_plots"

python3 make_convergence_summary.py

echo "Convergence analysis completed under $analysis_root"
