#!/usr/bin/env bash
# Reproduce the 5 kV RF-only convergence matrix used in results.md.
#
# Each historical/generated pair runs concurrently. Completed HDF5 files are
# skipped, so an interrupted study can be resumed by invoking this script again.

set -u
set -o pipefail

project_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$project_dir" || exit 1

output_root=outputs/convergence_5kv
input_file=1075200_init_dist.txt
input_sha256=3b591b8381cb7cf57a88d7db08d535f2de182720d9fcb3c335800d7b3283e7e1

if ! printf '%s  %s\n' "$input_sha256" "$input_file" | sha256sum --check --status; then
    echo "Historical input is missing or has the wrong SHA-256: $input_file" >&2
    exit 1
fi

validate_result() {
    local result_path=$1
    local expected_particles=$2
    local expected_source_kind=$3

    python3 - "$result_path" "$expected_particles" "$expected_source_kind" "$input_sha256" <<'PY'
from pathlib import Path
import sys

import h5py

path = Path(sys.argv[1])
expected_particles = int(sys.argv[2])
expected_source_kind = sys.argv[3]
historical_sha256 = sys.argv[4]

try:
    with h5py.File(path, "r") as h5:
        checks = {
            "particle count": int(h5.attrs["sample_particles"]) == expected_particles,
            "turn count": int(h5.attrs["n_turns"]) == 8083,
            "profile cadence": int(h5.attrs["record_every"]) == 10,
            "last recorded turn": int(h5["profiles/turn"][-1]) == 8083,
            "final coordinate length": len(h5["final_particles/zeta_m"]) == expected_particles,
            "5 kV RF start": float(h5["rf_program"].attrs["rf25_initial_v"]) == 5000.0,
            "source kind": str(h5.attrs["input_source_kind"]) == expected_source_kind,
            "source-tree hash": len(str(h5.attrs["project_source_sha256"])) == 64,
            "coordinate hash": len(str(h5.attrs["input_tracking_coordinate_sha256"])) == 64,
            "recorded invocation": bool(str(h5.attrs["invocation_command"])),
        }
        if expected_source_kind == "file":
            checks["historical input hash"] = (
                str(h5.attrs["input_sha256"]) == historical_sha256
            )
except (KeyError, OSError, TypeError, ValueError) as exc:
    print(f"{path}: cannot validate completed result: {exc}", file=sys.stderr)
    raise SystemExit(1)

failed = [name for name, passed in checks.items() if not passed]
if failed:
    print(f"{path}: failed completion checks: {', '.join(failed)}", file=sys.stderr)
    raise SystemExit(1)
PY
}

run_historical() {
    local particle_count=$1
    local seed=$2
    local output_dir="$output_root/historical_n${particle_count}_seed${seed}"

    if [[ -e "$output_dir/rf_only_results.h5" ]]; then
        if validate_result "$output_dir/rf_only_results.h5" "$particle_count" file; then
            echo "Skipping validated $output_dir"
            return 0
        fi
        echo "Refusing to overwrite an invalid existing result in $output_dir" >&2
        return 1
    fi

    mkdir -p "$output_dir"
    local max_particles=$particle_count
    if [[ $particle_count -eq 1075200 ]]; then
        max_particles=0
    fi

    python3 -u main.py run \
        --input "$input_file" \
        --max-particles "$max_particles" \
        --sampling-strata 168 \
        --seed "$seed" \
        --record-every 10 \
        --turns 8083 \
        --rf25-start-kv 5 \
        --no-figures \
        --output-dir "$output_dir" \
        >"$output_dir/run.log" 2>&1
    validate_result "$output_dir/rf_only_results.h5" "$particle_count" file
}

run_generated() {
    local particle_count=$1
    local seed=$2
    local output_dir="$output_root/generated_n${particle_count}_seed${seed}"
    local particles_per_microbunch=$((particle_count / 168))

    if [[ -e "$output_dir/rf_only_results.h5" ]]; then
        if validate_result "$output_dir/rf_only_results.h5" "$particle_count" generated; then
            echo "Skipping validated $output_dir"
            return 0
        fi
        echo "Refusing to overwrite an invalid existing result in $output_dir" >&2
        return 1
    fi

    if [[ $((particle_count % 168)) -ne 0 ]]; then
        echo "Generated particle count must be divisible by 168: $particle_count" >&2
        return 2
    fi

    mkdir -p "$output_dir"
    python3 -u main.py run \
        --particles-per-microbunch "$particles_per_microbunch" \
        --seed "$seed" \
        --record-every 10 \
        --turns 8083 \
        --rf25-start-kv 5 \
        --no-figures \
        --output-dir "$output_dir" \
        >"$output_dir/run.log" 2>&1
    validate_result "$output_dir/rf_only_results.h5" "$particle_count" generated
}

run_pair() {
    local particle_count=$1
    local seed=$2
    local historical_pid generated_pid historical_status generated_status

    echo "Starting matched pair: n=$particle_count seed=$seed"
    run_historical "$particle_count" "$seed" &
    historical_pid=$!
    run_generated "$particle_count" "$seed" &
    generated_pid=$!

    wait "$historical_pid"
    historical_status=$?
    wait "$generated_pid"
    generated_status=$?
    if [[ $historical_status -ne 0 || $generated_status -ne 0 ]]; then
        echo "Pair failed: n=$particle_count seed=$seed " \
             "(historical=$historical_status, generated=$generated_status)" >&2
        return 1
    fi
    echo "Completed matched pair: n=$particle_count seed=$seed"
}

# Seed 202208 supplies the particle-count ladder. The full-statistics pair is
# started first because it is the longest checkpoint.
for particle_count in 1075200 336000 100800 33600; do
    run_pair "$particle_count" 202208 || exit 1
done

# Four additional generated/reference-resampling pairs test seed stability at
# 100,800 particles, a compromise between tail statistics and local runtime.
for seed in 202209 202210 202211 202212; do
    run_pair 100800 "$seed" || exit 1
done

echo "All 5 kV convergence simulations completed under $output_root"
