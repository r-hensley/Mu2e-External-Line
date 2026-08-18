#!/usr/bin/env bash
# Download and verify the Git-ignored files used by the Prebys study.

set -euo pipefail

readonly SERVER_BASE="https://prebys.physics.ucdavis.edu/misc/AAAreadme/g4beamline_study"
readonly EXPORT_ARCHIVE="g4beamline_study_export.tgz"
readonly EXPORT_SHA256="7e5e4a30e5a85688a44154ab70fc5ded54f2ddfd8c53b2e62b3af3b09a8a7137"
readonly OLD_ARCHIVE="g4beamline_study_20170103.tgz"
readonly OLD_SHA256="7d6c5e6ee0c731be57e7f3b23ba6c164dcd0c5f3c4adf468a37a97e24f784350"

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly STUDY_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
readonly INPUT_DIR="${SCRIPT_DIR}/large_inputs"
readonly OUTPUT_DIR="${STUDY_DIR}/reference_results/large_outputs"
readonly TEMP_DIR="$(mktemp -d /tmp/prebys-g4bl-download-XXXXXX)"

PARTIAL_PATH=""

cleanup() {
    # Remove temporary downloads and any interrupted destination copy.

    if [[ -n "${PARTIAL_PATH}" ]]; then
        rm -f -- "${PARTIAL_PATH}"
    fi
    rm -rf -- "${TEMP_DIR}"
}
trap cleanup EXIT

require_command() {
    # Fail with a useful message when a required command is unavailable.

    local command_name="$1"
    if ! command -v "${command_name}" >/dev/null 2>&1; then
        echo "Missing required command: ${command_name}" >&2
        exit 1
    fi
}

sha256_of() {
    # Print only the SHA-256 digest for a local file.

    sha256sum -- "$1" | awk '{print $1}'
}

verify_file() {
    # Verify one file against an expected SHA-256 digest.

    local path="$1"
    local expected_sha256="$2"
    local actual_sha256

    actual_sha256="$(sha256_of "${path}")"
    if [[ "${actual_sha256}" != "${expected_sha256}" ]]; then
        echo "SHA-256 mismatch: ${path}" >&2
        echo "  expected: ${expected_sha256}" >&2
        echo "  actual:   ${actual_sha256}" >&2
        return 1
    fi
}

destination_is_ready() {
    # Return success for a correct destination and reject a conflicting one.

    local destination="$1"
    local expected_sha256="$2"

    if [[ ! -e "${destination}" ]]; then
        return 1
    fi
    if [[ ! -f "${destination}" || -L "${destination}" ]]; then
        echo "Refusing non-regular destination: ${destination}" >&2
        exit 1
    fi
    if ! verify_file "${destination}" "${expected_sha256}"; then
        echo "Refusing to overwrite the mismatched local file." >&2
        echo "Move it aside after investigating, then run this script again." >&2
        exit 1
    fi

    echo "Already verified: ${destination}"
    return 0
}

install_verified() {
    # Atomically install a verified file without overwriting existing data.

    local source="$1"
    local destination="$2"
    local expected_sha256="$3"

    verify_file "${source}" "${expected_sha256}"
    if destination_is_ready "${destination}" "${expected_sha256}"; then
        return 0
    fi

    PARTIAL_PATH="${destination}.partial.$$"
    cp -p -- "${source}" "${PARTIAL_PATH}"
    verify_file "${PARTIAL_PATH}" "${expected_sha256}"
    mv --no-clobber -- "${PARTIAL_PATH}" "${destination}"
    if [[ -e "${PARTIAL_PATH}" ]]; then
        rm -f -- "${PARTIAL_PATH}"
        PARTIAL_PATH=""
        destination_is_ready "${destination}" "${expected_sha256}"
        return 0
    fi
    PARTIAL_PATH=""
    echo "Installed: ${destination}"
}

download() {
    # Download a URL to a temporary file and preserve its server timestamp.

    local url="$1"
    local destination="$2"

    echo "Downloading: ${url}"
    curl --fail --location --remote-time \
        --retry 3 --retry-delay 2 --connect-timeout 20 \
        --output "${destination}" "${url}"
}

archive_member_destination() {
    # Return the curated destination for one selected export member.

    case "$1" in
        ExtrBeamDSC-magL12noDiff.dat|ExtrBeamTransToCMAG-D.dat)
            printf '%s/%s\n' "${INPUT_DIR}" "$1"
            ;;
        us_phase_space_col.pdf|us_phase_space_nocol.pdf)
            printf '%s/%s\n' "${OUTPUT_DIR}" "$1"
            ;;
        *)
            echo "Unrecognized selected archive member: $1" >&2
            exit 1
            ;;
    esac
}

extract_required_archive_files() {
    # Fetch both server archives, but extract only four required large files.

    local export_path="${TEMP_DIR}/${EXPORT_ARCHIVE}"
    local old_path="${TEMP_DIR}/${OLD_ARCHIVE}"
    local stage_dir="${TEMP_DIR}/archive-files"
    local member
    local basename
    local destination
    local expected_sha256
    local needs_export=0

    # The 2017 archive is checked for provenance but contains no ignored large
    # dependency. The later export supplies exactly the four members below.
    declare -A member_sha256=(
        [ExtrBeamDSC-magL12noDiff.dat]="0136b4eb07135ca7c0147fa2dc2d26b5c7728e79229327652e375f1e50f3d94b"
        [ExtrBeamTransToCMAG-D.dat]="cf7c836f923dd814f87ca7a6645edff3a02c0c9c05bcf021186002ed9e6b8149"
        [us_phase_space_col.pdf]="27d011577131a77d5e5e9066259773b2fbacbc4e8e9045f29c51ac026a28f2fa"
        [us_phase_space_nocol.pdf]="96112a39ea8db2427cef5e894cf3a88e38eed66210659dd1103a90c4779b64e9"
    )

    for basename in "${!member_sha256[@]}"; do
        destination="$(archive_member_destination "${basename}")"
        if ! destination_is_ready "${destination}" "${member_sha256[${basename}]}"; then
            needs_export=1
        fi
    done
    if (( needs_export == 0 )); then
        return 0
    fi

    download "${SERVER_BASE}/${EXPORT_ARCHIVE}" "${export_path}"
    verify_file "${export_path}" "${EXPORT_SHA256}"
    download "${SERVER_BASE}/${OLD_ARCHIVE}" "${old_path}"
    verify_file "${old_path}" "${OLD_SHA256}"
    echo "Verified both remote .tgz archives. The 2017 archive needs no extraction."

    mkdir -p -- "${stage_dir}"
    for basename in \
        ExtrBeamDSC-magL12noDiff.dat \
        ExtrBeamTransToCMAG-D.dat \
        us_phase_space_col.pdf \
        us_phase_space_nocol.pdf
    do
        member="g4beamline_study_export/${basename}"
        destination="$(archive_member_destination "${basename}")"
        expected_sha256="${member_sha256[${basename}]}"

        if destination_is_ready "${destination}" "${expected_sha256}"; then
            continue
        fi
        tar --extract --gzip --file "${export_path}" \
            --directory "${stage_dir}" --no-same-owner -- "${member}"
        install_verified "${stage_dir}/${member}" "${destination}" "${expected_sha256}"
    done
}

download_standalone_inputs() {
    # Download the five required ROOT inputs that are not in either tarball.

    local filename
    local expected_sha256
    local temporary_path
    declare -A input_sha256=(
        [exttracks_rot_1800000.root]="7df4c14fcb651c9259293b33ca9e5ae645d9ad8a6b61a13095af64c36f055edb"
        [exttracks_rot_e30_e40_1800000.root]="583286959bb1cdd94df638a62c52c10412e8e3f6a852d0cc9bb4a2fbf130cf74"
        [mu2e_downstream.root]="7e5cefcb2fad3f2eb615dbfa4e9d47a804a717ee541f6d2eedef96cd748ffc0a"
        [mu2e_upstream.root]="3853c89068ffd225a358b52309d87f812cf8d8034ecddd8fade34ff75f848968"
        [mu2e_upstream_gaussian.root]="50404b83ff69620be6bfdb06ab0ba8f4d4e403bddb38714461a1efd5f73c4196"
    )

    for filename in \
        exttracks_rot_1800000.root \
        exttracks_rot_e30_e40_1800000.root \
        mu2e_downstream.root \
        mu2e_upstream.root \
        mu2e_upstream_gaussian.root
    do
        expected_sha256="${input_sha256[${filename}]}"
        if destination_is_ready "${INPUT_DIR}/${filename}" "${expected_sha256}"; then
            continue
        fi
        temporary_path="${TEMP_DIR}/${filename}"
        download "${SERVER_BASE}/${filename}" "${temporary_path}"
        install_verified "${temporary_path}" "${INPUT_DIR}/${filename}" "${expected_sha256}"
    done
}

main() {
    # Restore all nine ignored large files and verify the resulting dataset.

    require_command awk
    require_command cp
    require_command curl
    require_command mktemp
    require_command mv
    require_command rm
    require_command sha256sum
    require_command tar

    mkdir -p -- "${INPUT_DIR}" "${OUTPUT_DIR}"
    extract_required_archive_files
    download_standalone_inputs
    echo "All required Git-ignored Prebys data files are present and verified."
}

main "$@"
