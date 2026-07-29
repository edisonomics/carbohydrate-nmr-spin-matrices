#!/usr/bin/env bash
# Load the local repository paths for the NMR spin-matrix workflows.
# Usage: source ./setup_paths.sh

_EDISON_SETUP_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

export EDISON_REPO_ROOT="${_EDISON_SETUP_DIR}"
export SPINACH_ROOT="${EDISON_REPO_ROOT}/lib/Spinach-2.10.1"
export EDISON_DATA_ROOT="${EDISON_REPO_ROOT}/data"
export EDISON_OUTPUT_ROOT="${EDISON_REPO_ROOT}/outputs"
export ALANINE_SRC="${EDISON_REPO_ROOT}/src/alanine"
export SUCROSE_SRC="${EDISON_REPO_ROOT}/src/sucrose"

if [ ! -d "${SPINACH_ROOT}" ]; then
    echo "Warning: SPINACH_ROOT does not exist: ${SPINACH_ROOT}" >&2
fi

mkdir -p "${EDISON_OUTPUT_ROOT}"

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
    echo "This script must be sourced so its exports persist in your shell."
    echo "Run: source ${EDISON_REPO_ROOT}/setup_paths.sh"
    exit 1
fi

echo "EDISON_REPO_ROOT=${EDISON_REPO_ROOT}"
echo "SPINACH_ROOT=${SPINACH_ROOT}"
echo "EDISON_DATA_ROOT=${EDISON_DATA_ROOT}"
echo "EDISON_OUTPUT_ROOT=${EDISON_OUTPUT_ROOT}"
echo "ALANINE_SRC=${ALANINE_SRC}"
echo "SUCROSE_SRC=${SUCROSE_SRC}"
