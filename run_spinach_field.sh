#!/usr/bin/env bash
# Student-facing launcher for Spinach carbohydrate fields.
# Usage: ./run_spinach_field.sh MOLECULE [FIELD_MHZ] [MATRIX_FILE]
# Candidate matrix for all fields: ./run_spinach_field.sh sucrose --matrix outputs/sucrose/sucrose_matrix_fit_600_900.txt

set -euo pipefail

repo_root="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

if [[ $# -lt 1 || $# -gt 3 ]]; then
    echo "Usage: ./run_spinach_field.sh MOLECULE [FIELD_MHZ] [MATRIX_FILE]" >&2
    echo "Example: ./run_spinach_field.sh glucose" >&2
    exit 2
fi

molecule="$1"
field="${2:-}"
matrix_override="${3:-}"
if [[ "$field" == "--matrix" ]]; then
    matrix_override="${3:-}"
    field=""
fi
if [[ -n "$matrix_override" && "$matrix_override" != /* ]]; then
    matrix_override="$repo_root/$matrix_override"
fi
spinach_root="${SPINACH_ROOT:-${repo_root}/lib/Spinach-2.10.1}"

if [[ ! -d "$spinach_root" ]]; then
    echo "Spinach was not found at: $spinach_root" >&2
    echo "Install/copy Spinach into lib/Spinach-2.10.1 or set SPINACH_ROOT." >&2
    exit 1
fi

matlab_bin="${MATLAB_BIN:-}"
if [[ -z "$matlab_bin" ]]; then
    matlab_bin="$(command -v matlab || true)"
fi
if [[ -z "$matlab_bin" ]]; then
    for candidate in /Applications/MATLAB_*.app/bin/matlab; do
        if [[ -x "$candidate" ]]; then
            matlab_bin="$candidate"
            break
        fi
    done
fi
if [[ -z "$matlab_bin" || ! -x "$matlab_bin" ]]; then
    echo "MATLAB was not found. Open MATLAB once, or set MATLAB_BIN to its executable." >&2
    exit 1
fi

export EDISON_REPO_ROOT="$repo_root"
export SPINACH_ROOT="$spinach_root"

if [[ -n "$field" ]]; then
    if [[ -n "$matrix_override" ]]; then
        matlab_call="run_carbohydrate_spinach_field('${molecule}','${field}','${matrix_override}');"
    else
        matlab_call="run_carbohydrate_spinach_field('${molecule}','${field}');"
    fi
else
    if [[ -n "$matrix_override" ]]; then
        matlab_call="run_carbohydrate_spinach_all('${molecule}','${matrix_override}');"
    else
        matlab_call="run_carbohydrate_spinach_all('${molecule}');"
    fi
fi

"$matlab_bin" -batch "addpath('${repo_root}/src/common'); ${matlab_call}"
