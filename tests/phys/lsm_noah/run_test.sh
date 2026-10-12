#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${1:-${SCRIPT_DIR}}"
WORK_DIR="$(cd "${WORK_DIR}" && pwd)"
cd "${WORK_DIR}"

FORT_BIN="${WORK_DIR}/driver_noah_fortran"
CPP_BIN="${WORK_DIR}/driver_noah_cpu"
INPUT_BIN="${WORK_DIR}/noah_input.bin"
FORT_OUT="${WORK_DIR}/noah_output_fortran.bin"
CPP_OUT="${WORK_DIR}/noah_output_cpp.bin"

echo "=== Running Noah LSM Differential Test ==="
python3 "${SCRIPT_DIR}/generate_inputs.py" "${INPUT_BIN}"
"${FORT_BIN}" "${INPUT_BIN}" "${FORT_OUT}" "${SCRIPT_DIR}"
"${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}" "${SCRIPT_DIR}"
python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPP_OUT}" --tol 1e-12
echo "=== Noah LSM Differential Test PASSED ==="
