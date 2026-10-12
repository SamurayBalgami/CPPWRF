#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${1:-${SCRIPT_DIR}}"
WORK_DIR="$(cd "${WORK_DIR}" && pwd)"
cd "${WORK_DIR}"

FORT_BIN="${WORK_DIR}/driver_kessler_fortran"
CPP_BIN="${WORK_DIR}/driver_kessler_cpu"
INPUT_BIN="${WORK_DIR}/kessler_input.bin"
FORT_OUT="${WORK_DIR}/kessler_output_fortran.bin"
CPP_OUT="${WORK_DIR}/kessler_output_cpp.bin"

echo "=== Running Kessler Microphysics Differential Test ==="
python3 "${SCRIPT_DIR}/generate_inputs.py" --ncols 2500 --nz 50 --dt 12.0 -o "${INPUT_BIN}"
"${FORT_BIN}" "${INPUT_BIN}" "${FORT_OUT}"
"${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}"
python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPP_OUT}" --tol 1e-12
echo "=== Kessler Differential Test PASSED ==="
