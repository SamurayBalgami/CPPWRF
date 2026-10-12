#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${1:-${SCRIPT_DIR}}"
WORK_DIR="$(cd "${WORK_DIR}" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${WORK_DIR}"

FORT_BIN="${WORK_DIR}/fortran_driver_kf"
CPP_BIN="${WORK_DIR}/driver_kf_cpu"
INPUT_BIN="${WORK_DIR}/input_kf.bin"
FORT_OUT="${WORK_DIR}/output_fortran.bin"
CPP_OUT="${WORK_DIR}/output_cpp.bin"

echo "=== Running Kain-Fritsch Cumulus Differential Test ==="
python3 "${SCRIPT_DIR}/generate_inputs.py"
"${REPO_ROOT}/tools/heavy.sh" "${FORT_BIN}"
"${REPO_ROOT}/tools/heavy.sh" "${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}"
python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPP_OUT}" --tol 1e-12
echo "=== Kain-Fritsch Cumulus Differential Test PASSED ==="
