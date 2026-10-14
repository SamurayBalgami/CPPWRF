#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${1:-${SCRIPT_DIR}}"
WORK_DIR="$(cd "${WORK_DIR}" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${WORK_DIR}"

FORT_BIN="${WORK_DIR}/fortran_driver_thompson"
CPP_BIN="${WORK_DIR}/driver_thompson_cpu"
INPUT_BIN="${WORK_DIR}/input_thompson.bin"
FORT_OUT="${WORK_DIR}/output_fortran.bin"
CPP_OUT="${WORK_DIR}/output_cpp.bin"

echo "=== Running Thompson Microphysics Differential Test ==="
python3 "${SCRIPT_DIR}/generate_inputs.py"
"${REPO_ROOT}/tools/heavy.sh" "${FORT_BIN}" "${INPUT_BIN}" "${FORT_OUT}"
if [ -x "${WORK_DIR}/test_tables" ]; then
  "${REPO_ROOT}/tools/heavy.sh" "${WORK_DIR}/test_tables" "${WORK_DIR}/tables_fortran.bin"
fi
"${REPO_ROOT}/tools/heavy.sh" "${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}"
python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPP_OUT}" --tol 1e-12
if [ -x "${WORK_DIR}/driver_thompson_cuda_dp" ]; then
  "${REPO_ROOT}/tools/heavy.sh" "${WORK_DIR}/driver_thompson_cuda_dp" "${INPUT_BIN}" "${WORK_DIR}/output_cpp_cuda_dp.bin"
  python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${WORK_DIR}/output_cpp_cuda_dp.bin" --tol 1e-12
fi
echo "=== Thompson Microphysics Differential Test PASSED ==="
