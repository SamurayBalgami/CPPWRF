#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${1:-${SCRIPT_DIR}}"
WORK_DIR="$(cd "${WORK_DIR}" && pwd)"

if [[ -f "${WORK_DIR}/driver_wsm6_fortran" ]]; then
  FORT_BIN="${WORK_DIR}/driver_wsm6_fortran"
  CPP_BIN="${WORK_DIR}/driver_wsm6_cpu"
  RUN_DIR="${WORK_DIR}"
else
  FORT_BIN="${WORK_DIR}/build/driver_wsm6_fortran"
  CPP_BIN="${WORK_DIR}/build/driver_wsm6_cpu"
  RUN_DIR="${WORK_DIR}/build"
  mkdir -p "${RUN_DIR}"
fi

cd "${RUN_DIR}"

REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
HEAVY="${REPO_ROOT}/tools/heavy.sh"

INPUT_BIN="${RUN_DIR}/wsm6_input.bin"
FORT_OUT="${RUN_DIR}/wsm6_output_fortran.bin"
CPP_OUT="${RUN_DIR}/wsm6_output_cpp.bin"

echo "=== Running WSM6 Microphysics Differential Test (Multi-loop dt=240s) ==="
python3 "${SCRIPT_DIR}/generate_inputs.py" --cols 2500 --nz 50 --dt 240.0 --out "${INPUT_BIN}"
"${HEAVY}" "${FORT_BIN}" "${INPUT_BIN}" "${FORT_OUT}"
"${HEAVY}" "${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}"
python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPP_OUT}" --tol 1e-12
echo "=== WSM6 Multi-Loop Differential Test PASSED ==="

echo "=== Running WSM6 Microphysics Differential Test (Single-loop dt=120s) ==="
INPUT_BIN_120="${RUN_DIR}/wsm6_input_120.bin"
FORT_OUT_120="${RUN_DIR}/wsm6_output_fortran_120.bin"
CPP_OUT_120="${RUN_DIR}/wsm6_output_cpp_120.bin"
python3 "${SCRIPT_DIR}/generate_inputs.py" --cols 2500 --nz 50 --dt 120.0 --out "${INPUT_BIN_120}"
"${HEAVY}" "${FORT_BIN}" "${INPUT_BIN_120}" "${FORT_OUT_120}"
"${HEAVY}" "${CPP_BIN}" "${INPUT_BIN_120}" "${CPP_OUT_120}"
python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT_120}" "${CPP_OUT_120}" --tol 1e-12
echo "=== WSM6 Single-Loop Differential Test PASSED ==="
