#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${1:-${SCRIPT_DIR}}"
WORK_DIR="$(cd "${WORK_DIR}" && pwd)"

if [[ -f "${WORK_DIR}/driver_ysu_fortran" ]]; then
  FORT_BIN="${WORK_DIR}/driver_ysu_fortran"
  CPP_BIN="${WORK_DIR}/driver_ysu_cpu"
  RUN_DIR="${WORK_DIR}"
else
  FORT_BIN="${WORK_DIR}/build/driver_ysu_fortran"
  CPP_BIN="${WORK_DIR}/build/driver_ysu_cpu"
  RUN_DIR="${WORK_DIR}/build"
  mkdir -p "${RUN_DIR}"
fi

cd "${RUN_DIR}"

REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
HEAVY="${REPO_ROOT}/tools/heavy.sh"

INPUT_BIN="${RUN_DIR}/ysu_input.bin"
FORT_OUT="${RUN_DIR}/ysu_output_fortran.bin"
CPP_OUT="${RUN_DIR}/ysu_output_cpp.bin"

echo "=== Running YSU PBL Differential Test (2500 cols x 50 levels) ==="
python3 "${SCRIPT_DIR}/generate_inputs.py" --cols 2500 --nz 50 --dt 60.0 --out "${INPUT_BIN}"
"${HEAVY}" "${FORT_BIN}" "${INPUT_BIN}" "${FORT_OUT}"
"${HEAVY}" "${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}"
python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPP_OUT}" --cols 2500 --nz 50 --tol 1e-12
echo "=== YSU Differential Test PASSED ==="
