#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${1:-${SCRIPT_DIR}}"
WORK_DIR="$(cd "${WORK_DIR}" && pwd)"

if [[ -f "${WORK_DIR}/driver_rrtmg_sw_fortran" ]]; then
  FORT_BIN="${WORK_DIR}/driver_rrtmg_sw_fortran"
  CPP_BIN="${WORK_DIR}/driver_rrtmg_sw_cpu"
  [[ -f "${CPP_BIN}" ]] || CPP_BIN="${WORK_DIR}/driver_rrtmg_sw_cpp"
  RUN_DIR="${WORK_DIR}"
else
  FORT_BIN="${WORK_DIR}/build/driver_rrtmg_sw_fortran"
  CPP_BIN="${WORK_DIR}/build/driver_rrtmg_sw_cpu"
  [[ -f "${CPP_BIN}" ]] || CPP_BIN="${WORK_DIR}/driver_rrtmg_sw_cpp"
  RUN_DIR="${WORK_DIR}/build"
  mkdir -p "${RUN_DIR}"
fi

cd "${RUN_DIR}"

REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
HEAVY="${REPO_ROOT}/tools/heavy.sh"

# Ensure data files are linked
ln -sf /home/ostriquetrum/wrf_ref/WRF/run/RRTMG_SW_DATA_DBL "${RUN_DIR}/RRTMG_SW_DATA"
ln -sf /home/ostriquetrum/wrf_ref/WRF/run/RRTMG_LW_DATA_DBL "${RUN_DIR}/RRTMG_LW_DATA"

INPUT_BIN="${RUN_DIR}/input_rrtmg_sw.bin"
FORT_OUT="${RUN_DIR}/output_fortran.bin"
CPP_OUT="${RUN_DIR}/output_cpp.bin"
TABLES_BIN="${RUN_DIR}/rrtmg_sw_tables.bin"

echo "=== Running RRTMG SW Differential Test (2500 cols) ==="
python3 "${SCRIPT_DIR}/generate_inputs.py" --cols 2500 --out "${INPUT_BIN}"
"${HEAVY}" "${FORT_BIN}"
"${HEAVY}" "${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}" "${TABLES_BIN}"
python3 "${SCRIPT_DIR}/compare.py" --fortran "${FORT_OUT}" --cpp "${CPP_OUT}" --input "${INPUT_BIN}" --tol 5e-12
echo "=== RRTMG SW Differential Test PASSED ==="
