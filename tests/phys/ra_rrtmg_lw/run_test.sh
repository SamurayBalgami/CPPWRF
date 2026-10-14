#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${1:-${SCRIPT_DIR}}"
WORK_DIR="$(cd "${WORK_DIR}" && pwd)"

if [[ -f "${WORK_DIR}/driver_rrtmg_lw_fortran" ]]; then
  FORT_BIN="${WORK_DIR}/driver_rrtmg_lw_fortran"
  CPP_BIN="${WORK_DIR}/driver_rrtmg_lw_cpu"
  [[ -f "${CPP_BIN}" ]] || CPP_BIN="${WORK_DIR}/driver_rrtmg_lw_cpp"
  RUN_DIR="${WORK_DIR}"
else
  FORT_BIN="${WORK_DIR}/build/driver_rrtmg_lw_fortran"
  CPP_BIN="${WORK_DIR}/build/driver_rrtmg_lw_cpu"
  [[ -f "${CPP_BIN}" ]] || CPP_BIN="${WORK_DIR}/build/driver_rrtmg_lw_cpp"
  RUN_DIR="${WORK_DIR}/build"
  mkdir -p "${RUN_DIR}"
fi

cd "${RUN_DIR}"

REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
HEAVY="${REPO_ROOT}/tools/heavy.sh"

# Ensure data file is linked
ln -sf /home/ostriquetrum/wrf_ref/WRF/run/RRTMG_LW_DATA_DBL "${RUN_DIR}/RRTMG_LW_DATA"

INPUT_BIN="${RUN_DIR}/input_rrtmg_lw.bin"
FORT_OUT="${RUN_DIR}/output_fortran.bin"
CPP_OUT="${RUN_DIR}/output_cpp.bin"
TABLES_BIN="${RUN_DIR}/rrtmg_lw_tables.bin"

if [[ ! -f "${TABLES_BIN}" ]]; then
  ln -sf "${SCRIPT_DIR}/rrtmg_lw_tables.bin" "${TABLES_BIN}"
fi

echo "=== Running RRTMG LW Differential Test (2500 cols) ==="
python3 "${SCRIPT_DIR}/generate_inputs.py" --cols 2500 --out "${INPUT_BIN}"
"${HEAVY}" "${FORT_BIN}" "${INPUT_BIN}" "${FORT_OUT}"
"${HEAVY}" "${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}" "${TABLES_BIN}"
python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPP_OUT}" 1e-6
echo "=== RRTMG LW Differential Test PASSED ==="
