#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Kain-Fritsch (cu_physics = 1) tests: WRF v4.8.0's kf_lutab tables and KF_eta_CPS (module_cu_kfeta.o of the -r8
# build, WRF_DIR) against the C++ tables (test_tables_kf) and the C++ KF (driver_kf_cpu of the CMake build),
# tolerance 0.
#   run_test.sh BUILD_DIR      (BUILD_DIR: the directory of driver_kf_cpu; the work files go there)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="$(cd "${1:?usage: run_test.sh BUILD_DIR}" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
. "${SCRIPT_DIR}/../wrf_fortran.sh"
wrf_need r8
cd "${WORK_DIR}"

FORT_BIN="${WORK_DIR}/fortran_driver_kf"
CPP_BIN="${WORK_DIR}/driver_kf_cpu"
INPUT_BIN="${WORK_DIR}/input_kf.bin"
FORT_OUT="${WORK_DIR}/output_fortran.bin"
CPP_OUT="${WORK_DIR}/output_cpp.bin"
wrf_fortran "${FORT_BIN}" "${SCRIPT_DIR}/fortran/driver_kf.f90"
wrf_fortran "${WORK_DIR}/dump_tables_kf" "${SCRIPT_DIR}/dump_tables_f.f90" "${SCRIPT_DIR}/dump_tables_c.c"

echo "=== Kain-Fritsch lookup tables (kf_lutab) ==="
rm -f tables_fortran.bin
"${WORK_DIR}/dump_tables_kf"
"${WORK_DIR}/test_tables_kf" tables_fortran.bin

# columns, levels, kfeta_trigger, moisture species (ice: QI/QS/QR, kessler: QR, none: mp_physics 0);
# 2499 x 51: the odd KL tail of WRF's vectorised loops
for cfg in "2500 50 1 ice" "2499 51 1 ice" "2500 50 3 ice" "2500 50 1 kessler" "2500 50 1 none"; do
  set -- $cfg
  echo "=== Kain-Fritsch Cumulus Differential Test ($1 cols x $2 levels, kfeta_trigger $3, species $4) ==="
  python3 "${SCRIPT_DIR}/generate_inputs.py" --ncols "$1" --nlevels "$2" --trigger "$3" --mp "$4"
  "${REPO_ROOT}/tools/heavy.sh" "${FORT_BIN}"
  "${REPO_ROOT}/tools/heavy.sh" "${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}"
  python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPP_OUT}" --tol 0
done
echo "=== Kain-Fritsch Cumulus Differential Test PASSED ==="
