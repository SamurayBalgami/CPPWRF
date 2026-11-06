#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# YSU PBL (bl_pbl_physics = 1) column test: WRF v4.8.0's ysu (module_bl_ysu.o, bl_ysu.o of the -r8 build,
# WRF_DIR) against the C++ YSU (driver_ysu_cpu of the CMake build), tolerance 0.
#   run_test.sh BUILD_DIR      (BUILD_DIR: the directory of driver_ysu_cpu; the work files go there)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="$(cd "${1:?usage: run_test.sh BUILD_DIR}" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
HEAVY="${REPO_ROOT}/tools/heavy.sh"
. "${SCRIPT_DIR}/../wrf_fortran.sh"
wrf_need r8
cd "${RUN_DIR}"

FORT_BIN="${RUN_DIR}/driver_ysu_fortran"
CPP_BIN="${RUN_DIR}/driver_ysu_cpu"
wrf_fortran "${FORT_BIN}" "${SCRIPT_DIR}/driver_ysu.f90"

INPUT_BIN="${RUN_DIR}/ysu_input.bin"
FORT_OUT="${RUN_DIR}/ysu_output_fortran.bin"
CPP_OUT="${RUN_DIR}/ysu_output_cpp.bin"

# columns, levels, ysu_topdown_pblmix (WRF default 1), flag_qc = flag_qi (0: no microphysics, as P1); 2499 x 41:
# the scalar tails of WRF's vectorised loops
for cfg in "2500 50 1 1" "2500 50 0 1" "2500 50 1 0" "2499 41 1 1"; do
  set -- $cfg
  echo "=== YSU Differential Test: $1 columns x $2 levels, ysu_topdown_pblmix $3, flag_qc/qi $4 ==="
  python3 "${SCRIPT_DIR}/generate_inputs.py" --cols "$1" --nz "$2" --dt 60.0 --topdown "$3" --flag-q "$4" --out "${INPUT_BIN}"
  "${HEAVY}" "${FORT_BIN}" "${INPUT_BIN}" "${FORT_OUT}"
  "${HEAVY}" "${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}"
  python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPP_OUT}" --cols "$1" --nz "$2" --tol 0
done
echo "=== YSU Differential Test PASSED ==="
