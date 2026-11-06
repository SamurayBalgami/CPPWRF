#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Revised MM5 surface layer (sf_sfclay_physics = 1) column test: WRF v4.8.0's sfclayrev (module_sf_sfclayrev.o,
# sf_sfclayrev.o of the -r8 build, WRF_DIR) against the C++ sfclayrev (driver_sfclayrev_cpu of the CMake build),
# tolerance 0.
#   run_test.sh BUILD_DIR      (BUILD_DIR: the directory of driver_sfclayrev_cpu; the work files go there)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="$(cd "${1:?usage: run_test.sh BUILD_DIR}" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
HEAVY="${REPO_ROOT}/tools/heavy.sh"
. "${SCRIPT_DIR}/../wrf_fortran.sh"
wrf_need r8
cd "${RUN_DIR}"

FORT_BIN="${RUN_DIR}/driver_sfclayrev_fortran"
CPP_BIN="${RUN_DIR}/driver_sfclayrev_cpu"
wrf_fortran "${FORT_BIN}" "${SCRIPT_DIR}/driver_sfclayrev.f90"

INPUT_BIN="${RUN_DIR}/sfclayrev_input.bin"
FORT_OUT="${RUN_DIR}/sfclayrev_output_fortran.bin"
CPP_OUT="${RUN_DIR}/sfclayrev_output_cpp.bin"

# columns, isfflx, isftcflx, iz0tlnd, kim_tofd (WRF 4.8 turbulent orographic form drag); 2499 columns: the odd i
# tail of WRF's vectorised loops
for cfg in "2500 1 0 0 0" "2500 1 0 0 1" "2500 0 0 0 0" "2500 1 1 1 0" "2500 1 2 2 0" "2499 1 0 0 0"; do
  set -- $cfg
  echo "=== SFCLAYREV Differential Test: $1 columns, isfflx $2, isftcflx $3, iz0tlnd $4, kim_tofd $5 ==="
  python3 "${SCRIPT_DIR}/generate_inputs.py" --cols "$1" --isfflx "$2" --isftcflx "$3" --iz0tlnd "$4" --out "${INPUT_BIN}"
  CPPWRF_TEST_KIM_TOFD=$5 "${HEAVY}" "${FORT_BIN}" "${INPUT_BIN}" "${FORT_OUT}"
  CPPWRF_TEST_KIM_TOFD=$5 "${HEAVY}" "${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}"
  python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPP_OUT}" --tol 0
done
echo "=== SFCLAYREV Differential Test PASSED ==="
