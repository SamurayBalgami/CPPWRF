#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Noah LSM (sf_surface_physics = 2) column test: WRF v4.8.0's SFLX with the tables of its SOIL_VEG_GEN_PARM
# (module_sf_noahlsm.o, module_sf_noahdrv.o of the -r8 build, WRF_DIR; WRF's run/ tables) against the C++
# NoahKernel (driver_noah_cpu of the CMake build), tolerance 0.
#   run_test.sh BUILD_DIR      (BUILD_DIR: the directory of driver_noah_cpu; the work files go there)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="$(cd "${1:?usage: run_test.sh BUILD_DIR}" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
HEAVY="${REPO_ROOT}/tools/heavy.sh"
. "${SCRIPT_DIR}/../wrf_fortran.sh"
wrf_need r8
cd "${WORK_DIR}"
for t in VEGPARM.TBL SOILPARM.TBL GENPARM.TBL; do cp -L "${W}/run/${t}" .; done

FORT_BIN="${WORK_DIR}/driver_noah_fortran"
CPP_BIN="${WORK_DIR}/driver_noah_cpu"
INPUT_BIN="${WORK_DIR}/noah_input.bin"
FORT_OUT="${WORK_DIR}/noah_output_fortran.bin"
CPP_OUT="${WORK_DIR}/noah_output_cpp.bin"
wrf_fortran "${FORT_BIN}" "${SCRIPT_DIR}/fortran/driver_noah.f90"

python3 "${SCRIPT_DIR}/generate_inputs.py" "${INPUT_BIN}"
for opt in 1 2; do  # opt_thcnd
  echo "=== Running Noah LSM Differential Test (opt_thcnd = ${opt}) ==="
  CPPWRF_TEST_OPT_THCND=${opt} "${HEAVY}" "${FORT_BIN}" "${INPUT_BIN}" "${FORT_OUT}"
  CPPWRF_TEST_OPT_THCND=${opt} "${HEAVY}" "${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}" "${WORK_DIR}"
  python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPP_OUT}" --tol 0
done
echo "=== Noah LSM Differential Test PASSED ==="
