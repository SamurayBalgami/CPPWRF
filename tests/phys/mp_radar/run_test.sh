#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# WSM6's radar reflectivity (do_radar_ref = 1) column test: WRF v4.8.0's mp_wsm6_init / radar_init and
# refl10cm_wsm6 (mp_wsm6.o, mp_radar.o of the -r8 -fdefault-double-8 build, WRF_R8D8_DIR: in the -r8 build their
# DOUBLE PRECISION is REAL*16) against the C++ radar_init + Wsm6ReflKernel (driver_radar_cpu of the CMake build),
# tolerance 0. Thompson's REFL_10CM is part of the Thompson column test.
#   run_test.sh BUILD_DIR      (BUILD_DIR: the directory of driver_radar_cpu; the work files go there)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="$(cd "${1:?usage: run_test.sh BUILD_DIR}" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
HEAVY="${REPO_ROOT}/tools/heavy.sh"
. "${SCRIPT_DIR}/../wrf_fortran.sh"
wrf_need r8d8
cd "${RUN_DIR}"

FORT_BIN="${RUN_DIR}/driver_radar_fortran"
CPP_BIN="${RUN_DIR}/driver_radar_cpu"
wrf_fortran "${FORT_BIN}" "${SCRIPT_DIR}/driver_radar.f90"
INPUT_BIN="${RUN_DIR}/radar_input.bin"
FORT_OUT="${RUN_DIR}/radar_output_fortran.bin"
CPP_OUT="${RUN_DIR}/radar_output_cpp.bin"

# hail_opt 0 and 1 (mp_wsm6_init's graupel / hail density); 41 levels: the scalar log10 of the odd last level
for cfg in "2500 50 0" "2500 50 1" "2499 41 0"; do
  set -- $cfg
  echo "=== WSM6 radar reflectivity: $1 columns, $2 levels, hail_opt $3 ==="
  python3 "${SCRIPT_DIR}/../mp_wsm6/generate_inputs.py" --cols "$1" --nz "$2" --dt 60.0 --hail "$3" --out "${INPUT_BIN}"
  "${HEAVY}" "${FORT_BIN}" "${INPUT_BIN}" "${FORT_OUT}"
  "${HEAVY}" "${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}"
  python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPP_OUT}" "$1" "$2"
done
echo "=== WSM6 radar reflectivity test PASSED ==="
