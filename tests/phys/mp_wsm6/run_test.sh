#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# WSM6 (mp_physics = 6) column test: WRF v4.8.0's wsm6 with the effective radii for RRTMG (module_mp_wsm6.o,
# mp_wsm6.o, mp_wsm6_effectRad.o of the -r8 build, WRF_DIR) against the C++ Wsm6Kernel + Wsm6EffRadKernel
# (driver_wsm6_cpu of the CMake build), tolerance 0.
#   run_test.sh BUILD_DIR      (BUILD_DIR: the directory of driver_wsm6_cpu; the work files go there)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="$(cd "${1:?usage: run_test.sh BUILD_DIR}" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
HEAVY="${REPO_ROOT}/tools/heavy.sh"
. "${SCRIPT_DIR}/../wrf_fortran.sh"
wrf_need r8
cd "${RUN_DIR}"

FORT_BIN="${RUN_DIR}/driver_wsm6_fortran"
CPP_BIN="${RUN_DIR}/driver_wsm6_cpu"
wrf_fortran "${FORT_BIN}" "${SCRIPT_DIR}/driver_wsm6.f90"

INPUT_BIN="${RUN_DIR}/wsm6_input.bin"
FORT_OUT="${RUN_DIR}/wsm6_output_fortran.bin"
CPP_OUT="${RUN_DIR}/wsm6_output_cpp.bin"

# dt 240 s: two microphysics loops (WRF's dtcldcr = 120 s), dt 120 / 60 s: one; hail_opt 1: the hail
# parameters of mp_wsm6_init; 2499 columns x 41 levels: the scalar tails of WRF's vectorised loops
for cfg in "2500 50 240.0 0" "2500 50 120.0 0" "2500 50 60.0 1" "2499 41 60.0 0"; do
  set -- $cfg
  echo "=== WSM6 Differential Test: $1 columns, $2 levels, dt $3 s, hail_opt $4 ==="
  python3 "${SCRIPT_DIR}/generate_inputs.py" --cols "$1" --nz "$2" --dt "$3" --hail "$4" --out "${INPUT_BIN}"
  "${HEAVY}" "${FORT_BIN}" "${INPUT_BIN}" "${FORT_OUT}"
  "${HEAVY}" "${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}"
  python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPP_OUT}" --tol 0
done
echo "=== WSM6 Differential Test PASSED ==="
