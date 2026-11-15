#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# WSM5 (mp_physics = 4) and WSM3 (mp_physics = 3) column test: WRF v4.8.0's wsm5init / wsm5 and wsm3init / wsm3
# (module_mp_wsm5.o, module_mp_wsm3.o of the -r8 build, WRF_DIR) against the C++ Wsm5Kernel (+ WSM6's effective radii)
# and Wsm3Kernel (driver_wsm_cpu of the CMake build), three calls in a row, tolerance 0 in every output and in
# wsm5init's / wsm3init's constants.
#   run_test.sh BUILD_DIR      (BUILD_DIR: the directory of driver_wsm_cpu; the work files go there)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="$(cd "${1:?usage: run_test.sh BUILD_DIR}" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
HEAVY="${REPO_ROOT}/tools/heavy.sh"
. "${SCRIPT_DIR}/../wrf_fortran.sh"
wrf_need r8
cd "${RUN_DIR}"

FORT_BIN="${RUN_DIR}/driver_wsm_fortran"
CPP_BIN="${RUN_DIR}/driver_wsm_cpu"
wrf_fortran "${FORT_BIN}" "${SCRIPT_DIR}/driver_wsm5.f90"

INPUT_BIN="${RUN_DIR}/wsm_input.bin"
FORT_OUT="${RUN_DIR}/wsm_output_fortran.bin"
CPP_OUT="${RUN_DIR}/wsm_output_cpp.bin"
# odd and even ni (the scalar tail of the vectorised i loops; ni = 1: no vector lanes), levels in the 48 / 128 tiers,
# dt 60 / 120 s (one microphysics loop), 180 s (nint(1.5) = 2 loops), 300 s (3 loops)
for scheme in 5 3; do
  for cfg in "41 25 44 60 7" "40 20 50 180 8" "1 30 60 120 9" "33 9 90 300 10" "64 12 35 240 11"; do
    set -- $cfg
    echo "=== WSM${scheme} Differential Test: $1 x $2 columns, $3 levels, dt $4 s ==="
    python3 "${SCRIPT_DIR}/generate_inputs.py" --ni "$1" --nj "$2" --nz "$3" --dt "$4" --steps 3 --scheme "$scheme" \
      --seed "$5" --out "${INPUT_BIN}"
    "${HEAVY}" "${FORT_BIN}" "${INPUT_BIN}" "${FORT_OUT}"
    "${HEAVY}" "${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}"
    python3 "${SCRIPT_DIR}/compare.py" "${INPUT_BIN}" "${FORT_OUT}" "${CPP_OUT}"
  done
done
echo "=== WSM5 / WSM3 Differential Test PASSED ==="
