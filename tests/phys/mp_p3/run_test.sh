#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# P3 (mp_physics = 50 / 51) column test: WRF v4.8.0's p3_init and mp_p3_wrapper_wrf (module_mp_p3.o of the -r8
# -fdefault-double-8 build, WRF_R8D8_DIR) against the C++ p3_init / run_p3 (driver_p3_cpu of the CMake build), several
# calls in a row (TH_OLD / QV_OLD carried), tolerance 0 in every output.
#   run_test.sh BUILD_DIR      (BUILD_DIR: the directory of driver_p3_cpu; the work files go there)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="$(cd "${1:?usage: run_test.sh BUILD_DIR}" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
HEAVY="${REPO_ROOT}/tools/heavy.sh"
. "${SCRIPT_DIR}/../wrf_fortran.sh"
wrf_need r8d8
cd "${RUN_DIR}"
# p3_init reads lookupTable_1 from '.' (WRF's run directory)
TABLE="$W/run/p3_lookupTable_1.dat-v5.4_2momI"
[ -f "$TABLE" ] || { echo "no $TABLE"; exit 1; }
ln -sf "$TABLE" p3_lookupTable_1.dat-v5.4_2momI

FORT_BIN="${RUN_DIR}/driver_p3_fortran"
CPP_BIN="${RUN_DIR}/driver_p3_cpu"
wrf_fortran "${FORT_BIN}" "${SCRIPT_DIR}/driver_p3.f90"

INPUT_BIN="${RUN_DIR}/p3_input.bin"
FORT_OUT="${RUN_DIR}/p3_output_fortran.bin"
CPP_OUT="${RUN_DIR}/p3_output_cpp.bin"
# columns odd / even (the scalar tail of the vectorised exner pow), levels in the 48 / 64 / 128 tiers, dt 30 ... 300 s
# (sedimentation sub-steps), the cold start (itimestep 1: saturation adjustment, no specified-Nc activation) and a
# later step with th_old / qv_old given, specified (50) and predicted (51) droplet number
for cfg in "1001 41 60 1 0 7" "1001 41 60 1 1 8" "1000 50 120 5 0 9" "999 60 30 5 1 10" "601 90 300 2 1 11" \
           "600 44 180 3 0 12"; do
  set -- $cfg
  echo "=== P3 Differential Test: $1 columns, $2 levels, dt $3 s, first itimestep $4, predicted Nc $5 ==="
  python3 "${SCRIPT_DIR}/generate_inputs.py" --cols "$1" --nz "$2" --dt "$3" --steps 3 --start-it "$4" \
    --predict-nc "$5" --seed "$6" --out "${INPUT_BIN}"
  "${HEAVY}" "${FORT_BIN}" "${INPUT_BIN}" "${FORT_OUT}" > /dev/null
  "${HEAVY}" "${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}"
  python3 "${SCRIPT_DIR}/compare.py" "${INPUT_BIN}" "${FORT_OUT}" "${CPP_OUT}"
done
echo "=== P3 Differential Test PASSED ==="
