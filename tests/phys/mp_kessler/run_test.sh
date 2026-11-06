#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Kessler (mp_physics = 1) column test: WRF v4.8.0's kessler (module_mp_kessler.o of the -r8 build, WRF_DIR)
# against the C++ KesslerKernel (driver_kessler_cpu of the CMake build), tolerance 0.
#   run_test.sh BUILD_DIR      (BUILD_DIR: the directory of driver_kessler_cpu; the work files go there)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="$(cd "${1:?usage: run_test.sh BUILD_DIR}" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
HEAVY="${REPO_ROOT}/tools/heavy.sh"
. "${SCRIPT_DIR}/../wrf_fortran.sh"
wrf_need r8
cd "${RUN_DIR}"

FORT_BIN="${RUN_DIR}/driver_kessler_fortran"
CPP_BIN="${RUN_DIR}/driver_kessler_cpu"
wrf_fortran "${FORT_BIN}" "${SCRIPT_DIR}/driver_kessler.f90"

INPUT_BIN="${RUN_DIR}/kessler_input.bin"
FORT_OUT="${RUN_DIR}/kessler_output_fortran.bin"
CPP_OUT="${RUN_DIR}/kessler_output_cpp.bin"

echo "=== Running Kessler Microphysics Differential Test ==="
# dt 12 s: one sedimentation step (nfall = 1) in every column; dt 60 s (the model's step of the real cases):
# nfall 1-5, the split sedimentation with its vectorised fall-speed loop; odd column and level counts: the
# trailing scalar iteration of WRF's vectorised loops
for cfg in "2500 50 12.0" "2499 41 12.0" "2500 50 60.0" "2499 41 60.0"; do
  set -- $cfg
  echo "--- $1 columns, $2 levels, dt $3 s"
  python3 "${SCRIPT_DIR}/generate_inputs.py" --ncols "$1" --nz "$2" --dt "$3" -o "${INPUT_BIN}"
  "${HEAVY}" "${FORT_BIN}" "${INPUT_BIN}" "${FORT_OUT}"
  "${HEAVY}" "${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}"
  python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPP_OUT}" --tol 0
done
echo "=== Kessler Differential Test PASSED ==="
