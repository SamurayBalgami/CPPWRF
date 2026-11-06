#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Column differential test for Noah-MP v5.1.2 (sf_surface_physics=4).
# Calls WRF v4.8.0's compiled NoahmpWRFinit and NoahmpWRFmain of the -r8 build that made the NMP1 references
# (WRF_DIR; tests/phys/wrf_fortran.sh) with WRF's NoahmpTable.TBL.
# Pass criterion: bit-identical (tol = 0.0) in CPU FP64 across all 150 state & diagnostic outputs.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
HEAVY="${REPO_ROOT}/tools/heavy.sh"
. "${SCRIPT_DIR}/../wrf_fortran.sh"
wrf_need r8
OUT_DIR="${1:-${PWD}/test_lsm_noahmp}"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

# NoahmpReadTable in WRF's libwrflib.a opens "NoahmpTable.TBL" in the current working directory
cp -L "${W}/run/NoahmpTable.TBL" "${OUT_DIR}/NoahmpTable.TBL"

echo "=== Building Fortran & C++ Noah-MP Differential Drivers (linking ${W}/main/libwrflib.a) ==="
wrf_fortran "${OUT_DIR}/driver_noahmp_fortran" "${SCRIPT_DIR}/driver_noahmp.f90"

"${HEAVY}" --light g++ -std=c++20 -DCPPWRF_DOUBLE -O3 -march=native -ffp-contract=off -fopenmp \
  -I"${REPO_ROOT}/src" \
  "${SCRIPT_DIR}/driver_noahmp.cpp" \
  "${REPO_ROOT}/src/phys/lsm_noahmp/noahmp.cpp" \
  "${REPO_ROOT}/src/phys/lsm_noahmp/noahmp_tables.cpp" \
  "${REPO_ROOT}/src/core/backend.cpp" \
  -lmvec -lm -o driver_noahmp_cpu_fp64

"${HEAVY}" --light g++ -std=c++20 -O3 -march=native -ffp-contract=off -fopenmp \
  -I"${REPO_ROOT}/src" \
  "${SCRIPT_DIR}/driver_noahmp.cpp" \
  "${REPO_ROOT}/src/phys/lsm_noahmp/noahmp.cpp" \
  "${REPO_ROOT}/src/phys/lsm_noahmp/noahmp_tables.cpp" \
  "${REPO_ROOT}/src/core/backend.cpp" \
  -lmvec -lm -o driver_noahmp_cpu_fp32

echo "=== 1. Noah-MP Initialization + 1-Step Differential Test (2500 cols) ==="
python3 "${SCRIPT_DIR}/generate_inputs.py" --cols 2500 --steps 1 --out "${OUT_DIR}/noahmp_in_1step.bin"
"${HEAVY}" --light ./driver_noahmp_fortran "${OUT_DIR}/noahmp_in_1step.bin" "${OUT_DIR}/noahmp_wrf_1step.bin"
"${HEAVY}" --light ./driver_noahmp_cpu_fp64 "${OUT_DIR}/noahmp_in_1step.bin" "${OUT_DIR}/noahmp_cpp_1step.bin" "${OUT_DIR}/NoahmpTable.TBL"
python3 "${SCRIPT_DIR}/compare.py" "${OUT_DIR}/noahmp_wrf_1step.bin" "${OUT_DIR}/noahmp_cpp_1step.bin" --tol 0.0

echo "=== 2. Noah-MP Initialization + 5-Step Multi-Step Evolution Test (2500 cols) ==="
python3 "${SCRIPT_DIR}/generate_inputs.py" --cols 2500 --steps 5 --out "${OUT_DIR}/noahmp_in_5step.bin"
"${HEAVY}" --light ./driver_noahmp_fortran "${OUT_DIR}/noahmp_in_5step.bin" "${OUT_DIR}/noahmp_wrf_5step.bin"
"${HEAVY}" --light ./driver_noahmp_cpu_fp64 "${OUT_DIR}/noahmp_in_5step.bin" "${OUT_DIR}/noahmp_cpp_5step.bin" "${OUT_DIR}/NoahmpTable.TBL"
python3 "${SCRIPT_DIR}/compare.py" "${OUT_DIR}/noahmp_wrf_5step.bin" "${OUT_DIR}/noahmp_cpp_5step.bin" --tol 0.0

echo "=== 3. Smoke-testing FP32 driver ==="
"${HEAVY}" --light ./driver_noahmp_cpu_fp32 "${OUT_DIR}/noahmp_in_1step.bin" /dev/null "${OUT_DIR}/NoahmpTable.TBL"

echo "=== ALL NOAH-MP COLUMN DIFFERENTIAL TESTS PASSED (BIT-IDENTICAL IN CPU FP64) ==="
