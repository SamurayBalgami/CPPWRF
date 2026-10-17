#!/usr/bin/env bash
# Column differential test for MYJ PBL (bl_pbl_physics=2) and Eta Surface Layer (sf_sfclay_physics=2)
# Calls WRF v4.6.1's own compiled objects in ${WRF_DIR}/main/libwrflib.a.
# Pass criterion: bit-identical (tol = 0.0) in CPU FP64.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
HEAVY="${REPO_ROOT}/tools/heavy.sh"
WRF_DIR="${WRF_DIR:?set WRF_DIR to a WRF v4.6.1 -r8d8 build (main/libwrflib.a)}"
OUT_DIR="${1:-${PWD}/test_bl_myj}"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

echo "=== Building Fortran & C++ Differential Drivers (linking ${WRF_DIR}/main/libwrflib.a) ==="
"${HEAVY}" --light gfortran -O2 -ffree-form -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 \
  -I"${WRF_DIR}/main" -I"${WRF_DIR}/phys" -I"${WRF_DIR}/frame" -I"${WRF_DIR}/share" \
  "${SCRIPT_DIR}/driver_sf_myjsfc.f90" "${WRF_DIR}/main/libwrflib.a" -lmvec -lm -o driver_sf_myjsfc_fortran

"${HEAVY}" --light gfortran -O2 -ffree-form -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 \
  -I"${WRF_DIR}/main" -I"${WRF_DIR}/phys" -I"${WRF_DIR}/frame" -I"${WRF_DIR}/share" \
  "${SCRIPT_DIR}/driver_bl_myj.f90" "${WRF_DIR}/main/libwrflib.a" -lmvec -lm -o driver_bl_myj_fortran

"${HEAVY}" --light g++ -std=c++20 -DCPPWRF_DOUBLE -O3 -march=native -ffp-contract=off -fopenmp \
  -I"${REPO_ROOT}/src" \
  "${SCRIPT_DIR}/driver_sf_myjsfc.cpp" \
  "${REPO_ROOT}/src/phys/sf_myjsfc/myjsfc.cpp" \
  "${REPO_ROOT}/src/core/backend.cpp" \
  -lmvec -lm -o driver_sf_myjsfc_cpu_fp64

"${HEAVY}" --light g++ -std=c++20 -O3 -march=native -ffp-contract=off -fopenmp \
  -I"${REPO_ROOT}/src" \
  "${SCRIPT_DIR}/driver_sf_myjsfc.cpp" \
  "${REPO_ROOT}/src/phys/sf_myjsfc/myjsfc.cpp" \
  "${REPO_ROOT}/src/core/backend.cpp" \
  -lmvec -lm -o driver_sf_myjsfc_cpu_fp32

"${HEAVY}" --light g++ -std=c++20 -DCPPWRF_DOUBLE -O3 -march=native -ffp-contract=off -fopenmp \
  -I"${REPO_ROOT}/src" \
  "${SCRIPT_DIR}/driver_bl_myj.cpp" \
  "${REPO_ROOT}/src/phys/bl_myj/myjpbl.cpp" \
  "${REPO_ROOT}/src/core/backend.cpp" \
  -lmvec -lm -o driver_bl_myj_cpu_fp64

"${HEAVY}" --light g++ -std=c++20 -O3 -march=native -ffp-contract=off -fopenmp \
  -I"${REPO_ROOT}/src" \
  "${SCRIPT_DIR}/driver_bl_myj.cpp" \
  "${REPO_ROOT}/src/phys/bl_myj/myjpbl.cpp" \
  "${REPO_ROOT}/src/core/backend.cpp" \
  -lmvec -lm -o driver_bl_myj_cpu_fp32

echo "=== 1. Eta Surface Layer (MYJSFC) Differential Test (itimestep=2, 2500 cols) ==="
python3 "${SCRIPT_DIR}/generate_inputs_sfc.py" --cols 2500 --itimestep 2 --out "${OUT_DIR}/sf_myjsfc_in_t2.bin"
"${HEAVY}" --light ./driver_sf_myjsfc_fortran "${OUT_DIR}/sf_myjsfc_in_t2.bin" "${OUT_DIR}/sf_myjsfc_wrf_t2.bin"
"${HEAVY}" --light ./driver_sf_myjsfc_cpu_fp64 "${OUT_DIR}/sf_myjsfc_in_t2.bin" "${OUT_DIR}/sf_myjsfc_cpp_t2.bin"
python3 "${SCRIPT_DIR}/compare_sfc.py" "${OUT_DIR}/sf_myjsfc_wrf_t2.bin" "${OUT_DIR}/sf_myjsfc_cpp_t2.bin" --tol 0.0

echo "=== 2. Eta Surface Layer (MYJSFC) Differential Test (itimestep=1, 2500 cols) ==="
python3 "${SCRIPT_DIR}/generate_inputs_sfc.py" --cols 2500 --itimestep 1 --out "${OUT_DIR}/sf_myjsfc_in_t1.bin"
"${HEAVY}" --light ./driver_sf_myjsfc_fortran "${OUT_DIR}/sf_myjsfc_in_t1.bin" "${OUT_DIR}/sf_myjsfc_wrf_t1.bin"
"${HEAVY}" --light ./driver_sf_myjsfc_cpu_fp64 "${OUT_DIR}/sf_myjsfc_in_t1.bin" "${OUT_DIR}/sf_myjsfc_cpp_t1.bin"
python3 "${SCRIPT_DIR}/compare_sfc.py" "${OUT_DIR}/sf_myjsfc_wrf_t1.bin" "${OUT_DIR}/sf_myjsfc_cpp_t1.bin" --tol 0.0

echo "=== 3. MYJ PBL (MYJPBL) Differential Test (2500 cols, multiple of 4) ==="
python3 "${SCRIPT_DIR}/generate_inputs.py" --cols 2500 --out "${OUT_DIR}/bl_myj_in_2500.bin"
"${HEAVY}" --light ./driver_bl_myj_fortran "${OUT_DIR}/bl_myj_in_2500.bin" "${OUT_DIR}/bl_myj_wrf_2500.bin"
"${HEAVY}" --light ./driver_bl_myj_cpu_fp64 "${OUT_DIR}/bl_myj_in_2500.bin" "${OUT_DIR}/bl_myj_cpp_2500.bin"
python3 "${SCRIPT_DIR}/compare.py" "${OUT_DIR}/bl_myj_wrf_2500.bin" "${OUT_DIR}/bl_myj_cpp_2500.bin" --tol 0.0

echo "=== 4. MYJ PBL (MYJPBL) Differential Test (2503 cols, VF=4 + 3-col epilogue) ==="
python3 "${SCRIPT_DIR}/generate_inputs.py" --cols 2503 --out "${OUT_DIR}/bl_myj_in_2503.bin"
"${HEAVY}" --light ./driver_bl_myj_fortran "${OUT_DIR}/bl_myj_in_2503.bin" "${OUT_DIR}/bl_myj_wrf_2503.bin"
"${HEAVY}" --light ./driver_bl_myj_cpu_fp64 "${OUT_DIR}/bl_myj_in_2503.bin" "${OUT_DIR}/bl_myj_cpp_2503.bin"
python3 "${SCRIPT_DIR}/compare.py" "${OUT_DIR}/bl_myj_wrf_2503.bin" "${OUT_DIR}/bl_myj_cpp_2503.bin" --tol 0.0

echo "=== 5. Smoke-testing FP32 drivers ==="
"${HEAVY}" --light ./driver_sf_myjsfc_cpu_fp32 "${OUT_DIR}/sf_myjsfc_in_t2.bin" /dev/null
"${HEAVY}" --light ./driver_bl_myj_cpu_fp32 "${OUT_DIR}/bl_myj_in_2500.bin" /dev/null

echo "=== ALL MYJ + MYJSFC COLUMN DIFFERENTIAL TESTS PASSED (BIT-IDENTICAL IN CPU FP64) ==="
