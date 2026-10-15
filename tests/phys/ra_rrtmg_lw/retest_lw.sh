#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
BUILD_DIR="${ROOT_DIR}/tests/phys/ra_rrtmg_lw/build"

cd "${ROOT_DIR}"

echo "=== 1. Compiling wrf_support and module_ra_rrtmg_lw.F with -DEM_CORE=1 ==="
gfortran -c -cpp -DEM_CORE=1 -DRWORDSIZE=8 -DIWORDSIZE=4 -ffree-form -fdefault-real-8 \
  -fconvert=big-endian -fallow-argument-mismatch -O3 \
  tests/phys/fortran/wrf_support.F90 \
  -o "${BUILD_DIR}/wrf_support.o" \
  -J"${BUILD_DIR}"

gfortran -c -cpp -DEM_CORE=1 -DRWORDSIZE=8 -DIWORDSIZE=4 -ffree-form -fdefault-real-8 \
  -fconvert=big-endian -fallow-argument-mismatch -O3 \
  -I"${BUILD_DIR}" \
  tests/phys/fortran/module_ra_rrtmg_lw.F \
  -o "${BUILD_DIR}/module_ra_rrtmg_lw.o" \
  -J"${BUILD_DIR}"

echo "=== 2. Linking driver_rrtmg_lw_fortran ==="
gfortran -fdefault-real-8 -fconvert=big-endian \
  tests/phys/ra_rrtmg_lw/driver_rrtmg_lw.f90 \
  "${BUILD_DIR}/module_ra_rrtmg_lw.o" \
  "${BUILD_DIR}/wrf_support.o" \
  -I"${BUILD_DIR}" \
  -o "${BUILD_DIR}/driver_rrtmg_lw_fortran" -O3

echo "=== 3. Compiling driver_rrtmg_lw_cpp ==="
g++ -std=c++17 -O3 -fopenmp \
  -I"${ROOT_DIR}/src" \
  "${ROOT_DIR}/tests/phys/ra_rrtmg_lw/driver_rrtmg_lw.cpp" \
  "${ROOT_DIR}/src/phys/ra_rrtmg_lw/rrtmg_lw.cpp" \
  "${ROOT_DIR}/src/phys/ra_rrtmg_lw/rrtmg_lw_tables.cpp" \
  -o "${BUILD_DIR}/driver_rrtmg_lw_cpp"

echo "=== 4. Running Fortran driver on input_rrtmg_lw.bin ==="
cd "${BUILD_DIR}"
./driver_rrtmg_lw_fortran input_rrtmg_lw.bin output_fortran.bin
cd "${ROOT_DIR}"

echo "=== 5. Running C++ driver on input_rrtmg_lw.bin ==="
"${BUILD_DIR}/driver_rrtmg_lw_cpp" "${BUILD_DIR}/input_rrtmg_lw.bin" "${BUILD_DIR}/output_cpp.bin" "${BUILD_DIR}/rrtmg_lw_tables.bin"

echo "=== 6. Differential Verification ==="
python3 tests/phys/ra_rrtmg_lw/compare.py "${BUILD_DIR}/output_fortran.bin" "${BUILD_DIR}/output_cpp.bin" 1e-6
