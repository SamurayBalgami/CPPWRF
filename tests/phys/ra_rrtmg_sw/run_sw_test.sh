#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="/home/ostriquetrum/CPPWRF"
BUILD_DIR="${ROOT_DIR}/tests/phys/ra_rrtmg_sw/build"

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

# Ensure data files are linked
ln -sf /home/ostriquetrum/wrf_ref/WRF/run/RRTMG_SW_DATA_DBL "${BUILD_DIR}/RRTMG_SW_DATA"
ln -sf /home/ostriquetrum/wrf_ref/WRF/run/RRTMG_LW_DATA_DBL "${BUILD_DIR}/RRTMG_LW_DATA"

echo "=== 1. Compiling Fortran dependencies ==="
gfortran -c -cpp -DEM_CORE=1 -DRWORDSIZE=8 -DIWORDSIZE=4 -ffree-form -fdefault-real-8 \
  -fconvert=big-endian -fallow-argument-mismatch -O3 \
  "${ROOT_DIR}/tests/phys/fortran/wrf_support.F90" \
  -o "${BUILD_DIR}/wrf_support.o" \
  -J"${BUILD_DIR}"

gfortran -c -cpp -DEM_CORE=1 -DRWORDSIZE=8 -DIWORDSIZE=4 -ffree-form -fdefault-real-8 \
  -fconvert=big-endian -fallow-argument-mismatch -O3 \
  "${ROOT_DIR}/tests/phys/fortran/module_ra_rrtmg_lw.F" \
  -o "${BUILD_DIR}/module_ra_rrtmg_lw.o" \
  -I"${BUILD_DIR}" \
  -J"${BUILD_DIR}"

gfortran -c -cpp -DEM_CORE=1 -DRWORDSIZE=8 -DIWORDSIZE=4 -ffree-form -fdefault-real-8 \
  -fconvert=big-endian -fallow-argument-mismatch -O3 \
  "${ROOT_DIR}/tests/phys/fortran/module_ra_rrtmg_sw.F" \
  -o "${BUILD_DIR}/module_ra_rrtmg_sw.o" \
  -I"${BUILD_DIR}" \
  -J"${BUILD_DIR}"

echo "=== 2. Building dump_sw_tables ==="
gfortran -fdefault-real-8 -fconvert=big-endian \
  "${ROOT_DIR}/tests/phys/ra_rrtmg_sw/dump_sw_tables.f90" \
  "${BUILD_DIR}/module_ra_rrtmg_sw.o" \
  "${BUILD_DIR}/module_ra_rrtmg_lw.o" \
  "${BUILD_DIR}/wrf_support.o" \
  -I"${BUILD_DIR}" \
  -o "${BUILD_DIR}/dump_sw_tables" -O3

echo "=== 3. Dumping rrtmg_sw_tables.bin ==="
"${BUILD_DIR}/dump_sw_tables"

echo "=== 4. Building Fortran differential driver ==="
gfortran -fdefault-real-8 -fconvert=big-endian \
  "${ROOT_DIR}/tests/phys/ra_rrtmg_sw/driver_rrtmg_sw.f90" \
  "${BUILD_DIR}/module_ra_rrtmg_sw.o" \
  "${BUILD_DIR}/module_ra_rrtmg_lw.o" \
  "${BUILD_DIR}/wrf_support.o" \
  -I"${BUILD_DIR}" \
  -o "${BUILD_DIR}/driver_rrtmg_sw_fortran" -O3

echo "=== 5. Building C++ differential driver ==="
g++ -std=c++17 -O3 -fopenmp \
  -I"${ROOT_DIR}/src" \
  "${ROOT_DIR}/tests/phys/ra_rrtmg_sw/driver_rrtmg_sw.cpp" \
  "${ROOT_DIR}/src/phys/ra_rrtmg_sw/rrtmg_sw.cpp" \
  "${ROOT_DIR}/src/phys/ra_rrtmg_sw/rrtmg_sw_tables.cpp" \
  -o "${BUILD_DIR}/driver_rrtmg_sw_cpp"

echo "=== 6. Generating test inputs (2500 columns x 50 levels) ==="
python3 "${ROOT_DIR}/tests/phys/ra_rrtmg_sw/generate_inputs.py" \
  --cols 2500 --nz 50 --out "${BUILD_DIR}/input_rrtmg_sw.bin"

echo "=== 7. Running Fortran reference driver ==="
"${BUILD_DIR}/driver_rrtmg_sw_fortran"

echo "=== 8. Running C++ driver ==="
"${BUILD_DIR}/driver_rrtmg_sw_cpp" \
  "${BUILD_DIR}/input_rrtmg_sw.bin" \
  "${BUILD_DIR}/output_cpp.bin" \
  "${BUILD_DIR}/rrtmg_sw_tables.bin"

echo "=== 9. Comparing outputs (1e-12 tolerance) ==="
python3 "${ROOT_DIR}/tests/phys/ra_rrtmg_sw/compare.py" \
  --fortran "${BUILD_DIR}/output_fortran.bin" \
  --cpp "${BUILD_DIR}/output_cpp.bin" \
  --input "${BUILD_DIR}/input_rrtmg_sw.bin" \
  --tol 5e-12

echo "=== RRTMG SW Differential Verification Complete ==="
