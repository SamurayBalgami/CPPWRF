#!/usr/bin/env bash
# Column differential test for Noah-MP v5.1.2 (sf_surface_physics=4).
# Calls WRF v4.8.0's compiled NoahmpWRFinit and NoahmpWRFmain in ${WRF_DIR}/main/libwrflib.a.
# Pass criterion: bit-identical (tol = 0.0) in CPU FP64 across all 150 state & diagnostic outputs.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
HEAVY="${REPO_ROOT}/tools/heavy.sh"
[ -n "${WRF_DIR:-}" ] || { echo "SKIP: set WRF_DIR to a WRF v4.8.0 -r8 or -r8d8 build (main/libwrflib.a)"; exit 77; }  # ctest: skipped
OUT_DIR="${1:-${PWD}/test_lsm_noahmp}"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

# NoahmpReadTable in WRF's libwrflib.a opens "NoahmpTable.TBL" in the current working directory
TBL_SRC="${WRF_DIR}/run/NoahmpTable.TBL"
if [[ ! -f "${TBL_SRC}" ]]; then
  TBL_SRC="${WRF_DIR}/phys/noahmp/parameters/NoahmpTable.TBL"
fi
ln -sf "${TBL_SRC}" "${OUT_DIR}/NoahmpTable.TBL"

echo "=== Building Fortran & C++ Noah-MP Differential Drivers (linking ${WRF_DIR}/main/libwrflib.a) ==="
NETCDF_LIB="${WRF_DIR}/../netcdf/lib"
"${HEAVY}" --light gfortran -O2 -ffree-form -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 \
  -I"${WRF_DIR}/main" -I"${WRF_DIR}/phys" -I"${WRF_DIR}/frame" -I"${WRF_DIR}/share" \
  "${SCRIPT_DIR}/driver_noahmp.f90" \
  "${WRF_DIR}/main/libwrflib.a" \
  "${WRF_DIR}/external/fftpack/fftpack5/libfftpack.a" \
  "${WRF_DIR}/external/io_grib1/libio_grib1.a" \
  "${WRF_DIR}/external/io_grib_share/libio_grib_share.a" \
  "${WRF_DIR}/external/io_int/libwrfio_int.a" \
  "${WRF_DIR}/external/esmf_time_f90/libesmf_time.a" \
  "${WRF_DIR}/external/RSL_LITE/librsl_lite.a" \
  "${WRF_DIR}/frame/module_internal_header_util.o" \
  "${WRF_DIR}/frame/pack_utils.o" \
  "${WRF_DIR}/external/io_netcdf/libwrfio_nf.a" \
  -L"${NETCDF_LIB}" -Wl,-rpath,"${NETCDF_LIB}" -lnetcdff -lnetcdf -lmvec -lm -o driver_noahmp_fortran

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
