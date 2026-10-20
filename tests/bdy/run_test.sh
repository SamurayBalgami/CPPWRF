#!/usr/bin/env bash
# ============================================================================
# Lateral boundary differential tests (C++ vs WRF)
# Builds and executes the differential tests comparing C++20 vs WRF Fortran
# ============================================================================
set -euo pipefail
[ -n "${WRF_DIR:-}" ] || { echo "SKIP: set WRF_DIR to a WRF v4.8.0 -r8 build"; exit 77; }  # ctest: skipped
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BUILD_DIR="${DIFF_TEST_BUILD_DIR:-${CPPWRF_SCRATCH:-${HOME}/cppwrf_scratch}/bdy_test}"

mkdir -p "${BUILD_DIR}"
cd "${REPO_ROOT}"

echo "================================================================="
echo "   Lateral boundary differential test (C++ vs WRF module_bc / module_bc_em)"
echo "================================================================="

# 1. Bridge (bdy_fortran_bridge.f90) + WRF's own module_bc.o and module_bc_em.o + error stubs
echo "[1/3] Compiling the Fortran bridge against ${WRF_DIR}..."
cp "${WRF_DIR}/share/module_bc.o" "${WRF_DIR}/dyn_em/module_bc_em.o" "${BUILD_DIR}/"
tools/heavy.sh gcc -c -O2 "${SCRIPT_DIR}/wrf_stubs.c" -o "${BUILD_DIR}/wrf_stubs.o"
tools/heavy.sh gfortran -c -fdefault-real-8 -fdefault-double-8 -O2 \
    -I${WRF_DIR}/frame \
    -I${WRF_DIR}/share \
    -I${WRF_DIR}/dyn_em \
    -I${WRF_DIR}/main \
    -I${WRF_DIR}/inc \
    "${SCRIPT_DIR}/bdy_fortran_bridge.f90" \
    -o "${BUILD_DIR}/bdy_fortran_bridge.o" -J"${BUILD_DIR}"
rm -f "${BUILD_DIR}/libwrfbdy_ref.a"
ar rcs "${BUILD_DIR}/libwrfbdy_ref.a" \
    "${BUILD_DIR}/module_bc.o" \
    "${BUILD_DIR}/module_bc_em.o" \
    "${BUILD_DIR}/wrf_stubs.o" \
    "${BUILD_DIR}/bdy_fortran_bridge.o"

# 2. Compile C++ differential test executable
echo "[2/3] Compiling C++ differential test executable..."
tools/heavy.sh g++ -std=c++20 -O2 \
    -I"${REPO_ROOT}/src" \
    -I/usr/include/hdf5/serial \
    "${SCRIPT_DIR}/test_bdy_differential.cpp" \
    "${BUILD_DIR}/libwrfbdy_ref.a" \
    -lgfortran -lm -lnetcdf \
    -o "${BUILD_DIR}/test_bdy_differential"

# 3. Run the test suite
echo "[3/3] Running lateral boundary differential test suite..."
tools/heavy.sh "${BUILD_DIR}/test_bdy_differential"
