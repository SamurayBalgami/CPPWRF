#!/usr/bin/env bash
# ============================================================================
# CPPWRF Task T4.3: Lateral Boundary Differential Test Suite Runner
# Builds and executes the differential tests comparing C++20 vs WRF Fortran
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"

mkdir -p "${BUILD_DIR}"
cd "${REPO_ROOT}"

echo "================================================================="
echo "   Building T4.3 Lateral Boundary Differential Test Suite"
echo "================================================================="

# 1. Compile Fortran bridge and static archive if not already built or if updated
if [ ! -f "${BUILD_DIR}/libwrfbdy_ref.a" ] || [ "${BUILD_DIR}/bdy_fortran_bridge.f90" -nt "${BUILD_DIR}/libwrfbdy_ref.a" ]; then
    echo "[1/3] Compiling Fortran reference objects and bridge..."
    tools/heavy.sh gfortran -c -fdefault-real-8 -fdefault-double-8 -O2 \
        -I/home/ostriquetrum/wrf_real/WRF-4.6.1-r8/frame \
        -I/home/ostriquetrum/wrf_real/WRF-4.6.1-r8/share \
        -I/home/ostriquetrum/wrf_real/WRF-4.6.1-r8/dyn_em \
        -I/home/ostriquetrum/wrf_real/WRF-4.6.1-r8/main \
        -I/home/ostriquetrum/wrf_real/WRF-4.6.1-r8/inc \
        -I"${BUILD_DIR}" "${BUILD_DIR}/bdy_fortran_bridge.f90" \
        -o "${BUILD_DIR}/bdy_fortran_bridge.o" -J"${BUILD_DIR}"

    ar rcs "${BUILD_DIR}/libwrfbdy_ref.a" \
        "${BUILD_DIR}/module_bc.o" \
        "${BUILD_DIR}/module_bc_em.o" \
        "${BUILD_DIR}/wrf_stubs.o" \
        "${BUILD_DIR}/bdy_fortran_bridge.o"
fi

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
