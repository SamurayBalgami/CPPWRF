#!/usr/bin/env bash
# Differential test of src/dyn/scalars.cpp (ScalarTransport, one RK stage) against WRF v4.8.0:
# scalar_bridge.f90 calls rk_update_scalar_pd, rk_scalar_tend and rk_update_scalar as solve_em does;
# periodic, specified, open and symmetric boundaries; tolerance 0 on the CPU in FP64.
set -euo pipefail
[ -n "${WRF_DIR:-}" ] || { echo "SKIP: set WRF_DIR to a WRF v4.8.0 -r8 build"; exit 77; }  # ctest: skipped
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
BUILD_DIR="${DIFF_TEST_BUILD_DIR:-${SCRIPT_DIR}/build}"
HEAVY="${REPO_ROOT}/tools/heavy.sh"
mkdir -p "${BUILD_DIR}"

echo "[1/3] Fortran bridge against ${WRF_DIR} (gfortran -fdefault-real-8, configure.wrf's flags)"
"${HEAVY}" gfortran -c -fdefault-real-8 -O2 -ftree-vectorize -funroll-loops \
  -I"${WRF_DIR}/frame" -I"${WRF_DIR}/share" -I"${WRF_DIR}/dyn_em" -I"${WRF_DIR}/main" -I"${WRF_DIR}/inc" \
  "${SCRIPT_DIR}/scalar_bridge.f90" -J"${BUILD_DIR}" -o "${BUILD_DIR}/scalar_bridge.o"

echo "[2/3] C++ test (g++ -std=c++20 -DCPPWRF_DOUBLE)"
"${HEAVY}" g++ -std=c++20 -DCPPWRF_DOUBLE -O3 -march=native -ffp-contract=off -fopenmp \
  -I"${REPO_ROOT}/src" \
  "${SCRIPT_DIR}/test_scalars.cpp" \
  "${REPO_ROOT}/src/dyn/scalars.cpp" "${REPO_ROOT}/src/dyn/state.cpp" "${REPO_ROOT}/src/dyn/lateral_bc.cpp" \
  "${REPO_ROOT}/src/dyn/large_step.cpp" "${REPO_ROOT}/src/dyn/diffusion_driver.cpp" \
  "${REPO_ROOT}/src/phys/driver/phy_prep.cpp" \
  "${REPO_ROOT}/src/core/halo.cpp" "${REPO_ROOT}/src/core/backend.cpp" \
  "${BUILD_DIR}/scalar_bridge.o" \
  "${WRF_DIR}/main/libwrflib.a" \
  "${WRF_DIR}/external/io_netcdf/libwrfio_nf.a" \
  "${WRF_DIR}/external/esmf_time_f90/libesmf_time.a" \
  "${WRF_DIR}/external/fftpack/fftpack5/libfftpack.a" \
  "${WRF_DIR}/external/io_int/libwrfio_int.a" \
  "${WRF_DIR}/external/io_int/module_internal_header_util.o" \
  "${WRF_DIR}/external/io_grib1/libio_grib1.a" \
  "${WRF_DIR}/external/io_grib_share/libio_grib_share.a" \
  "${WRF_DIR}/frame/pack_utils.o" \
  "${WRF_DIR}/external/RSL_LITE/librsl_lite.a" \
  -lnetcdff -lnetcdf -lgfortran -lmvec \
  -o "${BUILD_DIR}/test_scalars"

echo "[3/3] Running"
"${HEAVY}" "${BUILD_DIR}/test_scalars"
