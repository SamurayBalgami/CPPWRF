#!/usr/bin/env bash
# Differential test of src/dyn/lateral_bc.cpp (specified lateral boundaries) against WRF v4.8.0:
# lbc_bridge.f90 calls WRF's own relax_bdy_dry, spec_bdy_dry, spec_bdyupdate(_ph), zero_grad_bdy,
# relax/spec_bdy_scalar, flow_dep_bdy and spec_bdy_final as solve_em does; tolerance 0 on the CPU in FP64.
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
  "${SCRIPT_DIR}/lbc_bridge.f90" -J"${BUILD_DIR}" -o "${BUILD_DIR}/lbc_bridge.o"

echo "[2/3] C++ test (g++ -std=c++20 -DCPPWRF_DOUBLE)"
"${HEAVY}" g++ -std=c++20 -DCPPWRF_DOUBLE -O3 -march=native -ffp-contract=off -fopenmp \
  -I"${REPO_ROOT}/src" \
  "${SCRIPT_DIR}/test_lateral_bc.cpp" \
  "${REPO_ROOT}/src/dyn/lateral_bc.cpp" "${REPO_ROOT}/src/dyn/state.cpp" \
  "${REPO_ROOT}/src/core/halo.cpp" "${REPO_ROOT}/src/core/backend.cpp" \
  "${BUILD_DIR}/lbc_bridge.o" \
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
  -o "${BUILD_DIR}/test_lateral_bc"

echo "[3/3] Running"
"${HEAVY}" "${BUILD_DIR}/test_lateral_bc"
