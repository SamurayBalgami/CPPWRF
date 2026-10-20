#!/usr/bin/env bash
# Column differential test for the Noah ice surfaces: SFLX_GLACIAL (land ice) and SFLX_SEAICE
# (sea ice), called from WRF v4.8.0's own compiled objects ${WRF_DIR}/phys/module_sf_noahlsm.o,
# module_sf_noahlsm_glacial_only.o and module_sf_noah_seaice.o (their only outside reference,
# wrf_error_fatal3, is in stubs.f90).
# Pass criterion: bit-identical (tol = 0.0) in CPU FP64.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
HEAVY="${REPO_ROOT}/tools/heavy.sh"
WRF_DIR="${WRF_DIR:?set WRF_DIR to a WRF v4.8.0 -r8d8 build (phys/*.o)}"
OUT_DIR="${1:-${PWD}/test_lsm_noah_ice}"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

echo "=== Building Fortran & C++ Differential Drivers (linking WRF objects of ${WRF_DIR}/phys) ==="
"${HEAVY}" --light gfortran -O2 -ffree-form -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 \
  -I"${WRF_DIR}/phys" -J. \
  "${SCRIPT_DIR}/driver_noah_ice.f90" "${SCRIPT_DIR}/stubs.f90" "${WRF_DIR}/phys/module_sf_noahlsm.o" \
  "${WRF_DIR}/phys/module_sf_noahlsm_glacial_only.o" "${WRF_DIR}/phys/module_sf_noah_seaice.o" \
  -lmvec -lm -o driver_noah_ice_fortran

"${HEAVY}" --light g++ -std=c++20 -DCPPWRF_DOUBLE -O3 -march=native -ffp-contract=off -fopenmp \
  -I"${REPO_ROOT}/src" \
  "${SCRIPT_DIR}/driver_noah_ice.cpp" \
  "${REPO_ROOT}/src/phys/lsm_noah/noah.cpp" \
  "${REPO_ROOT}/src/phys/lsm_noah/noah_tables.cpp" \
  "${REPO_ROOT}/src/core/backend.cpp" \
  -lmvec -lm -o driver_noah_ice_cpu_fp64

for dt in 60 300 3600; do
  echo "=== Noah ice (SFLX_GLACIAL, SFLX_SEAICE) Differential Test (2500 cols, dt = ${dt} s) ==="
  python3 "${SCRIPT_DIR}/generate_inputs.py" --cols 2500 --dt "${dt}" --out "${OUT_DIR}/noah_ice_in_${dt}.bin"
  "${HEAVY}" --light ./driver_noah_ice_fortran "${OUT_DIR}/noah_ice_in_${dt}.bin" "${OUT_DIR}/noah_ice_wrf_${dt}.bin"
  "${HEAVY}" --light ./driver_noah_ice_cpu_fp64 "${OUT_DIR}/noah_ice_in_${dt}.bin" "${OUT_DIR}/noah_ice_cpp_${dt}.bin"
  python3 "${SCRIPT_DIR}/compare.py" "${OUT_DIR}/noah_ice_wrf_${dt}.bin" "${OUT_DIR}/noah_ice_cpp_${dt}.bin" --tol 0.0
done
echo "=== Noah ice Differential Test PASSED ==="
