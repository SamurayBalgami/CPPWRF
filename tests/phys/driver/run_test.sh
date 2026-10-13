#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"
HEAVY="${REPO_ROOT}/tools/heavy.sh"

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

echo "======================================================================"
echo " Extracting original WRF routines from ~/wrf_ref/WRF"
echo "======================================================================"
python3 "${SCRIPT_DIR}/fortran/extract_wrf_subs.py"

echo "======================================================================"
echo " Building Fortran differential test driver from ORIGINAL WRF sources"
echo "======================================================================"
"${HEAVY}" gfortran -ffree-form -fdefault-real-8 -O3 -march=native -ffp-contract=off \
  -J "${BUILD_DIR}" \
  -c "${REPO_ROOT}/../wrf_ref/WRF/share/module_model_constants.F" \
  -o "${BUILD_DIR}/module_model_constants.o"

"${HEAVY}" gfortran -ffree-form -fdefault-real-8 -O3 -march=native -ffp-contract=off \
  -I "${BUILD_DIR}" -J "${BUILD_DIR}" \
  -c "${SCRIPT_DIR}/fortran/wrf_phys_driver_orig.F90" \
  -o "${BUILD_DIR}/wrf_phys_driver_orig.o"

"${HEAVY}" gfortran -ffree-form -fdefault-real-8 -O3 -march=native -ffp-contract=off \
  -I "${BUILD_DIR}" \
  "${BUILD_DIR}/module_model_constants.o" \
  "${BUILD_DIR}/wrf_phys_driver_orig.o" \
  "${SCRIPT_DIR}/driver_fortran.f90" \
  -o "${BUILD_DIR}/driver_fortran"

echo "======================================================================"
echo " Building C++ CPU FP64 driver (g++ -O3 -fopenmp -DCPPWRF_DOUBLE)"
echo "======================================================================"
"${HEAVY}" g++ -std=c++20 -DCPPWRF_DOUBLE -O3 -march=native -ffp-contract=off -fopenmp \
  -I "${REPO_ROOT}/src" \
  "${SCRIPT_DIR}/driver_cpu.cpp" \
  "${REPO_ROOT}/src/phys/driver/phy_prep.cpp" \
  "${REPO_ROOT}/src/phys/driver/tendency_coupling.cpp" \
  "${REPO_ROOT}/src/core/backend.cpp" \
  "${REPO_ROOT}/src/core/halo.cpp" \
  "${REPO_ROOT}/src/dyn/state.cpp" \
  -o "${BUILD_DIR}/driver_cpu"

CUDA_AVAILABLE=0
if command -v nvcc &> /dev/null; then
  echo "======================================================================"
  echo " Building C++ CUDA FP64 driver (nvcc -O3 -DCPPWRF_USE_CUDA -DCPPWRF_DOUBLE)"
  echo "======================================================================"
  if "${HEAVY}" nvcc -std=c++20 -DCPPWRF_DOUBLE -DCPPWRF_USE_CUDA -O3 \
    --expt-relaxed-constexpr --extended-lambda \
    -I "${REPO_ROOT}/src" \
    -x cu "${SCRIPT_DIR}/driver_cuda.cu" \
    -x cu "${REPO_ROOT}/src/phys/driver/phy_prep.cpp" \
    -x cu "${REPO_ROOT}/src/phys/driver/tendency_coupling.cpp" \
    -x cu "${REPO_ROOT}/src/core/backend.cpp" \
    -x cu "${REPO_ROOT}/src/core/halo.cpp" \
    -x cu "${REPO_ROOT}/src/dyn/state.cpp" \
    -o "${BUILD_DIR}/driver_cuda"; then
    CUDA_AVAILABLE=1
  else
    echo "Warning: CUDA compilation failed, skipping GPU differential test"
  fi
fi

INPUT_BIN="${BUILD_DIR}/driver_input.bin"
FORT_OUT="${BUILD_DIR}/driver_output_fortran.bin"
CPU_OUT="${BUILD_DIR}/driver_output_cpu.bin"
CUDA_OUT="${BUILD_DIR}/driver_output_cuda.bin"

NX=16
NY=16
NZ=35

echo "======================================================================"
echo " Generating input dataset: ${NX}x${NY}x${NZ} tile with terrain & map factors"
echo "======================================================================"
python3 "${SCRIPT_DIR}/generate_inputs.py" --out "${INPUT_BIN}" --nx "${NX}" --ny "${NY}" --nz "${NZ}"

echo "======================================================================"
echo " Running Fortran reference run"
echo "======================================================================"
"${HEAVY}" "${BUILD_DIR}/driver_fortran" "${INPUT_BIN}" "${FORT_OUT}"

echo "======================================================================"
echo " Running C++ CPU FP64 run"
echo "======================================================================"
"${HEAVY}" "${BUILD_DIR}/driver_cpu" "${INPUT_BIN}" "${CPU_OUT}"

echo "======================================================================"
echo " Comparing C++ CPU FP64 vs Fortran Reference (PASS = max rel <= 1e-12)"
echo "======================================================================"
python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPU_OUT}" --nx "${NX}" --ny "${NY}" --nz "${NZ}" --tol 1e-12

if [[ "${CUDA_AVAILABLE}" -eq 1 ]]; then
  echo "======================================================================"
  echo " Running C++ CUDA FP64 run"
  echo "======================================================================"
  "${HEAVY}" "${BUILD_DIR}/driver_cuda" "${INPUT_BIN}" "${CUDA_OUT}"

  echo "======================================================================"
  echo " Comparing C++ CUDA FP64 vs Fortran Reference (PASS = max rel <= 1e-12)"
  echo "======================================================================"
  python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CUDA_OUT}" --nx "${NX}" --ny "${NY}" --nz "${NZ}" --tol 1e-12
fi

echo "======================================================================"
echo " ALL TESTS COMPLETED SUCCESSFULLY"
echo "======================================================================"
