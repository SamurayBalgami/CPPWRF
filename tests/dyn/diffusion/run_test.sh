#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Diffusion and damping kernels (km_opt 1-4, diff_6th, w_damp, damping layers): WRF v4.8.0's own module_diffusion_em,
# module_big_step_utilities_em, ... (the -r8 build, WRF_DIR; tests/phys/wrf_fortran.sh) called by driver_diffusion.f90
# against the C++ kernels (driver_diffusion.cpp), tolerance 0 in CPU FP64; with nvcc and a GPU also the CUDA FP64
# build (relative 1e-12).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
HEAVY="${REPO_ROOT}/tools/heavy.sh"
. "${REPO_ROOT}/tests/phys/wrf_fortran.sh"
wrf_need r8
BUILD_DIR="${DIFF_TEST_BUILD_DIR:-${SCRIPT_DIR}/build}/diffusion"
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

echo "== 1. Fortran driver (WRF objects of ${W})"
wrf_fortran "${BUILD_DIR}/driver_diffusion_f90" "${SCRIPT_DIR}/driver_diffusion.f90"

echo "== 2. C++ CPU FP64 driver"
"${HEAVY}" g++ -std=c++20 -DCPPWRF_DOUBLE -O3 -march=native -ffp-contract=off -fopenmp \
  -I "${REPO_ROOT}/src" "${SCRIPT_DIR}/driver_diffusion.cpp" "${REPO_ROOT}/src/core/backend.cpp" -lmvec \
  -o "${BUILD_DIR}/driver_diffusion_cpp"

INPUT_BIN="${BUILD_DIR}/diffusion_input.bin"
FORT_OUT="${BUILD_DIR}/diffusion_output_fortran.bin"
CPP_OUT="${BUILD_DIR}/diffusion_output_cpp.bin"
CUDA_OUT="${BUILD_DIR}/diffusion_output_cuda.bin"

echo "== 3. inputs (terrain, map factors, shear, steps)"
python3 "${SCRIPT_DIR}/generate_inputs.py" "${INPUT_BIN}"

echo "== 4. WRF and C++ CPU runs, comparison (tolerance 0)"
"${HEAVY}" "${BUILD_DIR}/driver_diffusion_f90" "${INPUT_BIN}" "${FORT_OUT}"
"${HEAVY}" "${BUILD_DIR}/driver_diffusion_cpp" "${INPUT_BIN}" "${CPP_OUT}"
python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPP_OUT}" --tol 0

if command -v nvcc > /dev/null && nvidia-smi > /dev/null 2>&1; then
  echo "== 5. C++ CUDA FP64 driver (relative 1e-12)"
  "${HEAVY}" nvcc -std=c++20 -DCPPWRF_DOUBLE -DCPPWRF_USE_CUDA -O3 --expt-relaxed-constexpr --extended-lambda \
    -I "${REPO_ROOT}/src" -x cu "${SCRIPT_DIR}/driver_diffusion_cuda.cu" -x cu "${REPO_ROOT}/src/core/backend.cpp" \
    -lmvec -o "${BUILD_DIR}/driver_diffusion_cuda"
  "${HEAVY}" "${BUILD_DIR}/driver_diffusion_cuda" "${INPUT_BIN}" "${CUDA_OUT}"
  python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CUDA_OUT}" --tol 1e-12
else
  echo "== 5. no nvcc / GPU: the CUDA build is not checked"
fi
