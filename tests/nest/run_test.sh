#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# ==============================================================================
# CPPWRF: differential (CPU) and CUDA test of the nest numerics: interpolation, boundary forcing, feedback and
# initialisation kernels against WRF v4.8.0's own compiled routines (the -r8 build, WRF_DIR; called through
# fortran_nest_bridge.f90, linked by tests/phys/wrf_fortran.sh) on the coarse fields of the R0 reference run
# (CPPWRF_WRF_RUNS). CPU FP64: tolerance 0; CUDA FP64 (with nvcc and a GPU): relative 1e-12.
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
. "${REPO_ROOT}/tests/phys/wrf_fortran.sh"
wrf_need r8
[ -n "${CPPWRF_WRF_RUNS:-}" ] || { echo "SKIP: set CPPWRF_WRF_RUNS (run_R0_r8)"; exit 77; }
OUT_DIR="${DIFF_TEST_BUILD_DIR:-${CPPWRF_SCRATCH:-${HOME}/cppwrf_scratch}}/nest_test"
HEAVY="${REPO_ROOT}/tools/heavy.sh"
mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

echo "=================================================================="
echo "CPPWRF nest numerics verification"
echo "=================================================================="

# 1. the coarse fields of the R0 reference run
DATASET="${OUT_DIR}/coarse_fields.bin"
python3 "${SCRIPT_DIR}/extract_coarse_fields.py" "${DATASET}"

# 2. CPU differential test: the C++ kernels and the bridge to WRF's routines in one program
echo "Compiling CPU differential test suite..."
CXXOPT=(-std=c++20 -DCPPWRF_DOUBLE -O3 -ffp-contract=off -fopenmp -I"${REPO_ROOT}/src" -I"${REPO_ROOT}")
"${HEAVY}" g++ "${CXXOPT[@]}" -c "${SCRIPT_DIR}/test_nest_differential.cpp" -o test_nest_differential.o
for s in dyn/lateral_bc core/halo core/backend; do
  "${HEAVY}" g++ "${CXXOPT[@]}" -c "${REPO_ROOT}/src/${s}.cpp" -o "$(basename "$s").o"
done
wrf_fortran "${OUT_DIR}/test_nest_differential" "${SCRIPT_DIR}/fortran_nest_bridge.f90" test_nest_differential.o \
  lateral_bc.o halo.o backend.o -fopenmp -lstdc++

echo "Executing CPU differential test suite..."
"${HEAVY}" "${OUT_DIR}/test_nest_differential" "${DATASET}"

# 3. CUDA FP64 device kernel test
if command -v nvcc &>/dev/null && nvidia-smi &>/dev/null; then
  echo "Compiling CUDA device objects..."
  NVOPT=(-std=c++20 -O3 -DCPPWRF_DOUBLE -DCPPWRF_USE_CUDA --expt-relaxed-constexpr --extended-lambda
         -I"${REPO_ROOT}/src" -I"${REPO_ROOT}")
  for s in dyn/lateral_bc core/halo core/backend; do
    "${HEAVY}" nvcc -x cu -c "${NVOPT[@]}" "${REPO_ROOT}/src/${s}.cpp" -o "$(basename "$s")_cuda.o"
  done
  "${HEAVY}" nvcc -c "${NVOPT[@]}" "${SCRIPT_DIR}/test_nest_cuda.cu" -o test_nest_cuda.o
  wrf_fortran "${OUT_DIR}/test_nest_cuda" "${SCRIPT_DIR}/fortran_nest_bridge.f90" test_nest_cuda.o \
    lateral_bc_cuda.o halo_cuda.o backend_cuda.o -lcudart -lstdc++
  echo "Executing CUDA FP64 device kernel test on GPU..."
  "${HEAVY}" "${OUT_DIR}/test_nest_cuda" "${DATASET}"
else
  echo "CUDA not detected; skipping GPU device test."
fi

echo "=================================================================="
echo "ALL CPU DIFFERENTIAL AND CUDA TESTS COMPLETED SUCCESSFULLY!"
echo "=================================================================="
