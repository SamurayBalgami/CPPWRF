#!/bin/bash
# ==============================================================================
# CPPWRF Automated Differential & CUDA Test Runner for Nest Numerics
# Tests all nesting interpolation, boundary forcing, feedback, and init kernels
# against unmodified WRF 4.6.1 Fortran (-fdefault-real-8)
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUT_DIR="${HOME}/cppwrf_scratch/nest_test"

export TMPDIR="${OUT_DIR}/tmp"
mkdir -p "${TMPDIR}"
mkdir -p "${OUT_DIR}"

echo "=================================================================="
echo "CPPWRF Nest Numerics Verification Test Suite"
echo "=================================================================="

# 1. Dataset Extraction
DATASET="${OUT_DIR}/coarse_fields.bin"
if [[ ! -f "${DATASET}" ]]; then
  echo "Extracting coarse dataset from WRF run..."
  "${REPO_ROOT}/tools/heavy.sh" --light python3 "${SCRIPT_DIR}/extract_coarse_fields.py" "${DATASET}"
fi

# 2. Compile Fortran Bridge against WRF 4.6.1 reference
WRF_ROOT="${HOME}/wrf_real/WRF-4.6.1-r8"
BRIDGE_OBJ="${OUT_DIR}/fortran_nest_bridge.o"
if [[ ! -f "${BRIDGE_OBJ}" || "${SCRIPT_DIR}/fortran_nest_bridge.f90" -nt "${BRIDGE_OBJ}" ]]; then
  echo "Compiling Fortran nest bridge against WRF 4.6.1 libraries..."
  "${REPO_ROOT}/tools/heavy.sh" gfortran -c -fdefault-real-8 -O3 -fPIC \
    -I"${WRF_ROOT}/inc" -I"${WRF_ROOT}/main" -I"${WRF_ROOT}/frame" \
    -I"${WRF_ROOT}/share" -I"${WRF_ROOT}/dyn_em" \
    -J"${OUT_DIR}" \
    "${SCRIPT_DIR}/fortran_nest_bridge.f90" -o "${BRIDGE_OBJ}"
fi

WRF_LIBS=(
  "${WRF_ROOT}/main/libwrflib.a"
  "${WRF_ROOT}/external/RSL_LITE/librsl_lite.a"
  "${WRF_ROOT}/frame/module_internal_header_util.o"
  "${WRF_ROOT}/frame/pack_utils.o"
  "${WRF_ROOT}/external/io_int/libwrfio_int.a"
  -L"${WRF_ROOT}/external/esmf_time_f90" -lesmf_time
  -L"${WRF_ROOT}/external/io_netcdf" -lwrfio_nf
  -L"${HOME}/wrf_real/netcdf/lib" -lnetcdff -lnetcdf
  "${WRF_ROOT}/external/fftpack/fftpack5/libfftpack.a"
  "${WRF_ROOT}/external/io_grib1/libio_grib1.a"
  "${WRF_ROOT}/external/io_grib_share/libio_grib_share.a"
  -lgfortran -lm -fopenmp
)

# 3. CPU Differential Test
echo "Compiling CPU differential test suite..."
"${REPO_ROOT}/tools/heavy.sh" g++ -std=c++20 -DCPPWRF_DOUBLE -O3 \
  -I"${REPO_ROOT}/src" -I"${REPO_ROOT}" \
  "${SCRIPT_DIR}/test_nest_differential.cpp" \
  "${REPO_ROOT}/src/dyn/lateral_bc.cpp" \
  "${REPO_ROOT}/src/core/halo.cpp" \
  "${REPO_ROOT}/src/core/backend.cpp" \
  "${BRIDGE_OBJ}" \
  "${WRF_LIBS[@]}" \
  -o "${OUT_DIR}/test_nest_differential"

echo "Executing CPU differential test suite..."
"${REPO_ROOT}/tools/heavy.sh" "${OUT_DIR}/test_nest_differential" "${DATASET}"

# 4. CUDA FP64 Device Kernel Test
if command -v nvcc &>/dev/null && nvidia-smi &>/dev/null; then
  echo "Compiling CUDA device objects..."
  "${REPO_ROOT}/tools/heavy.sh" nvcc -x cu -c -std=c++20 -O3 -DCPPWRF_DOUBLE -DCPPWRF_USE_CUDA \
    --expt-relaxed-constexpr --extended-lambda -I"${REPO_ROOT}/src" \
    "${REPO_ROOT}/src/dyn/lateral_bc.cpp" -o "${OUT_DIR}/lateral_bc_cuda.o"

  "${REPO_ROOT}/tools/heavy.sh" nvcc -x cu -c -std=c++20 -O3 -DCPPWRF_DOUBLE -DCPPWRF_USE_CUDA \
    --expt-relaxed-constexpr --extended-lambda -I"${REPO_ROOT}/src" \
    "${REPO_ROOT}/src/core/halo.cpp" -o "${OUT_DIR}/halo_cuda.o"

  "${REPO_ROOT}/tools/heavy.sh" nvcc -x cu -c -std=c++20 -O3 -DCPPWRF_DOUBLE -DCPPWRF_USE_CUDA \
    --expt-relaxed-constexpr --extended-lambda -I"${REPO_ROOT}/src" \
    "${REPO_ROOT}/src/core/backend.cpp" -o "${OUT_DIR}/backend_cuda.o"

  "${REPO_ROOT}/tools/heavy.sh" nvcc -c -std=c++20 -O3 -DCPPWRF_DOUBLE -DCPPWRF_USE_CUDA \
    --expt-relaxed-constexpr --extended-lambda -I"${REPO_ROOT}/src" -I"${REPO_ROOT}" \
    "${SCRIPT_DIR}/test_nest_cuda.cu" -o "${OUT_DIR}/test_nest_cuda.o"

  "${REPO_ROOT}/tools/heavy.sh" --light g++ -std=c++20 -DCPPWRF_DOUBLE \
    -I"${REPO_ROOT}/src" -I"${REPO_ROOT}" \
    "${OUT_DIR}/test_nest_cuda.o" \
    "${OUT_DIR}/lateral_bc_cuda.o" \
    "${OUT_DIR}/halo_cuda.o" \
    "${OUT_DIR}/backend_cuda.o" \
    "${BRIDGE_OBJ}" \
    "${WRF_LIBS[@]}" \
    -lcudart \
    -o "${OUT_DIR}/test_nest_cuda"

  echo "Executing CUDA FP64 device kernel test on GPU..."
  "${REPO_ROOT}/tools/heavy.sh" "${OUT_DIR}/test_nest_cuda" "${DATASET}"
else
  echo "CUDA not detected; skipping GPU device test."
fi

echo "=================================================================="
echo "ALL CPU DIFFERENTIAL AND CUDA TESTS COMPLETED SUCCESSFULLY!"
echo "=================================================================="
