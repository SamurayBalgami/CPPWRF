#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"
HEAVY="${REPO_ROOT}/tools/heavy.sh"
WRF_R8="/home/ostriquetrum/wrf_real/WRF-4.6.1-r8"

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

echo "======================================================================"
echo " 1. Building Fortran Differential Driver (gfortran -fdefault-real-8)"
echo "======================================================================"
"${HEAVY}" gfortran -std=gnu -ffree-form -fdefault-real-8 -O3 -ffp-contract=off \
  -I "${WRF_R8}/dyn_em" \
  -I "${WRF_R8}/frame" \
  -I "${WRF_R8}/share" \
  "${SCRIPT_DIR}/driver_diffusion.f90" \
  "${WRF_R8}/main/libwrflib.a" \
  "${WRF_R8}/external/io_netcdf/libwrfio_nf.a" \
  "${WRF_R8}/external/esmf_time_f90/libesmf_time.a" \
  "${WRF_R8}/external/fftpack/fftpack5/libfftpack.a" \
  "${WRF_R8}/external/io_int/libwrfio_int.a" \
  "${WRF_R8}/external/io_int/module_internal_header_util.o" \
  "${WRF_R8}/external/io_grib1/libio_grib1.a" \
  "${WRF_R8}/external/io_grib_share/libio_grib_share.a" \
  "${WRF_R8}/frame/pack_utils.o" \
  "${WRF_R8}/external/RSL_LITE/librsl_lite.a" \
  -lnetcdff -lnetcdf \
  -o "${BUILD_DIR}/driver_diffusion_f90"

echo "======================================================================"
echo " 2. Building C++ CPU Differential Driver (g++ -std=c++20 -DCPPWRF_DOUBLE)"
echo "======================================================================"
"${HEAVY}" g++ -std=c++20 -DCPPWRF_DOUBLE -O3 -march=native -ffp-contract=off -fopenmp \
  -I "${REPO_ROOT}/src" \
  "${SCRIPT_DIR}/driver_diffusion.cpp" \
  "${REPO_ROOT}/src/core/backend.cpp" \
  -o "${BUILD_DIR}/driver_diffusion_cpp"

echo "======================================================================"
echo " 3. Building C++ CUDA FP64 Driver (nvcc -std=c++20 -DCPPWRF_USE_CUDA -DCPPWRF_DOUBLE)"
echo "======================================================================"
"${HEAVY}" nvcc -std=c++20 -DCPPWRF_DOUBLE -DCPPWRF_USE_CUDA -O3 \
  --expt-relaxed-constexpr --extended-lambda \
  -I "${REPO_ROOT}/src" \
  -x cu "${SCRIPT_DIR}/driver_diffusion_cuda.cu" \
  -x cu "${REPO_ROOT}/src/core/backend.cpp" \
  -o "${BUILD_DIR}/driver_diffusion_cuda"

echo "======================================================================"
echo " 4. Generating Physical 3D Test Inputs (Terrain, Map Factors, Shear, Steps)"
echo "======================================================================"
INPUT_BIN="${BUILD_DIR}/diffusion_input.bin"
FORT_OUT="${BUILD_DIR}/diffusion_output_fortran.bin"
CPP_OUT="${BUILD_DIR}/diffusion_output_cpp.bin"
CUDA_OUT="${BUILD_DIR}/diffusion_output_cuda.bin"

python3 "${SCRIPT_DIR}/generate_inputs.py" "${INPUT_BIN}"

echo "======================================================================"
echo " 5. Running Fortran Reference Driver"
echo "======================================================================"
"${HEAVY}" "${BUILD_DIR}/driver_diffusion_f90" "${INPUT_BIN}" "${FORT_OUT}"

echo "======================================================================"
echo " 6. Running CPPWRF C++ CPU Driver"
echo "======================================================================"
"${HEAVY}" "${BUILD_DIR}/driver_diffusion_cpp" "${INPUT_BIN}" "${CPP_OUT}"

echo "======================================================================"
echo " 7. Running CPPWRF C++ CUDA Driver"
echo "======================================================================"
"${HEAVY}" "${BUILD_DIR}/driver_diffusion_cuda" "${INPUT_BIN}" "${CUDA_OUT}"

echo "======================================================================"
echo " 8. Comparing Fortran vs C++ CPU Outputs (Tolerance 1e-12)"
echo "======================================================================"
python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPP_OUT}" --tol 1e-12

echo "======================================================================"
echo " 9. Comparing Fortran vs C++ CUDA Outputs (Tolerance 1e-12)"
echo "======================================================================"
python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CUDA_OUT}" --tol 1e-12
