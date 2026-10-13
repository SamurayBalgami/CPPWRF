#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"
HEAVY="${REPO_ROOT}/tools/heavy.sh"
TBL_DIR="/home/ostriquetrum/wrf_real/wrf_run"
WRFINPUT="/home/ostriquetrum/wrf_real/wrf_run/wrfinput_d01_hybrid0"

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

echo "======================================================================"
echo " 1. Building Fortran table reference reader"
echo "======================================================================"
"${HEAVY}" gfortran -ffree-form -fdefault-real-8 -O3 -march=native \
  "${SCRIPT_DIR}/fortran/read_tables_ref.f90" \
  -o "${BUILD_DIR}/read_tables_ref"

echo "======================================================================"
echo " 2. Running Fortran table reference reader"
echo "======================================================================"
"${HEAVY}" "${BUILD_DIR}/read_tables_ref" "${TBL_DIR}" "${BUILD_DIR}/tables_ref.bin"

echo "======================================================================"
echo " 3. Building C++ table differential test"
echo "======================================================================"
"${HEAVY}" g++ -std=c++20 -DCPPWRF_DOUBLE -O3 -march=native -fopenmp \
  -I "${REPO_ROOT}/src" \
  "${SCRIPT_DIR}/test_init_tables.cpp" \
  "${REPO_ROOT}/src/phys/driver/landuse_tables.cpp" \
  "${REPO_ROOT}/src/phys/driver/physics_init.cpp" \
  "${REPO_ROOT}/src/phys/lsm_noah/noah_tables.cpp" \
  "${REPO_ROOT}/src/phys/cu_kf/kf_tables.cpp" \
  "${REPO_ROOT}/src/core/backend.cpp" \
  -o "${BUILD_DIR}/test_init_tables"

echo "======================================================================"
echo " 4. Running C++ table differential test"
echo "======================================================================"
"${HEAVY}" "${BUILD_DIR}/test_init_tables" "${TBL_DIR}" "${BUILD_DIR}/tables_ref.bin"

echo "======================================================================"
echo " 5. Extracting real-data wrfinput fields (${WRFINPUT})"
echo "======================================================================"
python3 "${SCRIPT_DIR}/extract_wrfinput.py" \
  --nc "${WRFINPUT}" \
  --out "${BUILD_DIR}/wrfinput_extracted.bin"

echo "======================================================================"
echo " 6. Building C++ real-data smoke test"
echo "======================================================================"
"${HEAVY}" g++ -std=c++20 -DCPPWRF_DOUBLE -O3 -march=native -fopenmp \
  -I "${REPO_ROOT}/src" \
  "${SCRIPT_DIR}/test_real_smoke.cpp" \
  "${REPO_ROOT}/src/phys/driver/landuse_tables.cpp" \
  "${REPO_ROOT}/src/phys/driver/physics_init.cpp" \
  "${REPO_ROOT}/src/phys/lsm_noah/noah_tables.cpp" \
  "${REPO_ROOT}/src/phys/cu_kf/kf_tables.cpp" \
  "${REPO_ROOT}/src/core/backend.cpp" \
  "${REPO_ROOT}/src/core/halo.cpp" \
  -o "${BUILD_DIR}/test_real_smoke"

echo "======================================================================"
echo " 7. Running C++ real-data smoke test"
echo "======================================================================"
"${HEAVY}" "${BUILD_DIR}/test_real_smoke" "${TBL_DIR}" "${BUILD_DIR}/wrfinput_extracted.bin"

echo "======================================================================"
echo " ALL PHYSICS INITIALIZATION TESTS COMPLETED SUCCESSFULLY!"
echo "======================================================================"
