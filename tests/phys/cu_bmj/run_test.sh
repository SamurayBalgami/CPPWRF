#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Column differential test for the Betts-Miller-Janjic cumulus scheme (cu_physics = 2): BMJINIT's tables and BMJDRV
# of WRF v4.8.0's own compiled module_cu_bmj.o (the -r8 build, WRF_DIR; tests/phys/wrf_fortran.sh) against
# cu_bmj::init_bmj_tables and cu_bmj::bmjdrv_column (src/phys/cu_bmj), on synthetic soundings (generate_inputs.py:
# deep, shallow, dry, cold columns, land and sea, CLDEFI in [EFIMN, 1]) with and without bmj_rad_feedback.
# Pass criterion: bit-identical (tol = 0.0) in CPU FP64; the FP32 driver is run as a smoke test.
#   run_test.sh [WORK_DIR]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
HEAVY="${REPO_ROOT}/tools/heavy.sh"
. "${SCRIPT_DIR}/../wrf_fortran.sh"
wrf_need r8
OUT_DIR="${1:-${PWD}/test_cu_bmj}"
mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

echo "=== Building the Fortran and C++ drivers (linking ${W}/main/libwrflib.a) ==="
wrf_fortran "${OUT_DIR}/driver_bmj_fortran" "${SCRIPT_DIR}/driver_bmj.f90"
for prec in fp64 fp32; do
  def=; [ "$prec" = fp64 ] && def=-DCPPWRF_DOUBLE
  "${HEAVY}" --light g++ -std=c++20 ${def} -O3 -march=native -ffp-contract=off -fopenmp -I"${REPO_ROOT}/src" \
    "${SCRIPT_DIR}/driver_bmj.cpp" "${REPO_ROOT}/src/core/backend.cpp" -lmvec -lm -o "driver_bmj_cpu_${prec}"
done

# columns, levels, stepcu, bmj_rad_feedback; odd column / level counts for the odd tails of WRF's vectorised loops
for cfg in "2000 44 1 0" "2001 45 3 1" "1999 60 1 1" "2002 33 6 0"; do
  set -- $cfg
  echo "=== BMJ differential test ($1 columns x $2 levels, stepcu $3, bmj_rad_feedback $4) ==="
  python3 "${SCRIPT_DIR}/generate_inputs.py" --cols "$1" --nz "$2" --stepcu "$3" --feedback "$4" --seed "$2" \
    --out bmj_in.bin
  "${HEAVY}" --light ./driver_bmj_fortran bmj_in.bin bmj_wrf.bin
  "${HEAVY}" --light ./driver_bmj_cpu_fp64 bmj_in.bin bmj_cpp.bin
  python3 "${SCRIPT_DIR}/compare.py" bmj_wrf.bin bmj_cpp.bin --input bmj_in.bin --tol 0.0
done

echo "=== FP32 smoke run ==="
"${HEAVY}" --light ./driver_bmj_cpu_fp32 bmj_in.bin bmj_fp32.bin
python3 "${SCRIPT_DIR}/compare.py" bmj_wrf.bin bmj_fp32.bin --input bmj_in.bin --tol 1e30 > /dev/null

echo "=== BMJ COLUMN DIFFERENTIAL TEST PASSED (BIT-IDENTICAL IN CPU FP64) ==="
