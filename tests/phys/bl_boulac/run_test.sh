#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# BouLac PBL (bl_pbl_physics = 8) column test: WRF v4.8.0's BOULAC (module_bl_boulac.o of the -r8 build, WRF_DIR)
# against the C++ BouLac (driver_boulac_cpu, built here with the C++ compiler), tolerance 0.
#   run_test.sh WORK_DIR
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
. "${SCRIPT_DIR}/../wrf_fortran.sh"
wrf_need r8
mkdir -p "${1:?usage: run_test.sh WORK_DIR}"
RUN_DIR="$(cd "$1" && pwd)"
cd "${RUN_DIR}"
wrf_fortran "${RUN_DIR}/driver_boulac_fortran" "${SCRIPT_DIR}/driver_boulac.f90"
g++ -std=c++20 -O3 -march=native -ffp-contract=off -DCPPWRF_DOUBLE -I"${REPO_ROOT}/src" \
  "${SCRIPT_DIR}/driver_boulac.cpp" -o "${RUN_DIR}/driver_boulac_cpu"
# columns, levels, dt, qc present; 44 levels: the K = 48 arrays, 60: kMaxLevels
for cfg in "3000 44 60 1" "3000 44 180 0" "2000 60 30 1" "1500 20 60 1"; do
  set -- $cfg
  echo "=== BouLac column test: $1 columns x $2 levels, dt $3, qc $4 ==="
  python3 "${SCRIPT_DIR}/generate_inputs.py" --cols "$1" --nz "$2" --dt "$3" --qc "$4" --seed "$2" --out boulac_in.bin
  ./driver_boulac_fortran boulac_in.bin boulac_out_fortran.bin
  ./driver_boulac_cpu boulac_in.bin boulac_out_cpp.bin
  python3 "${SCRIPT_DIR}/compare.py" boulac_out_fortran.bin boulac_out_cpp.bin --cols "$1" --nz "$2"
done
echo "=== BouLac column test PASSED ==="
