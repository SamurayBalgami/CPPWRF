#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Grell-Freitas-Li column replay: CPPWRF's cu_grell_freitas_li inputs/outputs of P8 (steps 1, 5, 10;
# ishallow = 1) through WRF v4.8.0's cu_grell_freitas_li (the -r8 -fdefault-double-8 build that made the P8
# references, WRF_R8D8_DIR; tests/phys/wrf_fortran.sh); every value of every column identical (same bits).
#   usage: run_test.sh [BUILD_DIR]   (needs CPPWRF_WRF_RUNS and WRF_R8D8_DIR, a WRF -r8 -double-8 build)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
BUILD_DIR=${1:-${REPO_ROOT}/build}
. "${SCRIPT_DIR}/../wrf_fortran.sh"
wrf_need r8d8
[ -n "${CPPWRF_WRF_RUNS:-}" ] || { echo "SKIP: set CPPWRF_WRF_RUNS"; exit 77; }
WORK=${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/gfl_test
rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"
wrf_fortran "$WORK/harness" "$SCRIPT_DIR/harness.f90"
ulimit -s unlimited
fail=0
for n in 1 5 10; do
  CPPWRF_GFL_DUMP=gfl$n.bin CPPWRF_GFL_DUMP_STEP=$n \
    "$BUILD_DIR/cppwrf" "$REPO_ROOT/configs/real_p8_r8.ini" --set run.run_time=$((60 * n)) --set output.dir="$WORK/out" > run$n.log 2>&1
  r=$(./harness gfl$n.bin 1 0 60 15000 $n | tail -1)
  echo "step $n: $r"
  # at least one convecting column, none differing
  echo "$r" | awk '{exit !($3 > 0 && $6 == 0)}' || fail=1
done
cd /; rm -rf "$WORK"
[ $fail = 0 ] && echo "PASS" || { echo "FAIL"; exit 1; }
