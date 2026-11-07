#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# New Tiedtke column replay: CPPWRF's cu_ntiedtke_run inputs/outputs of P7 (steps 2 and 4) through WRF
# v4.8.0's cu_ntiedtke_run (the -r8 -fdefault-double-8 build that made the P7 references, WRF_R8D8_DIR;
# tests/phys/wrf_fortran.sh); every output of every column identical (same bits).
#   usage: run_test.sh [BUILD_DIR]   (needs CPPWRF_WRF_RUNS and WRF_R8D8_DIR, a WRF -r8 -double-8 build)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
BUILD_DIR=${1:-${REPO_ROOT}/build}
. "${SCRIPT_DIR}/../wrf_fortran.sh"
wrf_need r8d8
[ -n "${CPPWRF_WRF_RUNS:-}" ] || { echo "SKIP: set CPPWRF_WRF_RUNS"; exit 77; }
WORK=${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/ntiedtke_test
rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"
wrf_fortran "$WORK/harness" "$SCRIPT_DIR/harness.f90"
fail=0
for n in 2 4; do
  CPPWRF_NTIEDTKE_DUMP=nt$n.bin CPPWRF_NTIEDTKE_DUMP_STEP=$n OMP_NUM_THREADS=1 \
    "$BUILD_DIR/cppwrf" "$REPO_ROOT/configs/real_p7_r8.ini" --set run.run_time=$((60 * n)) --set output.dir="$WORK/out" > run$n.log 2>&1
  r=$(./harness nt$n.bin | tail -1)
  echo "step $n: $r"
  echo "$r" | awk '{exit !($2 > 0 && $4 > 0 && $NF == 0)}' || fail=1  # columns > 0, some active, none mismatching
done
cd /; rm -rf "$WORK"
[ $fail = 0 ] && echo "PASS" || { echo "FAIL"; exit 1; }
