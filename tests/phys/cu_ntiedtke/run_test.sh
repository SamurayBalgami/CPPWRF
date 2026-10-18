#!/usr/bin/env bash
# New Tiedtke column replay: CPPWRF's cu_ntiedtke_run inputs/outputs of P7 (steps 2 and 4) through WRF
# v4.8.0's cu_ntiedtke_run (r8d8 objects); every column identical.
#   usage: run_test.sh [BUILD_DIR]   (needs CPPWRF_WRF_RUNS and WRF_R8D8_DIR, a WRF -r8 -double-8 build)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
BUILD_DIR=${1:-${REPO_ROOT}/build}
W=${WRF_R8D8_DIR:?set WRF_R8D8_DIR to a WRF v4.8.0 -r8 -double-8 build}
: "${CPPWRF_WRF_RUNS:?set CPPWRF_WRF_RUNS}"
WORK=${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/ntiedtke_test
rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"
gfortran -O2 -fdefault-real-8 -fdefault-double-8 -I"$W/phys" -o harness "$SCRIPT_DIR/harness.f90" \
  "$W/phys/physics_mmm/cu_ntiedtke.o" "$W/phys/ccpp_kind_types.o"
fail=0
for n in 2 4; do
  CPPWRF_NTIEDTKE_DUMP=nt$n.bin CPPWRF_NTIEDTKE_DUMP_STEP=$n OMP_NUM_THREADS=1 \
    "$BUILD_DIR/cppwrf" "$REPO_ROOT/configs/real_p7_r8.ini" --set run.run_time=$((60 * n)) --set output.dir="$WORK/out" > run$n.log 2>&1
  r=$(./harness nt$n.bin | tail -1)
  echo "step $n: $r"
  [ "$(echo "$r" | awk '{print $NF}')" = 0 ] || fail=1
done
cd /; rm -rf "$WORK"
[ $fail = 0 ] && echo "PASS" || { echo "FAIL"; exit 1; }
