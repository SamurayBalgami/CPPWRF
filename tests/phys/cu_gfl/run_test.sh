#!/usr/bin/env bash
# Grell-Freitas-Li column replay: CPPWRF's cu_grell_freitas_li inputs/outputs of P8 (steps 1, 5, 10;
# ishallow = 1) through WRF v4.8.0's cu_grell_freitas_li (r8d8 objects); every column identical.
#   usage: run_test.sh [BUILD_DIR]   (needs CPPWRF_WRF_RUNS and WRF_R8D8_DIR, a WRF -r8 -double-8 build)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
BUILD_DIR=${1:-${REPO_ROOT}/build}
W=${WRF_R8D8_DIR:?set WRF_R8D8_DIR to a WRF v4.8.0 -r8 -double-8 build}/phys
: "${CPPWRF_WRF_RUNS:?set CPPWRF_WRF_RUNS}"
WORK=${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/gfl_test
rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"
gfortran -O2 -fdefault-real-8 -fdefault-double-8 -I"$W" -o harness "$SCRIPT_DIR/harness.f90" \
  "$W/module_cu_gfl_wrf.o" "$W/module_cu_gfl_deep.o" "$W/module_cu_gfl_sh.o" \
  "$W/module_cu_gfl_common.o" "$W/ccpp_kind_types.o"
ulimit -s unlimited
fail=0
for n in 1 5 10; do
  CPPWRF_GFL_DUMP=gfl$n.bin CPPWRF_GFL_DUMP_STEP=$n \
    "$BUILD_DIR/cppwrf" "$REPO_ROOT/configs/real_p8_r8.ini" --set run.run_time=$((60 * n)) --set output.dir="$WORK/out" > run$n.log 2>&1
  r=$(./harness gfl$n.bin 1 0 60 15000 $n | tail -1)
  echo "step $n: $r"
  echo "$r" | grep -q "differing >0: *0 " || fail=1
done
cd /; rm -rf "$WORK"
[ $fail = 0 ] && echo "PASS" || { echo "FAIL"; exit 1; }
