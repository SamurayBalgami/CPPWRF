#!/usr/bin/env bash
# CPPWRF wrfout_d0N and wrfrst_d0N files against WRF v4.8.0's (CPU FP64 build, WRF real*8 runs):
# every variable, value (tolerance 0), dimension, attribute and global attribute must be identical.
#   usage: run_wrf_files_test.sh [BUILD_DIR] [case ...]
# References: ${CPPWRF_WRF_RUNS}/run_<C>_r8[/steps] (history every step) and
# ${CPPWRF_WRF_RUNS}/rst_<C>_r8 (tests/verify/make_restart_refs.sh: restart every minute).
# RQCBLTEN and RQIBLTEN are not compared: without qc/qi WRF leaves them uninitialised.
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
if [ -z "${CPPWRF_WRF_RUNS:-}" ] || [ ! -d "${CPPWRF_WRF_RUNS}" ]; then
  echo "SKIP: \${CPPWRF_WRF_RUNS} not set"; exit 0
fi
BUILD_DIR=${1:-${REPO_ROOT}/build}; shift || true
EXE=${BUILD_DIR}/cppwrf
CMP="python3 ${REPO_ROOT}/tests/init_real/compare_real.py --ignore RQCBLTEN,RQIBLTEN"
WORK=${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/wrf_files_test
R=${CPPWRF_WRF_RUNS}
cd "$REPO_ROOT"

# name | config | wrfout reference dir | restart reference dir (or -) | extra --set
CASES="
p1|configs/real_p1_r8.ini|$R/run_P1_r8/steps|$R/rst_P1_r8|
p2|configs/real_p2_r8.ini|$R/run_P2_r8/steps|-|
p3|configs/real_p3_r8.ini|$R/run_P3_r8/steps|$R/rst_P3_r8|
p4|configs/real_p4_r8.ini|$R/run_P4_r8/steps|$R/rst_P4_r8|
p5|configs/real_p5_r8.ini|$R/run_P5_r8|-|
p6|configs/real_p6_r8.ini|$R/run_P6_r8/steps|$R/rst_P6_r8|
p7|configs/real_p7_r8.ini|$R/run_P7_r8/steps|$R/rst_P7_r8|
p8|configs/real_p8_r8.ini|$R/run_P8_r8/steps|$R/rst_P8_r8|--set run.run_time=600
p9|configs/real_p9_r8.ini|$R/run_P9_r8/steps|$R/rst_P9_r8|
mp8|configs/real_mp8_r8.ini|$R/run_mp8_r8/steps|$R/rst_mp8_r8|
r1|configs/real_r1.ini|$R/run_R1_r8|-|
t1|configs/real_t1.ini|$R/run_T1_r8|-|--set output.wrf_namelist=$R/run_T1_r8/namelist.input
n1|configs/real_n1_steps.ini|$R/run_N1_r8/steps|-|
n2|configs/real_n2_steps.ini|$R/run_N2_r8/steps|$R/rst_N2_r8|
"
want=" $* "
fail=0 n=0
while IFS='|' read -r name cfg ref rst extra; do
  [ -z "$name" ] && continue
  [ $# -gt 0 ] && [[ "$want" != *" $name "* ]] && continue
  [ -d "$ref" ] || { echo "SKIP $name: $ref not found"; continue; }
  sets=(--set output.format=wrfout --set output.out_interval=1 --set output.dir="$WORK/$name")
  [ "$rst" != "-" ] && [ -d "$rst" ] && sets+=(--set restart.interval=60)
  # shellcheck disable=SC2206
  [ -n "$extra" ] && sets+=($extra)
  rm -rf "${WORK:?}/$name"; mkdir -p "$WORK/$name"
  if ! "$EXE" "$cfg" "${sets[@]}" > "$WORK/$name/run.log" 2>&1; then
    echo "FAIL $name: cppwrf failed"; tail -3 "$WORK/$name/run.log"; fail=1; continue
  fi
  ok=1
  for dir in "$ref" "$rst"; do
    [ "$dir" = "-" ] || [ ! -d "$dir" ] && continue
    for f in "$dir"/wrfout_d0* "$dir"/wrfrst_d0*; do
      [ -f "$f" ] || continue
      b=$(basename "$f"); dom=${b:7:3}
      ours=$WORK/$name/$b; [ -f "$ours" ] || ours=$WORK/$name/$dom/$b
      if [ ! -f "$ours" ]; then echo "FAIL $name: $b not written"; ok=0; continue; fi
      $CMP "$ours" "$f" > "$WORK/$name/cmp.log" 2>&1 || { ok=0; echo "FAIL $name $b"; head -20 "$WORK/$name/cmp.log"; }
      n=$((n + 1))
    done
  done
  if [ $ok = 1 ]; then echo "PASS $name"; rm -rf "${WORK:?}/$name"; else fail=1; fi
done <<< "$CASES"
echo "compared $n files"
exit $fail
