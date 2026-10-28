#!/usr/bin/env bash
# CPPWRF wrfout_d0N and wrfrst_d0N files against WRF v4.8.0's (CPU FP64 build, WRF real*8 runs):
# every variable, value (tolerance 0), dimension, attribute and global attribute must be identical.
#   usage: run_wrf_files_test.sh [BUILD_DIR] [case ...]
# CPPWRF_FILES_TRANSLATED=1: the CPPWRF configuration is made from the WRF run's namelist.input by
# tools/wrf2cppwrf.py instead of configs/ (the user's path from a WRF run directory).
# References: ${CPPWRF_WRF_RUNS}/run_<C>_r8[/steps] (history every step) and
# ${CPPWRF_WRF_RUNS}/rst_<C>_r8 (tests/verify/make_restart_refs.sh: restart every minute) and
# ${CPPWRF_WRF_RUNS}/rsm_<C>_r8 (make_resumed_refs.sh: WRF resumed from rst_<C>_r8 at minute 1; CPPWRF
# resumes from the same wrfrst file).
# RQCBLTEN / RQIBLTEN are not compared when WRF's file has no QCLOUD / QICE: without qc / qi WRF leaves them
# uninitialised. Neither are QZ0 and USTM of f1 (fractional sea ice with sfclayrev): WRF averages them with
# uninitialised open-water locals.
# MYNN's first step depends on a restart alarm ringing at its end (solve_em's restart_flag), so its history
# with restart_interval 1 differs from the plain run: mynn (no restart files) and mynnr (rst_/rsm_ only).
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
if [ -z "${CPPWRF_WRF_RUNS:-}" ] || [ ! -d "${CPPWRF_WRF_RUNS}" ]; then
  echo "SKIP: \${CPPWRF_WRF_RUNS} not set"; exit 77
fi
BUILD_DIR=${1:-${REPO_ROOT}/build}; shift || true
EXE=${BUILD_DIR}/cppwrf
CMP="python3 ${REPO_ROOT}/tests/init_real/compare_real.py --ignore"
# the variables not compared in WRF's file $1 (case $2): RQCBLTEN / RQIBLTEN without QCLOUD / QICE, f1's QZ0 / USTM
ignored() {
  python3 -c 'import sys, netCDF4
v = netCDF4.Dataset(sys.argv[1]).variables
print(",".join([x for x, q in (("RQCBLTEN", "QCLOUD"), ("RQIBLTEN", "QICE")) if q not in v]
               + (["QZ0", "USTM"] if sys.argv[2] == "f1" else []) or ["-"]))' "$1" "$2"
}
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
mynn|configs/real_mynn_r8.ini|$R/run_MYNN_r8/steps|-|
mynnr|configs/real_mynn_r8.ini|$R/rst_MYNN_r8|$R/rst_MYNN_r8|
r1|configs/real_r1.ini|$R/run_R1_r8|-|
t1|configs/real_t1.ini|$R/run_T1_r8|-|--set output.wrf_namelist=$R/run_T1_r8/namelist.input
n1|configs/real_n1_steps.ini|$R/run_N1_r8/steps|-|
n2|configs/real_n2_steps.ini|$R/run_N2_r8/steps|$R/rst_N2_r8|
n2s1|configs/real_n2_steps.ini|$R/run_N2S1_r8/steps|-|--set real.dir=$R/run_N2S1_r8/steps --set domains.smooth_option=1
n2s2|configs/real_n2_steps.ini|$R/run_N2S2_r8/steps|-|--set real.dir=$R/run_N2S2_r8/steps --set domains.smooth_option=2
n3|configs/real_n3_steps.ini|$R/run_N3_r8/steps|$R/rst_N3_r8|
nmp1|configs/real_nmp1_r8.ini|$R/run_NMP1_r8/steps|$R/rst_NMP1_r8|
v1|configs/real_v1_r8.ini|$R/run_V1_r8/steps|-|--set run.run_time=180
v2|configs/real_v2_r8.ini|$R/run_V2_r8/steps|-|--set run.run_time=180
v3|configs/real_v3_r8.ini|$R/run_V3_r8/steps|-|--set run.run_time=54
t6|configs/real_t6_r8.ini|$R/run_T6_r8/steps|$R/rst_T6_r8|
g2|configs/real_g2_r8.ini|$R/run_G2_r8/steps|$R/rst_G2_r8|
f1|configs/real_f1_r8.ini|$R/run_F1_r8/steps|$R/rst_F1_r8|
f2|configs/real_f2_r8.ini|$R/run_F2_r8/steps|$R/rst_F2_r8|
d0|configs/real_d0_r8.ini|$R/run_D0_r8/steps|$R/rst_D0_r8|
d0v|configs/real_d0v_r8.ini|$R/run_D0V_r8/steps|$R/rst_D0V_r8|
su1|configs/real_su1_r8.ini|$R/rst_SU1_r8|$R/rst_SU1_r8|--set run.run_time=180
hy1|configs/real_p1_r8.ini|$R/run_HY1_r8/steps|-|--set real.dir=$R/run_HY1_r8/steps --set dyn.hypsometric_opt=1
"
want=" $* "
fail=0 n=0
while IFS='|' read -r name cfg ref rst extra; do
  [ -z "$name" ] && continue
  [ $# -gt 0 ] && [[ "$want" != *" $name "* ]] && continue
  [ -d "$ref" ] || { echo "SKIP $name: $ref not found"; continue; }
  sets=(--set output.format=wrfout --set output.out_interval=1 --set output.dir="$WORK/$name")
  if { [ "$rst" != "-" ] && [ -d "$rst" ]; } || compgen -G "$ref/wrfrst_d0*" > /dev/null; then
    sets+=(--set restart.interval=60)
    # restart files only up to WRF's last one (they take ~0.5 GB each for the larger domains)
    last=$( (ls "$rst"/wrfrst_d01_* "$ref"/wrfrst_d01_* 2>/dev/null || true) | sed -E 's/.*_([0-9]{2}):([0-9]{2}):([0-9]{2})$/\1 \2 \3/' |
            awk '{t = $1 * 3600 + $2 * 60 + $3; if (t > m) m = t} END {print m + 0}')
    [ "$last" -gt 0 ] && sets+=(--set restart.end="$last")
  fi
  # shellcheck disable=SC2206
  [ -n "$extra" ] && sets+=($extra)
  rm -rf "${WORK:?}/$name"; mkdir -p "$WORK/$name"
  if [ "${CPPWRF_FILES_TRANSLATED:-0}" = 1 ]; then
    if ! python3 "$REPO_ROOT/tools/wrf2cppwrf.py" "$ref/namelist.input" -o "$WORK/$name/translated.ini" \
         --format wrfout > "$WORK/$name/translate.log" 2>&1; then
      echo "FAIL $name: wrf2cppwrf failed"; tail -3 "$WORK/$name/translate.log"; fail=1; continue
    fi
    cfg=$WORK/$name/translated.ini
  fi
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
      $CMP "$(ignored "$f" "$name")" "$ours" "$f" > "$WORK/$name/cmp.log" 2>&1 || { ok=0; echo "FAIL $name $b"; head -20 "$WORK/$name/cmp.log"; }
      n=$((n + 1))
    done
  done
  # resumed from WRF's restart file at minute 1 against WRF's resumed run (make_resumed_refs.sh)
  rsm=${rst/\/rst_/\/rsm_}
  if [ "$rst" != "-" ] && [ -d "$rsm" ]; then
    rdir=$WORK/$name/resumed; mkdir -p "$rdir"
    rsets=("${sets[@]}" --set output.dir="$rdir")
    for r in "$rst"/wrfrst_d01_*_00:01:00; do rsets+=(--set run.restart_file="$r"); done
    if ! "$EXE" "$cfg" "${rsets[@]}" > "$rdir/run.log" 2>&1; then
      echo "FAIL $name: resumed run failed"; tail -3 "$rdir/run.log"; ok=0
    else
      for f in "$rsm"/wrfout_d0*; do
        b=$(basename "$f"); ours=$rdir/$b; [ -f "$ours" ] || ours=$rdir/${b:7:3}/$b
        $CMP "$(ignored "$f" "$name")" "$ours" "$f" > "$rdir/cmp.log" 2>&1 || { ok=0; echo "FAIL $name resumed $b"; head -20 "$rdir/cmp.log"; }
        n=$((n + 1))
      done
    fi
  fi
  if [ $ok = 1 ]; then echo "PASS $name"; rm -rf "${WORK:?}/$name"; else fail=1; fi
done <<< "$CASES"
echo "compared $n files"
exit $fail
