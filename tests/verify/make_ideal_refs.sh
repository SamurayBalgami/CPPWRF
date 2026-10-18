#!/usr/bin/env bash
# WRF v4.8.0 real*8 references of the ideal cases: the single-precision reference runs in
# ${CPPWRF_WRF_REF}/run_<name> (namelist.input, input_sounding, tables) rerun with real*8 ideal.exe and
# wrf.exe (a WRF -r8 build without nesting, one per case: ${WRF_IDEAL_BIN}/<em_case>/{ideal,wrf}.exe)
# into ${CPPWRF_WRF_REF}/r8/run_<name>.
#   usage: make_ideal_refs.sh [name ...]   (default: every case below)
set -uo pipefail
: "${CPPWRF_WRF_REF:?set CPPWRF_WRF_REF}" "${WRF_IDEAL_BIN:?set WRF_IDEAL_BIN}"
declare -A CASE=([hill_p0]=em_hill2d_x [hill_p0_steps]=em_hill2d_x [straka_d1]=em_grav2d_x
                 [qss]=em_quarter_ss [qss_steps]=em_quarter_ss [squall2d_kessler]=em_squall2d_x)
names=("$@"); [ ${#names[@]} -eq 0 ] && names=(hill_p0 hill_p0_steps straka_d1 qss qss_steps squall2d_kessler)
for n in "${names[@]}"; do
  src=$CPPWRF_WRF_REF/run_$n dst=$CPPWRF_WRF_REF/r8/run_$n bin=$WRF_IDEAL_BIN/${CASE[$n]}
  [ -d "$src" ] && [ -x "$bin/wrf.exe" ] || { echo "$n: missing $src or $bin"; continue; }
  rm -rf "$dst"; mkdir -p "$dst"; cd "$dst" || exit 1
  for f in "$src"/*; do
    b=$(basename "$f")
    case $b in wrfout*|wrfrst*|wrfinput*|*.exe|*.log|namelist.output|*.npz) ;; *) ln -s "$(readlink -f "$f")" "$b" ;; esac
  done
  ln -s "$bin/ideal.exe" ideal.exe; ln -s "$bin/wrf.exe" wrf.exe
  ./ideal.exe > ideal.log 2>&1 && ./wrf.exe > wrf.log 2>&1
  echo "$n: $(grep -c 'SUCCESS COMPLETE' ideal.log wrf.log | tr '\n' ' ') $(ls | grep -c wrfout) wrfout"
done
