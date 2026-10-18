#!/usr/bin/env bash
# Reference wrfinput/wrfbdy files from WRF's real.exe for tests/init_real/run_real_test.sh.
# For every ${CPPWRF_WRF_RUNS}/run_<case> with a namelist.input and met_em files, runs real.exe of the
# -r8 build (cases named *_r8*) or the default single-precision build and keeps the result in
# ${CPPWRF_WRF_RUNS}/real480/<case>. A case directory without met_em files uses ${CPPWRF_WRF_RUNS}/wps_run;
# a run shorter than interval_seconds is extended to one interval (real_namelist.py).
#   usage: make_real_refs.sh WRF_SP_DIR WRF_R8_DIR [case ...]
set -uo pipefail
: "${CPPWRF_WRF_RUNS:?set CPPWRF_WRF_RUNS}"
SP=$1 R8=$2; shift 2
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HEAVY="${HERE}/../heavy.sh"
OUT=${CPPWRF_WRF_RUNS}/real480
cases=("$@")
[ ${#cases[@]} -eq 0 ] && for d in "${CPPWRF_WRF_RUNS}"/run_*/; do cases+=("$(basename "$d")"); done
for c in "${cases[@]}"; do
  src=${CPPWRF_WRF_RUNS}/$c
  [ -f "$src/namelist.input" ] || continue
  met=$src; ls "$met"/met_em.* >/dev/null 2>&1 || met=${CPPWRF_WRF_RUNS}/wps_run
  if [[ $c == *_r8* ]]; then W=$R8 P=double; else W=$SP P=float; fi
  w=$OUT/$c; rm -rf "$w"; mkdir -p "$w"; cd "$w" || exit 1
  python3 "$HERE/real_namelist.py" "$src/namelist.input" namelist.input
  for f in "$met"/met_em.*; do ln -s "$(readlink -f "$f")" .; done
  for f in "$W"/run/*; do
    case $(basename "$f") in namelist*|*.exe|wrfinput*|wrfbdy*|met_em*) ;; *) ln -s "$f" .;; esac
  done
  "$HEAVY" --light "$W/main/real.exe" > real.log 2>&1
  if grep -qs SUCCESS rsl.error.0000 real.log; then
    find . -maxdepth 1 -type l ! -name 'met_em.*' -delete
    rm -f rsl.* namelist.output
    echo $P > precision
    echo "$c: ok ($P)"
  else
    echo "$c: real.exe failed, skipped"; cd "$OUT" && rm -rf "$w"
  fi
done
