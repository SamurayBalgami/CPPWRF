#!/usr/bin/env bash
# WRF v4.8.0 runs resumed from a restart file, for the wrfrst-restart test of run_wrf_files_test.sh:
# ${CPPWRF_WRF_RUNS}/rsm_<C> = rst_<C> (make_restart_refs.sh) restarted at minute 1, history every minute.
#   usage: make_resumed_refs.sh P1_r8 ...
set -uo pipefail
: "${CPPWRF_WRF_RUNS:?set CPPWRF_WRF_RUNS}"
for c in "$@"; do
  src=${CPPWRF_WRF_RUNS}/rst_$c dst=${CPPWRF_WRF_RUNS}/rsm_$c
  [ -d "$src" ] || { echo "$c: no $src"; continue; }
  rm -rf "$dst"; mkdir -p "$dst"; cd "$dst" || exit 1
  for f in "$src"/*; do
    b=$(basename "$f")
    case $b in wrfout*|rsl*|wrfrst*|namelist.input|namelist.output|*.log) ;; *) ln -s "$(readlink -f "$f")" "$b" ;; esac
  done
  ln -s "$src/wrfrst_d01_2026-10-12_00:01:00" .
  [ -f "$src/wrfrst_d02_2026-10-12_00:01:00" ] && ln -s "$src/wrfrst_d02_2026-10-12_00:01:00" .
  rm=$(sed -nE 's/^ *run_minutes *= *([0-9]+).*/\1/p' "$src/namelist.input")  # the same end time
  sed -E -e 's/^( *restart *= *)\.false\./\1.true./' -e 's/^( *start_minute *= *)[0-9]+/\101/' \
         -e "s/^( *run_minutes *= *)[0-9]+/\\1$((rm - 1))/" -e 's/^( *restart_interval *= *)[0-9]+/\11000/' \
         "$src/namelist.input" > namelist.input
  grep -q start_minute namelist.input || sed -i 's/^\( *start_hour *=.*\)$/\1\n start_minute = 01, 01, 01,/' namelist.input
  # MALLOC_PERTURB_=128: see make_restart_refs.sh (GFL reads an uninitialised buffer)
  MALLOC_PERTURB_=128 "$(readlink -f "$src/wrf.exe")" > wrf.log 2>&1
  echo "$c: $(grep -c 'SUCCESS COMPLETE WRF' wrf.log rsl.error.0000 2>/dev/null | tr '\n' ' ') $(ls | grep -c wrfout) history files"
done
