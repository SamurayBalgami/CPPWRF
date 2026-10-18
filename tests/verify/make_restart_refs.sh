#!/usr/bin/env bash
# WRF v4.8.0 restart references for tests/verify/run_wrf_files_test.sh: reruns the 3-step run
# ${CPPWRF_WRF_RUNS}/run_<C>/steps with restart_interval = 1 (minute) in ${CPPWRF_WRF_RUNS}/rst_<C>.
#   usage: make_restart_refs.sh P1_r8 P3_r8 ...
set -uo pipefail
: "${CPPWRF_WRF_RUNS:?set CPPWRF_WRF_RUNS}"
for c in "$@"; do
  src=${CPPWRF_WRF_RUNS}/run_$c/steps dst=${CPPWRF_WRF_RUNS}/rst_$c
  [ -d "$src" ] || { echo "$c: no $src"; continue; }
  rm -rf "$dst"; mkdir -p "$dst"; cd "$dst" || exit 1
  # inputs and tables as links under their own names (the tables are looked up by name)
  for f in "$src"/*; do
    b=$(basename "$f")
    case $b in wrfout*|rsl*|wrfrst*|namelist.input|namelist.output|*.log) ;; *) ln -s "$(readlink -f "$f")" "$b" ;; esac
  done
  sed -E 's/^( *restart_interval *= *)[0-9]+/\11/' "$src/namelist.input" > namelist.input
  grep -q restart_interval namelist.input || sed -i 's/^ *history_interval/ restart_interval = 1,\n&/' namelist.input
  # MALLOC_PERTURB_=128: fills malloc'd buffers with 0x7f (finite float/double) so 0.0*cupclwm
  # in GFL's module_cu_gfl_wrf.F90 (cupclwm uninitialised when imid=0) does not hit 0.0*NaN
  # from negative 64-bit integers (0xffffffffffffff..) freed by HDF5 during wrfrst writes.
  MALLOC_PERTURB_=128 "$(readlink -f "$src/wrf.exe")" > wrf.log 2>&1
  echo "$c: $(grep -c 'SUCCESS COMPLETE WRF' wrf.log rsl.error.0000 2>/dev/null | tr '\n' ' ') $(ls | grep -c wrfrst) restart files"
done
