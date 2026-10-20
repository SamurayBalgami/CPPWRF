#!/usr/bin/env bash
# Data tar of tools/gcp/gcp_matrix.sh: the input files and tables of the given WRF run directories (links
# followed; outputs, executables and the large lookup tables no ported scheme reads left out) under runs/, and
# the WRF executables under wrf/.  usage: stage_matrix_data.sh OUT.tar WRF_R8_EXE WRF_R8D8_EXE RUN_DIR ...
set -euo pipefail
out=$1 r8=$2 r8d8=$3; shift 3
T=$(mktemp -d "${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/matrixdata.XXXX"); mkdir -p "$T/runs" "$T/wrf"
cp -L "$r8" "$T/wrf/wrf_r8.exe"; cp -L "$r8d8" "$T/wrf/wrf_r8d8.exe"
for d in "$@"; do
  n=$(basename "$d"); mkdir -p "$T/runs/$n"
  for f in "$d"/*; do
    b=$(basename "$f")
    case $b in wrfout*|wrfrst*|rsl.*|*.log|*.exe|steps|namelist.output|freezeH2O.dat|ishmael-*|p3_lookupTable*|snicar_*|aerosol*|CAM_ABS_DATA|CAM_AEROPT_DATA|create_p3*|*.nc) ;;
      *) cp -L "$f" "$T/runs/$n/" 2>/dev/null || true ;;
    esac
  done
done
tar cf "$out" -C "$T" runs wrf; rm -rf "$T"; ls -la "$out"
