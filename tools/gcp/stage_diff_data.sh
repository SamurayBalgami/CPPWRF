#!/usr/bin/env bash
# Packs the data of the differential tests (ctest -L diff, which tools/gcp/remote_suite.sh runs on the full VM suite)
# as a tar on stdout, for the bucket:
#   tools/gcp/stage_diff_data.sh | gcloud storage cp - gs://BUCKET/diff_data.tar
#   wrf_build/r8, wrf_build/r8d8   the pieces of the WRF v4.8.0 build trees the tests link (tests/phys/wrf_fortran.sh):
#                                  configure.wrf, main/libwrflib.a, module_wrf_top.o, the .mod files, inc/, the
#                                  external libraries, frame's two loose objects and the run/ tables (links followed)
#   wrf_runs_diff/run_R0_r8, run_P4_r8   the first frame of their wrfout_d01_2026-10-12_00:00:00 (tests/nest,
#                                  tests/phys/bl_myj read the input columns from it)
#   wrf_runs_diff/run_P7_r8/steps, run_P8_r8/steps   their wrfinput_d01, wrfbdy_d01, namelist.input and *.TBL
#                                  (the model runs of the New Tiedtke / Grell-Freitas-Li replays)
#   wrf_runs_diff/wrf_run          its wrfinput_d01 and wrfbdy_d01 (tests/io_netcdf/test_wrf_readers)
# tools/gcp/run_diff_tests.sh uses a run of the suite data when it is there and these otherwise.
# Environment: WRF_DIR, WRF_R8D8_DIR (default ~/wrf_real/WRF-4.8.0-r8, -r8d8), CPPWRF_WRF_RUNS (default ~/wrf_real480).
# The two wrfout frames (~90 MB) are written to a temporary directory under ${CPPWRF_SCRATCH:-~/cppwrf_scratch}.
set -euo pipefail
R8=${WRF_DIR:-$HOME/wrf_real/WRF-4.8.0-r8}
R8D8=${WRF_R8D8_DIR:-$HOME/wrf_real/WRF-4.8.0-r8d8}
RUNS=${CPPWRF_WRF_RUNS:-$HOME/wrf_real480}
TMP=$(mktemp -d "${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/stage_diff.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

for b in r8 r8d8; do
  W=$R8; [ $b = r8d8 ] && W=$R8D8
  ( cd "$W"
    { echo configure.wrf; echo main/libwrflib.a; echo main/module_wrf_top.o
      echo frame/module_internal_header_util.o; echo frame/pack_utils.o
      find main frame share phys dyn_em external -name '*.mod'
      find external -name 'lib*.a'
      find inc -type f
      for t in LANDUSE.TBL VEGPARM.TBL SOILPARM.TBL GENPARM.TBL NoahmpTable.TBL CCN_ACTIVATE.BIN; do echo "run/$t"; done
    } | sort -u > "$TMP/list_$b" )
  mkdir -p "$TMP/wrf_build"; ln -s "$W" "$TMP/wrf_build/$b"
done

for r in run_R0_r8 run_P4_r8; do
  mkdir -p "$TMP/wrf_runs_diff/$r"
  python3 -I - "$RUNS/$r/wrfout_d01_2026-10-12_00:00:00" "$TMP/wrf_runs_diff/$r/wrfout_d01_2026-10-12_00:00:00" <<'EOF'
import sys
import netCDF4 as nc
src, dst = nc.Dataset(sys.argv[1]), nc.Dataset(sys.argv[2], "w", format=nc.Dataset(sys.argv[1]).data_model)
src.set_auto_maskandscale(False)
for n, d in src.dimensions.items():
    dst.createDimension(n, None if d.isunlimited() else len(d))
dst.setncatts({k: src.getncattr(k) for k in src.ncattrs()})
for n, v in src.variables.items():
    w = dst.createVariable(n, v.dtype, v.dimensions)
    w.set_auto_maskandscale(False)
    w.setncatts({k: v.getncattr(k) for k in v.ncattrs()})
    w[:] = v[0:1] if v.dimensions and v.dimensions[0] == "Time" else v[:]
dst.close()
EOF
done

for r in run_P7_r8 run_P8_r8; do
  mkdir -p "$TMP/wrf_runs_diff/$r/steps"
  for f in "$RUNS/$r/steps"/{wrfinput_d01,wrfbdy_d01,namelist.input} "$RUNS/$r/steps"/*.TBL; do
    [ -e "$f" ] && ln -s "$(readlink -f "$f")" "$TMP/wrf_runs_diff/$r/steps/$(basename "$f")"
  done
done

mkdir -p "$TMP/wrf_runs_diff/wrf_run"
for f in wrfinput_d01 wrfbdy_d01; do ln -s "$(readlink -f "$RUNS/wrf_run/$f")" "$TMP/wrf_runs_diff/wrf_run/$f"; done

cd "$TMP"
{ sed 's|^|wrf_build/r8/|' list_r8; sed 's|^|wrf_build/r8d8/|' list_r8d8; find -L wrf_runs_diff -type f; } > files
tar -ch -T files -f -
