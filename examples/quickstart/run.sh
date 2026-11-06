#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# A real-data forecast from scratch: GFS analysis -> WPS -> cppwrf_real -> cppwrf -> wrfout (docs/USER_GUIDE.md,
# "End-to-end example"). The domains and the namelist are WRF's defaults (WPS 4.6's namelist.wps and WRF v4.8.0's
# test/em_real/namelist.input: a 15 km domain with a two-way 5 km nest, physics_suite CONUS) with the date of the
# GFS run, 3 forecast hours and WPS's low-resolution geography.
#
#   WPS_DIR=... WRF_DIR=... GEOG_DIR=... [CPPWRF_BUILD=build-cpu] examples/quickstart/run.sh WORK
#
# WPS_DIR   a built WPS (geogrid.exe, ungrib.exe, metgrid.exe and its geogrid/, metgrid/, ungrib/Variable_Tables/)
# WRF_DIR   WRF v4.8.0's source tree (its run/ directory has the physics tables; WRF itself is not needed)
# GEOG_DIR  WPS's low-resolution static data (geog_low_res_mandatory.tar.gz from the WPS download page, unpacked:
#           the directory WPS_GEOG_LOW_RES)
# CPPWRF_BUILD  the CPPWRF build directory with cppwrf and src/init_real/cppwrf_real (default: build-cpu in the
#           repository)
# GFS_DATE, GFS_CYCLE  the GFS run (default 20261012 00; NOAA's open data on AWS, two 0.25 degree files of
#           about 500 MB: the analysis and the 3 h forecast). Other dates: change the dates of both namelists too.
#
# WORK gets wps/ (the met_em files), run/ (namelist.input, wrfinput_d01, wrfinput_d02, wrfbdy_d01, run.ini, the
# tables) and out/d01/wrfout_d01_*, out/d02/wrfout_d02_* (a nested run writes each domain's files to output.dir/d0N;
# one file per hour, as WRF's frames_per_outfile = 1).
set -euo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EX=$REPO/examples/quickstart
WORK=${1:?usage: run.sh WORK}
: "${WPS_DIR:?set WPS_DIR}" "${WRF_DIR:?set WRF_DIR}" "${GEOG_DIR:?set GEOG_DIR}"
B=${CPPWRF_BUILD:-$REPO/build-cpu}
GFS_DATE=${GFS_DATE:-20261012} GFS_CYCLE=${GFS_CYCLE:-00}
mkdir -p "$WORK/wps" "$WORK/run" "$WORK/out"
WORK=$(cd "$WORK" && pwd) B=$(cd "$B" && pwd)

# 1. WPS: geogrid, ungrib (GFS through WPS's Vtable.GFS), metgrid
cd "$WORK/wps"
sed "s|GEOG_DIR|$GEOG_DIR|" "$EX/namelist.wps" > namelist.wps
ln -sfn "$WPS_DIR/geogrid" geogrid; ln -sfn "$WPS_DIR/metgrid" metgrid
ln -sf "$WPS_DIR/ungrib/Variable_Tables/Vtable.GFS" Vtable
L=ABCDEFGHIJ n=0
for h in 000 003; do
  [ -s "gfs.f$h" ] || curl -sf --retry 5 -o "gfs.f$h" \
    "https://noaa-gfs-bdp-pds.s3.amazonaws.com/gfs.$GFS_DATE/$GFS_CYCLE/atmos/gfs.t${GFS_CYCLE}z.pgrb2.0p25.f$h"
  ln -sf "gfs.f$h" "GRIBFILE.AA${L:$n:1}"; n=$((n + 1))
done
for p in geogrid ungrib metgrid; do
  "$WPS_DIR/$p.exe" > "$p.log" 2>&1
  grep -q "Successful completion of $p" "$p.log" || { echo "$p.exe failed: $WORK/wps/$p.log"; exit 1; }
done

# 2. cppwrf_real (real.exe): wrfinput_d01, wrfinput_d02, wrfbdy_d01. num_metgrid_levels and num_land_cat as in
# the met_em files (real.exe stops when they differ)
cd "$WORK/run"
read -r NLEV NCAT < <(python3 -c "
import netCDF4, glob
d = netCDF4.Dataset(sorted(glob.glob('$WORK/wps/met_em.d01.*'))[0])
print(len(d.dimensions['num_metgrid_levels']), d.getncattr('NUM_LAND_CAT'))")
sed -e "s/^ num_metgrid_levels .*/ num_metgrid_levels                  = $NLEV,/" \
    -e "s/^ num_land_cat .*/ num_land_cat                        = $NCAT,/" "$EX/namelist.input" > namelist.input
"$B/src/init_real/cppwrf_real" --namelist namelist.input --met-dir "$WORK/wps" --out-dir . > real.log 2>&1 ||
  { tail -5 real.log; echo "cppwrf_real failed: $WORK/run/real.log"; exit 1; }

# 3. the physics tables of WRF's run directory, the translated namelist and the model
for f in "$WRF_DIR"/run/*; do  # as for wrf.exe: everything but the namelist and the executables
  case ${f##*/} in namelist.input|*.exe|README*) ;; *) ln -sfn "$f" . ;; esac
done
python3 "$REPO/tools/wrf2cppwrf.py" namelist.input -o run.ini --format wrfout > translate.log 2>&1 ||
  { cat translate.log; exit 1; }
"$B/cppwrf" run.ini --set output.dir="$WORK/out" > cppwrf.log 2>&1 || { tail -5 cppwrf.log; echo "cppwrf failed: $WORK/run/cppwrf.log"; exit 1; }
ls "$WORK/out"/d0*/wrfout_d0*
