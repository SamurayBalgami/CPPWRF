#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# A weather forecast from the latest GFS run: GFS -> WPS -> cppwrf_real -> cppwrf -> maps and city forecasts.
# One domain centred on LAT / LON (default: Turkey), WRF v4.8.0's CONUS physics suite, hourly output.
#
#   WPS_DIR=... WRF_DIR=... GEOG_DIR=... CPPWRF_BUILD=... examples/forecast/forecast.sh WORK
#
# WPS_DIR       a built WPS (geogrid.exe, ungrib.exe, metgrid.exe and its geogrid/, metgrid/, ungrib/Variable_Tables/)
# WRF_DIR       WRF v4.8.0's source tree (its run/ directory has the physics tables)
# GEOG_DIR      WPS's low-resolution static data (the directory WPS_GEOG_LOW_RES)
# CPPWRF_BUILD  the CPPWRF build directory (its cppwrf and src/init_real/cppwrf_real)
# MODE          gpu-fp32 (default), gpu-fp64 or cpu-fp64
# LAT, LON      domain centre (default 39.0, 35.0)
# NX, NY        grid points west-east, south-north (default 220, 130)
# DX            grid spacing in m (default 12000)
# HOURS         forecast length in hours, a multiple of 3 (default 24)
# CYCLE         GFS run YYYYMMDDHH (default: the latest one on NOAA's open data whose last needed hour is there)
# CITIES        city forecasts, "Name:lat:lon,..." (default: a list of Turkish cities, those inside the domain)
#
# WORK gets gfs/ (the GFS files), wps/ (geo_em, met_em), run/ (wrfinput, wrfbdy, wrfout_d01_*, rsl.*) and maps/
# (the maps of every hour, forecast.gif, cities.png, cities.txt).
set -euo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
WORK=${1:?usage: forecast.sh WORK}
: "${WPS_DIR:?set WPS_DIR}" "${WRF_DIR:?set WRF_DIR}" "${GEOG_DIR:?set GEOG_DIR}" "${CPPWRF_BUILD:?set CPPWRF_BUILD}"
MODE=${MODE:-gpu-fp32} LAT=${LAT:-39.0} LON=${LON:-35.0} NX=${NX:-220} NY=${NY:-130} DX=${DX:-12000}
HOURS=${HOURS:-24} CYCLE=${CYCLE:-}
[ $((HOURS % 3)) = 0 ] || { echo "HOURS must be a multiple of 3"; exit 2; }
B=$(cd "$CPPWRF_BUILD" && pwd)
export CPPWRF_EXE=$B/cppwrf CPPWRF_REAL_EXE=$B/src/init_real/cppwrf_real
mkdir -p "$WORK"/{gfs,wps,run,maps}
WORK=$(cd "$WORK" && pwd)
GFS=https://noaa-gfs-bdp-pds.s3.amazonaws.com
t0=$SECONDS
step() { echo "[$(date -u +%H:%M:%S)] $*"; }

# 1. the GFS run: the latest whose forecast hour HOURS is published (NOAA posts a run about 3.5 to 5 h after its time)
if [ -z "$CYCLE" ]; then
  now=$(date -u +%s)
  for back in 0 6 12 18 24 30; do
    c=$(date -u -d "@$(( (now - 3600 * (3 + back)) / 21600 * 21600 ))" +%Y%m%d%H)
    if curl -sfI "$GFS/gfs.${c:0:8}/${c:8:2}/atmos/gfs.t${c:8:2}z.pgrb2.0p25.f$(printf %03d "$HOURS")" > /dev/null; then
      CYCLE=$c; break
    fi
  done
  [ -n "$CYCLE" ] || { echo "no complete GFS run found in the last 36 h"; exit 1; }
fi
D=${CYCLE:0:8} H=${CYCLE:8:2}
START=$(date -u -d "${D} ${H}:00" +%Y-%m-%d_%H:00:00)
END=$(date -u -d "${D} ${H}:00 ${HOURS} hours" +%Y-%m-%d_%H:00:00)
step "GFS run $D ${H}z, forecast $START -> $END, ${NX}x${NY} points at $((DX / 1000)) km around $LAT N $LON E, $MODE"

# 2. the GFS files (0.25 degree, about 500 MB each), four at a time
cd "$WORK/gfs"
for f in $(seq 0 3 "$HOURS"); do printf '%03d\n' "$f"; done |
  xargs -P 4 -I{} sh -c "[ -s gfs.$CYCLE.f{} ] || curl -sf --retry 5 -o gfs.$CYCLE.f{}.part '$GFS/gfs.$D/$H/atmos/gfs.t${H}z.pgrb2.0p25.f{}' && mv -f gfs.$CYCLE.f{}.part gfs.$CYCLE.f{} 2>/dev/null || [ -s gfs.$CYCLE.f{} ]"
step "GFS: $(ls gfs.$CYCLE.f* | wc -l) files, $(du -sh . | cut -f1)"

# 3. WPS: geogrid, ungrib, metgrid on one Lambert domain
cd "$WORK/wps"
rm -f met_em.d01.* FILE:* GRIBFILE.*
cat > namelist.wps <<EOF
&share
 wrf_core = 'ARW',
 max_dom = 1,
 start_date = '$START',
 end_date   = '$END',
 interval_seconds = 10800
/
&geogrid
 parent_id = 1, parent_grid_ratio = 1, i_parent_start = 1, j_parent_start = 1,
 e_we = $NX, e_sn = $NY,
 geog_data_res = 'lowres',
 dx = $DX, dy = $DX,
 map_proj = 'lambert',
 ref_lat = $LAT, ref_lon = $LON,
 truelat1 = $LAT, truelat2 = $LAT, stand_lon = $LON,
 geog_data_path = '$GEOG_DIR'
/
&ungrib
 out_format = 'WPS',
 prefix = 'FILE',
/
&metgrid
 fg_name = 'FILE'
/
EOF
ln -sfn "$WPS_DIR/geogrid" geogrid; ln -sfn "$WPS_DIR/metgrid" metgrid
ln -sf "$WPS_DIR/ungrib/Variable_Tables/Vtable.GFS" Vtable
L=ABCDEFGHIJKLMNOPQRSTUVWXYZ n=0
for f in "$WORK"/gfs/gfs.$CYCLE.f*; do
  ln -sf "$f" "GRIBFILE.A${L:$((n / 26)):1}${L:$((n % 26)):1}"; n=$((n + 1))
done
for p in geogrid ungrib metgrid; do
  "$WPS_DIR/$p.exe" > "$p.out" 2>&1 || true
  grep -q "Successful completion of $p" "$p.out" || { tail -20 "$p.out"; echo "$p.exe failed: $WORK/wps/$p.out"; exit 1; }
done
# WPS 4.6.0's metgrid built by gfortran 15 at -O writes met_em files without fields (the notebook builds it at -O0)
for f in met_em.d01.*; do
  [ "$(stat -c %s "$f")" -gt 1000000 ] || { echo "$f has no fields: build metgrid with -O0"; exit 1; }
done
step "WPS: $(ls met_em.d01.* | wc -l) met_em files"

# 4. cppwrf_real (real.exe): wrfinput_d01, wrfbdy_d01
cd "$WORK/run"
rm -f wrfout_d01_* rsl.* met_em.d01.*
read -r NLEV NCAT < <(python3 -c "
import netCDF4, glob
d = netCDF4.Dataset(sorted(glob.glob('$WORK/wps/met_em.d01.*'))[0])
print(len(d.dimensions['num_metgrid_levels']), d.getncattr('NUM_LAND_CAT'))")
DT=$(( DX * 6 / 1000 ))          # WRF's rule: time step 6 x dx (km) seconds
RADT=$(( DX / 1000 ))            # radiation every dx (km) minutes
s=(${START//[-_:]/ }) e=(${END//[-_:]/ })
cat > namelist.input <<EOF
 &time_control
 run_days = 0, run_hours = $HOURS, run_minutes = 0, run_seconds = 0,
 start_year = ${s[0]}, start_month = ${s[1]}, start_day = ${s[2]}, start_hour = ${s[3]},
 end_year = ${e[0]}, end_month = ${e[1]}, end_day = ${e[2]}, end_hour = ${e[3]},
 interval_seconds = 10800
 input_from_file = .true.,
 history_interval = 60,
 frames_per_outfile = 1,
 restart = .false.,
 restart_interval = 100000,
 io_form_history = 2
 io_form_restart = 2
 io_form_input = 2
 io_form_boundary = 2
 /
 &domains
 time_step = $DT,
 time_step_fract_num = 0,
 time_step_fract_den = 1,
 max_dom = 1,
 e_we = $NX,
 e_sn = $NY,
 e_vert = 48,
 dzbot = 30.
 dzstretch_s = 1.11
 dzstretch_u = 1.10
 p_top_requested = 5000,
 num_metgrid_levels = $NLEV,
 num_metgrid_soil_levels = 4,
 dx = $DX,
 dy = $DX,
 grid_id = 1,
 parent_id = 0,
 i_parent_start = 1,
 j_parent_start = 1,
 parent_grid_ratio = 1,
 parent_time_step_ratio = 1,
 feedback = 1,
 smooth_option = 0
 /
 &physics
 physics_suite = 'CONUS'
 mp_physics = -1,
 cu_physics = -1,
 ra_lw_physics = -1,
 ra_sw_physics = -1,
 bl_pbl_physics = -1,
 sf_sfclay_physics = -1,
 sf_surface_physics = -1,
 radt = $RADT,
 bldt = 0,
 cudt = 0,
 icloud = 1,
 num_land_cat = $NCAT,
 sf_urban_physics = 0,
 fractional_seaice = 1,
 /
 &fdda
 /
 &dynamics
 hybrid_opt = 2,
 w_damping = 1,
 diff_opt = 2,
 km_opt = 4,
 diff_6th_opt = 0,
 diff_6th_factor = 0.12,
 base_temp = 290.
 damp_opt = 3,
 zdamp = 5000.,
 dampcoef = 0.2,
 khdif = 0,
 kvdif = 0,
 non_hydrostatic = .true.,
 moist_adv_opt = 1,
 scalar_adv_opt = 1,
 gwd_opt = 1,
 /
 &bdy_control
 spec_bdy_width = 5,
 specified = .true.
 /
 &grib2
 /
 &namelist_quilt
 nio_tasks_per_group = 0,
 nio_groups = 1,
 /
EOF
ln -sf "$WORK"/wps/met_em.d01.* .
python3 "$REPO/tools/cppwrf-realexe" > real.log 2>&1 || { tail -20 real.log; tail -5 rsl.error.0000 2>/dev/null; exit 1; }
step "cppwrf_real: $(ls wrfinput_d01 wrfbdy_d01 | tr '\n' ' ')"

# 5. cppwrf (wrf.exe) with WRF's physics tables
for f in "$WRF_DIR"/run/*; do
  case ${f##*/} in namelist.input|*.exe|README*) ;; *) ln -sfn "$f" . ;; esac
done
t1=$SECONDS
python3 "$REPO/tools/cppwrf-wrf" --mode "$MODE" > wrf.log 2>&1 || { tail -20 wrf.log; tail -5 rsl.error.0000 2>/dev/null; exit 1; }
grep -q "SUCCESS COMPLETE WRF" rsl.out.0000
step "cppwrf: $HOURS h forecast in $((SECONDS - t1)) s, $(ls wrfout_d01_* | wc -l) output files"

# 6. maps and city forecasts
python3 "$REPO/examples/forecast/maps.py" "$WORK/run" "$WORK/maps" ${CITIES:+--cities "$CITIES"}
step "done in $((SECONDS - t0)) s: $WORK/maps"
