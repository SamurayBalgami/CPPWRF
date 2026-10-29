#!/usr/bin/env bash
# Startup script of the WPS VM (tools/gcp/gcp_wps.sh): downloads the GFS files, runs ungrib and metgrid, checks the
# met_em files that the WPS run directory also had, uploads the met_em files as one tar.
md() { curl -s -H "Metadata-Flavor: Google" "http://metadata.google.internal/computeMetadata/v1/instance/attributes/$1"; }
B=$(md bucket) SRC=$(md src) RES=$(md result) CYCLE=$(md cycle) HOURS=$(md hours) OBJ=$(md obj)
J=/opt/job; mkdir -p $J; exec > $J/job.log 2>&1
put() { gcloud storage cp "$@" "gs://$B/$RES/" || true; }
trap 'put $J/job.log $J/wps.out $J/ungrib.log $J/metgrid.log; date > $J/done; put $J/done; poweroff' EXIT
set -x
export DEBIAN_FRONTEND=noninteractive
apt-get -qq update && apt-get -qq install -y libgfortran5 libgomp1 libnetcdff7 libpng16-16 curl > /dev/null
gcloud storage cp "gs://$B/$SRC" $J/wps.tar.gz && tar xzf $J/wps.tar.gz -C $J
W=$J/wps; cd $W || exit 1
first=${HOURS%-*} last=${HOURS#*-}
iv=$(sed -nE 's/^ *interval_seconds *= *([0-9]+).*/\1/p' namelist.wps)
c0="${CYCLE:0:4}-${CYCLE:4:2}-${CYCLE:6:2} ${CYCLE:8:2}:00"
url=https://noaa-gfs-bdp-pds.s3.amazonaws.com/gfs.${CYCLE:0:8}/${CYCLE:8:2}/atmos
n=0 L=ABCDEFGHIJKLMNOPQRSTUVWXYZ
for ((h = first; h <= last; h += iv / 3600)); do
  f=gfs.t${CYCLE:8:2}z.pgrb2.0p25.f$(printf %03d $h)
  curl -sf --retry 5 -o $J/$f "$url/$f" || { echo "DOWNLOAD-FAILED $f" >> $J/wps.out; exit 0; }
  # link_grib.csh's names: GRIBFILE.AAA, AAB, ...
  ln -s $J/$f "GRIBFILE.A${L:$((n / 26)):1}${L:$((n % 26)):1}"
  n=$((n + 1))
done
s0=$(date -ud "$c0" +%s)
start=$(date -ud "@$((s0 + first * 3600))" +%Y-%m-%d_%H:00:00) end=$(date -ud "@$((s0 + last * 3600))" +%Y-%m-%d_%H:00:00)
nd=$(sed -nE 's/^ *max_dom *= *([0-9]+).*/\1/p' namelist.wps) sd= ed=
for ((i = 0; i < nd; i++)); do sd="$sd'$start', " ed="$ed'$end', "; done
sed -i -E -e "s/^ *start_date *=.*/ start_date = $sd/" -e "s/^ *end_date *=.*/ end_date   = $ed/" namelist.wps
cat namelist.wps
./ungrib.exe > $J/ungrib.log 2>&1
grep -q "Successful completion of ungrib" $J/ungrib.log || { echo "UNGRIB-FAILED" >> $J/wps.out; exit 0; }
rm -f $J/gfs.*
./metgrid.exe > $J/metgrid.log 2>&1
grep -q "Successful completion of metgrid" $J/metgrid.log || { echo "METGRID-FAILED" >> $J/wps.out; exit 0; }
for f in $J/check/met_em*; do
  [ -e "$f" ] || continue
  if cmp -s "$f" "$W/$(basename "$f")"; then echo "CHECK $(basename "$f") identical"; else echo "CHECK $(basename "$f") DIFFERS"; fi
done >> $J/wps.out
ls met_em* >> $J/wps.out
d=$(basename "$OBJ" .tar); mkdir -p $J/runs/$d && cp met_em* $J/runs/$d/ && tar cf $J/met_em.tar -C $J runs && gcloud storage cp $J/met_em.tar "gs://$B/$OBJ" && echo "STAGED gs://$B/$OBJ" >> $J/wps.out
exit 0
