#!/usr/bin/env bash
# met_em files from the NOAA GFS 0.25 degree forecast on a Google Compute Engine spot VM (WPS ungrib + metgrid).
#
#   CPPWRF_GCP_PROJECT=... CPPWRF_GCP_BUCKET=... tools/gcp/gcp_wps.sh WPS_RUN_DIR --cycle YYYYMMDDHH
#       --hours FIRST-LAST --obj OBJ [--machine TYPE] [--zones Z,...] [--out DIR]
#
# WPS_RUN_DIR: a WPS run directory (namelist.wps, geo_em.d0N.nc, Vtable, metgrid/METGRID.TBL, ungrib.exe and
# metgrid.exe of the build that made its met_em files). The VM downloads gfs.tHHz.pgrb2.0p25.fFFF of the cycle for
# FIRST..LAST every interval_seconds of namelist.wps from NOAA's open data (noaa-gfs-bdp-pds), runs ungrib and
# metgrid with the namelist's domains over that period, and uploads the met_em files as the bucket object OBJ (a
# tar of runs/NAME/met_em*, NAME = OBJ without .tar: a data tar of gcp_wrfref.sh, run.real=EXE:NAME). A met_em file
# of WPS_RUN_DIR at a time the VM also makes must come out identical (cmp; the line "CHECK ... identical" of
# OUT/wps.out): the VM's WPS is the one that made the directory's files. Machine: tools/gcp/gcp_pick.py's cheapest
# Intel spot VM with 4 vCPUs (ungrib and metgrid are serial) and its zone now (--machine, --zones as gcp_wrfref.sh).
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
. "$REPO/tools/gcp/pick_lib.sh"
: "${CPPWRF_GCP_PROJECT:?set CPPWRF_GCP_PROJECT}" "${CPPWRF_GCP_BUCKET:?set CPPWRF_GCP_BUCKET}"
P=$CPPWRF_GCP_PROJECT BK=$CPPWRF_GCP_BUCKET
WPS=$(readlink -f "$1"); shift
MACHINE= ZONES= CYCLE= HOURS= OBJ= OUT=
while [ $# -gt 0 ]; do
  case $1 in
    --cycle) CYCLE=$2; shift 2 ;; --hours) HOURS=$2; shift 2 ;; --obj) OBJ=$2; shift 2 ;;
    --machine) MACHINE=$2; shift 2 ;; --zones) ZONES=$2; shift 2 ;; --out) OUT=$2; shift 2 ;;
    *) sed -n '2,15p' "$0"; exit 2 ;;
  esac
done
[ -n "$CYCLE" ] && [ -n "$HOURS" ] && [ -n "$OBJ" ] || { sed -n '2,15p' "$0"; exit 2; }
STAMP=$(date +%Y%m%d-%H%M%S)-$(printf %04x $RANDOM) VM=cppwrf-wps-$STAMP RES=results/wps-$STAMP
OUT=${OUT:-${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/gcp/wps-$STAMP}; mkdir -p "$OUT"
# the WPS inputs (links followed) and the directory's met_em files of the period, to check the VM's output against
T=$(mktemp -d "${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/wps.XXXX"); mkdir -p "$T/wps/metgrid" "$T/check"
cp -L "$WPS"/{namelist.wps,Vtable,ungrib.exe,metgrid.exe} "$WPS"/geo_em.d0*.nc "$T/wps/"
cp -L "$WPS/metgrid/METGRID.TBL" "$T/wps/metgrid/"
first=${HOURS%-*} last=${HOURS#*-}
c0="${CYCLE:0:4}-${CYCLE:4:2}-${CYCLE:6:2} ${CYCLE:8:2}:00"
for f in "$WPS"/met_em.d0*.nc; do
  t=$(basename "$f" | sed -E 's/met_em\.d0.\.(.*)_(..):00:00\.nc/\1 \2:00/')
  h=$(( ($(date -ud "$t" +%s) - $(date -ud "$c0" +%s)) / 3600 ))
  [ "$h" -ge "$first" ] && [ "$h" -le "$last" ] && cp -L "$f" "$T/check/"
done
tar czf "$OUT/wps.tar.gz" -C "$T" wps check; rm -rf "$T"
gcloud storage cp "$OUT/wps.tar.gz" "gs://$BK/src/$STAMP.wps.tar.gz" --project "$P" -q
ZONE=
cleanup() { [ -n "$ZONE" ] && gcloud compute instances delete "$VM" --zone "$ZONE" --project "$P" -q 2>/dev/null || true; }
trap cleanup EXIT
while read -r -u 3 M z DISK; do
  if gcloud compute instances create "$VM" --project "$P" --zone "$z" --machine-type "$M" \
       --provisioning-model=SPOT --instance-termination-action=DELETE --max-run-duration=3h \
       --image-family=debian-13 --image-project=debian-cloud --boot-disk-size=60GB --boot-disk-type=$DISK \
       --scopes=storage-rw --metadata-from-file=startup-script="$REPO/tools/gcp/remote_wps.sh" \
       --metadata="^@^bucket=$BK@src=src/$STAMP.wps.tar.gz@result=$RES@cycle=$CYCLE@hours=$HOURS@obj=$OBJ" \
       > "$OUT/create.log" 2>&1; then ZONE=$z MACHINE=$M; break; fi
  why=$(grep -oE "stockout|ZONE_RESOURCE_POOL_EXHAUSTED|Quota '[A-Z0-9_]*'|does not exist|not available|not supported" "$OUT/create.log" | head -1)
  echo "[gcp_wps] $M $z: ${why:-error}"
  [ -n "$why" ] || { cat "$OUT/create.log"; exit 1; }  # not a capacity or quota problem
done 3< <(vm_candidates "$MACHINE" 4 1 "" "$ZONES")
[ -n "$ZONE" ] || { echo "[gcp_wps] no machine and zone available now"; cat "$OUT/create.log"; exit 1; }
echo "[gcp_wps] $VM ($MACHINE spot, $ZONE): results to gs://$BK/$RES"
gone=0
while :; do
  sleep 60
  gcloud storage ls "gs://$BK/$RES/done" --project "$P" > /dev/null 2>&1 && break
  st=$(gcloud compute instances describe "$VM" --zone "$ZONE" --project "$P" --format='value(status)' 2>/dev/null || echo GONE)
  case $st in RUNNING|STAGING|PROVISIONING) gone=0 ;; *) gone=$((gone + 1)) ;; esac
  [ $gone -lt 3 ] || { echo "[gcp_wps] VM $st (preempted?)"; break; }
done
gcloud storage cp "gs://$BK/$RES/*" "$OUT/" --project "$P" -q 2>/dev/null || true
gcloud storage rm -r "gs://$BK/$RES" "gs://$BK/src/$STAMP.wps.tar.gz" --project "$P" -q 2>/dev/null || true
rm -f "$OUT/wps.tar.gz"
echo "[gcp_wps] results in $OUT"; cat "$OUT/wps.out" 2>/dev/null || echo "no wps.out"
grep -q "^STAGED " "$OUT/wps.out" 2>/dev/null
