#!/usr/bin/env bash
# WRF-vs-CPPWRF namelist-variant matrix on a Google Compute Engine spot VM (Intel; see gcp_run.sh).
#
#   CPPWRF_GCP_PROJECT=... CPPWRF_GCP_BUCKET=... tools/gcp/gcp_matrix.sh VARIANTS.txt [--machine TYPE] [--data OBJ,...]
#
# VARIANTS.txt: one variant per line "NAME BASE_RUN WRF_EXE section.key=value ..." (BASE_RUN a run directory of
# the data tar's runs/, WRF_EXE a file of its wrf/: wrf_r8.exe or wrf_r8d8.exe). Each variant runs 3 minutes in
# WRF and CPPWRF (tools/gcp/matrix_variant.py); the result lines go to OUT/matrix.out. The data tar (default
# matrix_data.tar in the bucket) is made with tools/gcp/stage_matrix_data.sh. The source is the committed HEAD.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
: "${CPPWRF_GCP_PROJECT:?set CPPWRF_GCP_PROJECT}" "${CPPWRF_GCP_BUCKET:?set CPPWRF_GCP_BUCKET}"
P=$CPPWRF_GCP_PROJECT BK=$CPPWRF_GCP_BUCKET
VARS=$1; shift
MACHINE=c3-highcpu-44 DATA=matrix_data.tar OUT=
ZONES=us-central1-a,us-central1-b,us-central1-c,us-central1-f,us-east1-b,us-east1-c,us-east1-d,us-east4-a,us-east4-b,us-east4-c
while [ $# -gt 0 ]; do
  case $1 in --machine) MACHINE=$2; shift 2 ;; --data) DATA=$2; shift 2 ;; --out) OUT=$2; shift 2 ;; *) sed -n '2,10p' "$0"; exit 2 ;; esac
done
STAMP=$(date +%Y%m%d-%H%M%S)-$(printf %04x $RANDOM) VM=cppwrf-matrix-$STAMP RES=results/matrix-$STAMP
OUT=${OUT:-${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/gcp/matrix-$STAMP}; mkdir -p "$OUT"
git -C "$REPO" archive --format=tar.gz -o "$OUT/src.tar.gz" HEAD
gcloud storage cp "$OUT/src.tar.gz" "gs://$BK/src/$STAMP.tar.gz" --project "$P" -q
gcloud storage cp "$VARS" "gs://$BK/src/$STAMP.variants.txt" --project "$P" -q
ZONE=
cleanup() { [ -n "$ZONE" ] && gcloud compute instances delete "$VM" --zone "$ZONE" --project "$P" -q 2>/dev/null || true; }
trap cleanup EXIT
for z in ${ZONES//,/ }; do
  if gcloud compute instances create "$VM" --project "$P" --zone "$z" --machine-type "$MACHINE" \
       --provisioning-model=SPOT --instance-termination-action=DELETE --max-run-duration=3h \
       --image-family=debian-13 --image-project=debian-cloud --boot-disk-size=60GB --boot-disk-type=pd-balanced \
       --scopes=storage-rw --metadata-from-file=startup-script="$REPO/tools/gcp/remote_matrix.sh" \
       --metadata="^@^bucket=$BK@src=src/$STAMP.tar.gz@data=$DATA@result=$RES@variants=src/$STAMP.variants.txt" \
       > "$OUT/create.log" 2>&1; then ZONE=$z; break; fi
done
[ -n "$ZONE" ] || { echo "[gcp_matrix] no zone had $MACHINE"; cat "$OUT/create.log"; exit 1; }
echo "[gcp_matrix] $VM in $ZONE"
while :; do
  sleep 60
  gcloud storage cp "gs://$BK/$RES/matrix.out" "$OUT/" --project "$P" -q 2>/dev/null || true
  gcloud storage ls "gs://$BK/$RES/done" --project "$P" > /dev/null 2>&1 && break
  gcloud compute instances describe "$VM" --zone "$ZONE" --project "$P" > /dev/null 2>&1 || { echo "[gcp_matrix] VM gone (spot preemption?)"; break; }
done
gcloud storage cp "gs://$BK/$RES/*" "$OUT/" --project "$P" -q 2>/dev/null || true
gcloud storage rm -r "gs://$BK/$RES" "gs://$BK/src/$STAMP.tar.gz" "gs://$BK/src/$STAMP.variants.txt" --project "$P" -q 2>/dev/null || true
echo "[gcp_matrix] results in $OUT"; grep -v "^MATRIX_DONE" "$OUT/matrix.out" 2>/dev/null | cut -c1-300 || true
