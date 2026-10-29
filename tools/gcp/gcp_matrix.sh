#!/usr/bin/env bash
# WRF-vs-CPPWRF namelist-variant matrix on a Google Compute Engine spot VM (Intel; see gcp_run.sh).
#
#   CPPWRF_GCP_PROJECT=... CPPWRF_GCP_BUCKET=... tools/gcp/gcp_matrix.sh VARIANTS.txt [--machine TYPE] [--zones Z,...]
#       [--data OBJ,...]   (machine: tools/gcp/gcp_pick.py's cheapest Intel spot VM with a vCPU per variant line,
#       at least 8, and its zone now; --machine fixes the type, --zones limits the zones)
#
# VARIANTS.txt: one variant per line "NAME BASE_RUN WRF_EXE section.key=value ..." (BASE_RUN a run directory of
# the data tar's runs/, WRF_EXE a file of its wrf/: wrf_r8.exe or wrf_r8d8.exe). Each variant runs 3 minutes in
# WRF and CPPWRF (tools/gcp/matrix_variant.py); the result lines go to OUT/matrix.out. The data tar (default
# matrix_data.tar in the bucket) is made with tools/gcp/stage_matrix_data.sh. The source is the committed HEAD.
# Exit status 1 when the matrix is incomplete: a variant without a result line, a VM that ended early, a failed
# build or data download, a variant whose comparison failed (HARNESS-ERROR); differences are results (status 0).
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
. "$REPO/tools/gcp/pick_lib.sh"
: "${CPPWRF_GCP_PROJECT:?set CPPWRF_GCP_PROJECT}" "${CPPWRF_GCP_BUCKET:?set CPPWRF_GCP_BUCKET}"
P=$CPPWRF_GCP_PROJECT BK=$CPPWRF_GCP_BUCKET
VARS=$1; shift
MACHINE= DATA=matrix_data.tar OUT=
ZONES=
while [ $# -gt 0 ]; do
  case $1 in --machine) MACHINE=$2; shift 2 ;; --zones) ZONES=$2; shift 2 ;; --data) DATA=$2; shift 2 ;; --out) OUT=$2; shift 2 ;; *) sed -n '2,12p' "$0"; exit 2 ;; esac
done
# one vCPU per variant line (each runs WRF, then CPPWRF, single-threaded), at least 8 for the build
N=$(grep -cvE '^[[:space:]]*(#|$)' "$VARS" || true) NEED=$(( N > 8 ? N : 8 ))
STAMP=$(date +%Y%m%d-%H%M%S)-$(printf %04x $RANDOM) VM=cppwrf-matrix-$STAMP RES=results/matrix-$STAMP
OUT=${OUT:-${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/gcp/matrix-$STAMP}; mkdir -p "$OUT"
git -C "$REPO" archive --format=tar.gz -o "$OUT/src.tar.gz" HEAD
gcloud storage cp "$OUT/src.tar.gz" "gs://$BK/src/$STAMP.tar.gz" --project "$P" -q
gcloud storage cp "$VARS" "gs://$BK/src/$STAMP.variants.txt" --project "$P" -q
ZONE=
cleanup() { [ -n "$ZONE" ] && gcloud compute instances delete "$VM" --zone "$ZONE" --project "$P" -q 2>/dev/null || true; }
trap cleanup EXIT
while read -r -u 3 M z DISK; do
  if gcloud compute instances create "$VM" --project "$P" --zone "$z" --machine-type "$M" \
       --provisioning-model=SPOT --instance-termination-action=DELETE --max-run-duration=3h \
       --image-family=debian-13 --image-project=debian-cloud --boot-disk-size=60GB --boot-disk-type=$DISK \
       --scopes=storage-rw --metadata-from-file=startup-script="$REPO/tools/gcp/remote_matrix.sh" \
       --metadata="^@^bucket=$BK@src=src/$STAMP.tar.gz@data=$DATA@result=$RES@variants=src/$STAMP.variants.txt" \
       > "$OUT/create.log" 2>&1; then ZONE=$z MACHINE=$M; break; fi
  why=$(grep -oE "stockout|ZONE_RESOURCE_POOL_EXHAUSTED|Quota '[A-Z0-9_]*'|does not exist|not available|not supported" "$OUT/create.log" | head -1)
  echo "[gcp_matrix] $M $z: ${why:-error}"
  [ -n "$why" ] || { cat "$OUT/create.log"; exit 1; }  # not a capacity or quota problem
done 3< <(vm_candidates "$MACHINE" $NEED 1 "$DATA" "$ZONES")
[ -n "$ZONE" ] || { echo "[gcp_matrix] no machine and zone available now"; cat "$OUT/create.log"; exit 1; }
echo "[gcp_matrix] $VM ($MACHINE spot, $ZONE) for $N variants"
gone=0
while :; do
  sleep 60
  gcloud storage cp "gs://$BK/$RES/matrix.out" "$OUT/" --project "$P" -q 2>/dev/null || true
  gcloud storage ls "gs://$BK/$RES/done" --project "$P" > /dev/null 2>&1 && break
  st=$(gcloud compute instances describe "$VM" --zone "$ZONE" --project "$P" --format='value(status)' 2>/dev/null || echo GONE)
  case $st in RUNNING|STAGING|PROVISIONING) gone=0 ;; *) gone=$((gone + 1)) ;; esac
  # three checks in a row: one failed describe call (network) must not end a running job
  [ $gone -lt 3 ] || { echo "[gcp_matrix] VM $st (spot preemption?)"; break; }
done
gcloud storage cp "gs://$BK/$RES/*" "$OUT/" --project "$P" -q 2>/dev/null || true
gcloud storage rm -r "gs://$BK/$RES" "gs://$BK/src/$STAMP.tar.gz" "gs://$BK/src/$STAMP.variants.txt" --project "$P" -q 2>/dev/null || true
echo "[gcp_matrix] results in $OUT"; grep -v "^MATRIX_DONE" "$OUT/matrix.out" 2>/dev/null | cut -c1-300 || true
# completeness: every variant has a result line, the job ended, the build and the data were there
bad=0
grep -q "^MATRIX_DONE" "$OUT/matrix.out" 2>/dev/null || { echo "[gcp_matrix] INCOMPLETE: the job did not end (no MATRIX_DONE)"; bad=1; }
grep -E "^MATRIX_(BUILD|DATA)_FAILED|^[^ ]+ HARNESS-ERROR" "$OUT/matrix.out" 2>/dev/null && bad=1
missing=$(grep -v '^\s*#' "$VARS" | awk 'NF {print $1}' | while read -r n; do
  grep -q "^$n " "$OUT/matrix.out" 2>/dev/null || echo "$n"; done | tr '\n' ' ' || true)
[ -z "$missing" ] || { echo "[gcp_matrix] INCOMPLETE: no result for $missing"; bad=1; }
exit $bad
