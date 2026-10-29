#!/usr/bin/env bash
# Run wrf_suite.py's CPU FP64 cases on a Google Compute Engine spot VM and fetch the results.
#
#   CPPWRF_GCP_PROJECT=... CPPWRF_GCP_BUCKET=... tools/gcp/gcp_run.sh [--machine TYPE] [--zones Z1,Z2,...]
#       [--jobs N] [--threads N] [--cases NAME,...] [--data OBJECT] [--out DIR]
#
# The bucket holds the reference data once (tools/colab/stage_data.py's tar, uploaded with
# `gcloud storage cp`; --data names them, comma-separated, default suite_data.tar,ideal_data.tar); the source is the committed HEAD.
# Before the suite the VM runs tests/tools/test_wrf2cppwrf.py (the translator; OUT/unit.log); after the full suite
# (no --cases) tests/verify/run_wrf_files_test.sh and ideal_files_test.py (OUT/files.log, OUT/idealfiles.log).
# The VM is a spot VM (Debian 13: the same gcc and glibc as the reference workstation) that powers
# itself off at the end and is deleted here; --max-run-duration deletes it in any case. The zones are tried
# in order until one has the machine (spot capacity varies). Use Intel machines (c3, n2, n1): glibc's
# libmvec log differs on AMD, so AMD runs are not bit-identical to Intel-made WRF references.
# Exit status 0 only for a complete suite without FAIL or SKIP (wrf_suite.py's SUITE line and exit 0), passing
# translator tests and comparator self-test (tests/verify/selftest.py), files tests without FAIL or SKIP (full
# suite), and no data or build failure on the VM; a
# VM that ends early is INCOMPLETE (status 1).
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
: "${CPPWRF_GCP_PROJECT:?set CPPWRF_GCP_PROJECT}" "${CPPWRF_GCP_BUCKET:?set CPPWRF_GCP_BUCKET}"
P=$CPPWRF_GCP_PROJECT BK=$CPPWRF_GCP_BUCKET
MACHINE=c3-highcpu-176 ZONES=us-central1-a,us-central1-b,us-central1-c,us-central1-f,us-east1-b,us-east1-c,us-east1-d,us-east4-a,us-east4-b,us-east4-c,europe-west4-a,europe-west4-b,europe-west4-c JOBS=11 THREADS=16 CASES= DATA=suite_data.tar,ideal_data.tar OUT= MAXH=8
while [ $# -gt 0 ]; do
  case $1 in
    --machine) MACHINE=$2; shift 2 ;; --zones) ZONES=$2; shift 2 ;; --jobs) JOBS=$2; shift 2 ;;
    --threads) THREADS=$2; shift 2 ;; --cases) CASES=$2; shift 2 ;; --data) DATA=$2; shift 2 ;;
    --out) OUT=$2; shift 2 ;; *) sed -n '2,12p' "$0"; exit 2 ;;
  esac
done
STAMP=$(date +%Y%m%d-%H%M%S)-$(printf %04x $RANDOM) VM=cppwrf-suite-$STAMP RES=results/$STAMP
OUT=${OUT:-${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/gcp/$STAMP}
mkdir -p "$OUT"
[ -z "$(git -C "$REPO" status --porcelain --untracked-files=no)" ] || echo "note: uncommitted changes are not sent (HEAD only)"
git -C "$REPO" archive --format=tar.gz -o "$OUT/src.tar.gz" HEAD
gcloud storage cp "$OUT/src.tar.gz" "gs://$BK/src/$STAMP.tar.gz" --project "$P" -q
cleanup() { gcloud compute instances delete "$VM" --zone "$ZONE" --project "$P" -q 2>/dev/null || true; }
trap cleanup EXIT
ZONE=
for z in ${ZONES//,/ }; do
  if gcloud compute instances create "$VM" --project "$P" --zone "$z" --machine-type "$MACHINE" \
       --provisioning-model=SPOT --instance-termination-action=DELETE --max-run-duration="${MAXH}h" \
       --image-family=debian-13 --image-project=debian-cloud --boot-disk-size=200GB --boot-disk-type=pd-balanced \
       --scopes=storage-rw --metadata-from-file=startup-script="$REPO/tools/gcp/remote_suite.sh" \
       --metadata="^@^bucket=$BK@src=src/$STAMP.tar.gz@data=$DATA@result=$RES@jobs=$JOBS@threads=$THREADS@cases=$CASES" \
       > "$OUT/create.log" 2>&1; then ZONE=$z; break; fi
  why=$(grep -oE "stockout|ZONE_RESOURCE_POOL_EXHAUSTED|Quota '[A-Z0-9_]*'|does not exist|not available|not supported" "$OUT/create.log" | head -1)
  echo "[gcp_run] $z: ${why:-error}"
  [ -n "$why" ] || { cat "$OUT/create.log"; exit 1; }  # not a capacity or quota problem
done
[ -n "$ZONE" ] || { echo "[gcp_run] no zone has $MACHINE now"; cat "$OUT/create.log"; exit 1; }
echo "[gcp_run] $VM ($MACHINE spot, $ZONE): results to gs://$BK/$RES"
gone=0
while :; do
  sleep 60
  if gcloud storage ls "gs://$BK/$RES/done" --project "$P" > /dev/null 2>&1; then break; fi
  st=$(gcloud compute instances describe "$VM" --zone "$ZONE" --project "$P" --format='value(status)' 2>/dev/null || echo GONE)
  case $st in RUNNING|STAGING|PROVISIONING) gone=0 ;; *) gone=$((gone + 1)) ;; esac
  # three checks in a row: one failed describe call (network) must not end a running job
  [ $gone -lt 3 ] || { echo "[gcp_run] VM $st (preempted?)"; break; }
done
gcloud storage cp "gs://$BK/$RES/*" "$OUT/" --project "$P" -q || true
rm -f "$OUT/src.tar.gz"; gcloud storage rm "gs://$BK/src/$STAMP.tar.gz" --project "$P" -q || true
echo "[gcp_run] results in $OUT"
for k in PASS FAIL SKIP; do echo "$k $(grep -c "^$k " "$OUT/suite.log" 2>/dev/null || true)"; done
grep -E "^(FAIL|SKIP|ERROR) " "$OUT/suite.log" 2>/dev/null || true
bad=0
cat "$OUT/setup.out" 2>/dev/null | grep . && bad=1
grep -q "^exit 0$" "$OUT/unit.log" 2>/dev/null || { echo "[gcp_run] translator tests: $(tail -1 "$OUT/unit.log" 2>/dev/null || echo 'no unit.log')"; bad=1; }
grep -q "^exit 0$" "$OUT/selftest.log" 2>/dev/null || { echo "[gcp_run] comparator self-test: $(grep -E "WRONG|SELFTEST" "$OUT/selftest.log" 2>/dev/null || echo 'no selftest.log')"; bad=1; }
if ! grep -q "^SUITE " "$OUT/suite.log" 2>/dev/null; then echo "[gcp_run] INCOMPLETE: the suite did not end (no SUITE line)"; bad=1
else grep "^SUITE " "$OUT/suite.log"; grep -q "^exit 0$" "$OUT/suite.log" || bad=1; fi
if [ -z "$CASES" ]; then  # the full suite: the output files tests too (remote_suite.sh)
  for f in files idealfiles; do
    if [ ! -f "$OUT/$f.log" ]; then echo "[gcp_run] INCOMPLETE: no $f.log"; bad=1; continue; fi
    echo "$f: PASS $(grep -c '^PASS ' "$OUT/$f.log") FAIL $(grep -c '^FAIL ' "$OUT/$f.log") SKIP $(grep -c '^SKIP ' "$OUT/$f.log")"
    grep -E "^(FAIL|SKIP) " "$OUT/$f.log" | cut -c1-200 || true  # none: grep fails (set -e, pipefail)
    grep -E "^exit [1-9]" "$OUT/$f.log" && bad=1
  done
fi
exit $bad
