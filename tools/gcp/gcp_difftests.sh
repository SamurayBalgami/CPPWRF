#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The differential tests (ctest -L diff) and the other ctest tests on a Google Compute Engine spot VM, alone
# (the full suite of gcp_run.sh runs them too), optionally with deliberate faults that the tests must catch.
#
#   CPPWRF_GCP_PROJECT=... CPPWRF_GCP_BUCKET=... tools/gcp/gcp_difftests.sh [--faults DIR] [--machine TYPE]
#       [--zones Z,...] [--vcpus N] [--tests REGEX] [--data OBJ,...] [--out DIR]
#
# The source is the committed HEAD; --data (default diff_data.tar: tools/gcp/stage_diff_data.sh) names the data
# tars; --tests limits the ctest runs to the tests matching REGEX. --faults DIR: its *.patch files (patch -p1 on
# the source, a line "# tests: REGEX" naming the tests that must fail), each applied alone after the clean run
# (tools/gcp/run_diff_tests.sh). Machine: gcp_pick.py's
# cheapest Intel spot VM with --vcpus vCPUs (default 16: each fault rebuilds; --machine fixes the type, --zones
# limits the zones); remote_difftests.sh powers it off at the end and it is deleted here (--max-run-duration in
# any case).
# Exit status 0 only when diff.log and ctest.log end with "exit 0", every fault is DETECTED and no input or build
# failed on the VM; a VM that ends early is INCOMPLETE (status 1).
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
. "$REPO/tools/gcp/pick_lib.sh"
: "${CPPWRF_GCP_PROJECT:?set CPPWRF_GCP_PROJECT}" "${CPPWRF_GCP_BUCKET:?set CPPWRF_GCP_BUCKET}"
P=$CPPWRF_GCP_PROJECT BK=$CPPWRF_GCP_BUCKET
MACHINE= ZONES= DATA=diff_data.tar OUT= FAULTS= VCPUS=16 TESTS=
while [ $# -gt 0 ]; do
  case $1 in
    --machine) MACHINE=$2; shift 2 ;; --zones) ZONES=$2; shift 2 ;; --data) DATA=$2; shift 2 ;;
    --out) OUT=$2; shift 2 ;; --faults) FAULTS=$(cd "$2" && pwd); shift 2 ;; --vcpus) VCPUS=$2; shift 2 ;; --tests) TESTS=$2; shift 2 ;;
    *) sed -n "3,16p" "$0"; exit 2 ;;
  esac
done
STAMP=$(date +%Y%m%d-%H%M%S)-$(printf %04x $RANDOM) VM=cppwrf-diff-$STAMP RES=results/diff-$STAMP
OUT=${OUT:-${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/gcp/diff-$STAMP}; mkdir -p "$OUT"
[ -z "$(git -C "$REPO" status --porcelain --untracked-files=no)" ] || echo "note: uncommitted changes are not sent (HEAD only)"
git -C "$REPO" archive --format=tar.gz -o "$OUT/src.tar.gz" HEAD
gcloud storage cp "$OUT/src.tar.gz" "gs://$BK/src/$STAMP.tar.gz" --project "$P" -q
FOBJ=
if [ -n "$FAULTS" ]; then
  tar cf "$OUT/faults.tar" -C "$FAULTS" $(cd "$FAULTS" && ls *.patch)
  FOBJ=src/$STAMP.faults.tar; gcloud storage cp "$OUT/faults.tar" "gs://$BK/$FOBJ" --project "$P" -q
fi
ZONE=
cleanup() {  # also on an early exit: the VM and this job's bucket objects
  [ -n "$ZONE" ] && gcloud compute instances delete "$VM" --zone "$ZONE" --project "$P" -q 2>/dev/null
  gcloud storage rm -r "gs://$BK/$RES" "gs://$BK/src/$STAMP.tar.gz" ${FOBJ:+"gs://$BK/$FOBJ"} --project "$P" -q 2>/dev/null
  rm -f "$OUT/src.tar.gz" "$OUT/faults.tar"; true
}
trap cleanup EXIT
while read -r -u 3 M z DISK; do
  if gcloud compute instances create "$VM" --project "$P" --zone "$z" --machine-type "$M" \
       --provisioning-model=SPOT --instance-termination-action=DELETE --max-run-duration=6h \
       --image-family=debian-13 --image-project=debian-cloud --boot-disk-size=60GB --boot-disk-type=$DISK \
       --scopes=storage-rw --metadata-from-file=startup-script="$REPO/tools/gcp/remote_difftests.sh" \
       --metadata="^@^bucket=$BK@src=src/$STAMP.tar.gz@data=$DATA@faults=$FOBJ@result=$RES@tests=$TESTS" \
       > "$OUT/create.log" 2>&1; then ZONE=$z MACHINE=$M; break; fi
  why=$(grep -oE "stockout|ZONE_RESOURCE_POOL_EXHAUSTED|Quota '[A-Z0-9_]*'|does not exist|not available|not supported" "$OUT/create.log" | head -1)
  echo "[gcp_difftests] $M $z: ${why:-error}"
  [ -n "$why" ] || { cat "$OUT/create.log"; exit 1; }  # not a capacity or quota problem
done 3< <(vm_candidates "$MACHINE" "$VCPUS" 3 "$DATA" "$ZONES")
[ -n "$ZONE" ] || { echo "[gcp_difftests] no machine and zone available now"; cat "$OUT/create.log"; exit 1; }
echo "[gcp_difftests] $VM ($MACHINE spot, $ZONE): results to gs://$BK/$RES"
gone=0 bad=0
while :; do
  sleep 60
  gcloud storage ls "gs://$BK/$RES/done" --project "$P" > /dev/null 2>&1 && break
  st=$(gcloud compute instances describe "$VM" --zone "$ZONE" --project "$P" --format='value(status)' 2>/dev/null || echo GONE)
  case $st in RUNNING|STAGING|PROVISIONING) gone=0 ;; *) gone=$((gone + 1)) ;; esac
  # three checks in a row: one failed describe call (network) must not end a running job
  [ $gone -lt 3 ] || { echo "[gcp_difftests] VM $st (preempted?): INCOMPLETE"; bad=1; break; }
done
gcloud storage cp "gs://$BK/$RES/*" "$OUT/" --project "$P" -q 2>/dev/null || true
echo "[gcp_difftests] results in $OUT"
grep . "$OUT/setup.out" 2>/dev/null && bad=1
for f in diff ctest; do
  [ -f "$OUT/$f.log" ] || { echo "[gcp_difftests] INCOMPLETE: no $f.log"; bad=1; continue; }
  echo "== $f.log: $(grep -E 'tests passed' "$OUT/$f.log" || true)"
  grep -E '^\s*[0-9]+ - .*\((Failed|Timeout|SEGFAULT|Subprocess|Skipped|Not Run)|^SKIPPED ' "$OUT/$f.log" || true
  grep -qx "exit 0" "$OUT/$f.log" || bad=1
done
if [ -n "$FAULTS" ]; then
  cat "$OUT/faults.log" 2>/dev/null || { echo "[gcp_difftests] no faults.log"; bad=1; }
  n=$(ls "$FAULTS"/*.patch | wc -l)
  [ "$(grep -c ' DETECTED ' "$OUT/faults.log" 2>/dev/null || true)" = "$n" ] || bad=1
fi
exit $bad
