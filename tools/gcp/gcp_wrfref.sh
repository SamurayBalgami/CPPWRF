#!/usr/bin/env bash
# WRF reference runs for the suite on a Google Compute Engine spot VM (Intel; see gcp_run.sh).
#
#   CPPWRF_GCP_PROJECT=... CPPWRF_GCP_BUCKET=... tools/gcp/gcp_wrfref.sh REFS.txt [--machine TYPE] [--data OBJ,...]
#       [--files NAME,...] [--cases NAME,...] [--runs DIR] [--out DIR] [--check-only]
#
# REFS.txt: one reference per line "NAME BASE_RUN WRF_EXE section.key=value ..." as gcp_matrix.sh's variants
# (time_control.run_seconds / history_interval_s set the run). The VM makes run_NAME and run_NAME/steps (3 steps,
# history every step; tools/gcp/wrf_ref_run.py), rst_NAME and rsm_NAME (tests/verify/make_restart_refs.sh,
# make_resumed_refs.sh) and, with --files / --cases, runs run_wrf_files_test.sh / wrf_suite.py (exact case names)
# of the committed HEAD against them. The references are unpacked into --runs (default $CPPWRF_WRF_RUNS, an
# existing run_NAME stops the script), the links to the VM's inputs pointing to the files of DIR/BASE_RUN and
# wrf.exe to $WRF_R8 / $WRF_R8D8 (default ~/wrf_real/WRF-4.8.0-r8[d8]/main/wrf.exe). --check-only: only the tests
# on the VM (the references are neither uploaded nor unpacked).
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
: "${CPPWRF_GCP_PROJECT:?set CPPWRF_GCP_PROJECT}" "${CPPWRF_GCP_BUCKET:?set CPPWRF_GCP_BUCKET}"
P=$CPPWRF_GCP_PROJECT BK=$CPPWRF_GCP_BUCKET
REFS=$1; shift
MACHINE=c3-highcpu-22 DATA=matrix_data.tar OUT= FILES= CASES= FETCH=1 RUNS=${CPPWRF_WRF_RUNS:-$HOME/wrf_real480}
WRF_R8=${WRF_R8:-$HOME/wrf_real/WRF-4.8.0-r8/main/wrf.exe} WRF_R8D8=${WRF_R8D8:-$HOME/wrf_real/WRF-4.8.0-r8d8/main/wrf.exe}
ZONES=us-central1-a,us-central1-b,us-central1-c,us-central1-f,us-east1-b,us-east1-c,us-east1-d,us-east4-a,us-east4-b,us-east4-c
while [ $# -gt 0 ]; do
  case $1 in
    --machine) MACHINE=$2; shift 2 ;; --data) DATA=$2; shift 2 ;; --out) OUT=$2; shift 2 ;;
    --files) FILES=$2; shift 2 ;; --cases) CASES=$2; shift 2 ;; --runs) RUNS=$2; shift 2 ;;
    --check-only) FETCH=0; shift ;; *) sed -n '2,15p' "$0"; exit 2 ;;
  esac
done
names=$(grep -v '^\s*#' "$REFS" | awk 'NF {print $1}')
[ $FETCH = 0 ] || for n in $names; do [ ! -e "$RUNS/run_$n" ] || { echo "[gcp_wrfref] $RUNS/run_$n exists"; exit 1; }; done
STAMP=$(date +%Y%m%d-%H%M%S)-$(printf %04x $RANDOM) VM=cppwrf-wrfref-$STAMP RES=results/wrfref-$STAMP
OUT=${OUT:-${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/gcp/wrfref-$STAMP}; mkdir -p "$OUT"
[ -z "$(git -C "$REPO" status --porcelain --untracked-files=no)" ] || echo "note: uncommitted changes are not sent (HEAD only)"
git -C "$REPO" archive --format=tar.gz -o "$OUT/src.tar.gz" HEAD
gcloud storage cp "$OUT/src.tar.gz" "gs://$BK/src/$STAMP.tar.gz" --project "$P" -q
gcloud storage cp "$REFS" "gs://$BK/src/$STAMP.refs.txt" --project "$P" -q
ZONE=
cleanup() { [ -n "$ZONE" ] && gcloud compute instances delete "$VM" --zone "$ZONE" --project "$P" -q 2>/dev/null || true; }
trap cleanup EXIT
for z in ${ZONES//,/ }; do
  if gcloud compute instances create "$VM" --project "$P" --zone "$z" --machine-type "$MACHINE" \
       --provisioning-model=SPOT --instance-termination-action=DELETE --max-run-duration=4h \
       --image-family=debian-13 --image-project=debian-cloud --boot-disk-size=60GB --boot-disk-type=pd-balanced \
       --scopes=storage-rw --metadata-from-file=startup-script="$REPO/tools/gcp/remote_wrfref.sh" \
       --metadata="^@^bucket=$BK@src=src/$STAMP.tar.gz@data=$DATA@result=$RES@refs=src/$STAMP.refs.txt@files=$FILES@cases=$CASES@fetch=$FETCH" \
       > "$OUT/create.log" 2>&1; then ZONE=$z; break; fi
  why=$(grep -oE "stockout|ZONE_RESOURCE_POOL_EXHAUSTED|Quota '[A-Z_]*'|does not exist|not available|not supported" "$OUT/create.log" | head -1)
  echo "[gcp_wrfref] $z: ${why:-error}"
  [ -n "$why" ] || { cat "$OUT/create.log"; exit 1; }  # not a capacity or quota problem
done
[ -n "$ZONE" ] || { echo "[gcp_wrfref] no zone has $MACHINE now"; cat "$OUT/create.log"; exit 1; }
echo "[gcp_wrfref] $VM ($MACHINE spot, $ZONE): results to gs://$BK/$RES"
while :; do
  sleep 60
  gcloud storage cp "gs://$BK/$RES/refs.out" "$OUT/" --project "$P" -q 2>/dev/null || true
  gcloud storage ls "gs://$BK/$RES/done" --project "$P" > /dev/null 2>&1 && break
  st=$(gcloud compute instances describe "$VM" --zone "$ZONE" --project "$P" --format='value(status)' 2>/dev/null || echo GONE)
  case $st in RUNNING|STAGING|PROVISIONING) ;; *) echo "[gcp_wrfref] VM $st (preempted?)"; break ;; esac
done
for f in job.log refs.out cmake.log files.log suite.log; do
  gcloud storage cp "gs://$BK/$RES/$f" "$OUT/" --project "$P" -q 2>/dev/null || true
done
if [ $FETCH = 1 ] && gcloud storage ls "gs://$BK/$RES/refs.tar" --project "$P" > /dev/null 2>&1; then
  gcloud storage cat "gs://$BK/$RES/refs.tar" --project "$P" | tar xf - -C "$RUNS"
  # the VM's links -> the same local files (base run inputs, the other reference directories, wrf.exe)
  for n in $names; do
    for d in "$RUNS/run_$n" "$RUNS/rst_$n" "$RUNS/rsm_$n"; do
      [ -d "$d" ] || continue
      find "$d" -type l | while read -r l; do
        t=$(readlink "$l")
        case $t in
          /opt/job/data/runs/*) nt=$(readlink -f "$RUNS/${t#/opt/job/data/runs/}") ;;
          /opt/job/wrf_runs/*) nt=$RUNS/${t#/opt/job/wrf_runs/} ;;
          */wrf_r8d8.exe) nt=$WRF_R8D8 ;;
          */wrf_r8.exe) nt=$WRF_R8 ;;
          *) echo "[gcp_wrfref] link not mapped: $l -> $t"; continue ;;
        esac
        [ -e "$nt" ] || echo "[gcp_wrfref] missing locally: $l -> $nt"
        ln -sfn "$nt" "$l"
      done
    done
  done
  echo "[gcp_wrfref] references in $RUNS: $(for n in $names; do ls -d "$RUNS"/{run,rst,rsm}_"$n" 2>/dev/null; done | tr '\n' ' ')"
fi
gcloud storage rm -r "gs://$BK/$RES" "gs://$BK/src/$STAMP.tar.gz" "gs://$BK/src/$STAMP.refs.txt" --project "$P" -q 2>/dev/null || true
rm -f "$OUT/src.tar.gz"
echo "[gcp_wrfref] results in $OUT"; cat "$OUT/refs.out" 2>/dev/null || true
for f in files.log suite.log; do [ -f "$OUT/$f" ] && { echo "== $f"; grep -E "^(PASS|FAIL|SKIP)|compared" "$OUT/$f" || tail -5 "$OUT/$f"; }; done
exit 0
