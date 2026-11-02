#!/usr/bin/env bash
# The CPU FP64 suite of tools/gcp/gcp_run.sh on a Google Colab runtime with an Intel CPU (prepaid compute units
# instead of GCP VMs).
#
#   CPPWRF_GCP_BUCKET=... tools/colab/colab_suite.sh [--gpu L4|A100|T4] [--cases NAME,...] [--data OBJECT,...]
#       [--out DIR] [--jobs N] [--threads N] [--slow] [--hours H] [--rev REV] [--session NAME] [--keep]
#
# --gpu: the runtime type (default L4: Intel Cascade Lake, 12 vCPUs, 52 GB; the GPU is not used). G4 and the CPU-only
# runtime are AMD: AMD's glibc vector log / atan differ, so they cannot reproduce the Intel-made references; the job
# stops with NOT_INTEL on a runtime that is not Intel. --cases, --data, --slow as gcp_run.sh (the tars of
# tests/verify/vm_data.txt); --jobs suite jobs (default 6; with --cases at most one per case) with --threads each
# (default the runtime's vCPUs / jobs); --hours: the longest wait for the job (default 10; 20 with --slow); --rev: the
# commit to test (default HEAD; the tools are always this checkout's). --session: an existing runtime of an earlier
# --keep run (its Debian chroot and tars are reused); --keep: leave the runtime running (`colab stop -s NAME` later).
#
# On the runtime (tools/colab/remote_chroot_job.sh): a Debian 13 chroot with the GCP VMs' packages, the tars from Google
# Drive (cppwrf-drive:cppwrf-data, tools/colab/drive_sync.sh; a tar not there yet comes from the bucket with a
# read-only token, tools/colab/gcs_token.py, and is put on Drive), then tools/gcp/suite_body.sh as on the GCP VM.
# The runtime gets only the rclone section of cppwrf-drive and that token; both are deleted when the data are in
# place. Results: OUT (default ${CPPWRF_SCRATCH:-~/cppwrf_scratch}/colab/suite_<date>) with the suite's logs, job.log,
# versions.txt (packages, CPU, RAM) and data.txt (per tar: bytes, seconds, source). Exit status as gcp_run.sh's
# (tools/gcp/suite_status.sh): 0 only for a complete suite without FAIL or SKIP and passing tests.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
: "${CPPWRF_GCP_BUCKET:?set CPPWRF_GCP_BUCKET}"
BK=$CPPWRF_GCP_BUCKET
GPU=L4 CASES= DATA=suite_data.tar,ideal_data.tar OUT= JOBS=6 THREADS= SLOW= HOURS= REV=HEAD S= KEEP=0 NPROC=
while [ $# -gt 0 ]; do
  case $1 in
    --gpu) GPU=$2; shift 2 ;; --cases) CASES=$2; shift 2 ;; --data) DATA=$2; shift 2 ;; --out) OUT=$2; shift 2 ;;
    --jobs) JOBS=$2; shift 2 ;; --threads) THREADS=$2; shift 2 ;; --slow) SLOW=1; shift ;; --hours) HOURS=$2; shift 2 ;;
    --rev) REV=$2; shift 2 ;; --session) S=$2; shift 2 ;; --keep) KEEP=1; shift ;;
    *) sed -n '2,22p' "$0"; exit 2 ;;
  esac
done
case $GPU in L4|A100|T4) ;; *) echo "--gpu $GPU: only the Intel runtimes L4, A100, T4"; exit 2 ;; esac
HOURS=${HOURS:-$([ "$SLOW" = 1 ] && echo 20 || echo 10)}
if [ -n "$CASES" ]; then NC=$(tr ',' '\n' <<< "$CASES" | grep -c .); [ "$JOBS" -le "$NC" ] || JOBS=$NC; fi
STAMP=$(date +%Y%m%d_%H%M%S)
OUT=${OUT:-${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/colab/suite_$STAMP}
mkdir -p "$OUT"
[ "$REV" != HEAD ] || [ -z "$(git -C "$REPO" status --porcelain --untracked-files=no)" ] || echo "note: uncommitted changes are not sent (HEAD only)"
git -C "$REPO" archive --format=tar.gz -o "$OUT/src.tar.gz" "$REV"
git -C "$REPO" rev-parse --short "$REV" > "$OUT/commit.txt"
echo "[colab_suite] commit $(cat "$OUT/commit.txt"), $GPU, data $DATA${CASES:+, cases $CASES}"
TAG=colab_suite
. "$REPO/tools/colab/colab_lib.sh"
cj_start "$GPU"
THR=${THREADS:-$(( NPROC / JOBS > 0 ? NPROC / JOBS : 1 ))}
echo "[colab_suite] jobs $JOBS, threads $THR"
cj_launch "$REPO/tools/gcp/suite_body.sh" "$DATA" "$JOBS" "$THR" "$CASES" "${SLOW:-0}"
cj_wait "$HOURS" suite.log
cj_results
bash "$REPO/tools/gcp/suite_status.sh" "$OUT" "$CASES"
