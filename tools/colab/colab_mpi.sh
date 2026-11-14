#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Rank-layout invariance of the MPI build on a Colab runtime's CPU (tools/colab/mpi_body.sh in the Debian 13 chroot,
# as colab_suite.sh): suite cases run under mpirun with every layout of ranks (run.ranks, docs/MPI.md); each layout
# is compared with the WRF references at tolerance 0, as a run without MPI, and its output files byte for byte with
# those of the run with one rank.
#
#   CPPWRF_GCP_BUCKET=... tools/colab/colab_mpi.sh --cases NAME,... [--layouts 1,1x2,1x3,2x2,2x3,2x4,3x1,1x5,4x2]
#       [--mpi mpich|openmpi] [--gpu G4] [--jobs N] [--data OBJECT,...] [--out DIR] [--hours H] [--rev REV]
#       [--gpu-cases NAME,...] [--gpu-layouts 1,2x1,1x2,2x2] [--gpu-fp32 NAME,...] [--base REV]
#
# --gpu-cases: cuda_dp cases under mpirun with every --gpu-layouts layout on the runtime's GPU (tools/colab/mpi_gpu.sh).
# --gpu-fp32: cuda_dp cases run in FP32 by the --base revision's build (default origin/main) and by this revision's,
# as one process and as two ranks: every output file must be the base build's.
#       [--session NAME] [--keep]
#
# Data: the tars of the cases (wrf_suite.py --data-for) unless --data. --session / --keep as colab_suite.sh (a kept
# session keeps the build in /opt/mpi of its chroot). Results in OUT: suite_<layout>.log, sums_<layout>.txt,
# same.txt, mpi.txt, cmake.log; a summary is printed. Exit 0 when every layout's suite ends with exit 0 and every
# layout's files are those of one rank.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
: "${CPPWRF_GCP_BUCKET:?set CPPWRF_GCP_BUCKET}"
BK=$CPPWRF_GCP_BUCKET
GPU=G4 CASES= LAYOUTS=1,1x2,1x3,2x2,2x3,2x4,3x1,1x5,4x2 MPI=mpich DATA= OUT= JOBS=4 HOURS=4 REV=HEAD S= KEEP=0
NPROC= ANY_CPU=1 GCASES= GLAY=1,2x1,1x2,2x2 GSP= BASE=origin/main
while [ $# -gt 0 ]; do
  case $1 in
    --cases) CASES=$2; shift 2 ;; --layouts) LAYOUTS=$2; shift 2 ;; --mpi) MPI=$2; shift 2 ;; --gpu) GPU=$2; shift 2 ;;
    --jobs) JOBS=$2; shift 2 ;; --data) DATA=$2; shift 2 ;; --out) OUT=$2; shift 2 ;; --hours) HOURS=$2; shift 2 ;;
    --rev) REV=$2; shift 2 ;; --session) S=$2; shift 2 ;; --keep) KEEP=1; shift ;;
    --gpu-cases) GCASES=$2; shift 2 ;; --gpu-layouts) GLAY=$2; shift 2 ;;
    --gpu-fp32) GSP=$2; shift 2 ;; --base) BASE=$2; shift 2 ;;
    *) sed -n '3,15p' "$0"; exit 2 ;;
  esac
done
[ -n "$CASES$GCASES" ] || { echo "--cases or --gpu-cases"; exit 2; }
GPUJOB=0; [ -z "$GCASES" ] || GPUJOB=1
[ -n "$DATA" ] || DATA=$(python3 "$REPO/tests/verify/wrf_suite.py" --data-for "$(echo "$CASES,$GCASES" | sed 's/^,//; s/,$//')")
STAMP=$(date +%Y%m%d_%H%M%S)
OUT=${OUT:-${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/colab/mpi_$STAMP}
mkdir -p "$OUT"
[ "$REV" != HEAD ] || [ -z "$(git -C "$REPO" status --porcelain --untracked-files=no)" ] || echo "note: uncommitted changes are not sent (HEAD only)"
git -C "$REPO" archive --format=tar.gz -o "$OUT/src.tar.gz" "$REV"
git -C "$REPO" rev-parse --short "$REV" > "$OUT/commit.txt"
echo "[colab_mpi] commit $(cat "$OUT/commit.txt"), $GPU, $MPI, cases $CASES, layouts $LAYOUTS, data $DATA"
TAG=colab_mpi
. "$REPO/tools/colab/colab_lib.sh"
cj_start "$GPU"
[ "$GPUJOB" = 0 ] || cj_send "$REPO/tools/colab/mpi_gpu.sh" gpu.sh || exit 1
if [ -n "$GSP" ]; then
  git -C "$REPO" archive --format=tar.gz -o "$OUT/base_src.tar.gz" "$BASE" && cj_send "$OUT/base_src.tar.gz" in/base_src.tar.gz || exit 1
  rm -f "$OUT/base_src.tar.gz"
fi
cj_launch "$REPO/tools/colab/mpi_body.sh" "$DATA" "$CASES" "$LAYOUTS" "$JOBS" "$MPI" "$GPUJOB"
if [ "$GPUJOB" = 1 ]; then
  echo "import subprocess; subprocess.Popen(\"nohup bash /content/cpu/gpu.sh '$GCASES' '$GLAY' '$GSP' > /dev/null 2>&1 &\", shell=True); print('started')" |
    rx 60 | grep started > /dev/null || echo "[colab_mpi] the GPU side did not start"
fi
cj_wait "$HOURS" setup.out
cj_results
bad=0
grep . "$OUT/setup.out" 2> /dev/null && bad=1
cat "$OUT/mpi.txt" 2> /dev/null || true
for L in ${LAYOUTS//,/ }; do
  [ -n "$CASES" ] || break
  f=$OUT/suite_$L.log
  [ -f "$f" ] || { echo "layout $L: no log"; bad=1; continue; }
  echo "== layout $L: $(grep -E '^SUITE' "$f" | tail -1) $(tail -1 "$f")"
  grep -E '^(FAIL|SKIP) ' "$f" | cut -c1-400 || true
  grep -q '^exit 0' "$f" || bad=1
  [ "$L" = 1 ] || grep -q "^SAME $L " "$OUT/same.txt" 2> /dev/null || { grep "^DIFF $L:" "$OUT/same.txt" 2> /dev/null || echo "layout $L: not compared with one rank"; bad=1; }
done
grep '^SAME' "$OUT/same.txt" 2> /dev/null || true
if [ "$GPUJOB" = 1 ]; then
  cat "$OUT/gpu.txt" 2> /dev/null || echo "[colab_mpi] no gpu.txt"
  for L in ${GLAY//,/ }; do
    f=$OUT/gpu_suite_$L.log
    [ -f "$f" ] || { echo "GPU layout $L: no log"; bad=1; continue; }
    echo "== GPU FP64 layout $L: $(grep -E '^SUITE' "$f" | tail -1) $(tail -1 "$f")"
    grep -E '^(FAIL|SKIP) ' "$f" | cut -c1-400 || true
    grep -q '^exit 0' "$f" || bad=1
    [ "$L" = 1 ] || grep "^SAME $L " "$OUT/gpu_same.txt" 2> /dev/null || { grep "^DIFF $L:" "$OUT/gpu_same.txt" 2> /dev/null || echo "GPU layout $L: not compared with one rank"; bad=1; }
  done
  if [ -n "$GSP" ]; then
    cat "$OUT/gpu_fp32.txt" 2> /dev/null || { echo "[colab_mpi] no gpu_fp32.txt"; bad=1; }
    grep -q '^DIFF' "$OUT/gpu_fp32.txt" 2> /dev/null && bad=1
    grep -q '^SAME' "$OUT/gpu_fp32.txt" 2> /dev/null || bad=1
  fi
fi
exit $bad
