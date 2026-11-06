#!/usr/bin/env bash
# The CPU FP64 suite of tools/gcp/gcp_run.sh on a Google Colab runtime (prepaid compute units instead of GCP VMs).
#
#   CPPWRF_GCP_BUCKET=... tools/colab/colab_suite.sh [--gpu G4|CPU|L4|A100|T4] [--any-cpu] [--cases NAME,...]
#       [--data OBJECT,... | --auto-data] [--out DIR] [--jobs N] [--threads N] [--slow] [--hours H] [--rev REV]
#       [--session NAME] [--keep] [--diff-only [--tests REGEX] [--faults DIR]]
#
# The full suite runs on the G4 runtime's AMD CPU (AMD EPYC 9B45, 48 vCPUs, 176 GB; its GPU is not used):
#   tools/colab/colab_suite.sh --gpu G4 --any-cpu --jobs 12 --threads 4      (about 42 min; L4: 2 h 45 min)
# --gpu: the runtime type (default L4: Intel Cascade Lake, 12 vCPUs, 52 GB). G4 and the CPU-only runtime (--gpu CPU)
# are AMD and need --any-cpu: glibc libmvec's vector log / log10 differ there (rcpps), and the job stops with
# NOT_INTEL on a runtime that is not Intel unless --any-cpu. The suite gives Intel's bits on any x86-64 CPU: the
# model computes them as on Intel (CPPWRF_PORTABLE_LIBMVEC, src/core/libmvec_x86.hpp), and so does WRF's side of the
# differential tests (tests/phys/libmvec_intel.cpp, linked by tests/phys/wrf_fortran.sh). L4 (Intel) stays possible.
# --cases, --data, --slow as gcp_run.sh (the tars of
# tests/verify/vm_data.txt). Without --data (or with --auto-data) the tars come from the case registry: with --cases
# only those the cases read (DATA in tests/verify/cases/*.py: wrf_suite.py --data-for), else every tar of vm_data.txt
# (the full suite: the output files, ideal files and differential tests too) and with --slow the slow cases' tars.
# --jobs suite jobs (default 4: 6 ran out of the L4's 52 GB next to the differential
# tests' build; with --cases at most one per case) with --threads each
# (default the runtime's vCPUs / jobs); --hours: the longest wait for the job (default 10; 20 with --slow); --rev: the
# commit to test (default HEAD; the tools are always this checkout's). --session: an existing runtime of an earlier
# --keep run (its Debian chroot and tars are reused); --keep: leave the runtime running (`colab stop -s NAME` later).
# --diff-only: only the differential tests and the other ctest tests (tools/colab/diff_body.sh: run_diff_tests.sh, as
# the full suite runs them; data: the *diff*.tar of vm_data.txt unless --data), --tests: a ctest -R regex limiting
# them; --faults DIR: its *.patch files (patch -p1, a line "# tests: REGEX" naming the tests that must fail), each
# applied alone after the clean run (faults.log: DETECTED / MISSED). Exit 0 when diff.log and ctest.log end with
# "exit 0", setup.out is empty and every fault was DETECTED.
#
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
: "${CPPWRF_GCP_BUCKET:?set CPPWRF_GCP_BUCKET}"
BK=$CPPWRF_GCP_BUCKET
GPU=L4 CASES= DATA= OUT= JOBS=4 THREADS= SLOW= HOURS= REV=HEAD S= KEEP=0 NPROC= ANY_CPU=0 DIFF=0 TESTS= FAULTS=
while [ $# -gt 0 ]; do
  case $1 in
    --gpu) GPU=$2; shift 2 ;; --cases) CASES=$2; shift 2 ;; --data) DATA=$2; shift 2 ;; --auto-data) DATA=; shift ;;
    --out) OUT=$2; shift 2 ;;
    --jobs) JOBS=$2; shift 2 ;; --threads) THREADS=$2; shift 2 ;; --slow) SLOW=1; shift ;; --hours) HOURS=$2; shift 2 ;;
    --rev) REV=$2; shift 2 ;; --session) S=$2; shift 2 ;; --keep) KEEP=1; shift ;; --any-cpu) ANY_CPU=1; shift ;;
    --diff-only) DIFF=1; shift ;; --tests) TESTS=$2; shift 2 ;; --faults) FAULTS=$(cd "$2" && pwd); shift 2 ;;
    *) sed -n '2,29p' "$0"; exit 2 ;;
  esac
done
case $GPU in L4|A100|T4) ;; G4|CPU) [ "${ANY_CPU:-0}" = 1 ] || { echo "--gpu $GPU: an AMD runtime, needs --any-cpu"; exit 2; } ;;
  *) echo "--gpu $GPU: L4, A100, T4 (Intel), G4 or CPU (AMD, with --any-cpu)"; exit 2 ;; esac
HOURS=${HOURS:-$([ "$SLOW" = 1 ] && echo 20 || echo 10)}
[ "$DIFF" = 1 ] || [ -z "$TESTS$FAULTS" ] || { echo "--tests and --faults go with --diff-only"; exit 2; }
if [ "$DIFF" = 1 ] && [ -z "$DATA" ]; then  # the differential tests' tars
  DATA=$(grep -v '^#' "$REPO/tests/verify/vm_data.txt" | grep -E '^[A-Za-z0-9_]*diff[A-Za-z0-9_]*\.tar$' | paste -sd,)
elif [ -z "$DATA" ]; then  # the tars from the case registry
  if [ -n "$CASES" ]; then DATA=$(python3 "$REPO/tests/verify/wrf_suite.py" --data-for "$CASES") || exit 2
  else
    DATA=$(grep -v '^#' "$REPO/tests/verify/vm_data.txt" | grep . | paste -sd,)
    if [ "$SLOW" = 1 ]; then
      SC=$(python3 -c 'import sys; sys.path.insert(0, sys.argv[1]); import wrf_suite as w
print(",".join(c[0] for c in w.CASES if c[1] == "cpu_dp" and w.slow(c[3])))' "$REPO/tests/verify")
      for t in $(python3 "$REPO/tests/verify/wrf_suite.py" --data-for "$SC" | tr ',' ' '); do
        case ,$DATA, in *,$t,*) ;; *) DATA=$DATA,$t ;; esac
      done
    fi
  fi
fi
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
if [ "$DIFF" = 1 ]; then
  for f in ${FAULTS:+"$FAULTS"/*.patch}; do cj_send "$f" "in/fault_$(basename "$f")"; done
  cj_launch "$REPO/tools/colab/diff_body.sh" "$DATA" "$TESTS"
  cj_wait "$HOURS" diff.log
  cj_results
  bad=0
  grep . "$OUT/setup.out" 2> /dev/null && bad=1
  cat "$OUT/cpu.txt" 2> /dev/null || true
  for f in diff ctest; do
    [ -f "$OUT/$f.log" ] || { echo "[colab_suite] INCOMPLETE: no $f.log"; bad=1; continue; }
    echo "$f: $(grep -E 'tests passed|No tests were found' "$OUT/$f.log" || echo 'no ctest summary')"
    grep -E '^\s*[0-9]+ - .*\((Failed|Timeout|SEGFAULT|Subprocess|Skipped|Not Run)|^SKIPPED ' "$OUT/$f.log" || true
    grep -qx "exit 0" "$OUT/$f.log" || bad=1
  done
  if [ -n "$FAULTS" ]; then
    cat "$OUT/faults.log" 2> /dev/null || { echo "[colab_suite] no faults.log"; bad=1; }
    [ "$(grep -c ' DETECTED ' "$OUT/faults.log" 2> /dev/null)" = "$(ls "$FAULTS"/*.patch | wc -l)" ] || bad=1
  fi
  exit $bad
fi
THR=${THREADS:-$(( NPROC / JOBS > 0 ? NPROC / JOBS : 1 ))}
echo "[colab_suite] jobs $JOBS, threads $THR"
cj_launch "$REPO/tools/gcp/suite_body.sh" "$DATA" "$JOBS" "$THR" "$CASES" "${SLOW:-0}"
cj_wait "$HOURS" suite.log
cj_results
bash "$REPO/tools/gcp/suite_status.sh" "$OUT" "$CASES"
