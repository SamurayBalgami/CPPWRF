#!/usr/bin/env bash
# WRF v4.8.0 against CPPWRF on ONE Google Colab G4 runtime (AMD EPYC 9B45, 48 vCPUs, 176 GB + NVIDIA RTX PRO 6000
# Blackwell, 96 GB): WRF's MPI builds (real*4, real*8) on the CPU, CPPWRF's CUDA FP32 / FP64 builds on the GPU and its
# OpenMP FP64 build on the CPU, the same namelists and inputs, 1 simulated hour each.
#
#   CPPWRF_GCP_BUCKET=... tools/bench/colab_bench.sh [--cases R1,P4,P4R15,N3,L] [--layouts dm12,dm24,dm48,sm24x2,sm12x4]
#       [--wps DIR] [--out DIR] [--hours H] [--rev REV] [--session NAME] [--keep]
#
# Cases (WRF real*8 run directories of the suite's data tars, run for 1 h): R1 (dynamics only), P4 (full physics,
# radt 1 min), P4R15 (P4 with radt 15 min), N3 (two-way nest with P4 physics), L (a 9.1 times larger domain made on the
# runtime from GFS with WPS and real.exe, P4 physics with radt 15 min). --layouts: WRF's process layouts (dmN: N MPI
# ranks; smNxT: N ranks x T OpenMP threads of the dm+sm build), each tried with history off; the fastest runs again
# with history on; R1, P4 and L also run WRF real*8 in one process (the exact comparisons: WRF's result depends on its
# MPI layout, the suite's references are one-process runs). --wps: the WPS 4.6.0 directory whose geogrid / ungrib / metgrid executables (Debian 13 builds) and
# tables made the suite's met_em files (default ~/wrf_real/WPS-4.6.0; ~3 MB are uploaded). Everything else is
# downloaded on the runtime: WRF's release source (github.com/wrf-model/WRF), WPS's low-resolution geography (UCAR),
# the GFS files (NOAA open data); the data tars come from Drive / the bucket as for tools/colab/colab_suite.sh.
# --hours: the longest wait (default 5); --rev, --session, --keep as colab_suite.sh (a --keep session keeps the WRF
# builds, the met_em files and the CUDA builds for the next run with --session).
#
# The job: tools/bench/bench_body.sh in the Debian 13 chroot (tools/colab/remote_chroot_job.sh: gfortran 14, MPICH)
# and tools/bench/bench_gpu.sh on the runtime's own Ubuntu (CUDA); the runs are tools/bench/bench_run.py's.
# One session takes about 75 minutes of the G4 (12 compute units in November 2026); docs/PERF.md has its results.
# Results in OUT (default ${CPPWRF_SCRATCH:-~/cppwrf_scratch}/colab/bench_<date>): report.md (tools/bench/bench_report.py
# with the runtime's compute-unit rate from `colab usage`), results.txt (one JSON line per run and comparison),
# bench.log, logs/ (rsl.out.0000 of each WRF run, CPPWRF's logs, namelists, configurations, WRF's configure.wrf
# flags), usage.txt. Exit 1 when the job failed or a build / case did not run.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
: "${CPPWRF_GCP_BUCKET:?set CPPWRF_GCP_BUCKET}"
BK=$CPPWRF_GCP_BUCKET
CASES=P4,L,P4R15,N3,R1 LAYOUTS=dm12,dm24,dm48,sm24x2,sm12x4 WPS=$HOME/wrf_real/WPS-4.6.0 OUT= HOURS=5 REV=HEAD S= KEEP=0
NPROC= ANY_CPU=1
# R1: dycore_data.tar; P4: suite_data.tar; N3: n3_data.tar; the whole 1 h references of P4 and N3: v05_suite_refs.tar
DATA=dycore_data.tar,suite_data.tar,n3_data.tar,v05_suite_refs.tar
while [ $# -gt 0 ]; do
  case $1 in
    --cases) CASES=$2; shift 2 ;; --layouts) LAYOUTS=$2; shift 2 ;; --wps) WPS=$2; shift 2 ;; --out) OUT=$2; shift 2 ;;
    --hours) HOURS=$2; shift 2 ;; --rev) REV=$2; shift 2 ;; --session) S=$2; shift 2 ;; --keep) KEEP=1; shift ;;
    *) sed -n '2,27p' "$0"; exit 2 ;;
  esac
done
STAMP=$(date +%Y%m%d_%H%M%S)
OUT=${OUT:-${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/colab/bench_$STAMP}
mkdir -p "$OUT"
[ "$REV" != HEAD ] || [ -z "$(git -C "$REPO" status --porcelain --untracked-files=no)" ] || echo "note: uncommitted changes are not sent (HEAD only)"
git -C "$REPO" archive --format=tar.gz -o "$OUT/src.tar.gz" "$REV"
git -C "$REPO" rev-parse --short "$REV" > "$OUT/commit.txt"
T=$(mktemp -d "$OUT/wps.XXXX"); mkdir -p "$T/geogrid" "$T/metgrid"
cp -L "$WPS"/geogrid.exe "$WPS"/ungrib.exe "$WPS"/metgrid.exe "$T/"
cp -L "$WPS/geogrid/GEOGRID.TBL" "$T/geogrid/"; cp -L "$WPS/metgrid/METGRID.TBL" "$T/metgrid/"
cp -L "$WPS/ungrib/Variable_Tables/Vtable.GFS" "$T/Vtable"; cp "$REPO/tools/bench/namelist.wps" "$T/"
tar czf "$OUT/wps.tar.gz" -C "$T" .; rm -rf "$T"
echo "[colab_bench] commit $(cat "$OUT/commit.txt"), cases $CASES, WRF layouts $LAYOUTS, data $DATA"
TAG=colab_bench
. "$REPO/tools/colab/colab_lib.sh"
cj_start G4
cj_send "$OUT/wps.tar.gz" in/wps.tar.gz
cj_send "$REPO/tools/bench/bench_gpu.sh" gpu.sh
cj_launch "$REPO/tools/bench/bench_body.sh" "$DATA" "$CASES" "$LAYOUTS" "$STAMP"
echo "import subprocess; subprocess.Popen('nohup bash /content/cpu/gpu.sh $STAMP > /dev/null 2>&1 &', shell=True); print('started')" |
  rx 60 | grep started > /dev/null || echo "[colab_bench] the GPU side did not start"
bad=0
cj_wait "$HOURS" bench.log || bad=1
timeout 120 colab usage 2>&1 | grep -v '^\[colab\]' > "$OUT/usage.txt" || true
cj_results || bad=1
if timeout 1800 colab download -s "$S" /content/cpu/fetch/bench.tar.gz "$OUT/bench.tar.gz" > /dev/null 2>&1; then
  tar xzf "$OUT/bench.tar.gz" -C "$OUT" && rm -f "$OUT/bench.tar.gz"
else echo "[colab_bench] bench.tar.gz: download failed"; bad=1; fi
grep . "$OUT/setup.out" 2> /dev/null && bad=1
rate=$(sed -nE 's/.*Usage rate: *([0-9.]+).*/\1/p' "$OUT/usage.txt" | head -1)
if [ -f "$OUT/results.txt" ]; then
  python3 "$REPO/tools/bench/bench_report.py" "$OUT/results.txt" ${rate:+--rate "$rate"} > "$OUT/report.md"
  echo "[colab_bench] G4 rate ${rate:-?} compute units / hour; report: $OUT/report.md"
else echo "[colab_bench] no results.txt"; bad=1; fi
cat "$OUT/bench.log" 2> /dev/null | tail -40 || true
exit $bad
