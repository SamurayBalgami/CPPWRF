#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The quick start of README.md and docs/USER_GUIDE.md from a clean source tree (git archive of a commit) on ONE
# Google Colab G4 runtime: the CPU FP64 build on its AMD CPU and the CUDA FP32 / FP64 builds on its GPU, with the
# commands of the documentation, then the documented runs.
#
#   CPPWRF_GCP_BUCKET=... tools/colab/colab_quickstart.sh [--wps DIR] [--out DIR] [--rev REV] [--hours H] [--session NAME] [--keep]
#
# CPU side (tools/colab/quickstart_body.sh in the Debian 13 chroot of tools/colab/remote_chroot_job.sh): README's CPU
# build, `ctest`, the translator's tests, the analytic Straka case, WRF's em_quarter_ss with WRF's input_sounding,
# and examples/quickstart/run.sh (GFS -> WPS -> cppwrf_real -> cppwrf -> wrfout, two domains, 3 h) with the CPU build.
# GPU side (tools/colab/quickstart_gpu.sh on the runtime's Ubuntu with CUDA): README's two CUDA builds, then the same
# Straka case and the model step of the example (same wrfinput / wrfbdy / run.ini) in GPU FP32 and FP64. At the end the
# chroot compares the example's wrfout files: GPU FP64 against CPU FP64 at tolerance 0 (tests/verify/exact.py), GPU
# FP32 against CPU FP64 by relative L2.
# Downloaded on the runtime: WRF v4.8.0's release source (its run/ tables and test/em_quarter_ss/input_sounding; WRF is
# not built), WPS's low-resolution geography (UCAR), two GFS files (NOAA open data). --wps: the WPS 4.6.0 directory
# whose geogrid / ungrib / metgrid executables (Debian 13 builds) and tables are uploaded (default
# ~/wrf_real/WPS-4.6.0, ~3 MB). No data tars. --hours: the longest wait (default 2).
# Results in OUT (default ${CPPWRF_SCRATCH:-~/cppwrf_scratch}/colab/quickstart_<date>): quickstart.txt (one line per
# step: PASS / FAIL with times), quickstart.log, the build / run logs, gpu.txt. Exit 1 when a step failed.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
: "${CPPWRF_GCP_BUCKET:?set CPPWRF_GCP_BUCKET (tools/colab/colab_lib.sh sends a read-only token; no tar is read)}"
BK=$CPPWRF_GCP_BUCKET
WPS=$HOME/wrf_real/WPS-4.6.0 OUT= HOURS=2 REV=HEAD S= KEEP=0 NPROC= ANY_CPU=1
while [ $# -gt 0 ]; do
  case $1 in
    --wps) WPS=$2; shift 2 ;; --out) OUT=$2; shift 2 ;; --rev) REV=$2; shift 2 ;; --hours) HOURS=$2; shift 2 ;;
    --session) S=$2; shift 2 ;; --keep) KEEP=1; shift ;;
    *) sed -n '3,21p' "$0"; exit 2 ;;
  esac
done
STAMP=$(date +%Y%m%d_%H%M%S)
OUT=${OUT:-${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/colab/quickstart_$STAMP}
mkdir -p "$OUT"
[ "$REV" != HEAD ] || [ -z "$(git -C "$REPO" status --porcelain --untracked-files=no)" ] || echo "note: uncommitted changes are not sent (HEAD only)"
git -C "$REPO" archive --format=tar.gz -o "$OUT/src.tar.gz" "$REV"
git -C "$REPO" rev-parse --short "$REV" > "$OUT/commit.txt"
T=$(mktemp -d "$OUT/wps.XXXX"); mkdir -p "$T/geogrid" "$T/metgrid" "$T/ungrib/Variable_Tables"
cp -L "$WPS"/geogrid.exe "$WPS"/ungrib.exe "$WPS"/metgrid.exe "$T/"
cp -L "$WPS/geogrid/GEOGRID.TBL" "$T/geogrid/"; cp -L "$WPS/metgrid/METGRID.TBL" "$T/metgrid/"
cp -L "$WPS/ungrib/Variable_Tables/Vtable.GFS" "$T/ungrib/Variable_Tables/"
tar czf "$OUT/wps.tar.gz" -C "$T" .; rm -rf "$T"
echo "[colab_quickstart] commit $(cat "$OUT/commit.txt")"
TAG=colab_quickstart
. "$REPO/tools/colab/colab_lib.sh"
cj_start G4
cj_send "$OUT/wps.tar.gz" in/wps.tar.gz
cj_send "$REPO/tools/colab/quickstart_gpu.sh" qs_gpu.sh
cj_launch "$REPO/tools/colab/quickstart_body.sh" "" "$STAMP"
gpu_up() { echo 'import subprocess; print("GPU-UP" if subprocess.run(["pgrep", "-f", "qs_gpu.sh"]).returncode == 0 else "GPU-DOWN")' |
  rx 60 | grep -q GPU-UP; }
for i in 1 2 3; do  # colab exec can time out although the command ran: start only when it is not running
  gpu_up && break
  echo "import subprocess; subprocess.Popen('nohup bash /content/cpu/qs_gpu.sh $STAMP > /dev/null 2>&1 &', shell=True); print('started')" |
    rx 60 > /dev/null || true
  sleep 10
done
gpu_up || echo "[colab_quickstart] the GPU side did not start"
bad=0
cj_wait "$HOURS" quickstart.log || bad=1
timeout 120 colab usage 2>&1 | grep -v '^\[colab\]' > "$OUT/usage.txt" || true
cj_results || bad=1
rm -f "$OUT/wps.tar.gz"
grep . "$OUT/setup.out" 2> /dev/null && bad=1
cat "$OUT/quickstart.txt" 2> /dev/null || { echo "[colab_quickstart] no quickstart.txt"; bad=1; }
grep -q FAIL "$OUT/quickstart.txt" 2> /dev/null && bad=1
grep -q '^END' "$OUT/quickstart.txt" 2> /dev/null || bad=1
exit $bad
