#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The FP32 accuracy study (docs/FP32_ACCURACY.md) on ONE Google Colab G4 runtime: WRF v4.8.0 real*4 (and real*8 with
# MPI) built and run in the Debian 13 chroot (tools/bench/fp32_body.sh), CPPWRF's CUDA FP32 and FP64 builds on the
# GPU (tools/bench/fp32_gpu.sh), the statistics computed on the runtime (tools/verify/fp32_accuracy.py); only the
# statistics, the report and the figures come back.
#
#   CPPWRF_GCP_BUCKET=... tools/bench/colab_fp32.sh [--gpu G4|L4|A100] [--cases P4,P6,...] [--wps DIR] [--out DIR]
#       [--hours H] [--rev REV] [--session NAME] [--keep] [--no-suite] [--data TAR,...] [--extra]
#
# Cases: tools/bench/fp32_run.py's CASES (default all: P4, P6, P9, N3, V1, T1, SKV, MORR, ACM7 and L for 1 h; P4L
# 72 h and N3L 24 h from the inputs of the suite's long references). The data tars come from Drive (else the bucket)
# as for tools/colab/colab_suite.sh; the long cases' inputs are the first members of the bucket's long72_data.tar
# (run_P4_72h_r8) and long_data.tar (run_N3_24h_r8), read on the runtime with ranged reads (about 230 MB instead of
# 9 GB). --wps: as tools/bench/colab_bench.sh (the large domain L is made on the runtime). --hours: the longest wait
# (default 4); --gpu: the runtime (default G4; another one when no G4 is available). Results in OUT (default
# ${CPPWRF_SCRATCH:-~/cppwrf_scratch}/colab/fp32_<date>): fp1/ (runs.jsonl, stats/*.jsonl, report_rmse.md,
# report_max_abs.md, report_rel_l2.md, svg/, logs/, suite logs), usage.txt. --no-suite: without the suite cases.
# --data: the data tars instead of the default list (FP2's cases: --data fp2_data.tar); the long inputs are fetched
# only for P4L / N3L. --extra: also CPPWRF GPU FP32 without FMA contraction and CPU FP32 (FP2's cases need them).
# Exit 1 when the job failed.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
: "${CPPWRF_GCP_BUCKET:?set CPPWRF_GCP_BUCKET}"
BK=$CPPWRF_GCP_BUCKET
SUITE=1 EXTRA=0 GPU=G4 CASES=all WPS=$HOME/wrf_real/WPS-4.6.0 OUT= HOURS=4 REV=HEAD S= KEEP=0 NPROC= ANY_CPU=1
# P4, P9, V1, T1: suite_data.tar; P6: port_a06 / port_c01; N3: n3_data.tar; SKV: port_a11; MORR: port_c04; ACM7:
# port_c16; the whole 1 h references (the fp64-ref pair and the suite's cuda_dp cases): v05_suite_refs.tar;
# fp1_long.tar: the long cases' inputs (made on the runtime, below)
DATA=suite_data.tar,n3_data.tar,port_a06_data.tar,port_c01_data.tar,port_c16_data.tar,port_c04_data.tar,port_a11_data.tar,v05_suite_refs.tar,v05_files_links.tar,fp1_long.tar
while [ $# -gt 0 ]; do
  case $1 in
    --gpu) GPU=$2; shift 2 ;; --cases) CASES=$2; shift 2 ;; --wps) WPS=$2; shift 2 ;; --out) OUT=$2; shift 2 ;; --hours) HOURS=$2; shift 2 ;;
    --rev) REV=$2; shift 2 ;; --session) S=$2; shift 2 ;; --keep) KEEP=1; shift ;; --no-suite) SUITE=0; shift ;;
    --data) DATA=$2; shift 2 ;; --extra) EXTRA=1; shift ;;
    *) sed -n '3,23p' "$0"; exit 2 ;;
  esac
done
STAMP=$(date +%Y%m%d_%H%M%S)
OUT=${OUT:-${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/colab/fp32_$STAMP}
mkdir -p "$OUT"
[ "$REV" != HEAD ] || [ -z "$(git -C "$REPO" status --porcelain --untracked-files=no)" ] || echo "note: uncommitted changes are not sent (HEAD only)"
git -C "$REPO" archive --format=tar.gz -o "$OUT/src.tar.gz" "$REV"
git -C "$REPO" rev-parse --short "$REV" > "$OUT/commit.txt"
T=$(mktemp -d "$OUT/wps.XXXX"); mkdir -p "$T/geogrid" "$T/metgrid"
cp -L "$WPS"/geogrid.exe "$WPS"/ungrib.exe "$WPS"/metgrid.exe "$T/"
cp -L "$WPS/geogrid/GEOGRID.TBL" "$T/geogrid/"; cp -L "$WPS/metgrid/METGRID.TBL" "$T/metgrid/"
cp -L "$WPS/ungrib/Variable_Tables/Vtable.GFS" "$T/Vtable"; cp "$REPO/tools/bench/namelist.wps" "$T/"
tar czf "$OUT/wps.tar.gz" -C "$T" .; rm -rf "$T"
echo "[colab_fp32] commit $(cat "$OUT/commit.txt"), cases $CASES, data $DATA"
TAG=colab_fp32
. "$REPO/tools/colab/colab_lib.sh"
cj_start "$GPU"
cj_send "$OUT/wps.tar.gz" in/wps.tar.gz
cj_send "$REPO/tools/bench/fp32_gpu.sh" fp32_gpu.sh
# the long cases' inputs: the leading members of two bucket tars (namelist, tables, wrfinput, wrfbdy), up to the
# first wrfout member (offsets from a listing of their headers), joined into one tar that the job unpacks as a cached
# data tar
case ",$CASES," in *,P4L,*|*,N3L,*|,all,) LONG=1 ;; *) LONG=0 ;; esac
[ $LONG = 0 ] || cj_token || { echo "[colab_fp32] token upload failed"; exit 1; }
[ $LONG = 0 ] || rx 1800 <<EOF | grep LONG_OK > /dev/null || { echo "[colab_fp32] the long inputs were not fetched"; exit 1; }
import os, subprocess
T = "/content/tars"; os.makedirs(T, exist_ok=True)
with open(T + "/fp1_long.tar", "wb") as o:
    for obj, end in (("long72_data.tar", 122625536), ("long_data.tar", 106426368)):
        r = subprocess.run(["gcloud", "storage", "--access-token-file=/content/cpu/gcs.tok", "cat", "-r", f"0-{end - 1}",
                            "gs://$BK/" + obj], stdout=o)
        assert r.returncode == 0 and o.tell() % 512 == 0, obj
    o.write(b"\0" * 10240)
n = subprocess.run(["tar", "tf", T + "/fp1_long.tar"], capture_output=True, text=True).stdout.split()
assert any(x.endswith("run_P4_72h_r8/wrfbdy_d01") for x in n) and any(x.endswith("run_N3_24h_r8/wrfinput_d02") for x in n), n[-5:]
open(T + "/fp1_long.tar.ok", "w").close()
print("LONG_OK", os.path.getsize(T + "/fp1_long.tar"), len(n))
EOF
cj_launch "$REPO/tools/bench/fp32_body.sh" "$DATA" "$CASES" "$STAMP"
echo "import subprocess; subprocess.Popen('FP32_EXTRA_BUILDS=$EXTRA nohup bash /content/cpu/fp32_gpu.sh $STAMP $SUITE > /dev/null 2>&1 &', shell=True); print('started')" |
  rx 60 | grep started > /dev/null || echo "[colab_fp32] the GPU side did not start"
bad=0
cj_wait "$HOURS" ../fp1/out/fp1.log || bad=1
timeout 120 colab usage 2>&1 | grep -v '^\[colab\]' > "$OUT/usage.txt" || true
cj_results || bad=1
if timeout 1800 colab download -s "$S" /content/cpu/fetch/fp1.tar.gz "$OUT/fp1.tar.gz" > /dev/null 2>&1; then
  tar xzf "$OUT/fp1.tar.gz" -C "$OUT" && rm -f "$OUT/fp1.tar.gz"
else echo "[colab_fp32] fp1.tar.gz: download failed"; bad=1; fi
grep . "$OUT/setup.out" 2> /dev/null && bad=1
cat "$OUT/fp1/fp1.log" 2> /dev/null | tail -60 || true
exit $bad
