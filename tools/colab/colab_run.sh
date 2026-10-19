#!/usr/bin/env bash
# Run CPPWRF's GPU suite cases (and timings) on a Google Colab GPU VM with the `colab` CLI
# (google-colab-cli, authorised once with `colab usage`). The VM is always released at the end.
#
#   tools/colab/colab_run.sh [--gpu T4|L4|G4|A100|H100] [--cases NAME,...] [--bench NAME,...] [--keep] [--out DIR]
#
#   --cases  wrf_suite.py cases to run on the VM against the WRF references (default: the GPU core set below)
#   --bench  suite cases to run again without history output for the stepping time (default: P1 and D0)
#   --keep   leave the VM running (stop it later with `colab stop -s <session>`)
#   --out    local results directory (default: ${CPPWRF_SCRATCH:-~/cppwrf_scratch}/colab/<date>_<gpu>)
#
# The source is the committed HEAD (git archive); the data come from ${CPPWRF_WRF_RUNS} through
# tools/colab/stage_data.py (cached per case list). Results: summary.txt (GPU, build times, PASS/FAIL lines
# with the relative L2 differences, ns/cell/step) and the logs.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
: "${CPPWRF_WRF_RUNS:?set CPPWRF_WRF_RUNS}"
SCRATCH=${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}
GPU=L4 KEEP=0 OUT=
CASES=real_r1_fp64_cuda,real_p1_r8_fp64_cuda_1h,real_p1_fp32_1h,real_p4_r8_fp64_cuda_1h,real_p8_r8_fp64_cuda_1h,real_t6_r8_fp64_cuda_1h,real_d0_r8_fp64_cuda_1h,real_nmp1_r8_fp64_cuda_1h
BENCH=real_p1_r8_fp64_cuda_1h,real_p1_fp32_1h,real_d0_r8_fp64_cuda_1h
while [ $# -gt 0 ]; do
  case $1 in
    --gpu) GPU=$2; shift 2 ;;
    --cases) CASES=$2; shift 2 ;;
    --bench) BENCH=$2; shift 2 ;;
    --keep) KEEP=1; shift ;;
    --out) OUT=$2; shift 2 ;;
    *) sed -n '2,15p' "$0"; exit 2 ;;
  esac
done
case $GPU in T4|L4|G4|A100|H100) ;; *) echo "unknown --gpu $GPU (the CLI would silently take an A100)"; exit 2 ;; esac
STAMP=$(date +%Y%m%d_%H%M%S)
OUT=${OUT:-$SCRATCH/colab/${STAMP}_$GPU}
S=cppwrf-$(echo "$GPU" | tr A-Z a-z)-$STAMP
mkdir -p "$OUT" "$SCRATCH/colab/cache"

# ---- local: source and data
[ -z "$(git -C "$REPO" status --porcelain --untracked-files=no)" ] || echo "note: uncommitted changes are not sent (HEAD only)"
git -C "$REPO" archive --format=tar.gz -o "$OUT/src.tar.gz" HEAD
git -C "$REPO" rev-parse --short HEAD > "$OUT/commit.txt"
key=$(echo "$CASES,$BENCH" | tr ',' '\n' | sort -u | tr '\n' ',' | md5sum | cut -c1-12)
DATA=$SCRATCH/colab/cache/data_$key.tar
[ -f "$DATA" ] || python3 "$REPO/tools/colab/stage_data.py" --only "$CASES,$BENCH" "$DATA"

# ---- remote
colab usage | tail -3
echo "[colab_run] new session $S ($GPU)"
colab new -s "$S" --gpu "$GPU"
cleanup() {
  if [ "$KEEP" = 1 ]; then echo "[colab_run] --keep: session $S left running (colab stop -s $S)"
  else colab stop -s "$S" || echo "[colab_run] WARNING: could not stop $S, check: colab sessions"; fi
}
trap cleanup EXIT
for f in "$OUT/src.tar.gz:/content/src.tar.gz" "$DATA:/content/data.tar" "$OUT/commit.txt:/content/commit.txt" \
         "$REPO/tools/colab/remote_job.py:/content/remote_job.py"; do
  echo "[colab_run] upload ${f%%:*}"
  colab upload -s "$S" "${f%%:*}" "${f#*:}"
done
colab exec -s "$S" --timeout 120 --env CPPWRF_JOB_CASES="$CASES" --env CPPWRF_JOB_BENCH="$BENCH" <<'EOF'
import os, subprocess
env = dict(os.environ)
subprocess.Popen("nohup python3 /content/remote_job.py > /content/remote_job.log 2>&1 &", shell=True, env=env)
print("job started")
EOF
# poll: the job writes /content/job/done; show the summary as it grows
miss=0
while :; do
  sleep 60
  st=$(colab exec -s "$S" --timeout 60 <<'EOF' 2>/dev/null || true
import os
print("DONE" if os.path.exists("/content/job/done") else "RUNNING")
s = "/content/job/summary.txt"
print(open(s).read().splitlines()[-1] if os.path.exists(s) and open(s).read() else "(starting)")
EOF
)
  echo "[colab_run] $(date +%H:%M) $(echo "$st" | tr '\n' ' ')"
  case $st in *DONE*) break ;; esac
  case $st in
    *RUNNING*) miss=0 ;;
    *) miss=$((miss + 1)); echo "[colab_run] session not answering ($miss)"; colab sessions || true
       [ $miss -lt 5 ] || { echo "[colab_run] giving up"; exit 1; } ;;
  esac
done
colab download -s "$S" /content/results.tar.gz "$OUT/results.tar.gz"
tar xzf "$OUT/results.tar.gz" -C "$OUT"
rm -f "$OUT/src.tar.gz"
echo "[colab_run] results in $OUT/job"
cat "$OUT/job/summary.txt"
