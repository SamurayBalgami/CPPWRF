#!/usr/bin/env bash
# Run CPPWRF's GPU suite cases (and timings) on a Google Colab GPU VM with the `colab` CLI
# (google-colab-cli, authorised once with `colab usage`). The VM is always released at the end.
#
#   tools/colab/colab_run.sh [--gpu T4|L4|G4|A100|H100] (default G4: RTX PRO 6000, the fastest per credit)
#                            [--cases NAME,...] [--files NAME,...] [--bench NAME,...] [--keep] [--out DIR]
#
#   --cases  wrf_suite.py cases to run on the VM against the WRF references (default: the GPU core set below)
#   --files  run_wrf_files_test.sh cases to run with the CUDA FP64 build (files' structure and times; default none)
#   --bench  suite cases to run again without history output for the stepping time (default: P1 and D0)
#   --keep   leave the VM running (stop it later with `colab stop -s <session>`)
#   --out    local results directory (default: ${CPPWRF_SCRATCH:-~/cppwrf_scratch}/colab/<date>_<gpu>)
#   MAX_UPLOAD_MIN (environment, default 20): stop when the upload would take longer (slow regions)
#
# The source is the committed HEAD (git archive); the data come from ${CPPWRF_WRF_RUNS} through
# tools/colab/stage_data.py (cached per case list). Results: summary.txt (GPU, build times, PASS/FAIL lines
# with the relative L2 differences, ns/cell/step) and the logs.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
: "${CPPWRF_WRF_RUNS:?set CPPWRF_WRF_RUNS}"
SCRATCH=${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}
GPU=G4 KEEP=0 OUT= FILES= MAX_UPLOAD_MIN=${MAX_UPLOAD_MIN:-20}
CASES=real_r1_fp64_cuda,real_p1_r8_fp64_cuda_1h,real_p1_fp32_1h,real_p4_r8_fp64_cuda_1h,real_p8_r8_fp64_cuda_1h,real_t6_r8_fp64_cuda_1h,real_d0_r8_fp64_cuda_1h,real_nmp1_r8_fp64_cuda_1h
BENCH=real_p1_r8_fp64_cuda_1h,real_p1_fp32_1h,real_d0_r8_fp64_cuda_1h
while [ $# -gt 0 ]; do
  case $1 in
    --gpu) GPU=$2; shift 2 ;;
    --cases) CASES=$2; shift 2 ;;
    --files) FILES=$2; shift 2 ;;
    --bench) BENCH=$2; shift 2 ;;
    --keep) KEEP=1; shift ;;
    --out) OUT=$2; shift 2 ;;
    *) sed -n '2,17p' "$0"; exit 2 ;;
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
key=$( (echo "$CASES,$BENCH" | tr ',' '\n' | sort -u; echo "files:$FILES") | tr '\n' ',' | md5sum | cut -c1-12)
DATA=$SCRATCH/colab/cache/data_$key.tar
# GPU runs only: CPPWRF computes the Thompson tables on the VM instead of uploading WRF's 370 MB
[ -f "$DATA" ] || python3 "$REPO/tools/colab/stage_data.py" --only "$CASES,$BENCH" --files "$FILES" --no-thompson-tables "$DATA"

# ---- remote
timeout 120 colab usage | tail -3 || true
echo "[colab_run] new session $S ($GPU)"
cleanup() {
  if [ "$KEEP" = 1 ]; then echo "[colab_run] --keep: session $S left running (colab stop -s $S)"
  else timeout 300 colab stop -s "$S" || echo "[colab_run] WARNING: could not stop $S, check: colab sessions"; fi
}
trap cleanup EXIT
timeout 900 colab new -s "$S" --gpu "$GPU"
# uploads go through Jupyter's contents API, which refuses large files (and `colab upload` exits 0 anyway):
# send 32 MiB parts, join them on the VM and check the md5 (each part is cut just before it is sent: the disk
# never holds a second copy of the data)
upload() {
  local src=$1 dst=$2 parts p n=0 sum got
  parts=$(mktemp -d "$SCRATCH/colab/parts.XXXXXX")
  local total t0=$SECONDS
  total=$(( ($(stat -c %s "$src") + 33554431) / 33554432 ))
  [ "$total" -gt 0 ] || total=1
  while [ "$n" -lt "$total" ]; do
    p=$parts/p$(printf %04d "$n")
    dd if="$src" of="$p" bs=32M skip="$n" count=1 iflag=fullblock status=none
    if ! timeout 900 colab upload -s "$S" "$p" "$dst.part.$(basename "$p")" 2>&1 | grep -q "Uploaded"; then
      echo "[colab_run] upload of part $(basename "$p") of $src failed"; rm -rf "$parts"; return 1
    fi
    rm -f "$p"
    n=$((n + 1))
    # the VM bills while it waits: give up when this region's upload would take too long
    local eta=$(( (SECONDS - t0) * (total - n) / n / 60 ))
    [ "$total" -gt 1 ] && echo "[colab_run] $(basename "$src"): $n/$total parts, ~${eta} min left"
    if [ "$n" -ge 2 ] && [ "$eta" -gt "$MAX_UPLOAD_MIN" ]; then
      echo "[colab_run] upload too slow (~${eta} min left > MAX_UPLOAD_MIN=$MAX_UPLOAD_MIN), giving up"
      rm -rf "$parts"; return 1
    fi
  done
  rm -rf "$parts"
  sum=$(md5sum "$src" | cut -c1-32)
  got=$(printf '%s\n' "import glob, hashlib, os" \
    "parts = sorted(glob.glob('$dst.part.p*'))" \
    "o = open('$dst', 'wb')" \
    "for p in parts: o.write(open(p, 'rb').read()); os.remove(p)" \
    "o.close()" \
    "print('MD5', hashlib.md5(open('$dst', 'rb').read()).hexdigest())" | timeout 900 colab exec -s "$S" --timeout 600 || true)
  case $got in
    *"MD5 $sum"*) echo "[colab_run] uploaded $src ($n parts)" ;;
    *) echo "[colab_run] $dst: md5 mismatch"; return 1 ;;
  esac
}
upload "$OUT/src.tar.gz" /content/src.tar.gz || exit 1
upload "$DATA" /content/data.tar || exit 1
upload "$OUT/commit.txt" /content/commit.txt || exit 1
upload "$REPO/tools/colab/remote_job.py" /content/remote_job.py || exit 1
fenv=(); [ -z "$FILES" ] || fenv=(--env CPPWRF_JOB_FILES="$FILES")
timeout 300 colab exec -s "$S" --timeout 120 --env CPPWRF_JOB_CASES="$CASES" --env CPPWRF_JOB_BENCH="$BENCH" "${fenv[@]}" <<'EOF'
import os, subprocess
env = dict(os.environ)
subprocess.Popen("nohup python3 /content/remote_job.py > /content/remote_job.log 2>&1 &", shell=True, env=env)
print("job started")
EOF
# poll: the job writes /content/job/done (also when it fails); JOB-DIED: no done file and no job process
miss=0
while :; do
  sleep 60
  st=$(timeout 180 colab exec -s "$S" --timeout 60 <<'EOF' 2>/dev/null || true
import os, subprocess
alive = subprocess.run(["pgrep", "-f", "remote_job.py"], capture_output=True).returncode == 0
print("JOB-DONE" if os.path.exists("/content/job/done") else "JOB-RUNNING" if alive else "JOB-DIED")
s = "/content/job/summary.txt"
print(open(s).read().splitlines()[-1] if os.path.exists(s) and open(s).read() else "(starting)")
EOF
)
  echo "[colab_run] $(date +%H:%M) $(echo "$st" | grep -v '^\[colab\]' | tr '\n' ' ')"
  case $st in *JOB-DONE*|*JOB-DIED*) break ;; esac
  case $st in
    *JOB-RUNNING*) miss=0 ;;
    *) miss=$((miss + 1)); echo "[colab_run] session not answering ($miss)"; timeout 120 colab sessions || true
       [ $miss -lt 5 ] || { echo "[colab_run] giving up"; exit 1; } ;;
  esac
done
case $st in *JOB-DIED*) echo "[colab_run] the job ended without its done file (killed?)"; exit 1 ;; esac
timeout 1800 colab download -s "$S" /content/results.tar.gz "$OUT/results.tar.gz" || { echo "[colab_run] download failed"; exit 1; }
tar xzf "$OUT/results.tar.gz" -C "$OUT"
rm -f "$OUT/src.tar.gz"
echo "[colab_run] results in $OUT/job"
cat "$OUT/job/summary.txt"
if grep -qE "^suite (FAIL|SKIP|ERROR)|^files (FAIL|SKIP)|^job error" "$OUT/job/summary.txt" ||
   { [ -n "$FILES" ] && ! grep -q "^files compared " "$OUT/job/summary.txt"; } ||
   { [ -n "$CASES" ] && ! grep -q "^suite SUITE " "$OUT/job/summary.txt"; }; then
  echo "[colab_run] the suite did not pass"; exit 1
fi
