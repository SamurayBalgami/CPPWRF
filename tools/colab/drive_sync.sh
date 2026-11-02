#!/usr/bin/env bash
# Copy the suite's data tars from the bucket to Google Drive (folder cppwrf-data of the rclone remote
# cppwrf-drive), on a Colab CPU runtime: the data never pass through the workstation, and tools/colab/colab_suite.sh
# then reads them from Drive without GCS egress fees on every run.
#
#   CPPWRF_GCP_BUCKET=... tools/colab/drive_sync.sh [--par N] [--out DIR] [TAR ...]
#
#   TAR    bucket objects to copy (default: tests/verify/vm_data.txt's tars); a tar whose Drive copy has the
#          bucket object's size is left out, so a new tar of the data list costs only its own transfer
#   --par  tars copied at a time (default 3; each passes through the runtime's disk)
#   --out  log directory (default: ${CPPWRF_SCRATCH:-~/cppwrf_scratch}/colab/drive_sync_<date>)
#
# Once on the workstation (browser login; scope drive.file: rclone sees only the files it created):
#   rclone config create cppwrf-drive drive scope=drive.file
# The runtime gets only that remote's rclone section and a read-only token for the bucket
# (tools/colab/gcs_token.py, renewed every 5 min), both deleted at the end; the session is always stopped.
# At the end each Drive copy's MD5 is compared with the bucket object's. Exit status 0 only when every tar is on
# Drive with the bucket's size and MD5.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
: "${CPPWRF_GCP_BUCKET:?set CPPWRF_GCP_BUCKET}"
BK=$CPPWRF_GCP_BUCKET PAR=3 OUT= TARS=()
while [ $# -gt 0 ]; do
  case $1 in
    --par) PAR=$2; shift 2 ;; --out) OUT=$2; shift 2 ;;
    -*) sed -n '2,18p' "$0"; exit 2 ;; *) TARS+=("$1"); shift ;;
  esac
done
[ ${#TARS[@]} -gt 0 ] || mapfile -t TARS < <(grep -v '^#' "$REPO/tests/verify/vm_data.txt" | grep .)
STAMP=$(date +%Y%m%d_%H%M%S)
OUT=${OUT:-${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/colab/drive_sync_$STAMP}
mkdir -p "$OUT"
rclone listremotes | grep -qx "cppwrf-drive:" || { echo "no rclone remote cppwrf-drive (see the header)"; exit 1; }

# sizes in the bucket and on Drive: only the tars whose Drive copy is missing or has another size, largest first
gcloud storage ls -l "gs://$BK/" | awk '$1 ~ /^[0-9]+$/ {n = $3; sub(".*/", "", n); print n, $1}' > "$OUT/bucket.txt"
# (rclone 1.60 can hang listing a missing folder: create it first)
timeout 120 rclone mkdir cppwrf-drive:cppwrf-data
timeout 300 rclone lsf --format sp --separator ' ' cppwrf-drive:cppwrf-data | awk '{print $2, $1}' > "$OUT/drive.txt"
todo=()
for t in "${TARS[@]}"; do
  b=$(awk -v t="$t" '$1 == t {print $2}' "$OUT/bucket.txt")
  [ -n "$b" ] || { echo "[drive_sync] $t is not in gs://$BK"; exit 1; }
  [ "$(awk -v t="$t" '$1 == t {print $2}' "$OUT/drive.txt")" = "$b" ] || todo+=("$b $t")
done
mapfile -t todo < <(printf '%s\n' "${todo[@]}" | grep . | sort -rn | cut -d' ' -f2)
if [ ${#todo[@]} -gt 0 ]; then
  echo "[drive_sync] ${#todo[@]} tars to copy ($(awk 'NR == FNR {w[$1]; next} $1 in w {s += $2} END {printf "%.1f GB", s / 1e9}' \
    <(printf '%s\n' "${todo[@]}") "$OUT/bucket.txt"))"
  S=cppwrf-sync-$STAMP
  token() { (umask 077; python3 "$REPO/tools/colab/gcs_token.py" "$BK" > "$OUT/gcs.tok") &&
            timeout 300 colab upload -s "$S" "$OUT/gcs.tok" /content/sync/gcs.tok.new 2>&1 | grep -q Uploaded &&
            echo 'import os; os.replace("/content/sync/gcs.tok.new", "/content/sync/gcs.tok"); print("TOKEN")' |
              timeout 120 colab exec -s "$S" --timeout 60 2>/dev/null | grep -q TOKEN; local r=$?; rm -f "$OUT/gcs.tok"; return $r; }
  cleanup() {
    echo 'import subprocess; subprocess.run("rm -f /content/sync/gcs.tok* /content/sync/rclone.conf", shell=True); print("CLEAN")' |
      timeout 120 colab exec -s "$S" --timeout 60 2>/dev/null | grep -q CLEAN || echo "[drive_sync] WARNING: credentials not removed (the session stop deletes the VM)"
    timeout 300 colab stop -s "$S" || echo "[drive_sync] WARNING: could not stop $S, check: colab sessions"
    rm -f "$OUT/gcs.tok" "$OUT/rclone.conf"
  }
  trap cleanup EXIT
  timeout 900 colab new -s "$S"
  echo 'import os; os.makedirs("/content/sync", exist_ok=True); print("OK")' | timeout 120 colab exec -s "$S" --timeout 60
  (umask 077; rclone config show cppwrf-drive > "$OUT/rclone.conf")
  for f in "$OUT/rclone.conf" "$REPO/tools/colab/remote_drive_sync.sh"; do
    timeout 300 colab upload -s "$S" "$f" "/content/sync/$(basename "$f")" 2>&1 | grep -q Uploaded || { echo "[drive_sync] upload of $f failed"; exit 1; }
  done
  rm -f "$OUT/rclone.conf"
  token || { echo "[drive_sync] token upload failed"; exit 1; }
  timeout 600 colab exec -s "$S" --timeout 500 --env SYNC_ARGS="$BK $PAR ${todo[*]}" <<'EOF'
import os, subprocess
r = subprocess.run("chmod 600 /content/sync/rclone.conf /content/sync/gcs.tok; apt-get -qq update && "
                   "DEBIAN_FRONTEND=noninteractive apt-get -qq install -y rclone > /dev/null && rclone version | head -1",
                   shell=True, capture_output=True, text=True)
print(r.stdout, r.stderr[-500:])
subprocess.Popen("nohup bash /content/sync/remote_drive_sync.sh " + os.environ["SYNC_ARGS"] + " > /dev/null 2>&1 &", shell=True)
print("sync started")
EOF
  # a new token every 5 min: the downscoped token expires with the workstation's gcloud token, which gcloud hands
  # out again until it has a few minutes left (5 to 60 min); the runtime retries a 401 for 20 min
  seen=0
  while :; do
    sleep 300
    token || echo "[drive_sync] token renewal failed"
    log=$(echo 'print(open("/content/sync/sync.log").read() if __import__("os").path.exists("/content/sync/sync.log") else "")' |
          timeout 180 colab exec -s "$S" --timeout 60 2>/dev/null | grep -v '^\[colab\]' || true)
    printf '%s\n' "$log" | grep -E '^(TAR|FAILED) ' | tail -n +$((seen + 1)) | while read -r k t sz dl ul; do
      [ "$k" = TAR ] && echo "[drive_sync] $(date +%H:%M) $t $(( sz / 1000000 )) MB: download $dl s ($(( sz / 1000000 / (dl > 0 ? dl : 1) )) MB/s), Drive $ul s ($(( sz / 1000000 / (ul > 0 ? ul : 1) )) MB/s)" ||
        echo "[drive_sync] $(date +%H:%M) FAILED $t"
    done | tee -a "$OUT/sync.txt"
    seen=$(printf '%s\n' "$log" | grep -cE '^(TAR|FAILED) ' || true)
    printf '%s\n' "$log" | grep -qx done && break
  done
fi
# every tar on Drive with the bucket's size and MD5 (a composite object without MD5: size only)
timeout 600 rclone md5sum cppwrf-drive:cppwrf-data > "$OUT/drive_md5.txt"
bad=0
for t in "${TARS[@]}"; do
  want=$(gcloud storage objects describe "gs://$BK/$t" --format='value(md5_hash)' | python3 -c 'import base64, sys; s = sys.stdin.read().strip(); print(base64.b64decode(s).hex() if s else "")')
  got=$(awk -v t="$t" '$2 == t {print $1}' "$OUT/drive_md5.txt")
  sz=$(timeout 120 rclone lsf --format s "cppwrf-drive:cppwrf-data/$t" 2>/dev/null || true)
  if [ "$sz" != "$(awk -v t="$t" '$1 == t {print $2}' "$OUT/bucket.txt")" ] || { [ -n "$want" ] && [ "$want" != "$got" ]; }; then
    echo "[drive_sync] MISMATCH $t (size $sz, md5 $got, bucket md5 $want)"; bad=1
  fi
done
[ $bad = 0 ] && echo "[drive_sync] all ${#TARS[@]} tars on Drive (cppwrf-drive:cppwrf-data) with the bucket's size and MD5"
exit $bad
