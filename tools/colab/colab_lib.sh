# The workstation side of a job in the Debian 13 chroot of a Colab runtime with an Intel CPU
# (tools/colab/remote_chroot_job.sh), sourced by tools/colab/colab_suite.sh and colab_wrfref.sh. They set REPO, BK
# (the bucket), OUT (local results; OUT/src.tar.gz is the source to send), TAG (log prefix), S (session; empty: a new
# one) and KEEP (1: leave the session running), then call:
#   cj_start GPU          new session (or S), stop it on exit; fails unless the CPU is Intel; sets NPROC
#   cj_send FILE NAME     local FILE -> /content/cpu/NAME (in/NAME: an input of the job, copied to /opt/job)
#   cj_launch BODY DATA ARG...   the job: `bash BODY /opt/job ARG...` in the chroot after DATA's tars are unpacked
#   cj_wait HOURS LOG     poll every 5 min (LOG: the job's log under /opt/job shown with the progress), a new bucket
#                         token each time while the data are fetched (it expires with the workstation's gcloud token,
#                         5 to 60 min; the runtime retries a 401 for 20 min); fails when the job died or took longer
#   cj_results            OUT gets the job's logs (results.tar.gz), prints versions.txt and data.txt
# The runtime gets only the rclone section of the Drive remote cppwrf-drive and read-only bucket tokens
# (tools/colab/gcs_token.py); the token is deleted when the data are in place, both at the end and on exit.

rx() { timeout "$(( $1 + 120 ))" colab exec -s "$S" --timeout "$1" 2>/dev/null | grep -v '^\[colab\]'; }  # python on stdin

cj_cleanup() {
  rm -f "$OUT/gcs.tok" "$OUT/rclone.conf"
  echo 'import subprocess; subprocess.run("rm -f /content/cpu/gcs.tok* /content/cpu/rclone.conf", shell=True)' | rx 60 > /dev/null || true
  if [ "$KEEP" = 1 ]; then echo "[$TAG] --keep: session $S left running (colab stop -s $S)"
  else timeout 300 colab stop -s "$S" || echo "[$TAG] WARNING: could not stop $S, check: colab sessions"; fi
}

cj_start() {
  local gpu=$1 cpu
  rclone listremotes | grep -qx "cppwrf-drive:" || { echo "[$TAG] no rclone remote cppwrf-drive (tools/colab/drive_sync.sh)"; return 1; }
  trap cj_cleanup EXIT
  if [ -z "$S" ]; then
    S=cppwrf-cpu-$(echo "$gpu" | tr A-Z a-z)-$(date +%Y%m%d_%H%M%S)
    echo "[$TAG] new session $S ($gpu)"; timeout 900 colab new -s "$S" --gpu "$gpu"
  fi
  cpu=$(echo 'print(open("/proc/cpuinfo").read().split("model name")[1].split("\n")[0].strip(" \t:"), "|", __import__("os").cpu_count())' | rx 60 || true)
  echo "[$TAG] CPU: $cpu"
  case $cpu in *Intel*) ;; *) echo "[$TAG] not an Intel runtime: the references cannot be reproduced here"; return 1 ;; esac
  NPROC=${cpu##*| }
  echo 'import os; os.makedirs("/content/cpu/in", exist_ok=True); print("OK")' | rx 60 | grep -q OK || { echo "[$TAG] session not answering"; return 1; }
}

cj_send() {  # the contents API takes files up to ~100 MB (the source and the inputs are a few MB)
  timeout 900 colab upload -s "$S" "$1" "/content/cpu/$2" 2>&1 | grep -q Uploaded || { echo "[$TAG] upload of $1 failed"; return 1; }
}

cj_token() {  # a new read-only bucket token (at most 1 h)
  (umask 077; python3 "$REPO/tools/colab/gcs_token.py" "$BK" > "$OUT/gcs.tok") && cj_send "$OUT/gcs.tok" gcs.tok.new &&
    echo 'import os; os.replace("/content/cpu/gcs.tok.new", "/content/cpu/gcs.tok"); print("TOKEN")' | rx 60 | grep -q TOKEN
  local r=$?; rm -f "$OUT/gcs.tok"; return $r
}

cj_launch() {
  local body=$1 data=$2 t0=$SECONDS a args=; shift 2
  (umask 077; rclone config show cppwrf-drive > "$OUT/rclone.conf")
  cj_send "$OUT/src.tar.gz" src.tar.gz && cj_send "$body" in/body.sh &&
    cj_send "$REPO/tools/colab/remote_chroot_job.sh" remote_chroot_job.sh && cj_send "$OUT/rclone.conf" rclone.conf || return 1
  rm -f "$OUT/rclone.conf"
  cj_token || { echo "[$TAG] token upload failed"; return 1; }
  echo "[$TAG] uploaded $(( $(stat -c %s "$OUT/src.tar.gz") / 1024 )) KiB of source and tools in $((SECONDS - t0)) s"
  for a in "$BK" "$data" "$@"; do args="$args '${a//\'/}'"; done
  echo "import subprocess; subprocess.Popen(\"chmod 600 /content/cpu/rclone.conf /content/cpu/gcs.tok; nohup bash /content/cpu/remote_chroot_job.sh$args > /dev/null 2>&1 &\", shell=True); print('started')" |
    rx 60 | grep -q started || { echo "[$TAG] the job did not start"; return 1; }
}

cj_wait() {
  local hours=$1 log=$2 miss=0 end=$((SECONDS + $1 * 3600)) st=
  while :; do
    sleep 300
    st=$(rx 60 <<EOF || true
import os, subprocess
C = "/content/cpu"
alive = subprocess.run(["pgrep", "-f", "remote_chroot_job.sh"], capture_output=True).returncode == 0
print("JOB-DONE" if os.path.exists(C + "/done") else "JOB-RUNNING" if alive else "JOB-DIED")
print("DATA-DONE" if os.path.exists(C + "/data_done") else "DATA-RUNNING")
for p in (C + "/data.txt", "/content/deb13/opt/job/$log"):
    if os.path.exists(p):
        ls = open(p).read().splitlines()
        print(os.path.basename(p), len(ls), "lines;", ls[-1][:150] if ls else "")
EOF
)
    echo "[$TAG] $(date +%H:%M) $(echo "$st" | tr '\n' ' ')"
    case $st in *JOB-DONE*) return 0 ;; *JOB-DIED*) echo "[$TAG] the job ended without its done file (killed?)"; return 1 ;; esac
    case $st in *JOB-RUNNING*) miss=0 ;; *) miss=$((miss + 1)); [ $miss -lt 5 ] || { echo "[$TAG] session not answering, giving up"; return 1; } ;; esac
    case $st in *DATA-RUNNING*) cj_token || echo "[$TAG] token renewal failed" ;; esac
    [ $SECONDS -lt $end ] || { echo "[$TAG] the job did not end within $hours h"; return 1; }
  done
}

cj_results() {
  timeout 1800 colab download -s "$S" /content/cpu/results.tar.gz "$OUT/results.tar.gz" > /dev/null 2>&1 || { echo "[$TAG] download failed"; return 1; }
  tar xzf "$OUT/results.tar.gz" -C "$OUT" && rm -f "$OUT/results.tar.gz" "$OUT/src.tar.gz"
  echo "[$TAG] results in $OUT"; cat "$OUT/versions.txt" 2>/dev/null || true
  awk '$1 == "TAR" {printf "[%s] data %s: %.0f MB in %d s from %s\n", tag, $2, $3 / 1e6, $4, $5}
       $1 == "DRIVE" {printf "[%s] to Drive: %s (%.0f MB)\n", tag, $2, $3 / 1e6}' tag="$TAG" "$OUT/data.txt" 2>/dev/null || true
}
