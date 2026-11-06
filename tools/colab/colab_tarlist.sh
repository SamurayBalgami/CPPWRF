#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The file lists of data tars on Google Drive (cppwrf-drive:cppwrf-data, tools/colab/drive_sync.sh), read on a Colab
# runtime (the tars stream from Drive there; nothing is downloaded to this machine but the lists):
#
#   tools/colab/colab_tarlist.sh [--out DIR] [--session NAME] [--keep] TAR[,TAR...]   (default: tests/verify/vm_data.txt)
#
# DIR (default ${CPPWRF_SCRATCH:-~/cppwrf_scratch}/colab/tarlist_<date>) gets TAR.txt per tar (`tar tvf`, ending with
# "exit N") and job.log. tests/verify/case_data.py DIR checks the suite cases' data tars (DATA) against them.
# A CPU-only runtime (the CPU does not matter here); it is stopped at the end unless --keep.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT= S= KEEP=0 TARS=
while [ $# -gt 0 ]; do
  case $1 in
    --out) OUT=$2; shift 2 ;; --session) S=$2; shift 2 ;; --keep) KEEP=1; shift ;;
    -*) sed -n '3,11p' "$0"; exit 2 ;; *) TARS=$1; shift ;;
  esac
done
TARS=${TARS:-$(grep -v '^#' "$REPO/tests/verify/vm_data.txt" | grep . | paste -sd,)}
OUT=${OUT:-${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/colab/tarlist_$(date +%Y%m%d_%H%M%S)}
mkdir -p "$OUT"
TAG=colab_tarlist ANY_CPU=1 BK=
. "$REPO/tools/colab/colab_lib.sh"
cj_start CPU
(umask 077; rclone config show cppwrf-drive > "$OUT/rclone.conf")
cj_send "$OUT/rclone.conf" rclone.conf || exit 1
rm -f "$OUT/rclone.conf"
cat > "$OUT/tarlist_job.sh" <<EOF
C=/content/cpu
exec > \$C/job.log 2>&1
rm -rf \$C/ls \$C/ls.tar.gz \$C/done; mkdir -p \$C/ls
chmod 600 \$C/rclone.conf
command -v rclone > /dev/null || { apt-get -qq update && apt-get -qq install -y rclone > /dev/null; }
for t in ${TARS//,/ }; do
  while [ "\$(jobs -rp | wc -l)" -ge 3 ]; do sleep 5; done
  ( t0=\$SECONDS; rclone --config \$C/rclone.conf cat "cppwrf-drive:cppwrf-data/\$t" | tar tvf - > \$C/ls/\$t.txt 2>&1
    echo "exit \$? \$((SECONDS - t0)) s" >> \$C/ls/\$t.txt ) &
done
wait
rm -f \$C/rclone.conf
tar czf \$C/ls.tar.gz -C \$C/ls .
date > \$C/done
EOF
cj_send "$OUT/tarlist_job.sh" tarlist_job.sh || exit 1
echo 'import subprocess; subprocess.Popen("nohup bash /content/cpu/tarlist_job.sh > /dev/null 2>&1 &", shell=True); print("started")' |
  rx 60 | grep started > /dev/null || { echo "[$TAG] the job did not start"; exit 1; }
echo "[$TAG] listing ${TARS//,/ }"
for i in $(seq 72); do  # at most 6 h
  sleep 300
  st=$(echo 'import os, glob; print("DONE" if os.path.exists("/content/cpu/done") else "RUNNING", len(glob.glob("/content/cpu/ls/*.txt")))' | rx 60 | grep . || true)
  echo "[$TAG] $(date +%H:%M) $st"
  case $st in DONE*) break ;; esac
done
timeout 900 colab download -s "$S" /content/cpu/ls.tar.gz "$OUT/ls.tar.gz" > /dev/null 2>&1 || { echo "[$TAG] download failed"; exit 1; }
timeout 300 colab download -s "$S" /content/cpu/job.log "$OUT/job.log" > /dev/null 2>&1 || true
tar xzf "$OUT/ls.tar.gz" -C "$OUT" && rm -f "$OUT/ls.tar.gz"
grep -H "^exit" "$OUT"/*.txt
! grep -h "^exit" "$OUT"/*.txt | grep -qv "^exit 0 "
