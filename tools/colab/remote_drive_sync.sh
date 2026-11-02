#!/usr/bin/env bash
# The Colab side of tools/colab/drive_sync.sh: copies data tars from the bucket to the Google Drive folder
# cppwrf-drive:cppwrf-data, PAR tars at a time, each through the runtime's disk (gcloud storage cp with the
# downscoped token /content/sync/gcs.tok, which drive_sync.sh renews; rclone copyto; the local copy is deleted
# before the next tar). Log lines per tar: "TAR name bytes download_s upload_s" or "FAILED name"; then "done".
#   remote_drive_sync.sh BUCKET PAR TAR...
B=$1 PAR=$2; shift 2
D=/content/sync
cd $D || exit 1
exec >> $D/sync.log 2>&1
one() {
  local t=$1 try t0 t1 t2 sz
  # an expired token (it lives as long as the workstation's gcloud token, 5 to 60 min; drive_sync.sh sends a new one
  # every 5 min): 20 tries a minute apart; another error or a failed Drive upload: 3 tries
  local fails=0
  for try in $(seq 20); do
    t0=$(date +%s)
    if ! gcloud storage --access-token-file=$D/gcs.tok cp "gs://$B/$t" "$D/stage/$t" > "$D/cp_$t.log" 2>&1; then
      rm -f "$D/stage/$t"; grep -q "HTTPError 401" "$D/cp_$t.log" || fails=$((fails + 1))
      [ $fails -lt 3 ] || break; sleep 60; continue
    fi
    t1=$(date +%s)
    rclone --config $D/rclone.conf copyto "$D/stage/$t" "cppwrf-drive:cppwrf-data/$t" --drive-chunk-size 256M \
      --retries 5 > "$D/rc_$t.log" 2>&1 || { rm -f "$D/stage/$t"; fails=$((fails + 1)); [ $fails -lt 3 ] || break; sleep 60; continue; }
    t2=$(date +%s); sz=$(stat -c %s "$D/stage/$t"); rm -f "$D/stage/$t"
    echo "TAR $t $sz $((t1 - t0)) $((t2 - t1))"; return 0
  done
  echo "FAILED $t"
}
mkdir -p $D/stage
for t in "$@"; do
  while [ "$(jobs -rp | wc -l)" -ge "$PAR" ]; do sleep 10; done
  one "$t" &
done
wait
echo done
