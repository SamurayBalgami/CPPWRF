#!/usr/bin/env bash
# The Colab side of tools/colab/colab_suite.sh and colab_wrfref.sh (root, under nohup; /content/cpu/done marks the
# end): a job in a Debian 13 chroot, so that gcc, gfortran, glibc (libmvec), NetCDF and Python are those of the GCP VMs
# that made the references (the Debian 13 packages of tools/gcp/remote_suite.sh; the runtime's own Ubuntu is not used).
#
#   remote_chroot_job.sh BUCKET DATA ARG...
#
# /content/cpu holds src.tar.gz (git archive), in/ (in/body.sh, the job, and its input files), rclone.conf (the Drive
# remote cppwrf-drive, tools/colab/drive_sync.sh) and gcs.tok (a read-only token for the bucket,
# tools/colab/gcs_token.py, renewed by the caller). Per session (cached for the next run of a --keep session): the
# chroot /content/deb13 and the data tars in /content/tars. A tar comes from Drive (cppwrf-drive:cppwrf-data), else
# from the bucket, and is then put on Drive for the next runs. The tars are unpacked into /opt/job/data in DATA's order
# (a later one replaces an earlier one's files), in/ is copied to /opt/job, and `bash /opt/job/body.sh /opt/job ARG...`
# runs in the chroot, as the GCP VMs' jobs in /opt/job. Logs: /content/cpu/job.log, versions.txt (packages, CPU,
# RAM, disk), mem.log (used RAM and the largest processes every minute), data.txt ("TAR name bytes seconds source"),
# setup.out (not Intel, a missing tar, a failed setup; the job's own setup.out is appended); the job's *.log, *.out
# and *.txt and these in /content/cpu/results.tar.gz; the job's /opt/job/fetch/ (large results) moves to
# /content/cpu/fetch/, the files of /opt/job/drive/ go to the Drive folder (data.txt: "DRIVE name bytes").
# gcs.tok is deleted when the data are in place, rclone.conf at the end.
B=$1 DATA=$2; shift 2
C=/content/cpu R=/content/deb13 T=/content/tars
J=$R/opt/job
cd $C || exit 1
exec > $C/job.log 2>&1
rm -rf $C/done $C/data_done $C/results.tar.gz $C/fetch $C/res
: > $C/setup.out; : > $C/data.txt
finish() {
  rm -f $C/gcs.tok $C/rclone.conf
  mkdir -p $C/res
  ( cd $J 2> /dev/null && cp -f ./*.log ./*.out ./*.txt $C/res/ 2> /dev/null )
  cat $C/res/setup.out >> $C/setup.out 2> /dev/null; rm -f $C/res/setup.out
  [ -d $J/fetch ] && mv $J/fetch $C/fetch
  cp -f $C/*.log $C/*.out $C/*.txt $C/res/ 2> /dev/null
  tar czf $C/results.tar.gz -C $C/res .
  date > $C/done
}
trap finish EXIT
set -x
# the CPU: the references are Intel-made, and glibc's vector log / atan (rcpps) differ on AMD
grep -q GenuineIntel /proc/cpuinfo || { echo "NOT_INTEL $(grep -m1 'model name' /proc/cpuinfo)" >> $C/setup.out; exit 1; }
export DEBIAN_FRONTEND=noninteractive
command -v rclone > /dev/null && command -v debootstrap > /dev/null ||
  { apt-get -qq update && apt-get -qq install -y debootstrap rclone > /dev/null; } || { echo "SETUP_FAILED apt" >> $C/setup.out; exit 1; }
mnt() { for d in proc sys dev dev/pts; do mountpoint -q $R/$d || mount --bind /$d $R/$d; done; }
if [ ! -f $R/.ready ]; then  # once per session
  for d in dev/pts dev sys proc; do umount $R/$d 2> /dev/null; done
  rm -rf $R
  debootstrap trixie $R http://deb.debian.org/debian > $C/debootstrap.log 2>&1 || { echo "SETUP_FAILED debootstrap" >> $C/setup.out; exit 1; }
  mnt
  cp /etc/resolv.conf $R/etc/resolv.conf
  chroot $R env -i PATH=/usr/sbin:/usr/bin:/sbin:/bin DEBIAN_FRONTEND=noninteractive bash -c 'apt-get -qq update &&
    apt-get -qq install -y build-essential cmake gfortran libnetcdf-dev libnetcdff-dev python3-netcdf4 python3-numpy > /dev/null' \
    > $C/packages.log 2>&1 || { echo "SETUP_FAILED packages (packages.log)" >> $C/setup.out; exit 1; }
  touch $R/.ready
fi
mnt
chroot $R dpkg-query -W -f='${Package} ${Version}\n' libc6 libgcc-s1 gcc-14 g++-14 gfortran-14 cmake libnetcdf22 libnetcdff7 \
  python3.13 python3-numpy python3-netcdf4 > $C/versions.txt 2>&1
lscpu | grep -E '^(Vendor ID|Model name|Model|Stepping|CPU\(s\)):' >> $C/versions.txt
free -g | awk 'NR == 2 {print "RAM", $2, "GB"}' >> $C/versions.txt
rm -rf $J; mkdir -p $J/src $J/data $T
tar xzf $C/src.tar.gz -C $J/src || { echo "SRC_FAILED" >> $C/setup.out; exit 1; }
cp -r $C/in/. $J/
# the tars: three fetched at a time (Drive, else the bucket), unpacked one by one in the list's order
gcs() {  # an expired token (it lives as long as the workstation's gcloud token, 5 to 60 min): the caller sends a
  local i  # new one every 5 min, so retry for 20 min
  for i in $(seq 20); do
    gcloud storage --access-token-file=$C/gcs.tok cp "gs://$B/$1" "$T/$1" > "$C/try_$1.log" 2>&1; local r=$?
    cat "$C/try_$1.log" >> "$C/fetch_$1.log"; [ $r = 0 ] && return 0
    grep -q "HTTPError 401" "$C/try_$1.log" || return 1
    rm -f "$T/$1"; sleep 60
  done
  return 1
}
fetch() {
  local t=$1 t0=$SECONDS src=cache
  if [ ! -f "$T/$t.ok" ]; then
    rm -f "$T/$t"
    if rclone --config $C/rclone.conf copyto "cppwrf-drive:cppwrf-data/$t" "$T/$t" > "$C/fetch_$t.log" 2>&1 && [ -s "$T/$t" ]; then src=drive
    elif rm -f "$T/$t"; gcs "$t"; then
      src=bucket  # for the next runs: on Drive too, while the job runs (rclone.conf is deleted after it)
      touch "$T/$t.todrive"
    else rm -f "$T/$t"; touch "$T/$t.failed"; return; fi
    touch "$T/$t.ok"
  fi
  echo "TAR $t $(stat -c %s "$T/$t") $((SECONDS - t0)) $src" >> $C/data.txt
}
rm -f $T/*.failed
for t in ${DATA//,/ }; do
  while [ "$(jobs -rp | wc -l)" -ge 3 ]; do sleep 5; done
  fetch "$t" &
done
for t in ${DATA//,/ }; do
  while [ ! -e "$T/$t.ok" ] && [ ! -e "$T/$t.failed" ]; do sleep 5; done
  if [ -e "$T/$t.ok" ]; then tar xf "$T/$t" -C $J/data || echo "DATA_FAILED $t (unpack)" >> $C/setup.out
  else echo "DATA_FAILED $t" >> $C/setup.out; fi
done
wait
rm -f $C/gcs.tok
touch $C/data_done
for f in $T/*.todrive; do  # the tars that came from the bucket, to Drive next to the job
  [ -e "$f" ] || continue
  t=$(basename "$f" .todrive)
  rclone --config $C/rclone.conf copyto "$T/$t" "cppwrf-drive:cppwrf-data/$t" --drive-chunk-size 256M >> "$C/fetch_$t.log" 2>&1 &&
    { rm -f "$f"; echo "DRIVE $t $(stat -c %s "$T/$t")" >> $C/data.txt; } || echo "DRIVE_COPY_FAILED $t" >> "$C/fetch_$t.log"
done > /dev/null 2>&1 &
TODRIVE=$!
df -h / | awk 'NR == 2 {print "disk", $2, "used", $3}' >> $C/versions.txt
# memory every minute (the runtime has no swap: the kernel kills the largest process when RAM runs out)
(set +x; while sleep 60; do echo "$(date +%H:%M) $(free -m | awk 'NR == 2 {print $3 "/" $2 " MB"}')" \
   "$(ps -eo rss,comm --sort=-rss | awk 'NR > 1 && NR < 5 {printf "%s %d MB, ", $2, $1 / 1024}')"; done) > $C/mem.log 2>&1 &
MEM=$!
chroot $R env -i PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin HOME=/root \
  bash /opt/job/body.sh /opt/job "$@" > $C/body.log 2>&1
kill $MEM
wait $TODRIVE
for f in $J/drive/*; do  # the job's tars for the next runs (wrfref_body.sh STAGE)
  [ -f "$f" ] || continue
  if rclone --config $C/rclone.conf copyto "$f" "cppwrf-drive:cppwrf-data/$(basename "$f")" --drive-chunk-size 256M > $C/drive.log 2>&1
  then echo "DRIVE $(basename "$f") $(stat -c %s "$f")" >> $C/data.txt; else echo "DRIVE_FAILED $(basename "$f") (drive.log)" >> $C/setup.out; fi
done
