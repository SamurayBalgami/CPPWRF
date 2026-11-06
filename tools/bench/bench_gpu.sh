#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The GPU side of tools/bench/colab_bench.sh, on the Colab runtime itself (its Ubuntu with the CUDA toolkit; the
# Debian 13 chroot of tools/colab/remote_chroot_job.sh has no CUDA). Started next to the chroot job:
#
#   bench_gpu.sh STAMP
#
# Builds CPPWRF's CUDA FP32 (cuda) and FP64 (cuda_dp) binaries from /content/cpu/src.tar.gz (the release build of
# README.md for this GPU), waits until the chroot job (tools/bench/bench_body.sh) has written STAMP into
# /opt/wrfbench/gpu/body_started, links /opt/job and /opt/wrfbench to the chroot's (the run directories and
# configurations name chroot paths), writes ready ("cuda=EXIT cuda_dp=EXIT") and then runs each queue the chroot's
# tools/bench/bench_run.py hands over (queue_CASE.txt and go_CASE): per line "TAG BUILD RUNDIR INI OUTDIR REPS" REPS
# runs of `cppwrf INI --set output.dir=OUTDIR` in RUNDIR with CPPWRF_MEMINFO=1 (log RUNDIR/TAG.repN.log) while
# nvidia-smi samples the used device memory every 0.2 s; one JSON line per run in res_CASE.txt (exit, wall seconds,
# used MiB before the run and at its peak), then done_CASE. Ends with the chroot job (stop, or /content/cpu/done).
STAMP=$1
G=/content/bench_gpu R=/content/deb13 Q=/content/deb13/opt/wrfbench/gpu
mkdir -p $G
exec > $G/gpu.log 2>&1
set -x
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader > $G/gpu.txt
nvcc --version | tail -2 >> $G/gpu.txt
grep -m1 'model name' /proc/cpuinfo >> $G/gpu.txt
for i in $(seq 20); do  # the chroot job installs debootstrap at the same time: wait for apt's locks
  apt-get -o DPkg::Lock::Timeout=600 -qq update > $G/apt.log 2>&1 &&
    apt-get -o DPkg::Lock::Timeout=600 -qq install -y libnetcdf-dev >> $G/apt.log 2>&1 && break
  sleep 30
done
rm -rf $G/src && mkdir -p $G/src && tar xzf /content/cpu/src.tar.gz -C $G/src
st=
for b in cuda:OFF cuda_dp:ON; do
  n=${b%%:*}
  t0=$SECONDS
  cmake -S $G/src -B $G/build-$n -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=ON -DCPPWRF_DOUBLE=${b#*:} \
    > $G/build-$n.log 2>&1 && cmake --build $G/build-$n -j"$(nproc)" --target cppwrf >> $G/build-$n.log 2>&1
  st="$st $n=$? (${n}_build_s=$((SECONDS - t0)))"
done
echo "$st" >> $G/gpu.txt
until [ "$(cat $Q/body_started 2> /dev/null)" = "$STAMP" ]; do
  [ -e /content/cpu/done ] && exit 0
  sleep 10
done
ln -sfn $R/opt/job /opt/job
ln -sfn $R/opt/wrfbench /opt/wrfbench
cp $G/gpu.txt $Q/gpu.txt
for f in $G/build-*.log $G/apt.log; do tail -40 "$f" > "$Q/$(basename "$f")"; done
echo "$st" > $Q/ready
run_queue() {  # $1: case
  local c=$1 tag build rundir ini out reps k t0 t1 rc idle peak smi
  while read -r -u 3 tag build rundir ini out reps; do
    [ -n "$tag" ] || continue
    for k in $(seq "$reps"); do
      rm -rf "$out"
      idle=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
      nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -lms 200 > "$rundir/$tag.rep$k.mem" &
      smi=$!
      t0=$(date +%s.%N)
      (cd "$rundir" && CPPWRF_MEMINFO=1 OMP_NUM_THREADS=4 timeout 3600 "$G/build-$build/cppwrf" "$ini" \
        --set output.dir="$out" > "$rundir/$tag.rep$k.log" 2>&1)
      rc=$?
      t1=$(date +%s.%N)
      kill $smi; wait $smi 2> /dev/null
      peak=$(grep -E '^[0-9]+$' "$rundir/$tag.rep$k.mem" | sort -n | tail -1)
      echo "{\"tag\": \"$tag\", \"rep\": $k, \"exit\": $rc, \"wall\": $(awk "BEGIN {print $t1 - $t0}"), \"idle\": ${idle:-null}, \"peak\": ${peak:-null}}" >> "$Q/res_$c.txt"
    done
  done 3< "$Q/queue_$c.txt"
  touch "$Q/done_$c"
}
end=$((SECONDS + 8 * 3600))
while [ $SECONDS -lt $end ]; do
  for go in "$Q"/go_*; do
    [ -e "$go" ] || continue
    c=${go##*/go_}
    [ -e "$Q/done_$c" ] || run_queue "$c"
  done
  [ -e "$Q/stop" ] || [ -e /content/cpu/done ] && break
  sleep 5
done
cp $G/gpu.log $Q/gpu_host.log 2> /dev/null
touch $Q/gpu_exit
exit 0
