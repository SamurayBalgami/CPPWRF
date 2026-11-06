#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The GPU side of tools/colab/colab_quickstart.sh, on the Colab runtime itself (its Ubuntu with the CUDA toolkit; the
# Debian 13 chroot has no CUDA). Started next to the chroot job:
#
#   quickstart_gpu.sh STAMP
#
# Builds README.md's CUDA FP32 (build) and FP64 (build-gpu-dp) binaries from /content/cpu/src.tar.gz, runs the Straka
# case in each, waits until the chroot job (tools/colab/quickstart_body.sh) has written STAMP into
# /opt/qs/gpu/cpu_ready (the example's inputs are made), links /opt/qs to the chroot's (run.ini names chroot paths)
# and runs the example's run.ini in both builds (outputs /opt/qs/gpu/out_cuda, out_cuda_dp). One line per step in
# /opt/qs/gpu/gpu.txt ("PASS|FAIL step (seconds)"), then gpu_done.
STAMP=$1
H=/content/qs_gpu R=/content/deb13 G=/content/deb13/opt/qs/gpu
mkdir -p $H
exec > $H/gpu_host.log 2>&1
set -x
T=$H/gpu.txt; : > $T
step() {  # NAME COMMAND...
  local n=$1 t0=$SECONDS rc; shift
  "$@" > $H/gpu_$n.log 2>&1; rc=$?
  echo "$([ $rc = 0 ] && echo PASS || echo FAIL) $n ($((SECONDS - t0)) s)$([ $rc = 0 ] || echo " exit $rc, gpu_$n.log")" >> $T
  return $rc
}
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader > $H/gpu_info.log
nvcc --version | tail -2 >> $H/gpu_info.log
for i in $(seq 20); do  # the chroot job installs debootstrap at the same time: wait for apt's locks
  apt-get -o DPkg::Lock::Timeout=600 -qq update > $H/apt.log 2>&1 &&
    apt-get -o DPkg::Lock::Timeout=600 -qq install -y libnetcdf-dev >> $H/apt.log 2>&1 && break
  sleep 30
done
rm -rf $H/src && mkdir -p $H/src && tar xzf /content/cpu/src.tar.gz -C $H/src
cd $H/src || exit 1
NP=$(nproc)
# README.md, "Build": the GPU builds (all targets, as the README's commands)
step build_cuda_fp32 bash -c "cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=ON -DCPPWRF_DOUBLE=OFF &&
  cmake --build build -j $NP"
step build_cuda_fp64 bash -c "cmake -S . -B build-gpu-dp -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=ON -DCPPWRF_DOUBLE=ON &&
  cmake --build build-gpu-dp -j $NP"
step straka_cuda_fp32 build/cppwrf configs/straka.ini --set output.dir=$H/straka_cuda
step straka_cuda_fp64 build-gpu-dp/cppwrf configs/straka.ini --set output.dir=$H/straka_cuda_dp
end=$((SECONDS + 3 * 3600))
until [ "$(cat $G/cpu_ready 2> /dev/null)" = "$STAMP" ]; do
  [ -e /content/cpu/done ] || [ $SECONDS -ge $end ] && exit 0
  sleep 10
done
ln -sfn $R/opt/qs /opt/qs
cd /opt/qs/ex/run || exit 1
if [ -f run.ini ]; then
  step example_cuda_fp32 env OMP_NUM_THREADS=4 $H/src/build/cppwrf run.ini --set output.dir=$G/out_cuda
  step example_cuda_fp64 env OMP_NUM_THREADS=4 $H/src/build-gpu-dp/cppwrf run.ini --set output.dir=$G/out_cuda_dp
else echo "FAIL example_cuda (no run.ini: the CPU example failed)" >> $T; fi
for f in $H/gpu_*.log $H/apt.log; do tail -60 "$f" > "$G/$(basename "$f")"; done
grep -h "cppwrf: done" $H/gpu_example_*.log >> $G/gpu_host.log 2> /dev/null
cp $T $G/gpu.txt
touch $G/gpu_done
