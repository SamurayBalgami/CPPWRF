#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The GPU side of tools/bench/colab_fp32.sh, on the Colab runtime itself (its Ubuntu with the CUDA toolkit), started
# next to the chroot job (tools/bench/fp32_body.sh):
#
#   fp32_gpu.sh STAMP [SUITE]   (SUITE 0: no suite cases; FP32_EXTRA_BUILDS=1: also cuda_nofma and cpu)
#
# Builds CPPWRF's CUDA FP32 (cuda) and FP64 (cuda_dp) release binaries from /content/cpu/src.tar.gz, waits until the
# chroot job has prepared the cases (STAMP in /opt/fp1/prepared), links /opt/fp1, /opt/job and /opt/wrfbench to the
# chroot's, runs the GPU queue (tools/bench/fp32_run.py gpu), then the suite's tolerance-0 cuda_dp cases of the
# studied configurations (tests/verify/wrf_suite.py --cuda-dp: CPPWRF FP64 against the WRF real*8 references) and the
# FP32 cases whose tolerances this study sets, from a tar made by tools/colab/stage_data.py --only (as a GPU proof
# stages them); after the WRF runs (/opt/fp1/wrf_done) the statistics (fp32_run.py stats), then /opt/fp1/host_done.
STAMP=$1 SUITE=${2:-1}
G=/content/fp1_gpu R=/content/deb13 F=/content/deb13/opt/fp1
rm -rf $G && mkdir -p $G
exec > $G/gpu.log 2>&1
set -x
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader > $G/gpu.txt
nvcc --version | tail -2 >> $G/gpu.txt
grep -m1 'model name' /proc/cpuinfo >> $G/gpu.txt
for i in $(seq 20); do
  apt-get -o DPkg::Lock::Timeout=600 -qq update > $G/apt.log 2>&1 &&
    apt-get -o DPkg::Lock::Timeout=600 -qq install -y libnetcdf-dev >> $G/apt.log 2>&1 && break
  sleep 30
done
pip -q install netCDF4 > $G/pip.log 2>&1
mkdir -p $G/src && tar xzf /content/cpu/src.tar.gz -C $G/src
st=
# cuda (FP32), cuda_dp (FP64); with FP32_EXTRA_BUILDS=1 (the FP2 cases): cuda_nofma (FP32, no FMA contraction in any
# kernel) and cpu (CPU FP32: glibc's libm, no FMA)
B="cuda:ON:OFF cuda_dp:ON:ON"
[ "${FP32_EXTRA_BUILDS:-0}" = 1 ] && B="$B cuda_nofma:ON:OFF cpu:OFF:OFF"
for b in $B; do
  n=${b%%:*} r=${b#*:}
  t0=$SECONDS
  x=; [ "$n" = cuda_nofma ] && x=--fmad=false
  cmake -S $G/src -B $G/build-$n -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=${r%%:*} -DCPPWRF_DOUBLE=${r#*:} \
    ${x:+-DCMAKE_CUDA_FLAGS=$x} > $G/build-$n.log 2>&1 && cmake --build $G/build-$n -j"$(nproc)" --target cppwrf >> $G/build-$n.log 2>&1
  st="$st $n=$? (${n}_build_s=$((SECONDS - t0)))"
done
echo "builds:$st" >> $G/gpu.txt
until [ "$(cat $F/prepared 2> /dev/null)" = "$STAMP" ]; do
  [ -e /content/cpu/done ] && exit 0
  sleep 10
done
trap 'touch $F/host_done' EXIT  # the chroot job waits for it
ln -sfn $R/opt/fp1 /opt/fp1
ln -sfn $R/opt/job /opt/job
ln -sfn $R/opt/wrfbench /opt/wrfbench
O=/opt/fp1/out
cp $G/gpu.txt $O/gpu.txt
for f in $G/build-*.log; do tail -40 "$f" > "$O/$(basename "$f")"; done
echo "$(date +%H:%M:%S) GPU side: builds$st" >> $O/fp1.log
cd $G/src || exit 1
python3 tools/bench/fp32_run.py gpu $G
D=/opt/job/data
if [ "$SUITE" = 1 ]; then
# CPPWRF FP64 = WRF real*8 at tolerance 0 (every variable of every compared frame) for the studied configurations
CPPWRF_WRF_RUNS=$D/wrf_runs CPPWRF_WRF_REF=$D/wrf_ref OMP_NUM_THREADS=4 timeout 7200 python3 tests/verify/wrf_suite.py \
  --cuda-dp $G/build-cuda_dp/cppwrf --exact --allow-skip --out $G/suite_dp \
  --only real_p4_r8_fp64_cuda_1h,real_p6_r8_fp64_cuda_1h,real_p9_r8_fp64_cuda_1h,real_v1_r8_fp64_cuda_2h,real_t1_fp64_cuda,real_n3_r8_fp64_cuda_1h,real_skv_r8_fp64_cuda_1h,real_morr_r8_fp64_cuda_1h,real_acm7_r8_fp64_cuda_1h \
  > $O/suite_dp.log 2>&1
echo "$(date +%H:%M:%S) GPU side: suite cuda_dp exit $?: $(tail -1 $O/suite_dp.log)" >> $O/fp1.log
rm -rf $G/suite_dp
# the FP32 cases whose tolerances come from this study, from their staged references (tools/colab/stage_data.py)
FP32C=real_morr_fp32_1h,real_acm7_fp32_1h,real_skv_fp32_1h,real_p6_fp32_1h,real_p9_fp32_1h
CPPWRF_WRF_RUNS=$D/wrf_runs CPPWRF_WRF_REF=$D/wrf_ref python3 tools/colab/stage_data.py --only "$FP32C" \
  --no-thompson-tables $G/fp32stage.tar > $O/stage_fp32.log 2>&1
echo "stage exit $?" >> $O/stage_fp32.log
rm -rf $G/fp32data && mkdir -p $G/fp32data && tar xf $G/fp32stage.tar -C $G/fp32data && rm -f $G/fp32stage.tar
CPPWRF_WRF_RUNS=$G/fp32data/wrf_runs CPPWRF_WRF_REF=$G/fp32data/wrf_ref OMP_NUM_THREADS=4 timeout 3600 python3 \
  tests/verify/wrf_suite.py --cuda $G/build-cuda/cppwrf --only "$FP32C" --allow-skip --out $G/suite_sp > $O/suite_sp.log 2>&1
echo "$(date +%H:%M:%S) GPU side: suite FP32 exit $?: $(tail -1 $O/suite_sp.log)" >> $O/fp1.log
rm -rf $G/suite_sp $G/fp32data
fi
until [ -e /opt/fp1/wrf_done ]; do sleep 20; done
python3 tools/bench/fp32_run.py stats
cp $G/gpu.log $O/gpu_host.log 2> /dev/null
exit 0
