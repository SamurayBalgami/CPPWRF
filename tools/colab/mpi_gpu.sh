#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The GPU side of tools/colab/colab_mpi.sh on the Colab runtime itself (its Ubuntu with the CUDA toolkit), next to
# the chroot job (tools/colab/mpi_body.sh), whose data it uses (/opt/job is linked to the chroot's):
#
#   mpi_gpu.sh CASES LAYOUTS
#
# Builds CUDA FP64 with CPPWRF_MPI (the runtime's Open MPI package) from /content/cpu/src.tar.gz and runs CASES (cuda_dp
# suite cases) under mpirun with every layout of LAYOUTS (1: one rank), all ranks on the one GPU, against the WRF
# references (tolerance 0): gpu_suite_<layout>.log; the output files of every layout are compared byte for byte
# with those of the run with one rank (gpu_same.txt). Logs go to /opt/job (the chroot job's results); gpu_done ends
# the chroot job's wait.
CASES=$1 LAYOUTS=$2 FP32=${3:-}
# FP32 (cuda_dp suite cases, run here in FP32 for at most 30 simulated minutes): the CUDA FP32 build of the base
# revision (/content/cpu/in/base_src.tar.gz, without MPI) and this revision's with CPPWRF_MPI as one process and as
# 2 x 1 and 1 x 2 ranks on the GPU; every output file must be the base build's (gpu_fp32.txt).
G=/content/mpi_gpu R=/content/deb13 J=/content/deb13/opt/job
mkdir -p $G
exec > $G/gpu.log 2>&1
set -x
until [ -e /content/cpu/data_done ] || [ -e /content/cpu/done ]; do sleep 10; done
[ -e /content/cpu/done ] && exit 0
ln -sfn $R/opt/job /opt/job
finish() { cp $G/gpu.log $J/gpu_host.log 2> /dev/null; touch $J/gpu_done; }
trap finish EXIT
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader > $J/gpu.txt
nvcc --version | tail -2 >> $J/gpu.txt
for i in $(seq 20); do  # the chroot job may hold apt's locks
  apt-get -o DPkg::Lock::Timeout=600 -qq update > $G/apt.log 2>&1 &&
    apt-get -o DPkg::Lock::Timeout=600 -qq install -y libnetcdf-dev openmpi-bin libopenmpi-dev >> $G/apt.log 2>&1 && break
  sleep 30
done
# Open MPI: the runtime's Ubuntu MPICH package (4.2.0, built for PMIx) starts every rank as a world of one under its mpirun
MPIRUN="$(command -v mpirun.openmpi || command -v mpirun) --allow-run-as-root --oversubscribe"
{ $MPIRUN --version 2>&1 | head -3; } >> $J/gpu.txt
python3 -c 'import netCDF4' 2> /dev/null || pip install -q netCDF4 > $G/pip.log 2>&1
rm -rf $G/src && mkdir -p $G/src
tar xzf /content/cpu/src.tar.gz -C $G/src
t0=$SECONDS
cmake -S $G/src -B $G/b_dp -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=ON -DCPPWRF_DOUBLE=ON -DCPPWRF_MPI=ON \
  -DMPI_CXX_COMPILER="$(command -v mpicxx.openmpi || command -v mpicxx)" > $G/b_dp.log 2>&1 &&
  cmake --build $G/b_dp -j"$(nproc)" --target cppwrf >> $G/b_dp.log 2>&1
echo "build: exit $?, $((SECONDS - t0)) s" >> $J/gpu.txt
tail -40 $G/b_dp.log > $J/gpu_build.log
if [ -n "$FP32" ]; then
  rm -rf $G/base && mkdir -p $G/base && tar xzf /content/cpu/in/base_src.tar.gz -C $G/base
  ( cmake -S $G/base -B $G/b_base -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=ON -DCPPWRF_DOUBLE=OFF > $G/b_base.log 2>&1 &&
      cmake --build $G/b_base -j"$(nproc)" --target cppwrf >> $G/b_base.log 2>&1
    echo "build base FP32: exit $?" >> $J/gpu.txt ) &
  ( cmake -S $G/src -B $G/b_sp -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=ON -DCPPWRF_DOUBLE=OFF -DCPPWRF_MPI=ON \
      -DMPI_CXX_COMPILER="$(command -v mpicxx.openmpi || command -v mpicxx)" > $G/b_sp.log 2>&1 &&
      cmake --build $G/b_sp -j"$(nproc)" --target cppwrf >> $G/b_sp.log 2>&1
    echo "build FP32 with MPI: exit $?" >> $J/gpu.txt ) &
  wait
fi
[ -x $G/b_dp/cppwrf ] || exit 1
export CPPWRF_WRF_RUNS=/opt/job/data/wrf_runs CPPWRF_WRF_REF=/opt/job/data/wrf_ref
cd $G/src
: > $J/gpu_same.txt
for L in ${LAYOUTS//,/ }; do
  if [ "$L" = 1 ]; then N=1; else N=$(( ${L%x*} * ${L#*x} )); fi
  W=$G/run_$L
  { echo '#!/bin/sh'
    if [ "$L" = 1 ]; then echo "exec $MPIRUN -np 1 $G/b_dp/cppwrf \"\$@\""
    else
      echo "R=$L; case \"\$*\" in *2d_x*) R=${N}x1 ;; *2d_y*) R=1x${N} ;; esac"  # a 2-D (x-z) case: its ranks along x
      echo "CPPWRF_RANKS=\$R exec $MPIRUN -np $N $G/b_dp/cppwrf \"\$@\""
    fi; } > $W
  chmod +x $W
  t0=$SECONDS
  OMP_NUM_THREADS=4 python3 tests/verify/wrf_suite.py --cuda-dp $W --only "$CASES" --exact --keep --out $G/s_$L \
    > $J/gpu_suite_$L.log 2>&1
  echo "exit $? ($((SECONDS - t0)) s, $N ranks on one GPU)" >> $J/gpu_suite_$L.log
  (cd $G/s_$L 2> /dev/null && find . -type f ! -name '*.log' ! -name 'run.json' ! -name 'run.ini' ! -name 'diag.csv' ! -name 'rsl.*' \
    ! -name 'timing*' -print0 | sort -z | xargs -0 -r sha256sum) > $J/gpu_sums_$L.txt
  grep -rh "cppwrf: done\|host memory of the ranks\|device memory" $G/s_$L --include='*.log' 2> /dev/null | sort | uniq -c |
    head -40 > $J/gpu_runs_$L.txt
  rm -rf $G/s_$L*
  if [ "$L" != 1 ] && [ -s $J/gpu_sums_1.txt ]; then
    d=$(diff <(awk '{print $2, $1}' $J/gpu_sums_1.txt | sort) <(awk '{print $2, $1}' $J/gpu_sums_$L.txt | sort) | grep -c '^[<>]')
    if [ "$d" = 0 ] && [ -s $J/gpu_sums_$L.txt ]; then echo "SAME $L $(wc -l < $J/gpu_sums_$L.txt) files" >> $J/gpu_same.txt
    else echo "DIFF $L: $(diff <(awk '{print $2, $1}' $J/gpu_sums_1.txt | sort) <(awk '{print $2, $1}' $J/gpu_sums_$L.txt | sort) |
      grep '^[<>]' | awk '{print $2}' | sort -u | head -20 | tr '\n' ' ')" >> $J/gpu_same.txt; fi
  fi
done
if [ -n "$FP32" ]; then
  : > $J/gpu_fp32.txt
  sums() { (cd "$1" 2> /dev/null && find . -type f ! -name '*.log' ! -name 'run.json' ! -name 'diag.csv' -print0 | sort -z |
    xargs -0 -r sha256sum); }
  for name in ${FP32//,/ }; do
    args=$(python3 - "$name" <<'EOF'
import shlex, sys
sys.path.insert(0, "tests/verify")
import wrf_suite
c = {x[0]: x for x in wrf_suite.CASES}[sys.argv[1]]
s = c[3] + [f"run.run_time={int(min(c[5], 1800))}", "output.format=both"]
print(shlex.join([c[2]] + [a for x in s for a in ("--set", x)]))
EOF
)
    for r in base:$G/b_base/cppwrf one:$G/b_sp/cppwrf 2x1:$G/b_sp/cppwrf 1x2:$G/b_sp/cppwrf; do
      tag=${r%%:*} bin=${r#*:} pre=
      case $tag in 2x1|1x2) pre="env CPPWRF_RANKS=$tag $MPIRUN -np 2" ;; esac
      case "$tag $args" in "1x2 "*2d_x*) pre="env CPPWRF_RANKS=2x1 $MPIRUN -np 2" ;; esac
      rm -rf $G/f_$tag
      eval "OMP_NUM_THREADS=4 $pre $bin $args --set output.dir=$G/f_$tag" > $G/f_$tag.log 2>&1
      rc=$?
      sums $G/f_$tag > $G/f_$tag.sums
      if [ "$tag" = base ]; then echo "RUN $name base exit $rc, $(wc -l < $G/f_base.sums) files" >> $J/gpu_fp32.txt
      elif [ $rc = 0 ] && [ -s $G/f_base.sums ] && cmp -s $G/f_base.sums $G/f_$tag.sums; then
        echo "SAME $name $tag $(wc -l < $G/f_$tag.sums) files" >> $J/gpu_fp32.txt
      else echo "DIFF $name $tag exit $rc: $(diff $G/f_base.sums $G/f_$tag.sums | grep '^[<>]' | awk '{print $3}' | sort -u | head -8 | tr '\n' ' ') $(tail -1 $G/f_$tag.log | cut -c1-200)" >> $J/gpu_fp32.txt; fi
    done
  done
fi
exit 0
