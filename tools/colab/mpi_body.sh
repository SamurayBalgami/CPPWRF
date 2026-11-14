#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The chroot job of tools/colab/colab_mpi.sh: the CPU FP64 build with CPPWRF_MPI and suite cases run under mpirun
# with several rank layouts (docs/MPI.md), each against the WRF references at tolerance 0; the output files of every
# layout are also compared byte for byte with those of the run with one rank.
#
#   mpi_body.sh J CASES LAYOUTS JOBS MPI
#
# J holds src/ and data/ (tools/colab/remote_chroot_job.sh). LAYOUTS: PxQ,... (1: one rank). MPI: mpich or openmpi
# (the Debian packages, installed here). The build lives in /opt/mpi (kept by a --keep session). Logs in J:
# setup.out, cmake.log, mpi.txt (the MPI's version, the smoke test), suite_<layout>.log, sums_<layout>.txt (sha256
# of every output file), same.txt ("SAME layout files" / "DIFF layout: files").
J=$1 CASES=$2 LAYOUTS=$3 JOBS=$4 MPI=${5:-mpich} GPU=${6:-0}
# GPU=1: tools/colab/mpi_gpu.sh runs on the runtime next to this job, which ends when it has written gpu_done.
wait_gpu() {
  [ "$GPU" = 1 ] || return 0
  local end=$((SECONDS + 6 * 3600))
  until [ -e $J/gpu_done ] || [ $SECONDS -gt $end ]; do sleep 20; done
}
[ -n "$CASES" ] || { touch $J/setup.out; wait_gpu; exit 0; }
P=/opt/mpi
set -x
touch $J/setup.out
case $MPI in
  mpich) PK="mpich libmpich-dev"; MPIRUN="mpirun.mpich"; FLAGS= ;;
  openmpi) PK="openmpi-bin libopenmpi-dev"; MPIRUN="mpirun.openmpi"; FLAGS="--allow-run-as-root --oversubscribe --bind-to none" ;;
  *) echo "MPI_UNKNOWN $MPI" >> $J/setup.out; exit 1 ;;
esac
command -v $MPIRUN > /dev/null || DEBIAN_FRONTEND=noninteractive apt-get -qq install -y $PK > $J/mpi_packages.log 2>&1 ||
  { echo "SETUP_FAILED $PK (mpi_packages.log)" >> $J/setup.out; exit 1; }
{ $MPIRUN --version 2>&1 | head -5; echo "smoke: $($MPIRUN $FLAGS -np 2 sh -c 'echo rank' 2>&1 | tr '\n' ' ')"; } > $J/mpi.txt
# The transport: JOBS runs of 8 ranks at once must start (MPICH's UCX device ran out of shared-memory resources
# in the chroot before /dev/shm was mounted). The first setting that passes is used for every layout.
mkdir -p $P
# the chroot's /dev is the container's (64 MB): MPICH's shared memory needs a larger /dev/shm
mountpoint -q /dev/shm || mount -t tmpfs -o size=32g tmpfs /dev/shm || echo "no tmpfs on /dev/shm" >> $J/mpi.txt
cat > $P/hello.c <<'C'
#include <mpi.h>
#include <stdio.h>
int main(int argc, char** argv) {
  int r, n, s = 0;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &r);
  MPI_Comm_size(MPI_COMM_WORLD, &n);
  MPI_Allreduce(&r, &s, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (r == 0) printf("ok %d %d\n", n, s);
  MPI_Finalize();
  return 0;
}
C
case $MPI in mpich) mpicc.mpich -O1 -o $P/hello $P/hello.c ;; *) mpicc.openmpi -O1 -o $P/hello $P/hello.c ;; esac >> $J/mpi.txt 2>&1
{ df -h /dev/shm 2>&1 | tail -1; ulimit -a 2>&1 | grep -E 'locked|open files|processes'; } >> $J/mpi.txt
MPIENV=
for E in "" "UCX_TLS=tcp,self" "UCX_TLS=posix,self" "UCX_TLS=sysv,self"; do
  ok=0
  for q in $(seq 1 $JOBS); do env $E $MPIRUN $FLAGS -np 8 $P/hello > $P/hello_$q.out 2>&1 & done
  wait
  for q in $(seq 1 $JOBS); do grep -q '^ok 8 28$' $P/hello_$q.out && ok=$((ok + 1)); done
  echo "transport '${E:-default}': $ok of $JOBS runs of 8 ranks ($(grep -h 'init_worker' $P/hello_*.out | head -1 | cut -c1-160))" >> $J/mpi.txt
  if [ $ok = "$JOBS" ]; then MPIENV=$E; echo "transport used: '${E:-default}'" >> $J/mpi.txt; break; fi
  E=none
done
[ "$E" != none ] || { echo "SETUP_FAILED no MPI transport starts $JOBS runs of 8 ranks (mpi.txt)" >> $J/setup.out; exit 1; }
mkdir -p $P/src
python3 - $J/src $P/src <<'PY'
import filecmp, os, shutil, sys
src, dst = sys.argv[1], sys.argv[2]
keep = set()
for root, dirs, files in os.walk(src):
    rel = os.path.relpath(root, src)
    os.makedirs(os.path.join(dst, rel), exist_ok=True)
    for f in files:
        a, b = os.path.join(root, f), os.path.join(dst, rel, f)
        keep.add(os.path.normpath(b))
        if not os.path.exists(b) or not filecmp.cmp(a, b, shallow=False):
            shutil.copy(a, b)  # a new time stamp: rebuilt
for root, dirs, files in os.walk(dst):
    for f in files:
        b = os.path.normpath(os.path.join(root, f))
        if b not in keep:
            os.remove(b)
PY
[ -f $P/build/Makefile ] || cmake -S $P/src -B $P/build -DCMAKE_BUILD_TYPE=Release -DCPPWRF_DOUBLE=ON \
  -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_NATIVE=OFF -DCPPWRF_MPI=ON -DCMAKE_CXX_FLAGS="-march=skylake -mtune=skylake" \
  > $J/cmake.log 2>&1
cmake --build $P/build -j"$(nproc)" --target cppwrf >> $J/cmake.log 2>&1 ||
  { echo "BUILD_FAILED (cmake.log)" >> $J/setup.out; wait_gpu; exit 1; }
cd $P/src
export CPPWRF_WRF_RUNS=$J/data/wrf_runs CPPWRF_WRF_REF=$J/data/wrf_ref
NP=$(nproc)
: > $J/same.txt
for L in ${LAYOUTS//,/ }; do
  if [ "$L" = 1 ]; then N=1; else N=$(( ${L%x*} * ${L#*x} )); fi
  jobs=$(( NP / N < JOBS ? NP / N : JOBS )); [ $jobs -ge 1 ] || jobs=1
  thr=$(( NP / (N * jobs) )); [ $thr -ge 1 ] || thr=1
  # the suite's binary: the layout under mpirun (one rank: the plain program, also under mpirun)
  W=$P/run_$L
  { echo '#!/bin/sh'
    [ -z "$MPIENV" ] || echo "export $MPIENV"
    if [ "$L" = 1 ]; then echo "exec $MPIRUN $FLAGS -np 1 $P/build/cppwrf \"\$@\""
    else
      # a 2-D case has one patch across its plane: its ranks all along x (x-z) or y (y-z)
      echo "R=$L; case \"\$*\" in *2d_x*) R=${N}x1 ;; *2d_y*) R=1x${N} ;; esac"
      echo "CPPWRF_RANKS=\$R exec $MPIRUN $FLAGS -np $N $P/build/cppwrf \"\$@\""
    fi; } > $W
  chmod +x $W
  t0=$SECONDS
  OMP_NUM_THREADS=$thr python3 tests/verify/wrf_suite.py --cpu-dp $W --jobs $jobs --out $J/s_$L --only "$CASES" --exact \
    --keep > $J/suite_$L.log 2>&1
  echo "exit $? ($((SECONDS - t0)) s, $N ranks, $jobs jobs, $thr threads)" >> $J/suite_$L.log
  # every output file but the logs, the files that hold wall-clock times and the drop-in's run.ini (it names the
  # run's own directory)
  (cd $J/s_$L 2> /dev/null && find . -type f ! -name '*.log' ! -name 'run.json' ! -name 'run.ini' ! -name 'diag.csv' ! -name 'rsl.*' \
    ! -name 'timing*' -print0 | sort -z | xargs -0 -r sha256sum) > $J/sums_$L.txt
  rm -rf $J/s_$L*
  if [ "$L" != 1 ] && [ -s $J/sums_1.txt ]; then
    d=$(diff <(awk '{print $2, $1}' $J/sums_1.txt | sort) <(awk '{print $2, $1}' $J/sums_$L.txt | sort) | grep -c '^[<>]')
    if [ "$d" = 0 ] && [ -s $J/sums_$L.txt ]; then echo "SAME $L $(wc -l < $J/sums_$L.txt) files" >> $J/same.txt
    else echo "DIFF $L: $(diff <(awk '{print $2, $1}' $J/sums_1.txt | sort) <(awk '{print $2, $1}' $J/sums_$L.txt | sort) |
      grep '^[<>]' | awk '{print $2}' | sort -u | head -20 | tr '\n' ' ')" >> $J/same.txt; fi
  fi
done
wait_gpu
exit 0
