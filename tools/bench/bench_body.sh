#!/usr/bin/env bash
# The CPU side of tools/bench/colab_bench.sh in the Debian 13 chroot of tools/colab/remote_chroot_job.sh (gcc /
# gfortran 14, glibc 2.41, NetCDF of the GCP VMs that made the suite's references).
#
#   bench_body.sh J CASES LAYOUTS STAMP
#
# J: /opt/job (src/, data/wrf_runs/ from the data tars, wps.tar.gz). Builds, in /opt/wrfbench (kept for a later run in
# the same --keep session; /opt/job is made anew by every run):
# - WRF v4.8.0 from the official release source (github.com/wrf-model/WRF, downloaded here) with WRF's own configure
#   for GNU gfortran / gcc and MPICH: dmpar (option 34) and dm+sm (35), each real*4 and real*8 (configure -r8, the
#   suite's wrf_r8.exe): WRF-r4, WRF-r8, WRF-r4sm, WRF-r8sm, basic nesting, compile -j 12 em_real;
# - CPPWRF's CPU FP64 build (cpu-dp: README.md's release build, -march=native, CPPWRF_PORTABLE_LIBMVEC on) and the same
#   with CPPWRF_PORTABLE_LIBMVEC off (cpu-dp-hwlog: libmvec's log on this CPU, as the WRF builds here);
# - the large domain L (tools/bench/namelist.wps: 360 x 270 points at 15 km, 9.1 times P4's columns, same levels and
#   physics): WPS 4.6.0's geogrid / ungrib / metgrid (the executables that made the suite's met_em files, in
#   wps.tar.gz) on WPS's low-resolution mandatory geography (downloaded from UCAR here) and the GFS 0.25 degree analysis
#   and 6 h forecast of 2026-10-12 00 UTC (NOAA's open data, downloaded here), then real.exe of WRF-r8
#   (tools/bench/bench_run.py make-large).
# The CUDA builds and GPU runs are tools/bench/bench_gpu.sh's, outside the chroot (STAMP: this run's mark for it).
# Then tools/bench/bench_run.py runs each case of CASES (R1, P4, P4R15, N3, L) with the WRF layouts of LAYOUTS and
# writes J/results.txt, J/report.md (tools/bench/bench_report.py) and J/fetch/bench.tar.gz (the logs).
# Progress: J/bench.log. Nothing runs next to the timed runs (the builds and downloads end before them).
J=$1 CASES=$2 LAYOUTS=$3 STAMP=$4
W=/opt/wrfbench Q=/opt/wrfbench/gpu
mkdir -p $W $J/fetch/logs
touch $J/setup.out $J/results.txt
B() { echo "$(date +%H:%M:%S) $*" >> $J/bench.log; }
rm -rf $Q && mkdir -p $Q && echo "$STAMP" > $Q/body_started
for n in wrf.exe real.exe cppwrf mpiexec mpiexec.hydra hydra_pmi_proxy; do pkill -9 -x $n; done  # an earlier run's (--keep)
set -x
NP=$(nproc)
B "start: $NP vCPUs, cases $CASES, layouts $LAYOUTS"
export DEBIAN_FRONTEND=noninteractive
apt-get -qq update > $J/apt.log 2>&1 && apt-get -qq install -y mpich libmpich-dev tcsh m4 perl make curl ca-certificates \
  git file libpng16-16t64 > $J/apt.log 2>&1 || { echo "BENCH_SETUP_FAILED apt (apt.log)" >> $J/setup.out; B "apt failed"; }
lscpu > $J/lscpu.txt
{ gfortran --version | head -1; gcc --version | head -1; mpichversion 2> /dev/null | head -2; } > $J/compilers.txt 2>&1
# the chroot's /dev is the container's (64 MB): MPICH's shared memory (UCX) for 24 and more ranks needs a larger /dev/shm
mountpoint -q /dev/shm || mount -t tmpfs -o size=64g tmpfs /dev/shm || echo "BENCH_SETUP_FAILED /dev/shm" >> $J/setup.out
# MPI in the chroot: a hello on every vCPU (MPICH's libfabric may need its TCP provider without a network fabric)
printf '#include <mpi.h>\n#include <stdio.h>\nint main(int c, char **v) { int r; MPI_Init(&c, &v); MPI_Comm_rank(MPI_COMM_WORLD, &r);\nif (!r) printf("MPI OK\\n"); MPI_Finalize(); return 0; }\n' > $W/hello.c
mpicc -o $W/hello $W/hello.c > $J/mpi.log 2>&1
if ! timeout 120 mpiexec -n $NP $W/hello >> $J/mpi.log 2>&1; then
  export FI_PROVIDER=tcp
  timeout 120 mpiexec -n $NP $W/hello >> $J/mpi.log 2>&1 || echo "BENCH_SETUP_FAILED mpi (mpi.log)" >> $J/setup.out
fi
B "MPI hello: $(grep -q 'MPI OK' $J/mpi.log && echo ok || echo FAILED)${FI_PROVIDER:+ with FI_PROVIDER=$FI_PROVIDER}"
mkdir -p /opt/netcdf  # WRF's configure wants one NETCDF directory: Debian's files
ln -sfn /usr/bin /opt/netcdf/bin; ln -sfn /usr/include /opt/netcdf/include; ln -sfn /usr/lib/x86_64-linux-gnu /opt/netcdf/lib

wrf_src() {  # the release tarball (with its submodules: Noah-MP, MMM physics), else a git clone of the tag
  [ -s $W/wrf-src.tar.gz ] && return 0
  local u
  for u in https://github.com/wrf-model/WRF/releases/download/v4.8.0/v4.8.0.tar.gz \
           https://github.com/wrf-model/WRF/releases/download/v4.8.0/WRFV4.8.0.tar.gz; do
    if curl -sfL --retry 3 -o $W/src.part "$u" && tar tzf $W/src.part | grep -q 'phys/noahmp/src/'; then
      mv $W/src.part $W/wrf-src.tar.gz; echo "WRF source $u md5 $(md5sum < $W/wrf-src.tar.gz | cut -c1-32)" >> $J/wrf_config.txt
      return 0
    fi
  done
  rm -rf $W/src.part $W/WRF
  git clone -q --depth 1 --branch v4.8.0 --recurse-submodules --shallow-submodules https://github.com/wrf-model/WRF $W/WRF &&
    tar czf $W/wrf-src.tar.gz -C $W WRF && rm -rf $W/WRF && echo "WRF source git clone of tag v4.8.0" >> $J/wrf_config.txt
}
build_wrf() {  # NAME CONFIGURE_OPTION [-r8]
  local d=$W/WRF-$1
  [ -x $d/main/wrf.exe ] && [ -x $d/main/real.exe ] && return 0
  rm -rf $d && mkdir -p $d && tar xzf $W/wrf-src.tar.gz -C $d --strip-components=1 || return 1
  cd $d || return 1
  export NETCDF=/opt/netcdf WRFIO_NCD_LARGE_FILE_SUPPORT=1
  printf '%s\n1\n' "$2" | ./configure $3 > configure.log 2>&1
  [ -f configure.wrf ] || return 1
  ./compile -j 12 em_real > compile.log 2>&1
  [ -x main/wrf.exe ] && [ -x main/real.exe ]
}
wps() {  # met_em files of the large domain in $W/wps
  [ -f $W/wps/met_ok ] && return 0
  rm -rf $W/wps && mkdir -p $W/wps && cd $W/wps && tar xzf $J/wps.tar.gz || return 1
  if [ ! -d $W/geog/WPS_GEOG_LOW_RES ]; then
    mkdir -p $W/geog && curl -sfL --retry 3 https://www2.mmm.ucar.edu/wrf/src/wps_files/geog_low_res_mandatory.tar.gz |
      tar xzf - -C $W/geog || return 1
  fi
  local h L=ABCDEFGHIJ n=0
  for h in 000 006; do
    curl -sf --retry 5 -o gfs.f$h "https://noaa-gfs-bdp-pds.s3.amazonaws.com/gfs.20261012/00/atmos/gfs.t00z.pgrb2.0p25.f$h" || return 1
    ln -sf gfs.f$h GRIBFILE.AA${L:$n:1}; n=$((n + 1))
  done
  ./geogrid.exe > geogrid.log 2>&1; grep -q "Successful completion of geogrid" geogrid.log || return 1
  ./ungrib.exe > ungrib.log 2>&1; grep -q "Successful completion of ungrib" ungrib.log || return 1
  ./metgrid.exe > metgrid.log 2>&1; grep -q "Successful completion of metgrid" metgrid.log || return 1
  rm -f gfs.f* FILE:*
  touch met_ok
}
build_cpu() {  # NAME PORTABLE_LIBMVEC
  cmake -S $J/src -B $W/$1 -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_DOUBLE=ON \
    -DCPPWRF_PORTABLE_LIBMVEC=$2 > $J/build-$1.log 2>&1 && cmake --build $W/$1 -j 16 --target cppwrf >> $J/build-$1.log 2>&1
}

B "WRF source, WRF builds, WPS and CPPWRF CPU builds"
( wrf_src || exit 1
  build_wrf r4 34 "" > /dev/null 2>&1 & build_wrf r8 34 -r8 > /dev/null 2>&1 &
  build_wrf r4sm 35 "" > /dev/null 2>&1 & build_wrf r8sm 35 -r8 > /dev/null 2>&1 &
  wait ) > $J/wrf_build.log 2>&1 &
WRFB=$!
( wps ) > $J/wps.log 2>&1 &
WPSP=$!
t0=$SECONDS
build_cpu cpu-dp ON || echo "BENCH_BUILD_FAILED cpu-dp (build-cpu-dp.log)" >> $J/setup.out
build_cpu cpu-dp-hwlog OFF || echo "BENCH_BUILD_FAILED cpu-dp-hwlog (build-cpu-dp-hwlog.log)" >> $J/setup.out
B "CPPWRF CPU builds: $((SECONDS - t0)) s"
wait $WPSP; B "WPS: $(ls $W/wps/met_em.d01.* 2> /dev/null | wc -l) met_em files ($( [ -f $W/wps/met_ok ] && echo ok || echo FAILED))"
for f in geogrid ungrib metgrid; do tail -3 $W/wps/$f.log > $J/fetch/logs/wps_$f.log 2> /dev/null; done
wait $WRFB
for n in r4 r8 r4sm r8sm; do
  d=$W/WRF-$n
  if [ -x $d/main/wrf.exe ]; then B "WRF-$n built"; else
    B "WRF-$n FAILED"; echo "BENCH_BUILD_FAILED WRF-$n" >> $J/setup.out
    tail -30 $d/configure.log > $J/fetch/logs/WRF-$n.configure.log 2> /dev/null
    grep -B2 -A8 -m5 -E 'Error|error:' $d/compile.log > $J/fetch/logs/WRF-$n.compile_errors.log 2> /dev/null
  fi
  { echo "== WRF-$n"; grep -E '^(DMPARALLEL|OMPCPP|OMP |SFC|SCC|DM_FC|DM_CC|FCOPTIM|FCBASEOPTS_NO_G|PROMOTION|CFLAGS_LOCAL|FCCOMPAT|ARCH_LOCAL)\s' $d/configure.wrf 2> /dev/null; } >> $J/wrf_config.txt
done
B "builds done ($((SECONDS - t0)) s)"
if [ -f $W/wps/met_ok ] && [ -x $W/WRF-r8/main/real.exe ]; then
  [ -s $W/run_L/wrfbdy_d01 ] || python3 $J/src/tools/bench/bench_run.py make-large $J/data/wrf_runs/run_P4_r8 $W/wps
fi
cd $J/src || exit 1
for c in ${CASES//,/ }; do
  case $c in
    R1) python3 tools/bench/bench_run.py case R1 $J/data/wrf_runs/run_R1_r8 --serial --layouts "$LAYOUTS" ;;
    P4) python3 tools/bench/bench_run.py case P4 $J/data/wrf_runs/run_P4_r8 --ref --serial --layouts "$LAYOUTS" ;;
    P4R15) python3 tools/bench/bench_run.py case P4R15 $J/data/wrf_runs/run_P4_r8 --set physics.radt=15 --layouts "$LAYOUTS" ;;
    N3) python3 tools/bench/bench_run.py case N3 $J/data/wrf_runs/run_N3h_r8 --ref --layouts "$LAYOUTS" ;;
    L) if [ -s $W/run_L/wrfbdy_d01 ]; then
         python3 tools/bench/bench_run.py case L $W/run_L --set physics.radt=15 --serial --layouts "$LAYOUTS"
       else B "L: no inputs (WPS or real.exe failed)"; echo "BENCH_L_FAILED" >> $J/setup.out; fi ;;
    *) B "unknown case $c" ;;
  esac >> $J/run.log 2>&1
done
touch $Q/stop
for i in $(seq 12); do [ -e $Q/gpu_exit ] && break; sleep 5; done
cp $Q/gpu.txt $Q/gpu_host.log $Q/build-cuda*.log $Q/apt.log $J/fetch/logs/ 2> /dev/null
cp $Q/gpu.txt $J/gpu.txt 2> /dev/null
python3 tools/bench/bench_run.py report >> $J/run.log 2>&1
cp $J/results.txt $J/report.md $J/wrf_config.txt $J/compilers.txt $J/lscpu.txt $J/bench.log $J/fetch/logs/ 2> /dev/null
tar czf $J/fetch/bench.tar.gz -C $J/fetch logs
B "end"
exit 0
