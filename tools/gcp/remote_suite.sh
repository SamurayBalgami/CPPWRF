#!/usr/bin/env bash
# The VM side of tools/gcp/gcp_run.sh (startup script): installs the build tools, builds the CPU FP64 binary
# with the instruction set of the reference workstation (-march=skylake, as -march=native there), runs
# wrf_suite.py, copies the results to the bucket and powers the VM off.
# Metadata: bucket, src (object of the git archive), data (object of tools/colab/stage_data.py's tar),
#           result (prefix for the results), jobs, threads, cases (empty: all CPU cases).
md() { curl -s -H "Metadata-Flavor: Google" "http://metadata.google.internal/computeMetadata/v1/instance/attributes/$1"; }
B=$(md bucket) SRC=$(md src) DATA=$(md data) RES=$(md result) JOBS=$(md jobs) THR=$(md threads) CASES=$(md cases)
J=/opt/job; mkdir -p $J; exec > $J/job.log 2>&1
put() { gcloud storage cp "$@" "gs://$B/$RES/" || true; }
trap 'put $J/job.log; date > $J/done; put $J/done; poweroff' EXIT
set -x
export DEBIAN_FRONTEND=noninteractive
apt-get -qq update && apt-get -qq install -y build-essential cmake libnetcdf-dev python3-netcdf4 python3-numpy > /dev/null
gcloud storage cp "gs://$B/$SRC" $J/src.tar.gz && gcloud storage cp "gs://$B/$DATA" $J/data.tar
mkdir -p $J/src $J/data && tar xzf $J/src.tar.gz -C $J/src && tar xf $J/data.tar -C $J/data && rm $J/data.tar
cmake -S $J/src -B $J/build -DCMAKE_BUILD_TYPE=Release -DCPPWRF_DOUBLE=ON -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_NATIVE=OFF \
  -DCMAKE_CXX_FLAGS="-march=skylake -mtune=skylake" > $J/cmake.log 2>&1 && cmake --build $J/build -j"$(nproc)" --target cppwrf >> $J/cmake.log 2>&1
put $J/cmake.log
cd $J/src
only=(); [ -n "$CASES" ] && only=(--only "$CASES" --exact)
CPPWRF_WRF_RUNS=$J/data/wrf_runs CPPWRF_WRF_REF=$J/data/wrf_ref OMP_NUM_THREADS=$THR \
  python3 tests/verify/wrf_suite.py --cpu-dp $J/build/cppwrf --jobs "$JOBS" --out $J/suite "${only[@]}" > $J/suite.log 2>&1
echo "exit $?" >> $J/suite.log
put $J/suite.log
