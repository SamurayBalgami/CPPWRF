#!/usr/bin/env bash
# The VM side of tools/gcp/gcp_run.sh (startup script): installs the build tools, builds the CPU FP64 binary
# with the instruction set of the reference workstation (-march=skylake, as -march=native there), runs
# wrf_suite.py, copies the results to the bucket and powers the VM off.
# Metadata: bucket, src (object of the git archive), data (objects of tools/colab/stage_data.py's tars, comma-separated),
#           result (prefix for the results), jobs, threads, cases (empty: all CPU cases).
md() { curl -s -H "Metadata-Flavor: Google" "http://metadata.google.internal/computeMetadata/v1/instance/attributes/$1"; }
B=$(md bucket) SRC=$(md src) DATA=$(md data) RES=$(md result) JOBS=$(md jobs) THR=$(md threads) CASES=$(md cases)
J=/opt/job; mkdir -p $J; exec > $J/job.log 2>&1
put() { gcloud storage cp "$@" "gs://$B/$RES/" || true; }
trap 'put $J/job.log $J/setup.out; date > $J/done; put $J/done; poweroff' EXIT
set -x
export DEBIAN_FRONTEND=noninteractive
apt-get -qq update && apt-get -qq install -y build-essential cmake libnetcdf-dev python3-netcdf4 python3-numpy > /dev/null
gcloud storage cp "gs://$B/$SRC" $J/src.tar.gz
mkdir -p $J/src $J/data && tar xzf $J/src.tar.gz -C $J/src
# a data tar that did not arrive or a failed build: setup.out (gcp_run.sh then exits with 1)
touch $J/setup.out
for d in ${DATA//,/ }; do  # in the list's order: a later tar replaces an earlier one's files (whole references)
  ( set -o pipefail; gcloud storage cat "gs://$B/$d" | tar xf - -C $J/data ) || echo "DATA_FAILED $d" >> $J/setup.out
done
cmake -S $J/src -B $J/build -DCMAKE_BUILD_TYPE=Release -DCPPWRF_DOUBLE=ON -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_NATIVE=OFF \
  -DCMAKE_CXX_FLAGS="-march=skylake -mtune=skylake" > $J/cmake.log 2>&1 && cmake --build $J/build -j"$(nproc)" --target cppwrf >> $J/cmake.log 2>&1 \
  || echo "BUILD_FAILED (cmake.log)" >> $J/setup.out
put $J/cmake.log $J/setup.out
cd $J/src
# the translator's unit tests (tests/tools/test_wrf2cppwrf.py: namelists of the data's WRF runs, 2-step model runs)
CPPWRF_WRF_RUNS=$J/data/wrf_runs CPPWRF_WRF_REF=$J/data/wrf_ref CPPWRF_BIN=$J/build/cppwrf OMP_NUM_THREADS=1 \
  python3 tests/tools/test_wrf2cppwrf.py -v > $J/unit.log 2>&1; echo "exit $?" >> $J/unit.log
put $J/unit.log
# the comparators' own test: each deliberate fault (one ulp, NaN, -0, a missing variable or frame, ...) must fail
CPPWRF_SCRATCH=$J python3 tests/verify/selftest.py > $J/selftest.log 2>&1; echo "exit $?" >> $J/selftest.log
put $J/selftest.log
only=(); [ -n "$CASES" ] && only=(--only "$CASES" --exact)
(while sleep 60; do put $J/suite.log; done) &  # progress: the lines of the cases that ended
CPPWRF_WRF_RUNS=$J/data/wrf_runs CPPWRF_WRF_REF=$J/data/wrf_ref OMP_NUM_THREADS=$THR \
  python3 tests/verify/wrf_suite.py --cpu-dp $J/build/cppwrf --jobs "$JOBS" --out $J/suite "${only[@]}" > $J/suite.log 2>&1
echo "exit $?" >> $J/suite.log
put $J/suite.log
# the full suite also compares the output files with WRF's (wrfout / wrfrst / resumed runs, ideal cases), one
# run_wrf_files_test.sh per case in parallel; each prints PASS / FAIL / SKIP and its exit status
if [ -z "$CASES" ]; then
  sed -n '/^CASES="/,/^"/p' tests/verify/run_wrf_files_test.sh | grep '|' | cut -d'|' -f1 |
    xargs -P "$JOBS" -I{} sh -c "CPPWRF_WRF_RUNS=$J/data/wrf_runs CPPWRF_SCRATCH=$J OMP_NUM_THREADS=$THR \
      bash tests/verify/run_wrf_files_test.sh $J/build {} > $J/files_{}.log 2>&1; echo \"exit \$? {}\" >> $J/files_{}.log"
  cat $J/files_*.log > $J/files.log
  CPPWRF_WRF_REF=$J/data/wrf_ref CPPWRF_SCRATCH=$J OMP_NUM_THREADS=$THR python3 tests/verify/ideal_files_test.py $J/build \
    > $J/idealfiles.log 2>&1; echo "exit $?" >> $J/idealfiles.log
  put $J/files.log $J/idealfiles.log
fi
