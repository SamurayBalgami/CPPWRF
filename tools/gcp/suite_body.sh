#!/usr/bin/env bash
# The full CPU FP64 suite on a machine that holds the source and the unpacked data tars: tools/gcp/remote_suite.sh
# (GCP VM) and tools/colab/remote_cpu_suite.sh (Colab runtime, in a Debian 13 chroot) run it. Builds the CPU FP64
# binary with the instruction set of the reference workstation (-march=skylake, as -march=native there), runs the
# translator tests, the comparators' self-test and wrf_suite.py; the full suite (no cases) also the differential and
# unit tests (tools/gcp/run_diff_tests.sh, next to the suite) and the output files tests.
#
#   suite_body.sh J JOBS THREADS CASES SLOW [FILES]
#
# J holds src/ (git archive) and data/ (the tars unpacked in the data list's order); the logs go to J: setup.out
# (a failed build), cmake.log, unit.log, selftest.log, suite.log, files.log, idealfiles.log, diff.log, ctest.log.
# FILES=1: the output files tests also with CASES (tools/colab/colab_suite.sh --files).
# A shell function `put` (exported by the caller, optional) gets each finished log (the GCP VM copies it to the bucket).
J=$1 JOBS=$2 THR=$3 CASES=$4 SLOW=$5 FILES=${6:-0}
type put > /dev/null 2>&1 || put() { :; }
set -x
touch $J/setup.out
cmake -S $J/src -B $J/build -DCMAKE_BUILD_TYPE=Release -DCPPWRF_DOUBLE=ON -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_NATIVE=OFF \
  -DCMAKE_CXX_FLAGS="-march=skylake -mtune=skylake" > $J/cmake.log 2>&1 && cmake --build $J/build -j"$(nproc)" --target cppwrf >> $J/cmake.log 2>&1 \
  || echo "BUILD_FAILED (cmake.log)" >> $J/setup.out
put $J/cmake.log $J/setup.out
cd $J/src
# the differential tests against WRF's objects (diff_data.tar) and the ctest unit tests, in their own build, while
# the suite runs
if [ -z "$CASES" ]; then
  (OMP_NUM_THREADS=4 bash tools/gcp/run_diff_tests.sh $J/src $J/build_diff $J/data $J > $J/diff_setup.out 2>&1
   put $J/diff_cmake.log $J/diff.log $J/ctest.log $J/diff_setup.out) &
  DIFF_PID=$!
fi
# the translator's unit tests (tests/tools/test_wrf2cppwrf.py: namelists of the data's WRF runs, 2-step model runs)
CPPWRF_WRF_RUNS=$J/data/wrf_runs CPPWRF_WRF_REF=$J/data/wrf_ref CPPWRF_BIN=$J/build/cppwrf OMP_NUM_THREADS=1 \
  python3 tests/tools/test_wrf2cppwrf.py -v > $J/unit.log 2>&1; echo "exit $?" >> $J/unit.log
put $J/unit.log
# the comparators' own test: each deliberate fault (one ulp, NaN, -0, a missing variable or frame, ...) must fail
CPPWRF_SCRATCH=$J python3 tests/verify/selftest.py > $J/selftest.log 2>&1; echo "exit $?" >> $J/selftest.log
put $J/selftest.log
only=(); [ -n "$CASES" ] && only=(--only "$CASES" --exact)
[ "$SLOW" = 1 ] && only+=(--slow)  # the slow cases too (gcp_run.sh --slow)
(while sleep 60; do put $J/suite.log; done) &  # progress: the lines of the cases that ended
PROG_PID=$!
CPPWRF_WRF_RUNS=$J/data/wrf_runs CPPWRF_WRF_REF=$J/data/wrf_ref OMP_NUM_THREADS=$THR \
  python3 tests/verify/wrf_suite.py --cpu-dp $J/build/cppwrf --jobs "$JOBS" --out $J/suite "${only[@]}" > $J/suite.log 2>&1
echo "exit $?" >> $J/suite.log
kill $PROG_PID 2> /dev/null
put $J/suite.log
# the full suite also compares the output files with WRF's (wrfout / wrfrst / resumed runs, ideal cases), one
# run_wrf_files_test.sh per case in parallel; each prints PASS / FAIL / SKIP and its exit status
if [ -z "$CASES" ] || [ "$FILES" = 1 ]; then
  sed -n '/^CASES="/,/^"/p' tests/verify/run_wrf_files_test.sh | grep '|' | cut -d'|' -f1 |
    xargs -P "$JOBS" -I{} sh -c "CPPWRF_WRF_RUNS=$J/data/wrf_runs CPPWRF_SCRATCH=$J OMP_NUM_THREADS=$THR \
      bash tests/verify/run_wrf_files_test.sh $J/build {} > $J/files_{}.log 2>&1; echo \"exit \$? {}\" >> $J/files_{}.log"
  cat $J/files_*.log > $J/files.log
  CPPWRF_WRF_REF=$J/data/wrf_ref CPPWRF_SCRATCH=$J OMP_NUM_THREADS=$THR python3 tests/verify/ideal_files_test.py $J/build \
    > $J/idealfiles.log 2>&1; echo "exit $?" >> $J/idealfiles.log
  put $J/files.log $J/idealfiles.log
fi
if [ -z "$CASES" ]; then
  wait $DIFF_PID
  cat $J/diff_setup.out >> $J/setup.out; put $J/setup.out
fi
