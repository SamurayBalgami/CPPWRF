#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The job of tools/colab/colab_wrfref.sh in the Debian 13 chroot (tools/colab/remote_chroot_job.sh): the work of
# tools/gcp/remote_wrfref.sh without the bucket. Makes every reference of refs.txt (tools/gcp/wrf_ref_run.py: run_NAME
# and run_NAME/steps, in parallel) and its restart references (tests/verify/make_restart_refs.sh,
# make_resumed_refs.sh), lists the MD5 of every file WRF wrote (refs_md5.txt: "MD5 DIR/FILE"), packs them as
# fetch/refs.tar (links kept; FETCH=1), packs the suite's tar of STAGE (OBJ:CASE,...; tools/colab/stage_data.py --whole)
# as drive/OBJ (remote_chroot_job.sh puts it on Drive), then with FILES / CASES builds CPPWRF and runs
# run_wrf_files_test.sh / wrf_suite.py against them (their logs end with "exit N").
#
#   wrfref_body.sh J FILES CASES STAGE FETCH
#
# J: /opt/job with refs.txt, src/ and data/ (data/runs: the base runs, data/wrf: the WRF executables, as on the VM).
J=$1 FILES=$2 CASES=$3 STAGE=$4 FETCH=$5
W=$J/wrf_runs; mkdir -p $W
set -x
touch $J/setup.out
grep -v '^\s*#' $J/refs.txt | grep -v '^\s*$' | while read -r name base exe sets; do
  echo "python3 $J/src/tools/gcp/wrf_ref_run.py $J/data/runs $W $J/data/wrf/$exe $name $base $sets"
done > $J/cmds.txt
xargs -P "$(nproc)" -I{} sh -c '{}' < $J/cmds.txt > $J/refs.out 2>&1
names=$(grep -v '^\s*#' $J/refs.txt | awk 'NF {print $1}')
# shellcheck disable=SC2086
{ CPPWRF_WRF_RUNS=$W bash $J/src/tests/verify/make_restart_refs.sh $names
  CPPWRF_WRF_RUNS=$W bash $J/src/tests/verify/make_resumed_refs.sh $names; } >> $J/refs.out 2>&1
dirs=$(for n in $names; do for p in run rst rsm; do [ -d "$W/${p}_$n" ] && echo "${p}_$n"; done; done)
# shellcheck disable=SC2086
( cd $W && find $dirs -type f \( -name 'wrfout_*' -o -name 'wrfrst_*' -o -name 'auxhist*' -o -name 'wrfinput_*' \
    -o -name 'wrfbdy_*' -o -name 'wrflowinp_*' \) | sort | xargs -r md5sum | awk '{print "MD5", $2, $1}' ) > $J/refs_md5.txt
# shellcheck disable=SC2086
[ "$FETCH" = 0 ] || { mkdir -p $J/fetch && tar cf $J/fetch/refs.tar -C $W $dirs; }
if [ -n "$STAGE" ]; then  # the suite's tar of these references, to Drive
  mkdir -p $J/drive
  CPPWRF_WRF_RUNS=$W python3 $J/src/tools/colab/stage_data.py --whole --only "${STAGE#*:}" "$J/drive/${STAGE%%:*}" > $J/stage.log 2>&1 \
    && echo "STAGED ${STAGE%%:*}" >> $J/stage.log
fi
[ -n "$FILES$CASES" ] || exit 0
cmake -S $J/src -B $J/build -DCMAKE_BUILD_TYPE=Release -DCPPWRF_DOUBLE=ON -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_NATIVE=OFF \
  -DCMAKE_CXX_FLAGS="-march=skylake -mtune=skylake" > $J/cmake.log 2>&1 && cmake --build $J/build -j"$(nproc)" --target cppwrf >> $J/cmake.log 2>&1 \
  || echo "BUILD_FAILED (cmake.log)" >> $J/setup.out
cd $J/src
if [ -n "$FILES" ]; then
  # shellcheck disable=SC2086
  CPPWRF_WRF_RUNS=$W CPPWRF_SCRATCH=$J bash tests/verify/run_wrf_files_test.sh $J/build ${FILES//,/ } > $J/files.log 2>&1
  echo "exit $?" >> $J/files.log
fi
if [ -n "$CASES" ]; then
  CPPWRF_WRF_RUNS=$W OMP_NUM_THREADS=4 python3 tests/verify/wrf_suite.py --cpu-dp $J/build/cppwrf \
    --only "$CASES" --exact --jobs "$(( $(nproc) / 4 > 0 ? $(nproc) / 4 : 1 ))" --out $J/suite > $J/suite.log 2>&1
  echo "exit $?" >> $J/suite.log
fi
exit 0
