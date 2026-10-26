#!/usr/bin/env bash
# Startup script of the reference VM (tools/gcp/gcp_wrfref.sh): unpacks the WRF runs and wrf.exe, makes every
# reference of the list (tools/gcp/wrf_ref_run.py: run_NAME and run_NAME/steps, in parallel) and its restart
# references (tests/verify/make_restart_refs.sh, make_resumed_refs.sh), uploads them as refs.tar (links kept; not
# with fetch=0);
# with `files` / `cases` it then builds CPPWRF and runs run_wrf_files_test.sh / wrf_suite.py against them.
md() { curl -s -H "Metadata-Flavor: Google" "http://metadata.google.internal/computeMetadata/v1/instance/attributes/$1"; }
B=$(md bucket) SRC=$(md src) DATA=$(md data) RES=$(md result) REFS=$(md refs) FILES=$(md files) CASES=$(md cases) FETCH=$(md fetch)
J=/opt/job W=/opt/job/wrf_runs; mkdir -p $J $W; exec > $J/job.log 2>&1
put() { gcloud storage cp "$@" "gs://$B/$RES/" || true; }
trap 'put $J/job.log $J/refs.out $J/files.log $J/suite.log; date > $J/done; put $J/done; poweroff' EXIT
set -x
export DEBIAN_FRONTEND=noninteractive
apt-get -qq update && apt-get -qq install -y build-essential cmake libnetcdf-dev libnetcdff-dev python3-netcdf4 python3-numpy > /dev/null
gcloud storage cp "gs://$B/$SRC" $J/src.tar.gz && gcloud storage cp "gs://$B/$REFS" $J/refs.txt
mkdir -p $J/src $J/data && tar xzf $J/src.tar.gz -C $J/src
for d in ${DATA//,/ }; do gcloud storage cp "gs://$B/$d" $J/data.tar && tar xf $J/data.tar -C $J/data && rm $J/data.tar; done
(while sleep 60; do put $J/refs.out; done) &
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
[ "$FETCH" = 0 ] || { tar cf $J/refs.tar -C $W $dirs && put $J/refs.tar; }
[ -n "$FILES$CASES" ] || exit 0
cmake -S $J/src -B $J/build -DCMAKE_BUILD_TYPE=Release -DCPPWRF_DOUBLE=ON -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_NATIVE=OFF \
  -DCMAKE_CXX_FLAGS="-march=skylake -mtune=skylake" > $J/cmake.log 2>&1 && cmake --build $J/build -j"$(nproc)" --target cppwrf >> $J/cmake.log 2>&1
put $J/cmake.log
cd $J/src
# shellcheck disable=SC2086
[ -n "$FILES" ] && CPPWRF_WRF_RUNS=$W CPPWRF_SCRATCH=$J bash tests/verify/run_wrf_files_test.sh $J/build ${FILES//,/ } > $J/files.log 2>&1
[ -n "$CASES" ] && CPPWRF_WRF_RUNS=$W OMP_NUM_THREADS=4 python3 tests/verify/wrf_suite.py --cpu-dp $J/build/cppwrf \
  --only "$CASES" --exact --jobs "$(( $(nproc) / 4 ))" --out $J/suite > $J/suite.log 2>&1
exit 0
