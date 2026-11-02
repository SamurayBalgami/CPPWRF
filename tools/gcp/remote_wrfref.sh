#!/usr/bin/env bash
# Startup script of the reference VM (tools/gcp/gcp_wrfref.sh): unpacks the WRF runs and wrf.exe, makes every
# reference of the list (tools/gcp/wrf_ref_run.py: run_NAME and run_NAME/steps, in parallel) and its restart
# references (tests/verify/make_restart_refs.sh, make_resumed_refs.sh), uploads them as refs.tar (links kept; not
# with fetch=0);
# with `files` / `cases` it then builds CPPWRF and runs run_wrf_files_test.sh / wrf_suite.py against them (their logs
# end with "exit N"). Inputs that did not arrive and a failed build are lines of setup.out (gcp_wrfref.sh then fails).
# With `stage` (OBJ:CASE,...) it packs the references of these suite cases (tools/colab/stage_data.py --whole) and
# uploads the tar as the bucket object OBJ. With `cache` = 1 a reference whose line (hash) is in
# gs://BUCKET/wrfref_cache/ comes from there, and every reference that succeeds goes there (a spot preemption then does
# not lose the finished ones; run.segment's parts too).
md() { curl -s -H "Metadata-Flavor: Google" "http://metadata.google.internal/computeMetadata/v1/instance/attributes/$1"; }
B=$(md bucket) SRC=$(md src) DATA=$(md data) RES=$(md result) REFS=$(md refs) FILES=$(md files) CASES=$(md cases) FETCH=$(md fetch) STAGE=$(md stage) CACHE=$(md cache)
J=/opt/job W=/opt/job/wrf_runs; mkdir -p $J $W; exec > $J/job.log 2>&1
put() { gcloud storage cp "$@" "gs://$B/$RES/" || true; }
trap 'put $J/job.log $J/setup.out $J/refs.out $J/stage.log $J/files.log $J/suite.log; date > $J/done; put $J/done; poweroff' EXIT
set -x
export DEBIAN_FRONTEND=noninteractive
apt-get -qq update && apt-get -qq install -y build-essential cmake libnetcdf-dev libnetcdff-dev python3-netcdf4 python3-numpy > /dev/null
touch $J/setup.out
{ gcloud storage cp "gs://$B/$SRC" $J/src.tar.gz && gcloud storage cp "gs://$B/$REFS" $J/refs.txt; } || echo "SRC_FAILED" >> $J/setup.out
mkdir -p $J/src $J/data && { tar xzf $J/src.tar.gz -C $J/src || echo "SRC_FAILED" >> $J/setup.out; }
for d in ${DATA//,/ }; do
  { gcloud storage cp "gs://$B/$d" $J/data.tar && tar xf $J/data.tar -C $J/data; } || echo "DATA_FAILED $d" >> $J/setup.out
  rm -f $J/data.tar
done
(while sleep 60; do put $J/refs.out; done) &
cat > $J/cached.sh <<END
#!/usr/bin/env bash
# cached.sh KEY NAME COMMAND...: run_NAME from the bucket's cache, else COMMAND, cached when all its runs succeed
key=\$1 name=\$2; shift 2; obj=gs://$B/wrfref_cache/\$name-\$key.tar
if [ "$CACHE" = 1 ] && gcloud storage cp "\$obj" $J/cache_\$name.tar 2> /dev/null; then
  tar xf $J/cache_\$name.tar -C $W && rm $J/cache_\$name.tar && echo "\$name SUCCESS (from \$obj)"; exit 0
fi
[ "$CACHE" = 1 ] && export WRFREF_PARTS=gs://$B/wrfref_cache/\$name-\$key  # run.segment's parts
out=\$("\$@" 2>&1); echo "\$out"
if [ "$CACHE" = 1 ] && ! echo "\$out" | grep -v ' SUCCESS\$' | grep -q .; then
  tar cf $J/cache_\$name.tar -C $W run_\$name && gcloud storage cp $J/cache_\$name.tar "\$obj" -q > /dev/null 2>&1 && rm $J/cache_\$name.tar
fi
END
grep -v '^\s*#' $J/refs.txt | grep -v '^\s*$' | while read -r name base exe sets; do
  key=$(echo "$name $base $exe $sets" | md5sum | cut -c1-12)
  echo "bash $J/cached.sh $key $name python3 $J/src/tools/gcp/wrf_ref_run.py $J/data/runs $W $J/data/wrf/$exe $name $base $sets"
done > $J/cmds.txt
xargs -P "$(nproc)" -I{} sh -c '{}' < $J/cmds.txt > $J/refs.out 2>&1
names=$(grep -v '^\s*#' $J/refs.txt | awk 'NF {print $1}')
# shellcheck disable=SC2086
{ CPPWRF_WRF_RUNS=$W bash $J/src/tests/verify/make_restart_refs.sh $names
  CPPWRF_WRF_RUNS=$W bash $J/src/tests/verify/make_resumed_refs.sh $names; } >> $J/refs.out 2>&1
dirs=$(for n in $names; do for p in run rst rsm; do [ -d "$W/${p}_$n" ] && echo "${p}_$n"; done; done)
# shellcheck disable=SC2086
[ "$FETCH" = 0 ] || { tar cf $J/refs.tar -C $W $dirs && put $J/refs.tar; }
if [ -n "$STAGE" ]; then  # the suite's tar of these references, straight to the bucket
  { CPPWRF_WRF_RUNS=$W python3 $J/src/tools/colab/stage_data.py --whole --only "${STAGE#*:}" $J/stage.tar \
      && gcloud storage cp $J/stage.tar "gs://$B/${STAGE%%:*}" && echo "STAGED gs://$B/${STAGE%%:*}"; } > $J/stage.log 2>&1
  put $J/stage.log
fi
[ -n "$FILES$CASES" ] || exit 0
cmake -S $J/src -B $J/build -DCMAKE_BUILD_TYPE=Release -DCPPWRF_DOUBLE=ON -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_NATIVE=OFF \
  -DCMAKE_CXX_FLAGS="-march=skylake -mtune=skylake" > $J/cmake.log 2>&1 && cmake --build $J/build -j"$(nproc)" --target cppwrf >> $J/cmake.log 2>&1 \
  || echo "BUILD_FAILED (cmake.log)" >> $J/setup.out
put $J/cmake.log $J/setup.out
cd $J/src
if [ -n "$FILES" ]; then
  # shellcheck disable=SC2086
  CPPWRF_WRF_RUNS=$W CPPWRF_SCRATCH=$J bash tests/verify/run_wrf_files_test.sh $J/build ${FILES//,/ } > $J/files.log 2>&1
  echo "exit $?" >> $J/files.log
fi
if [ -n "$CASES" ]; then
  CPPWRF_WRF_RUNS=$W OMP_NUM_THREADS=4 python3 tests/verify/wrf_suite.py --cpu-dp $J/build/cppwrf \
    --only "$CASES" --exact --jobs "$(( $(nproc) / 4 ))" --out $J/suite > $J/suite.log 2>&1
  echo "exit $?" >> $J/suite.log
fi
exit 0
