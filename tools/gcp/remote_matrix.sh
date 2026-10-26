#!/usr/bin/env bash
# Startup script of the variant-matrix VM (tools/gcp/gcp_matrix.sh): builds CPPWRF, unpacks the WRF runs and
# wrf.exe, runs tools/gcp/matrix_variant.py for every line of the variant list in parallel, uploads the results.
md() { curl -s -H "Metadata-Flavor: Google" "http://metadata.google.internal/computeMetadata/v1/instance/attributes/$1"; }
B=$(md bucket) SRC=$(md src) DATA=$(md data) RES=$(md result) VARS=$(md variants)
J=/opt/job; mkdir -p $J; exec > $J/job.log 2>&1
put() { gcloud storage cp "$@" "gs://$B/$RES/" || true; }
trap 'put $J/job.log $J/matrix.out; date > $J/done; put $J/done; poweroff' EXIT
set -x
export DEBIAN_FRONTEND=noninteractive
apt-get -qq update && apt-get -qq install -y build-essential cmake libnetcdf-dev libnetcdff-dev python3-netcdf4 python3-numpy > /dev/null
gcloud storage cp "gs://$B/$SRC" $J/src.tar.gz && gcloud storage cp "gs://$B/$VARS" $J/variants.txt
mkdir -p $J/src $J/data $J/work && tar xzf $J/src.tar.gz -C $J/src
for d in ${DATA//,/ }; do gcloud storage cp "gs://$B/$d" $J/data.tar && tar xf $J/data.tar -C $J/data && rm $J/data.tar; done
cmake -S $J/src -B $J/build -DCMAKE_BUILD_TYPE=Release -DCPPWRF_DOUBLE=ON -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_NATIVE=OFF \
  -DCMAKE_CXX_FLAGS="-march=skylake -mtune=skylake" > $J/cmake.log 2>&1 && cmake --build $J/build -j"$(nproc)" --target cppwrf cppwrf_real >> $J/cmake.log 2>&1
put $J/cmake.log
(while sleep 60; do put $J/matrix.out; done) &
# a variant whose executable is real_*.exe compares WRF real.exe with cppwrf_real (tools/gcp/real_variant.py)
grep -v '^\s*#' $J/variants.txt | grep -v '^\s*$' | while read -r name base exe sets; do
  case $exe in
    real_*) echo "python3 $J/src/tools/gcp/real_variant.py $J/data/runs $J/work $J/data/wrf/$exe $J/build/src/init_real/cppwrf_real $name $base $sets" ;;
    *) echo "python3 $J/src/tools/gcp/matrix_variant.py $J/data/runs $J/work $J/data/wrf/$exe $J/build/cppwrf $name $base $sets" ;;
  esac
done > $J/cmds.txt
xargs -P "$(nproc)" -I{} sh -c '{}' < $J/cmds.txt > $J/matrix.out 2>&1
[ -d $J/work/keep ] && tar cf $J/keep.tar -C $J/work/keep . && put $J/keep.tar  # run.keep=1 references
echo "MATRIX_DONE" >> $J/matrix.out
