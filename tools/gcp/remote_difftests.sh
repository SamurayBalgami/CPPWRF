#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Startup script of tools/gcp/gcp_difftests.sh's VM: installs the build tools and gfortran, unpacks the source and
# the data tars, runs tools/gcp/run_diff_tests.sh (with the fault patches of the `faults` tar), copies the logs to
# the bucket and powers the VM off. Inputs that did not arrive are lines of setup.out.
# Metadata: bucket, src, data (comma-separated tars), faults (object of a tar of *.patch files, or empty), result,
#           tests (a ctest -R regex limiting the run, or empty).
md() { curl -s -H "Metadata-Flavor: Google" "http://metadata.google.internal/computeMetadata/v1/instance/attributes/$1"; }
B=$(md bucket) SRC=$(md src) DATA=$(md data) FAULTS=$(md faults) RES=$(md result) TESTS=$(md tests)
J=/opt/job; mkdir -p $J; exec > $J/job.log 2>&1
put() { gcloud storage cp "$@" "gs://$B/$RES/" || true; }
trap 'put $J/job.log $J/setup.out $J/diff_cmake.log $J/diff.log $J/ctest.log $J/faults.log $J/fault_*.log; date > $J/done; put $J/done; poweroff' EXIT
set -x
export DEBIAN_FRONTEND=noninteractive
apt-get -qq update && apt-get -qq install -y build-essential cmake gfortran libnetcdf-dev libnetcdff-dev \
  python3-netcdf4 python3-numpy patch > /dev/null
touch $J/setup.out
mkdir -p $J/src $J/data $J/faults
{ gcloud storage cp "gs://$B/$SRC" $J/src.tar.gz && tar xzf $J/src.tar.gz -C $J/src; } || echo "SRC_FAILED" >> $J/setup.out
for d in ${DATA//,/ }; do
  ( set -o pipefail; gcloud storage cat "gs://$B/$d" | tar xf - -C $J/data ) || echo "DATA_FAILED $d" >> $J/setup.out
done
[ -z "$FAULTS" ] || { gcloud storage cp "gs://$B/$FAULTS" $J/faults.tar && tar xf $J/faults.tar -C $J/faults; } ||
  echo "FAULTS_FAILED" >> $J/setup.out
DIFF_TESTS=$TESTS bash $J/src/tools/gcp/run_diff_tests.sh $J/src $J/build $J/data $J $(ls $J/faults/*.patch 2>/dev/null) >> $J/setup.out
exit 0
