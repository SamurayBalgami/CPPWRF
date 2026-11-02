#!/usr/bin/env bash
# The VM side of tools/gcp/gcp_run.sh (startup script): installs the build tools, fetches the source and the data,
# runs tools/gcp/suite_body.sh (the CPU FP64 build, translator tests, comparator self-test, wrf_suite.py; the full
# suite (no cases) also the differential, unit and output files tests), copies the results to the bucket and powers
# the VM off.
# Metadata: bucket, src (object of the git archive), data (objects of tools/colab/stage_data.py's tars, comma-separated),
#           result (prefix for the results), jobs, threads, cases (empty: all CPU cases).
md() { curl -s -H "Metadata-Flavor: Google" "http://metadata.google.internal/computeMetadata/v1/instance/attributes/$1"; }
B=$(md bucket) SRC=$(md src) DATA=$(md data) RES=$(md result) JOBS=$(md jobs) THR=$(md threads) CASES=$(md cases) SLOW=$(md slow)
J=/opt/job; mkdir -p $J; exec > $J/job.log 2>&1
put() { gcloud storage cp "$@" "gs://$B/$RES/" || true; }
trap 'put $J/job.log $J/setup.out; date > $J/done; put $J/done; poweroff' EXIT
set -x
export DEBIAN_FRONTEND=noninteractive
apt-get -qq update && apt-get -qq install -y build-essential cmake gfortran libnetcdf-dev libnetcdff-dev python3-netcdf4 python3-numpy > /dev/null
gcloud storage cp "gs://$B/$SRC" $J/src.tar.gz
mkdir -p $J/src $J/data && tar xzf $J/src.tar.gz -C $J/src
# a data tar that did not arrive or a failed build: setup.out (gcp_run.sh then exits with 1)
touch $J/setup.out
for d in ${DATA//,/ }; do  # in the list's order: a later tar replaces an earlier one's files (whole references)
  ( set -o pipefail; gcloud storage cat "gs://$B/$d" | tar xf - -C $J/data ) || echo "DATA_FAILED $d" >> $J/setup.out
done
# the build and the tests (tools/gcp/suite_body.sh, shared with the Colab runtime's tools/colab/remote_cpu_suite.sh)
export -f put; export B RES  # put() reads them
bash $J/src/tools/gcp/suite_body.sh $J "$JOBS" "$THR" "$CASES" "$SLOW"
