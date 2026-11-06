#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The differential tests alone (tools/colab/colab_suite.sh --diff-only) in the Debian 13 chroot of a Colab runtime:
# tools/gcp/run_diff_tests.sh as the full suite runs it (its own CPU FP64 build, ctest -L diff, then the other ctest
# tests), with the deliberate faults /opt/job/fault_*.patch (colab_suite.sh --faults).
#
#   diff_body.sh J TESTS
#
# J holds src/ (git archive) and data/ (diff_data.tar, port_c01_diff.tar unpacked); TESTS: a ctest -R regex limiting
# both ctest runs (empty: all). Logs in J: cpu.txt, diff_cmake.log, diff.log, ctest.log, faults.log, fault_*.log and
# setup.out (a failed build).
J=$1 TESTS=$2
set -x
touch $J/setup.out
lscpu | grep -E '^(Vendor ID|Model name|CPU\(s\)):' > $J/cpu.txt
cd $J/src || exit 1
DIFF_TESTS=$TESTS OMP_NUM_THREADS=4 bash tools/gcp/run_diff_tests.sh $J/src $J/build_diff $J/data $J \
  $(ls $J/fault_*.patch 2> /dev/null) > $J/diff_setup.out 2>&1
cat $J/diff_setup.out >> $J/setup.out
