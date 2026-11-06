#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The differential tests on a VM (tools/gcp/remote_suite.sh for the full suite, remote_difftests.sh alone): a CPU
# FP64 build with -DCPPWRF_DIFF_TESTS=ON, `ctest -L diff` against the WRF objects of tools/gcp/stage_diff_data.sh's
# tar, then the other ctest tests (unit, io_netcdf) except the two that need whole WRF runs (real_vs_wrf,
# wrf_files_vs_wrf: the suite's files tests and real_variant cover them) and the translator's (run on its own).
#
#   run_diff_tests.sh SRC BUILD DATA OUT [FAULT.patch ...]
#
# DATA: the unpacked data tars (DATA/wrf_build/{r8,r8d8}; the runs of DATA/wrf_runs where there, else
# DATA/wrf_runs_diff). OUT gets diff.log and ctest.log, each ending with "exit N" (ctest's status; a skipped test
# is a failure here: "SKIPPED" line). Each FAULT.patch (patch -p1 on SRC; a line "# tests: REGEX" names the ctest
# tests it must break) is applied after these, built, its tests run, and taken back: OUT/faults.log gets
# "FAULT NAME DETECTED" when one of them failed and "FAULT NAME MISSED" when all passed. DIFF_TESTS (a ctest -R
# regex) limits both runs to those tests.
set -uo pipefail
SRC=$(readlink -f "$1") BUILD=$2 DATA=$(readlink -f "$3") OUT=$(readlink -f "$4"); shift 4
mkdir -p "$BUILD" "$OUT"
cmake -S "$SRC" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DCPPWRF_DOUBLE=ON -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_NATIVE=OFF \
  -DCPPWRF_DIFF_TESTS=ON -DCMAKE_CXX_FLAGS="-march=skylake -mtune=skylake" > "$OUT/diff_cmake.log" 2>&1 &&
  cmake --build "$BUILD" -j"$(nproc)" >> "$OUT/diff_cmake.log" 2>&1 ||
  { echo "BUILD_FAILED (diff_cmake.log)"; echo "exit 1" > "$OUT/diff.log"; exit 1; }

# the WRF runs the tests read, file by file: the diff tar's files (those the tests read), then the suite's others
# (a suite run directory need not hold them, or holds reduced references: run_R0_r8, run_P4_r8 wrfout)
RUNS=$OUT/diff_runs; rm -rf "$RUNS"; mkdir -p "$RUNS"
for d in "$DATA"/wrf_runs/*; do
  [ -e "$DATA/wrf_runs_diff/$(basename "$d")" ] || ln -s "$d" "$RUNS/"
done
for d in "$DATA"/wrf_runs_diff/*; do
  [ -e "$d" ] || continue
  for s in "$d" "$DATA/wrf_runs/$(basename "$d")"; do
    [ -d "$s" ] || continue
    ( cd "$s" && find -L . -type d ) | while read -r p; do mkdir -p "$RUNS/$(basename "$d")/$p"; done
    ( cd "$s" && find -L . ! -type d ) | while read -r p; do
      [ -e "$RUNS/$(basename "$d")/$p" ] || ln -s "$s/$p" "$RUNS/$(basename "$d")/$p"
    done
  done
done
export WRF_DIR=$DATA/wrf_build/r8 WRF_R8D8_DIR=$DATA/wrf_build/r8d8 CPPWRF_WRF_RUNS=$RUNS
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4} OMP_WAIT_POLICY=passive

ctest_log() {  # LOG ctest-args...: the run, the skipped tests as SKIPPED lines, "exit N"
  local log=$1; shift
  ( cd "$BUILD" && ctest --output-on-failure --timeout 5400 "$@" ) > "$log" 2>&1; local s=$?
  # ctest: "  12 - test_phys_wsm6 (Skipped)" / "(Not Run)" under "The following tests did not run:"
  grep -E '^\s+[0-9]+ - .* \((Skipped|Not Run|Disabled)\)' "$log" | sed 's/^/SKIPPED /' >> "$log"
  grep -q '^SKIPPED ' "$log" && [ $s = 0 ] && s=1
  echo "exit $s" >> "$log"
  return $s
}
# serial: the column tests of one directory share its work files
sel=(); [ -n "${DIFF_TESTS:-}" ] && sel=(-R "$DIFF_TESTS")
ctest_log "$OUT/diff.log" -L diff "${sel[@]}"; st=$?
ctest_log "$OUT/ctest.log" -LE diff -E '^(real_vs_wrf|wrf_files_vs_wrf|wrf2cppwrf_translator)$' -j"$(nproc)" "${sel[@]}" || st=1

: > "$OUT/faults.log"
for f in "$@"; do
  name=$(basename "$f" .patch) re=$(sed -n 's/^# tests: //p' "$f" | head -1)
  if [ -z "$re" ] || ! patch -p1 -d "$SRC" -s --dry-run < "$f" > /dev/null 2>&1; then
    echo "FAULT $name NOT-APPLIED" >> "$OUT/faults.log"; continue
  fi
  patch -p1 -d "$SRC" -s < "$f"
  if cmake --build "$BUILD" -j"$(nproc)" > "$OUT/fault_$name.log" 2>&1; then
    ( cd "$BUILD" && ctest --output-on-failure --timeout 5400 -R "$re" ) >> "$OUT/fault_$name.log" 2>&1 && r=MISSED || r=DETECTED
    grep -qE '\((Skipped|Not Run)\)|No tests were found' "$OUT/fault_$name.log" && r=MISSED
  else r=BUILD-FAILED; fi
  echo "FAULT $name $r ($re)" >> "$OUT/faults.log"
  patch -p1 -R -d "$SRC" -s < "$f"
done
[ $# = 0 ] || cmake --build "$BUILD" -j"$(nproc)" > /dev/null 2>&1  # the sources as they were
exit $st
