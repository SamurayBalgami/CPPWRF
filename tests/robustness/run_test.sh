#!/usr/bin/env bash
# Deliberate faults the model must refuse, with their messages (ctest robustness; analytic cases, seconds, no data):
#  - a misspelled configuration key, and a key in the wrong section: "unknown configuration key";
#  - a column above the physics' level limit (nz = 129 with Kessler): an error in Release builds too;
#  - a resume from CPPWRF's own restart file with other options: refused (the same options resume, bit-identical);
#  - run.check_nonfinite_every_step: an unstable run stops at the first non-finite step (the step a diagnostic every
#    step finds), naming the first non-finite field and point.
#   run_test.sh CPPWRF_BINARY WORK_DIR      (from the repository root)
set -u
BIN=$1 W=$2
mkdir -p "$W"
fails=0
note() { echo "$1"; [ "${1:0:4}" = FAIL ] && fails=$((fails + 1)); }
# run NAME CONFIG ARGS...: the run's output in $W/NAME, its log in $W/NAME.log; returns the exit status
run() {
  local n=$1 c=$2; shift 2
  rm -rf "$W/$n"
  "$BIN" "$c" --set output.dir="$W/$n" "$@" > "$W/$n.log" 2>&1
}
SQ=(configs/squall2d.ini --set domain.nx=20 --set domain.nz=20 --set run.run_time=12 --set output.diag_interval=1
    --set output.out_interval=1)
expect_ok() {
  local n=$1; shift
  if run "$n" "$@"; then note "ok   $n"; else note "FAIL $n: exit $? ($(tail -1 "$W/$n.log"))"; fi
}
expect_fail() {  # NAME MESSAGE CONFIG ARGS...
  local n=$1 msg=$2; shift 2
  if run "$n" "$@"; then note "FAIL $n: ran (expected '$msg')"
  elif grep -qF -- "$msg" "$W/$n.log"; then note "ok   $n: $(grep -F -- "$msg" "$W/$n.log" | head -1)"
  else note "FAIL $n: no '$msg' in: $(tail -2 "$W/$n.log" | tr '\n' ' ')"; fi
}

expect_ok base "${SQ[@]}"
expect_fail misspelled_key "unknown configuration key 'phys.mp_physic'" "${SQ[@]}" --set phys.mp_physic=1
expect_fail wrong_section "unknown configuration key 'dyn.mp_physics'" "${SQ[@]}" --set dyn.mp_physics=1
expect_fail nz_limit "Kessler: nz = 129 exceeds the column limit 128" "${SQ[@]}" --set domain.nz=129

# CPPWRF's own restart file: written at steps 1 and 2 (restart.interval = dt)
if run rst_a "${SQ[@]}" --set restart.interval=6; then
  expect_ok rst_same "${SQ[@]}" --set run.restart_file="$W/rst_a/restart_000001.nc"
  if [ -f "$W/rst_same/snap_000002.bin" ] && cmp -s "$W/rst_a/snap_000002.bin" "$W/rst_same/snap_000002.bin"; then
    note "ok   rst_same: step 2 bit-identical"
  else note "FAIL rst_same: step 2 differs from the uninterrupted run"; fi
  expect_fail rst_options "written by a run with other options (dyn.khdif = 0 in the file, 10 here)" "${SQ[@]}" \
    --set run.restart_file="$W/rst_a/restart_000001.nc" --set dyn.khdif=10
elif grep -q "require NetCDF" "$W/rst_a.log"; then
  echo "skip restart checks: no NetCDF in this build"
else
  note "FAIL rst_a: exit ($(tail -1 "$W/rst_a.log"))"
fi

# an unstable dry bubble (dt far beyond the acoustic limit): the step of the first non-finite value
WB=(configs/warm_bubble.ini --set domain.nx=40 --set domain.nz=20 --set time.dt=100 --set run.run_time=5000
    --set output.out_interval=1000)
run nan_diag "${WB[@]}" --set output.diag_interval=1
k1=$(grep -o "non-finite values at step [0-9]*" "$W/nan_diag.log" | head -1 | grep -o "[0-9]*$")
run nan_every "${WB[@]}" --set output.diag_interval=1000 --set run.check_nonfinite_every_step=true
k2=$(grep -o "non-finite values at step [0-9]*" "$W/nan_every.log" | head -1 | grep -o "[0-9]*$")
if [ -z "$k1" ]; then note "FAIL nan_diag: the unstable run gave no non-finite value ($(tail -1 "$W/nan_diag.log"))"
elif [ "$k1" = 50 ]; then note "FAIL nan_diag: non-finite only at the last step (no early blow-up to test with)"
elif [ "$k2" != "$k1" ]; then note "FAIL nan_every: stopped at step '${k2}', the diagnostic every step at $k1"
elif grep -q "first: none" "$W/nan_every.log" || ! grep -q "first: .* at (i, j, k)" "$W/nan_every.log"; then
  note "FAIL nan_every: no field and point named: $(grep -m1 non-finite "$W/nan_every.log")"
else note "ok   nan_every: $(grep -m1 non-finite "$W/nan_every.log")"; fi

[ $fails = 0 ] && echo "PASS robustness" || echo "FAIL robustness ($fails)"
exit $((fails > 0))
