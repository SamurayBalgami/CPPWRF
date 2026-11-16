#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# one run: $1 checkout, $2 build dir, $3 config (relative to the checkout), $4 output dir, rest: --set ...
# CPPWRF_WRF_RUNS must point at the data (run_N3h_r8 of n3_data.tar)
C=$1; B=$2; CFG=$3; OUT=$4; shift 4
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$C"
/usr/bin/time -f "wall %e s, maxrss %M KB" "$B/cppwrf" "$CFG" --set output.dir="$OUT" "$@" > "$OUT/run.log" 2>&1
echo "rc=$? $(tail -1 $OUT/run.log)"; grep -E 'stepping' "$OUT/run.log"
