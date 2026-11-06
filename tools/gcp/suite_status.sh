#!/usr/bin/env bash
# The result of a CPU suite run (tools/gcp/gcp_run.sh, tools/colab/colab_suite.sh): prints the PASS / FAIL / SKIP counts
# and lines of OUT's logs and exits with 0 only for a complete suite without FAIL or SKIP (wrf_suite.py's SUITE line and
# exit 0), passing translator tests and comparator self-test, an empty setup.out and, for the full suite (CASES empty),
# files, ideal files, differential and ctest tests without FAIL or SKIP.
#
#   suite_status.sh OUT CASES [FILES]   (FILES=1: a targeted suite with the output files tests, colab_suite.sh --files)
set -euo pipefail
OUT=$1 CASES=$2 FILES=${3:-0}
for k in PASS FAIL SKIP; do echo "$k $(grep -c "^$k " "$OUT/suite.log" 2>/dev/null || true)"; done
grep -E "^(FAIL|SKIP|ERROR) " "$OUT/suite.log" 2>/dev/null || true
bad=0
cat "$OUT/setup.out" 2>/dev/null | grep . && bad=1
grep -q "^exit 0$" "$OUT/unit.log" 2>/dev/null || { echo "[gcp_run] translator tests: $(tail -1 "$OUT/unit.log" 2>/dev/null || echo 'no unit.log')"; bad=1; }
grep -q "^exit 0$" "$OUT/selftest.log" 2>/dev/null || { echo "[gcp_run] comparator self-test: $(grep -E "WRONG|SELFTEST" "$OUT/selftest.log" 2>/dev/null || echo 'no selftest.log')"; bad=1; }
if ! grep -q "^SUITE " "$OUT/suite.log" 2>/dev/null; then echo "[gcp_run] INCOMPLETE: the suite did not end (no SUITE line)"; bad=1
else grep "^SUITE " "$OUT/suite.log"; grep -q "^exit 0$" "$OUT/suite.log" || bad=1; fi
if [ -z "$CASES" ] || [ "$FILES" = 1 ]; then  # the full suite: the output files tests too (remote_suite.sh)
  for f in files idealfiles; do
    if [ ! -f "$OUT/$f.log" ]; then echo "[gcp_run] INCOMPLETE: no $f.log"; bad=1; continue; fi
    echo "$f: PASS $(grep -c '^PASS ' "$OUT/$f.log") FAIL $(grep -c '^FAIL ' "$OUT/$f.log") SKIP $(grep -c '^SKIP ' "$OUT/$f.log")"
    grep -E "^(FAIL|SKIP) " "$OUT/$f.log" | cut -c1-200 || true  # none: grep fails (set -e, pipefail)
    grep -E "^exit [1-9]" "$OUT/$f.log" && bad=1
  done
fi
if [ -z "$CASES" ]; then
  for f in diff ctest; do
    if [ ! -f "$OUT/$f.log" ]; then echo "[gcp_run] INCOMPLETE: no $f.log"; bad=1; continue; fi
    echo "$f: $(grep -E 'tests passed' "$OUT/$f.log" || echo 'no ctest summary')"
    grep -E '^\s*[0-9]+ - .*\((Failed|Timeout|SEGFAULT|Subprocess|Skipped|Not Run)|^SKIPPED ' "$OUT/$f.log" || true
    grep -qx "exit 0" "$OUT/$f.log" || bad=1
  done
fi
exit $bad
