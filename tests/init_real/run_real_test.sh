#!/usr/bin/env bash
# cppwrf_real against WRF real.exe: every reference case in ${CPPWRF_WRF_RUNS}/real
# (tools/init_real/make_real_refs.sh) must give identical wrfinput_d0N and wrfbdy_d01 files.
#   usage: run_real_test.sh [BUILD_DIR] [case ...]
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
if [ -z "${CPPWRF_WRF_RUNS:-}" ] || [ ! -d "${CPPWRF_WRF_RUNS}/real" ]; then
  echo "SKIP: \${CPPWRF_WRF_RUNS}/real not found"; exit 0
fi
BUILD_DIR=${1:-${REPO_ROOT}/build}; shift || true
EXE=${BUILD_DIR}/src/init_real/cppwrf_real
SCHEMA=${REPO_ROOT}/data/wrf_schema/wrf-4.8.0.tsv
WORK=${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/real_test
cases=("$@")
[ ${#cases[@]} -eq 0 ] && for d in "${CPPWRF_WRF_RUNS}"/real/*/; do cases+=("$(basename "$d")"); done
fail=0
for c in "${cases[@]}"; do
  ref=${CPPWRF_WRF_RUNS}/real/$c
  rm -rf "$WORK"; mkdir -p "$WORK"
  if ! "$EXE" --namelist "$ref/namelist.input" --met-dir "$ref" --out-dir "$WORK" \
       --precision "$(cat "$ref/precision")" --schema "$SCHEMA" > "$WORK/run.log" 2>&1; then
    echo "FAIL $c: cppwrf_real failed"; tail -3 "$WORK/run.log"; fail=1; continue
  fi
  ok=1
  for f in "$ref"/wrfinput_d0? "$ref"/wrfbdy_d01; do
    [ -f "$f" ] || continue
    python3 "${SCRIPT_DIR}/compare_real.py" "$WORK/$(basename "$f")" "$f" > "$WORK/cmp.log" 2>&1 \
      || { ok=0; cat "$WORK/cmp.log"; }
  done
  if [ $ok = 1 ]; then echo "PASS $c"; else echo "FAIL $c"; fail=1; fi
done
rm -rf "$WORK"
exit $fail
