#!/usr/bin/env bash
# Job slots for heavy work on a shared workstation (sized for 12 cores, 15 GB RAM, 1 GPU).
#   tools/heavy.sh <command...>              run a heavy job (build, long run, WRF, tests) in one
#                                            of 4 slots; waits while all slots are busy
#   tools/heavy.sh --light <command...>      quick job (python check, ncdump, one small compile,
#                                            a few model steps): no slot, runs at once with 2 threads
#   tools/heavy.sh --exclusive <command...>  take ALL slots (benchmarks, timings, GPU profiling).
#                                            Waits until every slot is free at the same moment and
#                                            never holds some slots while waiting; for its first
#                                            2 minutes of waiting new heavy jobs hold back.
# Defaults OMP_NUM_THREADS=3 (light: 2) and OMP_WAIT_POLICY=passive unless already set.
set -uo pipefail
LOCKDIR=/tmp/cppwrf_slots
NSLOTS=4
mkdir -p "$LOCKDIR"
[ $# -ge 1 ] || { echo "usage: tools/heavy.sh [--light|--exclusive] <command...>" >&2; exit 2; }
if [ "$1" = "--light" ]; then
  shift
  export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}" OMP_WAIT_POLICY="${OMP_WAIT_POLICY:-passive}"
  exec nice -n 5 "$@"
fi
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-3}" OMP_WAIT_POLICY="${OMP_WAIT_POLICY:-passive}"
# Nested call (a script run under heavy.sh that calls heavy.sh again): the outer call already
# holds a slot, so run directly.
if [ -n "${CPPWRF_IN_SLOT:-}" ]; then
  [ "$1" = "--exclusive" ] && shift
  exec "$@"
fi
export CPPWRF_IN_SLOT=1

if [ "$1" = "--exclusive" ]; then
  shift
  flag="$LOCKDIR/exclusive_waiting.$$"
  touch "$flag"
  trap 'rm -f "$flag"' EXIT
  while true; do
    got=()
    for i in $(seq 0 $((NSLOTS - 1))); do
      fd=$((10 + i))
      eval "exec $fd>\"$LOCKDIR/slot$i\""
      if flock -n "$fd"; then got+=("$fd"); else break; fi
    done
    [ "${#got[@]}" -eq "$NSLOTS" ] && break
    for fd in "${got[@]}"; do flock -u "$fd"; eval "exec $fd>&-"; done  # never hold a partial set
    sleep 5
  done
  rm -f "$flag"
  echo "[heavy] exclusive: $*" >&2
  "$@"
  exit $?
fi

# A waiting --exclusive job holds back new jobs for at most 2 minutes (from its request).
exclusive_blocks() {
  local now f
  now=$(date +%s)
  for f in "$LOCKDIR"/exclusive_waiting.*; do
    [ -e "$f" ] || continue
    if ! kill -0 "${f##*.}" 2>/dev/null; then rm -f "$f"; continue; fi  # stale flag
    [ $((now - $(stat -c %Y "$f" 2>/dev/null || echo 0))) -lt 120 ] && return 0
  done
  return 1
}
while true; do
  if ! exclusive_blocks; then
    for i in $(seq 0 $((NSLOTS - 1))); do
      exec 9>"$LOCKDIR/slot$i"
      if flock -n 9; then
        echo "[heavy] slot $i: $*" >&2
        "$@"
        exit $?
      fi
    done
  fi
  sleep 5
done
