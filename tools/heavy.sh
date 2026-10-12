#!/usr/bin/env bash
# Shared-machine job slots for all agents (12 cores, 15 GB RAM, 1 GPU).
#   tools/heavy.sh <command...>              run a heavy job (build, long run, WRF, tests) in one
#                                            of 3 slots; waits while all slots are busy
#   tools/heavy.sh --exclusive <command...>  take ALL slots (benchmarks, timings, GPU profiling);
#                                            new heavy jobs wait until it finishes
# Defaults OMP_NUM_THREADS=4 and OMP_WAIT_POLICY=passive unless already set.
# Wrap builds, model runs and test suites with this.
set -uo pipefail
LOCKDIR=/tmp/cppwrf_slots
NSLOTS=3
mkdir -p "$LOCKDIR"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}" OMP_WAIT_POLICY="${OMP_WAIT_POLICY:-passive}"
[ $# -ge 1 ] || { echo "usage: tools/heavy.sh [--exclusive] <command...>" >&2; exit 2; }

if [ "$1" = "--exclusive" ]; then
  shift
  touch "$LOCKDIR/exclusive_waiting"
  # take every slot in a fixed order (no deadlock between exclusive jobs)
  exec 3>"$LOCKDIR/slot0" 4>"$LOCKDIR/slot1" 5>"$LOCKDIR/slot2"
  flock 3; flock 4; flock 5
  rm -f "$LOCKDIR/exclusive_waiting"
  echo "[heavy] exclusive: $*" >&2
  "$@"
  exit $?
fi

while true; do
  if [ ! -e "$LOCKDIR/exclusive_waiting" ]; then
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
