#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The proof of src/core/libmvec_x86.hpp (glibc libmvec's vector log / log10 with Intel's rcpps) on a Google Colab
# runtime: tools/colab/mvec_body.sh in its Debian 13 chroot (glibc 2.41 as the GCP VMs). No data tars.
#
#   tools/colab/colab_mvec.sh [--gpu L4|G4|CPU|A100|T4] [--n N] [--out DIR] [--rev REV] [--session NAME] [--keep]
#
# --gpu: the runtime (default L4: Intel; G4 and CPU are AMD), any CPU vendor; --n: random inputs per function for
# tools/libmvec_proof.cpp (default 10^10; on Intel all must equal libmvec); --rev: the commit (default HEAD).
# Results in OUT (default ${CPPWRF_SCRATCH:-~/cppwrf_scratch}/colab/mvec_<date>): the logs of mvec_body.sh and
# job.log; exit 0 when the unit test passed, the table equals this CPU's rcpps (Intel) and the proof found no
# difference (Intel), and both deliberate faults made the unit test fail.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
GPU=L4 N=10000000000 OUT= REV=HEAD S= KEEP=0 NPROC= ANY_CPU=1 BK=${CPPWRF_GCP_BUCKET:-none}
while [ $# -gt 0 ]; do
  case $1 in
    --gpu) GPU=$2; shift 2 ;; --n) N=$2; shift 2 ;; --out) OUT=$2; shift 2 ;; --rev) REV=$2; shift 2 ;;
    --session) S=$2; shift 2 ;; --keep) KEEP=1; shift ;;
    *) sed -n '3,13p' "$0"; exit 2 ;;
  esac
done
OUT=${OUT:-${CPPWRF_SCRATCH:-$HOME/cppwrf_scratch}/colab/mvec_$(date +%Y%m%d_%H%M%S)}
mkdir -p "$OUT"
git -C "$REPO" archive --format=tar.gz -o "$OUT/src.tar.gz" "$REV"
echo "[colab_mvec] commit $(git -C "$REPO" rev-parse --short "$REV"), $GPU, n $N"
TAG=colab_mvec
. "$REPO/tools/colab/colab_lib.sh"
cj_token() { :; }  # no data: no bucket token
cj_start "$GPU"
cj_launch "$REPO/tools/colab/mvec_body.sh" "" "$N"
cj_wait 4 proof.log
cj_results
for f in cpu.txt unit.log check.log regen.log proof.log binary.log bug_table.log bug_order.log; do
  [ -f "$OUT/$f" ] && { echo "== $f"; grep -v -E '^\+' "$OUT/$f" | tail -25; }
done
ok=1
grep -q 'test_libmvec passed' "$OUT/unit.log" 2> /dev/null || ok=0
for b in table order; do grep -q '^exit 1$' "$OUT/bug_$b.log" 2> /dev/null && grep -q 'test_libmvec FAILED' "$OUT/bug_$b.log" || ok=0; done
if grep -q GenuineIntel "$OUT/cpu.txt" 2> /dev/null; then
  grep -q '^exit 0$' "$OUT/check.log" && grep -q 'libmvec_proof: identical' "$OUT/proof.log" && grep -q '^diff exit 0$' "$OUT/regen.log" || ok=0
fi
[ "$ok" = 1 ] && echo "[colab_mvec] PASS" || { echo "[colab_mvec] FAIL"; exit 1; }
