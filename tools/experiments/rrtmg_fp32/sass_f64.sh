#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# FP64 SASS instructions (DADD, DMUL, DFMA, DSETP, F2F to/from F64, DMNMX, MUFU *64H) per kernel of a CUDA binary;
# total and in the RRTMG kernels (names with rrtmg_lw_batch / rrtmg_sw_batch). $@: cppwrf binaries
for x in "$@"; do
  cuobjdump -sass "$x" 2>/dev/null | awk '/Function :/{fn=$3} /[[:space:]](DADD|DMUL|DFMA|DSETP|F2F\.F64|F2F\.F32\.F64|DMNMX|MUFU\.RCP64H|MUFU\.RSQ64H)/{n[fn]++} END{for(f in n) print n[f], f}' > /tmp/f64.$$
  echo "$x: f64 instructions $(awk '{s+=$1}END{print s+0}' /tmp/f64.$$), in RRTMG kernels $(grep -E 'rrtmg_(lw|sw)_batch' /tmp/f64.$$ | awk '{s+=$1}END{print s+0}')"
  grep -E 'rrtmg_(lw|sw)_batch' /tmp/f64.$$ | sort -nr | head -5 | cut -c1-140
done; rm -f /tmp/f64.$$
