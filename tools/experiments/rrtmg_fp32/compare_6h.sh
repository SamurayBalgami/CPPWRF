#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Statistics of the 6 h runs against the FP64 run: OUT_ROOT holds r64_6h, b32_6h, e32_6h, b32nf_6h, e32nf_6h
# (b = upstream, e = experiment, nf = built with --fmad=false); prints RMSE tables at 3 h and 6 h per domain.
D=${OUT_ROOT:?set OUT_ROOT}; T=$(dirname "$0")/../../verify/fp32_accuracy.py; O=$D/acc; rm -rf $O; mkdir -p $O
for dom in d01 d02; do
  for t in b32 e32 b32nf e32nf; do python3 -I $T compare $D/${t}_6h $D/r64_6h --domain $dom --pair $t --out $O/${t}_$dom.jsonl > /dev/null & done
  python3 -I $T compare $D/b32nf_6h $D/b32_6h --domain $dom --pair bb --out $O/bb_$dom.jsonl > /dev/null &
  wait
done
python3 -I - <<'PY'
import json, os
O = os.environ["OUT_ROOT"] + "/acc"
def L(t, dom): return {(x["var"], round(x["hours"])): x["rmse"] for x in map(json.loads, open(f"{O}/{t}_{dom}.jsonl")) if "rmse" in x}
V = ["T2","Q2","U10","V10","PSFC","TSK","SWDOWN","HFX","LH","PBLH","RAIN","T","QVAPOR","U","V","W","PH","QCLOUD","QRAIN","QICE","QSNOW","QGRAUP"]
for dom in ("d01", "d02"):
    d = {t: L(t, dom) for t in ("b32", "b32nf", "e32", "e32nf", "bb")}
    for h in (3, 6):
        print(f"\n### {dom} {h} h: RMSE vs FP64 | base(FMA) base(noFMA) | exp(FMA) exp(noFMA) | exp/base mean ratio | base FMA vs noFMA")
        rs = []
        for v in V:
            k = (v, h)
            if k not in d["b32"] or d["b32"][k] == 0: continue
            bm = (d["b32"][k] + d["b32nf"][k]) / 2; em = (d["e32"][k] + d["e32nf"][k]) / 2; r = em / bm; rs.append(r)
            print(f"| {v} | {d['b32'][k]:.3g} {d['b32nf'][k]:.3g} | {d['e32'][k]:.3g} {d['e32nf'][k]:.3g} | {r:.2f} | {d['bb'][k]:.3g} |")
        rs.sort(); print(f"median {rs[len(rs)//2]:.2f}  min {rs[0]:.2f}  max {rs[-1]:.2f}")
PY
