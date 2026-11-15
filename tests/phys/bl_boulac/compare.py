#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""compare.py WRF.bin CPP.bin --cols N --nz K: the BouLac outputs bit for bit (NaN equals only NaN, -0 differs from
+0); exit 1 on any difference."""
import argparse
import sys
import numpy as np

NAMES = ["tke", "dlk", "wu", "wv", "wt", "wq", "exch_h", "exch_m", "rublten", "rvblten", "rthblten", "rqvblten",
         "rqcblten"]
p = argparse.ArgumentParser()
p.add_argument("fortran")
p.add_argument("cpp")
p.add_argument("--cols", type=int, required=True)
p.add_argument("--nz", type=int, required=True)
a = p.parse_args()
f = np.fromfile(a.fortran, dtype="<f8")
c = np.fromfile(a.cpp, dtype="<f8")
n3 = a.cols * a.nz
if f.size != c.size or f.size != len(NAMES) * n3 + a.cols:
    print(f"ERROR: sizes {f.size} / {c.size}, expected {len(NAMES) * n3 + a.cols}")
    sys.exit(1)
ok = True
for idx, name in enumerate(NAMES + ["pblh"]):
    lo = idx * n3
    hi = lo + (n3 if idx < len(NAMES) else a.cols)
    x, y = f[lo:hi], c[lo:hi]
    diff = x.view("<u8") != y.view("<u8")
    nd = int(diff.sum())
    mx = float(np.nanmax(np.abs(x[diff] - y[diff]))) if nd else 0.0
    nonzero = int((x != 0).sum())
    print(f"{name:<10} nonzero {nonzero:>8}  differing {nd:>8}  max|d| {mx:.3e}  {'PASS' if nd == 0 else 'FAIL'}")
    ok = ok and nd == 0
print("BouLac column test", "PASSED" if ok else "FAILED")
sys.exit(0 if ok else 1)
