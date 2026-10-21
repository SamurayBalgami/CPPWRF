#!/usr/bin/env python3
"""compare.py WRF.bin CPP.bin NI [--tol 0]: per-field max |difference| of the MYNN surface-layer outputs."""
import sys
import numpy as np
names = "pblh chs chs2 cqs cqs2 cpm znt ust ustm zol mol rmol psim psih qgh hfx qfx lh flhc flqc qsfc u10 v10 th2 t2 q2 gz1oz0 wspd br ch ck cka cd cda".split()
ni = int(sys.argv[3])
tol = float(sys.argv[sys.argv.index("--tol") + 1]) if "--tol" in sys.argv else 0.0
a = np.fromfile(sys.argv[1], dtype=np.float64).reshape(len(names), ni)
b = np.fromfile(sys.argv[2], dtype=np.float64).reshape(len(names), ni)
bad = 0
for k, x, y in zip(names, a, b):
    d = np.abs(x - y)
    same_nan = np.isnan(x) & np.isnan(y)
    d[same_nan] = 0
    nd = int((d > tol).sum() + (np.isnan(x) ^ np.isnan(y)).sum())
    if nd:
        bad += 1
        i = int(np.nanargmax(d))
        print(f"DIFF {k:7s} {nd:6d} points, max {d[i]:.3e} at {i}: wrf {x[i]!r} cpp {y[i]!r}")
print(f"{'PASS' if bad == 0 else 'FAIL'}: {len(names) - bad}/{len(names)} fields identical")
sys.exit(1 if bad else 0)
