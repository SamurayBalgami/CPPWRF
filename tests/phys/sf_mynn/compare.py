#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""compare.py WRF.bin CPP.bin NI: the MYNN surface-layer outputs, every value identical as tests/verify/exact.py (NaN
equals only NaN, -0 differs from +0)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import colcmp  # noqa: E402

names = ("pblh chs chs2 cqs cqs2 cpm znt ust ustm zol mol rmol psim psih qgh hfx qfx lh flhc flqc qsfc u10 v10 th2 t2 "
         "q2 gz1oz0 wspd br ch ck cka cd cda").split()
ni = int(sys.argv[3])
for p in sys.argv[1:3]:
    if os.path.getsize(p) != len(names) * ni * 8:
        sys.exit(f"{p}: {os.path.getsize(p)} bytes, expected {len(names) * ni * 8}")
a = colcmp.read_exact(sys.argv[1], np.float64, len(names) * ni).reshape(len(names), ni)
b = colcmp.read_exact(sys.argv[2], np.float64, len(names) * ni).reshape(len(names), ni)
bad = sum(not colcmp.same(k, x, y, quiet=True) for k, x, y in zip(names, a, b))
print(f"{'PASS' if bad == 0 else 'FAIL'}: {len(names) - bad}/{len(names)} fields identical")
sys.exit(1 if bad else 0)
