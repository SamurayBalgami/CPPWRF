#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Compares the outputs of driver_p3.f90 and driver_p3.cpp bit for bit (tolerance 0: -0 differs from +0, NaN
equals only a NaN), field by field and step by step."""
import struct
import sys

import numpy as np

F3 = "th qv qc qr qnr qi qni qir qib qnc th_old qv_old refl_10cm re_cloud re_ice v_ice d_ice rho_ice".split()
F2 = "rainnc rainncv sr snownc snowncv".split()


def main():
    inp, a, b = sys.argv[1:4]
    n, nz, steps = struct.unpack("<iii", open(inp, "rb").read(12))
    x = np.fromfile(a, "<f8")
    y = np.fromfile(b, "<f8")
    per = len(F3) * n * nz + len(F2) * n
    if x.size != steps * per or y.size != steps * per:
        print(f"FAIL: sizes {x.size} {y.size}, expected {steps * per}")
        return 1
    bad = 0
    for s in range(steps):
        off = s * per
        for name, cnt in [(f, n * nz) for f in F3] + [(f, n) for f in F2]:
            u = x[off:off + cnt].view("<u8"); v = y[off:off + cnt].view("<u8")
            xu, yv = x[off:off + cnt], y[off:off + cnt]
            same = (u == v) | (np.isnan(xu) & np.isnan(yv))
            if not same.all():
                idx = np.flatnonzero(~same)
                d = np.nanmax(np.abs(xu[idx] - yv[idx])) if idx.size else 0
                i0 = idx[0]
                print(f"step {s + 1} {name}: {idx.size} differ, first index {i0} (i={i0 % n}, k={i0 // n}) "
                      f"wrf={xu[i0]!r} cpp={yv[i0]!r}, max |d| {d:g}")
                bad += 1
            off += cnt
    nonzero = {f: 0 for f in F3}
    print("PASS" if bad == 0 else f"FAIL: {bad} field/step pairs differ")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
