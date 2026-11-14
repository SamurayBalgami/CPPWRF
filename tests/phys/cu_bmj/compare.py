#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Compares the outputs of driver_bmj.f90 (WRF) and driver_bmj.cpp: every value bit for bit (--tol 0) or within a
relative tolerance; prints per field the differing points and the largest difference."""
import argparse
import sys

import numpy as np

JTB, ITB = 134, 76


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("wrf")
    ap.add_argument("cpp")
    ap.add_argument("--input", required=True)
    ap.add_argument("--tol", type=float, default=0.0)
    a = ap.parse_args()
    ncols, nz = np.fromfile(a.input, dtype="<i4", count=2)
    w = np.fromfile(a.wrf, dtype="<f8")
    c = np.fromfile(a.cpp, dtype="<f8")
    if w.size != c.size:
        print(f"FAIL: sizes differ ({w.size} vs {c.size})")
        return 1
    names = [("QS0", JTB), ("SQS", JTB), ("PTBL", JTB * ITB)]
    names += [(n, nz * ncols) for n in ("RTHCUTEN", "RQVCUTEN", "CCLDFRA", "QCCONV", "QICONV")]
    names += [(n, ncols) for n in ("RAINCV", "PRATEC", "CUTOP", "CUBOT", "CLDEFI", "CONVCLD")]
    off, bad = 0, 0
    for n, m in names:
        x, y = w[off:off + m], c[off:off + m]
        off += m
        if a.tol == 0:
            d = (x.view("<u8") != y.view("<u8")) & ~(np.isnan(x) & np.isnan(y))
        else:
            d = np.abs(x - y) > a.tol * np.maximum(np.abs(x), 1e-30)
        nonzero = int(np.count_nonzero(x))
        if d.any():
            bad += 1
            print(f"  {n}: {int(d.sum())} of {m} differ, max |d| {np.nanmax(np.abs(x - y)):.3e} (nonzero {nonzero})")
        else:
            print(f"  {n}: identical (nonzero {nonzero})")
    act = int(np.count_nonzero(w[off - 6 * ncols:off - 5 * ncols]))
    cub = w[off - 3 * ncols:off - 2 * ncols]
    print(f"  columns with convective rain: {act} of {ncols}; with a cloud (CUBOT <= nz): {int((cub <= nz).sum())}")
    if bad:
        print("FAIL")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
