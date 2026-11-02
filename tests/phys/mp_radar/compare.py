#!/usr/bin/env python3
"""Compares driver_radar.f90's and driver_radar_cpu's outputs bit for bit: the radar_init tables of mp_wsm6_init
and WSM6's REFL_10CM.   compare.py FORTRAN_OUT CPP_OUT N_COLS NZ"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import colcmp  # noqa: E402

NB = 50
PARTS = [("xxDs", NB), ("xdts", NB), ("xxDg", NB), ("xdtg", NB), ("simpson", NB + 1),
         ("k_w pi5 lamda4", 3), ("m_w_0 m_i_0", 4),
         ("xcre", 4), ("xcse", 4), ("xcge", 4), ("xcrg", 4), ("xcsg", 4), ("xcgg", 4),
         ("xam_r ... xogg2", 19)]


def main():
    fort, cpp, ncols, nz = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    f = np.fromfile(fort, dtype=np.float64)
    c = np.fromfile(cpp, dtype=np.float64)
    parts = PARTS + [("refl_10cm", ncols * nz)]
    n = sum(m for _, m in parts)
    if f.size != n or c.size != n:
        print(f"FAIL: {f.size} / {c.size} values, expected {n}")
        return 1
    ok, o = True, 0
    for name, m in parts:
        ok &= colcmp.same(name, f[o:o + m], c[o:o + m])
        o += m
    r = f[-ncols * nz:]
    print(f"refl_10cm: {np.sum(r > -35.)} values above -35 dBZ, max {r.max():.3f}")
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
