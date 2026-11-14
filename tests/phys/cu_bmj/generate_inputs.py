#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Columns for the BMJ (cu_physics = 2) differential test: soundings of deep convection (warm, moist, conditionally
unstable), shallow convection (a moist boundary layer under a dry inversion, a saturated layer under an inversion), dry and stable columns, cold columns
(freezing level near the ground), land and sea points and CLDEFI between EFIMN and 1.

Binary file (little endian): int32 ncols, nz, stepcu, feedback; float64 dt; then float64 arrays, the level index
fastest within a column: t, qv, p, pi, rho, dz8w (nz per column), p8w (nz + 1 per column), xland, cldefi (one per
column).
"""
import argparse
import struct

import numpy as np

R_D, CP, G = 287.0, 7.0 * 287.0 / 2.0, 9.81


def sounding(rng, kind, nz, ptop):
    psfc = rng.uniform(95000.0, 102000.0)
    if kind == "deep":
        tsfc, gam, rh0, rhtop, zinv = rng.uniform(296, 306), rng.uniform(6.0e-3, 7.5e-3), rng.uniform(.8, .97), .5, 0
    elif kind == "shallow":
        tsfc, gam, rh0, rhtop, zinv = rng.uniform(290, 300), 5.5e-3, rng.uniform(.8, .95), .15, rng.uniform(1200, 2500)
    elif kind == "dry":
        tsfc, gam, rh0, rhtop, zinv = rng.uniform(285, 312), 8.0e-3, rng.uniform(.05, .3), .05, 0
    elif kind == "cold":
        tsfc, gam, rh0, rhtop, zinv = rng.uniform(266, 276), rng.uniform(5.0e-3, 7.0e-3), rng.uniform(.7, .95), .4, 0
    elif kind == "wetshallow":  # a saturated cloud layer under an inversion (BMJ raises the shallow cloud top)
        tsfc, gam, rh0, rhtop, zinv = rng.uniform(288, 300), rng.uniform(6.0e-3, 7.5e-3), rng.uniform(.9, 1.0), .1, \
            rng.uniform(2500, 4000)
    else:  # mixed: anything in between
        tsfc, gam, rh0, rhtop, zinv = rng.uniform(280, 305), rng.uniform(4.5e-3, 8e-3), rng.uniform(.3, .97), \
            rng.uniform(.05, .7), rng.choice([0, rng.uniform(800, 3000)])
    # pressure at the w levels: hybrid-like eta from psfc to ptop
    eta = np.linspace(1.0, 0.0, nz + 1) ** rng.uniform(1.2, 1.6)
    p8w = ptop + eta * (psfc - ptop)
    p8w += rng.uniform(-30, 30, nz + 1) * (eta * (1 - eta))
    p = 0.5 * (p8w[:-1] + p8w[1:])
    # temperature by height (hydrostatic, dry)
    t = np.zeros(nz)
    z = 0.0
    tl = tsfc
    for k in range(nz):
        tv = tl
        dz = R_D * tv / G * np.log((p8w[k]) / p8w[k + 1])
        zm = z + 0.5 * dz
        if zinv and zinv < zm < zinv + 600:
            tk = tsfc - gam * zinv + 5.0 * (zm - zinv) / 600.0
        elif zinv and zm >= zinv + 600:
            tk = tsfc - gam * zinv + 5.0 - gam * (zm - zinv - 600)
        else:
            tk = tsfc - gam * zm
        t[k] = max(tk, 195.0 + rng.uniform(0, 15))
        tl = t[k]
        z += dz
    t += rng.normal(0, 0.3, nz)
    es = 611.2 * np.exp(17.67 * (t - 273.15) / (t - 29.65))
    qs = 0.622 * es / np.maximum(p - es, 1.0)
    frac = (psfc - p) / (psfc - ptop)
    rh = rh0 + (rhtop - rh0) * np.minimum(1.0, frac * 2.0)
    if kind == "wetshallow":
        rh = np.where(p > psfc - rng.uniform(15000, 30000), rh0 + rng.uniform(-.03, .03, nz), rhtop)
    if zinv and kind != "wetshallow":
        rh = np.where(p < psfc - 0.8 * (psfc - ptop) * 0.25, rh * 0.4, rh)
    rh = np.clip(rh + rng.normal(0, 0.03, nz), 0.0, 1.05)
    qv = np.maximum(rh * qs, 1.e-7)
    qv[p < 10000.0] = 2.e-6
    if rng.uniform() < 0.02:
        qv[:] = 0.0  # EPSQ floor
    pi = (p / 1.e5) ** (R_D / CP)
    rho = p / (R_D * t * (1.0 + 0.608 * qv))
    dz8w = (p8w[:-1] - p8w[1:]) / (rho * G) * rng.uniform(0.999, 1.001)
    return t, qv, p, pi, rho, dz8w, p8w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cols", type=int, default=2000)
    ap.add_argument("--nz", type=int, default=44)
    ap.add_argument("--stepcu", type=int, default=1)
    ap.add_argument("--dt", type=float, default=60.0)
    ap.add_argument("--feedback", type=int, default=0)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    kinds = ["deep", "deep", "shallow", "wetshallow", "wetshallow", "dry", "cold", "mixed", "mixed"]
    cols = []
    for c in range(a.cols):
        cols.append(sounding(rng, kinds[c % len(kinds)], a.nz, 5000.0))
    xland = np.where(rng.uniform(size=a.cols) < 0.4, 2.0, 1.0)
    cldefi = np.where(rng.uniform(size=a.cols) < 0.3, 0.6, rng.uniform(0.2, 1.0, a.cols))
    with open(a.out, "wb") as f:
        f.write(struct.pack("<4i", a.cols, a.nz, a.stepcu, a.feedback))
        f.write(struct.pack("<d", a.dt))
        for n in range(7):
            np.concatenate([c[n] for c in cols]).astype("<f8").tofile(f)
        xland.astype("<f8").tofile(f)
        cldefi.astype("<f8").tofile(f)


if __name__ == "__main__":
    main()
