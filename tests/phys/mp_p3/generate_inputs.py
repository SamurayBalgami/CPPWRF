#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Columns for the P3 (mp_physics 50 / 51) differential test: warm, mixed-phase and cold soundings near and above
saturation with cloud water, rain, ice (rime mass and volume across the rime-density limits), droplet number for
the predicted-Nc variant, zeros and values at the clipping thresholds (qsmall, 1e-8, 1e-12), strong updraft-like
th_old / qv_old (the 'A' term) or the cold start (th_old = qv_old = 0, itimestep 1)."""
import argparse
import struct

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cols", type=int, default=1001)
    ap.add_argument("--nz", type=int, default=41)
    ap.add_argument("--dt", type=float, default=60.0)
    ap.add_argument("--steps", type=int, default=3)
    ap.add_argument("--start-it", type=int, default=1)
    ap.add_argument("--predict-nc", type=int, default=0)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    n, nz = a.cols, a.nz
    zw = 18000.0 * np.linspace(0, 1, nz + 1) ** 1.25
    dz = np.diff(zw)
    zm = 0.5 * (zw[1:] + zw[:-1])
    th = np.zeros((nz, n)); qv = np.zeros_like(th); qc = np.zeros_like(th); qr = np.zeros_like(th)
    qnr = np.zeros_like(th); qi = np.zeros_like(th); qni = np.zeros_like(th); qir = np.zeros_like(th)
    qib = np.zeros_like(th); qnc = np.zeros_like(th); p = np.zeros_like(th); dzz = np.zeros_like(th)
    tho = np.zeros_like(th); qvo = np.zeros_like(th)

    def esat(t, ice):
        tc = t - 273.15
        return np.where(ice & (t < 273.15), 611.2 * np.exp(22.46 * tc / (tc + 272.62)),
                        611.2 * np.exp(17.62 * tc / (tc + 243.12)))

    for c in range(n):
        ts = rng.uniform(250.0, 306.0)
        t = np.maximum(ts - 0.0065 * zm, rng.uniform(195.0, 215.0))
        ps = rng.uniform(85000.0, 102000.0)
        pc = ps * np.exp(-zm / 7500.0)
        rh = np.clip(rng.uniform(0.6, 1.08, nz) + rng.normal(0, 0.03, nz), 0.05, 1.2)
        es = esat(t, np.ones(nz, bool))
        qs = 0.622 * es / (pc - es)
        q = rh * qs
        cold = t < 273.15
        cloudy = rh > rng.uniform(0.85, 1.0)
        lc = np.where(cloudy & (t > 235.0), 10 ** rng.uniform(-6.5, -2.7, nz), 0.0)
        lr = np.where(rng.random(nz) < 0.5, 10 ** rng.uniform(-7.5, -2.5, nz), 0.0)
        li = np.where(cold & (rng.random(nz) < 0.7), 10 ** rng.uniform(-8.5, -2.8, nz), 0.0)
        # thresholds and tiny amounts
        for arr in (lc, lr, li):
            m = rng.random(nz)
            arr[m < 0.04] = 1.0e-14 * rng.uniform(0.5, 2.0)
            arr[(m >= 0.04) & (m < 0.08)] = 10 ** rng.uniform(-13, -8)
        rhoa = pc / (287.0 * t)
        ni = li / (10 ** rng.uniform(-12.0, -6.0, nz))
        ni[li == 0] = 0.0
        frim = np.where(rng.random(nz) < 0.6, rng.uniform(0.0, 1.05, nz), 0.0)
        rim = li * frim
        rhor = rng.uniform(20.0, 950.0, nz)
        bim = np.where(rim > 0, rim / rhor, 0.0)
        bim[rng.random(nz) < 0.03] = 1.0e-16
        dr = 10 ** rng.uniform(-4.3, -2.7, nz)
        nrr = np.where(lr > 0, lr / (np.pi / 6 * 1000.0 * dr ** 3), 0.0)
        ncc = np.where(lc > 0, rng.uniform(20e6, 600e6, nz) / rhoa, 0.0)
        if c % 17 == 0:  # clear column: only vapour (the skipped columns)
            lc[:] = lr[:] = li[:] = ni[:] = rim[:] = bim[:] = nrr[:] = ncc[:] = 0.0
            q *= 0.5
        theta = t * (1.0e5 / pc) ** (287.15 / 1005.0)
        th[:, c] = theta; qv[:, c] = q; qc[:, c] = lc; qr[:, c] = lr; qnr[:, c] = nrr
        qi[:, c] = li; qni[:, c] = ni; qir[:, c] = rim; qib[:, c] = bim; qnc[:, c] = ncc
        p[:, c] = pc; dzz[:, c] = dz
        if a.start_it > 1:  # the previous step's values: cooling (ascent) or warming
            tho[:, c] = theta + rng.normal(0.0, 0.3, nz)
            qvo[:, c] = q * (1.0 + rng.normal(0.0, 0.02, nz))
    with open(a.out, "wb") as f:
        f.write(struct.pack("<iiiid", n, nz, a.steps, a.predict_nc, a.dt))
        f.write(struct.pack("<i", a.start_it))
        for arr in (th, qv, qc, qr, qnr, qi, qni, qir, qib, qnc, tho, qvo, p, dzz):
            f.write(np.ascontiguousarray(arr, dtype="<f8").tobytes())  # (k, i): i fastest = Fortran (i, k)


if __name__ == "__main__":
    main()
