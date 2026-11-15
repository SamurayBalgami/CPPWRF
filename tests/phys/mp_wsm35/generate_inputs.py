#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Slabs for the WSM5 / WSM3 differential test (ni x nj columns, nz levels): warm, mixed-phase and cold soundings
below and above saturation with cloud water, rain, cloud ice and snow at every height (snow and ice in warm layers:
melting; rain and cloud water in cold layers: freezing, below -40 C: homogeneous freezing), amounts at the thresholds
(qmin 1e-15, qcrmin 1e-9, qc0) and zeros, clear columns, updrafts and downdrafts (WSM3's freezing-level term), the
precipitation accumulations of an earlier step. WSM3 reads qc (qci) and qr (qrs) only."""
import argparse
import struct

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ni", type=int, default=41)
    ap.add_argument("--nj", type=int, default=25)
    ap.add_argument("--nz", type=int, default=44)
    ap.add_argument("--dt", type=float, default=60.0)
    ap.add_argument("--steps", type=int, default=3)
    ap.add_argument("--scheme", type=int, default=5)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    ni, nj, nz = a.ni, a.nj, a.nz
    shape = (nj, nz, ni)  # C order (j, k, i) = Fortran (i, k, j)
    f = {n: np.zeros(shape) for n in "th q qc qr qi qs den pii p delz w".split()}
    rain = rng.uniform(0.0, 20.0, (nj, ni)) * (rng.random((nj, ni)) < 0.5)
    snow = rng.uniform(0.0, 5.0, (nj, ni)) * (rng.random((nj, ni)) < 0.5)

    def esat(t, ice):
        tc = t - 273.15
        return np.where(ice & (t < 273.15), 611.2 * np.exp(22.46 * tc / (tc + 272.62)),
                        611.2 * np.exp(17.62 * tc / (tc + 243.12)))

    for j in range(nj):
        for i in range(ni):
            ztop = rng.uniform(14000.0, 20000.0)
            zw = ztop * np.linspace(0, 1, nz + 1) ** rng.uniform(1.0, 1.4)
            dz = np.diff(zw)
            zm = 0.5 * (zw[1:] + zw[:-1])
            ts = rng.uniform(255.0, 305.0)
            t = np.maximum(ts - rng.uniform(0.005, 0.0085) * zm, rng.uniform(195.0, 215.0))
            t += rng.normal(0.0, 0.4, nz)
            ps = rng.uniform(80000.0, 103000.0)
            pc = ps * np.exp(-zm / rng.uniform(7000.0, 8200.0))
            rh = np.clip(rng.uniform(0.5, 1.1, nz) + rng.normal(0, 0.04, nz), 0.02, 1.25)
            es = esat(t, np.ones(nz, bool))
            qsat = 0.622 * es / (pc - es)
            qv = rh * qsat
            cloudy = rh > rng.uniform(0.8, 1.0)
            lc = np.where(cloudy & (rng.random(nz) < 0.8), 10 ** rng.uniform(-6.5, -2.6, nz), 0.0)
            lr = np.where(rng.random(nz) < 0.5, 10 ** rng.uniform(-8.5, -2.4, nz), 0.0)
            li = np.where(rng.random(nz) < 0.5, 10 ** rng.uniform(-8.5, -3.0, nz), 0.0)
            ls = np.where(rng.random(nz) < 0.5, 10 ** rng.uniform(-8.5, -2.5, nz), 0.0)
            for arr in (lc, lr, li, ls):
                m = rng.random(nz)
                arr[m < 0.03] = 1.0e-15 * rng.uniform(0.5, 2.0)
                arr[(m >= 0.03) & (m < 0.06)] = 1.0e-9 * rng.uniform(0.5, 2.0)
                arr[(m >= 0.06) & (m < 0.08)] = 10 ** rng.uniform(-14, -9)
            if (i + 3 * j) % 13 == 0:  # clear column
                lc[:] = lr[:] = li[:] = ls[:] = 0.0
            if (i + j) % 11 == 0:  # heavy precipitation
                lr *= 10.0
                ls *= 10.0
            if a.scheme == 3:  # qci / qrs: cloud water or ice, rain or snow
                lc = lc + li
                lr = lr + ls
            pi_ = (pc / 1.0e5) ** (287.0 / 1004.5)
            f["th"][j, :, i] = t / pi_
            f["pii"][j, :, i] = pi_
            f["p"][j, :, i] = pc
            f["q"][j, :, i] = qv
            f["qc"][j, :, i] = lc
            f["qr"][j, :, i] = lr
            f["qi"][j, :, i] = li
            f["qs"][j, :, i] = ls
            f["den"][j, :, i] = pc / (287.0 * t * (1.0 + 0.608 * qv))
            f["delz"][j, :, i] = dz
            f["w"][j, :, i] = rng.normal(0.0, 1.5, nz) * (rng.random() < 0.8)
    with open(a.out, "wb") as o:
        o.write(struct.pack("<iiiiid", ni, nj, nz, a.steps, a.scheme, a.dt))
        for n in "th q qc qr qi qs den pii p delz w".split():
            o.write(np.ascontiguousarray(f[n], dtype="<f8").tobytes())
        for arr in (rain, snow):
            o.write(np.ascontiguousarray(arr, dtype="<f8").tobytes())


if __name__ == "__main__":
    main()
