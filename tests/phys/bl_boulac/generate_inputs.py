#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""generate_inputs.py --cols N --nz K --dt DT --out FILE [--qc 0|1] [--seed S]: random BouLac PBL columns (stable,
neutral and convective profiles, surface fluxes of both signs, TKE from the minimum to a few m2 s-2) in the stream
layout of driver_boulac.f90 / driver_boulac.cpp (int32 n, nz; float64 dt, cp, g; 3D arrays (i fastest) u, v, th, rho,
qv, qc, dz8w, tke; 2D ust, hfx, qfx)."""
import argparse
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--cols", type=int, default=2000)
p.add_argument("--nz", type=int, default=44)
p.add_argument("--dt", type=float, default=60.0)
p.add_argument("--qc", type=int, default=1)
p.add_argument("--seed", type=int, default=1)
p.add_argument("--out", required=True)
a = p.parse_args()
rng = np.random.default_rng(a.seed)
n, nz = a.cols, a.nz

dz = np.empty((nz, n))
dz0 = rng.uniform(15.0, 60.0, n)
for k in range(nz):
    dz[k] = dz0 * (1.0 + 0.09 * k) * rng.uniform(0.97, 1.03, n)
z = np.cumsum(dz, axis=0) - 0.5 * dz
zi = rng.uniform(100.0, 2500.0, n)  # mixed-layer depth of the profile
lapse = rng.uniform(0.003, 0.01, n)
th = np.empty((nz, n))
th0 = rng.uniform(280.0, 305.0, n)
kind = rng.integers(0, 4, n)  # 0 convective, 1 stable, 2 neutral, 3 noisy
for k in range(nz):
    above = np.maximum(z[k] - zi, 0.0)
    t = th0 + lapse * above + 3.0 * (z[k] > zi)
    t = np.where(kind == 1, th0 + 0.02 * np.minimum(z[k], 300.0) + lapse * z[k], t)
    t = np.where(kind == 2, th0 + 0.0 * z[k], t)
    t = np.where(kind == 3, th0 + lapse * z[k] + rng.normal(0.0, 0.4, n), t)
    t = np.where(kind == 0, t - 0.8 * np.exp(-z[k] / 50.0), t)
    th[k] = t
th[:, :8] = 300.0  # some columns with equal thetas (the bbb = 0 and th(izz) = th(iz) branches)
u = rng.normal(5.0, 4.0, n)[None, :] + np.log1p(z / 10.0) * rng.uniform(0.5, 2.0, n)[None, :] + \
    rng.normal(0.0, 0.5, (nz, n))
v = rng.normal(0.0, 4.0, n)[None, :] + rng.normal(0.0, 0.5, (nz, n))
rho = 1.2 * np.exp(-z / 8500.0) * rng.uniform(0.98, 1.02, n)[None, :]
qv = np.maximum(0.012 * np.exp(-z / 2500.0) * rng.uniform(0.3, 1.2, n)[None, :], 1e-6)
qc = np.where(rng.uniform(size=(nz, n)) < 0.15, rng.uniform(0.0, 5e-4, (nz, n)), 0.0) if a.qc else np.zeros((nz, n))
tke = np.where(rng.uniform(size=(nz, n)) < 0.2, 0.0001, rng.lognormal(-1.5, 1.5, (nz, n)))
tke[:, 8:12] = 0.0001  # first-step columns (boulacinit)
ust = rng.uniform(0.02, 1.0, n)
hfx = rng.normal(80.0, 150.0, n)
hfx[kind == 1] = -np.abs(hfx[kind == 1]) * 0.3
hfx[12:16] = 0.0
qfx = rng.normal(5e-5, 6e-5, n)

with open(a.out, "wb") as f:
    np.array([n, nz], dtype="<i4").tofile(f)
    np.array([a.dt, 7.0 * 287.0 / 2.0, 9.81], dtype="<f8").tofile(f)
    for x in (u, v, th, rho, qv, qc, dz, tke):
        np.ascontiguousarray(x, dtype="<f8").tofile(f)
    for x in (ust, hfx, qfx):
        np.ascontiguousarray(x, dtype="<f8").tofile(f)
