#!/usr/bin/env python3
"""Random PBL columns for the MYNN-EDMF differential test (record order: driver_bl_mynn.f90).

  generate_inputs.py --cols N --nz NZ --nsteps S --dt DT --dx DX [--seed 1] OUT.bin

Convective (mixed layer, superadiabatic surface layer), neutral and stable columns over land and water, with
and without boundary-layer clouds (qc), ice and snow aloft, strong resolved w in some columns, and surface fluxes
that cross the driver's limits (hfx > 1200, qfx > 9e-4).
"""
import argparse
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--cols", type=int, default=2000)
p.add_argument("--nz", type=int, default=50)
p.add_argument("--nsteps", type=int, default=1)
p.add_argument("--dt", type=float, default=60.0)
p.add_argument("--dx", type=float, default=12000.0)
p.add_argument("--seed", type=int, default=1)
p.add_argument("out")
a = p.parse_args()
r = np.random.default_rng(a.seed)
n, nz = a.cols, a.nz
g, rd, cp = 9.81, 287.0, 1004.5
kind = r.integers(0, 3, n)  # 0 convective, 1 neutral, 2 stable
xland = np.where(r.random(n) < 0.6, 1.0, 2.0)
dz0 = r.uniform(15.0, 60.0, n)
growth = r.uniform(1.06, 1.12, n)
dz = np.minimum(dz0[:, None] * growth[:, None] ** np.arange(nz)[None, :], 800.0)
zf = np.concatenate([np.zeros((n, 1)), np.cumsum(dz, axis=1)], axis=1)
zm = 0.5 * (zf[:, 1:] + zf[:, :-1])
ps = r.uniform(70000.0, 103000.0, n)
t1 = r.uniform(240.0, 305.0, n)
h = np.where(kind == 0, r.uniform(300.0, 2500.0, n), np.where(kind == 1, r.uniform(150.0, 800.0, n), r.uniform(50.0, 300.0, n)))
lapse = np.where(zm < 11000.0, 0.0065, 0.0)
T = t1[:, None] - np.cumsum(lapse * dz, axis=1) + lapse[:, :1] * dz[:, :1]
# boundary layer: mixed (dry adiabat) below h for convective, inversion for stable
in_bl = zm < h[:, None]
T = np.where((kind[:, None] == 0) & in_bl, t1[:, None] - 0.0098 * (zm - zm[:, :1]), T)
inv = np.where(kind[:, None] == 2, r.uniform(1.0, 8.0, n)[:, None] * np.clip(1.0 - zm / h[:, None], 0, 1), 0.0)
T = T + inv
sup = (kind[:, None] == 0) & (zm < 60.0)
T = T + np.where(sup, r.uniform(0.3, 2.5, n)[:, None] * (1.0 - zm / 60.0), 0.0)
T = np.maximum(T, 190.0)
# pressure (hydrostatic, dry approximation)
pf = np.empty((n, nz + 1))
pf[:, 0] = ps
for k in range(nz):
    pf[:, k + 1] = pf[:, k] * np.exp(-g * dz[:, k] / (rd * T[:, k]))
p3 = 0.5 * (pf[:, 1:] + pf[:, :-1])
ex = (p3 / 1e5) ** (rd / cp)
th = T / ex
es = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
qsat = 0.622 * es / np.maximum(p3 - es, 1.0)
rh = np.clip(r.uniform(0.3, 0.97, n)[:, None] - 0.5 * zm / 12000.0, 0.02, 1.0)
qv = np.maximum(rh * qsat, 1e-7)
qv = np.where((kind[:, None] == 0) & in_bl, qv[:, :1] * np.ones((1, nz)), qv)
qv = np.minimum(qv, 1.02 * qsat)
cloudy = r.random(n) < 0.35
near_top = (np.abs(zm - h[:, None]) < 300.0) & cloudy[:, None]
qc = np.where(near_top & (T > 245.0), r.uniform(1e-6, 8e-4, (n, nz)), 0.0)
qi = np.where((T < 260.0) & (r.random((n, nz)) < 0.3), r.uniform(1e-8, 2e-4, (n, nz)), 0.0)
qs = np.where((T < 270.0) & (r.random((n, nz)) < 0.3), r.uniform(1e-8, 5e-4, (n, nz)), 0.0)
rho = p3 / (rd * T * (1.0 + 0.608 * qv))
ws = r.choice([0.5, 3.0, 8.0, 16.0, 30.0], n) * r.uniform(0.5, 1.2, n)
ang = r.uniform(0, 2 * np.pi, n)
prof = np.log1p(zm / 0.1) / np.log1p(zm[:, :1] / 0.1)
u = ws[:, None] * np.cos(ang)[:, None] * prof ** 0.5 + r.normal(0, 0.5, (n, nz))
v = ws[:, None] * np.sin(ang)[:, None] * prof ** 0.5 + r.normal(0, 0.5, (n, nz))
wamp = np.where(r.random(n) < 0.1, r.uniform(1.0, 3.0, n), r.uniform(0.0, 0.3, n))
w = wamp[:, None] * r.normal(0, 1, (n, nz))
rthraten = r.normal(0, 2e-5, (n, nz))
ts = T[:, 0] + np.where(kind == 0, r.uniform(1.0, 12.0, n), np.where(kind == 1, r.uniform(-0.5, 0.5, n), r.uniform(-6.0, -0.5, n)))
qsfc = qv[:, 0] * r.uniform(0.9, 1.2, n)
ust = np.where(kind == 2, r.uniform(0.02, 0.3, n), r.uniform(0.1, 1.0, n))
ch = r.uniform(0.001, 0.05, n)
hfx = np.where(kind == 0, r.uniform(20.0, 600.0, n), np.where(kind == 1, r.uniform(-10.0, 30.0, n), r.uniform(-80.0, -1.0, n)))
hfx = np.where(r.random(n) < 0.02, r.uniform(1200.0, 1500.0, n), hfx)
qfx = np.where(kind == 2, r.uniform(-2e-5, 2e-5, n), r.uniform(0.0, 3e-4, n))
qfx = np.where(r.random(n) < 0.02, r.uniform(9e-4, 1.2e-3, n), qfx)
wspd = np.sqrt(u[:, 0] ** 2 + v[:, 0] ** 2) + 0.1
znt = np.where(xland == 1.0, r.uniform(0.01, 1.0, n), 1e-4)
uoce = np.where(r.random(n) < 0.2, r.uniform(-0.5, 0.5, n), 0.0)
voce = np.where(r.random(n) < 0.2, r.uniform(-0.5, 0.5, n), 0.0)
with open(a.out, "wb") as f:
    np.array([n, nz, a.nsteps], dtype=np.int32).tofile(f)
    np.array([a.dt, a.dx], dtype=np.float64).tofile(f)
    for x in (dz, u, v, w, th, T, p3, ex, rho, qv, qc, qi, qs, rthraten):
        np.asfortranarray(x).T.astype(np.float64).tofile(f)  # (i, k) Fortran order: i fastest
    for x in (xland, ts, qsfc, ps, ust, ch, hfx, qfx, wspd, znt, uoce, voce):
        np.asarray(x, dtype=np.float64).tofile(f)
