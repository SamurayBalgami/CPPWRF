#!/usr/bin/env python3
"""Random surface-layer columns for the MYNN surface-layer differential test (record order: fields.txt).

  generate_inputs.py --cols N --itimestep T --nsteps S [--isfflx 1] [--seed 1] OUT.bin

Land (no snow / snow >= 5 cm), open water, lakes and partial sea ice; stable, neutral and unstable cases;
calm to strong winds; the INOUT fields as a model state of an earlier step would have them (zeros at T = 1).
"""
import argparse
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--cols", type=int, default=2500)
p.add_argument("--itimestep", type=int, default=1)
p.add_argument("--nsteps", type=int, default=1)
p.add_argument("--isfflx", type=int, default=1)
p.add_argument("--seed", type=int, default=1)
p.add_argument("out")
a = p.parse_args()
r = np.random.default_rng(a.seed)
n = a.cols
kind = r.integers(0, 5, n)  # 0 land, 1 land with snow, 2 ocean, 3 lake, 4 partial sea ice (water point)
xland = np.where(kind <= 1, 1.0, 2.0)
snowh = np.where(kind == 1, r.uniform(0.05, 1.0, n), np.where(kind == 0, r.uniform(0.0, 0.049, n) * (r.random(n) < 0.3), 0.0))
xice = np.where(kind == 4, r.uniform(0.5, 0.99, n), 0.0)
lakemask = np.where(kind == 3, 1.0, 0.0)
psfc = r.uniform(70000.0, 103000.0, n)
dz1 = r.uniform(15.0, 80.0, n)
low = r.random(n) < 0.15  # za <= 7 m (u10 from the log profile or level 2)
dz1[low] = r.uniform(5.0, 14.0, low.sum())
dz2 = dz1 * r.uniform(1.0, 1.3, n)
t1 = r.uniform(240.0, 310.0, n)
tsk = t1 + r.normal(0.0, 4.0, n)
tsk = np.where(kind == 2, np.maximum(tsk, 271.4), tsk)
tsk = np.where((kind == 1) | (kind == 4), np.minimum(tsk, 273.15), tsk)
p1 = psfc - 1.2 * 9.81 * dz1 * 0.5
p2 = p1 - 1.2 * 9.81 * (dz1 + dz2) * 0.5
t2 = t1 - 0.0065 * (dz1 + dz2) * 0.5
qv1 = r.uniform(0.2, 0.9, n) * 0.6112 * np.exp(17.67 * (t1 - 273.15) / (t1 - 29.65)) * 0.622 / (p1 / 1000.0)
qv2 = qv1 * r.uniform(0.85, 1.0, n)
ws = r.choice([0.2, 2.0, 6.0, 15.0, 30.0], n) * r.uniform(0.3, 1.2, n)
ang = r.uniform(0, 2 * np.pi, n)
u1, v1 = ws * np.cos(ang), ws * np.sin(ang)
u2, v2 = u1 * r.uniform(1.0, 1.4, n), v1 * r.uniform(1.0, 1.4, n)
th1 = t1 * (1e5 / p1) ** (287.0 / 1004.5)
th2 = t2 * (1e5 / p2) ** (287.0 / 1004.5)
rho1 = p1 / (287.0 * t1 * (1 + 0.608 * qv1))
rho2 = p2 / (287.0 * t2 * (1 + 0.608 * qv2))
mavail = np.where(kind <= 1, r.uniform(0.05, 1.0, n), 1.0)
dx = r.choice([1000.0, 3000.0, 9000.0, 27000.0], n)
water_depth = r.uniform(5.0, 4000.0, n)
in3 = {"u": (u1, u2), "v": (v1, v2), "t": (t1, t2), "qv": (qv1, qv2), "p": (p1, p2), "dz8w": (dz1, dz2),
       "th": (th1, th2), "rho": (rho1, rho2)}
in2 = [psfc, mavail, xland, tsk, snowh, dx, xice, lakemask, water_depth]
names = "pblh chs chs2 cqs cqs2 cpm znt ust ustm zol mol rmol psim psih qgh hfx qfx lh flhc flqc qsfc u10 v10 th2 t2 q2 gz1oz0 wspd br ch ck cka cd cda".split()
io = {k: np.zeros(n) for k in names}
io["znt"] = np.where(kind <= 1, r.choice([0.01, 0.05, 0.1, 0.5, 1.0], n), 1e-4)
io["pblh"] = r.uniform(50.0, 3000.0, n)
if a.itimestep > 1:  # a state of a previous step
    io["ust"] = r.uniform(0.01, 1.5, n)
    io["ustm"] = io["ust"] * r.uniform(0.9, 1.1, n)
    io["hfx"] = r.uniform(-80.0, 400.0, n)
    io["qfx"] = r.uniform(-2e-5, 2e-4, n)
    io["lh"] = io["qfx"] * 2.5e6
    io["mol"] = r.uniform(-1.0, 1.0, n)
    io["zol"] = r.uniform(-5.0, 5.0, n)
    io["qsfc"] = qv1 * r.uniform(0.9, 1.3, n)
    for k in ("chs", "chs2", "cqs", "cqs2", "flhc", "flqc", "rmol", "psim", "psih", "gz1oz0", "wspd", "br", "ch", "cpm"):
        io[k] = r.uniform(0.0, 1.0, n)
else:
    io["ust"] = np.where(r.random(n) < 0.5, 0.0, r.uniform(0.01, 0.5, n))
with open(a.out, "wb") as f:
    np.array([n, a.itimestep, a.isfflx, a.nsteps], dtype=np.int32).tofile(f)
    for k in ("u", "v", "t", "qv", "p", "dz8w", "th", "rho"):
        for lev in in3[k]:
            np.asarray(lev, dtype=np.float64).tofile(f)
    for x in in2:
        np.asarray(x, dtype=np.float64).tofile(f)
    for k in names:
        io[k].astype(np.float64).tofile(f)
