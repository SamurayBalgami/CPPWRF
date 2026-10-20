#!/usr/bin/env python3
"""
Column inputs for the Noah ice differential test (SFLX_GLACIAL and SFLX_SEAICE).
Regimes: cold and dry, snowfall, freezing rain, rain, melting snowpack (T1 above freezing,
strong solar), snowpack below the schemes' lower bounds (reset), thin snowpack under strong
forcing (melts completely within a long step), day and night forcing; a fifth of the columns
have a rough snow-free surface (Z0BRD up to 0.05 m, partly buried snow in SNOWZ0).
Every column is run through both routines; WRF stops when SNEQV/SNOWH > 1, so the input
snow density stays in 0.1..0.5.

Layout (native endian): int32 ncols, nsoil; float64 dt, lvcoef; float64 sldpth[nsoil];
float64 [ncols] for each name in IN_VARS; float64 stc[nsoil][ncols].
"""

import argparse
import struct
import numpy as np

IN_VARS = ["ffrozp", "lwdn", "solnet", "sfcprs", "prcp", "sfctmp", "q2", "th2", "q2sat", "dqsdt2",
           "alb", "snoalb", "tbot", "z0brd", "emissi", "embrd", "t1", "snowh", "sneqv", "ch",
           "sncovr", "snotime1", "ribb", "sithick"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cols", type=int, default=2500)
    ap.add_argument("--dt", type=float, default=60.0)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", default="noah_ice_input.bin")
    a = ap.parse_args()

    rng = np.random.default_rng(a.seed)
    n, nsoil = a.cols, 4
    sldpth = np.array([0.1, 0.3, 0.6, 1.0])
    lvcoef = 0.5  # GENPARM.TBL LVCOEF_DATA

    regime = rng.integers(0, 7, size=n)  # 0 cold dry, 1 snowfall, 2 freezing rain, 3 rain, 4 melt, 5 thin snow, 6 melt-out
    u = lambda lo, hi: rng.uniform(lo, hi, size=n)

    sfctmp = np.where(regime == 0, u(235.0, 265.0), 0.0)
    sfctmp = np.where(regime == 1, u(245.0, 273.0), sfctmp)
    sfctmp = np.where(regime == 2, u(262.0, 274.0), sfctmp)
    sfctmp = np.where(regime == 3, u(273.5, 280.0), sfctmp)
    sfctmp = np.where(regime == 4, u(273.5, 286.0), sfctmp)
    sfctmp = np.where(regime == 5, u(255.0, 278.0), sfctmp)
    sfctmp = np.where(regime == 6, u(285.0, 300.0), sfctmp)

    t1 = sfctmp + u(-6.0, 3.0)
    t1 = np.where(regime == 2, np.minimum(t1, 273.10), t1)
    t1 = np.where(regime == 4, u(273.2, 278.0), t1)
    t1 = np.where(regime == 6, u(274.0, 280.0), t1)

    prcp = np.where(np.isin(regime, [1, 2, 3]), u(1.0e-5, 2.0e-3), 0.0)
    ffrozp = np.where(regime == 1, 1.0, 0.0)
    ffrozp = np.where(regime == 5, rng.integers(0, 2, size=n).astype(float), ffrozp)
    prcp = np.where(regime == 5, np.where(rng.random(n) < 0.5, u(1.0e-5, 1.0e-3), 0.0), prcp)

    day = rng.random(n) < 0.6
    solnet = np.where(day, u(10.0, 350.0), 0.0)
    solnet = np.where(regime == 4, u(200.0, 700.0), solnet)
    lwdn = u(150.0, 380.0)
    solnet = np.where(regime == 6, u(900.0, 1150.0), solnet)
    lwdn = np.where(regime == 6, u(380.0, 460.0), lwdn)

    sfcprs = u(60000.0, 102000.0)
    psfc = sfcprs + u(30.0, 120.0)
    th2 = sfctmp * (psfc / sfcprs) ** 0.2857

    a23m4, a4 = 17.67 * (273.15 - 29.65), 29.65
    es = 611.2 * np.exp(17.67 * (sfctmp - 273.15) / (sfctmp - a4))
    qs = 0.622 * es / (sfcprs - 0.378 * es)
    q2sat = qs / (1.0 + qs)
    q2 = q2sat * u(0.3, 1.05)
    dqsdt2 = q2sat * a23m4 / ((sfctmp - a4) ** 2)

    alb = u(0.55, 0.80)
    snoalb = u(0.60, 0.85)
    tbot = u(250.0, 271.4)
    z0brd = np.where(rng.random(n) < 0.2, u(0.015, 0.05), u(0.001, 0.01))
    emissi = u(0.95, 0.99)
    embrd = u(0.95, 0.99)
    ch = u(1.0e-4, 0.03)
    sncovr = u(0.5, 1.0)
    snotime1 = np.where(rng.random(n) < 0.3, 0.0, u(0.0, 2.0e5))
    ribb = u(-1.0, 1.0)
    sithick = u(0.5, 3.0)

    sneqv = u(0.12, 1.2)
    sneqv = np.where(regime == 4, u(0.005, 0.3), sneqv)
    sneqv = np.where(regime == 5, u(0.0, 0.12), sneqv)
    sneqv = np.where(regime == 6, u(0.0, 0.012), sneqv)
    snowh = sneqv / u(0.1, 0.5)
    snowh = np.where(sneqv == 0.0, 0.0, snowh)

    stc = np.empty((nsoil, n))
    top = np.minimum(t1, 273.1) + u(-4.0, 0.0)
    for k in range(nsoil):
        w = (k + 1) / (nsoil + 1)
        stc[k] = np.minimum((1.0 - w) * top + w * tbot + u(-0.5, 0.5), 273.14)

    vals = dict(ffrozp=ffrozp, lwdn=lwdn, solnet=solnet, sfcprs=sfcprs, prcp=prcp, sfctmp=sfctmp, q2=q2,
                th2=th2, q2sat=q2sat, dqsdt2=dqsdt2, alb=alb, snoalb=snoalb, tbot=tbot, z0brd=z0brd,
                emissi=emissi, embrd=embrd, t1=t1, snowh=snowh, sneqv=sneqv, ch=ch, sncovr=sncovr,
                snotime1=snotime1, ribb=ribb, sithick=sithick)
    with open(a.out, "wb") as f:
        f.write(struct.pack("=ii", n, nsoil))
        f.write(struct.pack("=dd", a.dt, lvcoef))
        f.write(sldpth.astype("=f8").tobytes())
        for name in IN_VARS:
            f.write(np.asarray(vals[name], dtype="=f8").tobytes())
        f.write(stc.astype("=f8").tobytes())
    print(f"wrote {a.out}: {n} columns, dt = {a.dt}, regimes {np.bincount(regime, minlength=7).tolist()}")


if __name__ == "__main__":
    main()
