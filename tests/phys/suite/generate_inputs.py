#!/usr/bin/env python3
"""Synthetic but physical real-data tile for the P1 physics-suite differential test.

Writes suite_input.bin: the fields a wrfinput_d01 would provide for a 24x20x44 domain with
terrain, open water, sea ice, land ice (glacier), snow cover, stable and unstable columns and
realistic soil states. The Fortran reference (fortran/suite_ref.F90) fills grid% with these
fields and calls WRF's own start_domain (phy_init included), so the state both codes start
from is the one WRF itself produces.

Record format (stream, little endian): 32-byte name, int32 nd, int32 n1,n2,n3, float64 data
in Fortran order (first index fastest). 3D fields are (i,k,j) as in WRF memory order.
"""
import struct
import sys

import numpy as np

NX, NY, NZ, NSOIL = 24, 20, 44, 4  # mass points
P_TOP = 5000.0
G, RD, RV, CP = 9.81, 287.0, 461.6, 7.0 * 287.0 / 2.0
P1000, T0 = 1.0e5, 300.0
P00, T00, TLP, TISO = 1.0e5, 290.0, 50.0, 200.0
ISWATER, ISICE, ISURBAN = 17, 15, 13


def eta_levels(n):
    # WRF-like stretched eta levels (dense near the surface)
    x = np.linspace(0.0, 1.0, n)
    znw = 1.0 - (0.62 * x + 0.38 * x ** 2.2)
    znw[0], znw[-1] = 1.0, 0.0
    return znw


def main(out):
    rng = np.random.RandomState(20261013)
    znw = eta_levels(NZ + 1)
    dnw = znw[1:] - znw[:-1]
    znu = 0.5 * (znw[1:] + znw[:-1])
    dn = np.zeros(NZ)
    fnm = np.zeros(NZ)
    fnp = np.zeros(NZ)
    rdn = np.zeros(NZ)
    for k in range(1, NZ):
        dn[k] = 0.5 * (dnw[k] + dnw[k - 1])
        rdn[k] = 1.0 / dn[k]
        fnp[k] = 0.5 * dnw[k] / dn[k]
        fnm[k] = 0.5 * dnw[k - 1] / dn[k]

    i = np.arange(NX)[:, None]
    j = np.arange(NY)[None, :]
    # ---- surface types (i,j) ----
    water = (i < 6) | ((i < 9) & (j < 4))
    seaice = water & (j >= 14)              # north-western sea is ice covered
    # terrain: a ridge rising to the east, with a peak carrying a glacier
    ht = np.where(water, 0.0, 1900.0 * np.exp(-(((i - 18) / 4.5) ** 2 + ((j - 11) / 5.0) ** 2))
                  + 25.0 * (i - 6).clip(0))
    glacier = (~water) & (ht > 1700.0)
    xland = np.where(water & ~seaice, 2.0, 1.0)
    landmask = np.where(water & ~seaice, 0.0, 1.0)
    xice = np.where(seaice, 1.0, 0.0)
    cats = np.array([1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 19, 20])
    ivgtyp = cats[rng.randint(0, len(cats), size=(NX, NY))]
    ivgtyp = np.where(water, ISWATER, ivgtyp)
    ivgtyp = np.where(seaice, ISICE, ivgtyp)
    ivgtyp = np.where(glacier, ISICE, ivgtyp)
    soils = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 17, 18, 19])
    isltyp = soils[rng.randint(0, len(soils), size=(NX, NY))]
    isltyp = np.where(water, 14, isltyp)
    isltyp = np.where(seaice | glacier, 16, isltyp)
    isltyp[10, 3] = 14                       # land point with water soil (SFLX reset to 7)
    vegfra = np.where(water | seaice | glacier, 0.0, rng.uniform(5.0, 95.0, size=(NX, NY)))
    vegfra = np.where(ivgtyp == 16, 1.0, vegfra)
    shdmax = np.where(vegfra > 0, np.minimum(vegfra + rng.uniform(0, 20, size=(NX, NY)), 99.0), 0.0)
    shdmin = np.where(vegfra > 0, np.maximum(vegfra - rng.uniform(0, 30, size=(NX, NY)), 0.0), 0.0)
    albbck = np.where(water, 0.08, rng.uniform(0.12, 0.30, size=(NX, NY)))
    snoalb = np.where(water & ~seaice, 0.08, rng.uniform(0.45, 0.75, size=(NX, NY)))
    # ---- skin / soil temperatures: warm sea, warm lowlands (unstable), cold snowy ridge ----
    tsk = 296.0 + 4.0 * np.sin(0.4 * i) * np.cos(0.3 * j) - 0.0065 * ht
    tsk = np.where(water, 291.5 + 0.1 * j, tsk)
    tsk = np.where(seaice, 262.0 + 0.5 * i, tsk)
    tsk = np.where(glacier, 258.0, tsk)
    sst = np.where(water, 291.5 + 0.1 * j, 0.0)
    sst = np.where(seaice, 271.4, sst)
    tmn = np.where(water & ~seaice, sst, 286.0 - 0.0065 * ht)
    tmn = np.where(seaice, 271.4, tmn)
    snowy = (~water) & ((ht > 1100.0) | ((j > 15) & (i > 12)))
    snow = np.where(snowy, rng.uniform(5.0, 120.0, size=(NX, NY)), 0.0)  # SWE kg/m2
    snow = np.where(glacier, 400.0, snow)
    snow = np.where(seaice, 20.0, snow)
    snowh = np.where(snow > 0, snow * 0.001 * rng.uniform(4.0, 8.0, size=(NX, NY)), 0.0)
    tsk = np.where(snowy, np.minimum(tsk, 272.0), tsk)
    canwat = np.where(water | seaice | glacier, 0.0, rng.uniform(0.0, 0.4, size=(NX, NY)))
    zs = np.array([0.05, 0.25, 0.70, 1.50])
    dzs = np.array([0.10, 0.30, 0.60, 1.00])
    tslb = np.zeros((NX, NSOIL, NY))
    smois = np.zeros((NX, NSOIL, NY))
    sh2o = np.zeros((NX, NSOIL, NY))
    for n in range(NSOIL):
        w = zs[n] / 1.5
        tslb[:, n, :] = (1 - w) * tsk + w * tmn
        smois[:, n, :] = rng.uniform(0.12, 0.42, size=(NX, NY))
    smois = np.where(water[:, None, :] | seaice[:, None, :] | glacier[:, None, :], 1.0, smois)
    frozen = tslb < 273.15
    sh2o = np.where(frozen, 0.6 * smois, smois)
    sh2o = np.where(water[:, None, :], 1.0, sh2o)
    tslb = np.where((water & ~seaice)[:, None, :], sst[:, None, :], tslb)

    # ---- base state (start_em, use_baseparam_fr_nml) and a moist atmosphere ----
    p_surf = P00 * np.exp(-T00 / TLP + np.sqrt((T00 / TLP) ** 2 - 2.0 * G * ht / TLP / RD))
    mub = p_surf - P_TOP
    mu2 = 120.0 * np.sin(0.3 * i + 0.2 * j)                          # dry column-mass perturbation
    mu = mub + mu2
    pw = znw[None, :, None] * mu[:, None, :] + P_TOP              # dry pressure at w levels
    ph = znu[None, :, None] * mu[:, None, :] + P_TOP
    # temperature: surface-based mixed/stable layers depending on tsk, tropopause aloft
    z_est = 7000.0 * np.log(p_surf[:, None, :] / ph)
    tair = tsk[:, None, :] - 0.0065 * z_est
    inv = (snowy | seaice | glacier)[:, None, :] & (z_est < 400.0)     # surface inversions
    tair = np.where(inv, tsk[:, None, :] + 6.0 * z_est / 400.0, tair)
    tair = np.maximum(tair, 213.0)
    theta = tair * (P1000 / ph) ** (RD / CP)
    rh = 0.85 - 0.6 * (1.0 - ph / p_surf[:, None, :]).clip(0, 1)
    es = 611.2 * np.exp(17.67 * (tair - 273.15) / (tair - 29.65))
    qv = (0.622 * rh * es / (ph - rh * es)).clip(1e-6, 0.02)
    qv = np.where(ph < 20000.0, 3e-6, qv)
    t2 = theta * (1.0 + RV / RD * qv) - T0                         # moist theta - t0 (use_theta_m=1)
    # geopotential: hypsometric integration of the full (moist) column
    phb_full = np.zeros((NX, NZ + 1, NY))
    phb_full[:, 0, :] = ht * G
    for k in range(NZ):
        tv = tair[:, k, :] * (1.0 + 0.608 * qv[:, k, :])
        phb_full[:, k + 1, :] = phb_full[:, k, :] + RD * tv * np.log(pw[:, k, :] / pw[:, k + 1, :])
    # base geopotential as start_em computes it (hypsometric_opt = 2)
    phb = np.zeros((NX, NZ + 1, NY))
    phb[:, 0, :] = ht * G
    for k in range(1, NZ + 1):
        pb_h = znu[k - 1] * mub + P_TOP
        temp = np.maximum(TISO, T00 + TLP * np.log(pb_h / P00))
        alb = (RD / P1000) * temp * (P00 / pb_h) ** (RD / CP) * (pb_h / P1000) ** (-1.0 / 1.4)
        pfu = znw[k] * mub + P_TOP
        pfd = znw[k - 1] * mub + P_TOP
        phb[:, k, :] = phb[:, k - 1, :] + alb * pb_h * np.log(pfd / pfu)
    ph2 = phb_full - phb
    # winds: veering profile, strong near the coast, weak (stable) over snow
    zc = 7000.0 * np.log(P00 / ph)
    ubase = 6.0 + 12.0 * (1 - np.exp(-zc / 3000.0))
    u = np.zeros((NX + 1, NZ, NY))
    v = np.zeros((NX, NZ, NY + 1))
    for k in range(NZ):
        uu = ubase[:, k, :] * (1.0 + 0.3 * np.sin(0.5 * i + 0.3 * j))
        vv = 0.5 * ubase[:, k, :] * np.cos(0.25 * i - 0.4 * j) + 2.0
        u[:NX, k, :] = uu
        u[NX, k, :] = uu[-1, :]
        v[:, k, :NY] = vv
        v[:, k, NY] = vv[:, -1]
    u[:, 0, :] *= 0.8
    v[:, 0, :] *= 0.8
    w = np.zeros((NX, NZ + 1, NY))
    w[:, 1:NZ, :] = 0.05 * np.sin(0.2 * i)[:, None, :] * np.ones((1, NZ - 1, 1))

    xlat = 38.0 + 0.135 * j + 0.0 * i
    xlong = 32.0 + 0.17 * i + 0.0 * j
    # radiation forcing kept fixed by the drivers (ra_*_physics = 0): daytime values
    glw = 300.0 + 30.0 * np.cos(0.2 * i) + 0.0 * j
    swdown = 550.0 + 150.0 * np.sin(0.3 * j) + 0.0 * i
    gsw = swdown * (1.0 - albbck)

    recs = []

    def add(name, a):
        a = np.asarray(a, dtype=np.float64)
        dims = list(a.shape) + [1] * (3 - a.ndim)
        recs.append((name, a.ndim, dims, a))

    for n, a in [("ZNW", znw), ("ZNU", znu), ("DNW", dnw), ("RDNW", 1.0 / dnw), ("DN", dn),
                 ("RDN", rdn), ("FNM", fnm), ("FNP", fnp), ("ZS", zs), ("DZS", dzs)]:
        add(n, a)
    for n, a in [("HT", ht), ("MUB", mub), ("MU_2", mu2), ("XLAND", xland), ("LANDMASK", landmask),
                 ("XICE", xice), ("IVGTYP", ivgtyp), ("ISLTYP", isltyp), ("LU_INDEX", ivgtyp),
                 ("VEGFRA", vegfra), ("SHDMAX", shdmax), ("SHDMIN", shdmin), ("ALBBCK", albbck),
                 ("SNOALB", snoalb), ("TMN", tmn), ("TSK", tsk), ("SST", sst), ("SNOW", snow),
                 ("SNOWH", snowh), ("CANWAT", canwat), ("XLAT", xlat), ("XLONG", xlong),
                 ("GLW", glw), ("SWDOWN", swdown), ("GSW", gsw)]:
        add(n, a)
    for n, a in [("U", u), ("V", v), ("W", w), ("T", t2), ("QVAPOR", qv), ("PH", ph2),
                 ("TSLB", tslb), ("SMOIS", smois), ("SH2O", sh2o)]:
        add(n, a)

    with open(out, "wb") as f:
        for name, nd, dims, a in recs:
            f.write(name.ljust(32).encode())
            f.write(struct.pack("<4i", nd, *dims))
            f.write(np.asfortranarray(a).tobytes(order="F"))
    print(f"wrote {out}: {len(recs)} fields, {NX}x{NY}x{NZ}, "
          f"water {int((water & ~seaice).sum())}, sea ice {int(seaice.sum())}, "
          f"glacier {int(glacier.sum())}, snow {int((snow > 0).sum())}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "suite_input.bin")
