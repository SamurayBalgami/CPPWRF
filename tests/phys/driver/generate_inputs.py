#!/usr/bin/env python3
"""
Generates random but physical inputs for WRF physics driver glue differential testing.
Includes terrain, non-unit map factors, realistic atmospheric vertical profile, and
realistic physics tendencies.
"""

import argparse
import struct
import numpy as np

def generate(out_file, nx=16, ny=16, nz=35, seed=42):
    np.random.seed(seed)

    n_moist = 6
    p_top = 5000.0  # 50 hPa
    use_theta_m = 1

    # 1. Vertical coordinates (eta levels znw)
    eta = np.linspace(0.0, 1.0, nz + 1)
    znw = 1.0 - eta**1.5
    znw[0] = 1.0
    znw[-1] = 0.0

    dnw = znw[1:] - znw[:-1]  # length nz (<0)
    znu = 0.5 * (znw[1:] + znw[:-1])  # length nz

    dn = np.zeros(nz + 1, dtype=np.float64)
    dn[1:nz] = 0.5 * (dnw[1:nz] + dnw[:nz-1])
    dn[0] = 0.0
    dn[nz] = 0.0

    fnm = np.zeros(nz + 1, dtype=np.float64)
    fnp = np.zeros(nz + 1, dtype=np.float64)
    for k in range(1, nz):
        fnm[k] = 0.5 * dnw[k - 1] / dn[k]
        fnp[k] = 0.5 * dnw[k] / dn[k]

    c1h = np.ones(nz, dtype=np.float64)
    c2h = np.zeros(nz, dtype=np.float64)
    c1f = np.ones(nz + 1, dtype=np.float64)
    c2f = np.zeros(nz + 1, dtype=np.float64)

    # 2. Horizontal 2D fields with terrain and non-unit map factors
    x = np.arange(nx, dtype=np.float64)
    y = np.arange(ny, dtype=np.float64)
    X, Y = np.meshgrid(x, y)  # shape (ny, nx)

    # Terrain: Bell-shaped mountain up to 1200 m
    ht = 1200.0 / (1.0 + ((X - nx / 2.0) / 3.0)**2 + ((Y - ny / 2.0) / 3.0)**2)

    # Map factors varying across domain (periodic)
    msftx = 1.0 + 0.06 * np.sin(2.0 * np.pi * X / nx) * np.cos(2.0 * np.pi * Y / ny)
    msfty = 1.0 + 0.04 * np.cos(2.0 * np.pi * X / nx) * np.sin(2.0 * np.pi * Y / ny)

    # U points: shape (ny, nx + 1)
    Xu, Yu = np.meshgrid(np.arange(nx + 1, dtype=np.float64), y)
    msfux = 1.0 + 0.06 * np.sin(2.0 * np.pi * Xu / nx) * np.cos(2.0 * np.pi * Yu / ny)
    msfuy = 1.0 + 0.04 * np.cos(2.0 * np.pi * Xu / nx) * np.sin(2.0 * np.pi * Yu / ny)

    # V points: shape (ny + 1, nx)
    Xv, Yv = np.meshgrid(x, np.arange(ny + 1, dtype=np.float64))
    msfvx = 1.0 + 0.06 * np.sin(2.0 * np.pi * Xv / nx) * np.cos(2.0 * np.pi * Yv / ny)
    msfvy = 1.0 + 0.04 * np.cos(2.0 * np.pi * Xv / nx) * np.sin(2.0 * np.pi * Yv / ny)
    msfvx_inv = 1.0 / msfvx

    # Surface base pressure from standard atmosphere
    p_sfc = 101325.0 * (1.0 - 0.0065 * ht / 288.15)**5.255
    mut = (p_sfc - p_top) + np.random.uniform(-300.0, 300.0, size=(ny, nx))

    # muu (ny, nx + 1) and muv (ny + 1, nx) with periodic wrapping
    muu = np.zeros((ny, nx + 1), dtype=np.float64)
    for i in range(nx + 1):
        im = (i - 1) % nx
        ip = i % nx
        muu[:, i] = 0.5 * (mut[:, im] + mut[:, ip])

    muv = np.zeros((ny + 1, nx), dtype=np.float64)
    for j in range(ny + 1):
        jm = (j - 1) % ny
        jp = j % ny
        muv[j, :] = 0.5 * (mut[jm, :] + mut[jp, :])

    # 3. 3D state fields
    g = 9.81
    z_nom = np.zeros((nz + 1, ny, nx), dtype=np.float64)
    for k in range(nz + 1):
        z_nom[k] = ht + (1.0 - znw[k]) * 16000.0
    phb = g * z_nom
    ph = np.random.uniform(-50.0, 50.0, size=(nz + 1, ny, nx))

    pb = np.zeros((nz, ny, nx), dtype=np.float64)
    for k in range(nz):
        pb[k] = znu[k] * mut + p_top
    p = np.random.uniform(-250.0, 250.0, size=(nz, ny, nx))

    t = np.zeros((nz, ny, nx), dtype=np.float64)
    for k in range(nz):
        th = 295.0 + (1.0 - znu[k]) * 125.0 + np.random.uniform(-2.0, 2.0, size=(ny, nx))
        t[k] = th - 300.0

    moist = np.zeros((n_moist, nz, ny, nx), dtype=np.float64)
    for k in range(nz):
        qv_scale = 0.014 * np.exp(-z_nom[k] / 3000.0)
        moist[0, k] = np.clip(qv_scale * (1.0 + np.random.uniform(-0.15, 0.15, size=(ny, nx))), 1e-6, 0.02)
        if 5 <= k <= 18:
            moist[1, k] = np.random.uniform(0.0, 1.5e-3, size=(ny, nx))
        if k <= 10:
            moist[2, k] = np.random.uniform(0.0, 2.0e-3, size=(ny, nx))
        if k >= 15:
            moist[3, k] = np.random.uniform(0.0, 4.0e-4, size=(ny, nx))
        if 12 <= k <= 25:
            moist[4, k] = np.random.uniform(0.0, 8.0e-4, size=(ny, nx))
        if 8 <= k <= 22:
            moist[5, k] = np.random.uniform(0.0, 6.0e-4, size=(ny, nx))

    alt = np.zeros((nz, ny, nx), dtype=np.float64)
    rd = 287.0
    rcp = 287.0 / 1004.5
    for k in range(nz):
        p_full = pb[k] + p[k]
        pi_full = (p_full / 100000.0)**rcp
        t_dry = (t[k] + 300.0) / (1.0 + (461.6 / 287.0) * moist[0, k])
        t_abs = t_dry * pi_full
        rho = p_full / (rd * t_abs * (1.0 + 0.608 * moist[0, k]))
        alt[k] = (1.0 + moist[0, k]) / rho

    # Velocities with periodic boundary face matching u[:, :, nx] = u[:, :, 0]
    u = np.random.uniform(-15.0, 25.0, size=(nz, ny, nx + 1))
    u[:, :, nx] = u[:, :, 0]

    v = np.random.uniform(-15.0, 15.0, size=(nz, ny + 1, nx))
    v[:, ny, :] = v[:, 0, :]

    # 4. Physics tendencies
    rthraten = np.random.uniform(-4.0e-5, 1.5e-5, size=(nz, ny, nx))
    rthblten = np.random.uniform(-6.0e-4, 6.0e-4, size=(nz, ny, nx))
    rublten = np.random.uniform(-3.0e-3, 3.0e-3, size=(nz, ny, nx))
    rvblten = np.random.uniform(-3.0e-3, 3.0e-3, size=(nz, ny, nx))
    rqvblten = np.random.uniform(-1.5e-6, 1.5e-6, size=(nz, ny, nx))
    rqcblten = np.random.uniform(-5.0e-7, 5.0e-7, size=(nz, ny, nx))
    rqiblten = np.random.uniform(-1.0e-7, 1.0e-7, size=(nz, ny, nx))

    rthcuten = np.random.uniform(-1.0e-3, 2.0e-3, size=(nz, ny, nx))
    rucuten = np.random.uniform(-1.0e-3, 1.0e-3, size=(nz, ny, nx))
    rvcuten = np.random.uniform(-1.0e-3, 1.0e-3, size=(nz, ny, nx))
    rqvcuten = np.random.uniform(-2.0e-6, 2.0e-6, size=(nz, ny, nx))
    rqccuten = np.random.uniform(-1.0e-6, 1.0e-6, size=(nz, ny, nx))
    rqrcuten = np.random.uniform(-1.0e-6, 1.0e-6, size=(nz, ny, nx))
    rqicuten = np.random.uniform(-5.0e-7, 5.0e-7, size=(nz, ny, nx))
    rqscuten = np.random.uniform(-5.0e-7, 5.0e-7, size=(nz, ny, nx))

    rthshten = np.random.uniform(-3.0e-4, 3.0e-4, size=(nz, ny, nx))
    rushten = np.random.uniform(-5.0e-4, 5.0e-4, size=(nz, ny, nx))
    rvshten = np.random.uniform(-5.0e-4, 5.0e-4, size=(nz, ny, nx))
    rqvshten = np.random.uniform(-5.0e-7, 5.0e-7, size=(nz, ny, nx))
    rqcshten = np.random.uniform(-2.0e-7, 2.0e-7, size=(nz, ny, nx))
    rqrshten = np.random.uniform(-2.0e-7, 2.0e-7, size=(nz, ny, nx))
    rqishten = np.random.uniform(-1.0e-7, 1.0e-7, size=(nz, ny, nx))
    rqsshten = np.random.uniform(-1.0e-7, 1.0e-7, size=(nz, ny, nx))
    rqgshten = np.random.uniform(-1.0e-7, 1.0e-7, size=(nz, ny, nx))

    # Write binary file
    with open(out_file, 'wb') as f:
        # Header
        f.write(struct.pack('4i', nx, ny, nz, n_moist))
        f.write(struct.pack('d', p_top))
        f.write(struct.pack('i', use_theta_m))

        # 1D vertical arrays
        znw.astype(np.float64).tofile(f)
        znu.astype(np.float64).tofile(f)
        dnw.astype(np.float64).tofile(f)
        dn.astype(np.float64).tofile(f)
        fnm.astype(np.float64).tofile(f)
        fnp.astype(np.float64).tofile(f)
        c1h.astype(np.float64).tofile(f)
        c2h.astype(np.float64).tofile(f)
        c1f.astype(np.float64).tofile(f)
        c2f.astype(np.float64).tofile(f)

        # 2D arrays
        mut.astype(np.float64).tofile(f)
        muu.astype(np.float64).tofile(f)
        muv.astype(np.float64).tofile(f)
        msftx.astype(np.float64).tofile(f)
        msfty.astype(np.float64).tofile(f)
        msfux.astype(np.float64).tofile(f)
        msfuy.astype(np.float64).tofile(f)
        msfvx.astype(np.float64).tofile(f)
        msfvy.astype(np.float64).tofile(f)
        msfvx_inv.astype(np.float64).tofile(f)
        ht.astype(np.float64).tofile(f)

        # 3D arrays
        u.astype(np.float64).tofile(f)
        v.astype(np.float64).tofile(f)
        t.astype(np.float64).tofile(f)
        p.astype(np.float64).tofile(f)
        pb.astype(np.float64).tofile(f)
        alt.astype(np.float64).tofile(f)
        ph.astype(np.float64).tofile(f)
        phb.astype(np.float64).tofile(f)

        for m in range(n_moist):
            moist[m].astype(np.float64).tofile(f)

        # Tendencies
        rthraten.astype(np.float64).tofile(f)
        rthblten.astype(np.float64).tofile(f)
        rublten.astype(np.float64).tofile(f)
        rvblten.astype(np.float64).tofile(f)
        rqvblten.astype(np.float64).tofile(f)
        rqcblten.astype(np.float64).tofile(f)
        rqiblten.astype(np.float64).tofile(f)

        rthcuten.astype(np.float64).tofile(f)
        rucuten.astype(np.float64).tofile(f)
        rvcuten.astype(np.float64).tofile(f)
        rqvcuten.astype(np.float64).tofile(f)
        rqccuten.astype(np.float64).tofile(f)
        rqrcuten.astype(np.float64).tofile(f)
        rqicuten.astype(np.float64).tofile(f)
        rqscuten.astype(np.float64).tofile(f)

        rthshten.astype(np.float64).tofile(f)
        rushten.astype(np.float64).tofile(f)
        rvshten.astype(np.float64).tofile(f)
        rqvshten.astype(np.float64).tofile(f)
        rqcshten.astype(np.float64).tofile(f)
        rqrshten.astype(np.float64).tofile(f)
        rqishten.astype(np.float64).tofile(f)
        rqsshten.astype(np.float64).tofile(f)
        rqgshten.astype(np.float64).tofile(f)

    print(f"Generated input binary: {out_file} (nx={nx}, ny={ny}, nz={nz}, n_moist={n_moist})")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', default='driver_input.bin')
    parser.add_argument('--nx', type=int, default=16)
    parser.add_argument('--ny', type=int, default=16)
    parser.add_argument('--nz', type=int, default=35)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()
    generate(args.out, nx=args.nx, ny=args.ny, nz=args.nz, seed=args.seed)
