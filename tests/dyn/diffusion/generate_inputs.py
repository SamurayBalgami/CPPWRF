#!/usr/bin/env python3
"""
generate_inputs.py
Generates physical, realistic 3D input fields for differential testing of
WRF v4.6.1 Diffusion and Damping kernels against CPPWRF.
"""

import sys
import os
import struct
import numpy as np

def generate_inputs(output_path, nx=16, ny=16, nz=15):
    np.random.seed(42)

    dx = 1000.0  # meters
    dy = 1000.0  # meters
    dt = 6.0     # seconds
    g = 9.81
    r_d = 287.0
    cp = 1004.5
    p1000mb = 100000.0

    # 1. Vertical coordinates (eta from 1 to 0)
    # znw: w levels (nz + 1)
    znw = np.linspace(1.0, 0.0, nz + 1)
    # znu: mass levels (nz)
    znu = 0.5 * (znw[:-1] + znw[1:])
    
    dn = np.zeros(nz)
    rdn = np.zeros(nz)
    for k in range(nz):
        dn[k] = znw[k + 1] - znw[k]
        rdn[k] = 1.0 / dn[k]

    dnw = np.zeros(nz)
    rdnw = np.zeros(nz)
    for k in range(1, nz):
        dnw[k] = znu[k] - znu[k - 1]
        rdnw[k] = 1.0 / dnw[k]
    dnw[0] = dnw[1]
    rdnw[0] = rdnw[1]

    fnm = np.zeros(nz)
    fnp = np.zeros(nz)
    for k in range(1, nz):
        fnm[k] = (znw[k] - znu[k - 1]) / (znu[k] - znu[k - 1])
        fnp[k] = (znu[k] - znw[k]) / (znu[k] - znu[k - 1])
    fnm[0] = 0.0
    fnp[0] = 0.0

    cf1 = 1.5
    cf2 = -0.5
    cf3 = 0.0

    # Hybrid coordinate coefficients
    # c1 = eta, c2 = (1 - c1)*(1 - eta)
    c1 = np.zeros(nz)
    c2 = np.zeros(nz)
    c3 = np.zeros(nz)
    c4 = np.zeros(nz)
    for k in range(nz):
        c1[k] = znu[k]**2
        c2[k] = 1.0 - c1[k]
        c3[k] = c1[k]
        c4[k] = c2[k]

    c1f = np.zeros(nz + 1)
    c2f = np.zeros(nz + 1)
    for k in range(nz + 1):
        c1f[k] = znw[k]**2
        c2f[k] = 1.0 - c1f[k]

    # 1D base state soundings
    u_base = 10.0 + 15.0 * (1.0 - znu)
    v_base = 2.0 + 5.0 * (1.0 - znu)
    t_base = 300.0 - 50.0 * (1.0 - znu)
    qv_base = 0.012 * np.exp(- (1.0 - znu) * 3.0)

    # 2. 2D horizontal fields
    # Meshgrid x, y
    x = np.arange(nx) * dx
    y = np.arange(ny) * dy
    X, Y = np.meshgrid(x, y, indexing='ij')

    # Topography: Smooth bell hill centered in domain, vanishing at boundaries
    h_topo = 300.0 * (1.0 - np.cos(2 * np.pi * X / (nx * dx))) * (1.0 - np.cos(2 * np.pi * Y / (ny * dy))) / 4.0

    # Map scale factors (non-unit, periodic across domain, anisotropic msfx != msfy)
    msftx = 1.02 + 0.04 * np.cos(2 * np.pi * Y / (ny * dy))
    msfty = 1.01 + 0.03 * np.cos(2 * np.pi * X / (nx * dx))

    # Staggered map factors:
    # msfux at u points (nx + 1, ny)
    xu = (np.arange(nx + 1) - 0.5) * dx
    Xu, Yu = np.meshgrid(xu, y, indexing='ij')
    msfux = 1.02 + 0.04 * np.cos(2 * np.pi * Yu / (ny * dy))
    msfuy = 1.01 + 0.03 * np.cos(2 * np.pi * Xu / (nx * dx))

    # msfvx at v points (nx, ny + 1)
    yv = (np.arange(ny + 1) - 0.5) * dy
    Xv, Yv = np.meshgrid(x, yv, indexing='ij')
    msfvx = 1.02 + 0.04 * np.cos(2 * np.pi * Yv / (ny * dy))
    msfvy = 1.01 + 0.03 * np.cos(2 * np.pi * Xv / (nx * dx))

    # Column dry mass mut, muu, muv
    p_sfc = 101325.0 - 100.0 * (h_topo / 8.0)
    p_top = 10000.0
    mut = p_sfc - p_top

    # muu (nx + 1, ny)
    muu = np.zeros((nx + 1, ny))
    muu[1:nx, :] = 0.5 * (mut[:-1, :] + mut[1:, :])
    muu[0, :] = mut[0, :]
    muu[nx, :] = mut[-1, :]

    # muv (nx, ny + 1)
    muv = np.zeros((nx, ny + 1))
    muv[:, 1:ny] = 0.5 * (mut[:, :-1] + mut[:, 1:])
    muv[:, 0] = mut[:, 0]
    muv[:, ny] = mut[:, -1]

    # 3. 3D fields
    # Memory order: Fortran (i, k, j) where i is fastest, k middle, j slowest
    # phb: base geopotential (nx, nz + 1, ny)
    phb = np.zeros((nx, nz + 1, ny))
    for j in range(ny):
        for i in range(nx):
            z_sfc = h_topo[i, j]
            for k in range(nz + 1):
                # Standard atmosphere height scaling
                z_lev = z_sfc + (1.0 - znw[k]) * 12000.0
                phb[i, k, j] = g * z_lev

    # ph: perturbation geopotential (nx, nz + 1, ny)
    ph = np.zeros((nx, nz + 1, ny))
    for j in range(ny):
        for k in range(nz + 1):
            for i in range(nx):
                ph[i, k, j] = 50.0 * np.sin(2 * np.pi * i / nx) * np.cos(2 * np.pi * j / ny) * (1.0 - znw[k])

    # u: staggered in x (nx + 1, nz, ny) with strong vertical shear
    u = np.zeros((nx + 1, nz, ny))
    for j in range(ny):
        for k in range(nz):
            for i in range(nx + 1):
                shear_u = 8.0 * np.sin(np.pi * znu[k]) * (1.0 + 0.4 * np.cos(2 * np.pi * znu[k]))
                u[i, k, j] = u_base[k] + shear_u * (1.0 + 0.5 * np.sin(2 * np.pi * j / ny) * np.cos(2 * np.pi * i / nx))

    # v: staggered in y (nx, nz, ny + 1) with strong vertical shear
    v = np.zeros((nx, nz, ny + 1))
    for j in range(ny + 1):
        for k in range(nz):
            for i in range(nx):
                shear_v = 7.0 * np.cos(np.pi * znu[k]) * (1.0 + 0.4 * np.sin(2 * np.pi * znu[k]))
                v[i, k, j] = v_base[k] + shear_v * (1.0 + 0.5 * np.cos(2 * np.pi * i / nx) * np.sin(2 * np.pi * j / ny))

    # w: staggered in z (nx, nz + 1, ny)
    w = np.zeros((nx, nz + 1, ny))
    for j in range(ny):
        for k in range(nz + 1):
            for i in range(nx):
                w[i, k, j] = 0.5 * np.sin(2 * np.pi * i / nx) * np.sin(2 * np.pi * j / ny) * np.sin(np.pi * znw[k])

    # Spikes in w and ww at center to trigger w_damping
    ww = np.zeros((nx, nz + 1, ny))
    for j in range(ny):
        for k in range(nz + 1):
            for i in range(nx):
                c1c2 = c1f[k] * mut[i, j] + c2f[k]
                ww[i, k, j] = w[i, k, j] * c1c2 / 10.0
    # Put strong vertical velocities at (nx//2, ny//2, nz//2)
    ic, jc, kc = nx // 2, ny // 2, nz // 2
    w[ic, kc, jc] = 25.0
    ww[ic, kc, jc] = 25.0 * (c1f[kc] * mut[ic, jc] + c2f[kc]) / (rdnw[kc] * dt) * 2.5
    w[ic+1, kc, jc] = -22.0
    ww[ic+1, kc, jc] = -22.0 * (c1f[kc] * mut[ic+1, jc] + c2f[kc]) / (rdnw[kc] * dt) * 2.2

    # theta: potential temperature (nx, nz, ny)
    theta = np.zeros((nx, nz, ny))
    thp = np.zeros((nx, nz, ny))  # perturbation theta
    for j in range(ny):
        for k in range(nz):
            for i in range(nx):
                thp[i, k, j] = 2.0 * np.sin(2 * np.pi * i / nx) * np.cos(2 * np.pi * j / ny)
                theta[i, k, j] = t_base[k] + thp[i, k, j]

    # Sharp step in interior box [4..11, 4..11] to activate diff_6th_opt=2 limiter
    for j in range(4, 12):
        for k in range(nz):
            for i in range(4, 12):
                theta[i, k, j] += 15.0
                thp[i, k, j] += 15.0
    for j in range(4, 12):
        for k in range(nz):
            for i in range(4, 13):
                u[i, k, j] += 10.0
    for j in range(4, 13):
        for k in range(nz):
            for i in range(4, 12):
                v[i, k, j] += 10.0

    # Temperature t (perturbation: t = theta - 300)
    t = theta - 300.0

    # Pressure p and pb (nx, nz, ny)
    pb = np.zeros((nx, nz, ny))
    p = np.zeros((nx, nz, ny))
    for j in range(ny):
        for k in range(nz):
            for i in range(nx):
                pb[i, k, j] = c1[k] * mut[i, j] + c2[k] * p_top
                p[i, k, j] = 20.0 * np.sin(2 * np.pi * i / nx)

    # Inverse density alt and density rho (nx, nz, ny)
    alt = np.zeros((nx, nz, ny))
    rho = np.zeros((nx, nz, ny))
    for j in range(ny):
        for k in range(nz):
            for i in range(nx):
                total_p = pb[i, k, j] + p[i, k, j]
                total_th = theta[i, k, j]
                temp = total_th * (total_p / p1000mb)**(r_d / cp)
                density = total_p / (r_d * temp)
                rho[i, k, j] = density
                alt[i, k, j] = 1.0 / density

    # Moisture qv (nx, nz, ny)
    moist = np.zeros((nx, nz, ny))
    for j in range(ny):
        for k in range(nz):
            for i in range(nx):
                moist[i, k, j] = qv_base[k] * (1.0 + 0.1 * np.sin(2 * np.pi * i / nx))
                if 4 <= i < 12 and 4 <= j < 12:
                    moist[i, k, j] += 0.005

    # TKE (nx, nz + 1, ny)
    tke = np.zeros((nx, nz + 1, ny))
    for j in range(ny):
        for k in range(nz + 1):
            for i in range(nx):
                tke[i, k, j] = 1.0 + 0.5 * np.cos(np.pi * i / nx) * np.sin(np.pi * j / ny)

    # Surface fluxes
    ust = np.full((nx, ny), 0.3)
    hfx = np.full((nx, ny), 20.0)
    qfx = np.full((nx, ny), 0.0001)

    # Write binary file
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, 'wb') as f:
        # Header: nx, ny, nz, dx, dy, dt (standard alignment without compiler padding)
        f.write(struct.pack('=iiiddd', nx, ny, nz, dx, dy, dt))

        # 1D arrays (nz)
        f.write(dn.astype(np.float64).tobytes())
        f.write(rdn.astype(np.float64).tobytes())
        f.write(dnw.astype(np.float64).tobytes())
        f.write(rdnw.astype(np.float64).tobytes())
        f.write(fnm.astype(np.float64).tobytes())
        f.write(fnp.astype(np.float64).tobytes())
        f.write(struct.pack('=ddd', cf1, cf2, cf3))

        f.write(c1.astype(np.float64).tobytes())
        f.write(c2.astype(np.float64).tobytes())
        f.write(c3.astype(np.float64).tobytes())
        f.write(c4.astype(np.float64).tobytes())
        f.write(c1f.astype(np.float64).tobytes())
        f.write(c2f.astype(np.float64).tobytes())

        f.write(u_base.astype(np.float64).tobytes())
        f.write(v_base.astype(np.float64).tobytes())
        f.write(t_base.astype(np.float64).tobytes())
        f.write(qv_base.astype(np.float64).tobytes())

        # 2D arrays (Fortran column-major order: i fastest, j slowest)
        f.write(mut.astype(np.float64).tobytes(order='F'))
        f.write(muu.astype(np.float64).tobytes(order='F'))
        f.write(muv.astype(np.float64).tobytes(order='F'))
        f.write(msftx.astype(np.float64).tobytes(order='F'))
        f.write(msfty.astype(np.float64).tobytes(order='F'))
        f.write(msfux.astype(np.float64).tobytes(order='F'))
        f.write(msfuy.astype(np.float64).tobytes(order='F'))
        f.write(msfvx.astype(np.float64).tobytes(order='F'))
        f.write(msfvy.astype(np.float64).tobytes(order='F'))
        f.write(ust.astype(np.float64).tobytes(order='F'))
        f.write(hfx.astype(np.float64).tobytes(order='F'))
        f.write(qfx.astype(np.float64).tobytes(order='F'))

        # 3D arrays (Fortran column-major: i fastest, k middle, j slowest)
        # Note: in Python, shape is (nx, nz, ny) with order='F' -> i fastest, then k, then j!
        f.write(phb.astype(np.float64).tobytes(order='F'))
        f.write(ph.astype(np.float64).tobytes(order='F'))
        f.write(u.astype(np.float64).tobytes(order='F'))
        f.write(v.astype(np.float64).tobytes(order='F'))
        f.write(w.astype(np.float64).tobytes(order='F'))
        f.write(ww.astype(np.float64).tobytes(order='F'))
        f.write(thp.astype(np.float64).tobytes(order='F'))
        f.write(theta.astype(np.float64).tobytes(order='F'))
        f.write(t.astype(np.float64).tobytes(order='F'))
        f.write(pb.astype(np.float64).tobytes(order='F'))
        f.write(p.astype(np.float64).tobytes(order='F'))
        f.write(alt.astype(np.float64).tobytes(order='F'))
        f.write(rho.astype(np.float64).tobytes(order='F'))
        f.write(moist.astype(np.float64).tobytes(order='F'))
        f.write(tke.astype(np.float64).tobytes(order='F'))

    print(f"Generated input dataset at {output_path} (nx={nx}, ny={ny}, nz={nz})")

if __name__ == '__main__':
    out_file = sys.argv[1] if len(sys.argv) > 1 else 'diffusion_input.bin'
    generate_inputs(out_file)
