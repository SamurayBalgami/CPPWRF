#!/usr/bin/env python3
import argparse
import struct
import numpy as np

def generate_ysu_inputs(num_cols=2500, nz=50, dt=60.0, out_path="ysu_input.bin", seed=42):
    np.random.seed(seed)

    # Physical constants (WRF share/module_model_constants.F)
    g = 9.81
    rd = 287.0
    rv = 461.6
    cp = 7.0 * rd / 2.0  # 1004.5
    rovcp = rd / cp
    rovg = rd / g
    ep1 = rv / rd - 1.0
    ep2 = rd / rv
    karman = 0.4
    xlv = 2.5e6
    p1000mb = 100000.0
    svp1 = 0.6112
    svp2 = 17.67
    svp3 = 29.65
    svpt0 = 273.15

    ysu_topdown_pblmix = 0
    flag_qc = 1
    flag_qi = 1

    # 3D arrays: (nz, num_cols) in python, flattening to (nz, num_cols)
    u3d = np.zeros((nz, num_cols), dtype=np.float64)
    v3d = np.zeros((nz, num_cols), dtype=np.float64)
    t3d = np.zeros((nz, num_cols), dtype=np.float64)
    qv3d = np.zeros((nz, num_cols), dtype=np.float64)
    qc3d = np.zeros((nz, num_cols), dtype=np.float64)
    qi3d = np.zeros((nz, num_cols), dtype=np.float64)
    p3d = np.zeros((nz, num_cols), dtype=np.float64)
    p3di = np.zeros((nz + 1, num_cols), dtype=np.float64)
    pi3d = np.zeros((nz, num_cols), dtype=np.float64)
    dz8w = np.zeros((nz, num_cols), dtype=np.float64)
    rthraten = np.zeros((nz, num_cols), dtype=np.float64)

    # 2D arrays: (num_cols,)
    psfcpa = np.zeros(num_cols, dtype=np.float64)
    znt = np.zeros(num_cols, dtype=np.float64)
    ust = np.zeros(num_cols, dtype=np.float64)
    psim = np.zeros(num_cols, dtype=np.float64)
    psih = np.zeros(num_cols, dtype=np.float64)
    xland = np.zeros(num_cols, dtype=np.float64)
    hfx = np.zeros(num_cols, dtype=np.float64)
    qfx = np.zeros(num_cols, dtype=np.float64)
    wspd = np.zeros(num_cols, dtype=np.float64)
    br = np.zeros(num_cols, dtype=np.float64)
    u10 = np.zeros(num_cols, dtype=np.float64)
    v10 = np.zeros(num_cols, dtype=np.float64)
    uoce = np.zeros(num_cols, dtype=np.float64)
    voce = np.zeros(num_cols, dtype=np.float64)

    for i in range(num_cols):
        psfcpa[i] = np.random.uniform(96000.0, 102000.0)
        p3di[0, i] = psfcpa[i]

        # Geometric height grid: near surface ~25m up to ~15km
        z_interfaces = np.zeros(nz + 1)
        z_interfaces[0] = 0.0
        cur_dz = np.random.uniform(20.0, 35.0)
        for k in range(nz):
            z_interfaces[k + 1] = z_interfaces[k] + cur_dz
            dz8w[k, i] = cur_dz
            cur_dz = min(cur_dz * np.random.uniform(1.05, 1.10), 600.0)

        # Regimes across columns
        # 0 to 600: Unstable convective PBL (hfx > 0, br < 0)
        # 600 to 1200: Stable nocturnal PBL (hfx < 0, br > 0)
        # 1200 to 1800: Marine boundary layer (xland = 2, clouds at top of PBL)
        # 1800 to 2200: Neutral / strong wind shear
        # 2200 to 2500: Deep convection / supercooled / edge cases

        z1 = 0.5 * dz8w[0, i]

        if i < 600:
            # Daytime convective
            xland[i] = 1.0
            t_sfc = np.random.uniform(295.0, 310.0)
            gamma = np.random.uniform(0.0055, 0.0075)
            hfx[i] = np.random.uniform(80.0, 450.0)
            qfx[i] = np.random.uniform(2e-5, 2.5e-4)
            znt[i] = np.random.uniform(0.02, 0.4)
            ust[i] = np.random.uniform(0.15, 0.8)
            br[i] = np.random.uniform(-5.0, -0.05)
            lnz = np.log((z1 + znt[i]) / znt[i])
            psim[i] = max(lnz - np.random.uniform(0.5, 3.0), 0.5)
            psih[i] = max(lnz - np.random.uniform(1.0, 4.0), 0.5)

        elif i < 1200:
            # Stable nocturnal
            xland[i] = 1.0
            t_sfc = np.random.uniform(265.0, 285.0)
            gamma = np.random.uniform(0.004, 0.0065)
            hfx[i] = np.random.uniform(-40.0, -2.0)
            qfx[i] = np.random.uniform(0.0, 3e-5)
            znt[i] = np.random.uniform(0.01, 0.3)
            ust[i] = np.random.uniform(0.05, 0.35)
            br[i] = np.random.uniform(0.05, 5.0)
            lnz = np.log((z1 + znt[i]) / znt[i])
            psim[i] = lnz + np.random.uniform(0.5, 4.0)
            psih[i] = lnz + np.random.uniform(0.5, 5.0)

        elif i < 1800:
            # Marine boundary layer
            xland[i] = 2.0
            t_sfc = np.random.uniform(280.0, 302.0)
            gamma = 0.0065
            hfx[i] = np.random.uniform(-10.0, 120.0)
            qfx[i] = np.random.uniform(1e-5, 2e-4)
            znt[i] = np.random.uniform(1e-4, 2e-3)
            ust[i] = np.random.uniform(0.1, 0.9)
            br[i] = np.random.uniform(-1.0, 0.5)
            lnz = np.log((z1 + znt[i]) / znt[i])
            psim[i] = max(lnz + np.random.uniform(-1.0, 2.0), 0.5)
            psih[i] = max(lnz + np.random.uniform(-1.0, 2.0), 0.5)

        elif i < 2200:
            # Neutral / strong winds
            xland[i] = 1.0 if np.random.rand() > 0.5 else 2.0
            t_sfc = np.random.uniform(280.0, 298.0)
            gamma = 0.0065
            hfx[i] = np.random.uniform(-5.0, 25.0)
            qfx[i] = np.random.uniform(0.0, 4e-5)
            znt[i] = np.random.uniform(0.005, 0.2)
            ust[i] = np.random.uniform(0.6, 2.2)
            br[i] = np.random.uniform(-0.02, 0.02)
            lnz = np.log((z1 + znt[i]) / znt[i])
            psim[i] = max(lnz, 0.5)
            psih[i] = max(lnz, 0.5)

        else:
            # Edge cases
            xland[i] = 1.0 if np.random.rand() > 0.5 else 2.0
            t_sfc = np.random.uniform(240.0, 315.0)
            gamma = np.random.uniform(0.003, 0.009)
            hfx[i] = np.random.uniform(-60.0, 500.0)
            qfx[i] = np.random.uniform(0.0, 3e-4)
            znt[i] = np.random.uniform(1e-4, 0.5)
            ust[i] = np.random.uniform(0.01, 2.5)
            br[i] = np.random.uniform(-20.0, 20.0)
            lnz = np.log((z1 + znt[i]) / znt[i])
            psim[i] = max(lnz + np.random.uniform(-2.0, 3.0), 0.5)
            psih[i] = max(lnz + np.random.uniform(-2.0, 3.0), 0.5)

        # Build vertical profiles for this column
        u_top = np.random.uniform(10.0, 40.0)
        v_top = np.random.uniform(-20.0, 20.0)

        for k in range(nz):
            z_mid = 0.5 * (z_interfaces[k] + z_interfaces[k + 1])
            # Hydrostatic pressure
            p3di[k + 1, i] = psfcpa[i] * np.exp(-g * z_interfaces[k + 1] / (rd * (t_sfc - 0.5 * gamma * z_interfaces[k + 1])))
            p3d[k, i] = 0.5 * (p3di[k, i] + p3di[k + 1, i])
            pi3d[k, i] = (p3d[k, i] / p1000mb) ** rovcp

            # Temperature with surface inversion if stable
            if i >= 600 and i < 1200 and z_mid < 400.0:
                # Strong inversion in lowest 400m
                t3d[k, i] = t_sfc + (z_mid / 400.0) * np.random.uniform(3.0, 8.0)
            else:
                t3d[k, i] = t_sfc - gamma * z_mid

            # Winds
            shear_fac = (z_mid / 10000.0) ** 0.4
            u3d[k, i] = u_top * shear_fac + np.random.uniform(-0.5, 0.5)
            v3d[k, i] = v_top * shear_fac + np.random.uniform(-0.5, 0.5)

            # Vapor mixing ratio
            es = svp1 * 1000.0 * np.exp(svp2 * (t3d[k, i] - svpt0) / (t3d[k, i] - svp3))
            qvs = ep2 * es / (p3d[k, i] - es)
            rh = np.random.uniform(0.1, 0.85)
            qv3d[k, i] = max(rh * qvs, 1e-6)

            # Clouds
            if z_mid > 800.0 and z_mid < 2500.0 and np.random.rand() > 0.6:
                qc3d[k, i] = np.random.uniform(0.05e-3, 1.2e-3)
            else:
                qc3d[k, i] = 0.0

            if t3d[k, i] < 265.0 and np.random.rand() > 0.7:
                qi3d[k, i] = np.random.uniform(0.01e-3, 0.5e-3)
            else:
                qi3d[k, i] = 0.0

            # Radiative cooling (K/s)
            rthraten[k, i] = np.random.uniform(-2.5e-5, -0.5e-5)

        wspd[i] = np.sqrt(u3d[0, i] ** 2 + v3d[0, i] ** 2)
        u10[i] = u3d[0, i] * 0.75
        v10[i] = v3d[0, i] * 0.75

    # Write binary file
    with open(out_path, "wb") as f:
        # Dimensions
        f.write(struct.pack("<ii", num_cols, nz))
        # Timestep
        f.write(struct.pack("<d", dt))
        # Flags
        f.write(struct.pack("<iii", ysu_topdown_pblmix, flag_qc, flag_qi))
        # Constants
        f.write(struct.pack("<dddddddddd", cp, g, rovcp, rd, rovg, ep1, ep2, karman, xlv, rv))

        # 3D arrays: write in Fortran column-major order (k changes slower: all cols for level k)
        # Specifically: for k in range(nz): write array[k, :]
        f.write(u3d.tobytes())
        f.write(v3d.tobytes())
        f.write(t3d.tobytes())
        f.write(qv3d.tobytes())
        f.write(qc3d.tobytes())
        f.write(qi3d.tobytes())
        f.write(p3d.tobytes())
        f.write(p3di.tobytes())
        f.write(pi3d.tobytes())
        f.write(dz8w.tobytes())
        f.write(rthraten.tobytes())

        # 2D arrays
        f.write(psfcpa.tobytes())
        f.write(znt.tobytes())
        f.write(ust.tobytes())
        f.write(psim.tobytes())
        f.write(psih.tobytes())
        f.write(xland.tobytes())
        f.write(hfx.tobytes())
        f.write(qfx.tobytes())
        f.write(wspd.tobytes())
        f.write(br.tobytes())
        f.write(u10.tobytes())
        f.write(v10.tobytes())
        f.write(uoce.tobytes())
        f.write(voce.tobytes())

    print(f"Generated YSU input: {num_cols} cols x {nz} levels into {out_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cols", type=int, default=2500)
    parser.add_argument("--nz", type=int, default=50)
    parser.add_argument("--dt", type=float, default=60.0)
    parser.add_argument("--out", type=str, default="ysu_input.bin")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    generate_ysu_inputs(args.cols, args.nz, args.dt, args.out, args.seed)
