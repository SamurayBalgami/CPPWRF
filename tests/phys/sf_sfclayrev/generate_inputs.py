#!/usr/bin/env python3
import argparse
import struct
import numpy as np

def generate_sfclayrev_inputs(num_cols=2500, nz=50, out_path="sfclayrev_input.bin", seed=42):
    np.random.seed(seed)

    # Physical constants (from phys_constants.hpp / WRF share/module_model_constants.F)
    g = 9.81
    r_d = 287.0
    r_v = 461.6
    cp = 7.0 * r_d / 2.0  # 1004.5
    rovcp = r_d / cp
    r = r_d
    xlv = 2.5e6
    svp1 = 0.6112
    svp2 = 17.67
    svp3 = 29.65
    svpt0 = 273.15
    ep1 = r_v / r_d - 1.0
    ep2 = r_d / r_v
    karman = 0.4
    p1000mb = 100000.0

    isfflx = 1
    shalwater_z0 = 0
    isftcflx = 0
    iz0tlnd = 0
    scm_force_flux = 0

    # Allocate 1D arrays of size num_cols
    u1d = np.zeros(num_cols, dtype=np.float64)
    v1d = np.zeros(num_cols, dtype=np.float64)
    t1d = np.zeros(num_cols, dtype=np.float64)
    qv1d = np.zeros(num_cols, dtype=np.float64)
    p1d = np.zeros(num_cols, dtype=np.float64)
    dz8w1d = np.zeros(num_cols, dtype=np.float64)
    psfcpa = np.zeros(num_cols, dtype=np.float64)
    tsk = np.zeros(num_cols, dtype=np.float64)
    xland = np.zeros(num_cols, dtype=np.float64)
    dx = np.zeros(num_cols, dtype=np.float64)
    mavail = np.zeros(num_cols, dtype=np.float64)
    pblh = np.zeros(num_cols, dtype=np.float64)
    lakemask = np.zeros(num_cols, dtype=np.float64)
    water_depth = np.zeros(num_cols, dtype=np.float64)

    regime = np.zeros(num_cols, dtype=np.float64)
    hfx = np.zeros(num_cols, dtype=np.float64)
    qfx = np.zeros(num_cols, dtype=np.float64)
    qsfc = np.zeros(num_cols, dtype=np.float64)
    mol = np.zeros(num_cols, dtype=np.float64)
    rmol = np.zeros(num_cols, dtype=np.float64)
    gz1oz0 = np.zeros(num_cols, dtype=np.float64)
    wspd = np.zeros(num_cols, dtype=np.float64)
    br = np.zeros(num_cols, dtype=np.float64)
    psim = np.zeros(num_cols, dtype=np.float64)
    psih = np.zeros(num_cols, dtype=np.float64)
    fm = np.zeros(num_cols, dtype=np.float64)
    fh = np.zeros(num_cols, dtype=np.float64)
    znt = np.zeros(num_cols, dtype=np.float64)
    zol = np.zeros(num_cols, dtype=np.float64)
    ust = np.zeros(num_cols, dtype=np.float64)
    cpm = np.zeros(num_cols, dtype=np.float64)
    chs2 = np.zeros(num_cols, dtype=np.float64)
    cqs2 = np.zeros(num_cols, dtype=np.float64)
    chs = np.zeros(num_cols, dtype=np.float64)
    flhc = np.zeros(num_cols, dtype=np.float64)
    flqc = np.zeros(num_cols, dtype=np.float64)
    qgh = np.zeros(num_cols, dtype=np.float64)
    ustm = np.zeros(num_cols, dtype=np.float64)

    # Segment columns into distinct regimes
    # 0 to 500: Unstable convective land (daytime)
    # 500 to 1000: Stable nocturnal land (nighttime)
    # 1000 to 1500: Marine boundary layer (xland = 2.0)
    # 1500 to 2000: Neutral / forced convection & extreme winds
    # 2000 to 2500: Edge cases (calm winds, ust < 0.001, large br, threshold boundaries)

    for i in range(num_cols):
        psfcpa[i] = np.random.uniform(95000.0, 102500.0)
        p1d[i] = psfcpa[i] - np.random.uniform(200.0, 800.0)
        dz8w1d[i] = np.random.uniform(20.0, 80.0)
        dx[i] = np.random.uniform(1000.0, 30000.0)
        mavail[i] = np.random.uniform(0.1, 1.0)
        pblh[i] = np.random.uniform(200.0, 2500.0)
        water_depth[i] = np.random.uniform(2.0, 200.0)

        # Default physical values
        znt[i] = np.random.uniform(0.001, 0.1)
        ust[i] = np.random.uniform(0.05, 0.5)
        hfx[i] = 0.0
        qfx[i] = 0.0
        mol[i] = 0.0

        if i < 500:
            # Unstable convective land
            xland[i] = 1.0
            lakemask[i] = 0.0
            tsk[i] = np.random.uniform(290.0, 315.0)
            t1d[i] = tsk[i] - np.random.uniform(2.0, 15.0)  # tsk > t1 -> unstable
            u1d[i] = np.random.uniform(-10.0, 10.0)
            v1d[i] = np.random.uniform(-10.0, 10.0)
            qv1d[i] = np.random.uniform(0.005, 0.020)
            znt[i] = np.random.uniform(0.01, 0.5)
            ust[i] = np.random.uniform(0.1, 0.8)
            hfx[i] = np.random.uniform(50.0, 400.0)
            qfx[i] = np.random.uniform(1e-5, 2e-4)
            mol[i] = np.random.uniform(-2.0, -0.01)

        elif i < 1000:
            # Stable nocturnal land
            xland[i] = 1.0
            lakemask[i] = 0.0
            tsk[i] = np.random.uniform(260.0, 290.0)
            t1d[i] = tsk[i] + np.random.uniform(1.0, 12.0)  # tsk < t1 -> stable
            u1d[i] = np.random.uniform(-8.0, 8.0)
            v1d[i] = np.random.uniform(-8.0, 8.0)
            qv1d[i] = np.random.uniform(0.001, 0.008)
            znt[i] = np.random.uniform(0.01, 0.4)
            ust[i] = np.random.uniform(0.05, 0.4)
            hfx[i] = np.random.uniform(-50.0, -1.0)
            qfx[i] = np.random.uniform(0.0, 5e-5)
            mol[i] = np.random.uniform(0.01, 1.5)

        elif i < 1500:
            # Marine boundary layer
            xland[i] = 2.0
            lakemask[i] = 0.0 if np.random.rand() > 0.3 else 1.0
            tsk[i] = np.random.uniform(280.0, 305.0)
            t1d[i] = tsk[i] + np.random.uniform(-4.0, 4.0)
            u1d[i] = np.random.uniform(-15.0, 15.0)
            v1d[i] = np.random.uniform(-15.0, 15.0)
            qv1d[i] = np.random.uniform(0.005, 0.022)
            znt[i] = np.random.uniform(1e-4, 2e-3)
            ust[i] = np.random.uniform(0.1, 1.2)
            hfx[i] = np.random.uniform(-30.0, 150.0)
            qfx[i] = np.random.uniform(1e-5, 3e-4)
            mol[i] = np.random.uniform(-1.0, 1.0)

        elif i < 2000:
            # Neutral / forced convection & high winds
            xland[i] = 1.0 if np.random.rand() > 0.5 else 2.0
            lakemask[i] = 0.0
            tsk[i] = np.random.uniform(275.0, 300.0)
            t1d[i] = tsk[i] + np.random.uniform(-0.1, 0.1)  # near neutral
            w_spd = np.random.uniform(15.0, 45.0)
            ang = np.random.uniform(0.0, 2.0 * np.pi)
            u1d[i] = w_spd * np.cos(ang)
            v1d[i] = w_spd * np.sin(ang)
            qv1d[i] = np.random.uniform(0.002, 0.015)
            znt[i] = np.random.uniform(0.001, 0.2)
            ust[i] = np.random.uniform(0.5, 2.0)
            hfx[i] = np.random.uniform(-10.0, 20.0)
            qfx[i] = np.random.uniform(0.0, 5e-5)
            mol[i] = 0.0

        else:
            # Edge cases
            xland[i] = 1.0 if np.random.rand() > 0.5 else 2.0
            lakemask[i] = 0.0
            sub = i - 2000
            if sub < 50:
                # Near calm winds
                u1d[i] = np.random.uniform(-0.02, 0.02)
                v1d[i] = np.random.uniform(-0.02, 0.02)
                tsk[i] = 295.0
                t1d[i] = 295.0
                ust[i] = 0.0005  # ust < 0.001
            elif sub < 100:
                # Strong stable br > 250
                u1d[i] = 0.05
                v1d[i] = 0.05
                tsk[i] = 260.0
                t1d[i] = 290.0  # large inversion
                ust[i] = 0.05
            elif sub < 150:
                # Extreme unstable br < -250
                u1d[i] = 0.05
                v1d[i] = 0.05
                tsk[i] = 315.0
                t1d[i] = 270.0  # huge lapse rate
                ust[i] = 0.05
            elif sub < 250:
                # Zero / tiny mixing ratios
                u1d[i] = np.random.uniform(1.0, 5.0)
                v1d[i] = np.random.uniform(1.0, 5.0)
                tsk[i] = 250.0
                t1d[i] = 250.0
                qv1d[i] = 1e-12
                ust[i] = 0.2
            else:
                # Random combinations
                u1d[i] = np.random.uniform(-25.0, 25.0)
                v1d[i] = np.random.uniform(-25.0, 25.0)
                tsk[i] = np.random.uniform(240.0, 320.0)
                t1d[i] = np.random.uniform(240.0, 320.0)
                qv1d[i] = np.random.uniform(1e-6, 0.025)
                znt[i] = np.random.uniform(1e-5, 1.0)
                ust[i] = np.random.uniform(0.0001, 2.5)

        ustm[i] = ust[i]
        cpm[i] = cp * (1.0 + 0.8 * qv1d[i])
        qsfc[i] = qv1d[i]

    # Pack into binary file
    with open(out_path, "wb") as f:
        # Dimensions
        f.write(struct.pack("<ii", num_cols, nz))
        # Option flags
        f.write(struct.pack("<iiiii", isfflx, shalwater_z0, isftcflx, iz0tlnd, scm_force_flux))
        # Constants
        f.write(struct.pack("<ddddddddddddd", cp, g, rovcp, r, xlv, svp1, svp2, svp3, svpt0, ep1, ep2, karman, p1000mb))

        # Array data
        arrays_to_write = [
            u1d, v1d, t1d, qv1d, p1d, dz8w1d, psfcpa, tsk, xland, dx, mavail, pblh, lakemask, water_depth,
            regime, hfx, qfx, qsfc, mol, rmol, gz1oz0, wspd, br, psim, psih, fm, fh,
            znt, zol, ust, cpm, chs2, cqs2, chs, flhc, flqc, qgh, ustm
        ]
        for arr in arrays_to_write:
            f.write(arr.tobytes())

    print(f"Generated {num_cols} columns into {out_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cols", type=int, default=2500)
    parser.add_argument("--nz", type=int, default=50)
    parser.add_argument("--out", type=str, default="sfclayrev_input.bin")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    generate_sfclayrev_inputs(args.cols, args.nz, args.out, args.seed)
