#!/usr/bin/env python3
"""
Generate diverse input regimes for Noah LSM differential testing.
Generates >= 2000 columns covering:
1. Bare soil
2. Diverse vegetation classes (crops, grass, forests, shrubs)
3. Snow cover (cold snow, melting snow)
4. Frozen soil (subfreezing soil profiles, phase change)
5. Water points (skipped)
6. Dry and saturated soils
7. Day and night forcing (zero vs high solar, downward LW, rain, freezing rain, snowfall)
"""

import struct
import numpy as np

def generate_noah_inputs(filename="noah_input.bin", n_cols=2500, seed=42):
    np.random.seed(seed)
    nsoil = 4
    dt = 60.0  # 60s timestep
    zlvl = 10.0  # 10m atmospheric forcing height
    sldpth = np.array([0.1, 0.3, 0.6, 1.0], dtype=np.float64)
    isurban = 1

    # Proportional regime allocation
    c1 = int(n_cols * 400 / 2500)
    c2 = int(n_cols * 1200 / 2500)
    c3 = int(n_cols * 1600 / 2500)
    c4 = int(n_cols * 2000 / 2500)
    c5 = int(n_cols * 2100 / 2500)
    c6 = int(n_cols * 2300 / 2500)

    vegtyp = np.zeros(n_cols, dtype=np.int32)
    soiltyp = np.zeros(n_cols, dtype=np.int32)
    slopetyp = np.random.randint(1, 10, size=n_cols, dtype=np.int32)

    shdfac = np.zeros(n_cols, dtype=np.float64)
    shdmin = np.full(n_cols, 0.05, dtype=np.float64)
    shdmax = np.full(n_cols, 0.85, dtype=np.float64)
    alb = np.zeros(n_cols, dtype=np.float64)
    snoalb = np.zeros(n_cols, dtype=np.float64)
    tbot = np.zeros(n_cols, dtype=np.float64)
    z0brd = np.zeros(n_cols, dtype=np.float64)
    embrd = np.zeros(n_cols, dtype=np.float64)

    lwdn = np.zeros(n_cols, dtype=np.float64)
    soldn = np.zeros(n_cols, dtype=np.float64)
    solnet = np.zeros(n_cols, dtype=np.float64)
    sfcprs = np.zeros(n_cols, dtype=np.float64)
    prcp = np.zeros(n_cols, dtype=np.float64)
    sfctmp = np.zeros(n_cols, dtype=np.float64)
    th2 = np.zeros(n_cols, dtype=np.float64)
    q2 = np.zeros(n_cols, dtype=np.float64)
    sfcspd = np.zeros(n_cols, dtype=np.float64)
    q2sat = np.zeros(n_cols, dtype=np.float64)
    dqsdt2 = np.zeros(n_cols, dtype=np.float64)
    ffrozp = np.zeros(n_cols, dtype=np.float64)

    cmc = np.zeros(n_cols, dtype=np.float64)
    t1 = np.zeros(n_cols, dtype=np.float64)
    snowh = np.zeros(n_cols, dtype=np.float64)
    sneqv = np.zeros(n_cols, dtype=np.float64)
    albedo = np.zeros(n_cols, dtype=np.float64)
    ch = np.zeros(n_cols, dtype=np.float64)
    cm = np.zeros(n_cols, dtype=np.float64)
    emissi = np.zeros(n_cols, dtype=np.float64)
    z0 = np.zeros(n_cols, dtype=np.float64)
    xlai = np.zeros(n_cols, dtype=np.float64)
    snotime1 = np.zeros(n_cols, dtype=np.float64)
    ribb = np.zeros(n_cols, dtype=np.float64)

    stc = np.zeros((nsoil, n_cols), dtype=np.float64)
    smc = np.zeros((nsoil, n_cols), dtype=np.float64)
    sh2o = np.zeros((nsoil, n_cols), dtype=np.float64)

    def calc_qsat(temp, pres):
        # Saturation vapor pressure (Tetens formula matching WRF)
        es = 611.2 * np.exp(17.67 * (temp - 273.15) / (temp - 29.65))
        qs = 0.622 * es / (pres - 0.378 * es)
        # dqs/dT
        desdt = es * 17.67 * (273.15 - 29.65) / ((temp - 29.65) ** 2)
        dqs = 0.622 * pres * desdt / ((pres - 0.378 * es) ** 2)
        return qs, dqs

    for c in range(n_cols):
        # Defaults
        sfcprs[c] = np.random.uniform(80000.0, 101325.0)
        sfcspd[c] = np.random.uniform(0.5, 20.0)
        ch[c] = np.random.uniform(0.005, 0.05)
        cm[c] = np.random.uniform(0.005, 0.05)
        ribb[c] = np.random.uniform(-0.1, 0.5)
        snoalb[c] = np.random.uniform(0.6, 0.85)
        z0brd[c] = np.random.uniform(0.01, 0.5)
        embrd[c] = np.random.uniform(0.88, 0.98)
        alb[c] = np.random.uniform(0.12, 0.25)
        albedo[c] = alb[c]
        emissi[c] = embrd[c]
        z0[c] = z0brd[c]
        xlai[c] = np.random.uniform(0.5, 4.0)

        # Regimes
        if c < c1:
            # 1. Bare soil
            vegtyp[c] = 19  # Barren or Sparsely Vegetated
            soiltyp[c] = np.random.randint(1, 13)
            shdfac[c] = np.random.uniform(0.0, 0.05)
            sfctmp[c] = np.random.uniform(280.0, 310.0)
            t1[c] = sfctmp[c] + np.random.uniform(-5.0, 10.0)
            tbot[c] = np.random.uniform(285.0, 295.0)
            for k in range(nsoil):
                stc[k, c] = sfctmp[c] + np.random.uniform(-2.0, 5.0)
                smc[k, c] = np.random.uniform(0.05, 0.35)
                sh2o[k, c] = smc[k, c]
            if np.random.rand() > 0.4:
                soldn[c] = np.random.uniform(200.0, 900.0)
                solnet[c] = soldn[c] * (1.0 - alb[c])
            else:
                soldn[c] = 0.0
                solnet[c] = 0.0
            lwdn[c] = np.random.uniform(250.0, 420.0)
            prcp[c] = np.random.uniform(0.0, 0.002) if np.random.rand() > 0.7 else 0.0

        elif c < c2:
            # 2. Vegetated land
            vegtyp[c] = np.random.choice([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18])
            soiltyp[c] = np.random.randint(1, 13)
            shdfac[c] = np.random.uniform(0.1, 0.9)
            sfctmp[c] = np.random.uniform(275.0, 305.0)
            t1[c] = sfctmp[c] + np.random.uniform(-3.0, 6.0)
            tbot[c] = np.random.uniform(280.0, 295.0)
            cmc[c] = np.random.uniform(0.0, 0.0005) if np.random.rand() > 0.5 else 0.0
            for k in range(nsoil):
                stc[k, c] = sfctmp[c] + np.random.uniform(-3.0, 3.0)
                smc[k, c] = np.random.uniform(0.10, 0.40)
                sh2o[k, c] = smc[k, c]
            if np.random.rand() > 0.3:
                soldn[c] = np.random.uniform(100.0, 950.0)
                solnet[c] = soldn[c] * (1.0 - alb[c])
            else:
                soldn[c] = 0.0
                solnet[c] = 0.0
            lwdn[c] = np.random.uniform(220.0, 400.0)
            prcp[c] = np.random.uniform(0.0, 0.004) if np.random.rand() > 0.6 else 0.0

        elif c < c3:
            # 3. Snow cover
            vegtyp[c] = np.random.choice([2, 5, 7, 8, 11, 14, 20, 21, 22])
            soiltyp[c] = np.random.randint(1, 13)
            shdfac[c] = np.random.uniform(0.1, 0.7)
            snowh[c] = np.random.uniform(0.02, 1.2)
            sneqv[c] = snowh[c] * np.random.uniform(0.1, 0.35)
            snotime1[c] = np.random.uniform(0.0, 200000.0)
            if c < c2 + (c3 - c2) // 2:
                # Subfreezing snow
                sfctmp[c] = np.random.uniform(240.0, 270.0)
                t1[c] = sfctmp[c] + np.random.uniform(-4.0, 2.0)
            else:
                # Melting snow
                sfctmp[c] = np.random.uniform(273.2, 280.0)
                t1[c] = 273.15 + np.random.uniform(-0.5, 2.0)
            tbot[c] = np.random.uniform(265.0, 280.0)
            for k in range(nsoil):
                stc[k, c] = t1[c] + np.random.uniform(-2.0, 4.0)
                smc[k, c] = np.random.uniform(0.15, 0.40)
                sh2o[k, c] = smc[k, c] * (0.8 if stc[k, c] > 273.15 else 0.2)
            soldn[c] = np.random.uniform(0.0, 600.0)
            solnet[c] = soldn[c] * (1.0 - 0.7)
            lwdn[c] = np.random.uniform(150.0, 320.0)
            if np.random.rand() > 0.5:
                prcp[c] = np.random.uniform(0.0001, 0.002)
                ffrozp[c] = 1.0 if sfctmp[c] <= 273.15 else 0.0

        elif c < c4:
            # 4. Frozen soil
            vegtyp[c] = np.random.choice([7, 8, 9, 20, 21, 22, 23])
            soiltyp[c] = np.random.randint(1, 13)
            shdfac[c] = np.random.uniform(0.05, 0.5)
            sfctmp[c] = np.random.uniform(235.0, 268.0)
            t1[c] = sfctmp[c] + np.random.uniform(-3.0, 3.0)
            tbot[c] = np.random.uniform(260.0, 272.0)
            for k in range(nsoil):
                stc[k, c] = np.random.uniform(245.0, 271.0)
                smc[k, c] = np.random.uniform(0.15, 0.42)
                sh2o[k, c] = smc[k, c] * np.random.uniform(0.05, 0.25)
            soldn[c] = np.random.uniform(0.0, 400.0)
            solnet[c] = soldn[c] * (1.0 - alb[c])
            lwdn[c] = np.random.uniform(140.0, 260.0)
            if np.random.rand() > 0.6:
                prcp[c] = np.random.uniform(0.0001, 0.001)
                ffrozp[c] = 1.0

        elif c < c5:
            # 5. Water points (skipped)
            vegtyp[c] = 16 if (c % 2 == 0) else np.random.randint(1, 10)
            soiltyp[c] = 14 if (c % 2 != 0) else np.random.randint(1, 10)
            sfctmp[c] = np.random.uniform(280.0, 300.0)
            t1[c] = sfctmp[c]
            tbot[c] = 285.0
            for k in range(nsoil):
                stc[k, c] = 285.0
                smc[k, c] = 0.3
                sh2o[k, c] = 0.3

        elif c < c6:
            # 6. Extreme dry and saturated soils
            vegtyp[c] = np.random.choice([1, 2, 7, 8, 19])
            soiltyp[c] = np.random.randint(1, 13)
            shdfac[c] = np.random.uniform(0.0, 0.8)
            sfctmp[c] = np.random.uniform(275.0, 310.0)
            t1[c] = sfctmp[c] + np.random.uniform(-4.0, 8.0)
            tbot[c] = 290.0
            is_dry = (c % 2 == 0)
            for k in range(nsoil):
                stc[k, c] = sfctmp[c]
                smc[k, c] = np.random.uniform(0.015, 0.04) if is_dry else np.random.uniform(0.44, 0.48)
                sh2o[k, c] = smc[k, c]
            soldn[c] = np.random.uniform(0.0, 1000.0)
            solnet[c] = soldn[c] * (1.0 - alb[c])
            lwdn[c] = np.random.uniform(200.0, 450.0)
            prcp[c] = np.random.uniform(0.0, 0.005) if not is_dry else 0.0

        else:
            # 7. Mixed precipitation, dewfall, freezing rain
            vegtyp[c] = np.random.randint(1, 15)
            soiltyp[c] = np.random.randint(1, 13)
            shdfac[c] = np.random.uniform(0.2, 0.8)
            tbot[c] = 280.0
            if c < c6 + (n_cols - c6) // 2:
                # Freezing rain: sfctmp > 273.15, t1 <= 273.15, prcp > 0, ffrozp = 0
                sfctmp[c] = np.random.uniform(273.5, 276.0)
                t1[c] = np.random.uniform(268.0, 272.5)
                prcp[c] = np.random.uniform(0.0005, 0.003)
                ffrozp[c] = 0.0
            else:
                # Dewfall / frostfall: negative ETP conditions
                sfctmp[c] = np.random.uniform(265.0, 285.0)
                t1[c] = sfctmp[c] - np.random.uniform(1.0, 5.0)
                prcp[c] = 0.0
                soldn[c] = 0.0
                solnet[c] = 0.0
                lwdn[c] = np.random.uniform(150.0, 250.0)
            for k in range(nsoil):
                stc[k, c] = t1[c]
                smc[k, c] = np.random.uniform(0.2, 0.35)
                sh2o[k, c] = smc[k, c] * (1.0 if stc[k, c] > 273.15 else 0.3)

        # Compute humidity and theta
        th2[c] = sfctmp[c] * ((100000.0 / sfcprs[c]) ** 0.286)
        qs, dqs = calc_qsat(sfctmp[c], sfcprs[c])
        q2sat[c] = qs
        dqsdt2[c] = dqs
        rh = np.random.uniform(0.2, 0.98)
        q2[c] = min(rh * qs, qs * 0.99)

    # Write binary input file
    with open(filename, "wb") as f:
        f.write(struct.pack("=iiddi", n_cols, nsoil, dt, zlvl, isurban))
        f.write(sldpth.tobytes())
        f.write(vegtyp.tobytes())
        f.write(soiltyp.tobytes())
        f.write(slopetyp.tobytes())
        f.write(shdfac.tobytes())
        f.write(shdmin.tobytes())
        f.write(shdmax.tobytes())
        f.write(alb.tobytes())
        f.write(snoalb.tobytes())
        f.write(tbot.tobytes())
        f.write(z0brd.tobytes())
        f.write(embrd.tobytes())
        f.write(lwdn.tobytes())
        f.write(soldn.tobytes())
        f.write(solnet.tobytes())
        f.write(sfcprs.tobytes())
        f.write(prcp.tobytes())
        f.write(sfctmp.tobytes())
        f.write(th2.tobytes())
        f.write(q2.tobytes())
        f.write(sfcspd.tobytes())
        f.write(q2sat.tobytes())
        f.write(dqsdt2.tobytes())
        f.write(ffrozp.tobytes())

        f.write(cmc.tobytes())
        f.write(t1.tobytes())
        f.write(snowh.tobytes())
        f.write(sneqv.tobytes())
        f.write(albedo.tobytes())
        f.write(ch.tobytes())
        f.write(cm.tobytes())
        f.write(emissi.tobytes())
        f.write(z0.tobytes())
        f.write(xlai.tobytes())
        f.write(snotime1.tobytes())
        f.write(ribb.tobytes())
        f.write(stc.tobytes(order="C"))
        f.write(smc.tobytes(order="C"))
        f.write(sh2o.tobytes(order="C"))

    print(f"Generated {n_cols} columns into {filename}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", nargs="?", default="noah_input.bin")
    parser.add_argument("-o", "--output", dest="out_fn", default=None)
    parser.add_argument("--ncols", type=int, default=2500)
    args = parser.parse_args()
    fn = args.out_fn if args.out_fn else args.filename
    generate_noah_inputs(fn, n_cols=args.ncols)
