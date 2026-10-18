#!/usr/bin/env python3
"""
Generate synthetic multi-regime column inputs for Noah-MP v5.1.2 differential testing.
Covers all MODIS_IGBP_MODIS_NOAH landuse categories (vegetated, cropland, wetland,
urban, barren, glacier/ice, water), soil textures (1..12), snow regimes (0, 1, 2, 3
layers, thin snow, melting snow), frozen/warm soil, daytime/nighttime radiation,
and convective/non-convective/snow/graupel/hail precipitation.
"""

import argparse
import struct
import numpy as np


def generate_inputs(filepath, n_cols=2500, n_steps=5):
    rng = np.random.default_rng(42)
    nsoil = 4
    yr = 2026
    dt = 30.0
    dx = 9000.0
    julian_start = 284.5
    xice_threshold = 0.02
    dzs = np.array([0.10, 0.30, 0.60, 1.00], dtype=np.float64)

    ivgtyp = np.zeros(n_cols, dtype=np.int32)
    isltyp = np.zeros(n_cols, dtype=np.int32)
    xlat = np.zeros(n_cols, dtype=np.float64)
    xlong = np.zeros(n_cols, dtype=np.float64)
    coszen = np.zeros(n_cols, dtype=np.float64)
    vegfra = np.zeros(n_cols, dtype=np.float64)
    shdmax = np.zeros(n_cols, dtype=np.float64)
    tmn = np.zeros(n_cols, dtype=np.float64)
    xland = np.ones(n_cols, dtype=np.float64)
    xice = np.zeros(n_cols, dtype=np.float64)

    swdown = np.zeros(n_cols, dtype=np.float64)
    glw = np.zeros(n_cols, dtype=np.float64)
    rainbl = np.zeros(n_cols, dtype=np.float64)
    sr = np.zeros(n_cols, dtype=np.float64)
    mp_rainc = np.zeros(n_cols, dtype=np.float64)
    mp_rainnc = np.zeros(n_cols, dtype=np.float64)
    mp_shcv = np.zeros(n_cols, dtype=np.float64)
    mp_snow = np.zeros(n_cols, dtype=np.float64)
    mp_graup = np.zeros(n_cols, dtype=np.float64)
    mp_hail = np.zeros(n_cols, dtype=np.float64)

    dz8w_1 = np.zeros(n_cols, dtype=np.float64)
    t_phy_1 = np.zeros(n_cols, dtype=np.float64)
    qv_curr_1 = np.zeros(n_cols, dtype=np.float64)
    u_phy_1 = np.zeros(n_cols, dtype=np.float64)
    v_phy_1 = np.zeros(n_cols, dtype=np.float64)
    p8w_1 = np.zeros(n_cols, dtype=np.float64)
    p8w_2 = np.zeros(n_cols, dtype=np.float64)

    tsk = np.zeros(n_cols, dtype=np.float64)
    snow = np.zeros(n_cols, dtype=np.float64)
    snowh = np.zeros(n_cols, dtype=np.float64)
    canwat = np.zeros(n_cols, dtype=np.float64)
    lai = np.zeros(n_cols, dtype=np.float64)

    tslb = np.zeros((nsoil, n_cols), dtype=np.float64)
    smois = np.zeros((nsoil, n_cols), dtype=np.float64)
    sh2o = np.zeros((nsoil, n_cols), dtype=np.float64)

    # Land vegetation types (excluding 17=water)
    land_veg_types = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20]

    for i in range(n_cols):
        regime = i % 10
        if regime == 0:
            # Warm daytime vegetated, no snow, no rain
            vt = land_veg_types[i % 12]  # 1..12
            st = (i % 12) + 1
            t_air = rng.uniform(285.0, 305.0)
            cz = rng.uniform(0.15, 0.95)
            sw = cz * rng.uniform(300.0, 950.0)
            swe = 0.0
            sdepth = 0.0
            prcp = 0.0
            s_ratio = 0.0
        elif regime == 1:
            # Nighttime vegetated with dew/frost conditions
            vt = land_veg_types[i % 12]
            st = (i % 12) + 1
            t_air = rng.uniform(268.0, 290.0)
            cz = rng.uniform(-0.8, 0.0)
            sw = 0.0
            swe = 0.0
            sdepth = 0.0
            prcp = 0.0
            s_ratio = 0.0
        elif regime == 2:
            # Heavy rain on vegetated / cropland
            vt = land_veg_types[i % 14]
            st = (i % 12) + 1
            t_air = rng.uniform(280.0, 298.0)
            cz = rng.uniform(0.05, 0.7)
            sw = cz * rng.uniform(100.0, 400.0)
            swe = 0.0
            sdepth = 0.0
            prcp = rng.uniform(0.05, 1.5)  # mm per step
            s_ratio = 0.0
        elif regime == 3:
            # Thin snow (0 layers: snowh < 0.025) and 1-layer snow (0.025 <= snowh <= 0.05)
            vt = land_veg_types[i % 12]
            st = (i % 12) + 1
            t_air = rng.uniform(260.0, 275.0)
            cz = rng.uniform(0.1, 0.6)
            sw = cz * rng.uniform(150.0, 500.0)
            sdepth = rng.uniform(0.005, 0.048)
            swe = sdepth * rng.uniform(150.0, 250.0)
            prcp = rng.uniform(0.0, 0.2)
            s_ratio = 1.0 if t_air < 272.0 else 0.3
        elif regime == 4:
            # 2-layer snow (0.05 < snowh <= 0.25) with snowfall
            vt = land_veg_types[i % 12]
            st = (i % 12) + 1
            t_air = rng.uniform(252.0, 271.0)
            cz = rng.uniform(0.0, 0.5)
            sw = cz * rng.uniform(100.0, 400.0)
            sdepth = rng.uniform(0.055, 0.24)
            swe = sdepth * rng.uniform(180.0, 280.0)
            prcp = rng.uniform(0.05, 0.6)
            s_ratio = rng.uniform(0.7, 1.0)
        elif regime == 5:
            # Deep 3-layer snow (snowh > 0.25), some melting (t_air > 273.15)
            vt = land_veg_types[i % 12]
            st = (i % 12) + 1
            t_air = rng.uniform(265.0, 281.0)
            cz = rng.uniform(0.2, 0.8)
            sw = cz * rng.uniform(300.0, 800.0)
            sdepth = rng.uniform(0.26, 0.90)
            swe = sdepth * rng.uniform(200.0, 320.0)
            prcp = rng.uniform(0.0, 0.3)
            s_ratio = 0.0 if t_air > 275.0 else 0.5
        elif regime == 6:
            # Barren (16) or Urban (13) or Tundra (18, 19, 20)
            vt = [13, 16, 18, 19, 20][(i // 10) % 5]
            st = (i % 12) + 1
            t_air = rng.uniform(265.0, 305.0)
            cz = rng.uniform(-0.3, 0.9)
            sw = max(0.0, cz) * rng.uniform(200.0, 900.0)
            sdepth = rng.uniform(0.0, 0.15) if t_air < 273.0 else 0.0
            swe = sdepth * 200.0
            prcp = rng.uniform(0.0, 0.4)
            s_ratio = 1.0 if t_air < 272.0 else 0.0
        elif regime == 7:
            # Glacier / Permanent Snow & Ice (IVGTYP = 15, XICE = 0)
            vt = 15
            st = (i % 12) + 1
            t_air = rng.uniform(248.0, 278.0)
            cz = rng.uniform(-0.2, 0.7)
            sw = max(0.0, cz) * rng.uniform(200.0, 750.0)
            # Test both thin/moderate snow on glacier and deep snow on glacier
            sdepth = rng.uniform(0.02, 0.80)
            swe = sdepth * rng.uniform(200.0, 350.0)
            prcp = rng.uniform(0.0, 0.5)
            s_ratio = 1.0 if t_air < 273.15 else 0.2
        elif regime == 8:
            # Open water (XLAND = 2.0, IVGTYP = 17)
            vt = 17
            st = 14
            xland[i] = 2.0
            xice[i] = 0.0
            t_air = rng.uniform(275.0, 300.0)
            cz = rng.uniform(0.0, 0.8)
            sw = cz * 500.0
            sdepth = 0.0
            swe = 0.0
            prcp = 0.0
            s_ratio = 0.0
        else:
            # Sea-ice (XLAND = 1.0, XICE = 1.0, IVGTYP = 15)
            vt = 15
            st = 14
            xland[i] = 1.0
            xice[i] = 1.0
            t_air = rng.uniform(250.0, 272.0)
            cz = rng.uniform(0.0, 0.5)
            sw = cz * 300.0
            sdepth = 0.10
            swe = 25.0
            prcp = 0.0
            s_ratio = 1.0

        ivgtyp[i] = vt
        isltyp[i] = st
        xlat[i] = rng.uniform(-60.0, 65.0)
        xlong[i] = rng.uniform(-180.0, 180.0)
        coszen[i] = cz
        vf = rng.uniform(5.0, 90.0) if vt not in (15, 16, 17) else 0.0
        vegfra[i] = vf
        shdmax[i] = min(98.0, vf + rng.uniform(5.0, 25.0)) if vt not in (15, 16, 17) else 0.0
        tmn[i] = t_air + rng.uniform(-3.0, 3.0)

        swdown[i] = sw
        glw[i] = rng.uniform(210.0, 410.0)
        rainbl[i] = prcp
        sr[i] = s_ratio
        if prcp > 0.0:
            f_conv = rng.uniform(0.0, 0.3)
            f_shcv = rng.uniform(0.0, 0.1)
            f_nc = max(0.0, 1.0 - f_conv - f_shcv - rng.uniform(0.0, 0.05))
            mp_rainc[i] = prcp * f_conv
            mp_shcv[i] = prcp * f_shcv
            mp_rainnc[i] = prcp * f_nc
            if s_ratio > 0.0:
                mp_snow[i] = prcp * s_ratio * 0.7
                mp_graup[i] = prcp * s_ratio * 0.2
                mp_hail[i] = prcp * s_ratio * 0.1

        dz8w_1[i] = rng.uniform(35.0, 65.0)
        t_phy_1[i] = t_air
        qv_curr_1[i] = rng.uniform(0.001, 0.014)
        u_phy_1[i] = rng.uniform(-12.0, 12.0)
        v_phy_1[i] = rng.uniform(-12.0, 12.0)
        p_sfc = rng.uniform(85000.0, 102000.0)
        p8w_1[i] = p_sfc
        p8w_2[i] = p_sfc - dz8w_1[i] * rng.uniform(9.5, 11.5)

        tsk[i] = t_air + rng.uniform(-2.0, 2.0)
        snow[i] = swe
        snowh[i] = sdepth
        canwat[i] = rng.uniform(0.0, 0.5)
        lai[i] = rng.uniform(0.2, 5.0) if vt not in (15, 16, 17) else 0.0

        for k in range(nsoil):
            tslb[k, i] = t_air + rng.uniform(-4.0, 4.0)
            sm = rng.uniform(0.12, 0.38)
            smois[k, i] = sm
            sh2o[k, i] = sm

    with open(filepath, "wb") as f:
        f.write(struct.pack("=iiii", n_cols, nsoil, n_steps, yr))
        f.write(struct.pack("=dddd", dt, dx, julian_start, xice_threshold))
        f.write(dzs.tobytes())

        f.write(ivgtyp.tobytes())
        f.write(isltyp.tobytes())

        for arr in [
            xlat, xlong, coszen, vegfra, shdmax, tmn, xland, xice,
            swdown, glw, rainbl, sr, mp_rainc, mp_rainnc, mp_shcv, mp_snow, mp_graup, mp_hail,
            dz8w_1, t_phy_1, qv_curr_1, u_phy_1, v_phy_1, p8w_1, p8w_2,
            tsk, snow, snowh, canwat, lai,
        ]:
            f.write(arr.tobytes())

        for arr3 in [tslb, smois, sh2o]:
            f.write(arr3.tobytes())

    print(f"Generated {n_cols} Noah-MP test columns ({n_steps} steps) -> {filepath}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cols", type=int, default=2500)
    parser.add_argument("--steps", type=int, default=5)
    parser.add_argument("--out", type=str, required=True)
    args = parser.parse_args()
    generate_inputs(args.out, args.cols, args.steps)
