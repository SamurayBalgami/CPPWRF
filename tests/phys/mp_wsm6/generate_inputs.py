#!/usr/bin/env python3
"""
Input generator for WSM6 6-class microphysics differential testing.
Generates N >= 2000 atmospheric columns covering:
- Tropical soundings (warm rain, melting layers)
- Mid-latitude soundings (mixed-phase, snow, graupel)
- Arctic / cold soundings (pure ice, aggregation)
- Convective core (high concentrations of all 6 classes: qc, qr, qi, qs, qg, qv)
- Freezing rain / supercooled rain (T < 0 C)
- Melting snow and graupel (T > 0 C)
- Edge cases: all zero hydrometeors, autoconversion thresholds, tiny values near qmin / qcrmin.
"""

import argparse
import struct
import numpy as np

def generate_wsm6_inputs(n_cols=2500, nz=50, dt=120.0, seed=42):
    np.random.seed(seed)

    # Physical constants matching WRF
    r_d = 287.0
    r_v = 461.6
    cp = 7.0 * r_d / 2.0  # 1004.5
    cpv = 4.0 * r_v       # 1846.4 in WRF or standard cpv = 1870.0; WRF uses 4.*r_v in constants or cpv=1870.
    # In WRF module_mp_wsm6: cpv is passed in from driver. In module_model_constants.F: cpv = 4.*r_v = 1846.4
    # Driver sets cpv = cpv
    cpv = 4.0 * r_v
    g = 9.81
    t0c = 273.15
    ep1 = r_v / r_d - 1.0
    ep2 = r_d / r_v
    qmin = 1.0e-15
    xls = 2.85e6
    xlv0 = 2.50e6
    xlf0 = xls - xlv0
    den0 = 1.28
    denr = 1000.0
    dens = 100.0
    cliq = 4190.0
    cice = 2106.0
    psat = 610.78
    hail_opt = 0

    p1000mb = 100000.0

    # Vertical coordinates: z from 0 to 20 km
    k_faces = np.linspace(0, 1, nz + 1)
    zw = 20000.0 * (k_faces ** 1.3)
    dz = zw[1:] - zw[:-1]
    z_mid = 0.5 * (zw[:-1] + zw[1:])

    z = np.tile(z_mid[:, np.newaxis], (1, n_cols))
    delz = np.tile(dz[:, np.newaxis], (1, n_cols))

    # Pressure profile
    H_scale = 8000.0
    p = p1000mb * np.exp(-z / H_scale)
    pii = (p / p1000mb) ** (r_d / cp)

    T = np.zeros((nz, n_cols), dtype=np.float64)

    # 5 atmospheric regimes
    c1 = n_cols // 5
    c2 = 2 * n_cols // 5
    c3 = 3 * n_cols // 5
    c4 = 4 * n_cols // 5

    # 1. Tropical (warm surface ~300 K, deep troposphere)
    Ts_trop = np.random.uniform(298.0, 305.0, c1)
    for idx, i in enumerate(range(0, c1)):
        T[:, i] = np.maximum(Ts_trop[idx] - 0.006 * z[:, i], 195.0)

    # 2. Mid-latitude (surface ~285 K, melting level around 2-3 km)
    Ts_mid = np.random.uniform(280.0, 292.0, c2 - c1)
    for idx, i in enumerate(range(c1, c2)):
        T[:, i] = np.maximum(Ts_mid[idx] - 0.0065 * z[:, i], 200.0)

    # 3. Cold / Arctic (surface ~250 K, below freezing throughout)
    Ts_arc = np.random.uniform(235.0, 265.0, c3 - c2)
    for idx, i in enumerate(range(c2, c3)):
        T[:, i] = np.maximum(Ts_arc[idx] - 0.0055 * z[:, i], 200.0)

    # 4. Convective / Mixed-phase storm (surface ~295 K, strong lapse rate)
    Ts_conv = np.random.uniform(290.0, 300.0, c4 - c3)
    for idx, i in enumerate(range(c3, c4)):
        T[:, i] = np.maximum(Ts_conv[idx] - 0.0072 * z[:, i], 195.0)

    # 5. Freezing rain / Inversion / Edge cases
    Ts_edge = np.random.uniform(265.0, 295.0, n_cols - c4)
    for idx, i in enumerate(range(c4, n_cols)):
        # Include temperature inversions in some columns
        if idx % 2 == 0:
            T_base = Ts_edge[idx] - 0.0065 * z[:, i]
            # Warm nose aloft
            warm_nose = 8.0 * np.exp(-((z[:, i] - 2500.0) / 1000.0) ** 2)
            T[:, i] = np.maximum(T_base + warm_nose, 200.0)
        else:
            T[:, i] = np.maximum(Ts_edge[idx] - 0.0065 * z[:, i], 195.0)

    theta = T / pii

    # Vapor saturation mixing ratio (approximate Tetens)
    svp1 = 0.6112
    svp2 = 17.67
    svp3 = 29.65
    es = 1000.0 * svp1 * np.exp(svp2 * (T - t0c) / (T - svp3))
    qvs = ep2 * es / (p - es)
    qvs = np.maximum(qvs, 1.0e-7)

    # Initialize hydrometeor fields: shape (nz, n_cols)
    qv = np.zeros((nz, n_cols), dtype=np.float64)
    qc = np.zeros((nz, n_cols), dtype=np.float64)
    qr = np.zeros((nz, n_cols), dtype=np.float64)
    qi = np.zeros((nz, n_cols), dtype=np.float64)
    qs = np.zeros((nz, n_cols), dtype=np.float64)
    qg = np.zeros((nz, n_cols), dtype=np.float64)

    # Populate hydrometeors by regime
    for i in range(n_cols):
        rh_target = np.random.uniform(0.3, 1.05, nz)
        qv[:, i] = np.maximum(rh_target * qvs[:, i], 1.0e-7)

        if i < c1:
            # Tropical: warm clouds, rain below 4.5 km, ice/snow aloft
            cloud_mask = (z[:, i] > 800.0) & (z[:, i] < 6000.0)
            qc[cloud_mask, i] = np.random.uniform(0.1e-3, 1.5e-3, np.sum(cloud_mask))
            rain_mask = (z[:, i] < 4500.0) & (T[:, i] > t0c)
            qr[rain_mask, i] = np.random.uniform(0.05e-3, 3.0e-3, np.sum(rain_mask))
            cirrus_mask = (z[:, i] > 8000.0) & (T[:, i] < 240.0)
            qi[cirrus_mask, i] = np.random.uniform(0.01e-3, 0.2e-3, np.sum(cirrus_mask))

        elif i < c2:
            # Mid-latitude: mixed phase
            c_mask = (T[:, i] > 260.0) & (T[:, i] < 285.0)
            qc[c_mask, i] = np.random.uniform(0.05e-3, 0.8e-3, np.sum(c_mask))
            r_mask = (T[:, i] > t0c) & (z[:, i] < 3500.0)
            qr[r_mask, i] = np.random.uniform(0.05e-3, 1.5e-3, np.sum(r_mask))
            s_mask = (T[:, i] <= t0c) & (T[:, i] > 245.0)
            qs[s_mask, i] = np.random.uniform(0.1e-3, 2.0e-3, np.sum(s_mask))
            i_mask = (T[:, i] <= 245.0)
            qi[i_mask, i] = np.random.uniform(0.02e-3, 0.3e-3, np.sum(i_mask))
            g_mask = (T[:, i] < 275.0) & (T[:, i] > 255.0)
            qg[g_mask, i] = np.random.uniform(0.05e-3, 1.0e-3, np.sum(g_mask))

        elif i < c3:
            # Cold / Arctic: almost pure ice, snow
            i_mask = (T[:, i] < t0c) & (z[:, i] > 1000.0)
            qi[i_mask, i] = np.random.uniform(0.01e-3, 0.4e-3, np.sum(i_mask))
            s_mask = (z[:, i] < 6000.0)
            qs[s_mask, i] = np.random.uniform(0.05e-3, 1.5e-3, np.sum(s_mask))

        elif i < c4:
            # Convective storm: high loadings of all 6 classes
            c_mask = (z[:, i] > 1000.0) & (z[:, i] < 8000.0)
            qc[c_mask, i] = np.random.uniform(0.2e-3, 2.5e-3, np.sum(c_mask))
            r_mask = (z[:, i] < 5000.0)
            qr[r_mask, i] = np.random.uniform(0.2e-3, 5.0e-3, np.sum(r_mask))
            s_mask = (z[:, i] > 3000.0) & (z[:, i] < 11000.0)
            qs[s_mask, i] = np.random.uniform(0.2e-3, 4.0e-3, np.sum(s_mask))
            g_mask = (z[:, i] > 2000.0) & (z[:, i] < 9000.0)
            qg[g_mask, i] = np.random.uniform(0.2e-3, 6.0e-3, np.sum(g_mask))
            i_mask = (z[:, i] > 7000.0)
            qi[i_mask, i] = np.random.uniform(0.05e-3, 0.8e-3, np.sum(i_mask))

        else:
            # Edge cases
            sub = i - c4
            if sub % 6 == 0:
                # Completely dry, zero hydrometeors
                pass
            elif sub % 6 == 1:
                # Near autoconversion threshold qc ~ 0.42e-3
                qc[:, i] = 0.419e-3 * np.random.uniform(0.9, 1.1, nz)
            elif sub % 6 == 2:
                # Near snow aggregation threshold qs ~ 0.6e-3
                qs[:, i] = 0.6e-3 * np.random.uniform(0.9, 1.1, nz)
            elif sub % 6 == 3:
                # Freezing rain: rain above freezing layer falling into subfreezing surface layer
                low_cold = T[:, i] < t0c
                qr[low_cold, i] = 1.0e-3
            elif sub % 6 == 4:
                # Pure graupel falling through melting layer
                low_warm = T[:, i] > t0c
                qg[low_warm, i] = 2.0e-3
            elif sub % 6 == 5:
                # Homogeneous freezing boundary: T ~ -40 C (233.15 K)
                t_40_mask = np.abs(T[:, i] - (t0c - 40.0)) < 3.0
                qc[t_40_mask, i] = 1.0e-3

    # Density rho = p / (r_d * T * (1 + 0.608 * qv))
    den = p / (r_d * T * (1.0 + ep1 * qv))

    # 2D surface arrays
    rain = np.random.uniform(0.0, 10.0, n_cols).astype(np.float64)
    rainncv = np.zeros(n_cols, dtype=np.float64)
    snow = np.random.uniform(0.0, 5.0, n_cols).astype(np.float64)
    snowncv = np.zeros(n_cols, dtype=np.float64)
    graupel = np.random.uniform(0.0, 2.0, n_cols).astype(np.float64)
    graupelncv = np.zeros(n_cols, dtype=np.float64)
    sr = np.zeros(n_cols, dtype=np.float64)

    return {
        'n_cols': n_cols, 'nz': nz, 'delt': dt,
        'g': g, 'cpd': cp, 'cpv': cpv, 'rd': r_d, 'rv': r_v,
        't0c': t0c, 'ep1': ep1, 'ep2': ep2, 'qmin': qmin,
        'xls': xls, 'xlv0': xlv0, 'xlf0': xlf0,
        'den0': den0, 'denr': denr, 'dens': dens,
        'cliq': cliq, 'cice': cice, 'psat': psat, 'hail_opt': hail_opt,
        'th': theta, 'q': qv, 'qc': qc, 'qr': qr, 'qi': qi, 'qs': qs, 'qg': qg,
        'den': den, 'pii': pii, 'p': p, 'delz': delz,
        'rain': rain, 'rainncv': rainncv,
        'snow': snow, 'snowncv': snowncv,
        'graupel': graupel, 'graupelncv': graupelncv,
        'sr': sr
    }

def write_binary(data, filename):
    with open(filename, 'wb') as f:
        # Header
        f.write(struct.pack('ii', data['n_cols'], data['nz']))
        f.write(struct.pack('ddddddddddddddddddd',
                            data['delt'], data['g'], data['cpd'], data['cpv'],
                            data['rd'], data['rv'], data['t0c'], data['ep1'],
                            data['ep2'], data['qmin'], data['xls'], data['xlv0'],
                            data['xlf0'], data['den0'], data['denr'], data['dens'],
                            data['cliq'], data['cice'], data['psat']))
        f.write(struct.pack('i', data['hail_opt']))

        # 3D arrays: shape (nz, n_cols)
        for key in ['th', 'q', 'qc', 'qr', 'qi', 'qs', 'qg', 'den', 'pii', 'p', 'delz']:
            f.write(data[key].astype(np.float64).tobytes())

        # 2D arrays: shape (n_cols,)
        for key in ['rain', 'rainncv', 'snow', 'snowncv', 'graupel', 'graupelncv', 'sr']:
            f.write(data[key].astype(np.float64).tobytes())

def main():
    parser = argparse.ArgumentParser(description='WSM6 input generator')
    parser.add_argument('--cols', type=int, default=2500, help='Number of atmospheric columns')
    parser.add_argument('--nz', type=int, default=50, help='Number of vertical levels')
    parser.add_argument('--dt', type=float, default=120.0, help='Timestep in seconds')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    parser.add_argument('--out', type=str, default='wsm6_input.bin', help='Output binary file')
    args = parser.parse_args()

    data = generate_wsm6_inputs(n_cols=args.cols, nz=args.nz, dt=args.dt, seed=args.seed)
    write_binary(data, args.out)
    print(f"Generated {args.cols} columns x {args.nz} levels WSM6 inputs to {args.out}")

if __name__ == '__main__':
    main()
