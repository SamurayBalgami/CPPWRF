#!/usr/bin/env python3
"""
Input generator for Kessler microphysics differential testing.
Generates N >= 2000 realistic atmospheric columns covering:
- Tropical soundings
- Mid-latitude soundings
- Dry / Arctic soundings
- Saturated layers (supersaturation => condensation)
- Supercooled layers (T < 273.15 K)
- Hydrometeors present / absent
- Exactly zero and near-zero mixing ratios
- Autoconversion threshold edge cases (qc ~ 0.001)
- High rain mixing ratios triggering Courant substepping (crmax > 0.75)
- Small dt and large dt support
"""

import argparse
import struct
import numpy as np

def generate_soundings(n_cols=2500, nz=50, dt=12.0, seed=42):
    np.random.seed(seed)

    # Physical constants matching WRF
    r_d = 287.0
    r_v = 461.6
    cp = 7.0 * r_d / 2.0  # 1004.5
    xlv = 2.5e6
    ep2 = r_d / r_v
    svp1 = 0.6112
    svp2 = 17.67
    svp3 = 29.65
    svpt0 = 273.15
    rhowater = 1000.0
    p1000mb = 100000.0

    # Vertical coordinates: z from 0 to 20 km with stretching
    k_faces = np.linspace(0, 1, nz + 1)
    zw = 20000.0 * (k_faces ** 1.3)  # Interface levels
    dz = zw[1:] - zw[:-1]            # nz layer thicknesses
    z_mid = 0.5 * (zw[:-1] + zw[1:]) # nz layer midpoints

    # 3D arrays: shape (nz, n_cols)
    z = np.tile(z_mid[:, np.newaxis], (1, n_cols))
    dz8w = np.tile(dz[:, np.newaxis], (1, n_cols))

    # Hydrostatic pressure profile: p ~ p0 * exp(-z / H)
    H_scale = 8000.0
    p = p1000mb * np.exp(-z / H_scale)
    pii = (p / p1000mb) ** (r_d / cp)

    # Temperature profile T (K) per column
    T = np.zeros((nz, n_cols), dtype=np.float64)

    # Partition columns dynamically into 5 regimes
    c1 = n_cols // 5
    c2 = 2 * n_cols // 5
    c3 = 3 * n_cols // 5
    c4 = 4 * n_cols // 5
    n_trop = c1
    n_mid = c2 - c1
    n_arc = c3 - c2
    n_conv = c4 - c3
    n_edge = n_cols - c4

    # 1. Tropical (0 .. c1)
    Ts_trop = np.random.uniform(298.0, 305.0, n_trop)
    for idx, i in enumerate(range(0, c1)):
        T_col = Ts_trop[idx] - 0.006 * z[:, i]
        T[:, i] = np.maximum(T_col, 200.0)

    # 2. Mid-latitude (c1 .. c2)
    Ts_mid = np.random.uniform(280.0, 295.0, n_mid)
    for idx, i in enumerate(range(c1, c2)):
        T_col = Ts_mid[idx] - 0.0065 * z[:, i]
        T[:, i] = np.maximum(T_col, 215.0)

    # 3. Dry / Arctic (c2 .. c3)
    Ts_arc = np.random.uniform(240.0, 270.0, n_arc)
    for idx, i in enumerate(range(c2, c3)):
        T_col = Ts_arc[idx] - 0.0055 * z[:, i]
        T[:, i] = np.maximum(T_col, 210.0)

    # 4. Supercooled / deep convective (c3 .. c4)
    Ts_conv = np.random.uniform(275.0, 292.0, n_conv)
    for idx, i in enumerate(range(c3, c4)):
        T_col = Ts_conv[idx] - 0.007 * z[:, i]
        T[:, i] = np.maximum(T_col, 205.0)

    # 5. Edge cases & varied conditions (c4 .. n_cols)
    Ts_edge = np.random.uniform(260.0, 300.0, n_edge)
    for idx, i in enumerate(range(c4, n_cols)):
        T_col = Ts_edge[idx] - 0.0065 * z[:, i]
        T[:, i] = np.maximum(T_col, 200.0)

    # Potential temperature theta = T / pii
    theta = T / pii

    # Air density rho = p / (r_d * T)
    rho = p / (r_d * T)

    # Saturation vapor pressure es and saturation mixing ratio qvs
    es = 1000.0 * svp1 * np.exp(svp2 * (T - svpt0) / (T - svp3))
    qvs = ep2 * es / (p - es)

    # Moisture initialization
    qv = np.zeros((nz, n_cols), dtype=np.float64)
    qc = np.zeros((nz, n_cols), dtype=np.float64)
    qr = np.zeros((nz, n_cols), dtype=np.float64)

    # 1. Tropical (0 .. c1)
    rh_trop = np.random.uniform(0.6, 0.95, (nz, n_trop))
    n_sat = min(200, n_trop // 2)
    rh_trop[5:20, :n_sat] = np.random.uniform(1.01, 1.08, (15, n_sat))
    qv[:, :c1] = rh_trop * qvs[:, :c1]
    qc[5:20, :n_sat] = np.random.uniform(0.0001, 0.0025, (15, n_sat))
    n_rain = min(150, n_sat)
    qr[0:15, :n_rain] = np.random.uniform(0.0005, 0.004, (15, n_rain))

    # 2. Mid-latitude (c1 .. c2)
    rh_mid = np.random.uniform(0.3, 0.85, (nz, n_mid))
    freezing_k = np.argmin(np.abs(T[:, c1:c2] - 273.15), axis=0)
    for idx in range(n_mid):
        col_idx = c1 + idx
        fk = freezing_k[idx]
        k_lo = max(0, fk - 4)
        k_hi = min(nz, fk + 8)
        rh_mid[k_lo:k_hi, idx] = np.random.uniform(1.0, 1.05)
        qc[k_lo:k_hi, col_idx] = np.random.uniform(0.0002, 0.0020)
        qr[:k_lo, col_idx] = np.random.uniform(0.0001, 0.002)
    qv[:, c1:c2] = rh_mid * qvs[:, c1:c2]

    # 3. Dry / Arctic (c2 .. c3)
    rh_arc = np.random.uniform(0.001, 0.20, (nz, n_arc))
    rh_arc[20:, :] = 0.0001
    qv[:, c2:c3] = rh_arc * qvs[:, c2:c3]
    n_evap = min(100, n_arc // 4)
    qr[10:25, c2:c2+n_evap] = np.random.uniform(0.0005, 0.003, (15, n_evap))

    # 4. Supercooled & heavy rain (c3 .. c4)
    rh_conv = np.random.uniform(0.5, 1.02, (nz, n_conv))
    qv[:, c3:c4] = rh_conv * qvs[:, c3:c4]
    n_hrain = min(250, n_conv // 2)
    qr[5:35, c3:c3+n_hrain] = np.random.uniform(0.004, 0.012, (30, n_hrain))
    qc[15:35, c3:c3+n_hrain] = np.random.uniform(0.001, 0.003, (20, n_hrain))

    # 5. Boundary edge cases (c4 .. n_cols)
    s0 = c4
    ds = n_edge // 5
    # (a) Zero hydrometeors
    qv[:, s0:s0+ds] = 0.5 * qvs[:, s0:s0+ds]
    # (b) Near-zero 1e-15
    qv[:, s0+ds:s0+2*ds] = 1e-15
    qc[:, s0+ds:s0+2*ds] = 1e-15
    qr[:, s0+ds:s0+2*ds] = 1e-15
    # (c) Autoconversion threshold qc ~ 0.001
    qv[:, s0+2*ds:s0+3*ds] = 0.8 * qvs[:, s0+2*ds:s0+3*ds]
    qc[:, s0+2*ds:s0+2*ds+ds//2] = 0.001
    qc[:, s0+2*ds+ds//2:s0+3*ds] = 0.00100001
    # (d) Extreme rain
    qr[10:30, s0+3*ds:s0+4*ds] = 0.015
    qv[:, s0+3*ds:s0+4*ds] = 0.9 * qvs[:, s0+3*ds:s0+4*ds]
    # (e) Mixed random
    qv[:, s0+4*ds:] = np.random.uniform(0.0, 1.1, (nz, n_cols - (s0+4*ds))) * qvs[:, s0+4*ds:]
    qc[:, s0+4*ds:] = np.random.uniform(0.0, 0.003, (nz, n_cols - (s0+4*ds)))
    qr[:, s0+4*ds:] = np.random.uniform(0.0, 0.005, (nz, n_cols - (s0+4*ds)))

    # 2D surface rain
    rainnc = np.random.uniform(0.0, 10.0, n_cols)
    rainncv = np.zeros(n_cols, dtype=np.float64)

    return {
        'n_cols': n_cols,
        'nz': nz,
        'dt': dt,
        'xlv': xlv,
        'cp': cp,
        'ep2': ep2,
        'svp1': svp1,
        'svp2': svp2,
        'svp3': svp3,
        'svpt0': svpt0,
        'rhowater': rhowater,
        'theta': theta,
        'qv': qv,
        'qc': qc,
        'qr': qr,
        'rho': rho,
        'pii': pii,
        'dz8w': dz8w,
        'z': z,
        'rainnc': rainnc,
        'rainncv': rainncv
    }

def save_binary(data, filepath):
    with open(filepath, 'wb') as f:
        # Header: n_cols, nz (int32), dt (float64)
        f.write(struct.pack('=iid', data['n_cols'], data['nz'], data['dt']))
        # Constants: 8 doubles
        f.write(struct.pack('=8d',
            data['xlv'], data['cp'], data['ep2'], data['svp1'],
            data['svp2'], data['svp3'], data['svpt0'], data['rhowater']))
        # 3D arrays: (nz, n_cols) contiguous in C order
        data['theta'].astype(np.float64).tofile(f)
        data['qv'].astype(np.float64).tofile(f)
        data['qc'].astype(np.float64).tofile(f)
        data['qr'].astype(np.float64).tofile(f)
        data['rho'].astype(np.float64).tofile(f)
        data['pii'].astype(np.float64).tofile(f)
        data['dz8w'].astype(np.float64).tofile(f)
        data['z'].astype(np.float64).tofile(f)
        # 2D arrays: n_cols doubles
        data['rainnc'].astype(np.float64).tofile(f)
        data['rainncv'].astype(np.float64).tofile(f)

def main():
    parser = argparse.ArgumentParser(description='Generate Kessler test inputs')
    parser.add_argument('--output', '-o', default='kessler_input.bin', help='Output binary path')
    parser.add_argument('--ncols', '-n', type=int, default=2500, help='Number of columns (>= 2000)')
    parser.add_argument('--nz', type=int, default=50, help='Number of vertical levels')
    parser.add_argument('--dt', type=float, default=12.0, help='Time step in seconds')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    args = parser.parse_args()

    print(f"Generating {args.ncols} columns x {args.nz} levels with dt={args.dt}s (seed={args.seed})...")
    data = generate_soundings(args.ncols, args.nz, args.dt, args.seed)
    save_binary(data, args.output)
    print(f"Saved test inputs to {args.output}")

if __name__ == '__main__':
    main()
