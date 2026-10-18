#!/usr/bin/env python3
"""
Generate realistic atmospheric sounding inputs for Thompson microphysics (mp_physics = 8)
differential testing against WRF Fortran reference.

Regimes covered (N >= 2000 columns, 50 vertical levels):
1. Warm Rain (T > 273 K, cloud droplets and rain, accretion, autoconversion)
2. Cold Convective / Mixed Phase (supercooled droplets, ice nucleation, graupel, snow)
3. Pristine Ice & Cirrus (cold upper levels, ice autoconversion, deposition/sublimation)
4. Snow Riming & Aggregation (subfreezing, snow collecting liquid & ice)
5. Graupel Riming & Collection (graupel-rain collection lookup tables)
6. Subsaturated Evaporation (dry subcloud layer, rain & snow evaporation/sublimation)
7. Sub-threshold Background (trace hydrometeors near/below R1 floor)
"""

import sys
import struct
import numpy as np

def generate_soundings(ncols=2500, nlevels=50, seed=42):
    np.random.seed(seed)

    dt = 60.0
    R_d = 287.04
    Cp = 1004.0
    g = 9.81
    p0 = 100000.0
    rcp = R_d / Cp

    # Vertical height grid from surface to ~16 km
    k_arr = np.linspace(0, 1, nlevels)
    dz_profile = 50.0 + 450.0 * (k_arr ** 1.3)
    z_faces = np.zeros(nlevels + 1)
    for k in range(nlevels):
        z_faces[k + 1] = z_faces[k] + dz_profile[k]
    z_mid = 0.5 * (z_faces[:-1] + z_faces[1:])

    # Base arrays: shape (nlevels, ncols)
    dz = np.tile(dz_profile[:, np.newaxis], (1, ncols)).astype(np.float64)
    dz *= (1.0 + 0.01 * np.random.uniform(-1, 1, (nlevels, ncols)))

    p = np.zeros((nlevels, ncols), dtype=np.float64)
    th = np.zeros((nlevels, ncols), dtype=np.float64)
    t = np.zeros((nlevels, ncols), dtype=np.float64)
    qv = np.zeros((nlevels, ncols), dtype=np.float64)
    qc = np.zeros((nlevels, ncols), dtype=np.float64)
    qr = np.zeros((nlevels, ncols), dtype=np.float64)
    qi = np.zeros((nlevels, ncols), dtype=np.float64)
    qs = np.zeros((nlevels, ncols), dtype=np.float64)
    qg = np.zeros((nlevels, ncols), dtype=np.float64)
    ni = np.zeros((nlevels, ncols), dtype=np.float64)
    nr = np.zeros((nlevels, ncols), dtype=np.float64)
    w = np.zeros((nlevels, ncols), dtype=np.float64)
    pii = np.zeros((nlevels, ncols), dtype=np.float64)

    # Surface 2D arrays
    rainnc = np.zeros(ncols, dtype=np.float64)
    rainncv = np.zeros(ncols, dtype=np.float64)
    snownc = np.zeros(ncols, dtype=np.float64)
    snowncv = np.zeros(ncols, dtype=np.float64)
    graupelnc = np.zeros(ncols, dtype=np.float64)
    graupelncv = np.zeros(ncols, dtype=np.float64)
    sr = np.zeros(ncols, dtype=np.float64)

    # Saturation vapor pressure formula (Tetens)
    def esat(tk):
        # Bolton 1980 / Tetens
        return 611.2 * np.exp(17.67 * (tk - 273.15) / (tk - 29.65))

    for col in range(ncols):
        # Surface pressure between 980 and 1020 hPa
        psfc = 101325.0 + np.random.uniform(-2000, 2000)

        if col < 500:
            # Regime 1: Warm Rain (entire troposphere warm or deep warm layer)
            tsfc = np.random.uniform(293.0, 305.0)
            gamma = np.random.uniform(0.0050, 0.0065)
            t_col = tsfc - gamma * z_mid
            # Ensure warm in lower 4 km
            t_col = np.maximum(t_col, 274.0 - 0.003 * z_mid)
            rh = np.random.uniform(0.70, 0.98, nlevels)
            w_val = np.random.uniform(-0.2, 1.5, nlevels)

            # Cloud and rain in warm layer
            warm_mask = (t_col > 275.0) & (z_mid < 4500.0)
            qc[warm_mask, col] = np.random.uniform(0.1e-3, 1.5e-3, np.sum(warm_mask))
            qr[warm_mask, col] = np.random.uniform(0.05e-3, 1.2e-3, np.sum(warm_mask))
            # Number concentration of rain: ~10^2 to 10^5 /kg
            nr[warm_mask, col] = np.random.uniform(500.0, 50000.0, np.sum(warm_mask))

        elif col < 1000:
            # Regime 2: Cold Convection / Mixed Phase
            tsfc = np.random.uniform(288.0, 298.0)
            gamma = np.random.uniform(0.0065, 0.0080)
            t_col = tsfc - gamma * z_mid
            rh = np.random.uniform(0.75, 0.99, nlevels)
            w_val = np.random.uniform(0.5, 8.0, nlevels) * np.sin(np.pi * z_mid / 14000.0)

            # Mixed-phase layer: 235 K < T < 273 K
            mixed_mask = (t_col >= 235.0) & (t_col <= 273.15)
            qc[mixed_mask, col] = np.random.uniform(0.05e-3, 0.8e-3, np.sum(mixed_mask))
            qi[mixed_mask, col] = np.random.uniform(0.01e-3, 0.3e-3, np.sum(mixed_mask))
            ni[mixed_mask, col] = np.random.uniform(1.0e3, 5.0e5, np.sum(mixed_mask))
            qs[mixed_mask, col] = np.random.uniform(0.05e-3, 1.0e-3, np.sum(mixed_mask))
            qg[mixed_mask, col] = np.random.uniform(0.1e-3, 2.0e-3, np.sum(mixed_mask))

            # Warm layer below 273 K: melting snow/graupel and rain
            rain_mask = (t_col > 273.15) & (z_mid < 3500.0)
            qr[rain_mask, col] = np.random.uniform(0.1e-3, 2.0e-3, np.sum(rain_mask))
            nr[rain_mask, col] = np.random.uniform(1.0e3, 1.0e5, np.sum(rain_mask))
            qs[rain_mask, col] = np.random.uniform(0.0, 0.2e-3, np.sum(rain_mask))
            qg[rain_mask, col] = np.random.uniform(0.0, 0.5e-3, np.sum(rain_mask))

        elif col < 1400:
            # Regime 3: Pristine Ice & Cirrus
            tsfc = np.random.uniform(270.0, 285.0)
            gamma = 0.0070
            t_col = tsfc - gamma * z_mid
            rh = np.random.uniform(0.30, 0.85, nlevels)
            w_val = np.random.uniform(-0.1, 0.5, nlevels)

            # High-altitude cold cirrus: T < 240 K
            cirrus_mask = (t_col < 240.0) & (z_mid > 6000.0)
            qi[cirrus_mask, col] = np.random.uniform(0.005e-3, 0.25e-3, np.sum(cirrus_mask))
            ni[cirrus_mask, col] = np.random.uniform(1.0e4, 2.0e6, np.sum(cirrus_mask))
            qs[cirrus_mask, col] = np.random.uniform(0.001e-3, 0.05e-3, np.sum(cirrus_mask))

        elif col < 1800:
            # Regime 4: Snow Riming & Aggregation
            tsfc = np.random.uniform(265.0, 275.0)
            gamma = 0.0060
            t_col = tsfc - gamma * z_mid
            rh = np.random.uniform(0.80, 0.98, nlevels)
            w_val = np.random.uniform(-0.2, 1.5, nlevels)

            # Snow-rich layer
            snow_mask = (t_col < 273.15) & (z_mid < 8000.0)
            qs[snow_mask, col] = np.random.uniform(0.1e-3, 1.8e-3, np.sum(snow_mask))
            # Supercooled droplets for riming
            qc[snow_mask & (t_col > 255.0), col] = np.random.uniform(0.02e-3, 0.4e-3, np.sum(snow_mask & (t_col > 255.0)))
            # Ice crystals for aggregation
            qi[snow_mask, col] = np.random.uniform(0.005e-3, 0.1e-3, np.sum(snow_mask))
            ni[snow_mask, col] = np.random.uniform(1.0e3, 2.0e5, np.sum(snow_mask))

        elif col < 2100:
            # Regime 5: Graupel Riming & Collection (active graupel-rain collision)
            tsfc = np.random.uniform(285.0, 295.0)
            gamma = 0.0075
            t_col = tsfc - gamma * z_mid
            rh = np.random.uniform(0.80, 0.98, nlevels)
            w_val = np.random.uniform(1.0, 6.0, nlevels)

            # Coexisting graupel and rain near melting level (265 K to 278 K)
            coll_mask = (t_col > 260.0) & (t_col < 280.0)
            qg[coll_mask, col] = np.random.uniform(0.2e-3, 2.5e-3, np.sum(coll_mask))
            qr[coll_mask, col] = np.random.uniform(0.1e-3, 1.5e-3, np.sum(coll_mask))
            nr[coll_mask, col] = np.random.uniform(1.0e3, 8.0e4, np.sum(coll_mask))
            qc[coll_mask & (t_col > 265.0), col] = np.random.uniform(0.05e-3, 0.5e-3, np.sum(coll_mask & (t_col > 265.0)))

        elif col < 2350:
            # Regime 6: Subsaturated Evaporation & Rain Sublimation
            tsfc = np.random.uniform(290.0, 305.0)
            gamma = 0.0085
            t_col = tsfc - gamma * z_mid
            # Dry subcloud layer
            rh = np.zeros(nlevels)
            rh[z_mid < 3000.0] = np.random.uniform(0.15, 0.45, np.sum(z_mid < 3000.0))
            rh[z_mid >= 3000.0] = np.random.uniform(0.60, 0.85, np.sum(z_mid >= 3000.0))
            w_val = np.random.uniform(-1.0, 0.2, nlevels)

            # Rain falling through dry subcloud layer
            sub_mask = (z_mid < 2500.0)
            qr[sub_mask, col] = np.random.uniform(0.1e-3, 1.0e-3, np.sum(sub_mask))
            nr[sub_mask, col] = np.random.uniform(500.0, 30000.0, np.sum(sub_mask))
            # Snow aloft
            snow_aloft = (z_mid >= 3000.0) & (t_col < 273.15)
            qs[snow_aloft, col] = np.random.uniform(0.05e-3, 0.5e-3, np.sum(snow_aloft))

        else:
            # Regime 7: Sub-threshold Background / Clear Sky
            tsfc = np.random.uniform(280.0, 295.0)
            gamma = 0.0065
            t_col = tsfc - gamma * z_mid
            rh = np.random.uniform(0.20, 0.70, nlevels)
            w_val = np.random.uniform(-0.1, 0.1, nlevels)

            # Trace amounts near R1 = 1e-12
            qc[:, col] = np.random.uniform(0.0, 5.0e-13, nlevels)
            qi[:, col] = np.random.uniform(0.0, 5.0e-13, nlevels)
            qr[:, col] = np.random.uniform(0.0, 5.0e-13, nlevels)
            qs[:, col] = np.random.uniform(0.0, 5.0e-13, nlevels)
            qg[:, col] = np.random.uniform(0.0, 5.0e-13, nlevels)

        # Hydrostatic pressure profile
        p_curr = psfc
        for k in range(nlevels):
            dz_k = dz[k, col]
            tk = t_col[k]
            # Hydrostatic equation: dp = - rho * g * dz = - (p / (R*T)) * g * dz
            p_mid = p_curr * np.exp(-g * dz_k * 0.5 / (R_d * tk))
            p[k, col] = p_mid
            p_curr = p_curr * np.exp(-g * dz_k / (R_d * tk))

            # Exner function
            pii_val = (p_mid / p0) ** rcp
            pii[k, col] = pii_val
            t[k, col] = tk
            th[k, col] = tk / pii_val

            # Water vapor qv from RH
            es = esat(tk)
            qvs = 0.622 * es / max(100.0, p_mid - es)
            qv[k, col] = max(1.0e-7, min(0.035, rh[k] * qvs))

        w[:, col] = w_val

    data = {
        "ncols": ncols,
        "nlevels": nlevels,
        "dt": dt,
        "th": th,
        "qv": qv,
        "qc": qc,
        "qr": qr,
        "qi": qi,
        "qs": qs,
        "qg": qg,
        "ni": ni,
        "nr": nr,
        "pii": pii,
        "p": p,
        "w": w,
        "dz": dz,
        "rainnc": rainnc,
        "rainncv": rainncv,
        "snownc": snownc,
        "snowncv": snowncv,
        "graupelnc": graupelnc,
        "graupelncv": graupelncv,
        "sr": sr
    }
    return data

def save_inputs(data, filename="input_thompson.bin"):
    with open(filename, "wb") as f:
        # Header: ncols (i32), nlevels (i32), dt (f64)
        header = struct.pack("=ii d", data["ncols"], data["nlevels"], data["dt"])
        f.write(header)

        # 3D arrays: (nlevels, ncols) in C order
        for key in ["th", "qv", "qc", "qr", "qi", "qs", "qg", "ni", "nr", "pii", "p", "w", "dz"]:
            f.write(data[key].astype(np.float64).tobytes(order="C"))

        # 2D arrays: (ncols,)
        for key in ["rainnc", "rainncv", "snownc", "snowncv", "graupelnc", "graupelncv", "sr"]:
            f.write(data[key].astype(np.float64).tobytes(order="C"))

    print(f"Saved {data['ncols']} columns x {data['nlevels']} levels to {filename}")

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncols", type=int, default=2500)
    ap.add_argument("--nlevels", type=int, default=50)
    args = ap.parse_args()
    data = generate_soundings(ncols=args.ncols, nlevels=args.nlevels)
    save_inputs(data)
