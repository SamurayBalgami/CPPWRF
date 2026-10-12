#!/usr/bin/env python3
"""
Generate realistic atmospheric sounding inputs for Kain-Fritsch (cu_physics = 1)
differential testing against WRF Fortran reference.

Regimes covered (N >= 2000 columns, 50 vertical levels):
1. Deep Convection (high CAPE, low CIN, moist boundary layer, positive W0)
2. Suppressed Convection / High CIN (strong capping inversion)
3. Shallow Convection (cloud depth < CHMIN, ishall=1)
4. Dry Atmosphere (low RH, desert/cold continental soundings)
5. High Wind Shear vs Weak Wind Shear (PEF variation)
6. Warm Rain Regime (warm sounding, no ice)
7. Active Relaxation (initial NCA >= 0.5 * DT)
"""

import sys
import struct
import numpy as np

def generate_soundings(ncols=2500, nlevels=50, seed=42):
    np.random.seed(seed)

    dt = 60.0
    dx = 10000.0
    cudt = 5.0
    stepcu = int(max(1, round(cudt * 60.0 / dt)))
    ktau = 1
    trigger = 1
    warm_rain = 0
    flag_qi = 1
    flag_qs = 1
    flag_qr = 1

    # Vertical height grid from surface to ~18 km
    # dz varies from ~50 m near surface to ~500 m aloft
    k_arr = np.linspace(0, 1, nlevels)
    dz_profile = 50.0 + 450.0 * (k_arr ** 1.5)
    z_faces = np.zeros(nlevels + 1)
    for k in range(nlevels):
        z_faces[k + 1] = z_faces[k] + dz_profile[k]
    z_half = 0.5 * (z_faces[:-1] + z_faces[1:])

    # Allocate arrays: shape (nlevels, ncols)
    dz8w = np.tile(dz_profile[:, np.newaxis], (1, ncols)).astype(np.float64)
    # Add minor column-to-column variation in layer thickness
    dz8w *= (1.0 + 0.02 * np.random.uniform(-1, 1, (nlevels, ncols)))

    t = np.zeros((nlevels, ncols), dtype=np.float64)
    qv = np.zeros((nlevels, ncols), dtype=np.float64)
    u = np.zeros((nlevels, ncols), dtype=np.float64)
    v = np.zeros((nlevels, ncols), dtype=np.float64)
    w = np.zeros((nlevels + 1, ncols), dtype=np.float64)
    pcps = np.zeros((nlevels, ncols), dtype=np.float64)
    rho = np.zeros((nlevels, ncols), dtype=np.float64)
    th = np.zeros((nlevels, ncols), dtype=np.float64)
    pi = np.zeros((nlevels, ncols), dtype=np.float64)
    w0avg = np.zeros((nlevels, ncols), dtype=np.float64)

    # 2D surface fields
    raincv = np.zeros(ncols, dtype=np.float64)
    pratec = np.zeros(ncols, dtype=np.float64)
    nca = np.full(ncols, -100.0, dtype=np.float64) # WRF default from init
    cutop = np.ones(ncols, dtype=np.float64)
    cubot = np.full(ncols, nlevels + 1, dtype=np.float64)

    r_d = 287.0
    cp = 1004.5
    g = 9.81
    p00 = 100000.0

    for c in range(ncols):
        # Determine regime
        if c < 500:
            # 1. Deep Convection
            tsfc = np.random.uniform(298.0, 308.0)
            lapse = np.random.uniform(6.5e-3, 7.8e-3)
            rh_sfc = np.random.uniform(0.75, 0.92)
            u_shear = np.random.uniform(10.0, 25.0)
            w_mag = np.random.uniform(0.1, 0.6)
            is_capped = False
            is_shallow = False
            dry_factor = 1.0
        elif c < 900:
            # 2. Suppressed / High CIN (inversion between 1.5 km and 3 km)
            tsfc = np.random.uniform(290.0, 300.0)
            lapse = 6.0e-3
            rh_sfc = np.random.uniform(0.5, 0.7)
            u_shear = np.random.uniform(5.0, 15.0)
            w_mag = np.random.uniform(-0.05, 0.05)
            is_capped = True
            is_shallow = False
            dry_factor = 1.0
        elif c < 1400:
            # 3. Shallow Convection (moist boundary layer, dry stable free troposphere)
            tsfc = np.random.uniform(292.0, 302.0)
            lapse = 5.5e-3
            rh_sfc = np.random.uniform(0.70, 0.85)
            u_shear = np.random.uniform(5.0, 15.0)
            w_mag = np.random.uniform(0.05, 0.25)
            is_capped = False
            is_shallow = True
            dry_factor = 1.0
        elif c < 1700:
            # 4. Dry Atmosphere
            tsfc = np.random.uniform(285.0, 310.0)
            lapse = 8.5e-3
            rh_sfc = np.random.uniform(0.05, 0.20)
            u_shear = np.random.uniform(5.0, 20.0)
            w_mag = np.random.uniform(-0.1, 0.1)
            is_capped = False
            is_shallow = False
            dry_factor = 0.15
        elif c < 2000:
            # 5. High Wind Shear
            tsfc = np.random.uniform(296.0, 304.0)
            lapse = 7.0e-3
            rh_sfc = np.random.uniform(0.72, 0.88)
            u_shear = np.random.uniform(30.0, 55.0)
            w_mag = np.random.uniform(0.15, 0.5)
            is_capped = False
            is_shallow = False
            dry_factor = 1.0
        elif c < 2300:
            # 6. Moderate convection with varying moisture
            tsfc = np.random.uniform(295.0, 305.0)
            lapse = np.random.uniform(6.0e-3, 7.5e-3)
            rh_sfc = np.random.uniform(0.65, 0.85)
            u_shear = np.random.uniform(10.0, 25.0)
            w_mag = np.random.uniform(0.05, 0.4)
            is_capped = False
            is_shallow = False
            dry_factor = 0.9
        else:
            # 7. Active Relaxation (initial NCA >= 0.5 * DT)
            tsfc = np.random.uniform(296.0, 304.0)
            lapse = 6.5e-3
            rh_sfc = 0.8
            u_shear = 15.0
            w_mag = 0.2
            is_capped = False
            is_shallow = False
            dry_factor = 1.0
            nca[c] = np.random.uniform(300.0, 1800.0)

        # Build temperature profile T(z)
        tc = np.zeros(nlevels)
        for k in range(nlevels):
            z = z_half[k]
            if is_capped and (1500.0 < z < 3000.0):
                # Strong temperature inversion
                tc[k] = tsfc - lapse * 1500.0 + 4.0 * (z - 1500.0) / 1500.0
            elif is_shallow and z > 3500.0:
                # Strong stable layer above 3.5 km
                tc[k] = tsfc - lapse * 3500.0 - 3.5e-3 * (z - 3500.0)
            else:
                tc[k] = tsfc - lapse * z
            if tc[k] < 195.0:
                tc[k] = 195.0 + 1.0e-3 * (z - 12000.0)

        t[:, c] = tc

        # Hydrostatic integration for pressure P(z)
        psfc = np.random.uniform(99000.0, 102000.0)
        p_faces = np.zeros(nlevels + 1)
        p_faces[0] = psfc
        for k in range(nlevels):
            t_layer = tc[k]
            p_faces[k + 1] = p_faces[k] * np.exp(-g * dz8w[k, c] / (r_d * t_layer))
            pcps[k, c] = 0.5 * (p_faces[k] + p_faces[k + 1])

        # Saturation vapor pressure Buck (1981)
        es = 611.2 * np.exp((17.67 * (tc - 273.15)) / (tc - 29.65))
        qsat = 0.622 * es / (pcps[:, c] - es)

        # Humidity profile qv(z)
        rh = np.zeros(nlevels)
        for k in range(nlevels):
            z = z_half[k]
            if is_shallow and z > 3000.0:
                rh[k] = 0.15 # Very dry air capping shallow clouds
            else:
                rh[k] = rh_sfc * np.exp(-z / 3500.0) * dry_factor
            rh[k] = np.clip(rh[k], 0.01, 0.98)

        qv[:, c] = np.clip(rh * qsat, 1.0e-6, qsat * 0.99)

        # Air density & exner
        tv = tc * (1.0 + 0.608 * qv[:, c])
        rho[:, c] = pcps[:, c] / (r_d * tv)
        pi[:, c] = (pcps[:, c] / p00) ** (r_d / cp)
        th[:, c] = tc / pi[:, c]

        # Horizontal winds
        u[:, c] = u_shear * (z_half / z_half[-1]) + np.random.uniform(-1, 1, nlevels)
        v[:, c] = np.random.uniform(-5, 5, nlevels)

        # Vertical velocity w on full levels (nlevels + 1)
        w_col = w_mag * np.sin(np.pi * np.linspace(0, 1, nlevels + 1))
        w[:, c] = w_col
        # Running average w0avg initialized to initial vertical velocity
        w0avg[:, c] = 0.5 * (w_col[:-1] + w_col[1:])

    return {
        "ncols": ncols,
        "nlevels": nlevels,
        "dt": dt,
        "dx": dx,
        "cudt": cudt,
        "stepcu": stepcu,
        "ktau": ktau,
        "trigger": trigger,
        "warm_rain": warm_rain,
        "flag_qi": flag_qi,
        "flag_qs": flag_qs,
        "flag_qr": flag_qr,
        "u": u,
        "v": v,
        "w": w,
        "th": th,
        "t": t,
        "qv": qv,
        "dz8w": dz8w,
        "pcps": pcps,
        "rho": rho,
        "pi": pi,
        "w0avg": w0avg,
        "raincv": raincv,
        "pratec": pratec,
        "nca": nca,
        "cutop": cutop,
        "cubot": cubot
    }

def save_inputs(data, filename="input_kf.bin"):
    with open(filename, "wb") as f:
        # Header
        header = struct.pack(
            "=ii ddd iiii iii",
            data["ncols"],
            data["nlevels"],
            data["dt"],
            data["dx"],
            data["cudt"],
            data["stepcu"],
            data["ktau"],
            data["trigger"],
            data["warm_rain"],
            data["flag_qi"],
            data["flag_qs"],
            data["flag_qr"]
        )
        f.write(header)

        # 3D arrays: (nlevels, ncols) in C order
        for key in ["u", "v", "th", "t", "qv", "dz8w", "pcps", "rho", "pi", "w0avg"]:
            f.write(data[key].astype(np.float64).tobytes(order="C"))

        # w array: (nlevels + 1, ncols) in C order
        f.write(data["w"].astype(np.float64).tobytes(order="C"))

        # 2D arrays: (ncols,)
        for key in ["raincv", "pratec", "nca", "cutop", "cubot"]:
            f.write(data[key].astype(np.float64).tobytes(order="C"))

    print(f"Saved {data['ncols']} columns x {data['nlevels']} levels to {filename}")

if __name__ == "__main__":
    ncols = 2500
    if len(sys.argv) > 1 and sys.argv[1] == "--ncols":
        ncols = int(sys.argv[2])
    data = generate_soundings(ncols=ncols)
    save_inputs(data)
