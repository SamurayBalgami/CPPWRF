#!/usr/bin/env python3
"""
Input generator for RRTMG Longwave (ra_lw_physics=4) differential verification.
Generates realistic, physically diverse atmospheric columns covering tropical,
midlatitude, arctic, supercooled, deep convective, and clear/cloudy edge cases.
"""
import argparse
import struct
import numpy as np

def generate_rrtmg_lw_inputs(num_cols=2500, nz=50, p_top=5000.0, out_path="input_rrtmg_lw.bin", seed=42):
    np.random.seed(seed)

    # Physical constants (WRF module_model_constants.F)
    g = 9.81
    rd = 287.0
    rv = 461.6
    cp = 1004.5
    rcp = rd / cp
    p1000mb = 100000.0

    icloud = 1
    cldovrlp = 2  # Maximum-random overlap (standard WRF)
    idcor = 0
    mp_physics = 6

    # 2D surface variables: (num_cols,)
    tsk = np.zeros(num_cols, dtype=np.float64)
    emiss = np.zeros(num_cols, dtype=np.float64)
    xland = np.zeros(num_cols, dtype=np.float64)
    xice = np.zeros(num_cols, dtype=np.float64)
    snow = np.zeros(num_cols, dtype=np.float64)
    xlat = np.zeros(num_cols, dtype=np.float64)

    # 3D layer variables: (nz, num_cols)
    p3d = np.zeros((nz, num_cols), dtype=np.float64)
    pi3d = np.zeros((nz, num_cols), dtype=np.float64)
    dz8w = np.zeros((nz, num_cols), dtype=np.float64)
    t3d = np.zeros((nz, num_cols), dtype=np.float64)
    rho3d = np.zeros((nz, num_cols), dtype=np.float64)
    qv3d = np.zeros((nz, num_cols), dtype=np.float64)
    qc3d = np.zeros((nz, num_cols), dtype=np.float64)
    qr3d = np.zeros((nz, num_cols), dtype=np.float64)
    qi3d = np.zeros((nz, num_cols), dtype=np.float64)
    qs3d = np.zeros((nz, num_cols), dtype=np.float64)
    qg3d = np.zeros((nz, num_cols), dtype=np.float64)
    cldfra3d = np.zeros((nz, num_cols), dtype=np.float64)

    # 3D interface variables: (nz+1, num_cols)
    p8w = np.zeros((nz + 1, num_cols), dtype=np.float64)
    t8w = np.zeros((nz + 1, num_cols), dtype=np.float64)

    # Standard sigma-like vertical coordinate distribution
    # Normalized eta from 1.0 (surface) to 0.0 (top)
    eta = np.linspace(1.0, 0.0, nz + 1)

    for i in range(num_cols):
        # Determine regime
        regime = i % 5  # 0: Tropical, 1: Mid-lat, 2: Dry/Arctic, 3: Deep convection, 4: Edge cases

        if regime == 0:
            # Tropical: warm surface, high moisture, marine boundary layer clouds
            p_sfc = np.random.uniform(100500.0, 101800.0)
            tsk[i] = np.random.uniform(298.0, 304.0)
            emiss[i] = np.random.uniform(0.96, 0.99)
            xland[i] = 2.0 if np.random.rand() > 0.3 else 1.0
            xice[i] = 0.0
            snow[i] = 0.0
            xlat[i] = np.random.uniform(-20.0, 20.0)
            t_sfc = tsk[i] - np.random.uniform(0.0, 2.0)
            lapse_rate = np.random.uniform(0.0055, 0.0065)
            q_sfc = np.random.uniform(0.015, 0.022)

        elif regime == 1:
            # Mid-latitude: moderate temps, mixed phase, mid/high clouds
            p_sfc = np.random.uniform(98000.0, 102500.0)
            tsk[i] = np.random.uniform(275.0, 295.0)
            emiss[i] = np.random.uniform(0.92, 0.98)
            xland[i] = 1.0 if np.random.rand() > 0.4 else 2.0
            xice[i] = 0.0
            snow[i] = np.random.uniform(0.0, 50.0) if tsk[i] < 280.0 else 0.0
            xlat[i] = np.random.uniform(30.0, 60.0) * (1.0 if np.random.rand() > 0.5 else -1.0)
            t_sfc = tsk[i] + np.random.uniform(-3.0, 3.0)
            lapse_rate = np.random.uniform(0.006, 0.0075)
            q_sfc = np.random.uniform(0.005, 0.012)

        elif regime == 2:
            # Dry / Arctic: cold surface, strong inversions, low water vapor, ice clouds
            p_sfc = np.random.uniform(97000.0, 103000.0)
            tsk[i] = np.random.uniform(230.0, 265.0)
            emiss[i] = np.random.uniform(0.95, 0.99)
            xland[i] = 1.0 if np.random.rand() > 0.5 else 2.0
            xice[i] = 1.0 if (xland[i] == 2.0 and tsk[i] < 271.0) else 0.0
            snow[i] = np.random.uniform(50.0, 500.0)
            xlat[i] = np.random.uniform(65.0, 85.0) * (1.0 if np.random.rand() > 0.5 else -1.0)
            t_sfc = tsk[i] + np.random.uniform(-5.0, 2.0)
            lapse_rate = np.random.uniform(0.003, 0.005)
            q_sfc = np.random.uniform(0.0001, 0.002)

        elif regime == 3:
            # Deep convection / Storm: thick high clouds, heavy precipitation
            p_sfc = np.random.uniform(99500.0, 101500.0)
            tsk[i] = np.random.uniform(295.0, 305.0)
            emiss[i] = np.random.uniform(0.94, 0.98)
            xland[i] = 1.0
            xice[i] = 0.0
            snow[i] = 0.0
            xlat[i] = np.random.uniform(-35.0, 35.0)
            t_sfc = tsk[i] - 1.0
            lapse_rate = np.random.uniform(0.0065, 0.008)
            q_sfc = np.random.uniform(0.014, 0.020)

        else:
            # Edge cases: completely clear sky (cldfra=0), fully overcast (cldfra=1),
            # very low humidity (1e-6), isothermal profile, or high tropopause
            p_sfc = np.random.uniform(95000.0, 103500.0)
            tsk[i] = np.random.uniform(240.0, 315.0)
            emiss[i] = np.random.uniform(0.85, 1.0)
            xland[i] = 1.0 if np.random.rand() > 0.5 else 2.0
            xice[i] = 1.0 if (tsk[i] < 270.0 and xland[i] == 2.0) else 0.0
            snow[i] = np.random.uniform(0.0, 1000.0)
            xlat[i] = np.random.uniform(-80.0, 80.0)
            t_sfc = tsk[i]
            lapse_rate = np.random.uniform(0.004, 0.007)
            q_sfc = np.random.uniform(1.0e-5, 0.015)

        # Construct pressure interfaces
        p8w[:, i] = p_top + (p_sfc - p_top) * (eta ** 1.3)
        p8w[-1, i] = p_top  # Ensure exact top match

        # Layer pressures (log-mean midpoint)
        for k in range(nz):
            p3d[k, i] = 0.5 * (p8w[k, i] + p8w[k + 1, i])
            pi3d[k, i] = (p3d[k, i] / p1000mb) ** rcp

        # Construct temperature profile
        z_curr = 0.0
        t8w[0, i] = t_sfc
        for k in range(nz):
            # Estimate dz from hypsometric equation
            dp = p8w[k, i] - p8w[k + 1, i]
            # Standard lapse rate with stratospheric isothermal layer above 12 km
            if z_curr < 11000.0:
                t_lay = max(195.0, t_sfc - lapse_rate * z_curr)
            elif z_curr < 20000.0:
                t_lay = 216.65
            else:
                t_lay = 216.65 + 0.001 * (z_curr - 20000.0)

            t3d[k, i] = t_lay
            rho3d[k, i] = p3d[k, i] / (rd * t_lay)
            dz = dp / (rho3d[k, i] * g)
            dz8w[k, i] = dz
            z_curr += dz
            t8w[k + 1, i] = max(195.0, t_lay - 0.5 * lapse_rate * dz)

        # Water vapor profile (exponential decrease with height)
        scale_height = np.random.uniform(1800.0, 2500.0)
        z_mid = 0.0
        for k in range(nz):
            z_mid += 0.5 * dz8w[k, i]
            qv = q_sfc * np.exp(-z_mid / scale_height)
            qv3d[k, i] = max(1.0e-8, min(0.04, qv))
            z_mid += 0.5 * dz8w[k, i]

        # Clouds and hydrometeors
        if regime == 0:
            # Marine boundary layer cloud at k = 2..5
            cld_k = np.random.randint(2, 6)
            cldfra3d[cld_k, i] = np.random.uniform(0.5, 0.95)
            qc3d[cld_k, i] = np.random.uniform(0.2e-3, 0.8e-3)
        elif regime == 1:
            # Mid-level cloud (k=8..15) and cirrus (k=25..35)
            if np.random.rand() > 0.2:
                k_mid = np.random.randint(8, 16)
                cldfra3d[k_mid, i] = np.random.uniform(0.4, 0.9)
                qc3d[k_mid, i] = np.random.uniform(0.1e-3, 0.5e-3)
                qi3d[k_mid, i] = np.random.uniform(0.05e-3, 0.2e-3)
            if np.random.rand() > 0.3:
                k_cir = np.random.randint(24, 34)
                cldfra3d[k_cir, i] = np.random.uniform(0.2, 0.8)
                qi3d[k_cir, i] = np.random.uniform(0.01e-3, 0.1e-3)
        elif regime == 2:
            # Ice cloud / Diamond dust
            if np.random.rand() > 0.3:
                k_ice = np.random.randint(4, 18)
                cldfra3d[k_ice, i] = np.random.uniform(0.2, 0.7)
                qi3d[k_ice, i] = np.random.uniform(0.005e-3, 0.05e-3)
                qs3d[k_ice, i] = np.random.uniform(0.01e-3, 0.05e-3)
        elif regime == 3:
            # Deep convective column across k=3..32
            k_top = np.random.randint(25, 38)
            for k in range(3, k_top):
                cldfra3d[k, i] = np.random.uniform(0.7, 1.0)
                if t3d[k, i] > 273.15:
                    qc3d[k, i] = np.random.uniform(0.4e-3, 1.5e-3)
                    qr3d[k, i] = np.random.uniform(0.2e-3, 1.0e-3)
                elif t3d[k, i] > 253.15:
                    qc3d[k, i] = np.random.uniform(0.1e-3, 0.5e-3)
                    qi3d[k, i] = np.random.uniform(0.2e-3, 1.0e-3)
                    qs3d[k, i] = np.random.uniform(0.1e-3, 0.8e-3)
                    qg3d[k, i] = np.random.uniform(0.05e-3, 0.5e-3)
                else:
                    qi3d[k, i] = np.random.uniform(0.1e-3, 0.8e-3)
                    qs3d[k, i] = np.random.uniform(0.05e-3, 0.4e-3)
        else:
            # Edge cases
            case_type = i % 4
            if case_type == 0:
                # Completely clear sky: exactly zero hydrometeors and zero cloud fraction
                pass
            elif case_type == 1:
                # Completely overcast (cldfra = 1.0) in 3 contiguous layers
                k_start = np.random.randint(10, 20)
                for k in range(k_start, k_start + 3):
                    cldfra3d[k, i] = 1.0
                    qc3d[k, i] = 0.5e-3
            elif case_type == 2:
                # Fractional cloud boundary threshold (0.01)
                k_t = np.random.randint(5, 25)
                cldfra3d[k_t, i] = 0.01
                qc3d[k_t, i] = 1.0e-5
            elif case_type == 3:
                # Ultra-cold cirrus near stratopause
                k_hi = np.random.randint(30, 42)
                cldfra3d[k_hi, i] = 0.6
                qi3d[k_hi, i] = 2.0e-5

    # Write binary file
    with open(out_path, "wb") as f:
        # Header: magic (4B), num_cols (4B), nz (4B), p_top (8B), icloud (4B), cldovrlp (4B), idcor (4B), mp_physics (4B)
        magic = 0x5252544D  # 'RRTM' in ASCII
        f.write(struct.pack("<3Id4I", magic, num_cols, nz, p_top, icloud, cldovrlp, idcor, mp_physics))

        # 2D surface variables
        f.write(tsk.tobytes())
        f.write(emiss.tobytes())
        f.write(xland.tobytes())
        f.write(xice.tobytes())
        f.write(snow.tobytes())
        f.write(xlat.tobytes())

        # 3D interface variables: (nz+1, num_cols)
        f.write(p8w.tobytes())
        f.write(t8w.tobytes())

        # 3D layer variables: (nz, num_cols)
        f.write(p3d.tobytes())
        f.write(pi3d.tobytes())
        f.write(dz8w.tobytes())
        f.write(t3d.tobytes())
        f.write(rho3d.tobytes())
        f.write(qv3d.tobytes())
        f.write(qc3d.tobytes())
        f.write(qr3d.tobytes())
        f.write(qi3d.tobytes())
        f.write(qs3d.tobytes())
        f.write(qg3d.tobytes())
        f.write(cldfra3d.tobytes())

    print(f"Generated {num_cols} atmospheric columns with {nz} levels to {out_path}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate RRTMG LW test inputs")
    parser.add_argument("--cols", type=int, default=2500, help="Number of columns")
    parser.add_argument("--nz", type=int, default=50, help="Number of vertical levels")
    parser.add_argument("--ptop", type=float, default=5000.0, help="Model top pressure (Pa)")
    parser.add_argument("--out", type=str, default="input_rrtmg_lw.bin", help="Output path")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    generate_rrtmg_lw_inputs(args.cols, args.nz, args.ptop, args.out, args.seed)
