#!/usr/bin/env python3
"""
Input generator for RRTMG Shortwave (ra_sw_physics=4) differential verification.
Generates realistic, physically diverse atmospheric columns covering tropical,
midlatitude, arctic, deep convective, dawn/dusk, daytime, and nighttime edge cases.
"""
import argparse
import struct
import numpy as np

def generate_rrtmg_sw_inputs(num_cols=2500, nz=50, p_top=5000.0, out_path="input_rrtmg_sw.bin", seed=42):
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

    # 2D solar and surface variables: (num_cols,)
    coszen = np.zeros(num_cols, dtype=np.float64)
    albedo = np.zeros(num_cols, dtype=np.float64)
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

    eta = np.linspace(1.0, 0.0, nz + 1)

    for i in range(num_cols):
        regime = i % 5  # 0: Tropical, 1: Mid-lat, 2: Dry/Arctic, 3: Deep convection, 4: Edge cases / Night

        # Solar zenith angle variation
        # 10% strictly nighttime (coszen <= 0), 10% dawn/dusk (0 < coszen < 0.15), 80% daytime (0.15 .. 1.0)
        u_zen = np.random.rand()
        if u_zen < 0.10:
            coszen[i] = np.random.uniform(-0.8, 0.0) # Night
        elif u_zen < 0.20:
            coszen[i] = np.random.uniform(0.01, 0.15) # Dawn/dusk
        else:
            coszen[i] = np.random.uniform(0.15, 1.0) # Daytime

        if regime == 0:
            # Tropical: warm, high moisture, marine boundary layer clouds
            p_sfc = np.random.uniform(100500.0, 101800.0)
            xland[i] = 2.0 if np.random.rand() > 0.3 else 1.0
            xice[i] = 0.0
            snow[i] = 0.0
            xlat[i] = np.random.uniform(-20.0, 20.0)
            albedo[i] = np.random.uniform(0.05, 0.08) if xland[i] == 2.0 else np.random.uniform(0.12, 0.20)
            t_sfc = np.random.uniform(298.0, 304.0)
            lapse_rate = np.random.uniform(0.0055, 0.0065)
            q_sfc = np.random.uniform(0.015, 0.022)

        elif regime == 1:
            # Mid-latitude: moderate temps, mixed phase, mid/high clouds
            p_sfc = np.random.uniform(98000.0, 102500.0)
            xland[i] = 1.0 if np.random.rand() > 0.4 else 2.0
            xice[i] = 0.0
            snow[i] = np.random.uniform(0.0, 50.0) if np.random.rand() > 0.7 else 0.0
            xlat[i] = np.random.uniform(30.0, 60.0) * (1.0 if np.random.rand() > 0.5 else -1.0)
            albedo[i] = np.random.uniform(0.15, 0.25) if snow[i] == 0.0 else np.random.uniform(0.4, 0.7)
            t_sfc = np.random.uniform(275.0, 295.0)
            lapse_rate = np.random.uniform(0.006, 0.0075)
            q_sfc = np.random.uniform(0.005, 0.012)

        elif regime == 2:
            # Dry / Arctic: cold surface, ice clouds, snow/ice high albedo
            p_sfc = np.random.uniform(97000.0, 103000.0)
            xland[i] = 1.0 if np.random.rand() > 0.5 else 2.0
            xice[i] = 1.0 if (xland[i] == 2.0) else 0.0
            snow[i] = np.random.uniform(50.0, 500.0)
            xlat[i] = np.random.uniform(65.0, 85.0) * (1.0 if np.random.rand() > 0.5 else -1.0)
            albedo[i] = np.random.uniform(0.65, 0.85)
            t_sfc = np.random.uniform(230.0, 265.0)
            lapse_rate = np.random.uniform(0.003, 0.005)
            q_sfc = np.random.uniform(0.0001, 0.002)

        elif regime == 3:
            # Deep convection: very thick high clouds, mixed liquid/ice/snow/graupel
            p_sfc = np.random.uniform(99500.0, 101500.0)
            xland[i] = 2.0 if np.random.rand() > 0.5 else 1.0
            xice[i] = 0.0
            snow[i] = 0.0
            xlat[i] = np.random.uniform(-15.0, 35.0)
            albedo[i] = np.random.uniform(0.06, 0.18)
            t_sfc = np.random.uniform(295.0, 305.0)
            lapse_rate = np.random.uniform(0.006, 0.007)
            q_sfc = np.random.uniform(0.014, 0.020)

        else:
            # Regime 4: Edge cases (clear sky, high zenith angle, super dry)
            p_sfc = np.random.uniform(95000.0, 103500.0)
            xland[i] = 1.0 if np.random.rand() > 0.5 else 2.0
            xice[i] = 1.0 if (xland[i] == 2.0 and np.random.rand() > 0.5) else 0.0
            snow[i] = np.random.uniform(0.0, 200.0) if np.random.rand() > 0.5 else 0.0
            xlat[i] = np.random.uniform(-80.0, 80.0)
            albedo[i] = np.random.uniform(0.05, 0.90)
            t_sfc = np.random.uniform(240.0, 310.0)
            lapse_rate = np.random.uniform(0.004, 0.008)
            q_sfc = np.random.uniform(0.00005, 0.018)

        # Build pressure profile p8w (interface pressures)
        p8w[:, i] = p_top + (p_sfc - p_top) * (eta ** 1.3)

        # Build temperature profile t8w
        z_approx = (rd * t_sfc / g) * np.log(p_sfc / p8w[:, i])
        t_prof = t_sfc - lapse_rate * z_approx
        t_prof = np.maximum(t_prof, 185.0)
        t8w[:, i] = t_prof

        # Layer midpoint pressures p3d and temperatures t3d
        for k in range(nz):
            p3d[k, i] = 0.5 * (p8w[k, i] + p8w[k + 1, i])
            pi3d[k, i] = (p3d[k, i] / p1000mb) ** rcp
            t3d[k, i] = 0.5 * (t8w[k, i] + t8w[k + 1, i])

            dp = p8w[k, i] - p8w[k + 1, i]
            tv = t3d[k, i] * (1.0 + 0.608 * q_sfc)
            rho3d[k, i] = p3d[k, i] / (rd * tv)
            dz8w[k, i] = dp / (rho3d[k, i] * g)

            # Water vapor profile
            p_ratio = p3d[k, i] / p_sfc
            qv3d[k, i] = q_sfc * (p_ratio ** 2.5)
            qv3d[k, i] = max(qv3d[k, i], 1.0e-7)

            # Hydrometeors and cloud fraction
            if regime == 0:
                # Boundary layer clouds below 850 hPa
                if p3d[k, i] > 85000.0 and np.random.rand() > 0.4:
                    cldfra3d[k, i] = np.random.uniform(0.2, 0.9)
                    qc3d[k, i] = np.random.uniform(1.0e-4, 8.0e-4)
            elif regime == 1:
                # Mid-level and high-level clouds
                if 40000.0 < p3d[k, i] < 75000.0 and np.random.rand() > 0.5:
                    cldfra3d[k, i] = np.random.uniform(0.3, 0.95)
                    if t3d[k, i] > 268.0:
                        qc3d[k, i] = np.random.uniform(5.0e-5, 5.0e-4)
                    else:
                        qi3d[k, i] = np.random.uniform(1.0e-5, 1.5e-4)
                        qs3d[k, i] = np.random.uniform(5.0e-6, 1.0e-4)
            elif regime == 2:
                # Thin Arctic ice clouds
                if 30000.0 < p3d[k, i] < 70000.0 and np.random.rand() > 0.6:
                    cldfra3d[k, i] = np.random.uniform(0.1, 0.8)
                    qi3d[k, i] = np.random.uniform(5.0e-6, 8.0e-5)
            elif regime == 3:
                # Deep convective column
                if p3d[k, i] > 20000.0:
                    cldfra3d[k, i] = np.random.uniform(0.6, 1.0)
                    if t3d[k, i] > 273.15:
                        qc3d[k, i] = np.random.uniform(2.0e-4, 1.5e-3)
                        qr3d[k, i] = np.random.uniform(1.0e-4, 8.0e-4)
                    else:
                        qi3d[k, i] = np.random.uniform(5.0e-5, 5.0e-4)
                        qs3d[k, i] = np.random.uniform(1.0e-4, 1.0e-3)
                        qg3d[k, i] = np.random.uniform(5.0e-5, 5.0e-4)
            elif regime == 4:
                # Clear sky or random
                if i % 2 == 0:
                    cldfra3d[k, i] = 0.0
                else:
                    cldfra3d[k, i] = np.random.uniform(0.0, 0.4) if np.random.rand() > 0.7 else 0.0

    # Write out binary file:
    # Header: num_cols (int32), nz (int32)
    # 2D variables: coszen, albedo, xland, xice, snow, xlat
    # 3D layer variables: p3d, pi3d, dz8w, t3d, rho3d, qv3d, qc3d, qr3d, qi3d, qs3d, qg3d, cldfra3d
    # 3D interface variables: p8w, t8w
    with open(out_path, "wb") as f:
        f.write(struct.pack("<ii", num_cols, nz))
        f.write(coszen.tobytes())
        f.write(albedo.tobytes())
        f.write(xland.tobytes())
        f.write(xice.tobytes())
        f.write(snow.tobytes())
        f.write(xlat.tobytes())
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
        f.write(p8w.tobytes())
        f.write(t8w.tobytes())

    print(f"Generated {num_cols} columns x {nz} levels of RRTMG-SW test inputs to {out_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate RRTMG-SW test inputs")
    parser.add_argument("--cols", type=int, default=2500, help="Number of columns")
    parser.add_argument("--nz", type=int, default=50, help="Number of vertical levels")
    parser.add_argument("--out", type=str, default="input_rrtmg_sw.bin", help="Output file path")
    args = parser.parse_args()
    generate_rrtmg_sw_inputs(num_cols=args.cols, nz=args.nz, out_path=args.out)
