#!/usr/bin/env python3
"""
Compares Fortran and C++ outputs for the WRF physics driver differential test.
Acceptance criterion: PASS if max relative difference <= tol (default 1e-12 for FP64)
or bit-exact across all fields.
"""

import argparse
import sys
import numpy as np

def compare_outputs(fort_file, cpp_file, nx=16, ny=16, nz=35, n_moist=6, tol=1e-12):
    # Definition of output fields and their shapes
    fields = [
        # phy_prep outputs
        ("th_phy", (nz, ny, nx)),
        ("th_phy_m_t0", (nz, ny, nx)),
        ("p_phy", (nz, ny, nx)),
        ("pi_phy", (nz, ny, nx)),
        ("t_phy", (nz, ny, nx)),
        ("rho", (nz, ny, nx)),
        ("u_phy", (nz, ny, nx)),
        ("v_phy", (nz, ny, nx)),
        ("z", (nz, ny, nx)),
        ("dz8w", (nz + 1, ny, nx)),
        ("z_at_w", (nz + 1, ny, nx)),
        ("p8w", (nz + 1, ny, nx)),
        ("t8w", (nz + 1, ny, nx)),
        ("p_hyd", (nz, ny, nx)),
        ("p_hyd_w", (nz + 1, ny, nx)),
        # Tendency coupling outputs
        ("ru_tendf", (nz, ny, nx + 1)),
        ("rv_tendf", (nz, ny + 1, nx)),
        ("t_tendf", (nz, ny, nx)),
        ("qv_tendf", (nz, ny, nx)),
        ("qc_tendf", (nz, ny, nx)),
        ("qr_tendf", (nz, ny, nx)),
        ("qi_tendf", (nz, ny, nx)),
        ("qs_tendf", (nz, ny, nx)),
        ("qg_tendf", (nz, ny, nx)),
        ("ru_tend", (nz, ny, nx + 1)),
        ("rv_tend", (nz, ny + 1, nx)),
        ("t_tend", (nz, ny, nx)),
    ]

    all_passed = True
    results = []

    FLOORS = {
        "th_phy": 1.0,
        "th_phy_m_t0": 1.0,
        "p_phy": 1.0,
        "pi_phy": 0.01,
        "t_phy": 1.0,
        "rho": 1e-4,
        "u_phy": 0.1,
        "v_phy": 0.1,
        "z": 1.0,
        "dz8w": 1.0,
        "z_at_w": 1.0,
        "p8w": 1.0,
        "t8w": 1.0,
        "p_hyd": 1.0,
        "p_hyd_w": 1.0,
        "ru_tendf": 1e-6,
        "rv_tendf": 1e-6,
        "t_tendf": 1e-6,
        "qv_tendf": 1e-12,
        "qc_tendf": 1e-12,
        "qr_tendf": 1e-12,
        "qi_tendf": 1e-12,
        "qs_tendf": 1e-12,
        "qg_tendf": 1e-12,
        "ru_tend": 1e-6,
        "rv_tend": 1e-6,
        "t_tend": 1e-6,
    }

    with open(fort_file, 'rb') as f_fort, open(cpp_file, 'rb') as f_cpp:
        for name, shape in fields:
            count = int(np.prod(shape))
            fort_data = np.fromfile(f_fort, dtype=np.float64, count=count)
            cpp_data = np.fromfile(f_cpp, dtype=np.float64, count=count)

            if len(fort_data) != count or len(cpp_data) != count:
                print(f"Error: Incomplete data for {name}. Expected {count}, got Fortran {len(fort_data)}, CPP {len(cpp_data)}")
                return False

            fort_arr = fort_data.reshape(shape)
            cpp_arr = cpp_data.reshape(shape)

            diff = np.abs(fort_arr - cpp_arr)
            max_abs = np.max(diff)

            # Denominator with floor to avoid div by zero / near-zero at crossings
            floor = FLOORS.get(name, 1e-6)
            denom = np.maximum(np.abs(fort_arr), floor)
            rel_diff = diff / denom
            max_rel = np.max(rel_diff)

            worst_idx = np.unravel_index(np.argmax(rel_diff), shape)
            worst_str = f"{worst_idx}"

            status = "PASS" if max_rel <= tol else "FAIL"
            if status == "FAIL":
                all_passed = False

            results.append((name, max_abs, max_rel, worst_str, status))

    # Print summary table
    print("\n" + "=" * 80)
    print(f"{'Field':<16} | {'Max Abs Diff':<14} | {'Max Rel Diff':<14} | {'Worst (k,j,i)':<14} | {'Status':<6}")
    print("-" * 80)
    for name, max_abs, max_rel, worst_idx, status in results:
        print(f"{name:<16} | {max_abs:<14.6e} | {max_rel:<14.6e} | {worst_idx:<14} | {status:<6}")
    print("=" * 80)

    if all_passed:
        print(f"\n>> ALL {len(fields)} FIELDS PASSED DIFFERENTIAL VERIFICATION (tolerance = {tol}) <<\n")
    else:
        print(f"\n>> VERIFICATION FAILED (some fields exceeded tolerance {tol}) <<\n")

    return all_passed

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('fortran_out', help='Fortran output binary')
    parser.add_argument('cpp_out', help='C++ output binary')
    parser.add_argument('--nx', type=int, default=16)
    parser.add_argument('--ny', type=int, default=16)
    parser.add_argument('--nz', type=int, default=35)
    parser.add_argument('--tol', type=float, default=1e-12)
    args = parser.parse_args()

    passed = compare_outputs(args.fortran_out, args.cpp_out, nx=args.nx, ny=args.ny, nz=args.nz, tol=args.tol)
    sys.exit(0 if passed else 1)
