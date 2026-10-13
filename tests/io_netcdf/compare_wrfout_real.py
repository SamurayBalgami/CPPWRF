#!/usr/bin/env python3
"""Compares CPPWRF wrfout against WRF reference wrfout for Real R0 case.
Checks variables: P, PB, PH, PHB, MU, MUB, T, ZNU, C1H, HGT, XLAT.
"""
import os
import sys
import numpy as np

try:
    import netCDF4 as nc
except ImportError:
    print("Warning: netCDF4 Python module not found; cannot run compare_wrfout_real.py")
    sys.exit(0)

def main():
    cpp_path = sys.argv[1] if len(sys.argv) > 1 else "tests/io_netcdf/test_out/real_r0_wrfout/wrfout_d01_2026-10-12_00:00:00"
    wrf_path = sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser("~/wrf_real/run_R0_r8/wrfout_d01_2026-10-12_00:00:00")

    if not os.path.exists(cpp_path):
        print(f"Error: CPPWRF wrfout file not found: {cpp_path}")
        return 1
    if not os.path.exists(wrf_path):
        print(f"Notice: WRF reference wrfout file not found: {wrf_path} (skipping)")
        return 0

    print("=========================================================================================")
    print("CPPWRF vs WRF Real R0 wrfout Field Comparison")
    print(f"  CPPWRF: {cpp_path}")
    print(f"  WRF ref:{wrf_path}")
    print("=========================================================================================")

    ds_cpp = nc.Dataset(cpp_path, "r")
    ds_wrf = nc.Dataset(wrf_path, "r")

    vars_to_check = [
        ("P", 1e-2, 1e-4),       # Perturbation pressure (Pa)
        ("PB", 1e-1, 1e-5),      # Base state pressure (Pa)
        ("PH", 1e-1, 1e-4),      # Perturbation geopotential (m2 s-2)
        ("PHB", 1.0, 1e-5),      # Base state geopotential (m2 s-2)
        ("MU", 1e-2, 1e-4),      # Perturbation dry air mass (Pa)
        ("MUB", 1e-1, 1e-5),     # Base state dry air mass (Pa)
        ("T", 1e-3, 1e-4),       # Perturbation potential temperature theta-300 (K)
        ("ZNU", 1e-5, 1e-5),     # Eta on half levels
        ("C1H", 1e-5, 1e-5),     # Hybrid coordinate C1H
        ("HGT", 1e-3, 1e-5),     # Terrain height (m)
        ("XLAT", 1e-4, 1e-5),    # Latitude (degrees)
    ]

    header = f"{'Variable':<10} | {'Max Abs Diff':<14} | {'Max Rel Diff':<14} | {'CPPWRF Range':<24} | {'WRF Range':<24} | {'Status'}"
    print(header)
    print("-" * len(header))

    all_passed = True

    for vname, abs_tol, rel_tol in vars_to_check:
        if vname not in ds_cpp.variables:
            print(f"{vname:<10} | MISSING IN CPPWRF")
            all_passed = False
            continue
        if vname not in ds_wrf.variables:
            print(f"{vname:<10} | MISSING IN WRF")
            all_passed = False
            continue

        a = np.array(ds_cpp.variables[vname][0], dtype=np.float64)
        b = np.array(ds_wrf.variables[vname][0], dtype=np.float64)

        if a.shape != b.shape:
            print(f"{vname:<10} | SHAPE MISMATCH: {a.shape} vs {b.shape}")
            all_passed = False
            continue

        diff = np.abs(a - b)
        max_abs = float(np.max(diff))
        scale = np.maximum(np.abs(b), 1.0)
        max_rel = float(np.max(diff / scale))

        a_min, a_max = float(np.min(a)), float(np.max(a))
        b_min, b_max = float(np.min(b)), float(np.max(b))

        passed = (max_abs <= abs_tol) or (max_rel <= rel_tol)
        status_str = "PASS" if passed else "FAIL"
        if not passed:
            all_passed = False

        print(f"{vname:<10} | {max_abs:<14.5e} | {max_rel:<14.5e} | [{a_min:9.3f}, {a_max:9.3f}] | [{b_min:9.3f}, {b_max:9.3f}] | {status_str}")

    print("=========================================================================================")
    if all_passed:
        print("ALL REAL R0 WRFOUT FIELDS MATCH WRF REFERENCE WITHIN FLOAT ROUNDOFF TOLERANCES!")
        return 0
    else:
        print("SOME FIELDS EXCEEDED TOLERANCE!")
        return 1

if __name__ == "__main__":
    sys.exit(main())
