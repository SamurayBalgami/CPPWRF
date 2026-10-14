#!/usr/bin/env python3
"""
Differential comparison script for RRTMG-SW (ra_sw_physics=4).
Compares output_fortran.bin vs output_cpp.bin.
Verifies against 1e-12 relative tolerance threshold.
"""
import argparse
import sys
import numpy as np

def compare_sw(fortran_file="output_fortran.bin", cpp_file="output_cpp.bin", input_file="input_rrtmg_sw.bin", tol=1e-12):
    print("=" * 80)
    print("RRTMG SHORTWAVE (SW) DIFFERENTIAL VERIFICATION REPORT")
    print(f"Fortran Reference: {fortran_file}")
    print(f"C++ Target:        {cpp_file}")
    print(f"Tolerance Target:  {tol}")
    print("=" * 80)

    # Read inputs to determine active counts
    with open(input_file, "rb") as f:
        header = np.fromfile(f, dtype=np.int32, count=2)
        num_cols, nz = int(header[0]), int(header[1])
        coszen = np.fromfile(f, dtype=np.float64, count=num_cols)
        albedo = np.fromfile(f, dtype=np.float64, count=num_cols)
        xland = np.fromfile(f, dtype=np.float64, count=num_cols)
        xice = np.fromfile(f, dtype=np.float64, count=num_cols)
        snow = np.fromfile(f, dtype=np.float64, count=num_cols)
        xlat = np.fromfile(f, dtype=np.float64, count=num_cols)
        # Skip layer inputs up to cldfra3d
        p3d = np.fromfile(f, dtype=np.float64, count=nz * num_cols)
        pi3d = np.fromfile(f, dtype=np.float64, count=nz * num_cols)
        dz8w = np.fromfile(f, dtype=np.float64, count=nz * num_cols)
        t3d = np.fromfile(f, dtype=np.float64, count=nz * num_cols)
        rho3d = np.fromfile(f, dtype=np.float64, count=nz * num_cols)
        qv3d = np.fromfile(f, dtype=np.float64, count=nz * num_cols)
        qc3d = np.fromfile(f, dtype=np.float64, count=nz * num_cols)
        qr3d = np.fromfile(f, dtype=np.float64, count=nz * num_cols)
        qi3d = np.fromfile(f, dtype=np.float64, count=nz * num_cols)
        qs3d = np.fromfile(f, dtype=np.float64, count=nz * num_cols)
        qg3d = np.fromfile(f, dtype=np.float64, count=nz * num_cols)
        cldfra3d = np.fromfile(f, dtype=np.float64, count=nz * num_cols).reshape((nz, num_cols))

    day_cols = np.sum(coszen > 0.0)
    night_cols = np.sum(coszen <= 0.0)
    cloudy_cols = np.sum(np.any(cldfra3d > 0.0, axis=0))
    clear_cols = num_cols - cloudy_cols

    print(f"Total Columns:     {num_cols}")
    print(f"Vertical Levels:   {nz}")
    print(f"Daytime Columns:   {day_cols} ({100.0 * day_cols / num_cols:.1f}%)")
    print(f"Nighttime Columns: {night_cols} ({100.0 * night_cols / num_cols:.1f}%)")
    print(f"Cloudy Columns:    {cloudy_cols} ({100.0 * cloudy_cols / num_cols:.1f}%)")
    print(f"Clear-Sky Columns: {clear_cols} ({100.0 * clear_cols / num_cols:.1f}%)")
    print("-" * 80)

    # Read Fortran outputs
    with open(fortran_file, "rb") as f:
        f_gsw = np.fromfile(f, dtype=np.float64, count=num_cols)
        f_swcf = np.fromfile(f, dtype=np.float64, count=num_cols)
        f_swupt = np.fromfile(f, dtype=np.float64, count=num_cols)
        f_swuptc = np.fromfile(f, dtype=np.float64, count=num_cols)
        f_swdnt = np.fromfile(f, dtype=np.float64, count=num_cols)
        f_swdntc = np.fromfile(f, dtype=np.float64, count=num_cols)
        f_swupb = np.fromfile(f, dtype=np.float64, count=num_cols)
        f_swupbc = np.fromfile(f, dtype=np.float64, count=num_cols)
        f_swdnb = np.fromfile(f, dtype=np.float64, count=num_cols)
        f_swdnbc = np.fromfile(f, dtype=np.float64, count=num_cols)
        f_swvisdir = np.fromfile(f, dtype=np.float64, count=num_cols)
        f_swvisdif = np.fromfile(f, dtype=np.float64, count=num_cols)
        f_swnirdir = np.fromfile(f, dtype=np.float64, count=num_cols)
        f_swnirdif = np.fromfile(f, dtype=np.float64, count=num_cols)
        f_rthratensw = np.fromfile(f, dtype=np.float64, count=nz * num_cols).reshape((nz, num_cols))
        f_rthratenswc = np.fromfile(f, dtype=np.float64, count=nz * num_cols).reshape((nz, num_cols))

    # Read C++ outputs
    with open(cpp_file, "rb") as f:
        c_gsw = np.fromfile(f, dtype=np.float64, count=num_cols)
        c_swcf = np.fromfile(f, dtype=np.float64, count=num_cols)
        c_swupt = np.fromfile(f, dtype=np.float64, count=num_cols)
        c_swuptc = np.fromfile(f, dtype=np.float64, count=num_cols)
        c_swdnt = np.fromfile(f, dtype=np.float64, count=num_cols)
        c_swdntc = np.fromfile(f, dtype=np.float64, count=num_cols)
        c_swupb = np.fromfile(f, dtype=np.float64, count=num_cols)
        c_swupbc = np.fromfile(f, dtype=np.float64, count=num_cols)
        c_swdnb = np.fromfile(f, dtype=np.float64, count=num_cols)
        c_swdnbc = np.fromfile(f, dtype=np.float64, count=num_cols)
        c_swvisdir = np.fromfile(f, dtype=np.float64, count=num_cols)
        c_swvisdif = np.fromfile(f, dtype=np.float64, count=num_cols)
        c_swnirdir = np.fromfile(f, dtype=np.float64, count=num_cols)
        c_swnirdif = np.fromfile(f, dtype=np.float64, count=num_cols)
        c_rthratensw = np.fromfile(f, dtype=np.float64, count=nz * num_cols).reshape((nz, num_cols))
        c_rthratenswc = np.fromfile(f, dtype=np.float64, count=nz * num_cols).reshape((nz, num_cols))

    variables = [
        ("gsw", f_gsw, c_gsw),
        ("swcf", f_swcf, c_swcf),
        ("swupt", f_swupt, c_swupt),
        ("swuptc", f_swuptc, c_swuptc),
        ("swdnt", f_swdnt, c_swdnt),
        ("swdntc", f_swdntc, c_swdntc),
        ("swupb", f_swupb, c_swupb),
        ("swupbc", f_swupbc, c_swupbc),
        ("swdnb", f_swdnb, c_swdnb),
        ("swdnbc", f_swdnbc, c_swdnbc),
        ("swvisdir", f_swvisdir, c_swvisdir),
        ("swvisdif", f_swvisdif, c_swvisdif),
        ("swnirdir", f_swnirdir, c_swnirdir),
        ("swnirdif", f_swnirdif, c_swnirdif),
        ("rthratensw", f_rthratensw, c_rthratensw),
        ("rthratenswc", f_rthratenswc, c_rthratenswc),
    ]

    all_pass = True
    print(f"{'Variable':<14} | {'Max Abs Err':<14} | {'Max Rel Err':<14} | {'Mean Abs Err':<14} | {'Status'}")
    print("-" * 80)

    for name, f_arr, c_arr in variables:
        abs_err = np.abs(f_arr - c_arr)
        max_abs = np.max(abs_err)
        mean_abs = np.mean(abs_err)

        denom = np.maximum(np.abs(f_arr), 1e-12)
        rel_err = abs_err / denom
        max_rel = np.max(rel_err)

        passed = (max_rel <= tol) or (max_abs <= 1e-12)
        if not passed:
            all_pass = False
        status_str = "PASS (bit/num)" if passed else "FAIL"
        print(f"{name:<14} | {max_abs:<14.6e} | {max_rel:<14.6e} | {mean_abs:<14.6e} | {status_str}")

    print("=" * 80)
    if all_pass:
        print(f"ALL VARIABLES PASSED DIFFERENTIAL VERIFICATION (tol <= {tol})!")
    else:
        print("WARNING: Some variables exceeded tolerance threshold.")
    print("=" * 80)

    return all_pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Differential test for RRTMG SW")
    parser.add_argument("--fortran", default="output_fortran.bin")
    parser.add_argument("--cpp", default="output_cpp.bin")
    parser.add_argument("--input", default="input_rrtmg_sw.bin")
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    success = compare_sw(args.fortran, args.cpp, args.input, args.tol)
    sys.exit(0 if success else 1)
