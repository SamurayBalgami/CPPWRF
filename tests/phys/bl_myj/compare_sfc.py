#!/usr/bin/env python3
import argparse
import sys
import numpy as np

VAR_NAMES = [
    "qsfc", "thz0", "qz0", "uz0", "vz0", "ustar", "znt", "pblh", "rmol", "akhs", "akms",
    "rib", "chs", "chs2", "cqs2", "hfx", "qfx", "flx_lh", "flhc", "flqc",
    "qgh", "cpm", "ct", "u10", "v10", "t02", "th02", "tshltr", "th10",
    "q02", "qshltr", "q10", "pshltr", "u10e", "v10e"
]

FLOORS = {
    "qsfc": 1e-6, "thz0": 1.0, "qz0": 1e-6, "uz0": 1e-3, "vz0": 1e-3,
    "ustar": 1e-3, "znt": 1e-5, "pblh": 1.0, "rmol": 1e-4, "akhs": 1e-4,
    "akms": 1e-4, "rib": 1e-3, "chs": 1e-4, "chs2": 1e-4, "cqs2": 1e-4,
    "hfx": 0.1, "qfx": 1e-6, "flx_lh": 1e-3, "flhc": 1e-3, "flqc": 1e-6,
    "qgh": 1e-6, "cpm": 1.0, "ct": 1e-3, "u10": 1e-3, "v10": 1e-3,
    "t02": 1.0, "th02": 1.0, "tshltr": 1.0, "th10": 1.0, "q02": 1e-6,
    "qshltr": 1e-6, "q10": 1e-6, "pshltr": 100.0, "u10e": 1e-3, "v10e": 1e-3
}

def compare_myjsfc_outputs(fortran_file, cpp_file, tol=0.0):
    with open(fortran_file, "rb") as ff, open(cpp_file, "rb") as cf:
        f_data = ff.read()
        c_data = cf.read()

    n_bytes = len(f_data)
    if n_bytes != len(c_data):
        print(f"ERROR: File size mismatch: Fortran={n_bytes} bytes, C++={len(c_data)} bytes")
        return False

    n_vars = len(VAR_NAMES)
    total_elements = n_bytes // 8
    n_cols = total_elements // n_vars

    if total_elements % n_vars != 0:
        print(f"ERROR: Total elements {total_elements} not divisible by {n_vars} variables")
        return False

    offset = 0
    all_pass = True

    print(f"Comparing Eta/MYJ surface layer outputs ({n_cols} columns, tol = {tol:e}):")
    print(f"{'Variable':<12} | {'Active Pts':<12} | {'Max Abs Diff':<14} | {'Max Rel Diff':<14} | {'Worst Col':<10} | {'Status':<6}")
    print("-" * 84)

    for name in VAR_NAMES:
        f_arr = np.frombuffer(f_data, dtype=np.float64, count=n_cols, offset=offset)
        c_arr = np.frombuffer(c_data, dtype=np.float64, count=n_cols, offset=offset)
        offset += n_cols * 8

        abs_diff = np.abs(f_arr - c_arr)
        floor = FLOORS.get(name, 1e-6)
        denom = np.maximum(np.abs(f_arr), floor)
        rel_diff = abs_diff / denom

        max_abs = np.max(abs_diff)
        max_rel = np.max(rel_diff)
        worst_col = int(np.argmax(abs_diff))

        active_pts = int(np.count_nonzero(np.abs(f_arr) > 1e-15))

        passed = (max_rel <= tol) or (max_abs <= tol)
        status = "PASS" if passed else "FAIL"
        if not passed:
            all_pass = False
            print(f"  [DIFF] {name} col={worst_col}: wrf={f_arr[worst_col]:.17e} cpp={c_arr[worst_col]:.17e} diff={abs_diff[worst_col]:.3e}")

        print(f"{name:<12} | {active_pts:<12} | {max_abs:<14.6e} | {max_rel:<14.6e} | {worst_col:<10} | {status:<6}")

    if all_pass:
        print("\nOVERALL VERDICT: PASS (all variables within tolerance)")
    else:
        print("\nOVERALL VERDICT: FAIL")

    return all_pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fortran_file", type=str)
    parser.add_argument("cpp_file", type=str)
    parser.add_argument("--tol", type=float, default=0.0)
    args = parser.parse_args()

    success = compare_myjsfc_outputs(args.fortran_file, args.cpp_file, args.tol)
    sys.exit(0 if success else 1)
