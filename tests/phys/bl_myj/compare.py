#!/usr/bin/env python3
"""
Differential validator for MYJ PBL scheme (bl_myj).
Compares unmodified Fortran WRF 4.8.0 outputs vs C++ outputs.
Checks 7 2D variables and 9 3D variables across all columns and levels.
"""

import argparse
import sys
import numpy as np

VARS_2D = [
    "qsfc", "thz0", "qz0", "ct", "pblh", "kpbl", "mixht"
]

VARS_3D = [
    "tke_myj", "exch_h", "el_myj", "rublten", "rvblten", "rthblten",
    "rqvblten", "rqcblten", "rqiblten"
]

FLOORS = {
    "qsfc": 1e-6,
    "thz0": 1.0,
    "qz0": 1e-6,
    "ct": 1e-3,
    "pblh": 1.0,
    "kpbl": 1.0,
    "mixht": 1.0,
    "tke_myj": 1e-3,
    "exch_h": 1e-4,
    "el_myj": 1e-3,
    "rublten": 1e-5,
    "rvblten": 1e-5,
    "rthblten": 1e-5,
    "rqvblten": 1e-8,
    "rqcblten": 1e-8,
    "rqiblten": 1e-8,
}

def compare_myjpbl_outputs(fortran_file, cpp_file, tol=0.0, nz=44):
    with open(fortran_file, "rb") as ff, open(cpp_file, "rb") as cf:
        f_data = ff.read()
        c_data = cf.read()

    n_bytes = len(f_data)
    if n_bytes != len(c_data):
        print(f"ERROR: File size mismatch: Fortran={n_bytes} bytes, C++={len(c_data)} bytes")
        return False

    total_elements = n_bytes // 8
    elements_per_col = len(VARS_2D) + len(VARS_3D) * nz
    if total_elements % elements_per_col != 0:
        print(f"ERROR: Total elements {total_elements} not divisible by {elements_per_col}")
        return False

    n_cols = total_elements // elements_per_col
    offset = 0
    all_pass = True

    print(f"Comparing MYJ PBL scheme outputs ({n_cols} columns, nz={nz}, tol = {tol:e}):")
    print(f"{'Variable':<12} | {'Shape':<10} | {'Active Pts':<12} | {'Max Abs Diff':<14} | {'Max Rel Diff':<14} | {'Status':<6}")
    print("-" * 80)

    for name in VARS_2D:
        f_arr = np.frombuffer(f_data, dtype=np.float64, count=n_cols, offset=offset)
        c_arr = np.frombuffer(c_data, dtype=np.float64, count=n_cols, offset=offset)
        offset += n_cols * 8

        abs_diff = np.abs(f_arr - c_arr)
        floor = FLOORS.get(name, 1e-6)
        denom = np.maximum(np.abs(f_arr), floor)
        rel_diff = abs_diff / denom

        max_abs = np.max(abs_diff)
        max_rel = np.max(rel_diff)
        active_pts = int(np.count_nonzero(np.abs(f_arr) > 1e-15))

        passed = (max_rel <= tol) or (max_abs <= tol)
        status = "PASS" if passed else "FAIL"
        if not passed:
            all_pass = False
            worst = int(np.argmax(abs_diff))
            print(f"  [DIFF] {name} col={worst}: wrf={f_arr[worst]:.17e} cpp={c_arr[worst]:.17e} diff={abs_diff[worst]:.3e}")

        shape_str = f"({n_cols})"
        print(f"{name:<12} | {shape_str:<10} | {active_pts:<12} | {max_abs:<14.6e} | {max_rel:<14.6e} | {status:<6}")

    n3d = n_cols * nz
    for name in VARS_3D:
        f_arr = np.frombuffer(f_data, dtype=np.float64, count=n3d, offset=offset)
        c_arr = np.frombuffer(c_data, dtype=np.float64, count=n3d, offset=offset)
        offset += n3d * 8

        abs_diff = np.abs(f_arr - c_arr)
        floor = FLOORS.get(name, 1e-6)
        denom = np.maximum(np.abs(f_arr), floor)
        rel_diff = abs_diff / denom

        max_abs = np.max(abs_diff)
        max_rel = np.max(rel_diff)
        active_pts = int(np.count_nonzero(np.abs(f_arr) > 1e-15))

        passed = (max_rel <= tol) or (max_abs <= tol)
        status = "PASS" if passed else "FAIL"
        if not passed:
            all_pass = False
            worst = int(np.argmax(abs_diff))
            k_w, c_w = divmod(worst, n_cols)
            print(f"  [DIFF] {name} col={c_w} k={k_w}: wrf={f_arr[worst]:.17e} cpp={c_arr[worst]:.17e} diff={abs_diff[worst]:.3e}")

        shape_str = f"({n_cols},{nz})"
        print(f"{name:<12} | {shape_str:<10} | {active_pts:<12} | {max_abs:<14.6e} | {max_rel:<14.6e} | {status:<6}")

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
    parser.add_argument("--nz", type=int, default=44)
    args = parser.parse_args()

    success = compare_myjpbl_outputs(args.fortran_file, args.cpp_file, args.tol, args.nz)
    sys.exit(0 if success else 1)
