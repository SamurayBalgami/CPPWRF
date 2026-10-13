#!/usr/bin/env python3
import argparse
import sys
import numpy as np

VAR_NAMES_3D = [
    "rublten", "rvblten", "rthblten", "rqvblten", "rqcblten", "rqiblten",
    "exch_h", "exch_m"
]

VAR_NAMES_2D = [
    "hpbl", "wstar", "delta", "kpbl"
]

FLOORS = {
    "rublten": 1e-6,
    "rvblten": 1e-6,
    "rthblten": 1e-6,
    "rqvblten": 1e-9,
    "rqcblten": 1e-9,
    "rqiblten": 1e-9,
    "exch_h": 1e-3,
    "exch_m": 1e-3,
    "hpbl": 10.0,
    "wstar": 0.01,
    "delta": 1.0,
    "kpbl": 1.0
}

def compare_ysu_outputs(fortran_file, cpp_file, n_cols=2500, nz=50, tol=1e-12):
    with open(fortran_file, "rb") as ff, open(cpp_file, "rb") as cf:
        f_data = ff.read()
        c_data = cf.read()

    n_bytes = len(f_data)
    if n_bytes != len(c_data):
        print(f"ERROR: File size mismatch: Fortran={n_bytes} bytes, C++={len(c_data)} bytes")
        return False

    n3d = n_cols * nz
    n2d = n_cols
    expected_bytes = (len(VAR_NAMES_3D) * n3d + len(VAR_NAMES_2D) * n2d) * 8
    if n_bytes != expected_bytes:
        print(f"ERROR: Expected {expected_bytes} bytes, got {n_bytes}")
        return False

    offset = 0
    all_pass = True

    print(f"Comparing YSU outputs ({n_cols} cols x {nz} levels, tol = {tol:e}):")
    print(f"{'Variable':<12} | {'Active Pts':<12} | {'Max Abs Diff':<14} | {'Max Rel Diff':<14} | {'Worst Col':<10} | {'Status':<6}")
    print("-" * 84)

    # 3D arrays
    for name in VAR_NAMES_3D:
        f_arr = np.frombuffer(f_data, dtype=np.float64, count=n3d, offset=offset)
        c_arr = np.frombuffer(c_data, dtype=np.float64, count=n3d, offset=offset)
        offset += n3d * 8

        abs_diff = np.abs(f_arr - c_arr)
        floor = FLOORS.get(name, 1e-6)
        denom = np.maximum(np.abs(f_arr), floor)
        rel_diff = abs_diff / denom

        max_abs = np.max(abs_diff)
        max_rel = np.max(rel_diff)
        worst_idx = np.argmax(rel_diff)
        worst_col = worst_idx % n_cols

        active_pts = np.count_nonzero(np.abs(f_arr) > 1e-15)

        passed = (max_rel <= tol) or (max_abs <= tol)
        status = "PASS" if passed else "FAIL"
        if not passed:
            all_pass = False

        print(f"{name:<12} | {active_pts:<12} | {max_abs:<14.6e} | {max_rel:<14.6e} | {worst_col:<10} | {status:<6}")

    # 2D arrays
    for name in VAR_NAMES_2D:
        f_arr = np.frombuffer(f_data, dtype=np.float64, count=n2d, offset=offset)
        c_arr = np.frombuffer(c_data, dtype=np.float64, count=n2d, offset=offset)
        offset += n2d * 8

        abs_diff = np.abs(f_arr - c_arr)
        floor = FLOORS.get(name, 1e-6)
        denom = np.maximum(np.abs(f_arr), floor)
        rel_diff = abs_diff / denom

        max_abs = np.max(abs_diff)
        max_rel = np.max(rel_diff)
        worst_col = np.argmax(rel_diff)

        active_pts = np.count_nonzero(np.abs(f_arr) > 1e-15)

        passed = (max_rel <= tol) or (max_abs <= tol)
        status = "PASS" if passed else "FAIL"
        if not passed:
            all_pass = False

        print(f"{name:<12} | {active_pts:<12} | {max_abs:<14.6e} | {max_rel:<14.6e} | {worst_col:<10} | {status:<6}")

    print("-" * 84)
    if all_pass:
        print("ALL YSU VARIABLES PASSED")
    else:
        print("SOME YSU VARIABLES FAILED")
    return all_pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fortran_out", type=str)
    parser.add_argument("cpp_out", type=str)
    parser.add_argument("--cols", type=int, default=2500)
    parser.add_argument("--nz", type=int, default=50)
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    ok = compare_ysu_outputs(args.fortran_out, args.cpp_out, args.cols, args.nz, args.tol)
    sys.exit(0 if ok else 1)
