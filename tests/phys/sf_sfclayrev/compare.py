#!/usr/bin/env python3
import argparse
import sys
import numpy as np

VAR_NAMES = [
    "regime", "hfx", "qfx", "lh", "qsfc", "mol", "rmol", "gz1oz0",
    "wspd", "br", "psim", "psih", "fm", "fh", "znt", "zol", "ust", "cpm",
    "chs2", "cqs2", "chs", "flhc", "flqc", "qgh", "ustm",
    "u10", "v10", "th2", "t2", "q2", "ck", "cka", "cd", "cda"
]

FLOORS = {
    "regime": 1.0,
    "hfx": 1e-3,
    "qfx": 1e-8,
    "lh": 1e-3,
    "qsfc": 1e-7,
    "mol": 1e-5,
    "rmol": 1e-5,
    "gz1oz0": 1e-3,
    "wspd": 0.1,
    "br": 1e-4,
    "psim": 1e-4,
    "psih": 1e-4,
    "fm": 1e-4,
    "fh": 1e-4,
    "znt": 1e-7,
    "zol": 1e-5,
    "ust": 1e-4,
    "cpm": 100.0,
    "chs2": 1e-6,
    "cqs2": 1e-6,
    "chs": 1e-6,
    "flhc": 1e-4,
    "flqc": 1e-6,
    "qgh": 1e-7,
    "ustm": 1e-4,
    "u10": 0.1,
    "v10": 0.1,
    "th2": 10.0,
    "t2": 10.0,
    "q2": 1e-7,
    "ck": 1e-6,
    "cka": 1e-6,
    "cd": 1e-6,
    "cda": 1e-6
}

def compare_outputs(fortran_file, cpp_file, tol=1e-12):
    with open(fortran_file, "rb") as ff, open(cpp_file, "rb") as cf:
        f_data = ff.read()
        c_data = cf.read()

    n_bytes = len(f_data)
    if n_bytes != len(c_data):
        print(f"ERROR: File size mismatch: Fortran={n_bytes} bytes, C++={len(c_data)} bytes")
        return False

    n_vars = len(VAR_NAMES)
    bytes_per_var = n_bytes // n_vars
    n_cols = bytes_per_var // 8

    print(f"Comparing {n_cols} columns across {n_vars} variables (tol = {tol:e}):")
    print(f"{'Variable':<12} | {'Active Pts':<10} | {'Max Abs Diff':<14} | {'Max Rel Diff':<14} | {'Worst Col':<10} | {'Status':<6}")
    print("-" * 80)

    all_pass = True
    f_arrs = np.frombuffer(f_data, dtype=np.float64).reshape((n_vars, n_cols))
    c_arrs = np.frombuffer(c_data, dtype=np.float64).reshape((n_vars, n_cols))

    for idx, name in enumerate(VAR_NAMES):
        f_val = f_arrs[idx]
        c_val = c_arrs[idx]

        abs_diff = np.abs(f_val - c_val)
        floor = FLOORS.get(name, 1e-6)
        denom = np.maximum(np.abs(f_val), floor)
        rel_diff = abs_diff / denom

        max_abs = np.max(abs_diff)
        max_rel = np.max(rel_diff)
        worst_col = np.argmax(rel_diff)

        # Active points: points where variable is non-zero (or physically active)
        active_pts = np.count_nonzero(np.abs(f_val) > 1e-15)

        passed = max_rel <= tol
        status = "PASS" if passed else "FAIL"
        if not passed:
            all_pass = False

        print(f"{name:<12} | {active_pts:<10} | {max_abs:<14.6e} | {max_rel:<14.6e} | {worst_col:<10} | {status:<6}")

    print("-" * 80)
    if all_pass:
        print("ALL VARIABLES PASSED")
    else:
        print("SOME VARIABLES FAILED")
    return all_pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fortran_out", type=str)
    parser.add_argument("cpp_out", type=str)
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()

    ok = compare_outputs(args.fortran_out, args.cpp_out, args.tol)
    sys.exit(0 if ok else 1)
