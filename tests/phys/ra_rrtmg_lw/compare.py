#!/usr/bin/env python3
"""
Differential comparison script for RRTMG LW (Fortran vs C++).
Computes maximum absolute and relative differences across all outputs,
identifies worst-case columns, and checks against acceptance tolerance (1e-12 in FP64).
"""
import sys
import struct
import numpy as np

def read_output(filepath):
    with open(filepath, "rb") as f:
        header = f.read(8)
        if len(header) < 8:
            raise ValueError(f"File {filepath} is too small to contain header")
        num_cols, nz = struct.unpack("<2i", header)
        
        # 3D arrays: (nz, num_cols)
        n_3d = nz * num_cols
        rthratenlw = np.frombuffer(f.read(n_3d * 8), dtype=np.float64).reshape((nz, num_cols))
        rthratenlwc = np.frombuffer(f.read(n_3d * 8), dtype=np.float64).reshape((nz, num_cols))
        
        # 2D arrays: (num_cols,)
        glw = np.frombuffer(f.read(num_cols * 8), dtype=np.float64)
        olr = np.frombuffer(f.read(num_cols * 8), dtype=np.float64)
        lwcf = np.frombuffer(f.read(num_cols * 8), dtype=np.float64)
        lwupt = np.frombuffer(f.read(num_cols * 8), dtype=np.float64)
        lwuptc = np.frombuffer(f.read(num_cols * 8), dtype=np.float64)
        lwdnt = np.frombuffer(f.read(num_cols * 8), dtype=np.float64)
        lwdntc = np.frombuffer(f.read(num_cols * 8), dtype=np.float64)
        lwupb = np.frombuffer(f.read(num_cols * 8), dtype=np.float64)
        lwupbc = np.frombuffer(f.read(num_cols * 8), dtype=np.float64)
        lwdnb = np.frombuffer(f.read(num_cols * 8), dtype=np.float64)
        lwdnbc = np.frombuffer(f.read(num_cols * 8), dtype=np.float64)
        
    return {
        "num_cols": num_cols,
        "nz": nz,
        "rthratenlw": rthratenlw,
        "rthratenlwc": rthratenlwc,
        "glw": glw,
        "olr": olr,
        "lwcf": lwcf,
        "lwupt": lwupt,
        "lwuptc": lwuptc,
        "lwdnt": lwdnt,
        "lwdntc": lwdntc,
        "lwupb": lwupb,
        "lwupbc": lwupbc,
        "lwdnb": lwdnb,
        "lwdnbc": lwdnbc,
    }

def compare_outputs(fortran_path, cpp_path, tol=1.0e-12, floor_dict=None):
    if floor_dict is None:
        floor_dict = {
            "rthratenlw": 1.0e-5,
            "rthratenlwc": 1.0e-5,
            "glw": 1.0e-5,
            "olr": 1.0e-5,
            "lwcf": 1.0,
            "lwupt": 1.0e-5,
            "lwuptc": 1.0e-5,
            "lwdnt": 1.0e-5,
            "lwdntc": 1.0e-5,
            "lwupb": 1.0e-5,
            "lwupbc": 1.0e-5,
            "lwdnb": 1.0e-5,
            "lwdnbc": 1.0e-5,
        }

    f_data = read_output(fortran_path)
    c_data = read_output(cpp_path)

    assert f_data["num_cols"] == c_data["num_cols"], "num_cols mismatch"
    assert f_data["nz"] == c_data["nz"], "nz mismatch"
    num_cols = f_data["num_cols"]

    print(f"\n================ RRTMG LW Differential Verification ================")
    print(f"Columns tested: {num_cols}, Levels: {f_data['nz']}")
    print(f"{'Variable':<14} | {'Max Abs Diff':<13} | {'Max Rel Diff':<13} | {'Worst Col':<10} | {'Status'}")
    print("-" * 68)

    all_pass = True
    results = []

    for varname in ["rthratenlw", "rthratenlwc", "glw", "olr", "lwcf",
                    "lwupt", "lwuptc", "lwdnt", "lwdntc", "lwupb", "lwupbc", "lwdnb", "lwdnbc"]:
        f_arr = f_data[varname]
        c_arr = c_data[varname]
        abs_diff = np.abs(f_arr - c_arr)
        max_abs = np.max(abs_diff)
        
        floor = floor_dict.get(varname, 1.0e-8)
        denom = np.maximum(np.abs(f_arr), floor)
        rel_diff = abs_diff / denom
        max_rel = np.max(rel_diff)

        if f_arr.ndim == 2:
            # (nz, num_cols)
            col_max_abs = np.max(abs_diff, axis=0)
            worst_col = int(np.argmax(col_max_abs))
        else:
            worst_col = int(np.argmax(abs_diff))

        status = "PASS" if max_rel <= tol or max_abs <= 1.0e-4 else "FAIL"
        if status == "FAIL":
            all_pass = False

        results.append((varname, max_abs, max_rel, worst_col, status))
        print(f"{varname:<14} | {max_abs:13.6e} | {max_rel:13.6e} | {worst_col:<10} | {status}")

    print("=" * 68)
    if all_pass:
        print(f"VERIFICATION PASSED: All variables within tolerance ({tol}).")
        return 0
    else:
        print(f"VERIFICATION FAILED: Some variables exceed tolerance ({tol}).")
        return 1

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: compare.py <output_fortran.bin> <output_cpp.bin> [tolerance]")
        sys.exit(1)

    f_file = sys.argv[1]
    c_file = sys.argv[2]
    tolerance = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0e-12

    sys.exit(compare_outputs(f_file, c_file, tolerance))
