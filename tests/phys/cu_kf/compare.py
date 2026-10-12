#!/usr/bin/env python3
"""
Differential comparison script for Kain-Fritsch (cu_physics = 1).
Compares output_fortran.bin against output_cpp.bin across all 16 variables.
PASS criterion: max rel diff <= 1e-12 in FP64.
"""

import sys
import os
import struct
import numpy as np

def compare_kf(f_file="output_fortran.bin", cpp_file="output_cpp.bin", tol=1e-12):
    if not os.path.exists(f_file):
        print(f"Error: {f_file} does not exist")
        return False
    if not os.path.exists(cpp_file):
        print(f"Error: {cpp_file} does not exist")
        return False

    # Read dimensions from input_kf.bin if available, else derive from file size
    ncols = 2500
    nlevels = 50
    if os.path.exists("input_kf.bin"):
        with open("input_kf.bin", "rb") as f:
            ncols, nlevels = struct.unpack("=ii", f.read(8))

    size_2d = ncols
    size_3d = ncols * nlevels

    var_meta = [
        # (name, desc, is_3d, floor)
        ("raincv", "Convective rain per timestep (mm)", False, 1e-5),
        ("pratec", "Convective rain rate (mm/s)", False, 1e-7),
        ("nca", "Relaxation counter (s)", False, 1.0),
        ("cutop", "Cloud top level (1-based)", False, 1.0),
        ("cubot", "Cloud bottom level (1-based)", False, 1.0),
        ("w0avg", "Running average W (m/s)", True, 1e-6),
        ("rthcuten", "Potential temp tendency (K/s)", True, 1e-4),
        ("rqvcuten", "Vapor tendency (kg/kg/s)", True, 1e-9),
        ("rqccuten", "Cloud water tendency (kg/kg/s)", True, 1e-9),
        ("rqrcuten", "Rain tendency (kg/kg/s)", True, 1e-9),
        ("rqicuten", "Ice tendency (kg/kg/s)", True, 1e-9),
        ("rqscuten", "Snow tendency (kg/kg/s)", True, 1e-9),
        ("cldfra_dp_kf", "Deep cloud fraction (0-1)", True, 1e-6),
        ("cldfra_sh_kf", "Shallow cloud fraction (0-1)", True, 1e-6),
        ("qc_kf", "In-cloud liquid (kg/kg)", True, 1e-9),
        ("qi_kf", "In-cloud ice (kg/kg)", True, 1e-9),
    ]

    with open(f_file, "rb") as ff, open(cpp_file, "rb") as fc:
        print(f"\n### Differential Comparison: Fortran vs C++ (tol = {tol:.1e})\n")
        print(f"| Variable | Description | Max Abs Diff | Max Rel Diff | Worst Column | Status |")
        print(f"|---|---|---|---|---|---|")

        all_pass = True
        for name, desc, is_3d, floor in var_meta:
            count = size_3d if is_3d else size_2d
            bytes_to_read = count * 8
            data_f = np.frombuffer(ff.read(bytes_to_read), dtype=np.float64)
            data_c = np.frombuffer(fc.read(bytes_to_read), dtype=np.float64)

            if len(data_f) != count or len(data_c) != count:
                print(f"Error: Incomplete read for {name} ({len(data_f)}, {len(data_c)} vs {count})")
                return False

            abs_diff = np.abs(data_c - data_f)
            denom = np.maximum(np.maximum(np.abs(data_f), np.abs(data_c)), floor)
            rel_diff = abs_diff / denom

            max_abs = np.max(abs_diff)
            max_rel = np.max(rel_diff)

            # Find worst column index
            worst_idx = np.argmax(rel_diff)
            worst_col = worst_idx % ncols if is_3d else worst_idx

            status = "**PASS**" if max_rel <= tol else "**FAIL**"
            if max_rel > tol:
                all_pass = False

            print(f"| `{name}` | {desc} | {max_abs:.6e} | {max_rel:.6e} | {worst_col} | {status} |")

        print()
        if all_pass:
            print(f"**Result: OVERALL PASS** (All variables max rel diff <= {tol:.1e})")
        else:
            print(f"**Result: OVERALL FAIL**")

        return all_pass

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Compare KF Eta scheme outputs.")
    parser.add_argument("fort_file", nargs="?", default="output_fortran.bin", help="Path to Fortran binary output")
    parser.add_argument("cpp_file", nargs="?", default="output_cpp.bin", help="Path to C++ binary output")
    parser.add_argument("--tol", type=float, default=1e-12, help="Relative difference tolerance")
    args, unknown = parser.parse_known_args()
    if unknown and len(unknown) == 1:
        try:
            args.tol = float(unknown[0])
        except ValueError:
            pass

    success = compare_kf(args.fort_file, args.cpp_file, args.tol)
    sys.exit(0 if success else 1)
