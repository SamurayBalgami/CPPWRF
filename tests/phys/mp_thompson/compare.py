#!/usr/bin/env python3
"""
Differential comparison script for Thompson microphysics (mp_physics = 8).
Compares output_fortran.bin against output_cpp.bin across all 19 variables.
PASS criterion: max rel diff <= 1e-12 in FP64.
"""

import sys
import os
import struct
import numpy as np

def compare_thompson(f_file="output_fortran.bin", cpp_file="output_cpp.bin", tol=1e-12):
    if not os.path.exists(f_file):
        print(f"Error: {f_file} does not exist")
        return False
    if not os.path.exists(cpp_file):
        print(f"Error: {cpp_file} does not exist")
        return False

    ncols = 2500
    nlevels = 50
    in_vars = {}
    if os.path.exists("input_thompson.bin"):
        with open("input_thompson.bin", "rb") as f:
            ncols, nlevels = struct.unpack("=ii", f.read(8))
            _dt = struct.unpack("=d", f.read(8))[0]
            size_2d = ncols
            size_3d = ncols * nlevels
            for name in ["th", "qv", "qc", "qr", "qi", "qs", "qg", "ni", "nr", "pii", "p", "w", "dz"]:
                in_vars[name] = np.frombuffer(f.read(size_3d * 8), dtype=np.float64)
            for name in ["rainnc", "rainncv", "snownc", "snowncv", "graupelnc", "graupelncv", "sr"]:
                in_vars[name] = np.frombuffer(f.read(size_2d * 8), dtype=np.float64)

    size_2d = ncols
    size_3d = ncols * nlevels

    var_meta = [
        # (name, desc, is_3d, floor)
        ("qv", "Water vapor mixing ratio (kg/kg)", True, 1e-9),
        ("qc", "Cloud water mixing ratio (kg/kg)", True, 1e-9),
        ("qr", "Rain water mixing ratio (kg/kg)", True, 1e-9),
        ("qi", "Cloud ice mixing ratio (kg/kg)", True, 1e-9),
        ("qs", "Snow mixing ratio (kg/kg)", True, 1e-9),
        ("qg", "Graupel mixing ratio (kg/kg)", True, 1e-9),
        ("ni", "Cloud ice number concentration (#/kg)", True, 1.0),
        ("nr", "Rain number concentration (#/kg)", True, 1.0),
        ("th", "Potential temperature (K)", True, 1.0),
        ("rainnc", "Accumulated total precipitation (mm)", False, 1e-4),
        ("rainncv", "Time-step total precipitation (mm)", False, 1e-4),
        ("snownc", "Accumulated snow+ice precipitation (mm)", False, 1e-4),
        ("snowncv", "Time-step snow+ice precipitation (mm)", False, 1e-4),
        ("graupelnc", "Accumulated graupel precipitation (mm)", False, 1e-4),
        ("graupelncv", "Time-step graupel precipitation (mm)", False, 1e-4),
        ("sr", "Snow ratio (fraction)", False, 1e-6),
        ("re_cloud", "Cloud water effective radius (m)", True, 1e-6),
        ("re_ice", "Cloud ice effective radius (m)", True, 1e-6),
        ("re_snow", "Snow effective radius (m)", True, 1e-6),
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
            if name in in_vars and len(in_vars[name]) == count:
                denom = np.maximum(denom, np.abs(in_vars[name]))
            rel_diff = abs_diff / denom

            max_abs = np.max(abs_diff)
            max_rel = np.max(rel_diff)

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
    parser = argparse.ArgumentParser(description="Compare Thompson microphysics outputs.")
    parser.add_argument("fort_file", nargs="?", default="output_fortran.bin", help="Path to Fortran binary output")
    parser.add_argument("cpp_file", nargs="?", default="output_cpp.bin", help="Path to C++ binary output")
    parser.add_argument("--tol", type=float, default=1e-12, help="Relative difference tolerance")
    args, unknown = parser.parse_known_args()
    if unknown and len(unknown) == 1:
        try:
            args.tol = float(unknown[0])
        except ValueError:
            pass

    ok = compare_thompson(args.fort_file, args.cpp_file, args.tol)
    sys.exit(0 if ok else 1)
