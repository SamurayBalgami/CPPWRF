#!/usr/bin/env python3
"""
Differential test comparison script for Kessler microphysics.
Compares Fortran vs C++ outputs on identical inputs.
Outputs max abs diff, max rel diff (with physical floor), worst column, and PASS/FAIL.
Acceptance criterion: max rel diff <= 1e-12 in FP64.
"""

import argparse
import struct
import sys
import numpy as np

def read_output(filepath):
    with open(filepath, 'rb') as f:
        header = f.read(8)
        if len(header) < 8:
            raise ValueError(f"Corrupt output header in {filepath}")
        n_cols, nz = struct.unpack('=ii', header)
        n3d = n_cols * nz
        n2d = n_cols

        t = np.fromfile(f, dtype=np.float64, count=n3d).reshape((nz, n_cols))
        qv = np.fromfile(f, dtype=np.float64, count=n3d).reshape((nz, n_cols))
        qc = np.fromfile(f, dtype=np.float64, count=n3d).reshape((nz, n_cols))
        qr = np.fromfile(f, dtype=np.float64, count=n3d).reshape((nz, n_cols))
        rainnc = np.fromfile(f, dtype=np.float64, count=n2d)
        rainncv = np.fromfile(f, dtype=np.float64, count=n2d)

        return {
            'n_cols': n_cols,
            'nz': nz,
            't': t,
            'qv': qv,
            'qc': qc,
            'qr': qr,
            'rainnc': rainnc,
            'rainncv': rainncv
        }

def compare_outputs(fortran_file, cpp_file, tol=1e-12):
    f_out = read_output(fortran_file)
    c_out = read_output(cpp_file)

    if f_out['n_cols'] != c_out['n_cols'] or f_out['nz'] != c_out['nz']:
        print(f"Dimension mismatch! Fortran: ({f_out['n_cols']}, {f_out['nz']}), C++: ({c_out['n_cols']}, {c_out['nz']})")
        return False, []

    n_cols = f_out['n_cols']
    nz = f_out['nz']

    # Floors for relative difference calculation
    floors = {
        't': 1.0,           # Potential temperature ~ 300 K
        'qv': 1e-12,        # Water vapor mixing ratio (kg/kg)
        'qc': 1e-12,        # Cloud water mixing ratio (kg/kg)
        'qr': 1e-12,        # Rain water mixing ratio (kg/kg)
        'rainnc': 1e-12,    # Accumulated rain (mm)
        'rainncv': 1e-12    # Step rain (mm)
    }

    variables = ['t', 'qv', 'qc', 'qr', 'rainnc', 'rainncv']
    all_pass = True
    results = []

    print(f"\nComparing Fortran ({fortran_file}) vs C++ ({cpp_file}):")
    print(f"Grid: {n_cols} columns x {nz} levels, Tolerance: {tol:.1e}")
    print("-" * 75)
    print(f"{'Variable':<10} | {'Max Abs Diff':<14} | {'Max Rel Diff':<14} | {'Worst Col':<10} | {'Status':<6}")
    print("-" * 75)

    for var in variables:
        f_arr = f_out[var]
        c_arr = c_out[var]
        floor = floors[var]

        abs_diff = np.abs(f_arr - c_arr)
        max_abs = np.max(abs_diff)

        scale = np.maximum(np.maximum(np.abs(f_arr), np.abs(c_arr)), floor)
        rel_diff = abs_diff / scale
        max_rel = np.max(rel_diff)

        # Find column of worst relative diff
        if f_arr.ndim == 2:
            # max over vertical levels per column
            col_rel = np.max(rel_diff, axis=0)
            worst_col = int(np.argmax(col_rel))
        else:
            worst_col = int(np.argmax(rel_diff))

        status = "PASS" if max_rel <= tol else "FAIL"
        if max_rel > tol:
            all_pass = False

        results.append({
            'var': var,
            'max_abs': max_abs,
            'max_rel': max_rel,
            'worst_col': worst_col,
            'status': status
        })

        print(f"{var:<10} | {max_abs:<14.6e} | {max_rel:<14.6e} | {worst_col:<10} | {status:<6}")

    print("-" * 75)
    overall = "ALL TESTS PASSED" if all_pass else "TESTS FAILED"
    print(f"Overall: {overall}\n")
    return all_pass, results

def main():
    parser = argparse.ArgumentParser(description='Compare Fortran and C++ microphysics outputs')
    parser.add_argument('fortran_out', nargs='?', default='kessler_output_fortran.bin', help='Fortran binary output')
    parser.add_argument('cpp_out', nargs='?', default='kessler_output_cpp.bin', help='C++ binary output')
    parser.add_argument('--tol', type=float, default=1e-12, help='Relative tolerance for PASS')
    args = parser.parse_args()

    passed, _ = compare_outputs(args.fortran_out, args.cpp_out, tol=args.tol)
    sys.exit(0 if passed else 1)

if __name__ == '__main__':
    main()
