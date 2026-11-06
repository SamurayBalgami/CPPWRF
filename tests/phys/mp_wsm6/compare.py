#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Differential test comparison script for WSM6 microphysics.
Compares Fortran vs C++ outputs on identical inputs.
Outputs max abs diff, max rel diff (with physical floor), worst column, and PASS/FAIL.
--tol 0 (CPU FP64): every value identical as tests/verify/exact.py (NaN equals only NaN, -0 differs from +0);
otherwise max rel diff <= tol.
"""
import os

import argparse
import struct
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import colcmp  # noqa: E402

def read_output(filepath):
    with open(filepath, 'rb') as f:
        header = f.read(8)
        if len(header) < 8:
            raise ValueError(f"Corrupt output header in {filepath}")
        n_cols, nz = struct.unpack('=ii', header)
        n3d = n_cols * nz
        n2d = n_cols

        th = np.fromfile(f, dtype=np.float64, count=n3d).reshape((nz, n_cols))
        q  = np.fromfile(f, dtype=np.float64, count=n3d).reshape((nz, n_cols))
        qc = np.fromfile(f, dtype=np.float64, count=n3d).reshape((nz, n_cols))
        qr = np.fromfile(f, dtype=np.float64, count=n3d).reshape((nz, n_cols))
        qi = np.fromfile(f, dtype=np.float64, count=n3d).reshape((nz, n_cols))
        qs = np.fromfile(f, dtype=np.float64, count=n3d).reshape((nz, n_cols))
        qg = np.fromfile(f, dtype=np.float64, count=n3d).reshape((nz, n_cols))

        rain       = np.fromfile(f, dtype=np.float64, count=n2d)
        rainncv    = np.fromfile(f, dtype=np.float64, count=n2d)
        snow       = np.fromfile(f, dtype=np.float64, count=n2d)
        snowncv    = np.fromfile(f, dtype=np.float64, count=n2d)
        graupel    = np.fromfile(f, dtype=np.float64, count=n2d)
        graupelncv = np.fromfile(f, dtype=np.float64, count=n2d)
        sr         = np.fromfile(f, dtype=np.float64, count=n2d)
        re = [np.fromfile(f, dtype=np.float64, count=n3d) for _ in range(3)]  # re_cloud, re_ice, re_snow
        if re[2].size != n3d or f.read(1):
            raise ValueError(f"{filepath}: not {10 * n3d + 7 * n2d} values after the header")
        re_cloud, re_ice, re_snow = (a.reshape((nz, n_cols)) for a in re)

        return {
            'n_cols': n_cols, 'nz': nz,
            'th': th, 'q': q, 'qc': qc, 'qr': qr, 'qi': qi, 'qs': qs, 'qg': qg,
            'rain': rain, 'rainncv': rainncv,
            'snow': snow, 'snowncv': snowncv,
            'graupel': graupel, 'graupelncv': graupelncv,
            'sr': sr,
            're_cloud': re_cloud, 're_ice': re_ice, 're_snow': re_snow
        }

def compare_outputs(fortran_file, cpp_file, tol=1e-12):
    f_out = read_output(fortran_file)
    c_out = read_output(cpp_file)

    if f_out['n_cols'] != c_out['n_cols'] or f_out['nz'] != c_out['nz']:
        print(f"Dimension mismatch! Fortran: ({f_out['n_cols']}, {f_out['nz']}), C++: ({c_out['n_cols']}, {c_out['nz']})")
        return False, []

    n_cols = f_out['n_cols']
    nz = f_out['nz']

    floors = {
        'th': 1.0,          # Potential temperature ~ 300 K
        'q': 1e-12,         # Water vapor mixing ratio (kg/kg)
        'qc': 1e-12,        # Cloud water mixing ratio (kg/kg)
        'qr': 1e-6,         # Rain water mixing ratio (kg/kg; physical floor for precipitating category)
        'qi': 1e-12,        # Cloud ice mixing ratio (kg/kg)
        'qs': 1e-6,         # Snow mixing ratio (kg/kg; physical floor for precipitating category)
        'qg': 1e-6,         # Graupel mixing ratio (kg/kg; physical floor for precipitating category)
        'rain': 1e-12,      # Accumulated rain (mm)
        'rainncv': 1e-12,   # Step rain (mm)
        'snow': 1e-12,      # Accumulated snow (mm)
        'snowncv': 1e-12,   # Step snow (mm)
        'graupel': 1e-12,   # Accumulated graupel (mm)
        'graupelncv': 1e-12,# Step graupel (mm)
        'sr': 1e-12,        # Snow ratio [0, 1]
        're_cloud': 1e-6, 're_ice': 1e-6, 're_snow': 1e-6  # effective radii (m)
    }

    variables = ['th', 'q', 'qc', 'qr', 'qi', 'qs', 'qg',
                 'rain', 'rainncv', 'snow', 'snowncv', 'graupel', 'graupelncv', 'sr', 're_cloud', 're_ice', 're_snow']
    all_pass = True
    results = []

    print(f"\nComparing Fortran ({fortran_file}) vs C++ ({cpp_file}):")
    print(f"Grid: {n_cols} columns x {nz} levels, Tolerance: {tol:.1e}")
    print("-" * 75)
    print(f"{'Variable':<12} | {'Max Abs Diff':<14} | {'Max Rel Diff':<14} | {'Worst Col':<10} | {'Status':<6}")
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

        if f_arr.ndim == 2:
            col_rel = np.max(rel_diff, axis=0)
            worst_col = int(np.argmax(col_rel))
        else:
            worst_col = int(np.argmax(rel_diff))

        ok = colcmp.same(var, f_arr, c_arr, quiet=True) if tol == 0 else max_rel <= tol
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_pass = False

        results.append({
            'var': var,
            'max_abs': max_abs,
            'max_rel': max_rel,
            'worst_col': worst_col,
            'status': status
        })

        print(f"{var:<12} | {max_abs:<14.6e} | {max_rel:<14.6e} | {worst_col:<10} | {status:<6}")

    print("-" * 75)
    overall = "PASS" if all_pass else "FAIL"
    print(f"Overall Result: {overall}")
    return all_pass, results

def main():
    parser = argparse.ArgumentParser(description='Compare Fortran and C++ WSM6 outputs')
    parser.add_argument('fortran_out', type=str, help='Fortran output binary')
    parser.add_argument('cpp_out', type=str, help='C++ output binary')
    parser.add_argument('--tol', type=float, default=1e-12, help='Relative tolerance (default: 1e-12)')
    args = parser.parse_args()

    passed, _ = compare_outputs(args.fortran_out, args.cpp_out, tol=args.tol)
    sys.exit(0 if passed else 1)

if __name__ == '__main__':
    main()
