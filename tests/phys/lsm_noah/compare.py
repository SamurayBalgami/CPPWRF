#!/usr/bin/env python3
"""
Differential test comparison script for Noah LSM.
Compares Fortran vs C++ outputs on identical inputs.
Outputs max abs diff, max rel diff (with physical floor), worst column, and PASS/FAIL.
Acceptance criterion: max rel diff <= 1e-12 in FP64.
"""

import argparse
import struct
import sys
import numpy as np

def read_noah_output(filepath):
    with open(filepath, 'rb') as f:
        header = f.read(8)
        if len(header) < 8:
            raise ValueError(f"Corrupt output header in {filepath}")
        n_cols, nsoil = struct.unpack('=ii', header)
        n2d = n_cols
        n3d = n_cols * nsoil

        data = {}
        data['n_cols'] = n_cols
        data['nsoil'] = nsoil

        # 2D diagnostic variables (27)
        var2d_names = [
            'eta', 'sheat', 'eta_kinematic', 'fdown', 'ec', 'edir',
            'ett', 'esnow', 'drip', 'dew', 'beta', 'etp', 'ssoil',
            'flx1', 'flx2', 'flx3', 'snomlt', 'sncovr', 'runoff1',
            'runoff2', 'runoff3', 'rc', 'pc', 'rsmin', 'soilw', 'soilm', 'q1'
        ]
        for name in var2d_names:
            data[name] = np.fromfile(f, dtype=np.float64, count=n2d)

        # Updated states (10)
        state_names = [
            'cmc', 't1', 'snowh', 'sneqv', 'albedo', 'emissi',
            'z0', 'xlai', 'snotime1', 'ribb'
        ]
        for name in state_names:
            data[name] = np.fromfile(f, dtype=np.float64, count=n2d)

        # 3D outputs (5)
        var3d_names = ['stc', 'smc', 'sh2o', 'et', 'smav']
        for name in var3d_names:
            data[name] = np.fromfile(f, dtype=np.float64, count=n3d).reshape((nsoil, n_cols))

        return data

def compare_outputs(fortran_file, cpp_file, tol=1e-12):
    f_out = read_noah_output(fortran_file)
    c_out = read_noah_output(cpp_file)

    if f_out['n_cols'] != c_out['n_cols'] or f_out['nsoil'] != c_out['nsoil']:
        print(f"Dimension mismatch! Fortran: ({f_out['n_cols']}, {f_out['nsoil']}), C++: ({c_out['n_cols']}, {c_out['nsoil']})")
        return False, []

    n_cols = f_out['n_cols']
    nsoil = f_out['nsoil']

    floors = {
        'eta': 1.0,            # W/m2
        'sheat': 1.0,          # W/m2
        'eta_kinematic': 1e-9, # kg/(m2 s)
        'fdown': 1.0,          # W/m2
        'ec': 1e-9,            # kg/(m2 s)
        'edir': 1e-9,          # kg/(m2 s)
        'ett': 1e-9,           # kg/(m2 s)
        'esnow': 1e-9,         # kg/(m2 s)
        'drip': 1e-12,         # m/s
        'dew': 1e-12,          # m/s
        'beta': 1e-6,          # fraction
        'etp': 1.0,            # W/m2
        'ssoil': 1.0,          # W/m2
        'flx1': 1.0,           # W/m2
        'flx2': 1.0,           # W/m2
        'flx3': 1.0,           # W/m2
        'snomlt': 1e-12,       # m/s
        'sncovr': 1e-6,        # fraction
        'runoff1': 1e-12,      # m/s
        'runoff2': 1e-12,      # m/s
        'runoff3': 1e-12,      # m/s
        'rc': 1.0,             # s/m
        'pc': 0.01,            # fraction
        'rsmin': 1.0,          # s/m
        'soilw': 1e-6,         # fraction
        'soilm': 1e-6,         # fraction
        'q1': 1e-9,            # kg/kg
        'cmc': 1e-12,          # m
        't1': 1.0,             # K
        'snowh': 1e-12,        # m
        'sneqv': 1e-12,        # m
        'albedo': 0.01,        # fraction
        'emissi': 0.01,        # fraction
        'z0': 1e-6,            # m
        'xlai': 0.01,          # fraction
        'snotime1': 1.0,       # s
        'ribb': 1e-6,          # Richardson number
        'stc': 1.0,            # K
        'smc': 1e-6,           # m3/m3
        'sh2o': 1e-6,          # m3/m3
        'et': 1e-12,           # m/s
        'smav': 1e-6           # fraction
    }

    variables = [
        'eta', 'sheat', 'eta_kinematic', 'fdown', 'ec', 'edir',
        'ett', 'esnow', 'drip', 'dew', 'beta', 'etp', 'ssoil',
        'flx1', 'flx2', 'flx3', 'snomlt', 'sncovr', 'runoff1',
        'runoff2', 'runoff3', 'rc', 'pc', 'rsmin', 'soilw', 'soilm', 'q1',
        'cmc', 't1', 'snowh', 'sneqv', 'albedo', 'emissi',
        'z0', 'xlai', 'snotime1', 'ribb',
        'stc', 'smc', 'sh2o', 'et', 'smav'
    ]

    all_pass = True
    results = []

    print(f"\nComparing Fortran ({fortran_file}) vs C++ ({cpp_file}):")
    print(f"Grid: {n_cols} columns x {nsoil} soil layers, Tolerance: {tol:.1e}")
    print("-" * 80)
    print(f"{'Variable':<15} | {'Max Abs Diff':<14} | {'Max Rel Diff':<14} | {'Worst Col':<10} | {'Status':<6}")
    print("-" * 80)

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

        print(f"{var:<15} | {max_abs:<14.6e} | {max_rel:<14.6e} | {worst_col:<10} | {status:<6}")

    print("-" * 80)
    overall = "PASSED" if all_pass else "FAILED"
    print(f"Overall Noah LSM Verification: {overall}")
    return all_pass, results

def main():
    parser = argparse.ArgumentParser(description="Differential test comparator for Noah LSM.")
    parser.add_argument("fortran_out", help="Fortran output binary file")
    parser.add_argument("cpp_out", help="C++ output binary file")
    parser.add_argument("--tol", type=float, default=1e-12, help="Relative difference tolerance (default: 1e-12)")
    args = parser.parse_args()

    passed, _ = compare_outputs(args.fortran_out, args.cpp_out, args.tol)
    sys.exit(0 if passed else 1)

if __name__ == "__main__":
    main()
