#!/usr/bin/env python3
"""
compare.py
Differential verification comparison between Fortran WRF v4.6.1 reference
and CPPWRF C++ diffusion and damping outputs.
"""

import sys
import struct
import numpy as np

def read_chunks(filename):
    chunks = {}
    with open(filename, 'rb') as f:
        while True:
            tag_bytes = f.read(32)
            if not tag_bytes or len(tag_bytes) < 32:
                break
            tag = tag_bytes.decode('ascii', errors='ignore').strip('\x00').strip()
            count_bytes = f.read(4)
            if not count_bytes:
                break
            count = struct.unpack('i', count_bytes)[0]
            data_bytes = f.read(count * 8)
            arr = np.frombuffer(data_bytes, dtype=np.float64)
            chunks[tag] = arr
    return chunks

def main():
    if len(sys.argv) < 3:
        print("Usage: compare.py <fortran_out.bin> <cpp_out.bin> [--tol TOL]")
        sys.exit(1)

    fort_file = sys.argv[1]
    cpp_file = sys.argv[2]
    tol = 1e-12

    if '--tol' in sys.argv:
        idx = sys.argv.index('--tol')
        if idx + 1 < len(sys.argv):
            tol = float(sys.argv[idx + 1])

    f_chunks = read_chunks(fort_file)
    c_chunks = read_chunks(cpp_file)

    all_passed = True
    print("=" * 95)
    print(f"{'Variable / Kernel Tag':<30} | {'Count':<7} | {'Max Ref':<11} | {'Max Abs Err':<12} | {'Max Rel Err':<12} | {'Status'}")
    print("-" * 95)

    for tag in f_chunks:
        if tag not in c_chunks:
            print(f"{tag:<30} | MISSING IN C++ OUTPUT")
            all_passed = False
            continue

        f_data = f_chunks[tag]
        c_data = c_chunks[tag]

        if len(f_data) != len(c_data):
            print(f"{tag:<30} | SIZE MISMATCH: Fortran {len(f_data)} vs C++ {len(c_data)}")
            all_passed = False
            continue

        abs_err = np.abs(c_data - f_data)
        max_abs = np.max(abs_err)
        max_ref = np.max(np.abs(f_data))
        denom = np.maximum(1.0, np.abs(f_data))
        rel_err = abs_err / denom
        max_rel = np.max(rel_err)

        passed = (max_rel <= tol)
        status = "PASS" if passed else "FAIL"
        if not passed:
            all_passed = False

        print(f"{tag:<30} | {len(f_data):<7} | {max_ref:<11.4e} | {max_abs:<12.4e} | {max_rel:<12.4e} | {status}")

    print("=" * 95)
    print("6th-Order Monotonic Flux Limiter Verification (opt1 vs opt2 on sharp steps):")
    limiter_active = True
    for tag_base in ['u_6th', 'v_6th', 'w_6th', 't_6th', 'm_6th']:
        tag1 = f"{tag_base}_opt1"
        tag2 = f"{tag_base}_opt2"
        if tag1 in f_chunks and tag2 in f_chunks and tag1 in c_chunks and tag2 in c_chunks:
            f_diff = np.max(np.abs(f_chunks[tag1] - f_chunks[tag2]))
            c_diff = np.max(np.abs(c_chunks[tag1] - c_chunks[tag2]))
            print(f"  {tag_base:<8}: Fortran max |opt1 - opt2| = {f_diff:.4e}, C++ max |opt1 - opt2| = {c_diff:.4e}")
            if f_diff <= 0.0 or c_diff <= 0.0:
                print(f"  --> FAIL: Limiter did not activate for {tag_base}!")
                limiter_active = False
            else:
                print(f"  --> PASS: opt1 != opt2 verified (limiter active)")
    if not limiter_active:
        all_passed = False

    print("=" * 95)
    if all_passed:
        print(f"ALL TESTS PASSED within relative tolerance {tol:1.0e}!")
        sys.exit(0)
    else:
        print(f"SOME TESTS FAILED relative tolerance {tol:1.0e}!")
        sys.exit(1)

if __name__ == '__main__':
    main()
