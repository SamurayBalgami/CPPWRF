#!/usr/bin/env python3
"""Compare a cppwrf_real output file with WRF real.exe's: every variable's values (tolerance 0),
the variable list and order, dimensions, types, variable attributes and global attributes.
TITLE is ignored (it names the program). Exit status 1 on any difference.
Also used for wrfout/wrfrst files (tests/verify/run_wrf_files_test.sh); --ignore VAR,... skips the
values of the listed variables.

  compare_real.py [--ignore VAR,...] cppwrf_file wrf_file"""
import sys

import numpy as np
from netCDF4 import Dataset


def values(a, b, skip=()):
    bad = 0
    for n, v in a.variables.items():
        if n not in b.variables or n in skip:
            continue
        w = b.variables[n]
        if v.dtype.kind in 'SU' or n == 'Times' or n.startswith('md___'):
            x = b''.join(np.array(v[:]).ravel())
            y = b''.join(np.array(w[:]).ravel())
            if x != y:
                print('%-14s text differs' % n)
                bad += 1
            continue
        x = np.array(v[:], dtype=np.float64)
        y = np.array(w[:], dtype=np.float64)
        if x.shape != y.shape:
            print('%-14s shape %s vs %s' % (n, x.shape, y.shape))
            bad += 1
            continue
        nd = int((x != y).sum())
        if nd:
            bad += 1
            print('%-14s ndiff=%d/%d maxabs=%.3e' % (n, nd, x.size, np.abs(x - y).max()))
    return bad


def metadata(a, b, ignore):
    bad = 0
    va, vb = list(a.variables), list(b.variables)
    if va != vb:
        bad += 1
        print('variable list/order differs; missing', [x for x in vb if x not in va],
              'extra', [x for x in va if x not in vb])
    dims = lambda f: [(d, 'U' if f.dimensions[d].isunlimited() else len(f.dimensions[d])) for d in f.dimensions]
    if dims(a) != dims(b):
        bad += 1
        print('dimensions', dims(a), dims(b))
    for n in vb:
        if n not in va:
            continue
        x, y = a.variables[n], b.variables[n]
        if x.dimensions != y.dimensions or x.dtype != y.dtype:
            bad += 1
            print('dims/type', n, x.dimensions, y.dimensions, x.dtype, y.dtype)
        ax = [(k, str(x.getncattr(k))) for k in x.ncattrs()]
        ay = [(k, str(y.getncattr(k))) for k in y.ncattrs()]
        if ax != ay:
            bad += 1
            print('attributes', n, ax, ay)
    ga = [(k, repr(a.getncattr(k)), type(a.getncattr(k)).__name__) for k in a.ncattrs() if k not in ignore]
    gb = [(k, repr(b.getncattr(k)), type(b.getncattr(k)).__name__) for k in b.ncattrs() if k not in ignore]
    if ga != gb:
        bad += 1
        for x, y in zip(ga, gb):
            if x != y:
                print('global attribute', x, y)
        if len(ga) != len(gb):
            print('global attribute count', len(ga), len(gb))
    return bad


def main():
    args = sys.argv[1:]
    skip = set()
    if len(args) == 4 and args[0] == '--ignore':
        skip = set(args[1].split(','))
        args = args[2:]
    if len(args) != 2:
        sys.exit('usage: compare_real.py [--ignore VAR,...] cppwrf_file wrf_file')
    a, b = Dataset(args[0]), Dataset(args[1])
    nv = values(a, b, skip)
    nm = metadata(a, b, {'TITLE'})
    print('%s: %d variables differ, %d metadata differences' % (args[1], nv, nm))
    sys.exit(1 if nv or nm else 0)


if __name__ == '__main__':
    main()
