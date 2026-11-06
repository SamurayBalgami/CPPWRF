# SPDX-License-Identifier: GPL-3.0-or-later
"""Tolerance-0 comparison of the differential tests' outputs (WRF's Fortran vs the C++), as tests/verify/exact.py:
NaN equals only NaN, -0 differs from +0. Used by the compare scripts of tests/phys, tests/dyn when --tol is 0."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verify"))
import exact  # noqa: E402,F401


def same(name, f, c, quiet=False):
    """True when the arrays f (WRF) and c (C++) hold the same values; prints the field's line otherwise."""
    f, c = np.asarray(f, dtype=np.float64), np.asarray(c, dtype=np.float64)
    if f.size == 0:
        print(f"{name}: FAIL (no values)")
        return False
    d = exact.describe(f, c)
    if d is None:
        if not quiet:
            print(f"{name}: identical ({f.size} values)")
        return True
    print(f"{name}: FAIL {d}")
    return False


def read_exact(path, dtype, count, offset=0):
    """`count` values of `dtype` from `path` at byte `offset`; an error when the file holds fewer."""
    a = np.fromfile(path, dtype=dtype, count=count, offset=offset)
    if a.size != count:
        raise ValueError(f"{path}: {a.size} of {count} values")
    return a
