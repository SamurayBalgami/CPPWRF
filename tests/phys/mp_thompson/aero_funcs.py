#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Samples for the aerosol functions of the Thompson scheme (aero_funcs.f90, aero_funcs.cpp) and their comparison
at tolerance 0.   aero_funcs.py gen FILE N  |  aero_funcs.py cmp FORTRAN_OUT CPP_OUT N"""
import sys
import numpy as np

if sys.argv[1] == "gen":
    n = int(sys.argv[3])
    r = np.random.RandomState(11)
    cols = [r.uniform(228.0, 315.0, n),                     # activ_ncloud: Tt
            np.where(r.uniform(size=n) < 0.1, r.uniform(-2.0, 0.01, n), 10.0 ** r.uniform(-2.5, 2.3, n)),  # Ww
            10.0 ** r.uniform(6.0, 10.3, n),                 # NCCN [#/m3]
            10.0 ** r.uniform(-4.5, -2.3, n),                # Eff_aero: D
            r.choice([0.04e-6, 0.8e-6], n),                  # Da
            r.uniform(1.2e-5, 1.9e-5, n),                    # visc
            r.uniform(0.25, 1.35, n),                        # rhoa
            r.uniform(210.0, 310.0, n),                      # Temp
            r.randint(0, 3, n).astype(float),                # species r, s, g
            r.uniform(-45.0, -0.01, n),                      # iceDeMott: tempc
            r.uniform(0.25, 1.4, n),                         # rho
            10.0 ** r.uniform(2.0, 9.0, n),                  # nifa [#/m3]
            r.uniform(195.0, 238.0, n),                      # iceKoop: temp
            r.uniform(1.0e-7, 2.0e-4, n),                    # qv
            np.zeros(n),                                     # qvs (below)
            10.0 ** r.uniform(6.0, 10.0, n),                 # naero
            np.full(n, 60.0)]                                # dt
    cols[14] = cols[13] / r.uniform(1.0, 2.0, n)  # qv / qvs from 1 to 2
    with open(sys.argv[2], "wb") as f:
        f.write(np.int32(n).tobytes())
        for c in cols:
            f.write(np.ascontiguousarray(c, dtype=np.float64).tobytes())
else:
    n = int(sys.argv[4])
    a = np.fromfile(sys.argv[2], dtype=np.float64).reshape(4, n)
    b = np.fromfile(sys.argv[3], dtype=np.float64).reshape(4, n)
    ok = True
    for i, name in enumerate(("activ_ncloud", "Eff_aero", "iceDeMott", "iceKoop")):
        bad = ~((a[i] == b[i]) & (np.signbit(a[i]) == np.signbit(b[i])) | (np.isnan(a[i]) & np.isnan(b[i])))
        nb = int(bad.sum())
        ok &= nb == 0
        first = f" (first {a[i][bad][0]!r} vs {b[i][bad][0]!r} at {int(np.argmax(bad))})" if nb else ""
        print(f"{name}: {'PASS' if nb == 0 else 'FAIL'} {nb}/{n} differ{first}")
    sys.exit(0 if ok else 1)
