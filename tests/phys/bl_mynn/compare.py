#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""compare.py WRF.bin CPP.bin NI NZ: the MYNN-EDMF outputs, every value identical as tests/verify/exact.py (NaN equals
only NaN, -0 differs from +0)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import colcmp  # noqa: E402

n3 = ("qke el_pbl sh3d sm3d exch_h exch_m tsq qsq cov qke_adv rublten rvblten rthblten rqvblten rqcblten rqiblten "
      "qc_bl qi_bl cldfra_bl edmf_a edmf_w edmf_qt edmf_thl edmf_ent edmf_qc").split()
n2 = "pblh maxwidth maxmf ztop_plume excess_h excess_q kpbl".split()
ni, nz = int(sys.argv[3]), int(sys.argv[4])
size = (len(n3) * nz + len(n2)) * ni
A = colcmp.read_exact(sys.argv[1], np.float64, size)
B = colcmp.read_exact(sys.argv[2], np.float64, size)
for p in sys.argv[1:3]:
    if os.path.getsize(p) != size * 8:
        sys.exit(f"{p}: {os.path.getsize(p)} bytes, expected {size * 8}")
m3 = len(n3) * ni * nz
recs = [(k, A[i * ni * nz:(i + 1) * ni * nz], B[i * ni * nz:(i + 1) * ni * nz]) for i, k in enumerate(n3)]
recs += [(k, A[m3 + i * ni:m3 + (i + 1) * ni], B[m3 + i * ni:m3 + (i + 1) * ni]) for i, k in enumerate(n2)]
bad = sum(not colcmp.same(k, x, y, quiet=True) for k, x, y in recs)
print(f"{'PASS' if bad == 0 else 'FAIL'}: {len(recs) - bad}/{len(recs)} fields identical")
sys.exit(1 if bad else 0)
