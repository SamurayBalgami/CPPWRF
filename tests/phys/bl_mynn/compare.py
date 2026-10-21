#!/usr/bin/env python3
"""compare.py WRF.bin CPP.bin NI NZ [--tol 0]: per-field max |difference| of the MYNN-EDMF outputs."""
import sys
import numpy as np
n3 = ("qke el_pbl sh3d sm3d exch_h exch_m tsq qsq cov qke_adv rublten rvblten rthblten rqvblten rqcblten rqiblten "
      "qc_bl qi_bl cldfra_bl edmf_a edmf_w edmf_qt edmf_thl edmf_ent edmf_qc").split()
n2 = "pblh maxwidth maxmf ztop_plume excess_h excess_q kpbl".split()
ni, nz = int(sys.argv[3]), int(sys.argv[4])
tol = float(sys.argv[sys.argv.index("--tol") + 1]) if "--tol" in sys.argv else 0.0
A = np.fromfile(sys.argv[1], dtype=np.float64)
B = np.fromfile(sys.argv[2], dtype=np.float64)
m3 = len(n3) * ni * nz
recs = [(k, A[i * ni * nz:(i + 1) * ni * nz], B[i * ni * nz:(i + 1) * ni * nz]) for i, k in enumerate(n3)]
recs += [(k, A[m3 + i * ni:m3 + (i + 1) * ni], B[m3 + i * ni:m3 + (i + 1) * ni]) for i, k in enumerate(n2)]
bad = 0
for k, x, y in recs:
    d = np.abs(x - y)
    d[np.isnan(x) & np.isnan(y)] = 0
    nd = int((d > tol).sum() + (np.isnan(x) ^ np.isnan(y)).sum())
    if nd:
        bad += 1
        j = int(np.nanargmax(d))
        print(f"DIFF {k:10s} {nd:7d} points, max {d[j]:.3e} at col {j % ni} k {j // ni + 1}: wrf {x[j]!r} cpp {y[j]!r}")
print(f"{'PASS' if bad == 0 else 'FAIL'}: {len(recs) - bad}/{len(recs)} fields identical")
sys.exit(1 if bad else 0)
