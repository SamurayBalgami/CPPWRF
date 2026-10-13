#!/usr/bin/env python3
"""CPPWRF real-data run vs a WRF wrfout with the same input.

  compare_real.py <wrfout> <cppwrf_out_dir> [--interior N]

Per output time: max abs / relative L2 difference of U, V, W, PH (full), MU (full), P (full),
THM (theta_m), QVAPOR, over the whole domain and over the interior (N points away from the
lateral boundaries, default 6). Exit 0 always (a report, the tolerances are set by the caller).
"""
import os
import sys

import netCDF4 as nc
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import verify as V  # noqa: E402


def main(wrfout, d, inner):
    ds = nc.Dataset(wrfout)
    steps = V.list_steps(d)
    xt = np.array(ds.variables["XTIME"][:], dtype=float) * 60.0
    snaps = {V.read_snap(d, s)["_time"]: s for s in steps}
    print(f"{'t':>6} {'var':>6} {'max|d|':>11} {'relL2':>10} {'relL2_in':>10} {'max|wrf|':>10}")
    for it, t in enumerate(xt):
        if t not in snaps:
            continue
        sn = V.read_snap(d, snaps[t])
        g = lambda n: np.array(ds.variables[n][it], dtype=float)  # noqa: E731
        pairs = {
            "U": (g("U"), sn["u"]), "V": (g("V"), sn["v"]), "W": (g("W"), sn["w"]),
            "PH": (g("PH") + g("PHB"), sn["ph"]), "MU": (g("MU") + g("MUB"), sn["mu"][0]),
            "P": (g("P") + g("PB"), sn["p"]), "QV": (g("QVAPOR"), sn["qv"]),
        }
        if "THM" in ds.variables and "theta_m" in sn:
            pairs["THM"] = (g("THM") + 300.0, sn["theta_m"])
        for k, (a, b) in pairs.items():
            e = b - a
            sl = (Ellipsis, slice(inner, -inner), slice(inner, -inner))
            rl = np.sqrt(np.sum(e * e) / max(np.sum(a * a), 1e-300))
            ri = np.sqrt(np.sum(e[sl] ** 2) / max(np.sum(a[sl] ** 2), 1e-300))
            print(f"{t:6.0f} {k:>6} {np.abs(e).max():11.3e} {rl:10.2e} {ri:10.2e} {np.abs(a).max():10.3e}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    inner = int(sys.argv[sys.argv.index("--interior") + 1]) if "--interior" in sys.argv else 6
    main(sys.argv[1], sys.argv[2], inner)
