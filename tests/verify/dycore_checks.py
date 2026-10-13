#!/usr/bin/env python3
"""Dycore consistency checks that need no external reference.

  dycore_checks.py <cppwrf_binary> <work_dir>

1. map_factor: squall2d (moist, Kessler, PD advection) with a constant map factor m = 2 and
   dx = 2000 m must equal m = 1, dx = 1000 m: the map-factor (ADT) form is invariant, and with
   m = 2 every scaling is exact in floating point, so the runs agree bit for bit. emdiv = 0
   because WRF's external-mode filter uses the grid dx (not invariant in WRF either).
2. coriolis: inertial oscillation on an f-plane (configs/inertial.ini), u = u0 cos(ft),
   v = -u0 sin(ft); tolerance 1e-3 m/s over one period (RK3 phase error ~ (f dt)^3 per step).
Exit 0 = all pass.
"""
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import verify as V  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def run(exe, cfg, out, *sets):
    cmd = [exe, os.path.join(ROOT, cfg), "--set", f"output.dir={out}"]
    for s in sets:
        cmd += ["--set", s]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def main(exe, work):
    os.makedirs(work, exist_ok=True)
    rep = V.Report("dycore consistency")
    a, b = os.path.join(work, "map1"), os.path.join(work, "map2")
    common = ["time.time_step_sound=4", "dyn.emdiv=0", "run.run_time=1800"]
    run(exe, "configs/squall2d.ini", a, *common)
    run(exe, "configs/squall2d.ini", b, *common, "domain.dx=2000", "domain.map_factor=2")
    worst = 0.0
    for s in V.list_steps(a):
        sa, sb = V.read_snap(a, s), V.read_snap(b, s)
        for k in ("th_pert", "u", "w", "mu", "qv", "qc", "qr", "p"):
            worst = max(worst, float(np.abs(sa[k] - sb[k]).max()))
    rep.check("map_factor_invariance", worst == 0.0, f"max |m=2,dx=2km - m=1,dx=1km| = {worst:.3e} (must be 0)")
    c = os.path.join(work, "inertial")
    run(exe, "configs/inertial.ini", c)
    err = 0.0
    for s in V.list_steps(c):
        sn = V.read_snap(c, s)
        t, f, u0 = sn["_time"], 1e-3, 10.0
        err = max(err, abs(sn["u"].mean() - u0 * np.cos(f * t)), abs(sn["v"].mean() + u0 * np.sin(f * t)))
    rep.check("coriolis_inertial", err <= 1e-3, f"max |(u,v) - exact| = {err:.2e} m/s over one period (tol 1e-3)")
    rep.print()
    return 1 if rep.failed else 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1], sys.argv[2]))
