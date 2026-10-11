#!/usr/bin/env python3
"""Self-test for verify.py: writes synthetic outputs in the I/O-contract format and
checks that good data PASSes and broken data FAILs. Run: python3 tests/verify/selftest.py"""
import json
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import verify  # noqa: E402


def write_run(d, prec="float64", **extra):
    os.makedirs(d, exist_ok=True)
    run = {"case": "synthetic", "precision": prec, "backend": "cpu", "grid": {"dx": 100.0}}
    run.update(extra)
    with open(os.path.join(d, "run.json"), "w") as f:
        json.dump(run, f)


def write_diag(d, rows):
    cols = ["step", "time", "max_abs_u", "max_abs_v", "max_abs_w", "min_th_pert",
            "max_th_pert", "dry_mass", "theta_mass", "total_energy", "wall_s"]
    with open(os.path.join(d, "diag.csv"), "w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(f"{r.get(c, 0.0):.17g}" for c in cols) + "\n")


def write_snap(d, step, time, fields, dtype="float64"):
    dt = np.dtype("<f8" if dtype == "float64" else "<f4")
    meta = {"step": step, "time": time, "variables": []}
    off = 0
    with open(os.path.join(d, f"snap_{step:06d}.bin"), "wb") as f:
        for name, (stag, a) in fields.items():
            b = np.ascontiguousarray(a, dtype=dt).tobytes()
            f.write(b)
            meta["variables"].append({"name": name, "stagger": stag, "shape": list(a.shape),
                                      "dtype": dtype, "offset": off, "nbytes": len(b)})
            off += len(b)
    with open(os.path.join(d, f"snap_{step:06d}.json"), "w") as f:
        json.dump(meta, f)


def straka_like(nx=512, nz=64, dx=100.0, front=15540.0, tmin=-9.77):
    """Symmetric cold pool: theta' = tmin at the centre, crossing -1 K exactly at +-front."""
    x = (np.arange(nx) + 0.5) * dx
    xc = nx * dx / 2
    r = np.abs(x - xc) / front
    row = tmin * (1 - r ** 2)          # = -1 K where r^2 = 1 + 1/tmin
    row = np.where(row < 0, row, 0.0)
    # rescale so that -1 K is reached exactly at r = 1
    row = np.where(r < 1, -1 + (tmin + 1) * (1 - r ** 2), np.where(r < 1.2, -1 * (1.2 - r) / 0.2, 0.0))
    th = np.zeros((nz, 1, nx))
    th[0, 0, :] = row
    return th


def fields_for(th, nx, nz, u_sym=True):
    u = np.zeros((nz, 1, nx + 1))
    u[:, :, : (nx + 1) // 2] = 1.0
    u[:, :, (nx + 1) // 2 + 1:] = -1.0       # antisymmetric -> mirror symmetric flow
    if not u_sym:
        u[3, 0, 10] += 0.5
    w = np.zeros((nz + 1, 1, nx))
    return {"u": ("U", u), "w": ("W", w), "th_pert": ("M", th)}


def expect(name, rep, should_fail, failures):
    ok = rep.failed == should_fail
    print(f"[{'ok' if ok else 'WRONG'}] {name}: verify {'FAILED' if rep.failed else 'PASSED'}"
          f" (expected {'FAIL' if should_fail else 'PASS'})")
    if not ok:
        failures.append(name)


def main():
    failures = []
    tmp = tempfile.mkdtemp(prefix="cppwrf_selftest_")
    nx, nz = 512, 64
    diag_ok = [{"step": 0, "time": 0, "dry_mass": 1e10, "theta_mass": 3e12},
               {"step": 1800, "time": 900, "dry_mass": 1e10, "theta_mass": 3e12}]

    # 1. Good Straka-like output must PASS.
    d = os.path.join(tmp, "straka_good")
    write_run(d)
    write_diag(d, diag_ok)
    init = np.zeros((nz, 1, nx))
    init[30, 0, nx // 2 - 20: nx // 2 + 20] = -15.0
    write_snap(d, 0, 0.0, fields_for(init, nx, nz))
    write_snap(d, 1800, 900.0, fields_for(straka_like(), nx, nz))
    rep = verify.run_check("straka", d)
    expect("straka good", rep, False, failures)
    m = verify.straka_metrics(d)
    if abs(m["front"] - 15540.0) > 60.0:
        print(f"[WRONG] front metric {m['front']:.1f} m, expected 15540 +- 60")
        failures.append("front metric")

    # 2. Front too far -> FAIL.
    d = os.path.join(tmp, "straka_front_bad")
    write_run(d)
    write_diag(d, diag_ok)
    write_snap(d, 0, 0.0, fields_for(init, nx, nz))
    write_snap(d, 1800, 900.0, fields_for(straka_like(front=17000.0), nx, nz))
    expect("straka wrong front", verify.run_check("straka", d), True, failures)

    # 3. Asymmetric u -> FAIL.
    d = os.path.join(tmp, "straka_asym")
    write_run(d)
    write_diag(d, diag_ok)
    write_snap(d, 0, 0.0, fields_for(init, nx, nz))
    write_snap(d, 1800, 900.0, fields_for(straka_like(), nx, nz, u_sym=False))
    expect("straka asymmetric", verify.run_check("straka", d), True, failures)

    # 4. NaN in a snapshot -> FAIL.
    d = os.path.join(tmp, "nan")
    write_run(d)
    write_diag(d, diag_ok)
    bad = straka_like()
    bad[5, 0, 7] = np.nan
    write_snap(d, 0, 0.0, fields_for(init, nx, nz))
    write_snap(d, 1800, 900.0, fields_for(bad, nx, nz))
    expect("nan snapshot", verify.run_check("straka", d), True, failures)

    # 5. Mass drift -> FAIL.
    d = os.path.join(tmp, "drift")
    write_run(d)
    write_diag(d, [diag_ok[0], {"step": 1800, "time": 900, "dry_mass": 1.000001e10, "theta_mass": 3e12}])
    write_snap(d, 0, 0.0, fields_for(init, nx, nz))
    write_snap(d, 1800, 900.0, fields_for(straka_like(), nx, nz))
    expect("mass drift", verify.run_check("straka", d), True, failures)

    # 6. Rest state: quiet PASS, noisy FAIL; FP32 tolerance used.
    for name, w, fail in (("rest quiet", 1e-12, False), ("rest noisy", 1e-2, True)):
        d = os.path.join(tmp, name.replace(" ", "_"))
        write_run(d, prec="float32")
        rows = [dict(r, max_abs_w=w) for r in diag_ok]
        write_diag(d, rows)
        z = np.zeros((nz, 1, nx))
        write_snap(d, 0, 0.0, {"th_pert": ("M", z)}, "float32")
        expect(name, verify.run_check("rest_state", d), fail, failures)

    # 7. Compare identical vs perturbed.
    a, b = os.path.join(tmp, "straka_good"), os.path.join(tmp, "straka_asym")
    expect("compare identical", verify.run_compare(a, a, 1e-12), False, failures)
    expect("compare different", verify.run_compare(a, b, 1e-6), True, failures)

    print("\nSELFTEST", "FAILED: " + ", ".join(failures) if failures else "PASSED", f"(data in {tmp})")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
