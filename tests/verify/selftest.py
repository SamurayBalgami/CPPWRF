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
    tmp = tempfile.mkdtemp(prefix="cppwrf_selftest_", dir=os.environ.get("CPPWRF_SCRATCH"))
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
    c = os.path.join(tmp, "straka_no_w")  # a variable the other run does not write
    write_run(c)
    write_diag(c, diag_ok)
    f = fields_for(straka_like(), nx, nz)
    del f["w"]
    write_snap(c, 1800, 900.0, f)
    expect("compare missing variable", verify.run_compare(a, c, 1e-6), True, failures)

    exact_faults(tmp, failures)

    print("\nSELFTEST", "FAILED: " + ", ".join(failures) if failures else "PASSED", f"(data in {tmp})")
    return 1 if failures else 0


def write_wrfout(path, times, fault=None, drop=None):
    """A small wrfout-like file (Times, U, W with WRF's -0, T, QVAPOR; global and variable attributes);
    fault(name, array, frame) edits a value array before it is written, fault('(file)', dataset, None) the file;
    drop: a variable left out."""
    import netCDF4
    ds = netCDF4.Dataset(path, "w")
    ds.createDimension("Time", None)
    ds.createDimension("DateStrLen", 19)
    for n, k in (("bottom_top", 4), ("bottom_top_stag", 5), ("south_north", 3), ("west_east", 5), ("west_east_stag", 6)):
        ds.createDimension(n, k)
    ds.setncatts({"TITLE": "test", "DX": np.float32(1000.0), "MP_PHYSICS": np.int32(6)})
    rng = np.random.default_rng(1)
    vs = {"U": ("bottom_top", "south_north", "west_east_stag"), "W": ("bottom_top_stag", "south_north", "west_east"),
          "T": ("bottom_top", "south_north", "west_east"), "QVAPOR": ("bottom_top", "south_north", "west_east")}
    ts = ds.createVariable("Times", "S1", ("Time", "DateStrLen"))
    for n, dims in vs.items():
        if n != drop:
            ds.createVariable(n, "f8", ("Time",) + dims).setncatts({"MemoryOrder": "XYZ", "units": "m s-1"})
    for it, t in enumerate(times):
        ts[it] = netCDF4.stringtochar(np.array([t], "S19"))
        for n, dims in vs.items():
            a = rng.standard_normal([len(ds.dimensions[x]) for x in dims])
            if n == "W":
                a[0] = -0.0
            if fault:
                fault(n, a, it)
            if n != drop:
                ds[n][it] = a
    if fault:
        fault("(file)", ds, None)
    ds.close()


def exact_faults(tmp, failures):
    """Tolerance-0 comparison (exact.py, used by wrf_suite.py and tools/gcp/matrix_variant.py): identical files
    are identical, each deliberate fault in the CPPWRF-side copy is a difference."""
    import exact
    T = ["2026-10-12_00:00:00", "2026-10-12_00:01:00", "2026-10-12_00:02:00"]
    T4 = T + ["2026-10-12_00:03:00"]

    def at(name, f, it=1, idx=(1, 1, 2)):
        return lambda n, a, i: f(a, idx) if (n, i) == (name, it) else None

    def ulp(a, i):
        a[i] = np.nextafter(a[i], np.inf)

    def nan(a, i):
        a[i] = np.nan

    def nan_and_one(n, a, i):
        if (n, i) == ("QVAPOR", 2):
            a[0, 0, 0] = np.nan
            a[2, 2, 2] += 1.0

    def plus_zero(a, i):
        a[0] = 0.0

    def gatt(n, ds, i):
        if n == "(file)":
            ds.MP_PHYSICS = np.int32(8)

    def vatt(n, ds, i):
        if n == "(file)":
            ds["T"].units = "K"

    # name, CPPWRF fault, CPPWRF frames, CPPWRF variable left out, a difference expected, subset (suite references)
    cases = [("identical", None, T, None, False, False), ("one ulp in T", at("T", ulp), T, None, True, False),
             ("NaN for a number", at("U", nan), T, None, True, False),
             ("NaN and 1.0 elsewhere", nan_and_one, T, None, True, False),
             ("+0 for WRF's -0", at("W", plus_zero, 0), T, None, True, False),
             ("a frame missing", None, T[:2], None, True, False), ("an extra frame", None, T4, None, True, False),
             ("a variable missing", None, T, "QVAPOR", True, False),
             ("a global attribute", gatt, T, None, True, False), ("a variable attribute", vatt, T, None, True, False),
             ("subset: identical, CPPWRF longer", None, T4, None, False, True),
             ("subset: one ulp", at("T", ulp, 2), T, None, True, True),
             ("subset: final frame missing", None, T[:2], None, True, True)]
    for name, fault, times, drop, should_differ, subset in cases:
        d = os.path.join(tmp, "exact_" + "".join(c if c.isalnum() else "_" for c in name))
        os.makedirs(d, exist_ok=True)
        write_wrfout(os.path.join(d, "wrf"), T)
        write_wrfout(os.path.join(d, "cpp"), times, fault, drop)
        bad, n = exact.compare_stream([os.path.join(d, "cpp")], [os.path.join(d, "wrf")], subset=subset,
                                      until=T[-1] if subset else None)
        ok = bool(bad) == should_differ
        print(f"[{'ok' if ok else 'WRONG'}] exact {name}: {'differs ' + str(bad)[:160] if bad else 'identical'}"
              f" (expected {'a difference' if should_differ else 'identical'})")
        if not ok:
            failures.append("exact " + name)
    d = os.path.join(tmp, "exact_wrf_variable_missing")  # a variable WRF writes and CPPWRF does not, reversed
    os.makedirs(d, exist_ok=True)
    write_wrfout(os.path.join(d, "wrf"), T, drop="QVAPOR")
    write_wrfout(os.path.join(d, "cpp"), T)
    bad, n = exact.compare_stream([os.path.join(d, "cpp")], [os.path.join(d, "wrf")])
    print(f"[{'ok' if bad else 'WRONG'}] exact a variable only CPPWRF writes: {str(bad)[:160] if bad else 'identical'}")
    if not bad:
        failures.append("exact cppwrf-only variable")


if __name__ == "__main__":
    sys.exit(main())
