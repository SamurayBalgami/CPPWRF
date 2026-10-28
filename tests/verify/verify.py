#!/usr/bin/env python3
"""CPPWRF independent verification.

Black box: reads only the files the model writes (the output format in docs/USER_GUIDE.md):
run.json, diag.csv, snap_<step>.bin/.json. Deliberately does not import scripts/
or any CPPWRF code, so a bug in those cannot hide a bug in the model.

Usage:
  verify.py check <case> <out_dir>          # case: rest_state | straka | warm_bubble | igw | squall2d
  verify.py compare <out_dir_a> <out_dir_b> [--rtol X]
  verify.py convergence <dir_coarse> <dir_mid> <dir_fine>   # Straka at dx, dx/2, dx/4
  verify.py wrf <wrfout> <out_dir> [--rms-tol K]   # theta' vs a real WRF run on the same grid
Exit code 0 = all PASS (SKIP allowed), 1 = any FAIL, 2 = usage/input error.
"""
import csv
import glob
import json
import math
import os
import sys

import numpy as np

# ---------------------------------------------------------------------------
# Tolerances. Each one has a reason; change them only with a reason.
# ---------------------------------------------------------------------------
TOL = {
    # Flat terrain, base state in discrete hydrostatic balance: the perturbation
    # equations have an exact zero solution, so any motion is roundoff or a bug.
    "rest_max_w": {"float64": 1e-8, "float32": 1e-4},          # m/s over the run
    "rest_max_thp": {"float64": 1e-8, "float32": 1e-4},        # K
    # Flux form + closed/periodic boundaries conserve dry mass and theta-mass
    # to roundoff (TN ch. 3). Relative drift over the run:
    "mass_drift": {"float64": 1e-10, "float32": 1e-5},
    "theta_mass_drift": {"float64": 1e-9, "float32": 1e-5},
    # Mirror symmetry of symmetric cases (relative to max |field|). Upwind-biased
    # advection is symmetric in exact arithmetic; roundoff breaks it slowly in FP32.
    "symmetry": {"float64": 1e-6, "float32": 1e-2},
    # Straka et al. (1993), converged reference: min theta' = -9.77 K, front 15.54 km.
    # At dx = 100 m published models fall in about -9.6..-9.9 K and 15.3..15.6 km.
    # min theta' is INFO only: WRF's diff_opt=1 applies no vertical diffusion to the lowest level
    # (vflux(0) = vflux(1)), unlike Straka's spec, which warms the cold pool to about -8.6 K.
    # Fidelity to WRF is checked field by field with `verify.py wrf` instead.
    "straka_ref_min_thp": -9.77,     # K
    "straka_ref_front": 15540.0, "straka_tol_front": 400.0,     # m
    # Moist cases. Water: transport alone is conservative to roundoff; WRF's Kessler is not
    # exactly (its rain sedimentation uses the mass-level spacing), measured 2.7e-4 of the total
    # water after 1 h of deep convection (0.75% of the rain), so allow 1e-3.
    "water_budget_mp0": {"float64": 1e-12, "float32": 1e-5},
    "water_budget_mp": 1e-3,
    "wsm6_sym_time": 600.0,  # s, WSM6 squall2d: last snapshot checked for FP64 mirror symmetry
    # Straka vs WRF v4.8.0 (FP32 output, diff_opt=1, top_lid=.false., same eta levels, p_top,
    # bubble; configs/straka_wrf.ini): measured theta' RMS difference 5.6e-6..7.4e-6 K over 900 s
    # (field RMS ~1 K), i.e. FP32 roundoff of WRF. Tolerance 1e-4 K (14x margin).
    "wrf_rms_thp": 1e-4,     # K
    "wrf_grid_z": 1e-2,      # m, mass-level heights must match (WRF writes FP32: ~2e-3 m at 6 km)
}


# ---------------------------------------------------------------------------
# Readers (independent of scripts/cppwrf_io.py on purpose)
# ---------------------------------------------------------------------------
def read_run(d):
    with open(os.path.join(d, "run.json")) as f:
        return json.load(f)


def read_diag(d):
    path = os.path.join(d, "diag.csv")
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"{path}: no rows")
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def list_steps(d):
    steps = []
    for p in glob.glob(os.path.join(d, "snap_*.json")):
        name = os.path.basename(p)[5:-5]
        if name.isdigit():
            steps.append(int(name))
    return sorted(steps)


def read_snap(d, step):
    base = os.path.join(d, f"snap_{step:06d}")
    with open(base + ".json") as f:
        meta = json.load(f)
    out = {"_time": float(meta["time"]), "_step": int(meta["step"]), "_meta": meta}
    with open(base + ".bin", "rb") as f:
        for v in meta["variables"]:
            dt = np.dtype({"float32": "<f4", "float64": "<f8"}[v["dtype"]])
            shape = tuple(v["shape"])
            f.seek(int(v["offset"]))
            a = np.frombuffer(f.read(int(np.prod(shape)) * dt.itemsize), dtype=dt)
            if a.size != int(np.prod(shape)):
                raise ValueError(f"{base}.bin: truncated variable {v['name']}")
            out[v["name"]] = a.reshape(shape).astype(np.float64)
    return out


def find_key(obj, key):
    """First value for `key` anywhere in a nested dict/list."""
    if isinstance(obj, dict):
        if key in obj:
            return obj[key]
        for v in obj.values():
            r = find_key(v, key)
            if r is not None:
                return r
    elif isinstance(obj, list):
        for v in obj:
            r = find_key(v, key)
            if r is not None:
                return r
    return None


def precision(run):
    p = str(find_key(run, "precision") or "").lower()
    if p in ("float32", "fp32", "single", "float"):
        return "float32"
    if p in ("float64", "fp64", "double"):
        return "float64"
    raise ValueError(f"run.json: unknown precision {p!r}")


# ---------------------------------------------------------------------------
# Result bookkeeping
# ---------------------------------------------------------------------------
class Report:
    def __init__(self, title):
        self.title, self.rows = title, []

    def add(self, name, status, detail):
        self.rows.append((name, status, detail))

    def check(self, name, ok, detail):
        self.add(name, "PASS" if ok else "FAIL", detail)

    def print(self):
        print(f"== {self.title}")
        for n, s, d in self.rows:
            print(f"  {s:4}  {n:24} {d}")

    @property
    def failed(self):
        return any(s == "FAIL" for _, s, _ in self.rows)


# ---------------------------------------------------------------------------
# Generic checks
# ---------------------------------------------------------------------------
def check_finite(rep, d, steps):
    diag = read_diag(d)
    bad = [k for k, v in diag.items() if not np.all(np.isfinite(v))]
    rep.check("finite_diag", not bad, f"non-finite columns: {bad}" if bad else "ok")
    worst = None
    for s in steps:
        snap = read_snap(d, s)
        for k, v in snap.items():
            if not k.startswith("_") and not np.all(np.isfinite(v)):
                worst = (s, k)
                break
        if worst:
            break
    rep.check("finite_snapshots", worst is None,
              f"first non-finite: step {worst[0]} var {worst[1]}" if worst else f"{len(steps)} snapshots ok")


def check_conservation(rep, d, prec):
    diag = read_diag(d)
    # WRF's constant-K vertical diffusion (vertical_diffusion_3dmp) applies g^2/(mu*alt_k)
    # outside the eta derivative, so it does not conserve theta-mass exactly. Faithful
    # to WRF, so with kvdif > 0 the theta-mass drift is reported, not judged.
    kv = float(find_key(read_run(d), "kvdif") or 0)
    for col, tol_key in (("dry_mass", "mass_drift"), ("theta_mass", "theta_mass_drift")):
        if col not in diag:
            rep.add(col, "SKIP", "column missing in diag.csv")
            continue
        v = diag[col]
        drift = abs(v[-1] - v[0]) / abs(v[0]) if v[0] != 0 else math.inf
        mp = int(find_key(read_run(d), "mp_physics") or 0)
        if col == "theta_mass" and mp > 0:
            rep.add(f"{col}_drift", "INFO", f"rel drift {drift:.3e} (mp_physics={mp}: latent heating changes theta-mass)")
            continue
        if col == "theta_mass" and kv > 0:
            rep.add(f"{col}_drift", "INFO", f"rel drift {drift:.3e} (kvdif={kv:g}: WRF form, not conservative)")
            continue
        rep.check(f"{col}_drift", drift <= TOL[tol_key][prec],
                  f"rel drift {drift:.3e} (tol {TOL[tol_key][prec]:.0e})")


def mirror_error(a, axis_len_stagger_u=False):
    """Relative asymmetry about the domain centre along x (last axis)."""
    m = a[..., ::-1]
    if axis_len_stagger_u:      # u at faces: mirror maps u -> -u
        err = np.max(np.abs(a + m))
    else:
        err = np.max(np.abs(a - m))
    scale = np.max(np.abs(a))
    return err / scale if scale > 0 else 0.0


def check_symmetry(rep, d, steps, prec):
    first = read_snap(d, steps[0])
    init_err = mirror_error(first["th_pert"])
    if init_err > 1e-6:
        rep.add("symmetry", "SKIP", f"initial th_pert not centred (asym {init_err:.2e})")
        return
    last = read_snap(d, steps[-1])
    tol = TOL["symmetry"][prec]
    e_th = mirror_error(last["th_pert"])
    e_u = mirror_error(last["u"], axis_len_stagger_u=True)
    e_w = mirror_error(last["w"])
    worst = max(e_th, e_u, e_w)
    rep.check("symmetry", worst <= tol,
              f"th' {e_th:.2e}  u {e_u:.2e}  w {e_w:.2e}  (tol {tol:.0e}, t={last['_time']:g}s)")


# ---------------------------------------------------------------------------
# Case-specific checks
# ---------------------------------------------------------------------------
def straka_metrics(d):
    """min theta' and front distance from the cold-bubble centre at the final snapshot."""
    run = read_run(d)
    steps = list_steps(d)
    first, last = read_snap(d, steps[0]), read_snap(d, steps[-1])
    dx = float(find_key(run, "dx"))
    th0, th = first["th_pert"], last["th_pert"]
    nx = th.shape[-1]
    x = (np.arange(nx) + 0.5) * dx
    xc = find_key(run, "xc")
    if xc is None:              # centroid of the initial cold anomaly
        w = np.clip(-th0, 0, None).sum(axis=(0, 1))
        xc = float((w * x).sum() / w.sum())
    xc = float(xc)
    row = th[0, 0, :] if th.ndim == 3 else th[0, :]   # lowest mass level, j = 0
    fronts = []
    for side in (+1, -1):
        idx = np.where((row <= -1.0) & (side * (x - xc) >= 0))[0]
        if idx.size == 0:
            continue
        i = idx.max() if side > 0 else idx.min()
        j = i + side
        if 0 <= j < nx and row[j] > -1.0:    # linear interpolation to -1 K crossing
            xf = x[i] + (x[j] - x[i]) * (-1.0 - row[i]) / (row[j] - row[i])
        else:
            xf = x[i]
        fronts.append(abs(xf - xc))
    return {"min_thp": float(th.min()), "front": max(fronts) if fronts else float("nan"),
            "fronts": fronts, "time": last["_time"], "dx": dx}


def run_check(case, d):
    rep = Report(f"{case}: {d}")
    try:
        run = read_run(d)
        prec = precision(run)
        steps = list_steps(d)
        if not steps:
            raise ValueError("no snapshots")
    except Exception as e:  # noqa: BLE001 - report any input problem as FAIL
        rep.add("inputs", "FAIL", str(e))
        rep.print()
        return rep
    rep.add("info", "INFO", f"precision={prec} backend={find_key(run, 'backend')} "
                            f"snapshots={len(steps)} t_end={read_snap(d, steps[-1])['_time']:g}s")
    check_finite(rep, d, steps)
    if rep.failed:
        rep.print()
        return rep
    check_conservation(rep, d, prec)

    if case == "rest_state":
        diag = read_diag(d)
        mw = float(np.max(diag["max_abs_w"]))
        mt = float(max(np.max(np.abs(diag["min_th_pert"])), np.max(np.abs(diag["max_th_pert"]))))
        rep.check("rest_max_w", mw <= TOL["rest_max_w"][prec], f"{mw:.3e} m/s (tol {TOL['rest_max_w'][prec]:.0e})")
        rep.check("rest_max_thp", mt <= TOL["rest_max_thp"][prec], f"{mt:.3e} K (tol {TOL['rest_max_thp'][prec]:.0e})")
    elif case == "straka":
        check_symmetry(rep, d, steps, prec)
        m = straka_metrics(d)
        if abs(m["time"] - 900.0) > 1.0:
            rep.add("straka_time", "FAIL", f"last snapshot at t={m['time']:g}s, need 900s")
        dmin = m["min_thp"] - TOL["straka_ref_min_thp"]
        rep.add("straka_min_thp", "INFO",
                f"{m['min_thp']:.3f} K (Straka ref -9.77, diff {dmin:+.3f}; WRF lowest-level BC, see TOL)")
        dfr = m["front"] - TOL["straka_ref_front"]
        rep.check("straka_front", abs(dfr) <= TOL["straka_tol_front"],
                  f"{m['front'] / 1000:.3f} km (ref 15.54, diff {dfr / 1000:+.3f}); sides "
                  + ", ".join(f"{f / 1000:.3f}" for f in m["fronts"]))
    elif case == "warm_bubble":
        check_symmetry(rep, d, steps, prec)
        diag = read_diag(d)
        rep.add("warm_bubble_max_w", "INFO", f"{np.max(diag['max_abs_w']):.4f} m/s (cross-backend reference)")
    elif case == "igw":
        # Skamarock & Klemp (1994): the perturbation is advected by the mean wind U
        # and spreads symmetrically about x0 + U t (linear, Boussinesq limit).
        last = read_snap(d, steps[-1])
        first = read_snap(d, steps[0])
        U = find_key(run, "u_mean") or find_key(run, "U")
        dx = float(find_key(run, "dx"))
        th0, th = first["th_pert"], last["th_pert"]
        x = (np.arange(th.shape[-1]) + 0.5) * dx
        w0 = np.abs(th0).sum(axis=tuple(range(th0.ndim - 1)))
        w1 = np.abs(th).sum(axis=tuple(range(th.ndim - 1)))
        x0, x1 = (w0 * x).sum() / w0.sum(), (w1 * x).sum() / w1.sum()
        if U is None:
            rep.add("igw_translation", "SKIP", "mean wind (u_mean) not in run.json")
        else:
            expect = x0 + float(U) * (last["_time"] - first["_time"])
            Lx = th.shape[-1] * dx
            err = (x1 - expect + Lx / 2) % Lx - Lx / 2            # periodic domain
            rep.check("igw_translation", abs(err) <= 2 * dx,
                      f"centroid {x1 / 1000:.2f} km, expected {expect % Lx / 1000:.2f} km (err {err:.0f} m)")
        rep.add("igw_amplitude", "INFO", f"th' range [{th.min():.5f}, {th.max():.5f}] K at t={last['_time']:g}s")
    elif case == "squall2d":
        # Deep convection amplifies roundoff (chaos): mirror symmetry is a hard check in FP64
        # only; in FP32 it is reported.
        mp = int(find_key(run, "mp_physics") or 0)
        if prec == "float64" and mp == 6:
            # WSM6 has threshold branches (nucleation, freezing, qmin cut-offs) that turn the
            # dycore's ~1e-12 mirror roundoff into 1e-5 after the first ice forms (measured:
            # 1.3e-12 at 600 s, 3e-5 at 1200 s; column kernel is deterministic). Hard check
            # before the ice phase, INFO afterwards.
            early = [s for s in steps if read_snap(d, s)["_time"] <= TOL["wsm6_sym_time"]]
            check_symmetry(rep, d, early, prec)
            last = read_snap(d, steps[-1])
            rep.add("symmetry_late", "INFO", f"th' asym {mirror_error(last['th_pert']):.2e} at t={last['_time']:g}s")
        elif prec == "float64":
            check_symmetry(rep, d, steps, prec)
        else:
            last = read_snap(d, steps[-1])
            rep.add("symmetry", "INFO", f"FP32 th' asym {mirror_error(last['th_pert']):.2e} (chaotic growth)")
        mp = int(find_key(run, "mp_physics") or 0)
        eta = np.array(find_key(run, "eta_levels"), dtype=float)
        dnw = -np.diff(eta)
        first, last = read_snap(d, steps[0]), read_snap(d, steps[-1])
        species = [k for k in ("qv", "qc", "qr", "qi", "qs", "qg") if k in last]
        neg = min(float(last[k].min()) for k in species)
        rep.check("moist_nonnegative", neg >= 0.0, f"min over {species}: {neg:.3e}")
        def water(s):  # kg per m^2 of domain, summed over columns (dry-air mass weighting)
            mu = s["mu"][0]
            return sum((mu[None, :, :] * dnw[:, None, None] * s[k]).sum() for k in species) / 9.81
        W0, W1 = water(first), water(last)
        rain = float(last["rainnc"].sum()) if "rainnc" in last else 0.0
        err = (W1 + rain - W0) / W0
        tol = TOL["water_budget_mp"] if mp > 0 else TOL["water_budget_mp0"][prec]
        rep.check("water_budget", abs(err) <= tol,
                  f"(water + rain - W0)/W0 = {err:+.2e} (tol {tol:.0e}, mp_physics={mp}, rain {rain:.3e} kg/m2 summed)")
        diag = read_diag(d)
        rep.add("squall_max_w", "INFO", f"{np.max(diag['max_abs_w']):.2f} m/s")
    else:
        rep.add("case", "FAIL", f"unknown case {case}")
    rep.print()
    return rep


def run_compare(a, b, rtol):
    """Same case, different backend or precision: final-snapshot relative L2 differences."""
    rep = Report(f"compare {a} vs {b}")
    sa, sb = list_steps(a), list_steps(b)
    common = sorted(set(sa) & set(sb))
    if not common:
        rep.add("steps", "FAIL", "no common snapshot steps")
        rep.print()
        return rep
    A, B = read_snap(a, common[-1]), read_snap(b, common[-1])
    for k in sorted(set(A) ^ set(B)):  # a variable one run writes and the other does not
        rep.add(k, "FAIL", f"only in {a if k in A else b}")
    for k in sorted(A):
        if k.startswith("_") or k not in B:
            continue
        if A[k].shape != B[k].shape:
            rep.add(k, "FAIL", f"shape {A[k].shape} vs {B[k].shape}")
            continue
        scale = np.linalg.norm(A[k])
        if k in ("p", "theta", "ph", "mu", "z_m"):   # large means: compare deviations
            ref = A[k].mean()
            scale = np.linalg.norm(A[k] - ref)
        rel = np.linalg.norm(A[k] - B[k]) / scale if scale > 0 else float(np.linalg.norm(A[k] - B[k]))
        rep.check(k, rel <= rtol, f"rel L2 diff {rel:.3e} (tol {rtol:.0e}) step {common[-1]}")
    rep.print()
    return rep


def run_convergence(dirs):
    rep = Report("straka convergence")
    ms = [straka_metrics(d) for d in dirs]
    for m in ms:
        rep.add(f"dx={m['dx']:g}", "INFO", f"front {m['front'] / 1000:.3f} km, min th' {m['min_thp']:.3f} K")
    f = [m["front"] for m in ms]
    d1, d2 = abs(f[0] - f[1]), abs(f[1] - f[2])
    if d2 == 0:
        rep.add("order", "INFO", "identical fine-grid fronts; order undefined")
    else:
        ratio = ms[0]["dx"] / ms[1]["dx"]
        p = math.log(d1 / d2) / math.log(ratio) if d1 > 0 else float("nan")
        rep.check("front_converges", d2 < d1, f"|df| {d1:.0f} m -> {d2:.0f} m, observed order {p:.2f}")
    rep.print()
    return rep


def run_wrf(wrfout, d, rms_tol):
    """theta' of a CPPWRF x-z run against a WRF wrfout file at every common output time."""
    import netCDF4
    rep = Report(f"wrf {os.path.basename(wrfout)} vs {d}")
    ds = netCDF4.Dataset(wrfout)
    t_wrf = np.array(ds.variables["T"][:], dtype=np.float64)[:, :, 0, :]   # (time, z, x), j = 0
    xt = np.array(ds.variables["XTIME"][:], dtype=np.float64) * 60.0       # minutes -> s
    gz = np.array(ds.variables["PH"][0] + ds.variables["PHB"][0], dtype=np.float64)[:, 0, :] / 9.81
    z_wrf = 0.5 * (gz[:-1] + gz[1:])
    steps = list_steps(d)
    snaps = {read_snap(d, s)["_time"]: s for s in steps}
    first = read_snap(d, steps[0])
    zc = first["z_m"][:, 0, :]
    if zc.shape != z_wrf.shape:
        rep.add("grid", "FAIL", f"shape {zc.shape} vs WRF {z_wrf.shape}")
        rep.print()
        return rep
    dz = float(np.abs(zc - z_wrf).max())
    rep.check("grid_z", dz <= TOL["wrf_grid_z"], f"max |z_cpp - z_wrf| {dz:.2e} m")
    n = 0
    for it, tw in enumerate(xt):
        if tw not in snaps:
            continue
        thc = read_snap(d, snaps[tw])["th_pert"][:, 0, :]
        e = thc - t_wrf[it]
        rms, ref = float(np.sqrt(np.mean(e * e))), float(np.sqrt(np.mean(t_wrf[it] ** 2)))
        rep.check(f"thp_t{tw:g}s", rms <= rms_tol,
                  f"RMS diff {rms:.2e} K, max {np.abs(e).max():.2e} K (field RMS {ref:.3f} K, tol {rms_tol})")
        n += 1
    if n == 0:
        rep.add("times", "FAIL", "no common output times")
    rep.print()
    return rep


def main(argv):
    if len(argv) >= 3 and argv[0] == "check":
        return 1 if run_check(argv[1], argv[2]).failed else 0
    if len(argv) >= 3 and argv[0] == "compare":
        rtol = float(argv[argv.index("--rtol") + 1]) if "--rtol" in argv else 1e-3
        return 1 if run_compare(argv[1], argv[2], rtol).failed else 0
    if len(argv) == 4 and argv[0] == "convergence":
        return 1 if run_convergence(argv[1:]).failed else 0
    if len(argv) >= 3 and argv[0] == "wrf":
        tol = float(argv[argv.index("--rms-tol") + 1]) if "--rms-tol" in argv else TOL["wrf_rms_thp"]
        return 1 if run_wrf(argv[1], argv[2], tol).failed else 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
