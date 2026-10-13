#!/usr/bin/env python3
"""Compare WRF reference run against CPPWRF for the 2D hill flow (em_hill2d_x).

Validates:
1. Initial state (t = 0) against WRF wrfinput_d01 / ref_fields.npz for:
   MUB, PHB, PH, T, U, W, P, PB, AL
2. Simulation evolution at t = 600, 1200, 1800, 2400, 3000, 3600 s for:
   W, U, T' (max abs diff and RMS diff).
"""

import argparse
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cppwrf_io import list_snapshots, read_run, read_snapshot

try:
    import netCDF4 as nc
    HAVE_NETCDF = True
except ImportError:
    HAVE_NETCDF = False


def compare_t0(cppwrf_dir: str, ref_npz_path: str, wrfinput_path: str = None):
    print("=" * 90)
    print(f"VERIFICATION AT t = 0 s ({cppwrf_dir})")
    print("=" * 90)

    snap0 = read_snapshot(cppwrf_dir, 0)
    ref = np.load(ref_npz_path)

    # Reference fields
    u_wrf = ref["U"][0]
    w_wrf = ref["W"][0]
    t_wrf = ref["T"][0]  # theta - 300 in WRF
    ph_wrf = ref["PH"][0]
    mub_wrf = ref["MUB"]
    phb_wrf = ref["PHB"]

    # CPPWRF fields
    u_cpp = snap0["u"][:, 0, :]
    w_cpp = snap0["w"][:, 0, :]
    theta_cpp = snap0["theta"][:, 0, :]
    t_cpp = theta_cpp - 300.0
    p_full_cpp = snap0["p"][:, 0, :]
    p_pert_cpp = snap0["p_pert"][:, 0, :]
    ph_full_cpp = snap0["ph"][:, 0, :]
    mu_full_cpp = snap0["mu"][0, 0, :]

    # In CPPWRF at t=0:
    # mu_pert = 0, p_pert = 0, ph_pert = 0, al_pert = 0
    # so base fields equal full fields:
    mub_cpp = mu_full_cpp
    phb_cpp = ph_full_cpp
    pb_cpp = p_full_cpp
    p_wrf = np.zeros_like(p_full_cpp)

    results = []

    # 1. MUB
    d_mub = np.abs(mub_cpp - mub_wrf)
    results.append(("MUB", d_mub.max(), (d_mub / mub_wrf).max(), "Column dry air mass [Pa]"))

    # 2. PHB
    d_phb = np.abs(phb_cpp - phb_wrf)
    results.append(("PHB", d_phb.max(), (d_phb / np.maximum(phb_wrf, 1.0)).max(), "Base geopotential [m2 s-2]"))

    # 3. PH (perturbation geopotential)
    # In CPPWRF ph_pert is exactly 0; in WRF wrfout it has single-precision roundoff noise
    d_ph = np.abs(0.0 - ph_wrf)
    results.append(("PH", d_ph.max(), d_ph.max() / np.maximum(phb_wrf.max(), 1.0), "Perturbation geopotential [m2 s-2]"))

    # 4. T (perturbation potential temperature theta - 300)
    d_t = np.abs(t_cpp - t_wrf)
    results.append(("T", d_t.max(), (d_t / np.abs(theta_cpp)).max(), "Potential temperature pert theta-300 [K]"))

    # 5. U
    d_u = np.abs(u_cpp - u_wrf)
    results.append(("U", d_u.max(), (d_u / u_wrf).max(), "Zonal velocity [m s-1]"))

    # 6. W
    # wrfinput has W=0, but wrfout frame 0 has diagnostic surface boundary w = u * dh/dx
    d_w_wrfout = np.abs(w_cpp - w_wrf)
    results.append(("W (vs wrfout)", d_w_wrfout.max(), float("nan"), "Vertical velocity [m s-1] (wrfout has diag w_sfc)"))
    results.append(("W (vs wrfinput)", np.abs(w_cpp).max(), 0.0, "Vertical velocity [m s-1] (wrfinput W=0)"))

    # 7. P (perturbation pressure)
    # P is 0 in both wrfinput and CPPWRF
    results.append(("P", np.abs(p_pert_cpp).max(), 0.0, "Perturbation pressure [Pa]"))

    # 8. PB & AL (from wrfinput if available)
    if wrfinput_path and os.path.exists(wrfinput_path) and HAVE_NETCDF:
        ds_in = nc.Dataset(wrfinput_path)
        pb_wrf = ds_in.variables["PB"][0, :, 0, :].data
        al_wrf = ds_in.variables["AL"][0, :, 0, :].data
        alb_wrf = ds_in.variables["ALB"][0, :, 0, :].data

        d_pb = np.abs(pb_cpp - pb_wrf)
        results.append(("PB", d_pb.max(), (d_pb / pb_wrf).max(), "Base pressure [Pa]"))

        # In CPPWRF, AL is perturbation inverse density = 0.0
        d_al = np.abs(0.0 - al_wrf)
        results.append(("AL", d_al.max(), (d_al / alb_wrf).max(), "Perturbation inverse density [m3 kg-1]"))

        # ALB
        r_d = 287.0
        p1000mb = 100000.0
        cvpm = -717.5 / 1004.5
        alb_cpp = (r_d / p1000mb) * theta_cpp * ((pb_cpp / p1000mb) ** cvpm)
        d_alb = np.abs(alb_cpp - alb_wrf)
        results.append(("ALB", d_alb.max(), (d_alb / alb_wrf).max(), "Base inverse density [m3 kg-1]"))

    print(f"{'Field':<16} | {'Max Abs Diff':<16} | {'Max Rel Diff':<16} | {'Description'}")
    print("-" * 90)
    for name, mabs, mrel, desc in results:
        rel_str = f"{mrel:<16.4e}" if not np.isnan(mrel) else f"{'N/A':<16}"
        print(f"{name:<16} | {mabs:<16.4e} | {rel_str} | {desc}")
    print("=" * 90 + "\n")


def compare_simulation(cppwrf_dir: str, ref_npz_path: str):
    print("=" * 90)
    print(f"SIMULATION COMPARISON (t = 0..3600 s) ({cppwrf_dir} vs WRF)")
    print("=" * 90)

    ref = np.load(ref_npz_path)
    ref_times = ref["times"]
    n_frames = len(ref_times)

    run_meta = read_run(cppwrf_dir)
    dt = float(run_meta.get("grid", {}).get("dt", 20.0))
    time_step_sec = 600.0  # 10 min interval

    print(f"{'Time (s)':<10} | {'W Max Diff':<12} | {'W RMS Diff':<12} | {'U Max Diff':<12} | {'U RMS Diff':<12} | {'T Max Diff':<12} | {'T RMS Diff':<12}")
    print("-" * 90)

    stats = []
    for f in range(n_frames):
        sim_sec = int(round(f * time_step_sec))
        step = int(round(sim_sec / dt))

        snap = read_snapshot(cppwrf_dir, step)
        w_cpp = snap["w"][:, 0, :]
        u_cpp = snap["u"][:, 0, :]
        t_cpp = snap["theta"][:, 0, :] - 300.0

        w_wrf = ref["W"][f]
        u_wrf = ref["U"][f]
        t_wrf = ref["T"][f]

        dw = np.abs(w_cpp - w_wrf)
        du = np.abs(u_cpp - u_wrf)
        dt_val = np.abs(t_cpp - t_wrf)

        w_max, w_rms = float(dw.max()), float(np.sqrt(np.mean(dw**2)))
        u_max, u_rms = float(du.max()), float(np.sqrt(np.mean(du**2)))
        t_max, t_rms = float(dt_val.max()), float(np.sqrt(np.mean(dt_val**2)))

        stats.append({
            "time": sim_sec, "step": step,
            "w_max": w_max, "w_rms": w_rms,
            "u_max": u_max, "u_rms": u_rms,
            "t_max": t_max, "t_rms": t_rms
        })

        print(f"{sim_sec:<10} | {w_max:<12.4e} | {w_rms:<12.4e} | {u_max:<12.4e} | {u_rms:<12.4e} | {t_max:<12.4e} | {t_rms:<12.4e}")

    print("=" * 90 + "\n")
    return stats


def main():
    parser = argparse.ArgumentParser(description="Compare CPPWRF hill2d run against WRF em_hill2d_x reference")
    parser.add_argument("--cppwrf-dir", default="out/hill2d_p0_dp", help="CPPWRF run output directory")
    parser.add_argument("--ref", default="/home/ostriquetrum/wrf_ref/run_hill_p0/ref_fields.npz", help="Path to ref_fields.npz")
    parser.add_argument("--wrfinput", default="/home/ostriquetrum/wrf_ref/run_hill_p0/wrfinput_d01", help="Path to wrfinput_d01")
    args = parser.parse_args()

    compare_t0(args.cppwrf_dir, args.ref, args.wrfinput)
    compare_simulation(args.cppwrf_dir, args.ref)


if __name__ == "__main__":
    main()
