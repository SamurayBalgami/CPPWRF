#!/usr/bin/env python3
"""Compare WRF reference run against CPPWRF for the Straka density current."""

import argparse
import os
import sys
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cppwrf_io import list_snapshots, read_diag, read_run, read_snapshot


def compute_front_and_min(row_lowest, x_coords, centre_x):
    """Compute outermost x where row_lowest <= -1.0 K (right side of centre), linearly interpolated."""
    nx = len(x_coords)
    idx = np.where((row_lowest <= -1.0) & (x_coords >= centre_x))[0]
    if len(idx) > 0:
        i = idx.max()
        j = i + 1
        if 0 <= j < nx and row_lowest[j] > -1.0:
            xf = x_coords[i] + (x_coords[j] - x_coords[i]) * (-1.0 - row_lowest[i]) / (row_lowest[j] - row_lowest[i])
        else:
            xf = x_coords[i]
        front_km = abs(xf - centre_x) / 1000.0
    else:
        front_km = float("nan")
    return front_km


def main():
    parser = argparse.ArgumentParser(description="Compare WRF and CPPWRF Straka runs")
    parser.add_argument("wrfout", help="Path to wrfout NetCDF file")
    parser.add_argument("cppwrf_dir", help="Path to CPPWRF run directory")
    parser.add_argument("--wall-wrf", type=float, default=None, help="WRF wall time in seconds")
    parser.add_argument("--out", default="out/plots/wrf_vs_cppwrf_straka.png", help="Output PNG path")
    args = parser.parse_args()

    # 1. Process WRF output
    ds = nc.Dataset(args.wrfout)
    t_wrf = ds.variables["T"][:]  # perturbation potential temperature (theta - 300)
    dx_wrf = float(ds.DX)
    nt_wrf, nz_wrf, ny_wrf, nx_wrf = t_wrf.shape

    last_t_idx = -1
    th_pert_wrf = t_wrf[last_t_idx]  # shape (nz, ny, nx)

    x_wrf = (np.arange(nx_wrf) + 0.5) * dx_wrf
    centre_wrf = (nx_wrf * dx_wrf) / 2.0  # 25600 m (domain centre in em_grav2d_x)

    # Height at mass levels: 0.5 * ((PH + PHB)[:-1] + (PH + PHB)[1:]) / 9.81
    ph = ds.variables["PH"][last_t_idx]
    phb = ds.variables["PHB"][last_t_idx]
    gz = ph + phb
    z_stag_wrf = gz / 9.81
    z_m_wrf = 0.5 * (z_stag_wrf[:-1, 0, :] + z_stag_wrf[1:, 0, :])

    row_lowest_wrf = th_pert_wrf[0, 0, :]
    min_th_wrf = float(np.min(th_pert_wrf))
    front_km_wrf = compute_front_and_min(row_lowest_wrf, x_wrf, centre_wrf)

    # 2. Process CPPWRF output
    cppwrf_run = read_run(args.cppwrf_dir)
    cppwrf_diag = read_diag(args.cppwrf_dir)
    cppwrf_steps = list_snapshots(args.cppwrf_dir)
    cppwrf_last_step = cppwrf_steps[-1]
    cppwrf_snap = read_snapshot(args.cppwrf_dir, cppwrf_last_step)

    th_pert_cpp = cppwrf_snap["th_pert"]  # [nz, ny, nx]
    z_m_cpp = cppwrf_snap.get("z_m")
    dx_cpp = float(cppwrf_run.get("grid", {}).get("dx", 100.0))

    nz_cpp, ny_cpp, nx_cpp = th_pert_cpp.shape
    j_cpp = ny_cpp // 2
    th_slice_cpp = th_pert_cpp[:, j_cpp, :]
    z_slice_cpp = z_m_cpp[:, j_cpp, :] if z_m_cpp is not None else np.linspace(50, 6350, nz_cpp)[:, None] * np.ones((nz_cpp, nx_cpp))

    x_cpp = (np.arange(nx_cpp) + 0.5) * dx_cpp
    centre_cpp = (nx_cpp * dx_cpp) / 2.0

    row_lowest_cpp = th_slice_cpp[0, :]
    min_th_cpp = float(np.min(th_pert_cpp))
    front_km_cpp = compute_front_and_min(row_lowest_cpp, x_cpp, centre_cpp)

    # Wall time for CPPWRF from diag.csv
    wall_cpp = float(cppwrf_diag["wall_s"][-1]) if "wall_s" in cppwrf_diag and len(cppwrf_diag["wall_s"]) > 0 else float("nan")
    wall_wrf_str = f"{args.wall_wrf:.2f}" if args.wall_wrf is not None else "N/A"

    # Print table
    print("\n" + "=" * 65)
    print(f"{'model':<20} | {'dx (m)':<8} | {'front_km':<10} | {'min_th_pert':<12} | {'wall_s':<8}")
    print("-" * 65)
    print(f"{'WRF v4.6.1 (serial)':<20} | {dx_wrf:<8.0f} | {front_km_wrf:<10.3f} | {min_th_wrf:<12.3f} | {wall_wrf_str:<8}")
    print(f"{'CPPWRF':<20} | {dx_cpp:<8.0f} | {front_km_cpp:<10.3f} | {min_th_cpp:<12.3f} | {wall_cpp:<8.2f}")
    print("=" * 65 + "\n")

    # 3. Plotting: WRF on top, CPPWRF on bottom
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(10, 6.5), constrained_layout=True)
    levels = np.arange(-16.5, -0.4, 1.0)

    # WRF panel (crop right half: x >= centre_wrf, distance from centre)
    right_mask_wrf = (x_wrf >= centre_wrf)
    x_rel_km_wrf = (x_wrf[right_mask_wrf] - centre_wrf) / 1000.0
    th_crop_wrf = th_pert_wrf[:, 0, right_mask_wrf]
    z_crop_wrf = z_m_wrf[:, right_mask_wrf] / 1000.0
    X_wrf = np.broadcast_to(x_rel_km_wrf[None, :], z_crop_wrf.shape)

    cf1 = ax_top.contourf(X_wrf, z_crop_wrf, th_crop_wrf, levels=levels, cmap="coolwarm", extend="both")
    ax_top.contour(X_wrf, z_crop_wrf, th_crop_wrf, levels=[-1.0], colors="black", linewidths=1.0)
    ax_top.set_xlim(0, 18.0)
    ax_top.set_ylim(0, 4.6)
    ax_top.set_title(f"WRF v4.6.1 (em_grav2d_x) t=900s, dx={dx_wrf:g}m (front={front_km_wrf:.3f} km, min θ'={min_th_wrf:.2f} K)")
    ax_top.set_ylabel("z (km)")

    # CPPWRF panel (crop right half: x >= centre_cpp, distance from centre)
    right_mask_cpp = (x_cpp >= centre_cpp)
    x_rel_km_cpp = (x_cpp[right_mask_cpp] - centre_cpp) / 1000.0
    th_crop_cpp = th_slice_cpp[:, right_mask_cpp]
    z_crop_cpp = z_slice_cpp[:, right_mask_cpp] / 1000.0
    X_cpp = np.broadcast_to(x_rel_km_cpp[None, :], z_crop_cpp.shape)

    cf2 = ax_bot.contourf(X_cpp, z_crop_cpp, th_crop_cpp, levels=levels, cmap="coolwarm", extend="both")
    ax_bot.contour(X_cpp, z_crop_cpp, th_crop_cpp, levels=[-1.0], colors="black", linewidths=1.0)
    ax_bot.set_xlim(0, 18.0)
    ax_bot.set_ylim(0, 4.6)
    cpp_name = os.path.basename(os.path.normpath(args.cppwrf_dir))
    ax_bot.set_title(f"CPPWRF ({cpp_name}) t=900s, dx={dx_cpp:g}m (front={front_km_cpp:.3f} km, min θ'={min_th_cpp:.2f} K)")
    ax_bot.set_ylabel("z (km)")
    ax_bot.set_xlabel("x distance from bubble centre (km)")

    plt.colorbar(cf2, ax=[ax_top, ax_bot], label="θ' (K)", orientation="vertical", shrink=0.85)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    fig.savefig(args.out, dpi=150)
    plt.close(fig)
    print(f"Saved figure: {args.out}")


if __name__ == "__main__":
    main()
