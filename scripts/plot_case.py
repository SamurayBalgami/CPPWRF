#!/usr/bin/env python3
"""Plot x-z cross sections of th_pert, w, and u from CPPWRF model snapshots."""

import argparse
import json
import os
import sys
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cppwrf_io import list_snapshots, read_run, read_snapshot


def find_key(d, key):
    """Recursively search for key in nested dictionary."""
    if not isinstance(d, dict):
        return None
    if key in d:
        return d[key]
    for v in d.values():
        if isinstance(v, dict):
            res = find_key(v, key)
            if res is not None:
                return res
    return None


def get_case_and_dx(run_dir):
    run_dict = read_run(run_dir)
    case = run_dict.get("case", "unknown")
    dx = find_key(run_dict, "dx")
    if dx is None:
        raise ValueError(f"Could not find dx in {run_dir}/run.json")
    return case, float(dx)


def plot_single_snapshot(run_dir, step, out_path, case, dx):
    snap_vars = read_snapshot(run_dir, step)

    # Read time from json if available
    json_path = os.path.join(run_dir, f"snap_{step:06d}.json")
    time_val = float(step)
    if os.path.exists(json_path):
        try:
            with open(json_path, "r", encoding="utf-8") as f:
                time_val = float(json.load(f).get("time", step))
        except Exception:
            pass

    th_pert = snap_vars["th_pert"]  # [nz, ny, nx]
    u = snap_vars["u"]              # [nz, ny, nx+1]
    w = snap_vars["w"]              # [nz+1, ny, nx]
    z_m = snap_vars.get("z_m")       # [nz, ny, nx]

    nz, ny, nx = th_pert.shape
    j = ny // 2

    # Destagger u to mass points in x: 0.5 * (u[..., :-1] + u[..., 1:])
    u_mass = 0.5 * (u[:, :, :-1] + u[:, :, 1:])
    # Destagger w to mass levels in z: 0.5 * (w[:-1, ...] + w[1:, ...])
    w_mass = 0.5 * (w[:-1, :, :] + w[1:, :, :])

    th_slice = th_pert[:, j, :]
    u_slice = u_mass[:, j, :]
    w_slice = w_mass[:, j, :]

    if z_m is not None:
        z_slice = z_m[:, j, :]
    else:
        # Fallback if z_m not available
        z_slice = np.linspace(0, 6400, nz)[:, np.newaxis] * np.ones((nz, nx))

    x_km = (np.arange(nx) + 0.5) * dx / 1000.0
    z_km = z_slice / 1000.0
    X_km = np.broadcast_to(x_km[np.newaxis, :], z_km.shape)

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), constrained_layout=True)
    ax_th, ax_w, ax_u = axes

    # Case-specific levels for th_pert
    if case == "straka":
        levels = np.arange(-16.5, -0.4, 1.0)
        cf_th = ax_th.contourf(X_km, z_km, th_slice, levels=levels, cmap="coolwarm", extend="both")
        ax_th.contour(X_km, z_km, th_slice, levels=[-1.0], colors="black", linewidths=1.0)
        centre_km = (nx * dx / 2.0) / 1000.0
        max_x_km = (nx * dx) / 1000.0
        for ax in axes:
            ax.set_xlim(centre_km, max_x_km)
            ax.set_ylim(0, 4.6)
    elif case == "igw":
        levels = [val for val in np.arange(-0.0015, 0.0031, 0.0005) if abs(val) > 1e-7]
        cf_th = ax_th.contourf(X_km, z_km, th_slice, levels=levels, cmap="RdBu_r", extend="both")
    elif case == "warm_bubble":
        levels = np.arange(0.25, 2.26, 0.25)
        cf_th = ax_th.contourf(X_km, z_km, th_slice, levels=levels, cmap="plasma", extend="both")
    else:
        vmin, vmax = float(np.nanmin(th_slice)), float(np.nanmax(th_slice))
        if vmin < 0 < vmax:
            vmax_abs = max(abs(vmin), abs(vmax))
            levels = np.linspace(-vmax_abs, vmax_abs, 15)
            cf_th = ax_th.contourf(X_km, z_km, th_slice, levels=levels, cmap="coolwarm", extend="both")
        else:
            cf_th = ax_th.contourf(X_km, z_km, th_slice, 15, cmap="viridis", extend="both")

    # w levels (symmetric about 0)
    w_abs = max(abs(float(np.nanmin(w_slice))), abs(float(np.nanmax(w_slice))))
    w_levels = np.linspace(-w_abs, w_abs, 15) if w_abs > 1e-7 else np.linspace(-1, 1, 15)
    cf_w = ax_w.contourf(X_km, z_km, w_slice, levels=w_levels, cmap="RdBu_r", extend="both")

    # u levels (symmetric about 0)
    u_abs = max(abs(float(np.nanmin(u_slice))), abs(float(np.nanmax(u_slice))))
    u_levels = np.linspace(-u_abs, u_abs, 15) if u_abs > 1e-7 else np.linspace(-1, 1, 15)
    cf_u = ax_u.contourf(X_km, z_km, u_slice, levels=u_levels, cmap="RdBu_r", extend="both")

    plt.colorbar(cf_th, ax=ax_th, label="K")
    plt.colorbar(cf_w, ax=ax_w, label="m/s")
    plt.colorbar(cf_u, ax=ax_u, label="m/s")

    time_str = f"{time_val:g}"
    ax_th.set_title(f"{case} t={time_str} s th_pert")
    ax_w.set_title(f"{case} t={time_str} s w")
    ax_u.set_title(f"{case} t={time_str} s u")

    ax_th.set_ylabel("z (km)")
    ax_w.set_ylabel("z (km)")
    ax_u.set_ylabel("z (km)")
    ax_u.set_xlabel("x (km)")

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Plot model cross-sections for a simulation case")
    parser.add_argument("run_dir", help="Path to simulation output directory")
    parser.add_argument("--out", default=None, help="Output directory or file path")
    parser.add_argument("--steps", choices=["all", "last"], default="all", help="Snapshots to plot")
    args = parser.parse_args()

    case, dx = get_case_and_dx(args.run_dir)
    all_steps = list_snapshots(args.run_dir)
    if not all_steps:
        print(f"No snapshots found in {args.run_dir}", file=sys.stderr)
        sys.exit(1)

    steps = [all_steps[-1]] if args.steps == "last" else all_steps

    out_arg = args.out
    for step in steps:
        if out_arg and out_arg.endswith(".png"):
            if len(steps) == 1:
                out_path = out_arg
            else:
                base, ext = os.path.splitext(out_arg)
                out_path = f"{base}_snap_{step:06d}{ext}"
        else:
            out_dir = out_arg if out_arg else f"out/plots/{case}"
            out_path = os.path.join(out_dir, f"snap_{step:06d}.png")

        plot_single_snapshot(args.run_dir, step, out_path, case, dx)


if __name__ == "__main__":
    main()
