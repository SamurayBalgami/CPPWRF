#!/usr/bin/env python3
"""Plot Straka convergence comparison and front position vs resolution."""

import argparse
import os
import sys
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cppwrf_io import list_snapshots, read_run, read_snapshot


def find_key(d, key):
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


def compute_straka_front(th_pert, dx):
    nx = th_pert.shape[-1]
    x = (np.arange(nx) + 0.5) * dx
    xc = (nx * dx) / 2.0
    j = th_pert.shape[1] // 2 if th_pert.ndim == 3 else 0
    row = th_pert[0, j, :] if th_pert.ndim == 3 else th_pert[0, :]  # lowest model level k = 0

    # Outermost x (from the centre, right side) where th_pert <= -1 K, linearly interpolated
    idx = np.where((row <= -1.0) & (x >= xc))[0]
    if len(idx) > 0:
        i = idx.max()
        j_next = i + 1
        if 0 <= j_next < nx and row[j_next] > -1.0:
            xf = x[i] + (x[j_next] - x[i]) * (-1.0 - row[i]) / (row[j_next] - row[i])
        else:
            xf = x[i]
        front_km = abs(xf - xc) / 1000.0
    else:
        front_km = float("nan")

    min_th_pert = float(np.min(th_pert))
    return front_km, min_th_pert


def main():
    parser = argparse.ArgumentParser(description="Straka convergence comparison across resolutions")
    parser.add_argument("run_dirs", nargs="+", help="Paths to Straka run directories")
    parser.add_argument("--out", default="out/plots/straka_convergence.png", help="Output PNG path")
    args = parser.parse_args()

    results = []
    run_data = []

    for rdir in args.run_dirs:
        rdict = read_run(rdir)
        dx_val = find_key(rdict, "dx")
        if dx_val is None:
            print(f"Warning: dx not found in {rdir}/run.json, skipping", file=sys.stderr)
            continue
        dx = float(dx_val)

        snaps = list_snapshots(rdir)
        if not snaps:
            print(f"Warning: No snapshots found in {rdir}, skipping", file=sys.stderr)
            continue
        last_step = snaps[-1]

        snap_vars = read_snapshot(rdir, last_step)
        th_pert = snap_vars["th_pert"]
        z_m = snap_vars.get("z_m")

        front_km, min_th = compute_straka_front(th_pert, dx)
        results.append((dx, front_km, min_th))
        run_data.append((rdir, dx, th_pert, z_m))

    # Print table to stdout
    print(f"{'dx':<8} {'front_km':<10} {'min_th_pert':<12}")
    for dx, f_km, m_th in results:
        print(f"{dx:<8.0f} {f_km:<10.3f} {m_th:<12.3f}")

    # Plot
    N = len(run_data)
    if N == 0:
        print("No valid run directories to plot", file=sys.stderr)
        sys.exit(1)

    fig, axes = plt.subplots(1, N + 1, figsize=(4 * N + 4, 3.8), constrained_layout=True)
    if N + 1 == 1:
        axes = [axes]

    levels = np.arange(-16.5, -0.4, 1.0)
    cf_last = None

    for idx, (rdir, dx, th_pert, z_m) in enumerate(run_data):
        ax = axes[idx]
        nz = th_pert.shape[0]
        ny = th_pert.shape[1] if th_pert.ndim == 3 else 1
        nx = th_pert.shape[-1]
        j = ny // 2

        th_slice = th_pert[:, j, :] if th_pert.ndim == 3 else th_pert

        if z_m is not None:
            z_slice = z_m[:, j, :] if z_m.ndim == 3 else z_m
        else:
            z_slice = np.linspace(0, 6400, nz)[:, np.newaxis] * np.ones((nz, nx))

        x_km = (np.arange(nx) + 0.5) * dx / 1000.0
        z_km = z_slice / 1000.0
        X_km = np.broadcast_to(x_km[np.newaxis, :], z_km.shape)

        cf = ax.contourf(X_km, z_km, th_slice, levels=levels, cmap="coolwarm", extend="both")
        ax.contour(X_km, z_km, th_slice, levels=[-1.0], colors="black", linewidths=1.0)
        cf_last = cf

        centre_km = (nx * dx / 2.0) / 1000.0
        max_x_km = (nx * dx) / 1000.0
        ax.set_xlim(centre_km, max_x_km)
        ax.set_ylim(0, 4.6)
        ax.set_title(f"dx = {dx:g} m")
        ax.set_xlabel("x (km)")
        if idx == 0:
            ax.set_ylabel("z (km)")

    # Colorbar for contour plots
    if cf_last is not None:
        plt.colorbar(cf_last, ax=axes[:N], label="θ' (K)", orientation="vertical", shrink=0.8)

    # Front position vs dx panel
    ax_front = axes[N]
    dxs = [r[0] for r in results]
    fronts = [r[1] for r in results]

    ax_front.plot(dxs, fronts, "o-", color="navy", lw=2, markersize=6)
    ax_front.set_title("Front position vs dx")
    ax_front.set_xlabel("dx (m)")
    ax_front.set_ylabel("front (km)")
    ax_front.grid(True, linestyle="--", alpha=0.6)

    out_path = args.out
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
