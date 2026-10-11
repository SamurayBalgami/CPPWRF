#!/usr/bin/env python3
"""Compare two simulation runs at their last common snapshot."""

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


def main():
    parser = argparse.ArgumentParser(description="Compare two model runs at last common snapshot")
    parser.add_argument("run_a", help="Path to first run directory")
    parser.add_argument("run_b", help="Path to second run directory")
    parser.add_argument("--var", default="th_pert", help="Variable name to compare (default: th_pert)")
    parser.add_argument("--out", default="out/plots/compare.png", help="Output PNG path")
    args = parser.parse_args()

    snaps_a = set(list_snapshots(args.run_a))
    snaps_b = set(list_snapshots(args.run_b))
    common = sorted(snaps_a & snaps_b)

    if not common:
        print(f"Error: No common snapshots found between {args.run_a} and {args.run_b}", file=sys.stderr)
        sys.exit(2)

    step = common[-1]
    snap_a = read_snapshot(args.run_a, step)
    snap_b = read_snapshot(args.run_b, step)

    if args.var not in snap_a or args.var not in snap_b:
        print(f"Error: Variable '{args.var}' not present in both snapshots", file=sys.stderr)
        sys.exit(2)

    data_a = snap_a[args.var]
    data_b = snap_b[args.var]

    if data_a.shape != data_b.shape:
        print(
            f"Error: Shape mismatch for variable '{args.var}': "
            f"{args.run_a} has shape {data_a.shape} but {args.run_b} has shape {data_b.shape}",
            file=sys.stderr,
        )
        sys.exit(2)

    diff = data_a - data_b
    max_abs_diff = float(np.max(np.abs(diff)))
    l2_diff = float(np.sqrt(np.sum(diff**2)))
    l2_a = float(np.sqrt(np.sum(data_a**2)))
    rel_l2 = float(l2_diff / l2_a) if l2_a > 1e-14 else (0.0 if l2_diff == 0.0 else float("inf"))

    print(f"max|A-B|: {max_abs_diff:.6e}")
    print(f"Relative L2: {rel_l2:.6e}")

    # Process slices for plotting
    run_dict_a = read_run(args.run_a)
    dx = find_key(run_dict_a, "dx")
    dx = float(dx) if dx is not None else 100.0

    z_m = snap_a.get("z_m")
    if z_m is None:
        z_m = snap_b.get("z_m")

    # Handle staggering for u and w if plotting them
    if args.var == "u":
        plot_a = 0.5 * (data_a[:, :, :-1] + data_a[:, :, 1:])
        plot_b = 0.5 * (data_b[:, :, :-1] + data_b[:, :, 1:])
        plot_diff = 0.5 * (diff[:, :, :-1] + diff[:, :, 1:])
    elif args.var == "w":
        plot_a = 0.5 * (data_a[:-1, :, :] + data_a[1:, :, :])
        plot_b = 0.5 * (data_b[:-1, :, :] + data_b[1:, :, :])
        plot_diff = 0.5 * (diff[:-1, :, :] + diff[1:, :, :])
    else:
        plot_a = data_a
        plot_b = data_b
        plot_diff = diff

    nz = plot_a.shape[0]
    ny = plot_a.shape[1] if plot_a.ndim == 3 else 1
    nx = plot_a.shape[-1]
    j = ny // 2

    slice_a = plot_a[:, j, :] if plot_a.ndim == 3 else plot_a
    slice_b = plot_b[:, j, :] if plot_b.ndim == 3 else plot_b
    slice_diff = plot_diff[:, j, :] if plot_diff.ndim == 3 else plot_diff

    if z_m is not None:
        z_slice = z_m[:, j, :] if z_m.ndim == 3 else z_m
    else:
        z_slice = np.linspace(0, 6400, nz)[:, np.newaxis] * np.ones((nz, nx))

    x_km = (np.arange(nx) + 0.5) * dx / 1000.0
    z_km = z_slice / 1000.0
    X_km = np.broadcast_to(x_km[np.newaxis, :], z_km.shape)

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), constrained_layout=True)
    ax_a, ax_b, ax_diff = axes

    # Common scale for A and B
    vmin = min(float(np.nanmin(slice_a)), float(np.nanmin(slice_b)))
    vmax = max(float(np.nanmax(slice_a)), float(np.nanmax(slice_b)))
    if vmin == vmax:
        vmin -= 1.0
        vmax += 1.0
    levels_data = np.linspace(vmin, vmax, 21)

    cf_a = ax_a.contourf(X_km, z_km, slice_a, levels=levels_data, cmap="viridis", extend="both")
    cf_b = ax_b.contourf(X_km, z_km, slice_b, levels=levels_data, cmap="viridis", extend="both")

    # Symmetric scale for difference panel
    diff_abs = max(abs(float(np.nanmin(slice_diff))), abs(float(np.nanmax(slice_diff))))
    if diff_abs < 1e-15:
        diff_abs = 1e-15
    levels_diff = np.linspace(-diff_abs, diff_abs, 21)
    cf_diff = ax_diff.contourf(X_km, z_km, slice_diff, levels=levels_diff, cmap="RdBu_r", extend="both")

    plt.colorbar(cf_a, ax=ax_a)
    plt.colorbar(cf_b, ax=ax_b)
    plt.colorbar(cf_diff, ax=ax_diff)

    name_a = os.path.basename(os.path.normpath(args.run_a))
    name_b = os.path.basename(os.path.normpath(args.run_b))
    ax_a.set_title(f"A: {name_a} ({args.var}) t_step={step}")
    ax_b.set_title(f"B: {name_b} ({args.var}) t_step={step}")
    ax_diff.set_title(f"Difference (A - B)")

    ax_a.set_ylabel("z (km)")
    ax_b.set_ylabel("z (km)")
    ax_diff.set_ylabel("z (km)")
    ax_diff.set_xlabel("x (km)")

    out_path = args.out
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
