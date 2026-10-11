#!/usr/bin/env python3
"""Plot diagnostic time series from diag.csv across one or more simulation runs."""

import argparse
import os
import sys
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cppwrf_io import read_diag


def main():
    parser = argparse.ArgumentParser(description="Plot diagnostics time series from diag.csv")
    parser.add_argument("run_dirs", nargs="+", help="Paths to simulation run directories")
    parser.add_argument("--out", default="out/plots/diag.png", help="Output PNG path")
    args = parser.parse_args()

    fig, axes = plt.subplots(5, 1, figsize=(10, 13), sharex=True, constrained_layout=True)
    ax_w, ax_tmin, ax_tmax, ax_dm, ax_tm = axes

    for run_dir in args.run_dirs:
        diag = read_diag(run_dir)
        if not diag:
            print(f"Warning: No diagnostics found in {run_dir}, skipping", file=sys.stderr)
            continue

        label = os.path.basename(os.path.normpath(run_dir))
        t = diag["time"]

        ax_w.plot(t, diag["max_abs_w"], label=label, lw=1.5)
        ax_tmin.plot(t, diag["min_th_pert"], label=label, lw=1.5)
        ax_tmax.plot(t, diag["max_th_pert"], label=label, lw=1.5)

        # Relative drift: (x - x0)/x0
        dm = diag["dry_mass"]
        if len(dm) > 0 and dm[0] != 0:
            dm_drift = np.abs((dm - dm[0]) / dm[0])
            dm_drift = np.maximum(dm_drift, 1e-16)  # Guard zeros for log y-axis
            ax_dm.plot(t, dm_drift, label=label, lw=1.5)

        tm = diag["theta_mass"]
        if len(tm) > 0 and tm[0] != 0:
            tm_drift = np.abs((tm - tm[0]) / tm[0])
            tm_drift = np.maximum(tm_drift, 1e-16)  # Guard zeros for log y-axis
            ax_tm.plot(t, tm_drift, label=label, lw=1.5)

    ax_w.set_ylabel("max |w| (m/s)")
    ax_w.set_title("Maximum Vertical Velocity")
    ax_w.grid(True, linestyle="--", alpha=0.6)
    ax_w.legend(loc="best")

    ax_tmin.set_ylabel("min θ' (K)")
    ax_tmin.set_title("Minimum Potential Temperature Perturbation")
    ax_tmin.grid(True, linestyle="--", alpha=0.6)

    ax_tmax.set_ylabel("max θ' (K)")
    ax_tmax.set_title("Maximum Potential Temperature Perturbation")
    ax_tmax.grid(True, linestyle="--", alpha=0.6)

    ax_dm.set_ylabel("|Δ dry_mass / dry_mass_0|")
    ax_dm.set_yscale("log")
    ax_dm.set_title("Dry Mass Relative Drift")
    ax_dm.grid(True, which="both", linestyle="--", alpha=0.6)

    ax_tm.set_ylabel("|Δ θ_mass / θ_mass_0|")
    ax_tm.set_yscale("log")
    ax_tm.set_title("Theta Mass Relative Drift")
    ax_tm.grid(True, which="both", linestyle="--", alpha=0.6)
    ax_tm.set_xlabel("time (s)")

    out_path = args.out
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
