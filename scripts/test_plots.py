#!/usr/bin/env python3
"""Automated tests for plotting scripts in scripts/."""

import json
import os
import subprocess
import sys
import tempfile
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_straka_convergence import compute_straka_front


def write_synthetic_run(target_dir, case_name, nx, ny, nz, dx, dt=1.0, front_target=None):
    os.makedirs(target_dir, exist_ok=True)

    # 1. run.json
    run_dict = {
        "case": case_name,
        "grid": {
            "nx": nx,
            "ny": ny,
            "nz": nz,
            "dx": dx,
            "dy": dx,
            "dt": dt,
            "ztop": 6400.0,
        },
    }
    with open(os.path.join(target_dir, "run.json"), "w", encoding="utf-8") as f:
        json.dump(run_dict, f)

    # 2. diag.csv
    diag_content = (
        "step,time,max_abs_u,max_abs_v,max_abs_w,min_th_pert,max_th_pert,dry_mass,theta_mass,total_energy,wall_s\n"
        "0,0,0,0,0,-15,0,1000000,2000000,3000000,0.01\n"
        "10,10,2.5,0,1.8,-12,0.2,1000001,2000002,3000001,0.05\n"
    )
    with open(os.path.join(target_dir, "diag.csv"), "w", encoding="utf-8") as f:
        f.write(diag_content)

    # 3. Snapshots (step 0 and step 10)
    xc = (nx * dx) / 2.0
    x = (np.arange(nx) + 0.5) * dx

    for step in [0, 10]:
        u = np.zeros((nz, ny, nx + 1), dtype=np.float32)
        w = np.zeros((nz + 1, ny, nx), dtype=np.float32)
        th = np.zeros((nz, ny, nx), dtype=np.float32)
        z_levels = np.linspace(50.0, 6350.0, nz, dtype=np.float32)
        z_m = np.broadcast_to(z_levels[:, None, None], (nz, ny, nx)).copy()

        if case_name == "straka":
            target = front_target if front_target is not None else (nx * dx / 4.0)
            # Linear crossing at xc +- target: th = -1.0 - 0.005 * (target - abs(x - xc))
            dist = np.abs(x - xc)
            # slope = 0.002 K / m
            row_th = -1.0 - 0.002 * (target - dist)
            th[0, :, :] = np.broadcast_to(row_th[None, :], (ny, nx))
            th[1:, :, :] = -5.0 * np.exp(-((dist[None, None, :] / target) ** 2))
        elif case_name == "warm_bubble":
            th[:] = 1.0
        elif case_name == "igw":
            th[:] = 0.001
        else:
            th[:] = 0.5

        # Pack into .bin
        u_b = u.tobytes()
        w_b = w.tobytes()
        th_b = th.tobytes()
        z_b = z_m.tobytes()
        bin_data = u_b + w_b + th_b + z_b

        meta = {
            "step": step,
            "time": float(step * dt),
            "variables": [
                {"name": "u", "stagger": "U", "shape": list(u.shape), "dtype": "float32", "offset": 0, "nbytes": len(u_b)},
                {"name": "w", "stagger": "W", "shape": list(w.shape), "dtype": "float32", "offset": len(u_b), "nbytes": len(w_b)},
                {"name": "th_pert", "stagger": "M", "shape": list(th.shape), "dtype": "float32", "offset": len(u_b) + len(w_b), "nbytes": len(th_b)},
                {"name": "z_m", "stagger": "M", "shape": list(z_m.shape), "dtype": "float32", "offset": len(u_b) + len(w_b) + len(th_b), "nbytes": len(z_b)},
            ],
        }

        with open(os.path.join(target_dir, f"snap_{step:06d}.bin"), "wb") as f:
            f.write(bin_data)
        with open(os.path.join(target_dir, f"snap_{step:06d}.json"), "w", encoding="utf-8") as f:
            json.dump(meta, f)


def check_png_valid(path):
    assert os.path.exists(path), f"File {path} does not exist"
    assert os.path.getsize(path) > 1000, f"File {path} is too small ({os.path.getsize(path)} bytes)"


def main():
    print("Running test_plots.py...")
    scripts_dir = os.path.dirname(os.path.abspath(__file__))

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Create test runs
        d_straka = os.path.join(tmp_dir, "run_straka")
        d_warm = os.path.join(tmp_dir, "run_warm")
        d_igw = os.path.join(tmp_dir, "run_igw")
        d_3d = os.path.join(tmp_dir, "run_3d")
        d_mismatch = os.path.join(tmp_dir, "run_mismatch")

        known_front = 1250.0  # meters from centre
        write_synthetic_run(d_straka, "straka", nx=64, ny=1, nz=16, dx=100.0, front_target=known_front)
        write_synthetic_run(d_warm, "warm_bubble", nx=64, ny=1, nz=16, dx=100.0)
        write_synthetic_run(d_igw, "igw", nx=64, ny=1, nz=16, dx=1000.0)
        write_synthetic_run(d_3d, "case_3d", nx=32, ny=4, nz=16, dx=100.0)
        write_synthetic_run(d_mismatch, "straka", nx=40, ny=1, nz=16, dx=100.0)

        # 1. Test plot_case.py
        print("Testing plot_case.py...")
        out_case_dir = os.path.join(tmp_dir, "plots_case")
        res = subprocess.run([
            sys.executable,
            os.path.join(scripts_dir, "plot_case.py"),
            d_straka,
            "--out", out_case_dir,
            "--steps", "all"
        ], capture_output=True, text=True)
        assert res.returncode == 0, f"plot_case.py failed:\n{res.stderr}"
        check_png_valid(os.path.join(out_case_dir, "snap_000000.png"))
        check_png_valid(os.path.join(out_case_dir, "snap_000010.png"))

        # Test single output file with --steps last
        single_png = os.path.join(tmp_dir, "single_case.png")
        res2 = subprocess.run([
            sys.executable,
            os.path.join(scripts_dir, "plot_case.py"),
            d_3d,
            "--out", single_png,
            "--steps", "last"
        ], capture_output=True, text=True)
        assert res2.returncode == 0, f"plot_case.py (3D) failed:\n{res2.stderr}"
        check_png_valid(single_png)

        # Test warm_bubble and igw cases
        for d in [d_warm, d_igw]:
            res = subprocess.run([
                sys.executable,
                os.path.join(scripts_dir, "plot_case.py"),
                d,
                "--out", os.path.join(tmp_dir, f"out_{os.path.basename(d)}.png"),
                "--steps", "last"
            ], capture_output=True, text=True)
            assert res.returncode == 0, f"plot_case.py failed on {d}:\n{res.stderr}"

        # 2. Test plot_diag.py
        print("Testing plot_diag.py...")
        out_diag = os.path.join(tmp_dir, "diag.png")
        res_diag = subprocess.run([
            sys.executable,
            os.path.join(scripts_dir, "plot_diag.py"),
            d_straka, d_warm, d_3d,
            "--out", out_diag
        ], capture_output=True, text=True)
        assert res_diag.returncode == 0, f"plot_diag.py failed:\n{res_diag.stderr}"
        check_png_valid(out_diag)

        # 3. Test plot_compare.py
        print("Testing plot_compare.py...")
        out_comp = os.path.join(tmp_dir, "comp.png")
        res_comp = subprocess.run([
            sys.executable,
            os.path.join(scripts_dir, "plot_compare.py"),
            d_straka, d_straka,
            "--out", out_comp
        ], capture_output=True, text=True)
        assert res_comp.returncode == 0, f"plot_compare.py failed:\n{res_comp.stderr}"
        assert "max|A-B|: 0.000000e+00" in res_comp.stdout, f"Expected 0 diff, got {res_comp.stdout}"
        check_png_valid(out_comp)

        # Test plot_compare failure on shape mismatch (must exit 2)
        res_comp_fail = subprocess.run([
            sys.executable,
            os.path.join(scripts_dir, "plot_compare.py"),
            d_straka, d_mismatch,
            "--out", os.path.join(tmp_dir, "comp_fail.png")
        ], capture_output=True, text=True)
        assert res_comp_fail.returncode == 2, f"Expected exit code 2, got {res_comp_fail.returncode}"

        # 4. Test plot_straka_convergence.py
        print("Testing plot_straka_convergence.py...")
        out_conv = os.path.join(tmp_dir, "conv.png")
        res_conv = subprocess.run([
            sys.executable,
            os.path.join(scripts_dir, "plot_straka_convergence.py"),
            d_straka,
            "--out", out_conv
        ], capture_output=True, text=True)
        assert res_conv.returncode == 0, f"plot_straka_convergence.py failed:\n{res_conv.stderr}"
        check_png_valid(out_conv)

        # Verify table in stdout
        assert "dx" in res_conv.stdout and "front_km" in res_conv.stdout and "min_th_pert" in res_conv.stdout

        # Verify front accuracy on synthetic field
        from cppwrf_io import read_snapshot as r_snap
        snap10 = r_snap(d_straka, 10)
        f_km, _ = compute_straka_front(snap10["th_pert"], dx=100.0)
        computed_front_m = f_km * 1000.0
        err_m = abs(computed_front_m - known_front)
        assert err_m <= 0.5 * 100.0, f"Front error {err_m} m exceeds half grid cell (50 m)"
        print(f"Known front: {known_front} m, computed: {computed_front_m} m, error: {err_m} m (PASSED)")

    print("ALL PLOT TESTS PASSED.")


if __name__ == "__main__":
    main()
