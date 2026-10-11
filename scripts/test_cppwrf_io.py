#!/usr/bin/env python3
"""Unit tests for cppwrf_io Python reader and cross-language verification."""

import argparse
import json
import os
import subprocess
import sys
import tempfile
import numpy as np

# Ensure scripts directory is in sys.path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cppwrf_io import list_snapshots, read_diag, read_run, read_snapshot


def test_python_roundtrip():
    print("Testing Python-side snapshot generation and read_snapshot...")
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Create run.json
        run_data = {"case": "straka", "grid": {"nx": 12, "ny": 8, "nz": 5}}
        with open(os.path.join(tmp_dir, "run.json"), "w", encoding="utf-8") as f:
            json.dump(run_data, f)

        # Create diag.csv
        with open(os.path.join(tmp_dir, "diag.csv"), "w", encoding="utf-8") as f:
            f.write("step,time,ke,pe\n")
            f.write("0,0.0,10.5,20.5\n")
            f.write("1,0.1,10.4,20.6\n")

        # Generate test arrays
        nz, ny, nx = 5, 8, 12
        rng = np.random.default_rng(42)
        u = rng.standard_normal((nz, ny, nx)).astype(np.float32)
        temp = rng.standard_normal((nz, ny, nx)).astype(np.float64)
        psfc = rng.standard_normal((1, ny, nx)).astype(np.float32)

        # Pack binary back to back
        u_bytes = u.tobytes()
        temp_bytes = temp.tobytes()
        psfc_bytes = psfc.tobytes()

        bin_data = u_bytes + temp_bytes + psfc_bytes

        meta = {
            "step": 42,
            "time": 10.5,
            "variables": [
                {
                    "name": "u",
                    "stagger": "U",
                    "shape": [nz, ny, nx],
                    "dtype": "float32",
                    "offset": 0,
                    "nbytes": len(u_bytes),
                },
                {
                    "name": "temp",
                    "stagger": "M",
                    "shape": [nz, ny, nx],
                    "dtype": "float64",
                    "offset": len(u_bytes),
                    "nbytes": len(temp_bytes),
                },
                {
                    "name": "psfc",
                    "stagger": "M",
                    "shape": [1, ny, nx],
                    "dtype": "float32",
                    "offset": len(u_bytes) + len(temp_bytes),
                    "nbytes": len(psfc_bytes),
                },
            ],
        }

        with open(os.path.join(tmp_dir, "snap_000042.bin"), "wb") as f:
            f.write(bin_data)
        with open(os.path.join(tmp_dir, "snap_000042.json"), "w", encoding="utf-8") as f:
            json.dump(meta, f)

        # Test list_snapshots
        snaps = list_snapshots(tmp_dir)
        assert snaps == [42], f"Expected [42], got {snaps}"

        # Test read_run
        read_run_data = read_run(tmp_dir)
        assert read_run_data == run_data, f"Run data mismatch: {read_run_data}"

        # Test read_diag
        diag = read_diag(tmp_dir)
        assert "step" in diag and "time" in diag and "ke" in diag and "pe" in diag
        np.testing.assert_allclose(diag["step"], np.array([0.0, 1.0]))
        np.testing.assert_allclose(diag["time"], np.array([0.0, 0.1]))
        np.testing.assert_allclose(diag["ke"], np.array([10.5, 10.4]))
        np.testing.assert_allclose(diag["pe"], np.array([20.5, 20.6]))

        # Test read_snapshot
        snap_vars = read_snapshot(tmp_dir, 42)
        assert "u" in snap_vars
        assert "temp" in snap_vars
        assert "psfc" in snap_vars

        assert snap_vars["u"].dtype == np.float32
        assert snap_vars["temp"].dtype == np.float64
        assert snap_vars["psfc"].dtype == np.float32

        np.testing.assert_array_equal(snap_vars["u"], u)
        np.testing.assert_array_equal(snap_vars["temp"], temp)
        np.testing.assert_array_equal(snap_vars["psfc"], psfc)

    print("Python roundtrip tests PASSED.")


def test_cross_language(cross_writer_bin: str):
    print(f"Testing C++ writer -> Python reader cross-language roundtrip using {cross_writer_bin}...")
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Run C++ writer executable
        res = subprocess.run([cross_writer_bin, tmp_dir], capture_output=True, text=True)
        if res.returncode != 0:
            print(f"C++ writer failed:\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")
            sys.exit(1)

        # Verify run.json
        run_data = read_run(tmp_dir)
        assert run_data["case"] == "straka"
        assert run_data["nx"] == 4
        assert run_data["ny"] == 3
        assert run_data["nz"] == 2

        # Verify diag.csv
        diag = read_diag(tmp_dir)
        assert "step" in diag and "time" in diag and "energy" in diag and "temp" in diag
        np.testing.assert_allclose(diag["step"], [0.0, 100.0, 200.0])
        np.testing.assert_allclose(diag["time"], [0.0, 12.5, 25.0])
        np.testing.assert_allclose(diag["energy"], [1000.5, 999.8, 998.2])
        np.testing.assert_allclose(diag["temp"], [300.0, 300.2, 300.4])

        # Verify snapshot listing
        snaps = list_snapshots(tmp_dir)
        assert snaps == [100, 200], f"Expected [100, 200], got {snaps}"

        # Verify each snapshot data
        nx, ny, nz = 4, 3, 2
        for step, time_val in [(100, 12.5), (200, 25.0)]:
            vars_dict = read_snapshot(tmp_dir, step)
            assert "u" in vars_dict
            assert "temp" in vars_dict
            assert "psfc" in vars_dict

            u = vars_dict["u"]
            temp = vars_dict["temp"]
            psfc = vars_dict["psfc"]

            assert u.shape == (nz, ny, nx)
            assert temp.shape == (nz, ny, nx)
            assert psfc.shape == (1, ny, nx)

            assert u.dtype == np.float32
            assert temp.dtype == np.float64
            assert psfc.dtype == np.float32

            # Analytical formulas from test_cross_writer.cpp:
            # u[k, j, i] = step + i*1.0f + j*10.0f + k*100.0f
            # temp[k, j, i] = step*2.0 + i*0.1 + j*0.01 + k*0.001
            # psfc[0, j, i] = step*10.0f + i + j*5.0f
            expected_u = np.zeros((nz, ny, nx), dtype=np.float32)
            expected_temp = np.zeros((nz, ny, nx), dtype=np.float64)
            expected_psfc = np.zeros((1, ny, nx), dtype=np.float32)

            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        expected_u[k, j, i] = float(step + i * 1.0 + j * 10.0 + k * 100.0)
                        expected_temp[k, j, i] = float(step * 2.0 + i * 0.1 + j * 0.01 + k * 0.001)
            for j in range(ny):
                for i in range(nx):
                    expected_psfc[0, j, i] = float(step * 10.0 + i + j * 5.0)

            np.testing.assert_array_equal(u, expected_u)
            np.testing.assert_array_equal(temp, expected_temp)
            np.testing.assert_array_equal(psfc, expected_psfc)

    print("Cross-language tests PASSED.")


def main():
    parser = argparse.ArgumentParser(description="Test cppwrf_io reader")
    parser.add_argument("--cross-writer", type=str, default=None, help="Path to test_io_cross_writer executable")
    args = parser.parse_args()

    test_python_roundtrip()

    cross_writer = args.cross_writer
    if not cross_writer:
        # Search candidate locations
        candidates = [
            "/tmp/cppwrf_io_build/test_io_cross_writer",
            "/tmp/cppwrf_io_build/src_io/test_io_cross_writer",
            os.path.join(os.path.dirname(__file__), "../build/test_io_cross_writer"),
        ]
        for c in candidates:
            if os.path.isfile(c) and os.access(c, os.X_OK):
                cross_writer = c
                break

    if cross_writer and os.path.isfile(cross_writer):
        test_cross_language(cross_writer)
    else:
        print("Note: cross-writer binary not provided or found; skipping cross-language test.")

    print("ALL PYTHON IO TESTS PASSED.")


if __name__ == "__main__":
    main()
