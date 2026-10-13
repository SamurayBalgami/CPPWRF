#!/usr/bin/env python3
"""Bitwise restart acceptance test suite for Task T6.2.
Tests straka (600 s straight vs 300 s + restart + 300 s)
and squall2d (600 s straight vs 300 s + restart + 300 s).
Validates:
1. cmp of final snap_*.bin (100% bitwise byte-for-byte match).
2. columns 1-10 of diag.csv are character-for-character identical.
Runs for:
- FP32 CPU
- FP64 CPU
- FP32 CUDA
- FP64 CUDA
"""

import os
import subprocess
import sys

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT_ROOT = os.path.join(BASE_DIR, "tests", "io_netcdf", "test_out", "acceptance")

def run_cmd(cmd):
    full_cmd = [os.path.join(BASE_DIR, "tools", "heavy.sh"), "--exclusive"] + cmd
    print(f">> Executing: {' '.join(full_cmd)}")
    res = subprocess.run(full_cmd, cwd=BASE_DIR, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if res.returncode != 0:
        print(f"Error output:\n{res.stdout}")
        raise RuntimeError(f"Command failed with exit code {res.returncode}")
    return res.stdout

def compare_diag(f_cont, f_rst, start_step):
    with open(f_cont) as f:
        c_lines = [l.strip().split(",") for l in f if l.strip()]
    with open(f_rst) as f:
        r_lines = [l.strip().split(",") for l in f if l.strip()]

    # Header check
    assert c_lines[0][:10] == r_lines[0][:10], f"Header mismatch: {c_lines[0][:10]} vs {r_lines[0][:10]}"

    c_sub = [row[:10] for row in c_lines[1:] if float(row[0]) >= start_step]
    r_sub = [row[:10] for row in r_lines[1:] if float(row[0]) >= start_step]

    assert len(c_sub) == len(r_sub), f"Row count mismatch for steps >= {start_step}: {len(c_sub)} vs {len(r_sub)}"
    for i, (cr, rr) in enumerate(zip(c_sub, r_sub)):
        assert cr == rr, f"Mismatch at step {cr[0]}:\n  cont: {cr}\n  rst:  {rr}"
    print(f"   PASS: diag.csv columns 1-10 identical across {len(c_sub)} rows (steps {start_step} to {c_sub[-1][0]})")

def compare_snap_bin(f_cont, f_rst):
    res = subprocess.run(["cmp", f_cont, f_rst], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert res.returncode == 0, f"cmp failed between {f_cont} and {f_rst}: {res.stderr.decode()}"
    size = os.path.getsize(f_cont)
    print(f"   PASS: cmp {f_cont} and {f_rst} bitwise identical ({size} bytes, exit code 0)")

def run_test(case_name, config_file, bdir, prec, dt, t_split, t_final):
    tag = f"{case_name}_{bdir}_{prec}"
    print(f"\n==================================================")
    print(f"Testing {tag}")
    print(f"==================================================")

    out_cont = os.path.join(OUT_ROOT, f"{tag}_cont")
    out_p1 = os.path.join(OUT_ROOT, f"{tag}_p1")
    out_p2 = os.path.join(OUT_ROOT, f"{tag}_p2")

    for d in (out_cont, out_p1, out_p2):
        subprocess.run(["rm", "-rf", d])
        os.makedirs(d, exist_ok=True)

    bin_path = os.path.join(BASE_DIR, bdir, "cppwrf")
    step_split = int(round(t_split / dt))
    step_final = int(round(t_final / dt))

    # 1. Continuous run to t_final
    print(f"1. Continuous run to {t_final} s (steps 0 to {step_final})...")
    run_cmd([bin_path, config_file,
             f"--set=run.run_time={t_final}",
             f"--set=output.dir={out_cont}",
             f"--set=output.out_interval={step_final}",
             f"--set=output.diag_interval=10"])

    # 2. Split segment 1 to t_split with restart output at t_split
    print(f"2. Split segment 1 to {t_split} s (steps 0 to {step_split}, restart.interval={t_split})...")
    run_cmd([bin_path, config_file,
             f"--set=run.run_time={t_split}",
             f"--set=output.dir={out_p1}",
             f"--set=output.out_interval={step_split}",
             f"--set=output.diag_interval=10",
             f"--set=restart.interval={t_split}"])

    rst_file = os.path.join(out_p1, f"restart_{step_split:06d}.nc")
    assert os.path.exists(rst_file), f"Restart file {rst_file} not found!"

    # 3. Split segment 2 restarting from t_split to t_final
    print(f"3. Split segment 2 restarting from {t_split} s to {t_final} s (steps {step_split} to {step_final})...")
    run_cmd([bin_path, config_file,
             f"--set=run.run_time={t_final}",
             f"--set=run.restart_file={rst_file}",
             f"--set=output.dir={out_p2}",
             f"--set=output.out_interval={step_final}",
             f"--set=output.diag_interval=10"])

    # 4. Compare final snapshot .bin
    snap_cont = os.path.join(out_cont, f"snap_{step_final:06d}.bin")
    snap_p2 = os.path.join(out_p2, f"snap_{step_final:06d}.bin")
    compare_snap_bin(snap_cont, snap_p2)

    # 5. Compare diag.csv columns 1-10
    diag_cont = os.path.join(out_cont, "diag.csv")
    diag_p2 = os.path.join(out_p2, "diag.csv")
    compare_diag(diag_cont, diag_p2, step_split)

    print(f"VERIFICATION SUCCESSFUL for {tag}!")

def main():
    os.makedirs(OUT_ROOT, exist_ok=True)

    # Test cases configurations:
    # straka: dt = 1.0, split = 300 s, final = 600 s
    # squall2d: dt = 6.0, split = 300 s, final = 600 s

    configs = [
        ("straka", "configs/straka.ini", "build-cpu", "fp32", 1.0, 300.0, 600.0),
        ("straka", "configs/straka.ini", "build-cpu-dp", "fp64", 1.0, 300.0, 600.0),
        ("straka", "configs/straka.ini", "build-cuda", "fp32", 1.0, 300.0, 600.0),
        ("straka", "configs/straka.ini", "build-cuda-dp", "fp64", 1.0, 300.0, 600.0),
        ("squall2d", "configs/squall2d.ini", "build-cpu", "fp32", 6.0, 300.0, 600.0),
        ("squall2d", "configs/squall2d.ini", "build-cpu-dp", "fp64", 6.0, 300.0, 600.0),
        ("squall2d", "configs/squall2d.ini", "build-cuda", "fp32", 6.0, 300.0, 600.0),
        ("squall2d", "configs/squall2d.ini", "build-cuda-dp", "fp64", 6.0, 300.0, 600.0),
    ]

    for c in configs:
        run_test(*c)

    print("\n==================================================")
    print("ALL RESTART ACCEPTANCE TESTS PASSED BITWISE!")
    print("==================================================")

if __name__ == "__main__":
    main()
