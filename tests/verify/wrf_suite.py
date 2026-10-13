#!/usr/bin/env python3
"""WRF-identity regression suite: every CPPWRF configuration that has a WRF
reference run, with the tolerance that the current code meets. Run after every change to the
dycore, boundaries, diffusion or physics coupling.

  wrf_suite.py --cpu-dp <build>/cppwrf --cuda <build>/cppwrf [--cuda-dp <build>/cppwrf] [--out DIR] [--only NAME]

Exit status 0 = all PASS. Cases whose reference or binary is missing are SKIPped (reported).
"""
import argparse
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
WR = os.path.expanduser("~/wrf_real")
WREF = os.path.expanduser("~/wrf_ref")

# name, binary key, config, extra --set, WRF wrfout, time [s], {var: max rel L2 (whole domain)}
CASES = [
    ("real_r1_fp64_cpu", "cpu_dp", "configs/real_r1.ini", [], f"{WR}/run_R1_r8/wrfout_d01_2026-10-12_00:00:00", 300,
     {"U": 1e-12, "W": 1e-10, "PH": 1e-13, "THM": 1e-13, "QV": 1e-13, "MU": 1e-13}),
    ("real_r1_fp64_cuda", "cuda_dp", "configs/real_r1.ini", [], f"{WR}/run_R1_r8/wrfout_d01_2026-10-12_00:00:00", 300,
     {"U": 1e-12, "W": 1e-10, "PH": 1e-13, "THM": 1e-13, "QV": 1e-13, "MU": 1e-13}),
    ("real_r0_fp32_1h", "cuda", "configs/real_r0.ini", [], f"{WR}/run_R0/wrfout_d01_2026-10-12_00:00:00", 3600,
     {"U": 3e-4, "W": 1e-2, "THM": 2e-6, "QV": 5e-5}),
    ("real_r2_fp32_1h", "cuda", "configs/real_r2.ini", [], f"{WR}/run_R2/wrfout_d01_2026-10-12_00:00:00", 3600,
     {"U": 3e-4, "W": 1e-2, "THM": 2e-6, "QV": 5e-5}),
    ("real_r3_fp32_1h", "cuda", "configs/real_r3.ini", [], f"{WR}/run_R3/wrfout_d01_2026-10-12_00:00:00", 3600,
     {"U": 3e-4, "W": 1e-2, "THM": 2e-6, "QV": 5e-5}),
    ("wdamp_w1_fp64", "cpu_dp", "configs/real_r3.ini",
     ["real.dir=" + f"{WR}/run_W1_r8", "time.dt=450", "run.run_time=4500", "output.out_interval=1",
      "output.diag_interval=1"], f"{WR}/run_W1_r8/wrfout_d01_2026-10-12_00:00:00", 4500,
     {"U": 1e-11, "W": 1e-10, "THM": 1e-13}),
]


def run(binary, cfg, sets, out):
    shutil.rmtree(out, ignore_errors=True)
    cmd = [binary, os.path.join(REPO, cfg), "--set", f"output.dir={out}"]
    for s in sets:
        cmd += ["--set", s]
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "3"))
    r = subprocess.run(cmd, cwd=REPO, env=env, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def compare(wrfout, out, t):
    r = subprocess.run([sys.executable, os.path.join(REPO, "tests/verify/compare_real.py"), wrfout, out],
                       capture_output=True, text=True)
    res = {}
    for line in r.stdout.splitlines():
        f = line.split()
        if len(f) >= 5 and f[0].replace(".", "").isdigit() and abs(float(f[0]) - t) < 1e-6:
            res[f[1]] = float(f[3])
    return res


def straka(binary, out):
    wo = f"{WREF}/run_straka_d1/wrfout_d01_0001-01-01_00:00:00"
    if not os.path.exists(wo):
        return None, "reference missing"
    rc, log = run(binary, "configs/straka_wrf.ini", [], out)
    if rc != 0:
        return False, log[-300:]
    r = subprocess.run([sys.executable, os.path.join(REPO, "tests/verify/verify.py"), "wrf", wo, out],
                       capture_output=True, text=True)
    return ("FAIL" not in r.stdout and r.returncode == 0), r.stdout.strip().splitlines()[-1] if r.stdout else r.stderr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpu-dp")
    ap.add_argument("--cuda")
    ap.add_argument("--cuda-dp")
    ap.add_argument("--out", default=os.path.expanduser("~/cppwrf_scratch/wrf_suite"))
    ap.add_argument("--only")
    a = ap.parse_args()
    bins = {"cpu_dp": a.cpu_dp, "cuda": a.cuda, "cuda_dp": a.cuda_dp}
    ok = True
    for name, key, cfg, sets, wo, t, tol in CASES:
        if a.only and a.only not in name:
            continue
        b = bins.get(key)
        if not b or not os.path.exists(b) or not os.path.exists(wo):
            print(f"SKIP  {name:22s} (binary or reference missing)")
            continue
        out = os.path.join(a.out, name)
        rc, log = run(b, cfg, sets, out)
        if rc != 0:
            print(f"FAIL  {name:22s} run failed: {log[-200:]}")
            ok = False
            continue
        res = compare(wo, out, t)
        bad = [f"{v} {res.get(v, float('nan')):.2e}>{tl:.0e}" for v, tl in tol.items()
               if not (res.get(v, float("inf")) <= tl)]
        txt = " ".join(f"{v}={res.get(v, float('nan')):.1e}" for v in tol)
        print(f"{'PASS' if not bad else 'FAIL'}  {name:22s} t={t:g}s {txt}" + (f"  [{'; '.join(bad)}]" if bad else ""))
        ok &= not bad
    if (not a.only or "straka" in a.only) and a.cpu_dp and os.path.exists(a.cpu_dp):
        good, msg = straka(a.cpu_dp, os.path.join(a.out, "straka"))
        if good is None:
            print(f"SKIP  straka_wrf             ({msg})")
        else:
            print(f"{'PASS' if good else 'FAIL'}  straka_wrf             {msg}")
            ok &= good
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
