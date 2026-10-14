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
    # P1 column physics (sfclayrev + Noah + YSU): 3 steps at roundoff in FP64; over 1 h the land
    # surface amplifies roundoff (CPU vs CUDA FP64 differ by the same amount), FP32 as the dycore.
    ("real_p1_r8_fp64_cpu", "cpu_dp", "configs/real_p1_r8.ini", [], f"{WR}/run_P1_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 1e-13, "W": 1e-10, "PH": 1e-14, "THM": 1e-14, "QV": 1e-14, "MU": 1e-14}),
    ("real_p1_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P1_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 3e-9, "W": 1e-7, "THM": 1e-9, "QV": 1e-8}),
    ("real_p1_fp32_1h", "cuda", "configs/real_p1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P1", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P1/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 3e-4, "W": 1e-2, "THM": 2e-6, "QV": 1e-4}),
    # P1 + WSM6 + diff_opt 1/km_opt 4 (run P5_r8): step 1 at roundoff. From step 2 on
    # WSM6 amplifies input roundoff ~1e8 (CPU vs CUDA FP64 QICE 7e-5): only step 1 is checked.
    ("real_p5_r8_fp64_cpu", "cpu_dp", "configs/real_p5_r8.ini", ["run.run_time=60"],
     f"{WR}/run_P5_r8/wrfout_d01_2026-10-12_00:00:00", 60,
     {"U": 1e-13, "W": 1e-10, "THM": 1e-14, "QV": 1e-14, "MU": 1e-14, "QCLOUD": 1e-10, "QICE": 1e-10}),
    # P2 (+Kain-Fritsch cu_physics=1): 3 steps at roundoff in FP64; at step 6 (t=360s) a 2-layer
    # saturated USL at (j=77,i=8) with WKL<1e-4 has TLCL-TENV == 0 algebraically (1-ULP trigger flip).
    ("real_p2_r8_fp64_cpu", "cpu_dp", "configs/real_p2_r8.ini", [], f"{WR}/run_P2_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 1e-13, "W": 1e-10, "PH": 1e-14, "THM": 1e-14, "QV": 1e-14, "MU": 1e-14}),
    ("real_p2_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p2_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P2_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P2_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 2e-4, "W": 5e-3, "THM": 3e-6, "QV": 1e-4}),
    ("real_p2_fp32_1h", "cuda", "configs/real_p2_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P2", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P2/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 2e-3, "W": 8e-2, "THM": 3e-5, "QV": 8e-4}),
    # P3 = P1 + RRTMG LW/SW every step (radt 1, o3input 2, ghg_input 1); tables in ~/wrf_real/rrtmg_tables
    ("real_p3_r8_fp64_cpu", "cpu_dp", "configs/real_p3_r8.ini", [], f"{WR}/run_P3_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 1e-13, "W": 1e-10, "PH": 1e-14, "THM": 1e-14, "QV": 1e-14, "MU": 1e-14}),
    ("real_p3_fp32_1h", "cuda", "configs/real_p3_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P3", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P3/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 3e-4, "W": 1e-2, "THM": 2e-6, "QV": 1e-4}),
    # P4 (KF + WSM6 + RRTMG with clouds, effective radii): step 1 at roundoff except one KF column where
    # a cancelling 1e-25 qr at cloud top has the opposite sign of WRF's and WSM6 turns it into a 1e-7
    # QRAIN difference (docs/PHYSICS_SUITE.md section 9); cloudy RRTMG is checked column by column.
    ("real_p4_r8_fp64_cpu", "cpu_dp", "configs/real_p4_r8.ini", ["run.run_time=60"],
     f"{WR}/run_P4_r8/steps/wrfout_d01_2026-10-12_00:00:00", 60,
     {"U": 1e-13, "W": 1e-10, "PH": 1e-14, "THM": 1e-11, "QV": 1e-9, "MU": 1e-14}),
    ("real_p4_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p4_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P4_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P4_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 6e-4, "W": 1e-2, "THM": 8e-6, "QV": 5e-4}),
    # Nesting: one-way (N1) and two-way (N2, feedback changes d01 by 5e-3) vs WRF r8
    ("real_n1_r8_fp64_cpu", "cpu_dp", "configs/real_n1_steps.ini", [],
     {d: f"{WR}/run_N1_r8/steps/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 180,
     {"d01.U": 1e-13, "d01.THM": 1e-14, "d02.U": 1e-12, "d02.THM": 1e-14, "d02.QV": 1e-13}),
    ("real_n2_r8_fp64_cpu", "cpu_dp", "configs/real_n2_steps.ini", [],
     {d: f"{WR}/run_N2_r8/steps/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 180,
     {"d01.U": 1e-13, "d01.THM": 1e-14, "d02.U": 1e-12, "d02.THM": 1e-14, "d02.QV": 1e-13}),
    ("real_n2_r8_fp64_cuda_1h", "cuda_dp", "configs/real_n2.ini", ["real.dir=" + f"{WR}/run_N2_r8"],
     {d: f"{WR}/run_N2_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {"d01.U": 1e-11, "d01.THM": 1e-13, "d02.U": 1e-11, "d02.THM": 1e-13}),
    ("real_n2_fp32_1h", "cuda", "configs/real_n2.ini", [],
     {d: f"{WR}/run_N2/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {"d01.U": 5e-4, "d01.THM": 3e-6, "d02.U": 6e-4, "d02.THM": 3e-6}),
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


def restart_p1(binary, out):
    """Physics restart: 3 steps straight vs restarted after step 1 must be bit-identical
    (dycore, surface/soil state, R*BLTEN, lagged t_1)."""
    if not os.path.exists(f"{WR}/run_P1_r8/steps/wrfinput_d01"):
        return None, "input missing"
    a, b = os.path.join(out, "straight"), os.path.join(out, "restarted")
    rc, log = run(binary, "configs/real_p1_r8.ini", ["restart.interval=60"], a)
    if rc != 0:
        return False, log[-300:]
    rc, log = run(binary, "configs/real_p1_r8.ini", [f"run.restart_file={a}/restart_000001.nc"], b)
    if rc != 0:
        return False, log[-300:]
    same = open(f"{a}/snap_000003.bin", "rb").read() == open(f"{b}/snap_000003.bin", "rb").read()
    return same, "step 3 bit-identical" if same else "step 3 differs after restart"


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
        refs = wo if isinstance(wo, dict) else {"": wo}  # nested runs: {domain: wrfout}
        if not b or not os.path.exists(b) or not all(os.path.exists(w) for w in refs.values()):
            print(f"SKIP  {name:22s} (binary or reference missing)")
            continue
        out = os.path.join(a.out, name)
        rc, log = run(b, cfg, sets, out)
        if rc != 0:
            print(f"FAIL  {name:22s} run failed: {log[-200:]}")
            ok = False
            continue
        res = {}
        for dom, w in refs.items():  # nested: out/<domain>, tolerance keys "<domain>.<var>"
            for v, x in compare(w, os.path.join(out, dom) if dom else out, t).items():
                res[f"{dom}.{v}" if dom else v] = x
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
    if (not a.only or "restart" in a.only) and a.cpu_dp and os.path.exists(a.cpu_dp):
        good, msg = restart_p1(a.cpu_dp, os.path.join(a.out, "restart_p1"))
        if good is None:
            print(f"SKIP  restart_p1_fp64        ({msg})")
        else:
            print(f"{'PASS' if good else 'FAIL'}  restart_p1_fp64        {msg}")
            ok &= good
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
