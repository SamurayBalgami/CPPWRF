#!/usr/bin/env python3
"""WRF-identity regression suite: every CPPWRF configuration that has a WRF
reference run, with the tolerance that the current code meets. Run after every change to the
dycore, boundaries, diffusion or physics coupling.

  wrf_suite.py --cpu-dp <build>/cppwrf --cuda <build>/cppwrf [--cuda-dp <build>/cppwrf] [--out DIR] [--only NAME[,NAME]]
               [--exact] [--jobs N]   (N CPU cases at once, each with OMP_NUM_THREADS threads) [--allow-skip] [--slow]

Cases with tolerance 0 (CPU FP64) run with WRF history output and their wrfout files are compared with WRF's at
tolerance 0 (exact.py: NaN equals only NaN, -0 differs from +0): every variable of the reference file, every
reference frame up to the case's time (the final one must be there), types and variable attributes, and the
global attributes when the reference's namelist.input is there (CPPWRF's history then reads it). The other
cases compare the model state with the reference's by relative L2 at the case's time (compare_real.py).
With the reference's namelist.input a tolerance-0 case FAILs when its configuration (config file and --set) sets
an option otherwise than tools/wrf2cppwrf.py's translation of that namelist.
Only the cases of the given binaries are selected; a given binary that does not exist is an error (exit 2).
Exit status 0 = every selected case PASSed; a case whose reference is missing is a SKIP and fails the suite unless
--allow-skip. The last line counts the results.
Slow cases (suite.slow in their --set list: the 24 h and 72 h runs) are left out unless --slow is given or --only --exact
names them.
"""
import argparse
import filecmp
import glob
import os
import re
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
# WRF reference runs (wrfinput/wrfbdy and wrfout of each case); the configs read them as ${CPPWRF_WRF_RUNS}
WR = os.environ.setdefault("CPPWRF_WRF_RUNS", os.path.expanduser("~/wrf_real480"))
os.environ.setdefault("CPPWRF_WRF_REF", os.path.expanduser("~/wrf_ref480"))
WREF = os.environ["CPPWRF_WRF_REF"]

# name, binary key, config, extra --set, WRF wrfout, time [s], {var: max rel L2 (whole domain)}; all 0.0: compared
# exactly in every variable of the wrfout file (the list names the variables that matter most for the case)
CASES = [
    # CPU FP64 real-data dycore: bit-identical to WRF (every field, every step; libmvec lanes, WRF's
    # summation orders and halo schedule reproduced)
    ("real_r1_fp64_cpu", "cpu_dp", "configs/real_r1.ini", [], f"{WR}/run_R1_r8/wrfout_d01_2026-10-12_00:00:00", 300,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    # R1D6 = R1 + diff_opt = 1, km_opt = 4, damp_opt = 1 (cal_dampkm), diff_6th_opt = 2 (u, v, w, theta, qv)
    ("real_r1d6_fp64_cpu", "cpu_dp", "configs/real_r1_d6.ini", ["real.dir=" + f"{WR}/run_R1D6_r8"],
     f"{WR}/run_R1D6_r8/wrfout_d01_2026-10-12_00:00:00", 600,
     {"U": 0.0, "V": 0.0, "W": 0.0, "PH": 0.0, "MU": 0.0, "P": 0.0, "QV": 0.0, "THM": 0.0}),
    ("real_r1d6_fp64_cuda", "cuda_dp", "configs/real_r1_d6.ini", ["real.dir=" + f"{WR}/run_R1D6_r8"],
     f"{WR}/run_R1D6_r8/wrfout_d01_2026-10-12_00:00:00", 600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_r1_fp64_cuda", "cuda_dp", "configs/real_r1.ini", [], f"{WR}/run_R1_r8/wrfout_d01_2026-10-12_00:00:00", 300,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    # T1 = R1 + diff_opt = 2, km_opt = 2 (1.5-order TKE, full metric stress tensor diffusion):
    # bit-identical on CPU FP64 across all fields including TKE
    ("real_t1_fp64_cpu", "cpu_dp", "configs/real_t1.ini", ["real.dir=" + f"{WR}/run_T1_r8"],
     f"{WR}/run_T1_r8/wrfout_d01_2026-10-12_00:00:00", 300,
     {"U": 0.0, "V": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0, "TKE": 0.0}),
    ("real_t1_fp64_cuda", "cuda_dp", "configs/real_t1.ini", ["real.dir=" + f"{WR}/run_T1_r8"], f"{WR}/run_T1_r8/wrfout_d01_2026-10-12_00:00:00", 300,
     {"U": 0.0, "V": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0, "TKE": 0.0}),
    ("real_r0_fp32_1h", "cuda", "configs/real_r0.ini", [], f"{WR}/run_R0/wrfout_d01_2026-10-12_00:00:00", 3600,
     {"U": 3e-4, "W": 1e-2, "THM": 2e-6, "QV": 5e-5}),
    ("real_r2_fp32_1h", "cuda", "configs/real_r2.ini", [], f"{WR}/run_R2/wrfout_d01_2026-10-12_00:00:00", 3600,
     {"U": 3e-4, "W": 1e-2, "THM": 2e-6, "QV": 5e-5}),
    ("real_r3_fp32_1h", "cuda", "configs/real_r3.ini", [], f"{WR}/run_R3/wrfout_d01_2026-10-12_00:00:00", 3600,
     {"U": 3e-4, "W": 1e-2, "THM": 2e-6, "QV": 5e-5}),
    # P1 column physics (sfclayrev + Noah + YSU): 3 steps at roundoff in FP64; over 1 h the land
    # surface amplifies roundoff (CPU vs CUDA FP64 differ by the same amount), FP32 as the dycore.
    ("real_p1_r8_fp64_cpu", "cpu_dp", "configs/real_p1_r8.ini", [], f"{WR}/run_P1_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_p1_r8_fp64_cpu_1h", "cpu_dp", "configs/real_p1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P1_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # P1T2 = P1 with opt_thcnd = 2 (Noah McCumber-Pielke conductivity on soil types 3 and 4)
    ("real_p1t2_r8_fp64_cpu_1h", "cpu_dp", "configs/real_p1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P1T2_r8", "phys.opt_thcnd=2", "output.out_interval=10",
      "output.diag_interval=10"],
     f"{WR}/run_P1T2_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_p1t2_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P1T2_r8", "phys.opt_thcnd=2", "output.out_interval=10",
      "output.diag_interval=10"],
     f"{WR}/run_P1T2_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # A01 = P1 with usemonalb = .true. (landuse_init keeps the wrfinput ALBBCK, SFLX too) and rdlai2d = .true. (Noah's
    # LAI from wrfinput); A01V = V2 with both (722 snow-covered land points: ALBEDO = SNOALB in landuse_init).
    # The A01 references leave out RQCBLTEN / RQIBLTEN (WRF's uninitialised arrays without qc / qi).
    ("real_a01_r8_fp64_cpu", "cpu_dp", "configs/real_a01_r8.ini", [],
     f"{WR}/run_A01_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_a01_r8_fp64_cpu_1h", "cpu_dp", "configs/real_a01_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P1_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_A01_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_a01_r8_fp64_cuda_1h", "cuda_dp", "configs/real_a01_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P1_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_A01_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_a01v_r8_fp64_cpu_1h", "cpu_dp", "configs/real_v2_r8.ini",
     ["run.run_time=3600", "output.out_interval=10", "output.diag_interval=10", "phys.usemonalb=true",
      "phys.rdlai2d=true"],
     f"{WR}/run_A01V_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {"U": 0.0, "V": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QICE": 0.0, "QSNOW": 0.0}),
    ("real_a01v_r8_fp64_cuda_1h", "cuda_dp", "configs/real_v2_r8.ini",
     ["run.run_time=3600", "output.out_interval=10", "output.diag_interval=10", "phys.usemonalb=true",
      "phys.rdlai2d=true"],
     f"{WR}/run_A01V_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_p2_r8_fp64_cpu_1h", "cpu_dp", "configs/real_p2_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P2_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P2_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_p7_r8_fp64_cpu_1h", "cpu_dp", "configs/real_p7_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P7_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P7_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_p8_r8_fp64_cpu_1h", "cpu_dp", "configs/real_p8_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P8_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P8_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_p1_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P1_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_p1_fp32_1h", "cuda", "configs/real_p1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P1", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P1/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 3e-4, "W": 1e-2, "THM": 2e-6, "QV": 1e-4}),
    # P5 = P1 + WSM6 + diff_opt 1/km_opt 4: bit-identical for 3 steps and 1 h
    ("real_p5_r8_fp64_cpu", "cpu_dp", "configs/real_p5_r8.ini", [],
     f"{WR}/run_P5_r8/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0, "QCLOUD": 0.0, "QICE": 0.0}),
    ("real_p5_r8_fp64_cpu_1h", "cpu_dp", "configs/real_p5_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P5_r8_1h", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P5_r8_1h/wrfout_d01_2026-10-12_00:00:00", 3600,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0, "QICE": 0.0, "QSNOW": 0.0}),
    ("real_p5_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p5_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P5_r8_1h", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P5_r8_1h/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # P2 (+Kain-Fritsch cu_physics=1): bit-identical (3 steps here, 1 h above)
    ("real_p2_r8_fp64_cpu", "cpu_dp", "configs/real_p2_r8.ini", [], f"{WR}/run_P2_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_p2_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p2_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P2_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P2_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_p2_fp32_1h", "cuda", "configs/real_p2_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P2", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P2/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 2e-3, "W": 8e-2, "THM": 3e-5, "QV": 8e-4}),
    # P3 = P1 + RRTMG LW/SW every step (radt 1, o3input 2, ghg_input 1); tables in ${CPPWRF_WRF_RUNS}/rrtmg_tables
    ("real_p3_r8_fp64_cpu", "cpu_dp", "configs/real_p3_r8.ini", [], f"{WR}/run_P3_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_p3_fp32_1h", "cuda", "configs/real_p3_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P3", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P3/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 3e-4, "W": 1e-2, "THM": 2e-6, "QV": 1e-4}),
    ("real_p3_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p3_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P3_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P3_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # P4 (KF + WSM6 + RRTMG with clouds, effective radii + diffusion): bit-identical for 3 steps and 1 h
    # (only the COSZEN diagnostic differs: radiation_driver's vector cos/sin)
    ("real_p4_r8_fp64_cpu", "cpu_dp", "configs/real_p4_r8.ini", [],
     f"{WR}/run_P4_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_p4_r8_fp64_cpu_1h", "cpu_dp", "configs/real_p4_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P4_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P4_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0, "QICE": 0.0}),
    ("real_p4_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p4_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P4_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P4_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_p4_fp32_1h", "cuda", "configs/real_p4_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P4", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P4/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 3e-3, "W": 7e-2, "THM": 3e-5, "QV": 2e-3}),
    # A05 = P4 with cu_rad_feedback (Kain-Fritsch's CLDFRA_DP / CLDFRA_SH and QC_CU / QI_CU into RRTMG, icloud_cu 2),
    # A05C = P4 with cldovrlp = 0 (RRTMG without clouds), A05G = P8 + WSM6 + RRTMG with cu_rad_feedback (the
    # Grell-Freitas QC_CU / QI_CU into qc / qi, icloud_cu 1): bit-identical for 3 steps and 1 h. CUDA FP64 tolerances:
    # 3x the relative L2 on Colab G4 (A05: U 6.4e-4, W 2.5e-2, THM 9.9e-6, QV 5.8e-4; A05C: 1.6e-5, 5.5e-4, 2.7e-7,
    # 1.9e-5; A05G: 5.4e-4, 1.3e-2, 7.7e-6, 2.7e-4: the KF / GF thresholds amplify the GPU's roundoff as in P4)
    *[c for tag, ref, gtol in (
        ("a05", "A05", {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
        ("a05c", "A05C", {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
        ("a05g", "A05G", {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0})) for c in (
        (f"real_{tag}_r8_fp64_cpu", "cpu_dp", f"configs/real_{tag}_r8.ini", [],
         f"{WR}/run_{ref}_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
         {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
        (f"real_{tag}_r8_fp64_cpu_1h", "cpu_dp", f"configs/real_{tag}_r8.ini",
         ["run.run_time=3600", "real.dir=" + f"{WR}/run_{ref}_r8", "output.out_interval=10", "output.diag_interval=10"],
         f"{WR}/run_{ref}_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
         {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0, "QICE": 0.0}),
        (f"real_{tag}_r8_fp64_cuda_1h", "cuda_dp", f"configs/real_{tag}_r8.ini",
         ["run.run_time=3600", "real.dir=" + f"{WR}/run_{ref}_r8", "output.out_interval=10", "output.diag_interval=10"],
         f"{WR}/run_{ref}_r8/wrfout_d01_2026-10-12_00:00:00", 3600, gtol))],
    # B01 = P4 with the adaptive time step (use_adaptive_time_step: dt from the CFL numbers of each step, stepping to the
    # history times; 1 h: steps of 60 to 120 s): bit-identical for 3 minutes (history every minute) and 1 h
    ("real_b01_r8_fp64_cpu", "cpu_dp", "configs/real_b01_r8.ini", [],
     f"{WR}/run_B01_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_b01_r8_fp64_cpu_1h", "cpu_dp", "configs/real_b01_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P4_r8", "time.history_interval_s=600"],
     f"{WR}/run_B01_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0, "QICE": 0.0}),
    ("real_b01_r8_fp64_cuda_1h", "cuda_dp", "configs/real_b01_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P4_r8", "time.history_interval_s=600"],
     f"{WR}/run_B01_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),  # G4: U 2.2e-4 W 6.2e-3 THM 4.5e-6 QV 2.0e-4
    # A08 = P4 with bucket_mm 0.05, bucket_J 1e4 and prec_acc_dt 2 min (diagnostic_output_calc: the buckets of RAINC /
    # RAINNC and the LW accumulations tip at the first step, I_RAINC ... I_ACLWDNBC; PREC_ACC_C / NC, SNOW_ACC_NC)
    ("real_a08_r8_fp64_cpu", "cpu_dp", "configs/real_a08_r8.ini", [],
     f"{WR}/run_A08_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "RAINC": 0.0, "RAINNC": 0.0, "I_RAINC": 0.0, "I_RAINNC": 0.0,
      "ACLWUPT": 0.0, "I_ACLWUPT": 0.0, "PREC_ACC_C": 0.0, "PREC_ACC_NC": 0.0, "SNOW_ACC_NC": 0.0}),
    ("real_a08_r8_fp64_cpu_1h", "cpu_dp", "configs/real_a08_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P4_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_A08_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "RAINC": 0.0, "RAINNC": 0.0, "I_RAINC": 0.0, "I_RAINNC": 0.0,
      "ACLWUPT": 0.0, "I_ACLWUPT": 0.0, "PREC_ACC_C": 0.0, "PREC_ACC_NC": 0.0, "SNOW_ACC_NC": 0.0}),
    ("real_a08_r8_fp64_cuda_1h", "cuda_dp", "configs/real_a08_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P4_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_A08_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),  # G4: U 1.8e-4 W 4.6e-3 THM 3.1e-6 QV 1.0e-4
    # A09 = P4 with output_diagnostics = 1, auxhist3_interval = 3 min (clwrf_output_calc): the wrfout files and the
    # auxhist3 files (T2MIN ... RAINNCVSTD) exact; the 1 h run also restarted at 30 min (wrfrst with the statistics)
    ("real_a09_r8_fp64_cpu", "cpu_dp", "configs/real_a09_r8.ini", ["output.auxhist3=true"],
     f"{WR}/run_A09_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "T2": 0.0, "Q2": 0.0, "U10": 0.0, "V10": 0.0, "RAINNCV": 0.0}),
    ("real_a09_r8_fp64_cpu_1h", "cpu_dp", "configs/real_a09_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_A09_r8", "output.out_interval=10", "output.diag_interval=10",
      "output.auxhist3=true"],
     f"{WR}/run_A09_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "T2": 0.0, "Q2": 0.0, "U10": 0.0, "V10": 0.0, "RAINNCV": 0.0}),
    ("real_a09_r8_fp64_cpu_1h_rst", "cpu_dp", "configs/real_a09_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_A09_r8", "output.out_interval=10", "output.diag_interval=10",
      "output.auxhist3=true", "suite.restart_at=1800"],
     f"{WR}/run_A09_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "T2": 0.0, "Q2": 0.0, "U10": 0.0, "V10": 0.0, "RAINNCV": 0.0}),
    ("real_a09_r8_fp64_cuda_1h", "cuda_dp", "configs/real_a09_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_A09_r8", "output.out_interval=10", "output.diag_interval=10",
      "output.auxhist3=true"],
     f"{WR}/run_A09_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # C01 = P6 with the aerosol-aware Thompson (mp_physics = 28, aer_init_opt = 0; WRF r8d8): bit-identical for
    # 3 minutes (history every minute) and 1 h
    ("real_c01_r8_fp64_cpu", "cpu_dp", "configs/real_c01_r8.ini", [],
     f"{WR}/run_C01_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {v: 0.0 for v in ("U", "W", "PH", "MU", "THM", "QV", "QCLOUD", "QNCLOUD", "QNWFA", "QNIFA", "QNWFA2D")}),
    ("real_c01_r8_fp64_cpu_1h", "cpu_dp", "configs/real_c01_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_C01_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {v: 0.0 for v in ("U", "W", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QNICE", "QNRAIN", "QNCLOUD", "QNWFA",
                       "QNIFA", "RAINNC")}),
    ("real_c01_r8_fp64_cuda_1h", "cuda_dp", "configs/real_c01_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_C01_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),  # G4: U 3.6e-7 W 7.3e-5 THM 1.1e-8 QV 1.6e-7
    # RADW / RADT = P4 (WSM6) / P6 (Thompson) with do_radar_ref = 1 (REFL_10CM at the history times; WRF r8d8):
    # bit-identical for 3 minutes (history every minute) and 1 h
    ("real_radw_r8_fp64_cpu", "cpu_dp", "configs/real_radw_r8.ini", [],
     f"{WR}/run_RADW_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {v: 0.0 for v in ("U", "W", "THM", "QV", "QRAIN", "QSNOW", "QGRAUP", "REFL_10CM")}),
    ("real_radw_r8_fp64_cpu_1h", "cpu_dp", "configs/real_radw_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P4_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_RADW_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {v: 0.0 for v in ("U", "W", "THM", "QV", "QRAIN", "QSNOW", "QGRAUP", "REFL_10CM")}),
    ("real_radw_r8_fp64_cuda_1h", "cuda_dp", "configs/real_radw_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P4_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_RADW_r8/wrfout_d01_2026-10-12_00:00:00", 3600,  # G4: U 1.8e-4 W 4.6e-3 THM 3.1e-6 QV 1.0e-4 REFL_10CM 3.4e-3
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "REFL_10CM": 0.0}),
    ("real_radt_r8_fp64_cpu", "cpu_dp", "configs/real_radt_r8.ini", [],
     f"{WR}/run_RADT_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {v: 0.0 for v in ("U", "W", "THM", "QV", "QRAIN", "QSNOW", "QGRAUP", "REFL_10CM")}),
    ("real_radt_r8_fp64_cpu_1h", "cpu_dp", "configs/real_radt_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_RADT_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {v: 0.0 for v in ("U", "W", "THM", "QV", "QRAIN", "QSNOW", "QGRAUP", "QNRAIN", "REFL_10CM")}),
    ("real_radt_r8_fp64_cuda_1h", "cuda_dp", "configs/real_radt_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_RADT_r8/wrfout_d01_2026-10-12_00:00:00", 3600,  # G4: U 3.2e-5 W 9.0e-4 THM 4.0e-7 QV 8.5e-6 REFL_10CM 3.9e-4
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "REFL_10CM": 0.0}),
    # KESSKF = P4 with Kessler (mp_physics = 1): Kain-Fritsch rain into Kessler's sedimentation (nfall > 1),
    # bit-identical for 3 steps and 1 h
    ("real_kesskf_r8_fp64_cpu", "cpu_dp", "configs/real_kesskf_r8.ini", [],
     f"{WR}/run_KESSKF_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "RAINNC")}),
    ("real_kesskf_r8_fp64_cpu_1h", "cpu_dp", "configs/real_kesskf_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P4_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_KESSKF_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "RAINNC")}),
    ("real_kesskf_r8_fp64_cuda_1h", "cuda_dp", "configs/real_kesskf_r8.ini",  # G4: U 9.2e-5 W 1.2e-3 THM 1.2e-6 QV 6.2e-5
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P4_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_KESSKF_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # DUD = V1 (Bay of Bengal, 00 UTC: sunrise over the domain) with the Dudhia shortwave (ra_sw_physics = 1) and
    # RRTMG longwave: bit-identical for 3 steps and 1 h (P4's domain is in the night for the first hours)
    ("real_dud_r8_fp64_cpu", "cpu_dp", "configs/real_dud_r8.ini", [],
     f"{WR}/run_DUD_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "THM", "QV", "QCLOUD", "GSW", "SWDOWN", "GLW")}),
    ("real_dud_r8_fp64_cpu_1h", "cpu_dp", "configs/real_dud_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_V1_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_DUD_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "GSW", "SWDOWN", "GLW")}),
    ("real_dud_r8_fp64_cuda_1h", "cuda_dp", "configs/real_dud_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_V1_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_DUD_r8/wrfout_d01_2026-10-12_00:00:00", 3600,  # G4: U 5.0e-4 W 4.0e-2 THM 4.4e-6 QV 9.9e-5
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # RRTM: V1 with the RRTM longwave (module_ra_rrtm.F) and RRTMG shortwave: 3 steps and 1 h bit-identical
    ("real_rrtm_r8_fp64_cpu", "cpu_dp", "configs/real_rrtm_r8.ini", [],
     f"{WR}/run_RRTM_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "THM", "QV", "QCLOUD", "GSW", "SWDOWN", "GLW", "OLR")}),
    ("real_rrtm_r8_fp64_cpu_1h", "cpu_dp", "configs/real_rrtm_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_V1_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_RRTM_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "GSW", "GLW", "OLR")}),
    ("real_rrtm_r8_fp64_cuda_1h", "cuda_dp", "configs/real_rrtm_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_V1_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_RRTM_r8/wrfout_d01_2026-10-12_00:00:00", 3600,  # G4: U 1.5e-4 W 9.8e-3 THM 1.6e-6 QV 3.6e-5
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # P6 = P5 with Thompson (mp_physics = 8; WRF r8 built with -fdefault-double-8, see PHYSICS_SUITE.md
    # section 11): bit-identical for 3 steps and 1 h
    ("real_p6_r8_fp64_cpu", "cpu_dp", "configs/real_p6_r8.ini", [], f"{WR}/run_P6_r8/steps/wrfout_d01_2026-10-12_00:00:00",
     180, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0}),
    ("real_p6_r8_fp64_cpu_1h", "cpu_dp", "configs/real_p6_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P6_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0, "QICE": 0.0, "QNICE": 0.0,
      "QNRAIN": 0.0}),
    ("real_p6_fp32_1h", "cuda", "configs/real_p6_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P6/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 5e-4, "W": 1e-2, "THM": 1e-5, "QV": 5e-4}),
    ("real_p6_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p6_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P6_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # P7 = P2 with New Tiedtke (cu_physics = 16) instead of KF; reference WRF-4.8.0-r8d8 (r8 + double-8).
    # The scheme has thresholds that roundoff flips (one column at step 4 jumps 1.3e-4 K under 1e-15 input noise).
    ("real_p7_r8_fp64_cpu", "cpu_dp", "configs/real_p7_r8.ini", [],
     f"{WR}/run_P7_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_p7_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p7_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P7_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P7_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_p7_fp32_1h", "cuda", "configs/real_p7_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P7", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P7/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 2e-3, "W": 5e-2, "THM": 3e-5, "QV": 8e-4}),
    # P8 = P2 with Grell-Freitas-Li (cu_physics = 3, GFLSCHEME, ishallow = 1); reference WRF-4.8.0-r8d8.
    # Column replay of cu_grell_freitas_li through WRF's objects is bit-identical (tests/phys/cu_gfl).
    ("real_p8_r8_fp64_cpu_10min", "cpu_dp", "configs/real_p8_r8.ini", ["run.run_time=600"],
     f"{WR}/run_P8_r8/steps/wrfout_d01_2026-10-12_00:00:00", 600,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    # P8 CUDA: Grell-Freitas thresholds amplify the GPU's last-bit differences, which depend on the architecture
    # (deterministic per GPU): QV 5.4e-9 on L4 / RTX 2060, 1.5e-7 on G4 (Blackwell: nvcc's sm_120 code contracts
    # other products of the surface layer and Noah into FMAs; its compute_89 PTX or --fmad=false gives the L4 result)
    ("real_p8_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p8_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P8_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P8_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_p8_fp32_1h", "cuda", "configs/real_p8_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 2e-3, "W": 5e-2, "THM": 3e-5, "QV": 8e-4}),
    # P9 = P1 with MYJ PBL (bl_pbl_physics = 2) + Eta surface layer (sf_sfclay_physics = 2);
    # reference WRF real*8. Bit-identical on CPU FP64 for both 3 steps and 1 h.
    ("real_p9_r8_fp64_cpu", "cpu_dp", "configs/real_p9_r8.ini", [],
     f"{WR}/run_P9_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_p9_r8_fp64_cpu_1h", "cpu_dp", "configs/real_p9_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P9_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P9_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_p9_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p9_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P9_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P9_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_p9_fp32_1h", "cuda", "configs/real_p9_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P9", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P9/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 6e-4, "W": 1e-2, "THM": 1e-5, "QV": 2e-3}),
    # MYNN = P6 with the MYNN surface layer (sf_sfclay_physics = 5), the MYNN-EDMF PBL (bl_pbl_physics = 5,
    # icloud_bl = 1) and RRTMG every 10 min; reference WRF-4.8.0-r8d8. Bit-identical for 3 steps and 1 h.
    ("real_mynn_r8_fp64_cpu", "cpu_dp", "configs/real_mynn_r8.ini", [],
     f"{WR}/run_MYNN_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynn_r8_fp64_cpu_1h", "cpu_dp", "configs/real_mynn_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_MYNN_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNN_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynn_r8_fp64_cuda_1h", "cuda_dp", "configs/real_mynn_r8.ini",  # L4: U 2.6e-4 W 5.8e-3 THM 5.5e-6 QV 1.1e-4
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_MYNN_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNN_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # MYNNSF1 / MYNNSF2 = MYNN with the MYNN-EDMF PBL on the revised MM5 (sf_sfclay_physics = 1) / Eta (2) surface
    # layer instead of the MYNN one, from the P6 inputs; reference WRF-4.8.0-r8d8. Bit-identical for 3 steps and 1 h.
    ("real_mynnsf1_r8_fp64_cpu", "cpu_dp", "configs/real_mynnsf1_r8.ini", [],
     f"{WR}/run_MYNNSF1_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynnsf1_r8_fp64_cpu_1h", "cpu_dp", "configs/real_mynnsf1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNNSF1_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynnsf1_r8_fp64_cuda_1h", "cuda_dp", "configs/real_mynnsf1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNNSF1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),  # L4: U 2.5e-4 W 5.3e-3 THM 5.3e-6 QV 2.9e-4, G4: U 2.7e-4
    ("real_mynnsf2_r8_fp64_cpu", "cpu_dp", "configs/real_mynnsf2_r8.ini", [],
     f"{WR}/run_MYNNSF2_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynnsf2_r8_fp64_cpu_1h", "cpu_dp", "configs/real_mynnsf2_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNNSF2_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynnsf2_r8_fp64_cuda_1h", "cuda_dp", "configs/real_mynnsf2_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNNSF2_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),  # L4: U 2.6e-4 W 5.4e-3 THM 5.6e-6 QV 1.1e-4
    # MYNNX1 / MYNNX2 = MYNN (P6 inputs, MYNN-EDMF PBL and MYNN surface layer, RRTMG radt 10) with the non-default
    # MYNN options: X1 bl_mynn_tkeadvect, cloudpdf 1, mixlength 2, closure 3.0, edmf_tke 1, mixqt 1, output 1,
    # tke_budget 1, sf_mynn_sfcflux_land 2 / _water 3, grav_settling 2; X2 bl_mynn_edmf 0, cloudmix 0, cloudpdf 0,
    # mixlength 0, closure 2.5, ess 0, mixscalars 1, sf_mynn_sfcflux_land 3 / _water 4, grav_settling 1.
    # Reference WRF-4.8.0-r8d8; bit-identical for 3 steps and 1 h. CUDA FP64 L4 / G4: X1 U 3.2e-4 / 3.1e-4,
    # W 8.6e-3 / 5.9e-3, THM 9.0e-6 / 5.4e-6, QV 1.4e-4 / 1.3e-4; X2 U 2.1e-4, W 6.3e-3 / 7.3e-3, THM 4.4e-6 / 4.8e-6,
    # QV 6.4e-5 / 6.7e-5.
    ("real_mynnx1_r8_fp64_cpu", "cpu_dp", "configs/real_mynnx1_r8.ini", [],
     f"{WR}/run_MYNNX1_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynnx1_r8_fp64_cpu_1h", "cpu_dp", "configs/real_mynnx1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNNX1_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynnx1_r8_fp64_cuda_1h", "cuda_dp", "configs/real_mynnx1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNNX1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_mynnx2_r8_fp64_cpu", "cpu_dp", "configs/real_mynnx2_r8.ini", [],
     f"{WR}/run_MYNNX2_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynnx2_r8_fp64_cpu_1h", "cpu_dp", "configs/real_mynnx2_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNNX2_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                       "QNICE", "QNRAIN", "RAINNC")}),
    ("real_mynnx2_r8_fp64_cuda_1h", "cuda_dp", "configs/real_mynnx2_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P6_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_MYNNX2_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # V1-V3: other regions and map projections from the global GFS of 2026-10-12 00 UTC (WPS on the low
    # resolution global geography): V1 tropics (Bay of Bengal, Mercator, 15 km, P4 physics), V2 Arctic
    # (Svalbard, polar stereographic, 15 km, P4 physics: sea ice, glacier land ice, polar night), V3 Alps
    # (Lambert, 3 km, no cumulus, Thompson + MYJ + Eta + Noah-MP + RRTMG radt 5). Checked over 6 h (V1, V2)
    # and 2 h (V3); the suite runs 2 h and 30 min.
    ("real_v1_r8_fp64_cpu_2h", "cpu_dp", "configs/real_v1_r8.ini",
     ["run.run_time=7200", "output.out_interval=60", "output.diag_interval=60"],
     f"{WR}/run_V1_r8/wrfout_d01_2026-10-12_00:00:00", 7200,
     {"U": 0.0, "V": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0, "QICE": 0.0}),
    ("real_v2_r8_fp64_cpu_2h", "cpu_dp", "configs/real_v2_r8.ini",
     ["run.run_time=7200", "output.out_interval=60", "output.diag_interval=60"],
     f"{WR}/run_V2_r8/wrfout_d01_2026-10-12_00:00:00", 7200,
     {"U": 0.0, "V": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QICE": 0.0, "QSNOW": 0.0}),
    ("real_v3_r8_fp64_cpu_30min", "cpu_dp", "configs/real_v3_r8.ini",
     ["run.run_time=1800", "output.out_interval=100", "output.diag_interval=100"],
     f"{WR}/run_V3_r8/wrfout_d01_2026-10-12_00:00:00", 1800,
     {"U": 0.0, "V": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0, "QICE": 0.0,
      "QSNOW": 0.0, "QGRAUP": 0.0, "QNICE": 0.0, "QNRAIN": 0.0}),
    ("real_v1_r8_fp64_cuda_2h", "cuda_dp", "configs/real_v1_r8.ini",
     ["run.run_time=7200", "output.out_interval=60", "output.diag_interval=60"],
     f"{WR}/run_V1_r8/wrfout_d01_2026-10-12_00:00:00", 7200, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_v2_r8_fp64_cuda_2h", "cuda_dp", "configs/real_v2_r8.ini",
     ["run.run_time=7200", "output.out_interval=60", "output.diag_interval=60"],
     f"{WR}/run_V2_r8/wrfout_d01_2026-10-12_00:00:00", 7200, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_v3_r8_fp64_cuda_30min", "cuda_dp", "configs/real_v3_r8.ini",
     ["run.run_time=1800", "output.out_interval=100", "output.diag_interval=100"],
     f"{WR}/run_V3_r8/wrfout_d01_2026-10-12_00:00:00", 1800, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # The options of WRF's default namelist (test/em_real, physics_suite CONUS): T6 = P2 with Tiedtke
    # (cu_physics = 6), G2 = P2 with GWDO (gwd_opt = 1; ELVMAX = HGT + 3 VAR in the wrfinput, WPS 4.6 writes
    # none), F1 / F2 = V2 with fractional sea ice (real.exe run again with fractional_seaice = 1) and
    # sfclayrev + YSU / Eta + MYJ. Bit-identical in every wrfout variable over 3 h.
    ("real_t6_r8_fp64_cpu", "cpu_dp", "configs/real_t6_r8.ini", [],
     f"{WR}/run_T6_r8/steps/wrfout_d01_2026-10-12_00:00:00", 600,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_t6_r8_fp64_cpu_3h", "cpu_dp", "configs/real_t6_r8.ini",
     ["run.run_time=10800", "real.dir=" + f"{WR}/run_T6_r8", "output.out_interval=30", "output.diag_interval=30"],
     f"{WR}/run_T6_r8/wrfout_d01_2026-10-12_00:00:00", 10800, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_t6_r8_fp64_cuda_1h", "cuda_dp", "configs/real_t6_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_T6_r8", "output.out_interval=30", "output.diag_interval=30"],
     f"{WR}/run_T6_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_g2_r8_fp64_cpu", "cpu_dp", "configs/real_g2_r8.ini", [],
     f"{WR}/run_G2_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_g2_r8_fp64_cpu_3h", "cpu_dp", "configs/real_g2_r8.ini",
     ["run.run_time=10800", "real.dir=" + f"{WR}/run_G2_r8", "output.out_interval=30", "output.diag_interval=30"],
     f"{WR}/run_G2_r8/wrfout_d01_2026-10-12_00:00:00", 10800, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_g2_r8_fp64_cuda_1h", "cuda_dp", "configs/real_g2_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_G2_r8", "output.out_interval=30", "output.diag_interval=30"],
     f"{WR}/run_G2_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_f1_r8_fp64_cpu", "cpu_dp", "configs/real_f1_r8.ini", [],
     f"{WR}/run_F1_r8/steps/wrfout_d01_2026-10-12_00:00:00", 300,
     {"U": 0.0, "V": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QICE": 0.0, "QSNOW": 0.0}),
    ("real_f1_r8_fp64_cpu_3h", "cpu_dp", "configs/real_f1_r8.ini",
     ["run.run_time=10800", "real.dir=" + f"{WR}/run_F1_r8", "output.out_interval=60", "output.diag_interval=60"],
     f"{WR}/run_F1_r8/wrfout_d01_2026-10-12_00:00:00", 10800,
     {"U": 0.0, "V": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QICE": 0.0, "QSNOW": 0.0}),
    ("real_f2_r8_fp64_cpu", "cpu_dp", "configs/real_f2_r8.ini", [],
     f"{WR}/run_F2_r8/steps/wrfout_d01_2026-10-12_00:00:00", 300,
     {"U": 0.0, "V": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QICE": 0.0, "QSNOW": 0.0}),
    ("real_f2_r8_fp64_cpu_3h", "cpu_dp", "configs/real_f2_r8.ini",
     ["run.run_time=10800", "real.dir=" + f"{WR}/run_F2_r8", "output.out_interval=60", "output.diag_interval=60"],
     f"{WR}/run_F2_r8/wrfout_d01_2026-10-12_00:00:00", 10800,
     {"U": 0.0, "V": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QICE": 0.0, "QSNOW": 0.0}),
    ("real_f1_r8_fp64_cuda_1h", "cuda_dp", "configs/real_f1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_F1_r8", "output.out_interval=60", "output.diag_interval=60"],
     f"{WR}/run_F1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_f2_r8_fp64_cuda_1h", "cuda_dp", "configs/real_f2_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_F2_r8", "output.out_interval=60", "output.diag_interval=60"],
     f"{WR}/run_F2_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # D0 / D0V: WRF's default real-data namelist (physics_suite CONUS, fractional sea ice, gwd_opt 1, diff_opt 2 /
    # km_opt 4) on the P2 and V2 domains; Thompson references from WRF -r8 -double-8.
    ("real_d0_r8_fp64_cpu", "cpu_dp", "configs/real_d0_r8.ini", [],
     f"{WR}/run_D0_r8/steps/wrfout_d01_2026-10-12_00:00:00", 1200,
     {"U": 0.0, "V": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0, "QICE": 0.0,
      "QSNOW": 0.0, "QGRAUP": 0.0, "QNICE": 0.0, "QNRAIN": 0.0}),
    ("real_d0_r8_fp64_cpu_6h", "cpu_dp", "configs/real_d0_r8.ini",
     ["run.run_time=21600", "real.dir=" + f"{WR}/run_D0_r8", "output.out_interval=60", "output.diag_interval=60"],
     f"{WR}/run_D0_r8/wrfout_d01_2026-10-12_00:00:00", 21600,
     {"U": 0.0, "V": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0, "QICE": 0.0,
      "QSNOW": 0.0, "QGRAUP": 0.0, "QNICE": 0.0, "QNRAIN": 0.0}),
    ("real_d0_r8_fp64_cuda_1h", "cuda_dp", "configs/real_d0_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_D0_r8", "output.out_interval=60", "output.diag_interval=60"],
     f"{WR}/run_D0_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_d0v_r8_fp64_cpu", "cpu_dp", "configs/real_d0v_r8.ini", [],
     f"{WR}/run_D0V_r8/steps/wrfout_d01_2026-10-12_00:00:00", 1200,
     {"U": 0.0, "V": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0, "QICE": 0.0,
      "QSNOW": 0.0, "QGRAUP": 0.0, "QNICE": 0.0, "QNRAIN": 0.0}),
    ("real_d0v_r8_fp64_cpu_6h", "cpu_dp", "configs/real_d0v_r8.ini",
     ["run.run_time=21600", "real.dir=" + f"{WR}/run_D0V_r8", "output.out_interval=60", "output.diag_interval=60"],
     f"{WR}/run_D0V_r8/wrfout_d01_2026-10-12_00:00:00", 21600,
     {"U": 0.0, "V": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0, "QICE": 0.0,
      "QSNOW": 0.0, "QGRAUP": 0.0, "QNICE": 0.0, "QNRAIN": 0.0}),
    ("real_d0v_r8_fp64_cuda_1h", "cuda_dp", "configs/real_d0v_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_D0V_r8", "output.out_interval=60", "output.diag_interval=60"],
     f"{WR}/run_D0V_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # SU1 = D0V with sst_update = 1 (wrflowinp_d01 every 6 h; the 7 h run includes the 6 h update).
    ("real_su1_r8_fp64_cpu", "cpu_dp", "configs/real_su1_r8.ini", [],
     f"{WR}/run_SU1_r8/steps/wrfout_d01_2026-10-12_00:00:00", 1200,
     {"U": 0.0, "V": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0, "QICE": 0.0,
      "QSNOW": 0.0, "QGRAUP": 0.0, "QNICE": 0.0, "QNRAIN": 0.0}),
    # SU2 = SU1's steps with a changed SEAICE in wrflowinp at the start (water to ice, ice to water, new fractions).
    ("real_su2_r8_fp64_cpu", "cpu_dp", "configs/real_su1_r8.ini", ["real.dir=" + f"{WR}/run_SU2_r8/steps"],
     f"{WR}/run_SU2_r8/steps/wrfout_d01_2026-10-12_00:00:00", 1200,
     {"U": 0.0, "V": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0, "QICE": 0.0,
      "QSNOW": 0.0, "QGRAUP": 0.0, "QNICE": 0.0, "QNRAIN": 0.0}),
    ("real_su2_r8_fp64_cuda", "cuda_dp", "configs/real_su1_r8.ini", ["real.dir=" + f"{WR}/run_SU2_r8/steps"],
     f"{WR}/run_SU2_r8/steps/wrfout_d01_2026-10-12_00:00:00", 1200, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_su1_r8_fp64_cpu_7h", "cpu_dp", "configs/real_su1_r8.ini",
     ["run.run_time=25200", "real.dir=" + f"{WR}/run_SU1_r8", "output.out_interval=30", "output.diag_interval=30"],
     f"{WR}/run_SU1_r8/wrfout_d01_2026-10-12_00:00:00", 25200,
     {"U": 0.0, "V": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0, "QICE": 0.0,
      "QSNOW": 0.0, "QGRAUP": 0.0, "QNICE": 0.0, "QNRAIN": 0.0}),
    ("real_su1_r8_fp64_cuda_1h", "cuda_dp", "configs/real_su1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_SU1_r8", "output.out_interval=30", "output.diag_interval=30"],
     f"{WR}/run_SU1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # HY1 = P1 with hypsometric_opt = 1 (real.exe and wrf.exe namelists; wrf.exe takes the namelist value)
    ("real_hy1_r8_fp64_cpu", "cpu_dp", "configs/real_p1_r8.ini", ["real.dir=" + f"{WR}/run_HY1_r8/steps",
                                                                  "dyn.hypsometric_opt=1"],
     f"{WR}/run_HY1_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "V": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_hy1_r8_fp64_cpu_1h", "cpu_dp", "configs/real_p1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_HY1_r8", "dyn.hypsometric_opt=1", "output.out_interval=10",
      "output.diag_interval=10"],
     f"{WR}/run_HY1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_hy1_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_HY1_r8", "dyn.hypsometric_opt=1", "output.out_interval=10",
      "output.diag_interval=10"],
     f"{WR}/run_HY1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # SLAB1 / SLAB4 = P1 / P4 with the 5-layer thermal diffusion land surface (sf_surface_physics = 1), inputs from
    # real.exe with the slab's 5 soil layers
    ("real_slab1_r8_fp64_cpu", "cpu_dp", "configs/real_slab1_r8.ini", [],
     f"{WR}/run_SLAB1_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "V": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_slab1_r8_fp64_cpu_1h", "cpu_dp", "configs/real_slab1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_SLAB1_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_SLAB1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_slab1_r8_fp64_cuda_1h", "cuda_dp", "configs/real_slab1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_SLAB1_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_SLAB1_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),  # G4: U 4.4e-13 W 7.9e-12 THM 7.8e-16 QV 2.8e-14
    ("real_slab4_r8_fp64_cpu_1h", "cpu_dp", "configs/real_slab4_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_SLAB4_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_SLAB4_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0, "QICE": 0.0}),
    ("real_slab4_r8_fp64_cuda_1h", "cuda_dp", "configs/real_slab4_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_SLAB4_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_SLAB4_r8/wrfout_d01_2026-10-12_00:00:00", 3600,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),  # G4: U 3.4e-4 W 4.8e-3 THM 5.3e-6 QV 2.9e-4
    # NMP1 = P1 with Noah-MP v5.1.2 (sf_surface_physics = 4) against WRF v4.8.0 real*8 and single precision.
    ("real_nmp1_r8_fp64_cpu", "cpu_dp", "configs/real_nmp1_r8.ini", [],
     f"{WR}/run_NMP1_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_nmp1_r8_fp64_cpu_1h", "cpu_dp", "configs/real_nmp1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_NMP1_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_NMP1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_nmp1_r8_fp64_cuda_1h", "cuda_dp", "configs/real_nmp1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_NMP1_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_NMP1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_nmp1_fp32_1h", "cuda", "configs/real_nmp1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_NMP1", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_NMP1/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 3e-4, "W": 1e-2, "THM": 2e-6, "QV": 1e-4}),
    # Nesting: one-way (N1) and two-way (N2, feedback) vs WRF r8, bit-identical on the CPU (both domains;
    # N2 for 1 h)
    ("real_n1_r8_fp64_cpu", "cpu_dp", "configs/real_n1_steps.ini", ["real.dir=" + f"{WR}/run_N1_r8/steps"],
     {d: f"{WR}/run_N1_r8/steps/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 180,
     {"d01.U": 0.0, "d01.THM": 0.0, "d01.QV": 0.0, "d02.U": 0.0, "d02.THM": 0.0, "d02.QV": 0.0}),
    ("real_n1_r8_fp64_cuda", "cuda_dp", "configs/real_n1_steps.ini", ["real.dir=" + f"{WR}/run_N1_r8/steps"],
     {d: f"{WR}/run_N1_r8/steps/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 180,
     {}),
    # N1 for 1 h in FP32 against WRF's single-precision run
    ("real_n1_fp32_1h", "cuda", "configs/real_n1.ini", [],
     {d: f"{WR}/run_N1/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {f"{d}.{x}": t for d in ("d01", "d02") for x, t in (("U", 9e-4), ("W", 1.5e-2), ("THM", 3e-6), ("QV", 1e-4))}),
    ("real_n2_r8_fp64_cpu", "cpu_dp", "configs/real_n2_steps.ini", ["real.dir=" + f"{WR}/run_N2_r8/steps"],
     {d: f"{WR}/run_N2_r8/steps/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 180,
     {"d01.U": 0.0, "d01.THM": 0.0, "d01.QV": 0.0, "d02.U": 0.0, "d02.THM": 0.0, "d02.QV": 0.0}),
    ("real_n2_r8_fp64_cpu_1h", "cpu_dp", "configs/real_n2.ini", ["real.dir=" + f"{WR}/run_N2_r8"],
     {d: f"{WR}/run_N2_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {"d01.U": 0.0, "d01.W": 0.0, "d01.THM": 0.0, "d02.U": 0.0, "d02.W": 0.0, "d02.THM": 0.0, "d02.QV": 0.0}),
    # N2S1 = N2 with smooth_option = 1 (sm121) of the fed-back parent fields
    ("real_n2s1_r8_fp64_cpu", "cpu_dp", "configs/real_n2_steps.ini",
     ["real.dir=" + f"{WR}/run_N2S1_r8/steps", "domains.smooth_option=1"],
     {d: f"{WR}/run_N2S1_r8/steps/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 180,
     {"d01.U": 0.0, "d01.THM": 0.0, "d01.QV": 0.0, "d02.U": 0.0, "d02.THM": 0.0, "d02.QV": 0.0}),
    ("real_n2s1_r8_fp64_cpu_1h", "cpu_dp", "configs/real_n2.ini",
     ["real.dir=" + f"{WR}/run_N2S1_r8", "domains.smooth_option=1"],
     {d: f"{WR}/run_N2S1_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {"d01.U": 0.0, "d01.W": 0.0, "d01.THM": 0.0, "d02.U": 0.0, "d02.W": 0.0, "d02.THM": 0.0, "d02.QV": 0.0}),
    ("real_n2s1_r8_fp64_cuda_1h", "cuda_dp", "configs/real_n2.ini",
     ["real.dir=" + f"{WR}/run_N2S1_r8", "domains.smooth_option=1"],
     {d: f"{WR}/run_N2S1_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {}),
    # N2S2 = N2 with smooth_option = 2 (smdsm, WRF default) of the fed-back parent fields
    ("real_n2s2_r8_fp64_cpu", "cpu_dp", "configs/real_n2_steps.ini",
     ["real.dir=" + f"{WR}/run_N2S2_r8/steps", "domains.smooth_option=2"],
     {d: f"{WR}/run_N2S2_r8/steps/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 180,
     {"d01.U": 0.0, "d01.THM": 0.0, "d01.QV": 0.0, "d02.U": 0.0, "d02.THM": 0.0, "d02.QV": 0.0}),
    ("real_n2s2_r8_fp64_cpu_1h", "cpu_dp", "configs/real_n2.ini",
     ["real.dir=" + f"{WR}/run_N2S2_r8", "domains.smooth_option=2"],
     {d: f"{WR}/run_N2S2_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {"d01.U": 0.0, "d01.W": 0.0, "d01.THM": 0.0, "d02.U": 0.0, "d02.W": 0.0, "d02.THM": 0.0, "d02.QV": 0.0}),
    ("real_n2s2_r8_fp64_cuda_1h", "cuda_dp", "configs/real_n2.ini",
     ["real.dir=" + f"{WR}/run_N2S2_r8", "domains.smooth_option=2"],
     {d: f"{WR}/run_N2S2_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {}),
    # N3 = two-way nest with the P4 physics on both domains (physics feedback, nest ozone from the parent)
    ("real_n3_r8_fp64_cpu", "cpu_dp", "configs/real_n3_steps.ini", [],
     {d: f"{WR}/run_N3s_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 180,
     {f"{d}.{x}": 0.0 for d in ("d01", "d02") for x in ("U", "W", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP")}),
    ("real_n3_r8_fp64_cpu_1h", "cpu_dp", "configs/real_n3.ini", [],
     {d: f"{WR}/run_N3h_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {f"{d}.{x}": 0.0 for d in ("d01", "d02") for x in ("U", "W", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP")}),
    # Long runs (24 h and 72 h): P4, P6, P9 and N3 from real.exe inputs (wrfbdy every 6 h; 24-72 h: GFS through
    # tools/gcp/gcp_wps.sh), history every 3 h, restarted in the middle from CPPWRF's own wrfrst files
    # (suite.restart_at); every 3-hourly frame and the wrfrst files at the restart equal WRF's uninterrupted run
    # (tools/gcp/gcp_wrfref.sh; references: long_data.tar, 72 h long72_data.tar, long72n3_data.tar). All are slow
    # (suite.slow: wrf_suite.py --slow, gcp_run.sh --slow; the 24 h ones made the VM suite 65 instead of 43 min)
    *[(f"real_{c.lower()}_r8_fp64_cpu_{h}h", "cpu_dp", f"configs/real_{c.lower()}_r8.ini",
       [f"run.run_time={h * 3600}", "real.dir=" + f"{WR}/run_{c}_{h}h_r8", "output.out_interval=180",
        "output.diag_interval=180", f"suite.restart_at={h * 1800}"] + ["suite.slow"],
       f"{WR}/run_{c}_{h}h_r8/wrfout_d01_2026-10-12_00:00:00", h * 3600, {v: 0.0 for v in vs})
      for h in (24, 72)
      for c, vs in (("P4", ("U", "W", "THM", "QV", "QCLOUD", "QRAIN", "QICE")),
                    ("P6", ("U", "W", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QNICE", "QNRAIN")),
                    ("P9", ("U", "W", "THM", "QV")))],
    *[(f"real_n3_r8_fp64_cpu_{h}h", "cpu_dp", "configs/real_n3.ini",
       [f"run.run_time={h * 3600}", "real.dir=" + f"{WR}/run_N3_{h}h_r8", "phys.table_dir=" + f"{WR}/run_N3_{h}h_r8",
        "output.out_interval=180", "output.diag_interval=180", "d02.out_interval=540",
        f"suite.restart_at={h * 1800}"] + ["suite.slow"],
       {d: f"{WR}/run_N3_{h}h_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, h * 3600,
       {f"{d}.{x}": 0.0 for d in ("d01", "d02") for x in ("U", "W", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW",
                                                         "QGRAUP")})
      for h in (24, 72)],
    # N3 CUDA (needs > 6 GB; measured on L4 and G4): the P4 physics' roundoff growth, so P4's tolerances
    ("real_n3_r8_fp64_cuda_1h", "cuda_dp", "configs/real_n3.ini", [],
     {d: f"{WR}/run_N3h_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {}),
    ("real_n2_r8_fp64_cuda_1h", "cuda_dp", "configs/real_n2.ini", ["real.dir=" + f"{WR}/run_N2_r8"],
     {d: f"{WR}/run_N2_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {"d01.U": 0.0, "d01.THM": 0.0, "d02.U": 0.0, "d02.THM": 0.0}),
    ("real_n2_fp32_1h", "cuda", "configs/real_n2.ini", [],
     {d: f"{WR}/run_N2/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {"d01.U": 5e-4, "d01.THM": 3e-6, "d02.U": 6e-4, "d02.THM": 3e-6}),
    # WRF ideal cases (port of ideal.exe, run.case = em_*) vs WRF r8 ideal runs: bit-identical on the CPU
    ("ideal_hill2d_x_fp64_cpu", "cpu_dp", "configs/wrf_hill2d_x.ini", [],
     f"{WREF}/r8/run_hill_p0/wrfout_d01_0001-01-01_00:00:00", 3600,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "MU": 0.0, "P": 0.0, "THM": 0.0}),
    ("ideal_grav2d_x_fp64_cpu", "cpu_dp", "configs/wrf_grav2d_x.ini", [],
     f"{WREF}/r8/run_straka_d1/wrfout_d01_0001-01-01_00:00:00", 900,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "MU": 0.0, "P": 0.0, "THM": 0.0}),
    ("ideal_quarter_ss_fp64_cpu", "cpu_dp", "configs/wrf_quarter_ss.ini",
     ["run.run_time=7200", "output.out_interval=50", "output.diag_interval=50"],
     f"{WREF}/r8/run_qss/wrfout_d01_0001-01-01_00:00:00", 7200,
     {"U": 0.0, "V": 0.0, "W": 0.0, "PH": 0.0, "MU": 0.0, "P": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0,
      "QRAIN": 0.0, "RAINNC": 0.0}),
    # em_quarter_ss with WRF's default open_xs/xe/ys/ye (3-D open lateral boundaries, corners included)
    ("ideal_quarter_ss_open_fp64_cpu", "cpu_dp", "configs/wrf_quarter_ss.ini",
     ["domain.bc_x=open", "domain.bc_y=open", "run.run_time=7200", "output.out_interval=10",
      "output.diag_interval=10"],
     f"{WREF}/r8/run_qss_open/wrfout_d01_0001-01-01_00:00:00", 7200,
     {"U": 0.0, "V": 0.0, "W": 0.0, "PH": 0.0, "MU": 0.0, "P": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0,
      "QRAIN": 0.0, "RAINNC": 0.0}),
    # quarter_ss with diff_opt = 2 (periodic x and y): km_opt = 2, and km_opt = 2 with isfflx = 0 and a
    # constant drag coefficient and heat flux (tke_drag_coefficient, tke_heat_flux)
    ("ideal_quarter_ss_k2_fp64_cpu", "cpu_dp", "configs/wrf_quarter_ss.ini",
     ["dyn.diff_opt=2", "dyn.km_opt=2", "run.run_time=600", "output.out_interval=10", "output.diag_interval=10"],
     f"{WREF}/r8/run_qss_k2/wrfout_d01_0001-01-01_00:00:00", 600,
     {"U": 0.0, "V": 0.0, "W": 0.0, "PH": 0.0, "MU": 0.0, "P": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0,
      "QRAIN": 0.0, "RAINNC": 0.0}),
    ("ideal_quarter_ss_k2_isf0_fp64_cpu", "cpu_dp", "configs/wrf_quarter_ss.ini",
     ["dyn.diff_opt=2", "dyn.km_opt=2", "phys.isfflx=0", "dyn.tke_heat_flux=0.05", "dyn.tke_drag_coefficient=0.0015",
      "run.run_time=600", "output.out_interval=10", "output.diag_interval=10"],
     f"{WREF}/r8/run_qss_k2_isf0/wrfout_d01_0001-01-01_00:00:00", 600,
     {"U": 0.0, "V": 0.0, "W": 0.0, "PH": 0.0, "MU": 0.0, "P": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0,
      "QRAIN": 0.0, "RAINNC": 0.0}),
    ("ideal_squall2d_x_fp64_cpu", "cpu_dp", "configs/wrf_squall2d_x.ini", ["output.out_interval=20"],
     f"{WREF}/r8/run_squall2d_kessler_10min/wrfout_d01_0001-01-01_00:00:00", 600,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "MU": 0.0, "P": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0}),
    ("ideal_squall2d_y_fp64_cpu", "cpu_dp", "configs/wrf_squall2d_y.ini", ["output.out_interval=20"],
     f"{WREF}/r8/run_squall2d_ys_10min/wrfout_d01_0001-01-01_00:00:00", 600,
     {"U": 0.0, "V": 0.0, "W": 0.0, "PH": 0.0, "MU": 0.0, "P": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0,
      "QRAIN": 0.0}),
    # symmetric lateral boundaries (WRF symmetric_xs/xe/ys/ye, rigid walls) with soundings without wind:
    # quarter_ss on all four sides, the Straka density current and squall2d_x (TKE) in x, long enough for the
    # storm / cold pools to reach the walls
    ("ideal_quarter_ss_sym_fp64_cpu", "cpu_dp", "configs/wrf_quarter_ss_sym.ini",
     ["run.run_time=3600", "output.out_interval=50", "output.diag_interval=50"],
     f"{WREF}/r8/run_qss_sym/wrfout_d01_0001-01-01_00:00:00", 3600,
     {"U": 0.0, "V": 0.0, "W": 0.0, "PH": 0.0, "MU": 0.0, "P": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0,
      "QRAIN": 0.0, "RAINNC": 0.0}),
    ("ideal_grav2d_x_sym_fp64_cpu", "cpu_dp", "configs/wrf_grav2d_x_sym.ini", ["run.run_time=2400", "output.out_interval=600"],
     f"{WREF}/r8/run_grav2d_x_sym/wrfout_d01_0001-01-01_00:00:00", 2400,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "MU": 0.0, "P": 0.0, "THM": 0.0}),
    ("ideal_squall2d_x_sym_fp64_cpu", "cpu_dp", "configs/wrf_squall2d_x_sym.ini",
     ["run.run_time=2400", "output.out_interval=200"],
     f"{WREF}/r8/run_squall2d_x_sym/wrfout_d01_0001-01-01_00:00:00", 2400,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "MU": 0.0, "P": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0}),
    ("ideal_squall2d_x_fp64_cuda", "cuda_dp", "configs/wrf_squall2d_x.ini", ["output.out_interval=20"],
     f"{WREF}/r8/run_squall2d_kessler_10min/wrfout_d01_0001-01-01_00:00:00", 600,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # FP32 against WRF's single-precision run: WRF r8 itself differs from it by U 8e-5, THM 5e-7
    ("ideal_squall2d_x_fp32", "cuda", "configs/wrf_squall2d_x.ini", ["output.out_interval=200"],
     f"{WREF}/run_squall2d_kessler/wrfout_d01_0001-01-01_00:00:00", 600,
     {"U": 3e-4, "W": 1e-3, "THM": 2e-6, "QV": 1e-4}),
    ("ideal_quarter_ss_fp64_cuda", "cuda_dp", "configs/wrf_quarter_ss.ini",
     ["run.run_time=7200", "output.out_interval=50", "output.diag_interval=50"],
     f"{WREF}/r8/run_qss/wrfout_d01_0001-01-01_00:00:00", 7200,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    # the other ideal cases in CUDA FP64 (the CPU cases' configurations and references)
    ("ideal_hill2d_x_fp64_cuda", "cuda_dp", "configs/wrf_hill2d_x.ini", [],
     f"{WREF}/r8/run_hill_p0/wrfout_d01_0001-01-01_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0}),
    ("ideal_grav2d_x_fp64_cuda", "cuda_dp", "configs/wrf_grav2d_x.ini", [],
     f"{WREF}/r8/run_straka_d1/wrfout_d01_0001-01-01_00:00:00", 900, {"U": 0.0, "W": 0.0, "THM": 0.0}),
    ("ideal_quarter_ss_open_fp64_cuda", "cuda_dp", "configs/wrf_quarter_ss.ini",
     ["domain.bc_x=open", "domain.bc_y=open", "run.run_time=7200", "output.out_interval=10",
      "output.diag_interval=10"],
     f"{WREF}/r8/run_qss_open/wrfout_d01_0001-01-01_00:00:00", 7200,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("ideal_quarter_ss_k2_fp64_cuda", "cuda_dp", "configs/wrf_quarter_ss.ini",
     ["dyn.diff_opt=2", "dyn.km_opt=2", "run.run_time=600", "output.out_interval=10", "output.diag_interval=10"],
     f"{WREF}/r8/run_qss_k2/wrfout_d01_0001-01-01_00:00:00", 600,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("ideal_quarter_ss_k2_isf0_fp64_cuda", "cuda_dp", "configs/wrf_quarter_ss.ini",
     ["dyn.diff_opt=2", "dyn.km_opt=2", "phys.isfflx=0", "dyn.tke_heat_flux=0.05", "dyn.tke_drag_coefficient=0.0015",
      "run.run_time=600", "output.out_interval=10", "output.diag_interval=10"],
     f"{WREF}/r8/run_qss_k2_isf0/wrfout_d01_0001-01-01_00:00:00", 600,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("ideal_squall2d_y_fp64_cuda", "cuda_dp", "configs/wrf_squall2d_y.ini", ["output.out_interval=20"],
     f"{WREF}/r8/run_squall2d_ys_10min/wrfout_d01_0001-01-01_00:00:00", 600,
     {"V": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("ideal_quarter_ss_sym_fp64_cuda", "cuda_dp", "configs/wrf_quarter_ss_sym.ini",
     ["run.run_time=3600", "output.out_interval=50", "output.diag_interval=50"],
     f"{WREF}/r8/run_qss_sym/wrfout_d01_0001-01-01_00:00:00", 3600,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("ideal_grav2d_x_sym_fp64_cuda", "cuda_dp", "configs/wrf_grav2d_x_sym.ini",
     ["run.run_time=2400", "output.out_interval=600"],
     f"{WREF}/r8/run_grav2d_x_sym/wrfout_d01_0001-01-01_00:00:00", 2400, {"U": 0.0, "W": 0.0, "THM": 0.0}),
    ("ideal_squall2d_x_sym_fp64_cuda", "cuda_dp", "configs/wrf_squall2d_x_sym.ini",
     ["run.run_time=2400", "output.out_interval=200"],
     f"{WREF}/r8/run_squall2d_x_sym/wrfout_d01_0001-01-01_00:00:00", 2400,
     {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("wdamp_w1_fp64", "cpu_dp", "configs/real_r3.ini",
     ["real.dir=" + f"{WR}/run_W1_r8", "time.dt=450", "run.run_time=4500", "output.out_interval=1",
      "output.diag_interval=1"], f"{WR}/run_W1_r8/wrfout_d01_2026-10-12_00:00:00", 4500,
     {"U": 0.0, "W": 0.0, "THM": 0.0}),
    ("wdamp_w1_fp64_cuda", "cuda_dp", "configs/real_r3.ini",
     ["real.dir=" + f"{WR}/run_W1_r8", "time.dt=450", "run.run_time=4500", "output.out_interval=1",
      "output.diag_interval=1"], f"{WR}/run_W1_r8/wrfout_d01_2026-10-12_00:00:00", 4500,
      {"U": 0.0, "W": 0.0, "THM": 0.0}),
]


# restart_<tag>_fp64: 3 steps straight vs restarted after step 1 (restart_p1)
RESTARTS = (("p1", "configs/real_p1_r8.ini", "run_P1_r8"), ("a01", "configs/real_a01_r8.ini", "run_P1_r8"),
            ("p7", "configs/real_p7_r8.ini", "run_P7_r8"),
            ("p8", "configs/real_p8_r8.ini", "run_P8_r8"), ("p9", "configs/real_p9_r8.ini", "run_P9_r8"),
            ("nmp1", "configs/real_nmp1_r8.ini", "run_NMP1_r8"), ("mynn", "configs/real_mynn_r8.ini", "run_MYNN_r8"),
            ("mynnsf1", "configs/real_mynnsf1_r8.ini", "run_P6_r8"), ("mynnsf2", "configs/real_mynnsf2_r8.ini", "run_P6_r8"),
            ("mynnx1", "configs/real_mynnx1_r8.ini", "run_P6_r8"), ("mynnx2", "configs/real_mynnx2_r8.ini", "run_P6_r8"),
            ("t1", "configs/real_t1.ini", "run_R1_r8"),
            ("a05", "configs/real_a05_r8.ini", "run_A05_r8"), ("a05g", "configs/real_a05g_r8.ini", "run_A05G_r8"),
            ("slab1", "configs/real_slab1_r8.ini", "run_SLAB1_r8"))
# restart_<tag>_fp64_cuda: the same in CUDA FP64 (the GPU run is deterministic: restarted = straight bit for bit)
SPECIAL = (["straka_wrf", "straka_wrf_cuda", "restart_p1_fp32_cuda"] + [f"restart_{t}_fp64" for t, _, _ in RESTARTS]
           + [f"restart_{t}_fp64_cuda" for t, _, _ in RESTARTS])


def run(binary, cfg, sets, out):
    shutil.rmtree(out, ignore_errors=True)
    cmd = [binary, os.path.join(REPO, cfg), "--set", f"output.dir={out}"]
    for s in sets:
        cmd += ["--set", s]
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "3"))
    r = subprocess.run(cmd, cwd=REPO, env=env, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def config_vs_namelist(cfg, sets, nml):
    """Options of the case's configuration (cfg and its --set) that differ from tools/wrf2cppwrf.py's translation of
    the namelist.input of the WRF run that made the reference (options both set; paths by their real path; the run
    length, output and restart settings aside), or the translation's errors. A case with output.auxhist3=true (the
    translation of io_form_auxhist3 = 2, which needs WRF history output) is translated with --format wrfout."""
    sys.path.insert(0, REPO)
    from tools.wrf2cppwrf import parse_namelist, WrfTranslator
    tr = WrfTranslator(parse_namelist(open(nml).read()), namelist_path=nml,
                       output_format="wrfout" if "output.auxhist3=true" in sets else "native")
    gen = tr.translate()
    if tr.errors:
        return [f"{nml} does not translate: {tr.errors}"]
    conf, sec = {}, None
    for line in open(os.path.join(REPO, cfg)):
        line = line.split("#")[0].strip()
        if re.fullmatch(r"\[.+\]", line):
            sec = line[1:-1]
        elif sec and "=" in line:
            conf.setdefault(sec, {})[line.split("=", 1)[0].strip()] = line.split("=", 1)[1].strip()
    for x in sets:
        k, v = x.split("=", 1)
        conf.setdefault(k.split(".", 1)[0], {})[k.split(".", 1)[1]] = v

    def norm(v):
        v = os.path.expandvars(str(v)).strip().strip("'\"").replace(" ", "")
        if v.startswith("/"):
            return os.path.realpath(v)
        try:
            return repr(float(v.replace("d", "e").replace("D", "e")))
        except ValueError:
            return {".true.": "true", ".false.": "false", "t": "true", "f": "false"}.get(v.lower(), v.lower())
    # a reference directory without its own wrfinput (tools/gcp/gcp_wrfref.sh links the base run's) does not name the
    # input directory
    skip = {("run", "run_time")} | (set() if os.path.exists(os.path.join(os.path.dirname(nml), "wrfinput_d01"))
                                    else {("real", "dir")})
    def same(a, b):  # files by their content (a reference directory's copy of a linked input_sounding)
        a, b = norm(a), norm(b)
        return a == b or (a.startswith("/") and os.path.isfile(a) and os.path.isfile(b) and filecmp.cmp(a, b, False))
    return [f"{s}.{k}: namelist {v} configuration {conf[s][k]}" for s, kv in gen.items() if s not in ("output", "restart")
            for k, v in kv.items() if (s, k) not in skip and k in conf.get(s, {}) and not same(v, conf[s][k])]


def restart_at(sets):
    """The time [s] of a case's restart (suite.restart_at=S in its --set list: the run stops there, writes wrfrst
    files and resumes from them), or None."""
    return next((float(x.split("=", 1)[1]) for x in sets if x.startswith("suite.restart_at=")), None)


def slow(sets):
    """A slow case (suite.slow in its --set list): run only with --slow or when --only --exact names it."""
    return "suite.slow" in sets


def exact_case(refs, out, t, sets):
    """Tolerance-0 comparison of the wrfout files in `out` (out/<domain> for nests; a list of directories: their
    files together, a later directory's frame replacing an earlier one's) with the references {domain: wrfout}:
    (differences, summary)."""
    import netCDF4
    import exact
    bad, frames, nvars = [], 0, 0
    for dom, w in refs.items():
        d = dom or re.search(r"wrfout_(d\d\d)_", os.path.basename(w)).group(1)
        with netCDF4.Dataset(w) as ds:
            xt = [float(x) * 60.0 for x in ds["XTIME"][:]]
            final = [exact.text(ds["Times"], i).decode() for i, x in enumerate(xt) if abs(x - t) < 1e-6]
            nvars += sum(1 for v in ds.variables.values() if v.dimensions[:1] == ("Time",) and v.dtype.kind != "S")
        if not final:
            bad.append(f"{d}: the reference has no frame at t={t:g}s")
            continue
        mine = [f for o in (out if isinstance(out, list) else [out])
                for f in glob.glob(os.path.join(o, f"wrfout_{d}_*")) + glob.glob(os.path.join(o, d, f"wrfout_{d}_*"))]
        dif, n = exact.compare_stream(mine, [w], subset=True, until=final[0],
                                      global_atts=os.path.exists(os.path.join(os.path.dirname(w), "namelist.input")))
        with netCDF4.Dataset(w) as ds:  # WRF leaves RQCBLTEN / RQIBLTEN uninitialised without qc / qi
            for v, q in (("RQCBLTEN", "QCLOUD"), ("RQIBLTEN", "QICE")):
                if q not in ds.variables:
                    dif.pop(v, None)
        frames += n
        bad += [f"{d}.{k}: {v}" for k, v in dif.items()]
        # output_diagnostics: the auxhist3 files next to the reference wrfout, compared as the wrfout files
        aux = sorted(glob.glob(os.path.join(os.path.dirname(w), f"auxhist3_{d}_*")))
        if aux:
            mine = [f for o in (out if isinstance(out, list) else [out])
                    for f in glob.glob(os.path.join(o, f"auxhist3_{d}_*")) + glob.glob(os.path.join(o, d, f"auxhist3_{d}_*"))]
            dif, n = exact.compare_stream(mine, aux, subset=True, until=final[0],
                                          global_atts=os.path.exists(os.path.join(os.path.dirname(w), "namelist.input")))
            if not n:
                dif.setdefault("(times)", "no auxhist3 frame in common")
            frames += n
            bad += [f"{d}.auxhist3.{k}: {v}" for k, v in dif.items()]
    return bad, f"{frames} frames, {nvars} variables"


def restarted_case(b, cfg, sets, refs, t, mid, out):
    """A run to `t` s restarted at `mid` s from its own wrfrst files: the first run (to one step past mid, wrfrst files
    at mid only: written within the run, as WRF's, not at its last step, where WRF's history alarm is still ringing)
    and the resumed one (mid to t) together against WRF's uninterrupted run at tolerance 0 (every reference frame),
    and the wrfrst files at mid against WRF's (the reference run's restart_interval writes them there):
    (differences, summary). As run_wrf_files_test.sh, RQCBLTEN / RQIBLTEN are not compared when WRF's restart file
    has no QCLOUD / QICE (WRF leaves them uninitialised)."""
    import netCDF4
    import exact
    first, resumed = os.path.join(out, "a_first"), os.path.join(out, "b_resumed")  # sorted: resumed frames win
    w0 = next(iter(refs.values()))
    with netCDF4.Dataset(w0) as ds:
        at = [exact.text(ds["Times"], i).decode() for i, x in enumerate(ds["XTIME"][:]) if abs(float(x) * 60.0 - mid) < 1e-6]
    if not at:
        return [f"{w0} has no frame at the restart time {mid:g}s"], ""
    nl = os.path.join(os.path.dirname(w0), "namelist.input")
    dt = float(re.search(r"^\s*time_step\s*=\s*(\d+)", open(nl).read(), re.M).group(1))
    rc, log = run(b, cfg, sets + [f"run.run_time={mid + dt:g}", f"restart.interval={mid:g}", f"restart.end={mid:g}"],
                  first)
    if rc != 0:
        return [f"run to the restart failed: {log[-200:]}"], ""
    # a nested run writes each domain's files in first/<domain>; the resumed run reads every domain's wrfrst file
    # from the directory of d01's (as WRF's run directory)
    rst = os.path.join(out, "rst")
    os.makedirs(rst, exist_ok=True)
    for pat in (f"wrfrst_d??_{at[0]}", f"d??/wrfrst_d??_{at[0]}"):
        for f in glob.glob(os.path.join(first, pat)):
            os.symlink(f, os.path.join(rst, os.path.basename(f)))
    rc, log = run(b, cfg, sets + [f"run.restart_file={rst}/wrfrst_d01_{at[0]}"], resumed)
    if rc != 0:
        return [f"resumed run failed: {log[-200:]}"], ""
    bad, what = exact_case(refs, [first, resumed], t, sets)
    for dom, w in refs.items():
        d = dom or re.search(r"wrfout_(d\d\d)_", os.path.basename(w)).group(1)
        wr = os.path.join(os.path.dirname(w), f"wrfrst_{d}_{at[0]}")
        mine = glob.glob(os.path.join(first, f"wrfrst_{d}_{at[0]}")) + glob.glob(os.path.join(first, d, f"wrfrst_{d}_{at[0]}"))
        if not os.path.exists(wr) or not mine:
            bad.append(f"{d}: no wrfrst_{d}_{at[0]} from {'WRF' if not os.path.exists(wr) else 'CPPWRF'}")
            continue
        with netCDF4.Dataset(wr) as ds:
            noise = [x for x, q in (("RQCBLTEN", "QCLOUD"), ("RQIBLTEN", "QICE")) if q not in ds.variables]
        dif, _ = exact.compare_stream(mine, [wr])
        bad += [f"{d}.wrfrst.{k}: {v}" for k, v in dif.items() if k not in noise]
    return bad, f"{what}, restarted at {mid:g}s (wrfrst identical)" if not bad else what


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


def restart_p1(binary, out, cfg="configs/real_p1_r8.ini", run_dir="run_P1_r8"):
    """Physics/TKE restart: 3 steps straight vs restarted after step 1 must be bit-identical
    (dycore, TKE, surface/soil state, R*BLTEN, lagged t_1; for P7 also R*CUTEN and RTHFTEN/RQVFTEN)."""
    if not os.path.exists(f"{WR}/{run_dir}/steps/wrfinput_d01") and not os.path.exists(f"{WR}/{run_dir}/wrfinput_d01"):
        return None, "input missing"
    a, b = os.path.join(out, "straight"), os.path.join(out, "restarted")
    rc, log = run(binary, cfg, ["restart.interval=60"], a)
    if rc != 0:
        return False, log[-300:]
    rc, log = run(binary, cfg, [f"run.restart_file={a}/restart_000001.nc"], b)
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
    ap.add_argument("--only", help="comma-separated substrings of case names")
    ap.add_argument("--exact", action="store_true", help="--only names are exact case names")
    ap.add_argument("--jobs", type=int, default=1, help="CPU cases run at the same time (GPU cases: one at a time)")
    ap.add_argument("--keep", action="store_true", help="keep the outputs of passing cases (default: deleted)")
    ap.add_argument("--allow-skip", action="store_true", help="cases without their reference do not fail the suite")
    ap.add_argument("--slow", action="store_true", help="also the slow cases (suite.slow: the 24 h and 72 h runs)")
    a = ap.parse_args()
    only = a.only.split(",") if a.only else None

    def wanted(name):
        return not only or (name in only if a.exact else any(o in name for o in only))

    unknown = [o for o in only or [] if not any(o == n if a.exact else o in n for n in [c[0] for c in CASES] + SPECIAL)]
    if unknown:  # a misspelt case would otherwise make an empty suite that passes
        print(f"ERROR no case named {', '.join(unknown)}", flush=True)
        sys.exit(2)

    def done(path, good):  # the outputs of passing cases are deleted (they take tens of GB)
        if good and not a.keep:
            shutil.rmtree(path, ignore_errors=True)

    def case(name, b, cfg, sets, refs, t, tol):
        out = os.path.join(a.out, name)
        mid = restart_at(sets)
        sets = [x for x in sets if not x.startswith("suite.")]
        exact = all(x == 0.0 for x in tol.values())
        if exact:  # WRF history output; with the reference's namelist.input CPPWRF writes WRF's global attributes
            nl = os.path.join(os.path.dirname(next(iter(refs.values()))), "namelist.input")
            if os.path.exists(nl):  # the configuration is the WRF run's that made the reference
                dif = config_vs_namelist(cfg, sets, nl)
                if dif:
                    return False, f"FAIL  {name:22s} configuration differs from {nl}: {'; '.join(dif)}"
            sets = sets + ["output.format=wrfout"] + ([f"output.wrf_namelist={nl}"] if os.path.exists(nl) else [])
        if exact and mid:
            bad, what = restarted_case(b, cfg, sets, refs, t, mid, out)
            done(out, not bad)
            return not bad, (f"{'PASS' if not bad else 'FAIL'}  {name:22s} t={t:g}s exact: {what}"
                             + (f"  [{'; '.join(x[:300] for x in bad[:6])}{' ...' if len(bad) > 6 else ''}]" if bad else ""))
        rc, log = run(b, cfg, sets, out)
        if rc != 0:
            return False, f"FAIL  {name:22s} run failed: {log[-200:]}"
        if exact:
            bad, what = exact_case(refs, out, t, sets)
            done(out, not bad)
            return not bad, (f"{'PASS' if not bad else 'FAIL'}  {name:22s} t={t:g}s exact: {what}"
                             + (f"  [{'; '.join(x[:300] for x in bad[:6])}{' ...' if len(bad) > 6 else ''}]" if bad else ""))
        res = {}
        for dom, w in refs.items():  # nested: out/<domain>, tolerance keys "<domain>.<var>"
            for v, x in compare(w, os.path.join(out, dom) if dom else out, t).items():
                res[f"{dom}.{v}" if dom else v] = x
        bad = [f"{v} {res.get(v, float('nan')):.2e}>{tl:.0e}" for v, tl in tol.items()
               if not (res.get(v, float("inf")) <= tl)]
        txt = " ".join(f"{v}={res.get(v, float('nan')):.1e}" for v in tol)
        done(out, not bad)
        return not bad, f"{'PASS' if not bad else 'FAIL'}  {name:22s} t={t:g}s {txt}" + (f"  [{'; '.join(bad)}]" if bad else "")

    def special(name, fn, path):
        good, msg = fn()
        if good is None:
            return None, f"SKIP  {name:22s} ({msg})"
        done(path, good)
        return good, f"{'PASS' if good else 'FAIL'}  {name:22s} {msg}"

    bins = {"cpu_dp": a.cpu_dp, "cuda": a.cuda, "cuda_dp": a.cuda_dp}
    for k, b in bins.items():  # a binary that was asked for and is not there (a failed build) is no SKIP
        if b and not os.path.exists(b):
            print(f"ERROR --{k.replace('_', '-')} {b}: no such file (build failed?)", flush=True)
            sys.exit(2)
    cpu, gpu = [], []
    count = {"PASS": 0, "FAIL": 0, "SKIP": 0}
    for name, key, cfg, sets, wo, t, tol in CASES:
        if not wanted(name) or not bins.get(key):  # cases of binaries not given are not selected
            continue
        if slow(sets) and not a.slow and not (a.exact and name in only):  # slow cases only on request
            continue
        b = bins[key]
        refs = wo if isinstance(wo, dict) else {"": wo}  # nested runs: {domain: wrfout}
        missing = [w for w in refs.values() if not os.path.exists(w)]
        if missing:
            print(f"SKIP  {name:22s} (reference missing: {', '.join(missing)})", flush=True)
            count["SKIP"] += 1
            continue
        (cpu if key == "cpu_dp" else gpu).append((t, name, lambda n=name, b=b, c=cfg, s=sets, r=refs, t=t, tl=tol: case(n, b, c, s, r, t, tl)))
    if a.cpu_dp and os.path.exists(a.cpu_dp):
        if not only or wanted("straka_wrf"):
            p = os.path.join(a.out, "straka")
            cpu.append((0, "straka_wrf", lambda p=p: special("straka_wrf", lambda: straka(a.cpu_dp, p), p)))
        for tag, cfg, rd in RESTARTS:
            rname = f"restart_{tag}_fp64"
            if not wanted(rname):
                continue
            p = os.path.join(a.out, f"restart_{tag}")
            cpu.append((0, rname, lambda r=rname, c=cfg, d=rd, p=p: special(r, lambda: restart_p1(a.cpu_dp, p, c, d), p)))
    if a.cuda and os.path.exists(a.cuda) and wanted("restart_p1_fp32_cuda"):  # single precision, GPU
        p = os.path.join(a.out, "restart_p1_fp32_cuda")
        gpu.append((0, "restart_p1_fp32_cuda", lambda p=p: special("restart_p1_fp32_cuda", lambda: restart_p1(a.cuda, p), p)))
    if a.cuda_dp and os.path.exists(a.cuda_dp):
        if wanted("straka_wrf_cuda"):
            p = os.path.join(a.out, "straka_cuda")
            gpu.append((0, "straka_wrf_cuda", lambda p=p: special("straka_wrf_cuda", lambda: straka(a.cuda_dp, p), p)))
        for tag, cfg, rd in RESTARTS:
            rname = f"restart_{tag}_fp64_cuda"
            if not wanted(rname):
                continue
            p = os.path.join(a.out, f"restart_{tag}_cuda")
            gpu.append((0, rname, lambda r=rname, c=cfg, d=rd, p=p: special(r, lambda: restart_p1(a.cuda_dp, p, c, d), p)))
    # longest simulated time first; each task prints its line as soon as it ends; GPU cases share the one GPU
    cpu = sorted(cpu, key=lambda x: -x[0])
    from concurrent.futures import ThreadPoolExecutor, as_completed

    def timed(fn):  # the case's line ends with its wall time
        t0 = time.monotonic()
        good, line = fn()
        return good, f"{line}  ({time.monotonic() - t0:.0f} s)"
    with ThreadPoolExecutor(max(1, a.jobs)) as pc, ThreadPoolExecutor(1) as pg:
        futs = {pc.submit(timed, f): n for _, n, f in cpu} | {pg.submit(timed, f): n for _, n, f in gpu}
        for f in as_completed(futs):
            try:
                good, line = f.result()
            except Exception as e:  # noqa: BLE001 - a case the harness cannot finish fails, the others go on
                good, line = False, f"FAIL  {futs[f]:22s} harness error: {type(e).__name__}: {e}"
            print(line, flush=True)
            count["SKIP" if good is None else "PASS" if good else "FAIL"] += 1
    print(f"SUITE PASS {count['PASS']} FAIL {count['FAIL']} SKIP {count['SKIP']}", flush=True)
    sys.exit(0 if not count["FAIL"] and (a.allow_skip or not count["SKIP"]) else 1)


if __name__ == "__main__":
    main()
