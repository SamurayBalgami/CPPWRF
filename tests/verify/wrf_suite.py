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
# WRF reference runs (wrfinput/wrfbdy and wrfout of each case); the configs read them as ${CPPWRF_WRF_RUNS}
WR = os.environ.setdefault("CPPWRF_WRF_RUNS", os.path.expanduser("~/wrf_real480"))
os.environ.setdefault("CPPWRF_WRF_REF", os.path.expanduser("~/wrf_ref480"))
WREF = os.environ["CPPWRF_WRF_REF"]

# name, binary key, config, extra --set, WRF wrfout, time [s], {var: max rel L2 (whole domain)}
CASES = [
    # CPU FP64 real-data dycore: bit-identical to WRF (every field, every step; libmvec lanes, WRF's
    # summation orders and halo schedule reproduced)
    ("real_r1_fp64_cpu", "cpu_dp", "configs/real_r1.ini", [], f"{WR}/run_R1_r8/wrfout_d01_2026-10-12_00:00:00", 300,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_r1_fp64_cuda", "cuda_dp", "configs/real_r1.ini", [], f"{WR}/run_R1_r8/wrfout_d01_2026-10-12_00:00:00", 300,
     {"U": 1e-12, "W": 1e-10, "PH": 1e-13, "THM": 1e-13, "QV": 1e-13, "MU": 1e-13}),
    # T1 = R1 + diff_opt = 2, km_opt = 2 (1.5-order TKE, full metric stress tensor diffusion):
    # bit-identical on CPU FP64 across all fields including TKE
    ("real_t1_fp64_cpu", "cpu_dp", "configs/real_t1.ini", [], f"{WR}/run_T1_r8/wrfout_d01_2026-10-12_00:00:00", 300,
     {"U": 0.0, "V": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0, "TKE": 0.0}),
    ("real_t1_fp64_cuda", "cuda_dp", "configs/real_t1.ini", [], f"{WR}/run_T1_r8/wrfout_d01_2026-10-12_00:00:00", 300,
     {"U": 1e-12, "V": 1e-12, "W": 1e-10, "PH": 1e-13, "THM": 1e-13, "QV": 1e-13, "MU": 1e-13, "TKE": 1e-10}),
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
     f"{WR}/run_P1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 3e-9, "W": 1e-7, "THM": 1e-9, "QV": 1e-8}),
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
    # P2 (+Kain-Fritsch cu_physics=1): bit-identical (3 steps here, 1 h above)
    ("real_p2_r8_fp64_cpu", "cpu_dp", "configs/real_p2_r8.ini", [], f"{WR}/run_P2_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_p2_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p2_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P2_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P2_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 2e-4, "W": 5e-3, "THM": 3e-6, "QV": 2e-4}),
    ("real_p2_fp32_1h", "cuda", "configs/real_p2_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P2", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P2/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 2e-3, "W": 8e-2, "THM": 3e-5, "QV": 8e-4}),
    # P3 = P1 + RRTMG LW/SW every step (radt 1, o3input 2, ghg_input 1); tables in ${CPPWRF_WRF_RUNS}/rrtmg_tables
    ("real_p3_r8_fp64_cpu", "cpu_dp", "configs/real_p3_r8.ini", [], f"{WR}/run_P3_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_p3_fp32_1h", "cuda", "configs/real_p3_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P3", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P3/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 3e-4, "W": 1e-2, "THM": 2e-6, "QV": 1e-4}),
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
     f"{WR}/run_P4_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 6e-4, "W": 1e-2, "THM": 8e-6, "QV": 5e-4}),
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
    # P7 = P2 with New Tiedtke (cu_physics = 16) instead of KF; reference WRF-4.8.0-r8d8 (r8 + double-8).
    # The scheme has thresholds that roundoff flips (one column at step 4 jumps 1.3e-4 K under 1e-15 input noise).
    ("real_p7_r8_fp64_cpu", "cpu_dp", "configs/real_p7_r8.ini", [],
     f"{WR}/run_P7_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_p7_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p7_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P7_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P7_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 1e-5, "W": 2e-4, "THM": 3e-7, "QV": 3e-5}),
    ("real_p7_fp32_1h", "cuda", "configs/real_p7_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P7", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P7/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 2e-3, "W": 5e-2, "THM": 3e-5, "QV": 8e-4}),
    # P8 = P2 with Grell-Freitas-Li (cu_physics = 3, GFLSCHEME, ishallow = 1); reference WRF-4.8.0-r8d8.
    # Column replay of cu_grell_freitas_li through WRF's objects is bit-identical (tests/phys/cu_gfl).
    ("real_p8_r8_fp64_cpu_10min", "cpu_dp", "configs/real_p8_r8.ini", ["run.run_time=600"],
     f"{WR}/run_P8_r8/steps/wrfout_d01_2026-10-12_00:00:00", 600,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_p8_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p8_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P8_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P8_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 1e-7, "W": 1e-6, "THM": 1e-9, "QV": 1e-7}),
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
     f"{WR}/run_P9_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 3e-9, "W": 1e-7, "THM": 1e-9, "QV": 1e-8}),
    ("real_p9_fp32_1h", "cuda", "configs/real_p9_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_P9", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_P9/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 6e-4, "W": 1e-2, "THM": 1e-5, "QV": 2e-3}),
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
    # NMP1 = P1 with Noah-MP v5.1.2 (sf_surface_physics = 4) against WRF v4.8.0 real*8 and single precision.
    ("real_nmp1_r8_fp64_cpu", "cpu_dp", "configs/real_nmp1_r8.ini", [],
     f"{WR}/run_NMP1_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_nmp1_r8_fp64_cpu_1h", "cpu_dp", "configs/real_nmp1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_NMP1_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_NMP1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_nmp1_r8_fp64_cuda_1h", "cuda_dp", "configs/real_nmp1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_NMP1_r8", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_NMP1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 3e-9, "W": 1e-7, "THM": 1e-9, "QV": 1e-8}),
    ("real_nmp1_fp32_1h", "cuda", "configs/real_nmp1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_NMP1", "output.out_interval=10", "output.diag_interval=10"],
     f"{WR}/run_NMP1/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 3e-4, "W": 1e-2, "THM": 2e-6, "QV": 1e-4}),
    # Nesting: one-way (N1) and two-way (N2, feedback) vs WRF r8, bit-identical on the CPU (both domains;
    # N2 for 1 h)
    ("real_n1_r8_fp64_cpu", "cpu_dp", "configs/real_n1_steps.ini", ["real.dir=" + f"{WR}/run_N1_r8/steps"],
     {d: f"{WR}/run_N1_r8/steps/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 180,
     {"d01.U": 0.0, "d01.THM": 0.0, "d01.QV": 0.0, "d02.U": 0.0, "d02.THM": 0.0, "d02.QV": 0.0}),
    ("real_n2_r8_fp64_cpu", "cpu_dp", "configs/real_n2_steps.ini", [],
     {d: f"{WR}/run_N2_r8/steps/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 180,
     {"d01.U": 0.0, "d01.THM": 0.0, "d01.QV": 0.0, "d02.U": 0.0, "d02.THM": 0.0, "d02.QV": 0.0}),
    ("real_n2_r8_fp64_cpu_1h", "cpu_dp", "configs/real_n2.ini", ["real.dir=" + f"{WR}/run_N2_r8"],
     {d: f"{WR}/run_N2_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {"d01.U": 0.0, "d01.W": 0.0, "d01.THM": 0.0, "d02.U": 0.0, "d02.W": 0.0, "d02.THM": 0.0, "d02.QV": 0.0}),
    ("real_n2_r8_fp64_cuda_1h", "cuda_dp", "configs/real_n2.ini", ["real.dir=" + f"{WR}/run_N2_r8"],
     {d: f"{WR}/run_N2_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {"d01.U": 1e-11, "d01.THM": 1e-13, "d02.U": 1e-11, "d02.THM": 1e-13}),
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
    ("ideal_squall2d_x_fp64_cpu", "cpu_dp", "configs/wrf_squall2d_x.ini", ["output.out_interval=20"],
     f"{WREF}/r8/run_squall2d_kessler_10min/wrfout_d01_0001-01-01_00:00:00", 600,
     {"U": 0.0, "W": 0.0, "PH": 0.0, "MU": 0.0, "P": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0}),
    ("ideal_squall2d_y_fp64_cpu", "cpu_dp", "configs/wrf_squall2d_y.ini", ["output.out_interval=20"],
     f"{WREF}/r8/run_squall2d_y_10min/wrfout_d01_0001-01-01_00:00:00", 600,
     {"U": 0.0, "V": 0.0, "W": 0.0, "PH": 0.0, "MU": 0.0, "P": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0,
      "QRAIN": 0.0}),
    ("ideal_squall2d_x_fp64_cuda", "cuda_dp", "configs/wrf_squall2d_x.ini", ["output.out_interval=20"],
     f"{WREF}/r8/run_squall2d_kessler_10min/wrfout_d01_0001-01-01_00:00:00", 600,
     {"U": 1e-11, "W": 1e-10, "THM": 1e-13, "QV": 1e-12}),
    # FP32 against WRF's single-precision run: WRF r8 itself differs from it by U 8e-5, THM 5e-7
    ("ideal_squall2d_x_fp32", "cuda", "configs/wrf_squall2d_x.ini", ["output.out_interval=200"],
     f"{WREF}/run_squall2d_kessler/wrfout_d01_0001-01-01_00:00:00", 600,
     {"U": 3e-4, "W": 1e-3, "THM": 2e-6, "QV": 1e-4}),
    ("ideal_quarter_ss_fp64_cuda", "cuda_dp", "configs/wrf_quarter_ss.ini",
     ["run.run_time=7200", "output.out_interval=50", "output.diag_interval=50"],
     f"{WREF}/r8/run_qss/wrfout_d01_0001-01-01_00:00:00", 7200,
     {"U": 1e-12, "W": 1e-11, "THM": 1e-14, "QV": 1e-12}),
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
    ap.add_argument("--only")
    ap.add_argument("--keep", action="store_true", help="keep the outputs of passing cases (default: deleted)")
    a = ap.parse_args()

    def done(path, good):  # the outputs of passing cases are deleted (they take tens of GB)
        if good and not a.keep:
            shutil.rmtree(path, ignore_errors=True)

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
        done(out, not bad)
    if (not a.only or "straka" in a.only) and a.cpu_dp and os.path.exists(a.cpu_dp):
        good, msg = straka(a.cpu_dp, os.path.join(a.out, "straka"))
        if good is None:
            print(f"SKIP  straka_wrf             ({msg})")
        else:
            print(f"{'PASS' if good else 'FAIL'}  straka_wrf             {msg}")
            ok &= good
            done(os.path.join(a.out, "straka"), good)
    if a.cpu_dp and os.path.exists(a.cpu_dp):
        for tag, cfg, rd in (("p1", "configs/real_p1_r8.ini", "run_P1_r8"),
                             ("p7", "configs/real_p7_r8.ini", "run_P7_r8"),
                             ("p8", "configs/real_p8_r8.ini", "run_P8_r8"),
                             ("p9", "configs/real_p9_r8.ini", "run_P9_r8"),
                             ("nmp1", "configs/real_nmp1_r8.ini", "run_NMP1_r8"),
                             ("t1", "configs/real_t1.ini", "run_R1_r8")):
            rname = f"restart_{tag}_fp64"
            if a.only and a.only not in rname:
                continue
            good, msg = restart_p1(a.cpu_dp, os.path.join(a.out, f"restart_{tag}"), cfg, rd)
            if good is None:
                print(f"SKIP  {rname:22s} ({msg})")
            else:
                print(f"{'PASS' if good else 'FAIL'}  {rname:22s} {msg}")
                ok &= good
                done(os.path.join(a.out, f"restart_{tag}"), good)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
