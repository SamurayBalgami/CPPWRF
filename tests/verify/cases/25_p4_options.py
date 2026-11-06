# SPDX-License-Identifier: GPL-3.0-or-later
"""Options on the P4 / V1 physics: A05 / A05C / A05G (cu_rad_feedback, cldovrlp 0), B01 (adaptive time
step), A08 (buckets, prec_acc_dt), A09 (output_diagnostics), C01 (aerosol-aware Thompson), RADW / RADT
(do_radar_ref), KESSKF (Kessler + Kain-Fritsch), DUD (Dudhia SW), RRTM (RRTM LW).

CASES: see tests/verify/cases/__init__.py."""
from cases.common import WR  # noqa: F401

CASES = [
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
]

# the rows of docs/STATUS.md's table of the CPU FP64 cases (wrf_suite.py --status-table): case, configuration,
# duration
STATUS = [
    ('A05 / A05C / A05G',
     'P4 with cu_rad_feedback (Kain-Fritsch cloud fraction and water in RRTMG) / P4 with cldovrlp 0 (RRTMG'
     ' without clouds) / P8 + WSM6 + RRTMG with cu_rad_feedback (Grell-Freitas QC_CU / QI_CU)',
     '1 h'),
    ('A08',
     'P4 with bucket_mm 0.05, bucket_J 1e4 and prec_acc_dt 2 min (RAINC / RAINNC and the LW accumulations '
     'tip at the first step: I_RAINC, I_ACLW*; PREC_ACC_C, PREC_ACC_NC, SNOW_ACC_NC)',
     '3 min, 1 h'),
    ('A09',
     'P4 with output_diagnostics = 1, auxhist3_interval 3 min (the auxhist3 files T2MIN ... RAINNCVSTD '
     'compared as wrfout)',
     '3 min, 1 h, 1 h restarted at 30 min'),
    ('C01',
     "P6 with the aerosol-aware Thompson (mp_physics 28, aer_init_opt 0: thompson_init's QNWFA / QNIFA "
     'profiles and QNWFA2D emission)',
     '1 h'),
    ('RADW / RADT',
     "P4 (WSM6) / P6 (Thompson) with do_radar_ref 1 (REFL_10CM at the history times: radar_init's bins, "
     'Rayleigh scattering of rain, snow and graupel, melting snow and graupel as Maxwell-Garnett water-ice'
     ' spheroids)',
     '3 min, 1 h'),
    ('KESSKF', "P4 with Kessler instead of WSM6 (Kain-Fritsch rain in Kessler's split sedimentation steps)", '1 h'),
    ('DUD',
     'V1 with the Dudhia shortwave and RRTMG longwave (00 UTC: sunrise over the Bay of Bengal, the whole '
     'domain sunlit within the hour)',
     '1 h'),
    ('RRTM', 'V1 with the RRTM longwave and RRTMG shortwave', '3 steps, 1 h'),
]

# the data tars (tests/verify/vm_data.txt) each case reads, in that file's order (tests/verify/case_data.py
# derives them from the tars' file lists; tools/colab/colab_suite.sh --cases fetches only these)
DATA = {
    'real_a05_r8_fp64_cpu': 'port_a05_data.tar',
    'real_a05_r8_fp64_cpu_1h': 'port_a05_data.tar',
    'real_a05_r8_fp64_cuda_1h': 'port_a05_data.tar',
    'real_a05c_r8_fp64_cpu': 'port_a05_data.tar',
    'real_a05c_r8_fp64_cpu_1h': 'port_a05_data.tar',
    'real_a05c_r8_fp64_cuda_1h': 'port_a05_data.tar',
    'real_a05g_r8_fp64_cpu': 'port_a05_data.tar',
    'real_a05g_r8_fp64_cpu_1h': 'port_a05_data.tar',
    'real_a05g_r8_fp64_cuda_1h': 'port_a05_data.tar',
    'real_b01_r8_fp64_cpu': 'suite_data.tar,port_b01_data.tar,v05_files_links.tar',
    'real_b01_r8_fp64_cpu_1h': 'suite_data.tar,port_b01_data.tar',
    'real_b01_r8_fp64_cuda_1h': 'suite_data.tar,port_b01_data.tar',
    'real_a08_r8_fp64_cpu': 'suite_data.tar,port_a08_data.tar,v05_files_links.tar',
    'real_a08_r8_fp64_cpu_1h': 'suite_data.tar,port_a08_data.tar',
    'real_a08_r8_fp64_cuda_1h': 'suite_data.tar,port_a08_data.tar',
    'real_a09_r8_fp64_cpu': 'port_a09_data.tar',
    'real_a09_r8_fp64_cpu_1h': 'port_a09_data.tar',
    'real_a09_r8_fp64_cpu_1h_rst': 'port_a09_data.tar',
    'real_a09_r8_fp64_cuda_1h': 'port_a09_data.tar',
    'real_c01_r8_fp64_cpu': 'port_a06_data.tar,port_c01_data.tar',
    'real_c01_r8_fp64_cpu_1h': 'port_a06_data.tar,port_c01_data.tar',
    'real_c01_r8_fp64_cuda_1h': 'port_a06_data.tar,port_c01_data.tar',
    'real_radw_r8_fp64_cpu': 'suite_data.tar,port_a10_data.tar,v05_files_links.tar',
    'real_radw_r8_fp64_cpu_1h': 'suite_data.tar,port_a10_data.tar',
    'real_radw_r8_fp64_cuda_1h': 'suite_data.tar,port_a10_data.tar',
    'real_radt_r8_fp64_cpu': 'port_a06_data.tar,port_c01_data.tar,port_a10_data.tar',
    'real_radt_r8_fp64_cpu_1h': 'port_a06_data.tar,port_c01_data.tar,port_a10_data.tar',
    'real_radt_r8_fp64_cuda_1h': 'port_a06_data.tar,port_c01_data.tar,port_a10_data.tar',
    'real_kesskf_r8_fp64_cpu': 'v05_iofields.tar,v05_files_data.tar,v05_files_links.tar',
    'real_kesskf_r8_fp64_cpu_1h': 'suite_data.tar,v05_suite_refs.tar,v05_iofields.tar',
    'real_kesskf_r8_fp64_cuda_1h': 'suite_data.tar,v05_suite_refs.tar,v05_iofields.tar',
    'real_dud_r8_fp64_cpu': 'suite_data.tar,port_c25_data.tar,v05_files_links.tar',
    'real_dud_r8_fp64_cpu_1h': 'suite_data.tar,port_c25_data.tar,v05_files_links.tar',
    'real_dud_r8_fp64_cuda_1h': 'suite_data.tar,port_c25_data.tar,v05_files_links.tar',
    'real_rrtm_r8_fp64_cpu': 'suite_data.tar,port_c25b_data.tar,v05_files_links.tar',
    'real_rrtm_r8_fp64_cpu_1h': 'suite_data.tar,port_c25b_data.tar,v05_files_links.tar',
    'real_rrtm_r8_fp64_cuda_1h': 'suite_data.tar,port_c25b_data.tar,v05_files_links.tar',
}
