"""P5 (WSM6), P2 (Kain-Fritsch), P3 (RRTMG), P4 (Kain-Fritsch + WSM6 + RRTMG).

CASES: see tests/verify/cases/__init__.py."""
from cases.common import WR  # noqa: F401

CASES = [
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
]

# the rows of docs/STATUS.md's table of the CPU FP64 cases (wrf_suite.py --status-table): case, configuration,
# duration
STATUS = [
    ('P5', 'P1 + WSM6 + diff_opt 1 / km_opt 4', '1 h'),
    ('P2', 'P1 + Kain-Fritsch', '1 h'),
    ('P3', 'P1 + RRTMG LW/SW (radt 1 min)', '3 steps'),
    ('P4', 'P2 + WSM6 + RRTMG (cloudy, effective radii) + diff_opt 1 / km_opt 4', '1 h'),
]
