"""Land surface: SLAB1 / SLAB4 (5-layer thermal diffusion), NMP1 (Noah-MP).

CASES: see tests/verify/cases/__init__.py."""
from cases.common import WR  # noqa: F401

CASES = [
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
]

# the rows of docs/STATUS.md's table of the CPU FP64 cases (wrf_suite.py --status-table): case, configuration,
# duration
STATUS = [
    ('SLAB1 / SLAB4', 'P1 / P4 with the 5-layer thermal diffusion land surface (sf_surface_physics 1; inputs from real.exe)', '1 h'),
    ('NMP1', 'P1 with Noah-MP', '1 h'),
]
