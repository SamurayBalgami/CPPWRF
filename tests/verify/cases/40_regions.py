"""V1-V3: other regions and map projections.

CASES: see tests/verify/cases/__init__.py."""
from cases.common import WR  # noqa: F401

CASES = [
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
]

# the rows of docs/STATUS.md's table of the CPU FP64 cases (wrf_suite.py --status-table): case, configuration,
# duration
STATUS = [
    ('V1', 'tropics (Bay of Bengal), Mercator, 15 km, P4 physics', '6 h'),
    ('V2', 'Arctic (Svalbard), polar stereographic, 15 km, P4 physics (sea ice, glacier land ice)', '6 h'),
    ('V3', 'Alps, Lambert, 3 km, Thompson + MYJ + Eta + Noah-MP + RRTMG (radt 5), no cumulus', '2 h'),
]
