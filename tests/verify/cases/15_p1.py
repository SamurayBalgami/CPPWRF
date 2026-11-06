# SPDX-License-Identifier: GPL-3.0-or-later
"""The P1 column physics (sfclayrev + Noah + YSU) and its variants P1T2, A01, A01V; the 1 h CPU runs of P2,
P7, P8.

CASES: see tests/verify/cases/__init__.py."""
from cases.common import WR  # noqa: F401

CASES = [
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
]

# the rows of docs/STATUS.md's table of the CPU FP64 cases (wrf_suite.py --status-table): case, configuration,
# duration
STATUS = [
    ('P1', 'sfclayrev + Noah + YSU', '1 h'),
    ('A01 / A01V',
     'P1 / V2 with usemonalb and rdlai2d (wrfinput ALBBCK and LAI; V2: snow-covered points start with '
     'ALBEDO = SNOALB)',
     '1 h'),
]

# the data tars (tests/verify/vm_data.txt) each case reads, in that file's order (tests/verify/case_data.py
# derives them from the tars' file lists; tools/colab/colab_suite.sh --cases fetches only these)
DATA = {
    'real_p1_r8_fp64_cpu': 'v05_files_data.tar,v05_files_links.tar',
    'real_p1_r8_fp64_cpu_1h': 'suite_data.tar,v05_suite_refs.tar',
    'real_p1t2_r8_fp64_cpu_1h': 'p1t2_data.tar,v05_suite_refs.tar',
    'real_p1t2_r8_fp64_cuda_1h': 'p1t2_data.tar,v05_suite_refs.tar',
    'real_a01_r8_fp64_cpu': 'suite_data.tar,port_a01_data.tar,v05_files_links.tar',
    'real_a01_r8_fp64_cpu_1h': 'suite_data.tar,port_a01_data.tar',
    'real_a01_r8_fp64_cuda_1h': 'suite_data.tar,port_a01_data.tar',
    'real_a01v_r8_fp64_cpu_1h': 'suite_data.tar,port_a01_data.tar,v05_files_links.tar',
    'real_a01v_r8_fp64_cuda_1h': 'suite_data.tar,port_a01_data.tar,v05_files_links.tar',
    'real_p2_r8_fp64_cpu_1h': 'suite_data.tar,v05_suite_refs.tar',
    'real_p7_r8_fp64_cpu_1h': 'suite_data.tar,v05_suite_refs.tar',
    'real_p8_r8_fp64_cpu_1h': 'suite_data.tar,v05_suite_refs.tar',
    'real_p1_r8_fp64_cuda_1h': 'suite_data.tar,v05_suite_refs.tar',
    'real_p1_fp32_1h': 'suite_data.tar',
}
