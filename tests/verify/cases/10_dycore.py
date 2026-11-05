"""Real-data dynamics: R1 (dynamics only), R1D6 (diff_opt 1 / km_opt 4, damp_opt 1, 6th-order
diffusion), T1 (diff_opt 2 / km_opt 2, TKE), the FP32 runs R0, R2, R3.

CASES: see tests/verify/cases/__init__.py."""
from cases.common import WR  # noqa: F401

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
]

# the rows of docs/STATUS.md's table of the CPU FP64 cases (wrf_suite.py --status-table): case, configuration,
# duration
STATUS = [
    ('R1', 'dynamics only', '5 steps'),
    ('T1', 'R1 with diff_opt 2 / km_opt 2 (TKE)', '5 steps'),
    ('R1D6', 'R1 with diff_opt 1 / km_opt 4, damp_opt 1, diff_6th_opt 2', '10 min'),
]

# the data tars (tests/verify/vm_data.txt) each case reads, in that file's order (tests/verify/case_data.py
# derives them from the tars' file lists; tools/colab/colab_suite.sh --cases fetches only these)
DATA = {
    'real_r1_fp64_cpu': 'dycore_data.tar,v05_suite_refs.tar',
    'real_r1d6_fp64_cpu': 'dycore_data.tar,v05_suite_refs.tar',
    'real_r1d6_fp64_cuda': 'dycore_data.tar,v05_suite_refs.tar',
    'real_r1_fp64_cuda': 'dycore_data.tar,v05_suite_refs.tar',
    'real_t1_fp64_cpu': 'suite_data.tar,v05_suite_refs.tar',
    'real_t1_fp64_cuda': 'suite_data.tar,v05_suite_refs.tar',
    'real_r0_fp32_1h': 'suite_data.tar',
    'real_r2_fp32_1h': 'suite_data.tar',
    'real_r3_fp32_1h': 'suite_data.tar',
}
