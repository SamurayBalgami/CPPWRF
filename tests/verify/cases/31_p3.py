# SPDX-License-Identifier: GPL-3.0-or-later
"""The P3 microphysics with one ice category (phys/module_mp_p3.F): P3MP (mp_physics = 50, constant droplet number)
and P3NC (mp_physics = 51, predicted droplet number QNCLOUD) on the P4 physics with Kain-Fritsch, each against its WRF
v4.8.0 real*8 run of the -fdefault-double-8 build (tools/colab/colab_wrfref.sh, the inputs of run_P4_r8).

CASES: see tests/verify/cases/__init__.py."""
from cases.common import WR  # noqa: F401

_D = "2026-10-12_00:00:00"
_T = {"U": 0.0, "V": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0, "QICE": 0.0, "QNICE": 0.0,
      "QNRAIN": 0.0, "QIR": 0.0, "QIB": 0.0, "REFL_10CM": 0.0, "RAINNC": 0.0, "SNOWNC": 0.0}
_TNC = dict(_T, QNCLOUD=0.0)
_1H = ["run.run_time=3600", "output.out_interval=10", "output.diag_interval=10"]

CASES = [
    # P3MP = P4 with P3 (50): 3 steps (history every step), 1 h, 1 h restarted at 30 min (wrfrst identical)
    ("real_p3mp_r8_fp64_cpu", "cpu_dp", "configs/real_p3mp_r8.ini", [],
     f"{WR}/run_P3MP_r8/steps/wrfout_d01_{_D}", 180, _T),
    ("real_p3mp_r8_fp64_cpu_1h", "cpu_dp", "configs/real_p3mp_r8.ini", ["real.dir=" + f"{WR}/run_P3MP_r8"] + _1H,
     f"{WR}/run_P3MP_r8/wrfout_d01_{_D}", 3600, _T),
    ("real_p3mp_r8_fp64_cpu_1h_rst", "cpu_dp", "configs/real_p3mp_r8.ini",
     ["real.dir=" + f"{WR}/run_P3MP_r8"] + _1H + ["suite.restart_at=1800"],
     f"{WR}/run_P3MP_r8/wrfout_d01_{_D}", 3600, _T),
    # GPU: FP64 exact (G4 2026-11-15); FP32 measured on G4 2026-11-15: U 4.3e-3, W 1.3e-1, THM 5.6e-5, QV 2.0e-3
    # (P4's KF columns, docs/FP32_ACCURACY.md); tolerances about 3x
    ("real_p3mp_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p3mp_r8.ini", ["real.dir=" + f"{WR}/run_P3MP_r8"] + _1H,
     f"{WR}/run_P3MP_r8/wrfout_d01_{_D}", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0, "QNICE": 0.0}),
    ("real_p3mp_fp32_1h", "cuda", "configs/real_p3mp_r8.ini", ["real.dir=" + f"{WR}/run_P3MP_r8"] + _1H,
     f"{WR}/run_P3MP_r8/wrfout_d01_{_D}", 3600, {"U": 1.3e-2, "W": 0.4, "THM": 1.7e-4, "QV": 6e-3}),
    # P3NC = P4 with P3 (51): 3 steps, 1 h
    ("real_p3nc_r8_fp64_cpu", "cpu_dp", "configs/real_p3nc_r8.ini", [],
     f"{WR}/run_P3NC_r8/steps/wrfout_d01_{_D}", 180, _TNC),
    ("real_p3nc_r8_fp64_cpu_1h", "cpu_dp", "configs/real_p3nc_r8.ini", ["real.dir=" + f"{WR}/run_P3NC_r8"] + _1H,
     f"{WR}/run_P3NC_r8/wrfout_d01_{_D}", 3600, _TNC),
]

# the rows of docs/STATUS.md's table of the CPU FP64 cases (wrf_suite.py --status-table): case, configuration,
# duration
STATUS = [
    ('P3MP', 'P4 with P3, one ice category (mp_physics 50) and KF (restarted at 30 min)', '1 h'),
    ('P3NC', 'P4 with P3 with predicted droplet number (mp_physics 51) and KF', '1 h'),
]

# the data tars (tests/verify/vm_data.txt) each case reads, in that file's order (tests/verify/case_data.py
# derives them from the tars' file lists; tools/colab/colab_suite.sh --cases fetches only these)
DATA = {c[0]: 'port_c09_data.tar' for c in CASES}
