# SPDX-License-Identifier: GPL-3.0-or-later
"""The BouLac PBL (bl_pbl_physics = 8,
phys/module_bl_boulac.F): BL1 / BL4 / BL9 = P1 / P4 / P9 with BouLac (over the revised MM5 surface layer, P4 with
WSM6 / KF / RRTMG, and over the Eta surface layer), each against its WRF v4.8.0 real*8 run (wrfref workflow of
Xassemblianist/cppwrf-data, tag c18-boulac: the inputs of run_P1_r8 / run_P4_r8 / run_P9_r8 with
physics.bl_pbl_physics = 8; BL1 with wrfrst at 30 min).

CASES: see tests/verify/cases/__init__.py."""
from cases.common import WR  # noqa: F401

_D = "2026-10-12_00:00:00"
_T = {"U": 0.0, "V": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0, "PBLH": 0.0, "HFX": 0.0, "UST": 0.0,
      "TKE_PBL": 0.0, "EL_PBL": 0.0}
_1H = ["run.run_time=3600", "output.out_interval=10", "output.diag_interval=10"]

CASES = []
for _n in ("1", "4", "9"):
    _cfg, _run = f"configs/real_bl{_n}_r8.ini", f"{WR}/run_BL{_n}_r8"
    CASES += [
        # 3 steps (history every step) and 1 h
        (f"real_bl{_n}_r8_fp64_cpu", "cpu_dp", _cfg, [], f"{_run}/steps/wrfout_d01_{_D}", 180, _T),
        (f"real_bl{_n}_r8_fp64_cpu_1h", "cpu_dp", _cfg, ["real.dir=" + _run] + _1H, f"{_run}/wrfout_d01_{_D}",
         3600, _T),
    ]
CASES += [
    # BL1 restarted at 30 min (its wrfrst identical to WRF's)
    ("real_bl1_r8_fp64_cpu_1h_rst", "cpu_dp", "configs/real_bl1_r8.ini",
     ["real.dir=" + f"{WR}/run_BL1_r8"] + _1H + ["suite.restart_at=1800"], f"{WR}/run_BL1_r8/wrfout_d01_{_D}", 3600, _T),
    # CUDA FP64 and FP32 (tolerances to be measured on Colab by the merge train: about 3x the measured rel L2)
    ("real_bl4_r8_fp64_cuda_1h", "cuda_dp", "configs/real_bl4_r8.ini", ["real.dir=" + f"{WR}/run_BL4_r8"] + _1H,
     f"{WR}/run_BL4_r8/wrfout_d01_{_D}", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_bl4_fp32_1h", "cuda", "configs/real_bl4_r8.ini", ["real.dir=" + f"{WR}/run_BL4_r8"] + _1H,
     f"{WR}/run_BL4_r8/wrfout_d01_{_D}", 3600, {"U": 1e-2, "W": 2e-1, "THM": 1e-4, "QV": 1e-2}),
]

STATUS = [
    ('BL1', 'P1 with the BouLac PBL (bl_pbl_physics 8; restarted at 30 min)', '1 h'),
    ('BL4', 'P4 with the BouLac PBL', '1 h'),
    ('BL9', 'P9 with the BouLac PBL over the Eta surface layer', '1 h'),
]

DATA = {c[0]: 'suite_data.tar,port_c18_data.tar' for c in CASES}  # the inputs and tables: links into run_P1_r8 / P4 / P9
