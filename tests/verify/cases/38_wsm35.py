# SPDX-License-Identifier: GPL-3.0-or-later
"""WSM5 and WSM3 microphysics
(mp_physics = 4 / 3, phys/module_mp_wsm5.F / module_mp_wsm3.F): W5 / W3 = P4 (Kain-Fritsch, RRTMG with the schemes'
effective radii, YSU, Noah) with WSM5 / WSM3 instead of WSM6, each against its WRF v4.8.0 real*8 run (wrfref workflow
of Xassemblianist/cppwrf-data, tag c02-wsm35: the inputs of run_P4_r8 with physics.mp_physics = 4 / 3; wrfrst at
30 min).

CASES: see tests/verify/cases/__init__.py."""
from cases.common import WR  # noqa: F401

_D = "2026-10-12_00:00:00"
_T5 = {"U": 0.0, "V": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0,
       "QICE": 0.0, "QSNOW": 0.0, "RAINNC": 0.0, "SNOWNC": 0.0, "RE_CLOUD": 0.0, "RE_ICE": 0.0, "RE_SNOW": 0.0}
_T3 = {k: v for k, v in _T5.items() if k not in ("QICE", "QSNOW")}
_1H = ["run.run_time=3600", "output.out_interval=10", "output.diag_interval=10"]

CASES = []
for _n, _t in (("5", _T5), ("3", _T3)):
    _cfg, _run = f"configs/real_p4w{_n}_r8.ini", f"{WR}/run_W{_n}_r8"
    CASES += [
        # 3 steps (history every step), 1 h, and 1 h restarted at 30 min (its wrfrst identical to WRF's)
        (f"real_p4w{_n}_r8_fp64_cpu", "cpu_dp", _cfg, [], f"{_run}/steps/wrfout_d01_{_D}", 180, _t),
        (f"real_p4w{_n}_r8_fp64_cpu_1h", "cpu_dp", _cfg, ["real.dir=" + _run] + _1H, f"{_run}/wrfout_d01_{_D}",
         3600, _t),
        (f"real_p4w{_n}_r8_fp64_cpu_1h_rst", "cpu_dp", _cfg, ["real.dir=" + _run] + _1H + ["suite.restart_at=1800"],
         f"{_run}/wrfout_d01_{_D}", 3600, _t),
        # CUDA FP64 (exact on the G4, 2026-11-16) and FP32 (measured on the G4: W5 U 6.3e-4, W 1.7e-2, THM 8.4e-6,
        # QV 5.4e-4; W3 about the same; tolerances about 3x)
        (f"real_p4w{_n}_r8_fp64_cuda_1h", "cuda_dp", _cfg, ["real.dir=" + _run] + _1H, f"{_run}/wrfout_d01_{_D}",
         3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
        (f"real_p4w{_n}_fp32_1h", "cuda", _cfg, ["real.dir=" + _run] + _1H, f"{_run}/wrfout_d01_{_D}", 3600,
         {"U": 1.9e-3, "W": 5.1e-2, "THM": 2.6e-5, "QV": 1.6e-3}),
    ]

STATUS = [
    ('W5', 'P4 with WSM5 (mp_physics 4) instead of WSM6; restarted at 30 min', '1 h'),
    ('W3', 'P4 with WSM3 (mp_physics 3) instead of WSM6; restarted at 30 min', '1 h'),
]

DATA = {c[0]: 'suite_data.tar,port_c02_data.tar' for c in CASES}
