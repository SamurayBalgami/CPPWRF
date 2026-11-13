# SPDX-License-Identifier: GPL-3.0-or-later
"""Option values that CPPWRF accepted without a whole-model comparison with WRF before (the coverage audit's
"unreached WRF code paths"): the RRTMG McICA cloud overlaps cldovrlp 1, 3, 4, 5 and idcor 1 (OVL1, OVL3, OVL4I,
OVL5), icloud 0 (ICL0), diff_6th_opt 1 (D6O1), km_opt 3 (KM3 with diff_opt 2, KM3H with diff_opt 1), tke_adv_opt 0 with km_opt 2 (TKA0) on
the P4 physics;
moist_adv_opt 0 / scalar_adv_opt 0 (ADV0) and bl_mynn_cloudpdf -1 / -2 (MPDFM1, MPDFM2) on the MYNN physics. Each
against its WRF v4.8.0 real*8 run (tools/colab/colab_wrfref.sh, the inputs of run_P4_r8 / run_P6_r8): 3 steps with
history every step, and 1 h with history every 10 min, CPU FP64 and CUDA FP64, every variable at tolerance 0.
KM3 (diff_opt 2 / km_opt 3) is compared over its first two steps only: WRF itself blows up on this 15 km domain in
the third step (NaN everywhere after it, max|U| 256 m/s after the second; CPPWRF's fatal non-finite check stops
there), so km_opt 3 runs 1 h in KM3H with diff_opt 1.

The other MYNN values the audit listed (bl_mynn_edmf 1, bl_mynn_edmf_mom 1, bl_mynn_mixlength 1, bl_mynn_cloudmix 1,
bl_mynn_cloudpdf 2, bl_mynn_closure 2.6) are WRF's defaults: MYNN, MYNNSF1 / MYNNSF2 and (but for the values it
sets) MYNNX1 run them already (tests/verify/cases/35_mynn.py).

CASES: see tests/verify/cases/__init__.py."""
from cases.common import WR  # noqa: F401

_D = "2026-10-12_00:00:00"
_P4 = {"U": 0.0, "V": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0, "QCLOUD": 0.0, "QRAIN": 0.0,
       "QICE": 0.0, "CLDFRA": 0.0, "SWDOWN": 0.0, "GLW": 0.0}
_MY = {v: 0.0 for v in ("U", "V", "W", "PH", "MU", "P", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP",
                        "QNICE", "QNRAIN", "RAINNC", "QKE", "CLDFRA")}
_1H = ["run.run_time=3600", "output.out_interval=10", "output.diag_interval=10"]

# (tag, WRF run, tolerances)
_SETS = [("ovl1", "OVL1", _P4), ("ovl3", "OVL3", _P4), ("ovl4i", "OVL4I", _P4), ("ovl5", "OVL5", _P4),
         ("icl0", "ICL0", _P4), ("d6o1", "D6O1", _P4), ("km3h", "KM3H", _P4), ("tka0", "TKA0", _P4),
         ("adv0", "ADV0", _MY), ("mpdfm1", "MPDFM1", _MY), ("mpdfm2", "MPDFM2", _MY)]

# KM3: two steps (WRF's run is NaN from the third step on), CPU and CUDA FP64
CASES = [(f"real_km3_r8_fp64{b}", p, "configs/real_km3_r8.ini", [], f"{WR}/run_KM3_r8/steps/wrfout_d01_{_D}", 120, _P4)
         for b, p in (("_cpu", "cpu_dp"), ("_cuda", "cuda_dp"))]
CASES += [c for tag, ref, tol in _SETS for c in (
    (f"real_{tag}_r8_fp64_cpu", "cpu_dp", f"configs/real_{tag}_r8.ini", [],
     f"{WR}/run_{ref}_r8/steps/wrfout_d01_{_D}", 180, tol),
    (f"real_{tag}_r8_fp64_cpu_1h", "cpu_dp", f"configs/real_{tag}_r8.ini", ["real.dir=" + f"{WR}/run_{ref}_r8"] + _1H,
     f"{WR}/run_{ref}_r8/wrfout_d01_{_D}", 3600, tol),
    (f"real_{tag}_r8_fp64_cuda_1h", "cuda_dp", f"configs/real_{tag}_r8.ini", ["real.dir=" + f"{WR}/run_{ref}_r8"] + _1H,
     f"{WR}/run_{ref}_r8/wrfout_d01_{_D}", 3600, tol))]

# the rows of docs/STATUS.md's table of the CPU FP64 cases (wrf_suite.py --status-table): case, configuration,
# duration
STATUS = [
    ('OVL1', 'P4 with cldovrlp 1 (RRTMG McICA random overlap)', '1 h'),
    ('OVL3', 'P4 with cldovrlp 3 (maximum overlap)', '1 h'),
    ('OVL4I', 'P4 with cldovrlp 4 (exponential overlap) and idcor 1 (latitude-dependent decorrelation length)', '1 h'),
    ('OVL5', 'P4 with cldovrlp 5 (exponential-random overlap, constant decorrelation length)', '1 h'),
    ('ICL0', 'P4 with icloud 0 (no cloud fraction in RRTMG)', '1 h'),
    ('D6O1', 'P4 with diff_6th_opt 1 (sixth-order diffusion without the flux limiter)', '1 h'),
    ('KM3', 'P4 with diff_opt 2 / km_opt 3 (3-D Smagorinsky; WRF is NaN from the third step on)', '2 steps'),
    ('KM3H', 'P4 with diff_opt 1 / km_opt 3 (3-D Smagorinsky Km, horizontal diffusion along the model surfaces)', '1 h'),
    ('TKA0', 'P4 with diff_opt 2 / km_opt 2 (1.5-order TKE) and tke_adv_opt 0', '1 h'),
    ('ADV0', 'MYNN with moist_adv_opt 0 and scalar_adv_opt 0', '1 h'),
    ('MPDFM1', 'MYNN with bl_mynn_cloudpdf -1 (subgrid clouds zeroed)', '1 h'),
    ('MPDFM2', 'MYNN with bl_mynn_cloudpdf -2 (subgrid clouds zeroed)', '1 h'),
]

# the data tars (tests/verify/vm_data.txt) each case reads, in that file's order (tests/verify/case_data.py
# derives them from the tars' file lists; tools/colab/colab_suite.sh --cases fetches only these)
DATA = {c[0]: 'port_cov1b_data.tar' if '_km3h_' in c[0] else 'port_cov1_data.tar' for c in CASES}
