# SPDX-License-Identifier: GPL-3.0-or-later
"""Nesting: N1 (one-way), N2 / N2S1 / N2S2 (two-way, smoothers), N3 (two-way with the P4 physics).

CASES: see tests/verify/cases/__init__.py."""
from cases.common import WR  # noqa: F401

CASES = [
    # Nesting: one-way (N1) and two-way (N2, feedback) vs WRF r8, bit-identical on the CPU (both domains;
    # N2 for 1 h)
    ("real_n1_r8_fp64_cpu", "cpu_dp", "configs/real_n1_steps.ini", ["real.dir=" + f"{WR}/run_N1_r8/steps"],
     {d: f"{WR}/run_N1_r8/steps/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 180,
     {"d01.U": 0.0, "d01.THM": 0.0, "d01.QV": 0.0, "d02.U": 0.0, "d02.THM": 0.0, "d02.QV": 0.0}),
    ("real_n1_r8_fp64_cuda", "cuda_dp", "configs/real_n1_steps.ini", ["real.dir=" + f"{WR}/run_N1_r8/steps"],
     {d: f"{WR}/run_N1_r8/steps/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 180,
     {}),
    # N1 for 1 h in FP32 against WRF's single-precision run
    ("real_n1_fp32_1h", "cuda", "configs/real_n1.ini", [],
     {d: f"{WR}/run_N1/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {f"{d}.{x}": t for d in ("d01", "d02") for x, t in (("U", 9e-4), ("W", 1.5e-2), ("THM", 3e-6), ("QV", 1e-4))}),
    ("real_n2_r8_fp64_cpu", "cpu_dp", "configs/real_n2_steps.ini", ["real.dir=" + f"{WR}/run_N2_r8/steps"],
     {d: f"{WR}/run_N2_r8/steps/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 180,
     {"d01.U": 0.0, "d01.THM": 0.0, "d01.QV": 0.0, "d02.U": 0.0, "d02.THM": 0.0, "d02.QV": 0.0}),
    ("real_n2_r8_fp64_cpu_1h", "cpu_dp", "configs/real_n2.ini", ["real.dir=" + f"{WR}/run_N2_r8"],
     {d: f"{WR}/run_N2_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {"d01.U": 0.0, "d01.W": 0.0, "d01.THM": 0.0, "d02.U": 0.0, "d02.W": 0.0, "d02.THM": 0.0, "d02.QV": 0.0}),
    # N2S1 = N2 with smooth_option = 1 (sm121) of the fed-back parent fields
    ("real_n2s1_r8_fp64_cpu", "cpu_dp", "configs/real_n2_steps.ini",
     ["real.dir=" + f"{WR}/run_N2S1_r8/steps", "domains.smooth_option=1"],
     {d: f"{WR}/run_N2S1_r8/steps/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 180,
     {"d01.U": 0.0, "d01.THM": 0.0, "d01.QV": 0.0, "d02.U": 0.0, "d02.THM": 0.0, "d02.QV": 0.0}),
    ("real_n2s1_r8_fp64_cpu_1h", "cpu_dp", "configs/real_n2.ini",
     ["real.dir=" + f"{WR}/run_N2S1_r8", "domains.smooth_option=1"],
     {d: f"{WR}/run_N2S1_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {"d01.U": 0.0, "d01.W": 0.0, "d01.THM": 0.0, "d02.U": 0.0, "d02.W": 0.0, "d02.THM": 0.0, "d02.QV": 0.0}),
    ("real_n2s1_r8_fp64_cuda_1h", "cuda_dp", "configs/real_n2.ini",
     ["real.dir=" + f"{WR}/run_N2S1_r8", "domains.smooth_option=1"],
     {d: f"{WR}/run_N2S1_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {}),
    # N2S2 = N2 with smooth_option = 2 (smdsm, WRF default) of the fed-back parent fields
    ("real_n2s2_r8_fp64_cpu", "cpu_dp", "configs/real_n2_steps.ini",
     ["real.dir=" + f"{WR}/run_N2S2_r8/steps", "domains.smooth_option=2"],
     {d: f"{WR}/run_N2S2_r8/steps/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 180,
     {"d01.U": 0.0, "d01.THM": 0.0, "d01.QV": 0.0, "d02.U": 0.0, "d02.THM": 0.0, "d02.QV": 0.0}),
    ("real_n2s2_r8_fp64_cpu_1h", "cpu_dp", "configs/real_n2.ini",
     ["real.dir=" + f"{WR}/run_N2S2_r8", "domains.smooth_option=2"],
     {d: f"{WR}/run_N2S2_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {"d01.U": 0.0, "d01.W": 0.0, "d01.THM": 0.0, "d02.U": 0.0, "d02.W": 0.0, "d02.THM": 0.0, "d02.QV": 0.0}),
    ("real_n2s2_r8_fp64_cuda_1h", "cuda_dp", "configs/real_n2.ini",
     ["real.dir=" + f"{WR}/run_N2S2_r8", "domains.smooth_option=2"],
     {d: f"{WR}/run_N2S2_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {}),
    # N3 = two-way nest with the P4 physics on both domains (physics feedback, nest ozone from the parent)
    ("real_n3_r8_fp64_cpu", "cpu_dp", "configs/real_n3_steps.ini", [],
     {d: f"{WR}/run_N3s_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 180,
     {f"{d}.{x}": 0.0 for d in ("d01", "d02") for x in ("U", "W", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP")}),
    ("real_n3_r8_fp64_cpu_1h", "cpu_dp", "configs/real_n3.ini", [],
     {d: f"{WR}/run_N3h_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, 3600,
     {f"{d}.{x}": 0.0 for d in ("d01", "d02") for x in ("U", "W", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW", "QGRAUP")}),
]

# the rows of docs/STATUS.md's table of the CPU FP64 cases (wrf_suite.py --status-table): case, configuration,
# duration
STATUS = [
    ('N1', 'one-way nest (both domains)', '3 parent steps'),
    ('N2', 'two-way nest (both domains)', '1 h'),
    ('N2S1 / N2S2', 'N2 with smooth_option 1 / 2 (smoothing of the fed-back fields, ht included)', '1 h'),
    ('N3', 'two-way nest with the P4 physics on both domains (physics feedback, nest ozone from the parent)', '1 h'),
]

# the data tars (tests/verify/vm_data.txt) each case reads, in that file's order (tests/verify/case_data.py
# derives them from the tars' file lists; tools/colab/colab_suite.sh --cases fetches only these)
DATA = {
    'real_n1_r8_fp64_cpu': 'suite_data.tar,v05_files_data.tar',
    'real_n1_r8_fp64_cuda': 'suite_data.tar,v05_files_data.tar',
    'real_n1_fp32_1h': '',
    'real_n2_r8_fp64_cpu': 'suite_data.tar,v05_files_data.tar',
    'real_n2_r8_fp64_cpu_1h': 'v05_suite_refs.tar,v05_files_links.tar',
    'real_n2s1_r8_fp64_cpu': 'nest_data.tar,v05_files_data.tar,v05_files_links.tar',
    'real_n2s1_r8_fp64_cpu_1h': 'nest_data.tar,v05_suite_refs.tar',
    'real_n2s1_r8_fp64_cuda_1h': 'nest_data.tar,v05_suite_refs.tar',
    'real_n2s2_r8_fp64_cpu': 'nest_data.tar,v05_files_data.tar,v05_files_links.tar',
    'real_n2s2_r8_fp64_cpu_1h': 'nest_data.tar,v05_suite_refs.tar',
    'real_n2s2_r8_fp64_cuda_1h': 'nest_data.tar,v05_suite_refs.tar',
    'real_n3_r8_fp64_cpu': 'n3_data.tar,v05_suite_refs.tar',
    'real_n3_r8_fp64_cpu_1h': 'n3_data.tar,v05_suite_refs.tar',
}
