"""Long runs (24 h and 72 h, slow: wrf_suite.py --slow).

CASES: see tests/verify/cases/__init__.py."""
from cases.common import WR  # noqa: F401

CASES = [
    # Long runs (24 h and 72 h): P4, P6, P9 and N3 from real.exe inputs (wrfbdy every 6 h; 24-72 h: GFS through
    # tools/gcp/gcp_wps.sh), history every 3 h, restarted in the middle from CPPWRF's own wrfrst files
    # (suite.restart_at); every 3-hourly frame and the wrfrst files at the restart equal WRF's uninterrupted run
    # (tools/gcp/gcp_wrfref.sh; references: long_data.tar, 72 h long72_data.tar, long72n3_data.tar). All are slow
    # (suite.slow: wrf_suite.py --slow, gcp_run.sh --slow; the 24 h ones made the VM suite 65 instead of 43 min)
    *[(f"real_{c.lower()}_r8_fp64_cpu_{h}h", "cpu_dp", f"configs/real_{c.lower()}_r8.ini",
       [f"run.run_time={h * 3600}", "real.dir=" + f"{WR}/run_{c}_{h}h_r8", "output.out_interval=180",
        "output.diag_interval=180", f"suite.restart_at={h * 1800}"] + ["suite.slow"],
       f"{WR}/run_{c}_{h}h_r8/wrfout_d01_2026-10-12_00:00:00", h * 3600, {v: 0.0 for v in vs})
      for h in (24, 72)
      for c, vs in (("P4", ("U", "W", "THM", "QV", "QCLOUD", "QRAIN", "QICE")),
                    ("P6", ("U", "W", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QNICE", "QNRAIN")),
                    ("P9", ("U", "W", "THM", "QV")))],
    *[(f"real_n3_r8_fp64_cpu_{h}h", "cpu_dp", "configs/real_n3.ini",
       [f"run.run_time={h * 3600}", "real.dir=" + f"{WR}/run_N3_{h}h_r8", "phys.table_dir=" + f"{WR}/run_N3_{h}h_r8",
        "output.out_interval=180", "output.diag_interval=180", "d02.out_interval=540",
        f"suite.restart_at={h * 1800}"] + ["suite.slow"],
       {d: f"{WR}/run_N3_{h}h_r8/wrfout_{d}_2026-10-12_00:00:00" for d in ("d01", "d02")}, h * 3600,
       {f"{d}.{x}": 0.0 for d in ("d01", "d02") for x in ("U", "W", "THM", "QV", "QCLOUD", "QRAIN", "QICE", "QSNOW",
                                                         "QGRAUP")})
      for h in (24, 72)],
]

# the rows of docs/STATUS.md's table of the CPU FP64 cases (wrf_suite.py --status-table): case, configuration,
# duration
STATUS = [
    ('P4, P6, P9, N3 long runs',
     'inputs from real.exe for 24 h / 72 h (wrfbdy every 6 h; 24-72 h from the GFS forecast through WPS), '
     'history every 3 h, CPPWRF restarted in the middle from its own wrfrst files',
     '24 h; 72 h (P4, P6, P9)'),
]

# the data tars (tests/verify/vm_data.txt) each case reads, in that file's order (tests/verify/case_data.py
# derives them from the tars' file lists; tools/colab/colab_suite.sh --cases fetches only these); the long runs' tars
# are not in vm_data.txt (gcp_run.sh --slow) and were not listed: from the comment above, unchecked
DATA = {
    'real_p4_r8_fp64_cpu_24h': 'long_data.tar',
    'real_p6_r8_fp64_cpu_24h': 'long_data.tar',
    'real_p9_r8_fp64_cpu_24h': 'long_data.tar',
    'real_p4_r8_fp64_cpu_72h': 'long72_data.tar',
    'real_p6_r8_fp64_cpu_72h': 'long72_data.tar',
    'real_p9_r8_fp64_cpu_72h': 'long72_data.tar',
    'real_n3_r8_fp64_cpu_24h': 'long_data.tar',
    'real_n3_r8_fp64_cpu_72h': 'long72n3_data.tar',
}
