"""HY1: hypsometric_opt 1.

CASES: see tests/verify/cases/__init__.py."""
from cases.common import WR  # noqa: F401

CASES = [
    # HY1 = P1 with hypsometric_opt = 1 (real.exe and wrf.exe namelists; wrf.exe takes the namelist value)
    ("real_hy1_r8_fp64_cpu", "cpu_dp", "configs/real_p1_r8.ini", ["real.dir=" + f"{WR}/run_HY1_r8/steps",
                                                                  "dyn.hypsometric_opt=1"],
     f"{WR}/run_HY1_r8/steps/wrfout_d01_2026-10-12_00:00:00", 180,
     {"U": 0.0, "V": 0.0, "W": 0.0, "PH": 0.0, "THM": 0.0, "QV": 0.0, "MU": 0.0}),
    ("real_hy1_r8_fp64_cpu_1h", "cpu_dp", "configs/real_p1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_HY1_r8", "dyn.hypsometric_opt=1", "output.out_interval=10",
      "output.diag_interval=10"],
     f"{WR}/run_HY1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
    ("real_hy1_r8_fp64_cuda_1h", "cuda_dp", "configs/real_p1_r8.ini",
     ["run.run_time=3600", "real.dir=" + f"{WR}/run_HY1_r8", "dyn.hypsometric_opt=1", "output.out_interval=10",
      "output.diag_interval=10"],
     f"{WR}/run_HY1_r8/wrfout_d01_2026-10-12_00:00:00", 3600, {"U": 0.0, "W": 0.0, "THM": 0.0, "QV": 0.0}),
]

# the rows of docs/STATUS.md's table of the CPU FP64 cases (wrf_suite.py --status-table): case, configuration,
# duration
STATUS = [
    ('HY1', 'P1 with hypsometric_opt 1 (real.exe and wrf.exe)', '1 h'),
]

# the data tars (tests/verify/vm_data.txt) each case reads, in that file's order (tests/verify/case_data.py
# derives them from the tars' file lists; tools/colab/colab_suite.sh --cases fetches only these)
DATA = {
    'real_hy1_r8_fp64_cpu': 'hy1_data.tar,v05_files_data.tar,v05_files_links.tar',
    'real_hy1_r8_fp64_cpu_1h': 'hy1_data.tar,v05_suite_refs.tar',
    'real_hy1_r8_fp64_cuda_1h': 'hy1_data.tar,v05_suite_refs.tar',
}
