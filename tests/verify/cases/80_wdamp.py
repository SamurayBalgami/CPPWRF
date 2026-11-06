# SPDX-License-Identifier: GPL-3.0-or-later
"""W1: active w_damping (dt 450 s).

CASES: see tests/verify/cases/__init__.py."""
from cases.common import WR  # noqa: F401

CASES = [
    ("wdamp_w1_fp64", "cpu_dp", "configs/real_r3.ini",
     ["real.dir=" + f"{WR}/run_W1_r8", "time.dt=450", "run.run_time=4500", "output.out_interval=1",
      "output.diag_interval=1"], f"{WR}/run_W1_r8/wrfout_d01_2026-10-12_00:00:00", 4500,
     {"U": 0.0, "W": 0.0, "THM": 0.0}),
    ("wdamp_w1_fp64_cuda", "cuda_dp", "configs/real_r3.ini",
     ["real.dir=" + f"{WR}/run_W1_r8", "time.dt=450", "run.run_time=4500", "output.out_interval=1",
      "output.diag_interval=1"], f"{WR}/run_W1_r8/wrfout_d01_2026-10-12_00:00:00", 4500,
      {"U": 0.0, "W": 0.0, "THM": 0.0}),
]

# the rows of docs/STATUS.md's table of the CPU FP64 cases (wrf_suite.py --status-table): case, configuration,
# duration
STATUS = [
    ('W1', 'active w_damping (dt 450 s)', '10 steps'),
]

# the data tars (tests/verify/vm_data.txt) each case reads, in that file's order (tests/verify/case_data.py
# derives them from the tars' file lists; tools/colab/colab_suite.sh --cases fetches only these)
DATA = {
    'wdamp_w1_fp64': 'suite_data.tar,v05_suite_refs.tar',
    'wdamp_w1_fp64_cuda': 'suite_data.tar,v05_suite_refs.tar',
}
