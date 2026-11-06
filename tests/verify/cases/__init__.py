# SPDX-License-Identifier: GPL-3.0-or-later
"""The cases of tests/verify/wrf_suite.py, one file per group: NN_<group>.py (loaded in the order of their names, so
the case order is fixed; a new group takes a free number) defines CASES, a list of
  (name, binary key, config, extra --set, WRF wrfout, time [s], {var: max rel L2 (whole domain)})
(binary key cpu_dp, cuda_dp or cuda; WRF wrfout: a file, or {domain: file} for nests; tolerances all 0.0: compared
exactly in every variable of the wrfout file, the list names the variables that matter most for the case).
A group file may define STATUS, its rows of docs/STATUS.md's table of the CPU FP64 cases (wrf_suite.py
--status-table), and DATA, the data tars (tests/verify/vm_data.txt) each case reads, {name: "TAR,TAR"} in
vm_data.txt's order (tests/verify/case_data.py derives and checks them; tools/colab/colab_suite.sh --cases fetches
only those). restarts.py defines RESTARTS, the restart checks, and DATA for wrf_suite.py's special checks (Straka, the
restart checks). Case names are unique."""
import glob
import importlib.util
import os
import sys

_DIR = os.path.dirname(os.path.abspath(__file__))
if os.path.dirname(_DIR) not in sys.path:  # the group files import cases.common
    sys.path.insert(0, os.path.dirname(_DIR))


def _load(path):
    spec = importlib.util.spec_from_file_location(f"cases._{os.path.basename(path)[:-3]}", path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


GROUPS = [_load(p) for p in sorted(glob.glob(os.path.join(_DIR, "[0-9][0-9]_*.py")))]
CASES = [c for g in GROUPS for c in g.CASES]
STATUS = [r for g in GROUPS for r in getattr(g, "STATUS", [])]
_dup = sorted({c[0] for c in CASES if [x[0] for x in CASES].count(c[0]) > 1})
if _dup:
    raise SystemExit(f"tests/verify/cases: case names defined twice: {', '.join(_dup)}")
from cases.restarts import RESTARTS  # noqa: E402,F401
import cases.restarts as _restarts  # noqa: E402
# the data tars each case and special check reads ("TAR,TAR" in tests/verify/vm_data.txt's order)
DATA = {n: t for g in GROUPS + [_restarts] for n, t in getattr(g, "DATA", {}).items()}
