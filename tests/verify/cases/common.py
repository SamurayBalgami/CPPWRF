# SPDX-License-Identifier: GPL-3.0-or-later
"""The data roots the case definitions read: WRF reference runs (wrfinput/wrfbdy and wrfout of each case; the
configurations read them as ${CPPWRF_WRF_RUNS}) and the ideal-case references (${CPPWRF_WRF_REF})."""
import os

WR = os.environ.setdefault("CPPWRF_WRF_RUNS", os.path.expanduser("~/wrf_real480"))
os.environ.setdefault("CPPWRF_WRF_REF", os.path.expanduser("~/wrf_ref480"))
WREF = os.environ["CPPWRF_WRF_REF"]
