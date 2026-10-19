#!/usr/bin/env python3
"""Pack the WRF reference data that the selected wrf_suite.py cases read, for a run on another machine.

  stage_data.py --only NAME[,NAME] [--no-thompson-tables] OUT.tar

Each case's run directory (real.dir of its configuration or of its --set overrides) is written under its own
name, so the remote side sets CPPWRF_WRF_RUNS to the unpacked directory and runs the same suite. Links are
followed. Left out: other history and restart files, rsl/log files, met_em input, executables and the
large tables no ported scheme reads; the Thompson .dat tables only for mp_physics = 8 (and not with
--no-thompson-tables: CPPWRF then computes them, close to but not bit-identical with WRF's files). The reference
wrfout files keep only the frames and variables tests/verify/compare_real.py compares.
"""
import argparse
import fnmatch
import os
import re
import sys
import tarfile
import tempfile

import netCDF4 as nc
import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "tests", "verify"))
import wrf_suite  # noqa: E402

WR = os.environ["CPPWRF_WRF_RUNS"]
SKIP = ["wrfout*", "wrfrst*", "rsl.*", "met_em*", "*.exe", "*.log", "namelist.output", "steps", "README*",
        "p3_lookup*", "ishmael*", "snicar*", "CAM_ABS_DATA", "CAM_AEROPT_DATA", "aerosol*.formatted",
        "CLM_*", "BROADBAND*", "ETAMPNEW*", "CCN_ACTIVATE*", "kernels*", "*.f90*", "eclipse*", "HLC.TBL",
        "tr*", "co2_trans", "grib*", "bulk*", "capacity.asc", "coeff_*", "constants.asc", "masses.asc",
        "termvels.asc", "*.nc"]
THOMPSON = ["freezeH2O.dat", "qr_acr_qg_V4.dat", "qr_acr_qsV2.dat"]
COMPARED = ["Times", "XTIME", "U", "V", "W", "PH", "PHB", "MU", "MUB", "P", "PB", "QVAPOR", "QCLOUD", "QRAIN", "QICE",
            "QSNOW", "QGRAUP", "QNICE", "QNRAIN", "TKE", "RAINNC", "THM"]


def trimmed(src, times, dst):
    """src with only the frames at `times` (seconds) and the compared variables."""
    with nc.Dataset(src) as a, nc.Dataset(dst, "w", format="NETCDF4") as b:
        xt = np.array(a["XTIME"][:], dtype=float) * 60.0
        frames = [i for i, t in enumerate(xt) if any(abs(t - u) < 1e-6 for u in times)]
        for n, d in a.dimensions.items():
            b.createDimension(n, None if d.isunlimited() else len(d))
        b.setncatts({k: a.getncattr(k) for k in a.ncattrs()})
        for n in COMPARED:
            if n not in a.variables:
                continue
            v = a[n]
            w = b.createVariable(n, v.dtype, v.dimensions, zlib=True, complevel=1)
            w.setncatts({k: v.getncattr(k) for k in v.ncattrs()})
            w[:] = v[frames] if v.dimensions and v.dimensions[0] == "Time" else v[:]


def config_value(cfg, section, key):
    sec = None
    for line in open(os.path.join(REPO, cfg)):
        line = line.split("#")[0].strip()
        m = re.match(r"\[(.+)\]", line)
        if m:
            sec = m.group(1)
        elif sec == section and "=" in line and line.split("=", 1)[0].strip() == key:
            return line.split("=", 1)[1].strip()
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", required=True)
    ap.add_argument("--no-thompson-tables", action="store_true")
    ap.add_argument("out")
    a = ap.parse_args()
    want = a.only.split(",")
    dirs = {}  # run directory -> history files to keep, Thompson tables needed
    ref_times = {}  # reference wrfout -> compared times
    for name, key, cfg, sets, wo, t, tol in wrf_suite.CASES:
        if not any(o in name for o in want):
            continue
        rd = next((s.split("=", 1)[1] for s in sets if s.startswith("real.dir=")), None)
        rd = rd or config_value(cfg, "real", "dir")
        if not rd or "${CPPWRF_WRF_RUNS}" not in rd and not rd.startswith(WR):
            print(f"{name}: skipped (not a real-data case of CPPWRF_WRF_RUNS)")
            continue
        rd = rd.replace("${CPPWRF_WRF_RUNS}", WR)
        mp = next((s.split("=", 1)[1] for s in sets if s.startswith("phys.mp_physics=")), None)
        mp = mp or config_value(cfg, "phys", "mp_physics") or "0"
        refs = wo.values() if isinstance(wo, dict) else [wo]
        for r in refs:
            d = dirs.setdefault(os.path.dirname(r), [set(), False])
            d[0].add(os.path.basename(r))
            ref_times.setdefault(r, set()).add(float(t))
        d = dirs.setdefault(rd, [set(), False])
        d[1] = d[1] or (mp.strip() == "8" and not a.no_thompson_tables)
        print(f"{name}: {os.path.relpath(rd, WR)}")
    with tarfile.open(a.out, "w", dereference=True) as tf:
        for d, (keep, thompson) in sorted(dirs.items()):
            for f in sorted(os.listdir(d)):
                p = os.path.join(d, f)
                if os.path.isdir(p) or (f in THOMPSON and not thompson):
                    continue
                arc = os.path.join("wrf_runs", os.path.relpath(p, WR))
                if p in ref_times:
                    with tempfile.TemporaryDirectory() as tmp:
                        q = os.path.join(tmp, f)
                        trimmed(p, ref_times[p], q)
                        tf.add(q, arcname=arc)
                elif f in THOMPSON or not any(fnmatch.fnmatch(f, s) for s in SKIP):
                    tf.add(p, arcname=arc, recursive=False)
    print(f"{a.out}: {os.path.getsize(a.out) / 2**20:.0f} MiB, {len(dirs)} directories")


if __name__ == "__main__":
    main()
