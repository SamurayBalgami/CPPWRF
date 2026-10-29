#!/usr/bin/env python3
"""Pack the WRF reference data that the selected wrf_suite.py cases read, for a run on another machine.

  stage_data.py --only NAME[,NAME] [--files NAME[,NAME]] [--no-thompson-tables] [--whole] [--refs-only] [--no-steps]
                OUT.tar|-

Each case's run directory (real.dir of its configuration or of its --set overrides) is written under its own
name, so the remote side sets CPPWRF_WRF_RUNS to the unpacked directory and runs the same suite. Links are
followed. Left out: other history and restart files, rsl/log files, met_em input, executables and the
large tables no ported scheme reads; the Thompson .dat tables only for mp_physics = 8 (and not with
--no-thompson-tables: CPPWRF then computes them, close to but not bit-identical with WRF's files). The reference
wrfout files keep only the frames the cases compare, and every variable for a tolerance-0 case (wrf_suite.py
compares all of them) or else the variables tests/verify/compare_real.py compares. The ideal cases' references
come with their namelist.input (wrf_suite.py checks the configuration and the global attributes with it).
--whole: the reference files whole (every frame and variable: wrf_suite.py then compares every frame up to the case's
time; it needs no work, so OUT can be - for a tar on stdout, e.g. | gcloud storage cp - gs://BUCKET/X.tar, which
writes nothing locally); --refs-only: only the reference files (and the ideal namelists), not the run directories;
--no-steps: no reference in a steps/ directory (tools/gcp/stage_files_data.py packs those whole).
--files: the cases of tests/verify/run_wrf_files_test.sh run with CPPWRF_FILES_NO_VALUES=1 (a GPU build): their input
directories, and their wrfout / wrfrst references (run, restart and resumed directories) as skeletons, with every
variable, dimension, attribute and time but zero numbers; whole only the wrfrst files of minute 1, which the resumed
runs start from.
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
WREF = os.environ.get("CPPWRF_WRF_REF", os.path.expanduser("~/wrf_ref480"))
SKIP = ["wrfout*", "wrfrst*", "rsl.*", "met_em*", "*.exe", "*.log", "namelist.output", "steps", "README*",
        "p3_lookup*", "ishmael*", "snicar*", "CAM_ABS_DATA", "CAM_AEROPT_DATA", "aerosol*.formatted",
        "CLM_*", "BROADBAND*", "ETAMPNEW*", "CCN_ACTIVATE*", "kernels*", "*.f90*", "eclipse*", "HLC.TBL",
        "tr*", "co2_trans", "grib*", "bulk*", "capacity.asc", "coeff_*", "constants.asc", "masses.asc",
        "termvels.asc", "*.nc"]
THOMPSON = ["freezeH2O.dat", "qr_acr_qg_V4.dat", "qr_acr_qsV2.dat"]
COMPARED = ["Times", "XTIME", "U", "V", "W", "PH", "PHB", "MU", "MUB", "P", "PB", "QVAPOR", "QCLOUD", "QRAIN", "QICE",
            "QSNOW", "QGRAUP", "QNICE", "QNRAIN", "TKE", "RAINNC", "THM"]


def trimmed(src, times, dst, every_variable):
    """src with only the frames at `times` (seconds) and the compared variables (every_variable: all)."""
    with nc.Dataset(src) as a, nc.Dataset(dst, "w", format="NETCDF4") as b:
        xt = np.array(a["XTIME"][:], dtype=float) * 60.0
        frames = [i for i, t in enumerate(xt) if any(abs(t - u) < 1e-6 for u in times)]
        for n, d in a.dimensions.items():
            b.createDimension(n, None if d.isunlimited() else len(d))
        b.setncatts({k: a.getncattr(k) for k in a.ncattrs()})
        for n in (list(a.variables) if every_variable else COMPARED):
            if n not in a.variables:
                continue
            v = a[n]
            w = b.createVariable(n, v.dtype, v.dimensions, zlib=True, complevel=1)
            w.setncatts({k: v.getncattr(k) for k in v.ncattrs()})
            w[:] = v[frames] if v.dimensions and v.dimensions[0] == "Time" else v[:]


def skeleton(src, dst):
    """src with every dimension, variable, attribute and text value, its numbers zero (compressed to almost nothing)."""
    with nc.Dataset(src) as a, nc.Dataset(dst, "w", format="NETCDF4") as b:
        for n, d in a.dimensions.items():
            b.createDimension(n, None if d.isunlimited() else len(d))
        b.setncatts({k: a.getncattr(k) for k in a.ncattrs()})
        for n, v in a.variables.items():
            w = b.createVariable(n, v.dtype, v.dimensions, zlib=True, complevel=1)
            w.setncatts({k: v.getncattr(k) for k in v.ncattrs()})
            w[:] = v[:] if v.dtype.kind in "SU" else np.zeros(v.shape, v.dtype)


def files_cases():
    """{name: (config, wrfout reference dir, restart dir or None, extra --set)} of run_wrf_files_test.sh's CASES."""
    text = open(os.path.join(REPO, "tests", "verify", "run_wrf_files_test.sh")).read()
    out = {}
    for line in text.split('CASES="', 1)[1].split('\n"', 1)[0].strip().splitlines():
        name, cfg, ref, rst, extra = line.split("|")
        out[name] = (cfg, ref.replace("$R", WR), None if rst == "-" else rst.replace("$R", WR), extra.replace("$R", WR))
    return out


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
    ap.add_argument("--files", default="")
    ap.add_argument("--no-thompson-tables", action="store_true")
    ap.add_argument("--whole", action="store_true")
    ap.add_argument("--refs-only", action="store_true")
    ap.add_argument("--no-steps", action="store_true")
    ap.add_argument("out")
    a = ap.parse_args()
    log = (lambda *x: print(*x, file=sys.stderr)) if a.out == "-" else print  # stdout may be the tar
    want = [o for o in a.only.split(",") if o]
    dirs = {}  # run directory -> history files to keep, Thompson tables needed
    ref_times = {}  # reference wrfout -> compared times
    exact = set()  # reference wrfout files of tolerance-0 cases
    ideal = set()  # files under CPPWRF_WRF_REF: ideal-case soundings and reference wrfout files
    if any("straka" in o for o in want):
        ideal.add(os.path.join(WREF, "run_straka_d1", "wrfout_d01_0001-01-01_00:00:00"))
    for name, key, cfg, sets, wo, t, tol in wrf_suite.CASES:
        if not any(o in name for o in want):
            continue
        snd = config_value(cfg, "case", "sounding_file")
        refs = list(wo.values()) if isinstance(wo, dict) else [wo]
        if snd and "${CPPWRF_WRF_REF}" in snd or any(r.startswith(WREF) for r in refs):
            if snd:
                ideal.add(snd.replace("${CPPWRF_WRF_REF}", WREF))
            ideal.update(r for r in refs if r.startswith(WREF))
            ideal.update(os.path.join(os.path.dirname(r), "namelist.input") for r in refs if r.startswith(WREF))
            log(f"{name}: ideal ({len(refs)} reference files)")
            continue
        rd = next((s.split("=", 1)[1] for s in sets if s.startswith("real.dir=")), None)
        rd = rd or config_value(cfg, "real", "dir")
        if not rd or "${CPPWRF_WRF_RUNS}" not in rd and not rd.startswith(WR):
            log(f"{name}: skipped (not a real-data case of CPPWRF_WRF_RUNS)")
            continue
        rd = rd.replace("${CPPWRF_WRF_RUNS}", WR)
        mp = next((s.split("=", 1)[1] for s in sets if s.startswith("phys.mp_physics=")), None)
        mp = mp or config_value(cfg, "phys", "mp_physics") or "0"
        refs = [r for r in (wo.values() if isinstance(wo, dict) else [wo]) if not (a.no_steps and "/steps/" in r)]
        for r in refs:
            d = dirs.setdefault(os.path.dirname(r), [set(), False])
            d[0].add(os.path.basename(r))
            ref_times.setdefault(r, set()).add(float(t))
            if all(x == 0.0 for x in tol.values()):
                exact.add(r)
        d = dirs.setdefault(rd, [set(), False]) if not a.refs_only else [set(), False]
        d[1] = d[1] or (mp.strip() == "8" and not a.no_thompson_tables)
        log(f"{name}: {os.path.relpath(rd, WR)}")
    for tag, cfg, _ in wrf_suite.RESTARTS:  # restart_<tag>_*: the input directory of the configuration
        if not a.refs_only and any(o in n for n in wrf_suite.SPECIAL if n.startswith(f"restart_{tag}_") for o in want):
            rd = config_value(cfg, "real", "dir").replace("${CPPWRF_WRF_RUNS}", WR)
            d = dirs.setdefault(rd, [set(), False])
            d[1] = d[1] or ((config_value(cfg, "phys", "mp_physics") or "0").strip() == "8" and not a.no_thompson_tables)
            log(f"restart_{tag}: {os.path.relpath(rd, WR)}")
    skel, whole = {}, {}  # files test references: skeletons, whole files (path in the tar -> file)
    fcases = files_cases()
    for name in [f for f in a.files.split(",") if f]:
        cfg, ref, rst, extra = fcases[name]
        m = re.search(r"real\.dir=(\S+)", extra)
        rd = m.group(1) if m else config_value(cfg, "real", "dir").replace("${CPPWRF_WRF_RUNS}", WR)
        d = dirs.setdefault(rd, [set(), False])
        d[1] = d[1] or ((config_value(cfg, "phys", "mp_physics") or "0").strip() == "8" and not a.no_thompson_tables)
        rsm = rst.replace("/rst_", "/rsm_") if rst else None
        for x in (ref, rst, rsm):
            if not x or not os.path.isdir(x):
                continue
            for f in sorted(os.listdir(x)):
                p, arc = os.path.join(x, f), os.path.join("wrf_runs", os.path.relpath(os.path.join(x, f), WR))
                if re.fullmatch(r"wrfrst_d\d\d_.*_00:01:00", f) and x == rst:
                    whole[arc] = p
                elif re.fullmatch(r"(wrfout|wrfrst)_d\d\d_.*", f):
                    if p in ref_times:
                        sys.exit(f"{p}: a reference of the suite cases and of the files test {name}: stage them apart")
                    skel.setdefault(arc, p)
                elif f in ("namelist.input", "myoutfields.txt"):
                    whole[arc] = p
        log(f"files {name}: {os.path.relpath(rd, WR)}, {os.path.relpath(ref, WR)}"
            + (f", {os.path.relpath(rst, WR)}, {os.path.relpath(rsm, WR)}" if rst else ""))
    tf = tarfile.open(fileobj=sys.stdout.buffer, mode="w|", dereference=True) if a.out == "-" else \
        tarfile.open(a.out, "w", dereference=True)
    added = set()
    with tf:
        for d, (keep, thompson) in sorted(dirs.items()):
            for f in sorted(os.listdir(d)):
                p = os.path.join(d, f)
                if os.path.isdir(p) or (f in THOMPSON and not thompson) or \
                        (a.refs_only and p not in ref_times and f not in ("namelist.input", "myoutfields.txt")):
                    continue
                arc = os.path.join("wrf_runs", os.path.relpath(p, WR))
                added.add(arc)
                if p in ref_times and a.whole:
                    tf.add(p, arcname=arc, recursive=False)
                elif p in ref_times:
                    with tempfile.TemporaryDirectory() as tmp:
                        q = os.path.join(tmp, f)
                        trimmed(p, ref_times[p], q, p in exact)
                        tf.add(q, arcname=arc)
                elif f in THOMPSON or not any(fnmatch.fnmatch(f, s) for s in SKIP):
                    tf.add(p, arcname=arc, recursive=False)
                else:
                    added.discard(arc)
        for arc, p in sorted(whole.items()):
            if arc not in added:
                tf.add(p, arcname=arc, recursive=False)
        for arc, p in sorted(skel.items()):
            if arc in whole:
                continue
            with tempfile.TemporaryDirectory() as tmp:
                q = os.path.join(tmp, os.path.basename(p))
                skeleton(p, q)
                tf.add(q, arcname=arc)
        for p in sorted(ideal):
            if os.path.exists(p):
                tf.add(p, arcname=os.path.join("wrf_ref", os.path.relpath(p, WREF)), recursive=False)
            else:
                log(f"missing: {p}")
    if a.out != "-":
        print(f"{a.out}: {os.path.getsize(a.out) / 2**20:.0f} MiB, {len(dirs)} directories")


if __name__ == "__main__":
    main()
