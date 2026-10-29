#!/usr/bin/env python3
"""wrf_ref_run.py RUNS DEST WRFEXE NAME BASE SECTION.key=value ... [--prepare-only]

A WRF reference run in the layout of the suite references (tests/verify/wrf_suite.py): DEST/run_NAME with the
inputs and tables of RUNS/BASE (links), wrf.exe -> WRFEXE and BASE's namelist.input with the given changes (key
replaced or added in &SECTION; run_* and history_interval_s from time_control.run_seconds /
time_control.history_interval_s, one history file), and DEST/run_NAME/steps: the same namelist for 3 time steps
with history every step. Both run wrf.exe (in parallel) unless --prepare-only. Special keys: run.real=EXE:DIR,...
makes the inputs first: real.exe (EXE, next to WRFEXE) runs in run_NAME on the met_em files of RUNS/DIR,... (the
first directory holding a file name gives it) with the same namelist, and both runs read its wrfinput / wrfbdy
instead of BASE's; run.steps=0 leaves the 3-step run out (long references). run.segment=S runs wrf.exe in parts of
S seconds, each resumed from the previous part's wrfrst files (WRF's own restart; for runs longer than a spot VM
lasts), and merges each domain's history files into the first part's; with the environment WRFREF_PARTS (a bucket
prefix) the run directory goes there after each part and a new run continues after the last part found there.
Prints one line per run: NAME[/steps] SUCCESS | WRF-FAILED [error lines] (NAME REAL-FAILED [error lines] when
real.exe fails)."""
import datetime
import glob
import os
import re
import subprocess
import sys

args = [a for a in sys.argv[1:] if a != "--prepare-only"]
prepare_only = len(args) != len(sys.argv) - 1
runs, dest, wrfexe, name, base = args[:5]
sets = [a.split("=", 1) for a in args[5:]]
special = {k: v for k, v in sets if k in ("run.real", "run.steps", "run.segment")}
sets = [(k, v) for k, v in sets if k not in special]
src = os.path.join(runs, base)
INPUTS = r"(wrfinput.*|wrfbdy.*|wrflowinp.*)$"


def setkey(nl, sec, key, val):  # as tools/gcp/matrix_variant.py (keeping the line breaks)
    m = re.search(rf"&{sec}\b(.*?)\n\s*/", nl, re.S | re.I)
    if not m:
        return nl.rstrip("\n") + f"\n\n&{sec}\n {key} = {val},\n/\n"
    body = m.group(1)
    pat = re.compile(rf"^[ \t]*{key}\s*=.*$", re.M | re.I)
    line = f" {key} = {val},"
    body = pat.sub(line, body) if pat.search(body) else body + "\n" + line
    return nl[:m.start(1)] + body + nl[m.end(1):]


def wrf_failed(d, logname="wrf.log"):
    """None when wrf.exe ended with SUCCESS in d, else its error lines."""
    log = open(os.path.join(d, logname), errors="replace").read()
    if "SUCCESS COMPLETE WRF" in log:
        return None
    lines = log.strip().splitlines()
    return [l.strip() for l in lines if 'ERROR' in l or 'FATAL' in l][:4] or lines[-1:]


def nlint(nl, key, default=None):
    """The (first) integer value of a namelist key; `default` when the namelist does not set it (WRF's default)."""
    m = re.search(rf"^\s*{key}\s*=\s*(\d+)", nl, re.M | re.I)
    if m is None and default is None:
        raise KeyError(f"namelist.input sets no {key}")
    return int(m.group(1)) if m else default


def segmented(R, nl, seg, parts_prefix):
    """wrf.exe in R in parts of `seg` seconds (run.segment), the wrfrst files of the namelist's restart_interval
    kept; then each domain's history files appended to the first part's file. None or the error lines."""
    total, ri, ndom = nlint(nl, "run_seconds"), nlint(nl, "restart_interval") * 60, nlint(nl, "max_dom", 1)
    if total % seg or seg % 60:
        return [f"run_seconds {total} is not a multiple of run.segment {seg} (whole minutes)"]
    t0 = datetime.datetime(*(nlint(nl, f"start_{k}", 0 if k in ("minute", "second") else None)
                             for k in ("year", "month", "day", "hour", "minute", "second")))
    tag = lambda t: t.strftime("%Y-%m-%d_%H:%M:%S")
    first = 0
    if parts_prefix:  # the last part a previous VM kept
        ls = subprocess.run(["gcloud", "storage", "ls", f"{parts_prefix}-part*.tar"], capture_output=True, text=True)
        done = sorted(int(m.group(1)) for m in re.finditer(r"-part(\d+)\.tar", ls.stdout))
        if done:
            subprocess.run(f"gcloud storage cat {parts_prefix}-part{done[-1]}.tar | tar xf - -C {os.path.dirname(R)}",
                           shell=True, check=True)
            first = done[-1]
    for k in range(first, total // seg):
        t = t0 + datetime.timedelta(seconds=k * seg)
        n = nl
        for key, v in [("run_seconds", seg), ("restart", ".true." if k else ".false."), ("restart_interval", seg // 60)] \
                + [(f"start_{u}", getattr(t, u)) for u in ("year", "month", "day", "hour", "minute", "second")]:
            n = setkey(n, "time_control", key, ", ".join([str(v)] * (1 if key in ("run_seconds", "restart",
                                                                                  "restart_interval") else ndom)))
        open(os.path.join(R, "namelist.input"), "w").write(n)
        subprocess.run(f"./wrf.exe > wrf_part{k}.log 2>&1", shell=True, cwd=R, env=env)
        bad = wrf_failed(R, f"wrf_part{k}.log")
        if bad:
            return [f"part {k}"] + bad
        if k:  # the wrfrst files this part started from, unless at a restart time of the reference's namelist
            for f in glob.glob(os.path.join(R, f"wrfrst_d??_{tag(t)}")):
                if (k * seg) % ri:
                    os.remove(f)
        if parts_prefix:
            subprocess.run(f"tar cf - -C {os.path.dirname(R)} {os.path.basename(R)} | "
                           f"gcloud storage cp - {parts_prefix}-part{k + 1}.tar -q > /dev/null 2>&1", shell=True)
    for f in glob.glob(os.path.join(R, "wrfrst_d??_*")):  # the last part's restart_interval may differ
        tf = datetime.datetime.strptime(f[-19:], "%Y-%m-%d_%H:%M:%S")
        if (tf - t0).total_seconds() % ri:
            os.remove(f)
    import netCDF4
    for d in range(1, ndom + 1):  # one history file per domain, as an uninterrupted run with frames_per_outfile
        files = sorted(glob.glob(os.path.join(R, f"wrfout_d{d:02d}_*")))
        main = os.path.join(R, f"wrfout_d{d:02d}_{tag(t0)}")
        with netCDF4.Dataset(main, "a") as a:
            a.set_auto_maskandscale(False)
            for f in files:
                if f == main:
                    continue
                with netCDF4.Dataset(f) as b:
                    b.set_auto_maskandscale(False)
                    n0, m = len(a.dimensions["Time"]), len(b.dimensions["Time"])
                    for v in b.variables:
                        if b[v].dimensions[:1] == ("Time",):
                            a[v][n0:n0 + m] = b[v][:]
                os.remove(f)
    open(os.path.join(R, "namelist.input"), "w").write(nl)  # the reference's namelist (suite configuration check)
    return None


def prepare(R, nl, inputs=src):
    """R with BASE's tables and the input files of the directory `inputs` (none: real.exe makes them)."""
    subprocess.run(["rm", "-rf", R])
    os.makedirs(R)
    for f in os.listdir(src):
        p = os.path.join(src, f)
        if not os.path.isdir(p) and not re.match(r"(wrfout.*|wrfrst.*|rsl\..*|namelist\.(output|input)|.*\.log|.*\.exe)$", f) \
                and not re.match(INPUTS, f):
            os.symlink(os.path.realpath(p), os.path.join(R, f))
    for f in os.listdir(inputs) if inputs else []:
        if re.match(INPUTS, f):
            os.symlink(os.path.realpath(os.path.join(inputs, f)), os.path.join(R, f))
    os.symlink(os.path.realpath(wrfexe), os.path.join(R, "wrf.exe"))
    open(os.path.join(R, "namelist.input"), "w").write(nl)


nl = open(os.path.join(src, "namelist.input")).read()
for k, v in [("run_days", "0"), ("run_hours", "0"), ("run_minutes", "0"), ("frames_per_outfile", "10000"),
             ("restart", ".false.")]:
    nl = setkey(nl, "time_control", k, v)
nl = re.sub(r"^\s*history_interval\s*=.*\n", "", nl, flags=re.M)
for sk, v in sets:
    sec, k = sk.split(".", 1)
    nl = setkey(nl, sec, k, v)
dt = int(re.search(r"^\s*time_step\s*=\s*(\d+)", nl, re.M).group(1))
steps = nl
for k, v in [("run_seconds", str(3 * dt)), ("history_interval_s", f"{dt}, {dt}, {dt}")]:
    steps = setkey(steps, "time_control", k, v)
R = os.path.join(dest, f"run_{name}")
env = dict(os.environ, OMP_NUM_THREADS="1")
real = special.get("run.real")
prepare(R, nl, None if real else src)
if real:  # real.exe on the met_em files of the listed directories, in R with the reference's namelist
    exe, dirs = real.split(":", 1)
    os.symlink(os.path.realpath(os.path.join(os.path.dirname(wrfexe), exe)), os.path.join(R, "real.exe"))
    for d in dirs.split(","):
        for f in sorted(os.listdir(os.path.join(runs, d))):
            if f.startswith("met_em") and not os.path.lexists(os.path.join(R, f)):
                os.symlink(os.path.realpath(os.path.join(runs, d, f)), os.path.join(R, f))
    if not prepare_only:
        subprocess.run("./real.exe > real.log 2>&1", shell=True, cwd=R, env=env)
        log = "".join(open(os.path.join(R, f), errors="replace").read()
                      for f in ("real.log", "rsl.out.0000", "rsl.error.0000") if os.path.exists(os.path.join(R, f)))
        if "SUCCESS COMPLETE REAL_EM" not in log:
            lines = log.strip().splitlines()
            print(f"{name} REAL-FAILED {[l.strip() for l in lines if 'ERROR' in l or 'FATAL' in l][:4] or lines[-1:]}",
                  flush=True)
            sys.exit(0)
runs_now = [R]
if special.get("run.steps", "1") != "0":
    prepare(os.path.join(R, "steps"), steps, R if real else src)
    runs_now.append(os.path.join(R, "steps"))
if prepare_only:
    sys.exit(0)
seg = int(special.get("run.segment", "0"))
jobs = {d: subprocess.Popen("./wrf.exe > wrf.log 2>&1", shell=True, cwd=d, env=env) for d in runs_now
        if not (seg and d == R)}
if seg:
    bad = segmented(R, nl, seg, os.environ.get("WRFREF_PARTS"))
    print(f"{name} SUCCESS" if not bad else f"{name} WRF-FAILED {bad}", flush=True)
for d, p in jobs.items():
    p.wait()
    bad = wrf_failed(d)
    tag = name if d == R else f"{name}/steps"
    print(f"{tag} SUCCESS" if not bad else f"{tag} WRF-FAILED {bad}", flush=True)
