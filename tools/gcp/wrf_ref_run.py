#!/usr/bin/env python3
"""wrf_ref_run.py RUNS DEST WRFEXE NAME BASE SECTION.key=value ... [--prepare-only]

A WRF reference run in the layout of the suite references (tests/verify/wrf_suite.py): DEST/run_NAME with the
inputs and tables of RUNS/BASE (links), wrf.exe -> WRFEXE and BASE's namelist.input with the given changes (key
replaced or added in &SECTION; run_* and history_interval_s from time_control.run_seconds /
time_control.history_interval_s, one history file), and DEST/run_NAME/steps: the same namelist for 3 time steps
with history every step. Both run wrf.exe (in parallel) unless --prepare-only. Prints one line per run:
NAME[/steps] SUCCESS | WRF-FAILED [error lines]."""
import os
import re
import subprocess
import sys

args = [a for a in sys.argv[1:] if a != "--prepare-only"]
prepare_only = len(args) != len(sys.argv) - 1
runs, dest, wrfexe, name, base = args[:5]
sets = [a.split("=", 1) for a in args[5:]]
src = os.path.join(runs, base)


def setkey(nl, sec, key, val):  # as tools/gcp/matrix_variant.py (keeping the line breaks)
    m = re.search(rf"&{sec}\b(.*?)\n\s*/", nl, re.S | re.I)
    if not m:
        return nl.rstrip("\n") + f"\n\n&{sec}\n {key} = {val},\n/\n"
    body = m.group(1)
    pat = re.compile(rf"^[ \t]*{key}\s*=.*$", re.M | re.I)
    line = f" {key} = {val},"
    body = pat.sub(line, body) if pat.search(body) else body + "\n" + line
    return nl[:m.start(1)] + body + nl[m.end(1):]


def prepare(R, nl):
    subprocess.run(["rm", "-rf", R])
    os.makedirs(R)
    for f in os.listdir(src):
        p = os.path.join(src, f)
        if not os.path.isdir(p) and not re.match(r"(wrfout.*|wrfrst.*|rsl\..*|namelist\.(output|input)|.*\.log|.*\.exe)$", f):
            os.symlink(os.path.realpath(p), os.path.join(R, f))
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
prepare(R, nl)
prepare(os.path.join(R, "steps"), steps)
if prepare_only:
    sys.exit(0)
env = dict(os.environ, OMP_NUM_THREADS="1")
jobs = {d: subprocess.Popen("./wrf.exe > wrf.log 2>&1", shell=True, cwd=d, env=env)
        for d in (R, os.path.join(R, "steps"))}
for d, p in jobs.items():
    p.wait()
    log = open(os.path.join(d, "wrf.log"), errors="replace").read()
    tag = name if d == R else f"{name}/steps"
    if "SUCCESS COMPLETE WRF" in log:
        print(f"{tag} SUCCESS", flush=True)
    else:
        lines = log.strip().splitlines()
        print(f"{tag} WRF-FAILED {[l.strip() for l in lines if 'ERROR' in l or 'FATAL' in l][:4] or lines[-1:]}", flush=True)
