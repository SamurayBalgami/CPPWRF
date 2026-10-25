#!/usr/bin/env python3
"""matrix_variant.py RUNS WORK WRFEXE CPPWRF NAME BASE SECTION.key=value ...

One WRF-vs-CPPWRF comparison of a namelist variant: the WRF run RUNS/BASE (inputs and namelist.input) runs
3 minutes with history every minute and the given namelist changes (key replaced or added in &SECTION), in
WORK/NAME; CPPWRF runs the configuration tools/wrf2cppwrf.py makes from that namelist; every common variable
of every common frame (by Times) of each domain's wrfout is compared at tolerance 0. Prints one line:
NAME d0N frames N differing: none | {VAR: 'time npoints/maxdiff', ...}  (or NAME WRF-FAILED / CPPWRF-FAILED);
variables only WRF writes are listed as wrf-only-vars."""
import glob, os, re, subprocess, sys
import netCDF4
import numpy as np

runs, work, wrfexe, cppwrf, name, base = sys.argv[1:7]
sets = [a.split("=", 1) for a in sys.argv[7:]]
repo = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
R = os.path.join(work, name)
subprocess.run(["rm", "-rf", R])
os.makedirs(R)
src = os.path.join(runs, base)
for f in os.listdir(src):
    if not re.match(r"(wrfout.*|wrfrst.*|rsl\..*|namelist\.(output|input)|.*\.log|.*\.exe)$", f):
        os.symlink(os.path.realpath(os.path.join(src, f)), os.path.join(R, f))
os.symlink(wrfexe, os.path.join(R, "wrf.exe"))
nl = open(os.path.join(src, "namelist.input")).read()


def setkey(nl, sec, key, val):
    m = re.search(rf"&{sec}\b(.*?)\n\s*/", nl, re.S | re.I)
    if not m:  # a group the base namelist lacks (e.g. &stoch)
        return nl.rstrip("\n") + f"\n\n&{sec}\n {key} = {val},\n/\n"
    body = m.group(1)
    pat = re.compile(rf"^\s*{key}\s*=.*$", re.M | re.I)
    line = f" {key} = {val},"
    body = pat.sub(line, body) if pat.search(body) else body + "\n" + line
    return nl[:m.start(1)] + body + nl[m.end(1):]


for k, v in [("run_days", "0"), ("run_hours", "0"), ("run_minutes", "0"), ("run_seconds", "180"),
             ("history_interval_s", "60, 60, 60"), ("frames_per_outfile", "10000"), ("restart", ".false.")]:
    nl = setkey(nl, "time_control", k, v)
nl = re.sub(r"^\s*history_interval\s*=.*\n", "", nl, flags=re.M)
for sk, v in sets:
    sec, k = sk.split(".", 1)
    if (sec, k) == ("run", "remove"):  # run.remove=f1,f2: files of the base run left out (both models)
        for f in v.split(","):
            if os.path.lexists(os.path.join(R, f)):
                os.remove(os.path.join(R, f))
        continue
    nl = setkey(nl, sec, k, v)
open(os.path.join(R, "namelist.input"), "w").write(nl)
env = dict(os.environ, OMP_NUM_THREADS="1")
# ideal cases: WRF_EXE wrf_X.exe with an ideal_X.exe next to it -> ideal.exe makes wrfinput_d01 from the
# variant's namelist first (boundary and grid options act in ideal.exe too)
iexe = os.path.join(os.path.dirname(wrfexe), "ideal_" + os.path.basename(wrfexe)[len("wrf_"):])
if os.path.basename(wrfexe).startswith("wrf_") and os.path.exists(iexe):
    for f in glob.glob(os.path.join(R, "wrfinput_d0*")):
        os.remove(f)
    os.symlink(iexe, os.path.join(R, "ideal.exe"))
    subprocess.run("./ideal.exe > ideal.log 2>&1", shell=True, cwd=R, env=env)
    if not os.path.exists(os.path.join(R, "wrfinput_d01")):
        print(f"{name} IDEAL-FAILED {open(os.path.join(R, 'ideal.log'), errors='replace').read().strip().splitlines()[-1:]}")
        sys.exit(0)
try:  # a hanging run must not keep the VM
    subprocess.run("./wrf.exe > wrf.log 2>&1", shell=True, cwd=R, env=env, timeout=3600)
except subprocess.TimeoutExpired:
    open(os.path.join(R, "wrf.log"), "a").write("\nTIMEOUT after 3600 s\n")
if "SUCCESS COMPLETE WRF" not in open(os.path.join(R, "rsl.out.0000") if os.path.exists(os.path.join(R, "rsl.out.0000")) else os.path.join(R, "wrf.log"), errors="replace").read():
    log = open(os.path.join(R, 'wrf.log'), errors='replace').read().strip().splitlines()
    why = [l.strip() for l in log if "ERROR" in l or "FATAL" in l or "Fix " in l][:4] or log[-1:]
    print(f"{name} WRF-FAILED {why}")
    sys.exit(0)
out = os.path.join(R, "cppwrf")
t = subprocess.run([sys.executable, os.path.join(repo, "tools/wrf2cppwrf.py"), os.path.join(R, "namelist.input"),
                    "--real-dir", R, "--format", "wrfout", "-o", os.path.join(R, "case.ini")], capture_output=True, text=True)
if t.returncode != 0:
    print(f"{name} TRANSLATE-FAILED {t.stderr.strip().splitlines()[-1:]}")
    sys.exit(0)
try:
    c = subprocess.run([cppwrf, os.path.join(R, "case.ini"), "--set", f"output.dir={out}", "--set", "run.stop_on_nonfinite=false"], capture_output=True, text=True,
                       cwd=repo, env=env, timeout=3600)
except subprocess.TimeoutExpired:
    c = subprocess.CompletedProcess([], 124, "", "TIMEOUT after 3600 s")
if c.returncode != 0:
    print(f"{name} CPPWRF-FAILED {(c.stderr or c.stdout).strip().splitlines()[-1:]}")
    sys.exit(0)
unread = re.findall(r"unaccessed config key '([^']*)'", c.stderr)
if unread:  # a translated option the model never read
    print(f"{name} UNREAD-KEYS {unread}", flush=True)


def load(files):
    m = {}
    for f in sorted(files):
        ds = netCDF4.Dataset(f)
        for t in range(len(ds["Times"])):
            m[ds["Times"][t].tobytes().decode()] = (ds, t)
    return m


for kind in ("wrfout", "wrfrst"):  # history frames, then restart files (both streams' common variables)
    for d in ("d01", "d02", "d03"):
        w = glob.glob(os.path.join(R, f"{kind}_{d}_*"))
        cdir = os.path.join(out, d) if os.path.isdir(os.path.join(out, d)) else out
        cw = glob.glob(os.path.join(cdir, f"{kind}_{d}_*"))
        if not w and not cw:
            continue
        A, B = load(cw), load(w)
        bad, nframes, wrf_only = {}, {}, set()
        for k in sorted(set(A) & set(B)):
            (a, ta), (b, tb) = A[k], B[k]
            for v in b.variables:
                if v not in a.variables and b[v].dtype.kind != "S":
                    wrf_only.add(v)
                if v not in a.variables or b[v].dtype.kind == "S" or "Time" not in b[v].dimensions:
                    continue
                x = np.atleast_1d(np.asarray(a[v][ta], float))
                y = np.atleast_1d(np.asarray(b[v][tb], float))
                if x.shape != y.shape:
                    if v not in bad:  # WRF's values per leading index (e.g. soil layers)
                        yl = y.reshape(y.shape[0], -1) if y.ndim > 1 else y.reshape(-1, 1)
                        bad[v] = f"shape cppwrf {x.shape} wrf {y.shape} wrf min/max per level " + \
                                 " ".join(f"{r.min():.4g}/{r.max():.4g}" for r in yl)
                    continue
                dd = np.abs(x - y)
                dd[np.isnan(x) & np.isnan(y)] = 0
                if dd.max() > 0:
                    nframes[v] = nframes.get(v, 0) + 1
                if dd.max() > 0 and v not in bad:
                    im = int(np.nanargmax(dd))
                    bad[v] = f"{k} {(dd > 0).sum()}/{dd.max():.2e} (cppwrf {x.flat[im]:.6g} wrf {y.flat[im]:.6g} at {im})"
        only = ""
        if set(A) != set(B):  # times written by one model only
            only = f" only-cppwrf {sorted(set(A) - set(B))} only-wrf {sorted(set(B) - set(A))}"
        if wrf_only:  # variables WRF writes and CPPWRF does not
            only += f" wrf-only-vars {sorted(wrf_only)}"
        for v in bad:
            if nframes.get(v, 0) > 1:
                bad[v] += f" [{nframes[v]} frames]"
        tag = d if kind == "wrfout" else f"{d} {kind}"
        print(f"{name} {tag} frames {len(set(A) & set(B))}{only} differing: {bad if bad else 'none'}", flush=True)
subprocess.run(["rm", "-rf", out] + glob.glob(os.path.join(R, "wrfout_*")) + glob.glob(os.path.join(R, "wrfrst_*")))
