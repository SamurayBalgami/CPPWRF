#!/usr/bin/env python3
"""real_variant.py RUNS WORK REALEXE CPPWRF_REAL NAME BASE SECTION.key=value ...

One WRF real.exe vs cppwrf_real comparison of a namelist variant (the real-data counterpart of matrix_variant.py):
RUNS/BASE holds met_em files, namelist.input and `precision` (float / double, the WRF build's real kind). Both
programs run in WORK/NAME with the changed namelist; every file WRF writes (wrfinput_d0N, wrfbdy_d01,
wrflowinp_d0N, ...) is compared with tests/init_real/compare_real.py (values at tolerance 0, variables, dimensions,
attributes). Special keys: run.remove=FILE,... leaves base files out; run.delete=section.key,... deletes namelist
lines (the key then takes its default). Prints one line per variant:
NAME real files N differing: none | {FILE: 'n variables / m metadata', ...}  (or REAL-FAILED / CPPWRF-REAL-FAILED,
with the files only one program wrote)."""
import glob, os, re, subprocess, sys

runs, work, realexe, cppreal, name, base = sys.argv[1:7]
sets = [a.split("=", 1) for a in sys.argv[7:]]
repo = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
R = os.path.join(work, name)
subprocess.run(["rm", "-rf", R])
os.makedirs(R)
src = os.path.join(runs, base)
for f in os.listdir(src):
    if not re.match(r"(wrfinput.*|wrfbdy.*|wrflowinp.*|rsl\..*|namelist\.(output|input)|.*\.log|.*\.exe)$", f):
        os.symlink(os.path.realpath(os.path.join(src, f)), os.path.join(R, f))
os.symlink(realexe, os.path.join(R, "real.exe"))
nl = open(os.path.join(src, "namelist.input")).read()


def setkey(nl, sec, key, val):
    m = re.search(rf"&{sec}\b(.*?)\n\s*/", nl, re.S | re.I)
    if not m:
        raise SystemExit(f"{name} BAD-NAMELIST no &{sec}")
    body = m.group(1)
    pat = re.compile(rf"^\s*{key}\s*=.*$", re.M | re.I)
    line = f" {key} = {val},"
    body = pat.sub(line, body) if pat.search(body) else body + "\n" + line
    return nl[:m.start(1)] + body + nl[m.end(1):]


for sk, v in sets:
    sec, k = sk.split(".", 1)
    if (sec, k) == ("run", "remove"):
        for f in v.split(","):
            if os.path.lexists(os.path.join(R, f)):
                os.remove(os.path.join(R, f))
    elif (sec, k) == ("run", "delete"):
        for sk2 in v.split(","):
            s2, k2 = sk2.split(".", 1)
            m = re.search(rf"&{s2}\b(.*?)\n\s*/", nl, re.S | re.I)
            if m:
                body = re.sub(rf"^\s*{k2}\s*=.*\n?", "", m.group(1), flags=re.M | re.I)
                nl = nl[:m.start(1)] + body + nl[m.end(1):]
    else:
        nl = setkey(nl, sec, k, v)
open(os.path.join(R, "namelist.input"), "w").write(nl)

env = dict(os.environ, OMP_NUM_THREADS="1")
subprocess.run("./real.exe > real.log 2>&1", shell=True, cwd=R, env=env)
log = open(os.path.join(R, "rsl.out.0000") if os.path.exists(os.path.join(R, "rsl.out.0000"))
           else os.path.join(R, "real.log"), errors="replace").read()
wrf_ok = "SUCCESS COMPLETE REAL_EM" in log
wrf_files = sorted(os.path.basename(f) for f in glob.glob(os.path.join(R, "wrf*_d0*")))
out = os.path.join(R, "cppwrf")
os.makedirs(out)
prec = open(os.path.join(src, "precision")).read().strip() if os.path.exists(os.path.join(src, "precision")) else "double"
c = subprocess.run([cppreal, "--namelist", os.path.join(R, "namelist.input"), "--met-dir", R, "--out-dir", out,
                    "--precision", prec, "--schema", os.path.join(repo, "data/wrf_schema/wrf-4.8.0.tsv")],
                   capture_output=True, text=True, env=env)
cpp_ok = c.returncode == 0
cpp_files = sorted(os.listdir(out)) if cpp_ok else []
if not wrf_ok:
    why = [l.strip() for l in log.splitlines() if "ERROR" in l or "FATAL" in l or "Fix " in l][:4] or log.splitlines()[-1:]
    status = "also fails" if not cpp_ok else f"runs and writes {cpp_files}"
    print(f"{name} REAL-FAILED {why}; cppwrf_real {status}", flush=True)
elif not cpp_ok:
    print(f"{name} CPPWRF-REAL-FAILED {(c.stderr or c.stdout).strip().splitlines()[-1:]}", flush=True)
else:
    bad = {}
    for f in wrf_files:
        if f not in cpp_files:
            bad[f] = "not written by cppwrf_real"
            continue
        r = subprocess.run([sys.executable, os.path.join(repo, "tests/init_real/compare_real.py"),
                            os.path.join(out, f), os.path.join(R, f)], capture_output=True, text=True)
        if r.returncode != 0:
            lines = r.stdout.strip().splitlines()
            bad[f] = "; ".join(lines[:6]) + (f" ... ({len(lines)} lines)" if len(lines) > 6 else "")
    for f in cpp_files:
        if f not in wrf_files:
            bad[f] = "not written by real.exe"
    print(f"{name} real files {len(wrf_files)} differing: {bad if bad else 'none'}", flush=True)
subprocess.run(["rm", "-rf", out] + glob.glob(os.path.join(R, "wrf*_d0*")))
