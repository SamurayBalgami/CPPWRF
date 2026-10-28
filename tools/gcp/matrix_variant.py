#!/usr/bin/env python3
"""matrix_variant.py RUNS WORK WRFEXE CPPWRF NAME BASE SECTION.key=value ...

One WRF-vs-CPPWRF comparison of a namelist variant: the WRF run RUNS/BASE (inputs and namelist.input) runs
3 minutes with history every minute and the given namelist changes (key replaced or added in &SECTION), in
WORK/NAME; CPPWRF runs the configuration tools/wrf2cppwrf.py makes from that namelist; every variable of every frame (by
Times) of each domain's wrfout and wrfrst files is compared at tolerance 0 (tests/verify/exact.py: NaN equals only
NaN, -0 differs from +0), with the types, dimensions, variable and global attributes (TITLE aside). Prints one
line per domain and stream: NAME d0N [wrfrst] frames N differing: none | {VAR: 'time npoints/maxdiff ...', ...};
times, variables or metadata only one model has are differences too ('(times)', '(wrf-only-vars)',
'(cppwrf-only-vars)', '(metadata)'). Or NAME WRF-FAILED with CPPWRF's status on the same namelist,
TRANSLATE-FAILED / CPPWRF-FAILED / IDEAL-FAILED / NO-OUTPUT / HARNESS-ERROR (the script itself failed).
run.keep=1 copies WRF's namelist.input, inputs and wrfout files to WORK/keep/run_NAME (a suite reference);
run.delete=section.key,... deletes namelist lines (the option then takes its Registry default), as real_variant.py.
Input files of the base run read by both models can be edited: run.sub=FILE:OLD:NEW replaces every OLD by NEW
(percent-encoded: %20 a blank, %0A a newline), run.head=FILE:N keeps the first N lines, run.copy=FILE:OTHER
takes the base run's file OTHER. CPPWRF always runs when WRF fails (NAME WRF-FAILED [...]; cppwrf: RUNS or its failure: does CPPWRF reject
the input as WRF does); run.both=1 is accepted for older variant files.
run.atts=1 is accepted for older variant files (the attributes are always compared).
run.file=NAME:TEXT writes the file NAME (e.g. an iofields_filename list) into the run directory.
run.threads=N runs CPPWRF with OMP_NUM_THREADS=N (default 1; WRF always runs serially). Values are
percent-decoded keys and values (%20 blank, %0A new line, %27 quote, %23 #, %7C |), as the variant line passes through sh -c."""
import glob, os, re, subprocess, sys
import urllib.parse
import numpy as np


def harness_error(t, e, tb):  # a variant that cannot run or be compared is a result line, never a missing one
    print(f"{sys.argv[5] if len(sys.argv) > 5 else '?'} HARNESS-ERROR {t.__name__}: {e}", flush=True)
    sys.__excepthook__(t, e, tb)


sys.excepthook = harness_error
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tests", "verify"))
import exact  # noqa: E402

runs, work, wrfexe, cppwrf, name, base = sys.argv[1:7]
sets = [(urllib.parse.unquote(k), urllib.parse.unquote(v)) for k, v in (a.split("=", 1) for a in sys.argv[7:])]
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
keep = False
threads = "1"
for sk, v in sets:
    sec, k = sk.split(".", 1)
    if (sec, k) == ("run", "keep"):  # run.keep=1: WRF's run (namelist, inputs, wrfout) kept as a reference
        keep = v not in ("0", "false", ".false.")
        continue
    if (sec, k) == ("run", "threads"):  # run.threads=N: CPPWRF's OpenMP threads (the result must not change)
        threads = str(int(v))
        continue
    if (sec, k) == ("run", "both"):  # always on (see above)
        continue
    if (sec, k) in (("run", "sub"), ("run", "head"), ("run", "copy")):  # an edited copy of an input file
        f, arg = v.split(":", 1)
        p = os.path.join(R, f)
        text = open(os.path.join(src, arg) if k == "copy" else p, encoding="latin-1", newline="").read()
        if k == "head":
            text = "".join(text.splitlines(True)[:int(arg)])
        elif k == "sub":
            old, new = arg.split(":", 1)
            if old not in text:
                print(f"{name} BAD-VARIANT {old!r} not in {f}")
                sys.exit(0)
            text = text.replace(old, new)
        if os.path.lexists(p):
            os.remove(p)
        open(p, "w", encoding="latin-1", newline="").write(text)
        continue
    if (sec, k) == ("run", "atts"):  # always on (see above)
        continue
    if (sec, k) == ("run", "file"):  # run.file=NAME:TEXT: a file of the run directory (both models)
        fn, text = v.split(":", 1)
        if os.path.lexists(os.path.join(R, fn)):  # never write through a link into the base run
            os.remove(os.path.join(R, fn))
        open(os.path.join(R, fn), "w").write(text)
        continue
    if (sec, k) == ("run", "remove"):  # run.remove=f1,f2: files of the base run left out (both models)
        for f in v.split(","):
            if os.path.lexists(os.path.join(R, f)):
                os.remove(os.path.join(R, f))
        continue
    if (sec, k) == ("run", "delete"):  # run.delete=section.key,...: namelist lines left out (Registry default)
        for sk2 in v.split(","):
            s2, k2 = sk2.split(".", 1)
            m = re.search(rf"&{s2}\b(.*?)\n\s*/", nl, re.S | re.I)
            if m:
                body = re.sub(rf"^\s*{k2}\s*=.*\n?", "", m.group(1), flags=re.M | re.I)
                nl = nl[:m.start(1)] + body + nl[m.end(1):]
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
out = os.path.join(R, "cppwrf")


def run_cppwrf():
    """Translate the variant's namelist and run CPPWRF: (None, CompletedProcess) or (failure text, None)."""
    t = subprocess.run([sys.executable, os.path.join(repo, "tools/wrf2cppwrf.py"), os.path.join(R, "namelist.input"),
                        "--real-dir", R, "--format", "wrfout", "-o", os.path.join(R, "case.ini")], capture_output=True, text=True)
    if t.returncode != 0:
        return f"TRANSLATE-FAILED {t.stderr.strip().splitlines()[-1:]}", None
    try:
        c = subprocess.run([cppwrf, os.path.join(R, "case.ini"), "--set", f"output.dir={out}", "--set", "run.stop_on_nonfinite=false"], capture_output=True, text=True,
                           cwd=repo, env=dict(env, OMP_NUM_THREADS=threads), timeout=3600)
    except subprocess.TimeoutExpired:
        c = subprocess.CompletedProcess([], 124, "", "TIMEOUT after 3600 s")
    if c.returncode != 0:
        return f"CPPWRF-FAILED {(c.stderr or c.stdout).strip().splitlines()[-1:]}", None
    return None, c


if "SUCCESS COMPLETE WRF" not in open(os.path.join(R, "rsl.out.0000") if os.path.exists(os.path.join(R, "rsl.out.0000")) else os.path.join(R, "wrf.log"), errors="replace").read():
    log = open(os.path.join(R, 'wrf.log'), errors='replace').read().strip().splitlines()
    why = [l.strip() for l in log if "ERROR" in l or "FATAL" in l or "Fix " in l][:4] or log[-1:]
    fail, _ = run_cppwrf()  # CPPWRF must reject what WRF rejects
    print(f"{name} WRF-FAILED {why}; cppwrf: {fail if fail else 'RUNS'}")
    subprocess.run(["rm", "-rf", out])
    sys.exit(0)
if keep:  # WORK/keep/run_NAME: what the suite needs of the reference (remote_matrix.sh uploads keep.tar)
    K = os.path.join(work, "keep", f"run_{name}")
    os.makedirs(K, exist_ok=True)
    for f in ["namelist.input", "wrfinput_d01", "input_sounding"] + sorted(glob.glob(os.path.join(R, "wrfout_d0*"))):
        f = os.path.join(R, f)
        if os.path.exists(f):
            subprocess.run(["cp", "-L", f, K])
fail, c = run_cppwrf()
if fail:
    print(f"{name} {fail}")
    sys.exit(0)
unread = re.findall(r"unaccessed config key '([^']*)'", c.stderr)
if unread:  # a translated option the model never read
    print(f"{name} UNREAD-KEYS {unread}", flush=True)


doms = sorted({m for f in glob.glob(os.path.join(R, "wrf*_d*_*")) + glob.glob(os.path.join(out, "wrf*_d*_*")) +
               glob.glob(os.path.join(out, "d*", "wrf*_d*_*"))
               for m in re.findall(r"^wrf(?:out|rst)_(d\d\d)_", os.path.basename(f))})
if not doms:
    print(f"{name} NO-OUTPUT neither model wrote a wrfout or wrfrst file", flush=True)
for kind in ("wrfout", "wrfrst"):  # history frames, then restart files
    for d in doms:
        cdir = os.path.join(out, d) if os.path.isdir(os.path.join(out, d)) else out
        if not glob.glob(os.path.join(R, f"{kind}_{d}_*")) + glob.glob(os.path.join(cdir, f"{kind}_{d}_*")):
            continue
        bad, n = exact.compare_stream(glob.glob(os.path.join(cdir, f"{kind}_{d}_*")),
                                      glob.glob(os.path.join(R, f"{kind}_{d}_*")))
        tag = d if kind == "wrfout" else f"{d} {kind}"
        print(f"{name} {tag} frames {n} differing: {bad if bad else 'none'}", flush=True)
subprocess.run(["rm", "-rf", out] + glob.glob(os.path.join(R, "wrfout_*")) + glob.glob(os.path.join(R, "wrfrst_*")))
