#!/usr/bin/env python3
"""WRF-identity regression suite: every CPPWRF configuration that has a WRF
reference run, with the tolerance that the current code meets. Run after every change to the
dycore, boundaries, diffusion or physics coupling.

  wrf_suite.py --cpu-dp <build>/cppwrf --cuda <build>/cppwrf [--cuda-dp <build>/cppwrf] [--out DIR] [--only NAME[,NAME]]
               [--exact] [--jobs N]   (N CPU cases at once, each with OMP_NUM_THREADS threads) [--allow-skip] [--slow]
  wrf_suite.py --list   (every case with all its fields, then the special and restart checks)
  wrf_suite.py --status-table [FILE]   (docs/STATUS.md's table of the CPU FP64 cases; FILE: replace it there)
The cases are defined in tests/verify/cases/ (one file per group, see its __init__.py).

Cases with tolerance 0 (CPU FP64) run with WRF history output and their wrfout files are compared with WRF's at
tolerance 0 (exact.py: NaN equals only NaN, -0 differs from +0): every variable of the reference file, every
reference frame up to the case's time (the final one must be there), types and variable attributes, and the
global attributes when the reference's namelist.input is there (CPPWRF's history then reads it). The other
cases compare the model state with the reference's by relative L2 at the case's time (compare_real.py).
With the reference's namelist.input a tolerance-0 case FAILs when its configuration (config file and --set) sets
an option otherwise than tools/wrf2cppwrf.py's translation of that namelist.
Only the cases of the given binaries are selected; a given binary that does not exist is an error (exit 2).
Exit status 0 = every selected case PASSed; a case whose reference is missing is a SKIP and fails the suite unless
--allow-skip. The last line counts the results.
Slow cases (suite.slow in their --set list: the 24 h and 72 h runs) are left out unless --slow is given or --only --exact
names them.
"""
import argparse
import filecmp
import glob
import os
import re
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
# the cases, one file per group (tests/verify/cases/NN_<group>.py), and the restart checks (cases/restarts.py);
# WRF reference runs under ${CPPWRF_WRF_RUNS}, ideal-case references under ${CPPWRF_WRF_REF} (cases/common.py)
from cases import CASES, RESTARTS, STATUS  # noqa: E402
from cases.common import WR, WREF  # noqa: E402,F401

# the checks that are not in CASES: the Straka density current against WRF and the restart checks
SPECIAL = (["straka_wrf", "straka_wrf_cuda", "restart_p1_fp32_cuda"] + [f"restart_{t}_fp64" for t, _, _ in RESTARTS]
           + [f"restart_{t}_fp64_cuda" for t, _, _ in RESTARTS])


def run(binary, cfg, sets, out):
    shutil.rmtree(out, ignore_errors=True)
    cmd = [binary, os.path.join(REPO, cfg), "--set", f"output.dir={out}"]
    for s in sets:
        cmd += ["--set", s]
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "3"))
    r = subprocess.run(cmd, cwd=REPO, env=env, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def config_vs_namelist(cfg, sets, nml):
    """Options of the case's configuration (cfg and its --set) that differ from tools/wrf2cppwrf.py's translation of
    the namelist.input of the WRF run that made the reference (options both set; paths by their real path; the run
    length, output and restart settings aside), or the translation's errors. A case with output.auxhist3=true (the
    translation of io_form_auxhist3 = 2, which needs WRF history output) is translated with --format wrfout."""
    sys.path.insert(0, REPO)
    from tools.wrf2cppwrf import parse_namelist, WrfTranslator
    tr = WrfTranslator(parse_namelist(open(nml).read()), namelist_path=nml,
                       output_format="wrfout" if "output.auxhist3=true" in sets else "native")
    gen = tr.translate()
    if tr.errors:
        return [f"{nml} does not translate: {tr.errors}"]
    conf, sec = {}, None
    for line in open(os.path.join(REPO, cfg)):
        line = line.split("#")[0].strip()
        if re.fullmatch(r"\[.+\]", line):
            sec = line[1:-1]
        elif sec and "=" in line:
            conf.setdefault(sec, {})[line.split("=", 1)[0].strip()] = line.split("=", 1)[1].strip()
    for x in sets:
        k, v = x.split("=", 1)
        conf.setdefault(k.split(".", 1)[0], {})[k.split(".", 1)[1]] = v

    def norm(v):
        v = os.path.expandvars(str(v)).strip().strip("'\"").replace(" ", "")
        if v.startswith("/"):
            return os.path.realpath(v)
        try:
            return repr(float(v.replace("d", "e").replace("D", "e")))
        except ValueError:
            return {".true.": "true", ".false.": "false", "t": "true", "f": "false"}.get(v.lower(), v.lower())
    # a reference directory without its own wrfinput (tools/gcp/gcp_wrfref.sh links the base run's) does not name the
    # input directory
    skip = {("run", "run_time")} | (set() if os.path.exists(os.path.join(os.path.dirname(nml), "wrfinput_d01"))
                                    else {("real", "dir")})
    def same(a, b):  # files by their content (a reference directory's copy of a linked input_sounding)
        a, b = norm(a), norm(b)
        return a == b or (a.startswith("/") and os.path.isfile(a) and os.path.isfile(b) and filecmp.cmp(a, b, False))
    return [f"{s}.{k}: namelist {v} configuration {conf[s][k]}" for s, kv in gen.items() if s not in ("output", "restart")
            for k, v in kv.items() if (s, k) not in skip and k in conf.get(s, {}) and not same(v, conf[s][k])]


def restart_at(sets):
    """The time [s] of a case's restart (suite.restart_at=S in its --set list: the run stops there, writes wrfrst
    files and resumes from them), or None."""
    return next((float(x.split("=", 1)[1]) for x in sets if x.startswith("suite.restart_at=")), None)


def slow(sets):
    """A slow case (suite.slow in its --set list): run only with --slow or when --only --exact names it."""
    return "suite.slow" in sets


def exact_case(refs, out, t, sets):
    """Tolerance-0 comparison of the wrfout files in `out` (out/<domain> for nests; a list of directories: their
    files together, a later directory's frame replacing an earlier one's) with the references {domain: wrfout}:
    (differences, summary)."""
    import netCDF4
    import exact
    bad, frames, nvars = [], 0, 0
    for dom, w in refs.items():
        d = dom or re.search(r"wrfout_(d\d\d)_", os.path.basename(w)).group(1)
        with netCDF4.Dataset(w) as ds:
            xt = [float(x) * 60.0 for x in ds["XTIME"][:]]
            final = [exact.text(ds["Times"], i).decode() for i, x in enumerate(xt) if abs(x - t) < 1e-6]
            nvars += sum(1 for v in ds.variables.values() if v.dimensions[:1] == ("Time",) and v.dtype.kind != "S")
        if not final:
            bad.append(f"{d}: the reference has no frame at t={t:g}s")
            continue
        mine = [f for o in (out if isinstance(out, list) else [out])
                for f in glob.glob(os.path.join(o, f"wrfout_{d}_*")) + glob.glob(os.path.join(o, d, f"wrfout_{d}_*"))]
        dif, n = exact.compare_stream(mine, [w], subset=True, until=final[0],
                                      global_atts=os.path.exists(os.path.join(os.path.dirname(w), "namelist.input")))
        with netCDF4.Dataset(w) as ds:  # WRF leaves RQCBLTEN / RQIBLTEN uninitialised without qc / qi
            for v, q in (("RQCBLTEN", "QCLOUD"), ("RQIBLTEN", "QICE")):
                if q not in ds.variables:
                    dif.pop(v, None)
        frames += n
        bad += [f"{d}.{k}: {v}" for k, v in dif.items()]
        # output_diagnostics: the auxhist3 files next to the reference wrfout, compared as the wrfout files
        aux = sorted(glob.glob(os.path.join(os.path.dirname(w), f"auxhist3_{d}_*")))
        if aux:
            mine = [f for o in (out if isinstance(out, list) else [out])
                    for f in glob.glob(os.path.join(o, f"auxhist3_{d}_*")) + glob.glob(os.path.join(o, d, f"auxhist3_{d}_*"))]
            dif, n = exact.compare_stream(mine, aux, subset=True, until=final[0],
                                          global_atts=os.path.exists(os.path.join(os.path.dirname(w), "namelist.input")))
            if not n:
                dif.setdefault("(times)", "no auxhist3 frame in common")
            frames += n
            bad += [f"{d}.auxhist3.{k}: {v}" for k, v in dif.items()]
    return bad, f"{frames} frames, {nvars} variables"


def restarted_case(b, cfg, sets, refs, t, mid, out):
    """A run to `t` s restarted at `mid` s from its own wrfrst files: the first run (to one step past mid, wrfrst files
    at mid only: written within the run, as WRF's, not at its last step, where WRF's history alarm is still ringing)
    and the resumed one (mid to t) together against WRF's uninterrupted run at tolerance 0 (every reference frame),
    and the wrfrst files at mid against WRF's (the reference run's restart_interval writes them there):
    (differences, summary). As run_wrf_files_test.sh, RQCBLTEN / RQIBLTEN are not compared when WRF's restart file
    has no QCLOUD / QICE (WRF leaves them uninitialised)."""
    import netCDF4
    import exact
    first, resumed = os.path.join(out, "a_first"), os.path.join(out, "b_resumed")  # sorted: resumed frames win
    w0 = next(iter(refs.values()))
    with netCDF4.Dataset(w0) as ds:
        at = [exact.text(ds["Times"], i).decode() for i, x in enumerate(ds["XTIME"][:]) if abs(float(x) * 60.0 - mid) < 1e-6]
    if not at:
        return [f"{w0} has no frame at the restart time {mid:g}s"], ""
    nl = os.path.join(os.path.dirname(w0), "namelist.input")
    dt = float(re.search(r"^\s*time_step\s*=\s*(\d+)", open(nl).read(), re.M).group(1))
    rc, log = run(b, cfg, sets + [f"run.run_time={mid + dt:g}", f"restart.interval={mid:g}", f"restart.end={mid:g}"],
                  first)
    if rc != 0:
        return [f"run to the restart failed: {log[-200:]}"], ""
    # a nested run writes each domain's files in first/<domain>; the resumed run reads every domain's wrfrst file
    # from the directory of d01's (as WRF's run directory)
    rst = os.path.join(out, "rst")
    os.makedirs(rst, exist_ok=True)
    for pat in (f"wrfrst_d??_{at[0]}", f"d??/wrfrst_d??_{at[0]}"):
        for f in glob.glob(os.path.join(first, pat)):
            os.symlink(f, os.path.join(rst, os.path.basename(f)))
    rc, log = run(b, cfg, sets + [f"run.restart_file={rst}/wrfrst_d01_{at[0]}"], resumed)
    if rc != 0:
        return [f"resumed run failed: {log[-200:]}"], ""
    bad, what = exact_case(refs, [first, resumed], t, sets)
    for dom, w in refs.items():
        d = dom or re.search(r"wrfout_(d\d\d)_", os.path.basename(w)).group(1)
        wr = os.path.join(os.path.dirname(w), f"wrfrst_{d}_{at[0]}")
        mine = glob.glob(os.path.join(first, f"wrfrst_{d}_{at[0]}")) + glob.glob(os.path.join(first, d, f"wrfrst_{d}_{at[0]}"))
        if not os.path.exists(wr) or not mine:
            bad.append(f"{d}: no wrfrst_{d}_{at[0]} from {'WRF' if not os.path.exists(wr) else 'CPPWRF'}")
            continue
        with netCDF4.Dataset(wr) as ds:
            noise = [x for x, q in (("RQCBLTEN", "QCLOUD"), ("RQIBLTEN", "QICE")) if q not in ds.variables]
        dif, _ = exact.compare_stream(mine, [wr])
        bad += [f"{d}.wrfrst.{k}: {v}" for k, v in dif.items() if k not in noise]
    return bad, f"{what}, restarted at {mid:g}s (wrfrst identical)" if not bad else what


def compare(wrfout, out, t):
    r = subprocess.run([sys.executable, os.path.join(REPO, "tests/verify/compare_real.py"), wrfout, out],
                       capture_output=True, text=True)
    res = {}
    for line in r.stdout.splitlines():
        f = line.split()
        if len(f) >= 5 and f[0].replace(".", "").isdigit() and abs(float(f[0]) - t) < 1e-6:
            res[f[1]] = float(f[3])
    return res


def straka(binary, out):
    wo = f"{WREF}/run_straka_d1/wrfout_d01_0001-01-01_00:00:00"
    if not os.path.exists(wo):
        return None, "reference missing"
    rc, log = run(binary, "configs/straka_wrf.ini", [], out)
    if rc != 0:
        return False, log[-300:]
    r = subprocess.run([sys.executable, os.path.join(REPO, "tests/verify/verify.py"), "wrf", wo, out],
                       capture_output=True, text=True)
    return ("FAIL" not in r.stdout and r.returncode == 0), r.stdout.strip().splitlines()[-1] if r.stdout else r.stderr


def restart_p1(binary, out, cfg="configs/real_p1_r8.ini", run_dir="run_P1_r8"):
    """Physics/TKE restart: 3 steps straight vs restarted after step 1 must be bit-identical
    (dycore, TKE, surface/soil state, R*BLTEN, lagged t_1; for P7 also R*CUTEN and RTHFTEN/RQVFTEN)."""
    if not os.path.exists(f"{WR}/{run_dir}/steps/wrfinput_d01") and not os.path.exists(f"{WR}/{run_dir}/wrfinput_d01"):
        return None, "input missing"
    a, b = os.path.join(out, "straight"), os.path.join(out, "restarted")
    rc, log = run(binary, cfg, ["restart.interval=60"], a)
    if rc != 0:
        return False, log[-300:]
    rc, log = run(binary, cfg, [f"run.restart_file={a}/restart_000001.nc"], b)
    if rc != 0:
        return False, log[-300:]
    same = open(f"{a}/snap_000003.bin", "rb").read() == open(f"{b}/snap_000003.bin", "rb").read()
    return same, "step 3 bit-identical" if same else "step 3 differs after restart"


STATUS_MARK = "<!-- generated: tests/verify/wrf_suite.py --status-table docs/STATUS.md (rows: STATUS in tests/verify/cases) -->"


def status_table(path=None):
    """The markdown table of STATUS (the groups' rows); with `path`, it replaces the table after STATUS_MARK there."""
    table = ["| Case | Configuration | Duration |", "|---|---|---|"] + [f"| {c} | {d} | {t} |" for c, d, t in STATUS]
    if not path:
        print("\n".join(table))
        return
    lines = open(path).read().split("\n")
    i = lines.index(STATUS_MARK) + 1
    j = i
    while j < len(lines) and lines[j].startswith("|"):
        j += 1
    open(path, "w").write("\n".join(lines[:i] + table + lines[j:]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpu-dp")
    ap.add_argument("--cuda")
    ap.add_argument("--cuda-dp")
    ap.add_argument("--out", default=os.path.expanduser("~/cppwrf_scratch/wrf_suite"))
    ap.add_argument("--only", help="comma-separated substrings of case names")
    ap.add_argument("--exact", action="store_true", help="--only names are exact case names")
    ap.add_argument("--jobs", type=int, default=1, help="CPU cases run at the same time (GPU cases: one at a time)")
    ap.add_argument("--keep", action="store_true", help="keep the outputs of passing cases (default: deleted)")
    ap.add_argument("--allow-skip", action="store_true", help="cases without their reference do not fail the suite")
    ap.add_argument("--slow", action="store_true", help="also the slow cases (suite.slow: the 24 h and 72 h runs)")
    ap.add_argument("--list", action="store_true", help="print every case (all its fields) and the restart and special "
                    "checks, and exit")
    ap.add_argument("--status-table", nargs="?", const="", metavar="FILE",
                    help="print docs/STATUS.md's case table (FILE: replace it there) and exit")
    a = ap.parse_args()
    if a.status_table is not None:
        status_table(a.status_table)
        return
    if a.list:
        for c in CASES:
            print(" | ".join(repr(x) for x in c))
        for n in SPECIAL:
            print(n)
        for r in RESTARTS:
            print(" | ".join(repr(x) for x in r))
        return
    only = a.only.split(",") if a.only else None

    def wanted(name):
        return not only or (name in only if a.exact else any(o in name for o in only))

    unknown = [o for o in only or [] if not any(o == n if a.exact else o in n for n in [c[0] for c in CASES] + SPECIAL)]
    if unknown:  # a misspelt case would otherwise make an empty suite that passes
        print(f"ERROR no case named {', '.join(unknown)}", flush=True)
        sys.exit(2)

    def done(path, good):  # the outputs of passing cases are deleted (they take tens of GB)
        if good and not a.keep:
            shutil.rmtree(path, ignore_errors=True)

    def case(name, b, cfg, sets, refs, t, tol):
        out = os.path.join(a.out, name)
        mid = restart_at(sets)
        sets = [x for x in sets if not x.startswith("suite.")]
        exact = all(x == 0.0 for x in tol.values())
        if exact:  # WRF history output; with the reference's namelist.input CPPWRF writes WRF's global attributes
            nl = os.path.join(os.path.dirname(next(iter(refs.values()))), "namelist.input")
            if os.path.exists(nl):  # the configuration is the WRF run's that made the reference
                dif = config_vs_namelist(cfg, sets, nl)
                if dif:
                    return False, f"FAIL  {name:22s} configuration differs from {nl}: {'; '.join(dif)}"
            sets = sets + ["output.format=wrfout"] + ([f"output.wrf_namelist={nl}"] if os.path.exists(nl) else [])
        if exact and mid:
            bad, what = restarted_case(b, cfg, sets, refs, t, mid, out)
            done(out, not bad)
            return not bad, (f"{'PASS' if not bad else 'FAIL'}  {name:22s} t={t:g}s exact: {what}"
                             + (f"  [{'; '.join(x[:300] for x in bad[:6])}{' ...' if len(bad) > 6 else ''}]" if bad else ""))
        rc, log = run(b, cfg, sets, out)
        if rc != 0:
            return False, f"FAIL  {name:22s} run failed: {log[-200:]}"
        if exact:
            bad, what = exact_case(refs, out, t, sets)
            done(out, not bad)
            return not bad, (f"{'PASS' if not bad else 'FAIL'}  {name:22s} t={t:g}s exact: {what}"
                             + (f"  [{'; '.join(x[:300] for x in bad[:6])}{' ...' if len(bad) > 6 else ''}]" if bad else ""))
        res = {}
        for dom, w in refs.items():  # nested: out/<domain>, tolerance keys "<domain>.<var>"
            for v, x in compare(w, os.path.join(out, dom) if dom else out, t).items():
                res[f"{dom}.{v}" if dom else v] = x
        bad = [f"{v} {res.get(v, float('nan')):.2e}>{tl:.0e}" for v, tl in tol.items()
               if not (res.get(v, float("inf")) <= tl)]
        txt = " ".join(f"{v}={res.get(v, float('nan')):.1e}" for v in tol)
        done(out, not bad)
        return not bad, f"{'PASS' if not bad else 'FAIL'}  {name:22s} t={t:g}s {txt}" + (f"  [{'; '.join(bad)}]" if bad else "")

    def special(name, fn, path):
        good, msg = fn()
        if good is None:
            return None, f"SKIP  {name:22s} ({msg})"
        done(path, good)
        return good, f"{'PASS' if good else 'FAIL'}  {name:22s} {msg}"

    bins = {"cpu_dp": a.cpu_dp, "cuda": a.cuda, "cuda_dp": a.cuda_dp}
    for k, b in bins.items():  # a binary that was asked for and is not there (a failed build) is no SKIP
        if b and not os.path.exists(b):
            print(f"ERROR --{k.replace('_', '-')} {b}: no such file (build failed?)", flush=True)
            sys.exit(2)
    cpu, gpu = [], []
    count = {"PASS": 0, "FAIL": 0, "SKIP": 0}
    for name, key, cfg, sets, wo, t, tol in CASES:
        if not wanted(name) or not bins.get(key):  # cases of binaries not given are not selected
            continue
        if slow(sets) and not a.slow and not (a.exact and name in only):  # slow cases only on request
            continue
        b = bins[key]
        refs = wo if isinstance(wo, dict) else {"": wo}  # nested runs: {domain: wrfout}
        missing = [w for w in refs.values() if not os.path.exists(w)]
        if missing:
            print(f"SKIP  {name:22s} (reference missing: {', '.join(missing)})", flush=True)
            count["SKIP"] += 1
            continue
        (cpu if key == "cpu_dp" else gpu).append((t, name, lambda n=name, b=b, c=cfg, s=sets, r=refs, t=t, tl=tol: case(n, b, c, s, r, t, tl)))
    if a.cpu_dp and os.path.exists(a.cpu_dp):
        if not only or wanted("straka_wrf"):
            p = os.path.join(a.out, "straka")
            cpu.append((0, "straka_wrf", lambda p=p: special("straka_wrf", lambda: straka(a.cpu_dp, p), p)))
        for tag, cfg, rd in RESTARTS:
            rname = f"restart_{tag}_fp64"
            if not wanted(rname):
                continue
            p = os.path.join(a.out, f"restart_{tag}")
            cpu.append((0, rname, lambda r=rname, c=cfg, d=rd, p=p: special(r, lambda: restart_p1(a.cpu_dp, p, c, d), p)))
    if a.cuda and os.path.exists(a.cuda) and wanted("restart_p1_fp32_cuda"):  # single precision, GPU
        p = os.path.join(a.out, "restart_p1_fp32_cuda")
        gpu.append((0, "restart_p1_fp32_cuda", lambda p=p: special("restart_p1_fp32_cuda", lambda: restart_p1(a.cuda, p), p)))
    if a.cuda_dp and os.path.exists(a.cuda_dp):
        if wanted("straka_wrf_cuda"):
            p = os.path.join(a.out, "straka_cuda")
            gpu.append((0, "straka_wrf_cuda", lambda p=p: special("straka_wrf_cuda", lambda: straka(a.cuda_dp, p), p)))
        for tag, cfg, rd in RESTARTS:
            rname = f"restart_{tag}_fp64_cuda"
            if not wanted(rname):
                continue
            p = os.path.join(a.out, f"restart_{tag}_cuda")
            gpu.append((0, rname, lambda r=rname, c=cfg, d=rd, p=p: special(r, lambda: restart_p1(a.cuda_dp, p, c, d), p)))
    # longest simulated time first; each task prints its line as soon as it ends; GPU cases share the one GPU
    cpu = sorted(cpu, key=lambda x: -x[0])
    from concurrent.futures import ThreadPoolExecutor, as_completed

    def timed(fn):  # the case's line ends with its wall time
        t0 = time.monotonic()
        good, line = fn()
        return good, f"{line}  ({time.monotonic() - t0:.0f} s)"
    with ThreadPoolExecutor(max(1, a.jobs)) as pc, ThreadPoolExecutor(1) as pg:
        futs = {pc.submit(timed, f): n for _, n, f in cpu} | {pg.submit(timed, f): n for _, n, f in gpu}
        for f in as_completed(futs):
            try:
                good, line = f.result()
            except Exception as e:  # noqa: BLE001 - a case the harness cannot finish fails, the others go on
                good, line = False, f"FAIL  {futs[f]:22s} harness error: {type(e).__name__}: {e}"
            print(line, flush=True)
            count["SKIP" if good is None else "PASS" if good else "FAIL"] += 1
    print(f"SUITE PASS {count['PASS']} FAIL {count['FAIL']} SKIP {count['SKIP']}", flush=True)
    sys.exit(0 if not count["FAIL"] and (a.allow_skip or not count["SKIP"]) else 1)


if __name__ == "__main__":
    main()
