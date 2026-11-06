#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""The runs of tools/bench/colab_bench.sh in the Debian 13 chroot of the Colab runtime (tools/bench/bench_body.sh).

  bench_run.py make-large BASE MET_DIR       the large domain's WRF run directory /opt/wrfbench/run_L: BASE's
                                             namelist and tables, real.exe (WRF real*8 build) on MET_DIR's met_em files
  bench_run.py case NAME BASE [--ref] [--serial] [--set SECTION.key=VALUE ...] [--layouts L,...]
  bench_run.py report                        /opt/job/report.md from /opt/job/results.txt (tools/bench/bench_report.py)

A case runs BASE's namelist for 1 simulated hour (run_hours = 1, restart off, all frames in one file per domain) with
the case's changes, once with history every 10 minutes ("on") and once with history_interval = 1440 ("off": WRF and
CPPWRF write only the initial frame):
- WRF v4.8.0 real*4 and real*8 (builds in /opt/wrfbench/WRF-r4, WRF-r8; dm+sm builds WRF-r4sm, WRF-r8sm) with each
  process layout of --layouts (dmN: N MPI ranks; smNxT: N ranks with T OpenMP threads each) and history off, then
  with the fastest layout (smallest stepping time) and history on;
- CPPWRF from tools/wrf2cppwrf.py's translation of the same namelist: the CUDA FP32 and FP64 builds on the GPU
  (tools/bench/bench_gpu.sh runs them outside the chroot: the queue /opt/wrfbench/gpu/queue_NAME.txt, then
  go_NAME; it answers with res_NAME.txt and done_NAME), the CPU FP64 build with OpenMP on every vCPU (cpu-dp: the
  default build, CPPWRF_PORTABLE_LIBMVEC on) and, history on only, the CPU FP64 build with this CPU's own libmvec
  log (cpu-dp-hwlog: CPPWRF_PORTABLE_LIBMVEC off, as the WRF executables built here).
Then the outputs are compared (tests/verify/exact.py): CPU FP64 against BASE's wrfout files (--ref: WRF real*8
references of the suite, made on Intel), GPU FP64 against CPU FP64, the hwlog build and the default build against
WRF real*8 of this runtime (tolerance 0), and the relative L2 differences at 1 h of GPU FP32 against WRF real*4,
GPU FP32 against WRF real*8 and WRF real*4 against WRF real*8. With --ref also WRF real*8 of this runtime against the
references; --serial adds a WRF real*8 run in one process (dm1, history on) to these comparisons.
Each run and comparison is a JSON line of /opt/job/results.txt; the logs (rsl.out.0000, CPPWRF's output, namelists)
go to /opt/job/fetch/logs/NAME/, the case's run directories are deleted at the end (the wrfout files are large).

Times: WRF's stepping is the sum of rsl.out.0000's "Timing for main" lines of domain 1 (a nest's steps inside its
parent's) minus the "Timing for Writing" lines among them (history output in med_before_solve_io); its initialisation
the "Timing for processing wrfinput file" and "lateral boundary" lines. CPPWRF's stepping is the first "cppwrf: done"
line (the parent domain's, which contains the nests' steps; output excluded). Wall: the whole process (mpiexec /
cppwrf), measured here.
"""
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import time

W = "/opt/wrfbench"
J = "/opt/job"
SRC = J + "/src"
RES = J + "/results.txt"
LOGS = J + "/fetch/logs"
Q = W + "/gpu"
NP = os.cpu_count() or 1
sys.path.insert(0, SRC + "/tests/verify")
L2_VARS = ("U", "V", "W", "T", "QVAPOR", "PH", "MU", "PSFC")


def log(msg):
    with open(J + "/bench.log", "a") as f:
        f.write(time.strftime("%H:%M:%S ") + msg + "\n")


def emit(rec):
    with open(RES, "a") as f:
        f.write(json.dumps(rec) + "\n")


LOAD = []


def sh(cmd, cwd, env=None, timeout=3600, out="run.log"):
    """Runs cmd (shell) in cwd with output to cwd/out: (exit status, wall seconds); LOAD[0]: the processes running just
    before the run, this one aside (/proc/stat procs_running: another process on the CPUs would show here)."""
    LOAD[:] = [max(int(x.split()[1]) for x in open("/proc/stat") if x.startswith("procs_running")) - 1]
    t0 = time.time()
    try:
        with open(os.path.join(cwd, out), "w") as f:
            rc = subprocess.run(cmd, shell=True, cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT,
                                timeout=timeout).returncode
    except subprocess.TimeoutExpired:
        rc = 124
    return rc, time.time() - t0


# ---- namelist (as tools/gcp/matrix_variant.py)
def setkey(nl, sec, key, val):
    m = re.search(rf"&{sec}\b(.*?)\n\s*/", nl, re.S | re.I)
    if not m:
        return nl.rstrip("\n") + f"\n\n&{sec}\n {key} = {val},\n/\n"
    body = m.group(1)
    pat = re.compile(rf"^\s*{key}\s*=.*$", re.M | re.I)
    line = f" {key} = {val},"
    body = pat.sub(line, body) if pat.search(body) else body + "\n" + line
    return nl[:m.start(1)] + body + nl[m.end(1):]


def vals(nl, key):
    m = re.search(rf"^\s*{key}\s*=\s*([^\n]*)$", nl, re.M | re.I)
    return [x.strip().strip("'\"") for x in m.group(1).split(",") if x.strip()] if m else []


def geti(nl, key, d=0, default=0):
    v = vals(nl, key)
    return int(float(v[d])) if len(v) > d else default


def grid(nl):
    """[cells, steps per hour] of each domain: (e_we - 1) (e_sn - 1) (e_vert - 1) and 3600 / dt (the parent's dt over
    parent_time_step_ratio)."""
    nd = geti(nl, "max_dom", 0, 1)
    dt0 = geti(nl, "time_step") + geti(nl, "time_step_fract_num") / max(1, geti(nl, "time_step_fract_den", 0, 1))
    dts, out = [], []
    for d in range(nd):
        p = geti(nl, "parent_id", d, 0) - 1
        dt = dt0 if d == 0 else dts[p] / geti(nl, "parent_time_step_ratio", d, 1)
        dts.append(dt)
        out.append([(geti(nl, "e_we", d) - 1) * (geti(nl, "e_sn", d) - 1) * (geti(nl, "e_vert", d) - 1),
                    round(3600 / dt)])
    return out


def make_dir(case, name, base, sets, hist):
    """W/runs/CASE/NAME: BASE's inputs (links) and its namelist for 1 h with history on / off and the case's sets."""
    R = f"{W}/runs/{case}/{name}"
    shutil.rmtree(R, ignore_errors=True)
    os.makedirs(R)
    for f in os.listdir(base):
        p = os.path.join(base, f)
        if os.path.isdir(p) or re.match(r"(wrfout.*|wrfrst.*|auxhist.*|rsl\..*|namelist\.(output|input)|.*\.log|.*\.exe"
                                        r"|met_em.*|FILE:.*|GRIBFILE.*)$", f):
            continue
        os.symlink(os.path.realpath(p), os.path.join(R, f))
    nl = open(os.path.join(base, "namelist.input")).read()
    nd = geti(nl, "max_dom", 0, 1)
    for k in ("history_interval", "history_interval_d", "history_interval_h", "history_interval_m",
              "history_interval_s"):
        nl = re.sub(rf"^\s*{k}\s*=.*\n", "", nl, flags=re.M | re.I)
    per = lambda v: ", ".join([str(v)] * nd)  # noqa: E731
    for k, v in (("run_days", "0"), ("run_hours", "1"), ("run_minutes", "0"), ("run_seconds", "0"),
                 ("restart", ".false."), ("restart_interval", "1440"), ("frames_per_outfile", per(1000)),
                 ("history_interval", per(10 if hist == "on" else 1440)), ("debug_level", "0")):
        nl = setkey(nl, "time_control", k, v)
    for s in sets:
        sk, v = s.split("=", 1)
        sec, k = sk.split(".", 1)
        nl = setkey(nl, sec, k, v)
    open(os.path.join(R, "namelist.input"), "w").write(nl)
    return R, nl


def save(case, R, files, prefix):
    os.makedirs(f"{LOGS}/{case}", exist_ok=True)
    for f in files:
        if os.path.exists(os.path.join(R, f)):
            shutil.copy(os.path.join(R, f), f"{LOGS}/{case}/{prefix}.{f}")


# ---- WRF
CORES = None
BIND = True


def topology():
    global CORES, BIND
    out = subprocess.run("lscpu -p=CORE,SOCKET | grep -v '^#' | sort -u | wc -l", shell=True, capture_output=True,
                         text=True).stdout.strip()
    CORES = int(out) if out.isdigit() and int(out) > 0 else NP
    BIND = subprocess.run("mpiexec -n 2 -bind-to core:1 true", shell=True, capture_output=True,
                          timeout=120).returncode == 0


def bind(n):
    """Each rank on CORES / n physical cores (its OpenMP threads too); more ranks than cores: one hardware thread."""
    if not BIND:
        return ""
    return f"-bind-to core:{CORES // n}" if n <= CORES else "-bind-to hwthread"


def wrf_timings(text):
    lines = text.splitlines()
    main = []
    for i, l in enumerate(lines):
        m = re.search(r"Timing for main: time \S+ on domain\s+(\d+):\s+([0-9.]+) elapsed", l)
        if m and m.group(1) == "1":
            main.append((i, float(m.group(2))))
    last = main[-1][0] if main else -1
    writes = [(i, float(m.group(1))) for i, l in enumerate(lines)
              if (m := re.search(r"Timing for Writing .* for domain\s+\d+:\s+([0-9.]+) elapsed", l))]
    w_in = sum(x for i, x in writes if i < last)
    init = sum(float(m.group(1)) for l in lines
               if (m := re.search(r"Timing for processing (?:wrfinput file|lateral boundary).*:\s+([0-9.]+) elapsed", l)))
    msum = sum(x for _, x in main)
    return dict(steps=len(main), main=msum, write_in=w_in, write=sum(x for _, x in writes), stepping=msum - w_in,
                init=init)


def wrf_run(case, base, sets, prec, layout, hist):
    m = re.match(r"(dm|sm)(\d+)(?:x(\d+))?$", layout)
    n, t = int(m.group(2)), int(m.group(3) or 1)
    exe = f"{W}/WRF-{prec}{'sm' if m.group(1) == 'sm' else ''}/main/wrf.exe"
    rec = dict(kind="run", case=case, model=f"wrf_{prec}", layout=layout, hist=hist, ok=False)
    if not os.path.exists(exe):
        rec["error"] = "no executable"
        emit(rec)
        return rec, None
    name = f"wrf{prec}_{layout}_{hist}"
    R, nl = make_dir(case, name, base, sets, hist)
    os.symlink(exe, R + "/wrf.exe")
    if prec == "r4":  # the suite's run directories are real*8 runs: their RRTMG_LW_DATA ... are WRF's *_DATA_DBL files
        run = f"{W}/WRF-r4/run"
        for f in os.listdir(run):
            if os.path.exists(f"{run}/{f}_DBL") and os.path.lexists(f"{R}/{f}"):
                os.remove(f"{R}/{f}")
                os.symlink(os.path.realpath(f"{run}/{f}"), f"{R}/{f}")
    env = dict(os.environ, OMP_NUM_THREADS=str(t), OMP_STACKSIZE="512M")
    cmd = f"ulimit -s unlimited; exec mpiexec -n {n} {bind(n)} ./wrf.exe"
    rc, wall = sh(cmd, R, env, out="wrf.log")
    text = open(R + "/rsl.out.0000", errors="replace").read() if os.path.exists(R + "/rsl.out.0000") else ""
    rec.update(ok=rc == 0 and "SUCCESS COMPLETE WRF" in text, exit=rc, wall=wall, cmd=cmd, threads=t, ranks=n, load=LOAD[0],
               grid=grid(nl), **wrf_timings(text))
    if not rec["ok"]:
        err = open(R + "/rsl.error.0000", errors="replace").read() if os.path.exists(R + "/rsl.error.0000") else ""
        rec["error"] = [x.strip() for x in (err + open(R + "/wrf.log", errors="replace").read()).splitlines()
                        if re.search(r"ERROR|FATAL|signal|[Ee]rror", x)][:6]
    save(case, R, ["rsl.out.0000", "wrf.log", "namelist.input"] + ([] if rec["ok"] else ["rsl.error.0000"]), name)
    emit(rec)
    log(f"{case} {name}: {'ok' if rec['ok'] else 'FAILED'} wall {wall:.1f} s, stepping {rec['stepping']:.2f} s, "
        f"{rec['steps']} steps ({LOAD[0]} other processes running before)")
    return rec, R


# ---- CPPWRF
def cpp_timings(text):
    done = re.findall(r"cppwrf: done, (\d+) steps, stepping ([0-9.]+) s, ([0-9.]+) ns/cell/step", text)
    mem = re.findall(r"cppwrf: device memory ([^:]*): ([0-9]+) MiB used", text)
    r = dict(domains=[[int(k), float(s), float(n)] for k, s, n in done], meminfo={w: int(v) for w, v in mem})
    if done:
        r.update(steps=int(done[0][0]), stepping=float(done[0][1]))
    return r


def translate(R):
    t = subprocess.run([sys.executable, SRC + "/tools/wrf2cppwrf.py", R + "/namelist.input", "--real-dir", R,
                        "--format", "wrfout", "-o", R + "/case.ini"], capture_output=True, text=True)
    open(R + "/translate.log", "w").write(t.stdout + t.stderr)
    return t.returncode == 0


def cpu_run(case, R, nl, hist, model, build):
    out = f"{R}/out_{model}"
    exe = f"{W}/{build}/cppwrf"
    rec = dict(kind="run", case=case, model=model, layout=f"omp{NP}", hist=hist, ok=False, grid=grid(nl))
    if not os.path.exists(exe):
        rec["error"] = "no executable"
        emit(rec)
        return rec
    env = dict(os.environ, OMP_NUM_THREADS=str(NP))
    rc, wall = sh(f"exec {exe} {R}/case.ini --set output.dir={out}", R, env, out=f"{model}.log")
    text = open(f"{R}/{model}.log", errors="replace").read()
    rec.update(exit=rc, wall=wall, threads=NP, load=LOAD[0], **cpp_timings(text))
    rec["ok"] = rc == 0 and "stepping" in rec
    if not rec["ok"]:
        rec["error"] = text.strip().splitlines()[-4:]
    save(case, R, [f"{model}.log"], f"cpp_{hist}")
    emit(rec)
    log(f"{case} {model} {hist}: {'ok' if rec['ok'] else 'FAILED'} wall {wall:.1f} s, stepping {rec.get('stepping', 0):.2f} s")
    return rec


def gpu_wait_ready(timeout):
    t0 = time.time()
    while time.time() - t0 < timeout:
        if os.path.exists(Q + "/ready"):
            return open(Q + "/ready").read().strip()
        time.sleep(10)
    return "timeout"


def gpu_runs(case, jobs, timeout=7200):
    """jobs: [(tag, build, R, nl, hist, reps)] run by tools/bench/bench_gpu.sh; one record per job (the repetition
    with the smallest stepping time)."""
    ready = gpu_wait_ready(3600)
    recs = []
    with open(f"{Q}/queue_{case}.txt", "w") as f:
        for tag, build, R, nl, hist, reps in jobs:
            f.write(f"{tag} {build} {R} {R}/case.ini {R}/out_{tag} {reps}\n")
    if "cuda=0" not in ready and "cuda_dp=0" not in ready:
        log(f"{case}: no GPU build ({ready})")
    open(f"{Q}/go_{case}", "w").write("go\n")
    t0 = time.time()
    while not os.path.exists(f"{Q}/done_{case}") and time.time() - t0 < timeout:
        time.sleep(5)
    res = {}
    if os.path.exists(f"{Q}/res_{case}.txt"):
        for line in open(f"{Q}/res_{case}.txt"):
            try:
                r = json.loads(line)
            except ValueError:
                continue
            res.setdefault(r["tag"], []).append(r)
    for tag, build, R, nl, hist, reps in jobs:
        model = tag.rsplit("_", 1)[0]
        best, steppings = None, []
        for r in res.get(tag, []):
            log_f = f"{R}/{tag}.rep{r['rep']}.log"
            text = open(log_f, errors="replace").read() if os.path.exists(log_f) else ""
            rec = dict(kind="run", case=case, model=model, layout="gpu", hist=hist, grid=grid(nl), exit=r["exit"],
                       wall=r["wall"], mem_idle=r.get("idle"), mem_peak=r.get("peak"), rep=r["rep"], reps=reps,
                       **cpp_timings(text))
            rec["ok"] = r["exit"] == 0 and "stepping" in rec
            if not rec["ok"]:
                rec["error"] = text.strip().splitlines()[-4:]
            steppings.append(rec.get("stepping"))
            if best is None or (rec["ok"] and (not best["ok"] or rec["stepping"] < best["stepping"])):
                best = rec
            save(case, R, [f"{tag}.rep{r['rep']}.log"], f"cpp_{hist}")
        if best is None:
            best = dict(kind="run", case=case, model=model, layout="gpu", hist=hist, ok=False, error="no result")
        best["rep_stepping"] = steppings
        emit(best)
        log(f"{case} {tag}: {'ok' if best['ok'] else 'FAILED'} wall {best.get('wall', 0):.1f} s, "
            f"stepping {best.get('stepping', 0):.2f} s, GPU memory peak {best.get('mem_peak')} MiB")
        recs.append(best)
    return recs


# ---- comparisons
def files(d, dom):
    """dom's wrfout files in d (WRF's run directory) or d/dom (CPPWRF's nested runs: output.dir/d01, d02, ...)."""
    return sorted(glob.glob(f"{d}/wrfout_{dom}_*") + glob.glob(f"{d}/{dom}/wrfout_{dom}_*"))


def doms(nl):
    return [f"d{i:02d}" for i in range(1, geti(nl, "max_dom", 0, 1) + 1)]


def cmp_exact(case, what, A, B, nl, subset=False):
    import exact
    rec = dict(kind="exact", case=case, what=what, domains={})
    for d in doms(nl):
        fa, fb = files(A, d), files(B, d)
        if not fa or not fb:
            rec["domains"][d] = dict(error=f"no files ({len(fa)} / {len(fb)})")
            continue
        try:
            bad, n = exact.compare_stream(fa, fb, subset=subset, global_atts=False)
            rec["domains"][d] = dict(frames=n, differ={k: str(v)[:200] for k, v in bad.items()})
        except Exception as e:  # noqa: BLE001 - a comparison that fails is a result line
            rec["domains"][d] = dict(error=f"{type(e).__name__}: {e}")
    emit(rec)
    log(f"{case} {what}: " + "; ".join(f"{d} " + (v.get("error") or (f"{v['frames']} frames, identical" if not v["differ"]
                                                                 else f"{v['frames']} frames, differ {sorted(v['differ'])[:12]}"))
                                        for d, v in rec["domains"].items()))


def cmp_l2(case, what, A, B, nl):
    """Relative L2 difference ||a - b|| / ||b|| of L2_VARS at the last common frame (B is the reference)."""
    import numpy as np
    import exact
    rec = dict(kind="l2", case=case, what=what, domains={})
    for d in doms(nl):
        opened = []
        try:
            fa, fb = exact.frames(files(A, d), opened), exact.frames(files(B, d), opened)
            common = sorted(set(fa) & set(fb))
            if not common:
                rec["domains"][d] = dict(error="no common frame")
                continue
            t = common[-1]
            (a, ia), (b, ib) = fa[t], fb[t]
            out = {}
            for v in L2_VARS:
                if v in a.variables and v in b.variables:
                    x = np.asarray(a[v][ia], dtype=np.float64)
                    y = np.asarray(b[v][ib], dtype=np.float64)
                    ny = float(np.sqrt(np.sum(y * y)))
                    out[v] = float(np.sqrt(np.sum((x - y) ** 2))) / ny if ny > 0 else None
            rec["domains"][d] = dict(time=t, l2=out)
        except Exception as e:  # noqa: BLE001
            rec["domains"][d] = dict(error=f"{type(e).__name__}: {e}")
        finally:
            for ds in opened:
                ds.close()
    emit(rec)
    log(f"{case} {what}: " + "; ".join(f"{d} " + (v.get("error") or " ".join(
        f"{k} {x:.2e}" for k, x in v["l2"].items() if x is not None)) for d, v in rec["domains"].items()))


# ---- one case
def case(name, base, sets, ref, layouts, serial=False):
    log(f"== case {name}: {base} {' '.join(sets)}")
    shutil.rmtree(f"{W}/runs/{name}", ignore_errors=True)
    best = {}
    for prec in ("r4", "r8"):  # WRF: the layouts with history off, the fastest with history on
        good = []
        for lay in layouts:
            rec, R = wrf_run(name, base, sets, prec, lay, "off")
            if R:
                shutil.rmtree(R, ignore_errors=True)
            if rec.get("ok"):
                good.append((rec["stepping"], lay))
        if good:
            best[prec] = min(good)[1]
            emit(dict(kind="best", case=name, model=f"wrf_{prec}", layout=best[prec]))
            log(f"{name} WRF {prec}: fastest layout {best[prec]}")
    wrf_on = {}
    for prec, lay in best.items():
        rec, R = wrf_run(name, base, sets, prec, lay, "on")
        if rec.get("ok"):
            wrf_on[prec] = R
    if serial and "r8" in best:  # WRF real*8 in one process (as the suite's references): its layout's effect on the bits
        rec, R = wrf_run(name, base, sets, "r8", "dm1", "on")
        if rec.get("ok"):
            wrf_on["r8_dm1"] = R
    # CPPWRF: one translated directory per history setting
    dirs = {}
    for hist in ("off", "on"):
        R, nl = make_dir(name, f"cpp_{hist}", base, sets, hist)
        if translate(R):
            dirs[hist] = (R, nl)
        else:
            log(f"{name} {hist}: wrf2cppwrf.py FAILED: {open(R + '/translate.log').read().strip().splitlines()[-1:]}")
            emit(dict(kind="run", case=name, model="translate", hist=hist, ok=False,
                      error=open(R + "/translate.log").read().strip().splitlines()[-3:]))
        save(name, R, ["translate.log", "case.ini", "namelist.input"], f"cpp_{hist}")
    jobs = []
    for hist, (R, nl) in dirs.items():
        reps = 2 if hist == "off" else 1  # the GPU's first run of a binary pays the module load
        jobs += [(f"gpu_fp32_{hist}", "cuda", R, nl, hist, reps), (f"gpu_fp64_{hist}", "cuda_dp", R, nl, hist, reps)]
    if jobs:
        gpu_runs(name, jobs)
    for hist, (R, nl) in dirs.items():
        cpu_run(name, R, nl, hist, "cpu_fp64", "cpu-dp")
    if "on" in dirs:
        R, nl = dirs["on"]
        cpu_run(name, R, nl, "on", "cpu_fp64_hwlog", "cpu-dp-hwlog")
        if ref:
            cmp_exact(name, "cpu_fp64 vs WRF r8 reference (suite, Intel)", R + "/out_cpu_fp64", base, nl, subset=True)
            for k in ("r8", "r8_dm1"):
                if k in wrf_on:
                    cmp_exact(name, f"WRF r8 {'dm1' if k == 'r8_dm1' else best['r8']} (this runtime) vs WRF r8 reference "
                              "(suite, Intel)", wrf_on[k], base, nl, subset=True)
        if "r8_dm1" in wrf_on:
            cmp_exact(name, "cpu_fp64_hwlog vs WRF r8 dm1 (this runtime)", R + "/out_cpu_fp64_hwlog", wrf_on["r8_dm1"], nl)
        cmp_exact(name, "gpu_fp64 vs cpu_fp64", R + "/out_gpu_fp64_on", R + "/out_cpu_fp64", nl)
        if "r8" in wrf_on:
            cmp_exact(name, "cpu_fp64_hwlog vs WRF r8 (this runtime)", R + "/out_cpu_fp64_hwlog", wrf_on["r8"], nl)
            cmp_exact(name, "cpu_fp64 vs WRF r8 (this runtime)", R + "/out_cpu_fp64", wrf_on["r8"], nl)
            cmp_l2(name, "gpu_fp32 vs WRF r8", R + "/out_gpu_fp32_on", wrf_on["r8"], nl)
        if "r4" in wrf_on:
            cmp_l2(name, "gpu_fp32 vs WRF r4", R + "/out_gpu_fp32_on", wrf_on["r4"], nl)
            if "r8" in wrf_on:
                cmp_l2(name, "WRF r4 vs WRF r8", wrf_on["r4"], wrf_on["r8"], nl)
    shutil.rmtree(f"{W}/runs/{name}", ignore_errors=True)
    log(f"== case {name} done")


def make_large(base, met):
    """W/run_L: BASE's tables and namelist with the met_em files' domain size; real.exe (real*8) makes the inputs."""
    import netCDF4
    R = W + "/run_L"
    shutil.rmtree(R, ignore_errors=True)
    os.makedirs(R)
    for f in os.listdir(base):
        p = os.path.join(base, f)
        if os.path.isdir(p) or re.match(r"(wrfout.*|wrfrst.*|wrfinput.*|wrfbdy.*|wrflowinp.*|rsl\..*|namelist\..*|"
                                        r".*\.log|.*\.exe)$", f):
            continue
        os.symlink(os.path.realpath(p), os.path.join(R, f))
    mets = sorted(glob.glob(met + "/met_em.d01.*"))
    for f in mets:
        os.symlink(f, os.path.join(R, os.path.basename(f)))
    ds = netCDF4.Dataset(mets[0])
    we, sn = ds.getncattr("WEST-EAST_GRID_DIMENSION"), ds.getncattr("SOUTH-NORTH_GRID_DIMENSION")
    nlev = len(ds.dimensions["num_metgrid_levels"])
    nsoil = int(ds.getncattr("NUM_METGRID_SOIL_LEVELS"))
    ds.close()
    nl = open(os.path.join(base, "namelist.input")).read()
    for k, v in (("e_we", we), ("e_sn", sn), ("num_metgrid_levels", nlev), ("num_metgrid_soil_levels", nsoil)):
        nl = setkey(nl, "domains", k, v)
    for k, v in (("start_hour", "00"), ("end_day", "12"), ("end_hour", "06"), ("interval_seconds", "21600")):
        nl = setkey(nl, "time_control", k, v)
    open(R + "/namelist.input", "w").write(nl)
    os.symlink(W + "/WRF-r8/main/real.exe", R + "/real.exe")
    topology()
    rc, wall = sh(f"exec mpiexec -n 16 {bind(16)} ./real.exe", R, out="real.log")
    text = "".join(open(f, errors="replace").read() for f in glob.glob(R + "/rsl.error.0000"))
    ok = rc == 0 and "SUCCESS COMPLETE REAL_EM" in text and os.path.exists(R + "/wrfbdy_d01")
    os.makedirs(LOGS + "/L", exist_ok=True)
    for f in ("namelist.input", "rsl.error.0000", "real.log"):
        if os.path.exists(os.path.join(R, f)):
            shutil.copy(os.path.join(R, f), f"{LOGS}/L/real.{f}")
    for f in glob.glob(R + "/rsl.*") + glob.glob(R + "/met_em.*") + [R + "/real.exe"]:
        os.remove(f)
    emit(dict(kind="real", case="L", ok=ok, wall=wall, e_we=int(we), e_sn=int(sn), levels=int(nlev)))
    log(f"L: real.exe {'ok' if ok else 'FAILED'} ({wall:.0f} s), {we} x {sn}, {nlev} metgrid levels")
    return 0 if ok else 1


def main():
    a = sys.argv[1:]
    if a and a[0] == "make-large":
        sys.exit(make_large(a[1], a[2]))
    if a and a[0] == "report":
        sys.path.insert(0, SRC + "/tools/bench")
        import bench_report
        open(J + "/report.md", "w").write(bench_report.report(RES, logs=LOGS))
        return
    if not a or a[0] != "case":
        print(__doc__)
        sys.exit(2)
    name, base = a[1], a[2]
    ref, sets, layouts, serial = False, [], "dm12,dm24,dm48,sm24x2,sm12x4", False
    i = 3
    while i < len(a):
        if a[i] == "--ref":
            ref = True
            i += 1
        elif a[i] == "--serial":
            serial = True
            i += 1
        elif a[i] == "--set":
            sets.append(a[i + 1])
            i += 2
        elif a[i] == "--layouts":
            layouts = a[i + 1]
            i += 2
        else:
            raise SystemExit(f"unknown argument {a[i]}")
    topology()
    log(f"topology: {NP} vCPUs, {CORES} cores, binding {'on' if BIND else 'off'}")
    case(name, base, sets, ref, [x for x in layouts.split(",") if x], serial)


if __name__ == "__main__":
    main()
