#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""The runs of tools/bench/colab_fp32.sh (the FP32 accuracy study, docs/FP32_ACCURACY.md) on a Colab G4 runtime.

  fp32_run.py prepare [CASE,...]  (chroot) the run directories /opt/fp1/runs/CASE/{wrf_*,cpp} with their namelists,
                                  CPPWRF's configuration (tools/wrf2cppwrf.py) and the GPU queue /opt/fp1/gpu_queue.txt
  fp32_run.py wrf                 (chroot) the WRF runs of the prepared cases, several at a time (at most SLOTS
                                  processes: FP1_SLOTS, default 2/3 of the vCPUs, at most 32, or all of 16 or fewer;
                                  an MPI run takes at most half of them, a third of 16 or fewer; the real*8 runs only
                                  with 32 vCPUs or more), one JSON line each in /opt/fp1/out/runs.jsonl;
                                  /opt/fp1/wrf_done at the end
  fp32_run.py gpu BUILD_DIR       (runtime) the GPU queue: BUILD_DIR/build-cuda (FP32) and build-cuda_dp (FP64)
  fp32_run.py stats               (runtime or chroot) tools/verify/fp32_accuracy.py compare for every pair of every
                                  case (/opt/fp1/out/stats/*.jsonl), then report (/opt/fp1/out/report_*.md, svg/)

A case is a WRF real*8 run directory of the suite's data tars (BASE: its inputs, tables and namelist) and a length.
Every model runs BASE's namelist with that length, history every HIST minutes, restart off:
- WRF v4.8.0 real*4 (the bench's dmpar build /opt/wrfbench/WRF-r4) in one process; for the long cases, which one
  process cannot finish within the session, with MPI ranks (and P4L also 24 h in one process, to show what the layout
  changes); WRF real*8 (WRF-r8) with MPI ranks on the long cases ("r8n": a real*8 run whose last bits differ from the
  suite's one-process Intel references, through the layout and this CPU's libmvec log);
- CPPWRF's CUDA FP32 and FP64 release builds from tools/wrf2cppwrf.py's translation of the same namelist. The FP64
  run is the real*8 reference of the comparisons: it is bit-identical to WRF real*8 (one process, Intel), which the
  suite checks at tolerance 0 (the fp64-ref pair repeats it against the reference files of BASE where they exist).
The real*4 runs read WRF's single-precision RRTMG tables (as tools/bench/bench_run.py) and compute their own Thompson
tables (the run directories' .dat files are real*8 runs' files).
"""
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import time

F = "/opt/fp1"
J = "/opt/job"
W = "/opt/wrfbench"
SRC = J + "/src"
WR = J + "/data/wrf_runs"
OUT = F + "/out"
NCPU = os.cpu_count() or 4
SLOTS = int(os.environ.get("FP1_SLOTS", str(min(32, NCPU * 2 // 3) if NCPU > 16 else NCPU)))
CAP = SLOTS // 2 if NCPU > 16 else max(1, SLOTS // 3)  # the ranks of one MPI run
sys.path.insert(0, SRC + "/tools/bench")

# name: (BASE, hours, history minutes, extra namelist settings, WRF runs [(tag, precision, ranks, hours or None)])
R4 = [("r4", "r4", 1, None)]
CASES = {
    "P4": ("run_P4_r8", 1, 10, [], R4 + [("r4m", "r4", 8, None)]),
    "P6": ("run_P6_r8", 1, 10, [], R4),
    "P9": ("run_P9_r8", 1, 10, [], R4),
    "N3": ("run_N3h_r8", 1, 10, [], R4 + [("r4m", "r4", 12, None)]),
    "V1": ("run_V1_r8", 1, 10, [], R4),
    "T1": ("run_T1_r8", 1, 10, [], R4),
    "SKV": ("run_SKV_r8", 1, 10, [], R4),
    "MORR": ("run_MORR_r8", 1, 10, [], R4),
    "ACM7": ("run_ACM7_r8", 1, 10, [], R4),
    "L": (W + "/run_L", 1, 10, ["physics.radt=15"], R4),
    "P4L": ("run_P4_72h_r8", 72, 60, [], [("r4", "r4", 8, None), ("r4s", "r4", 1, 24), ("r8n", "r8", 8, None)]),
    "N3L": ("run_N3_24h_r8", 24, 60, [], [("r4", "r4", 12, None), ("r8n", "r8", 12, None)]),
}
# the FP32 outliers of the 1 h study (ACM7, MORR, SKV) and their control P4 again, from the real*8 real.exe inputs of
# fp2_data.tar (24 h boundaries; wrfinput identical to run_P4_r8's / run_V1_r8's): "S" = 1 h with an ensemble of
# equally valid single-precision runs (WRF real*4 in one process and with 2, 4 and 8 MPI ranks; CPPWRF GPU FP32, GPU
# FP32 without FMA contraction, CPU FP32), "L" = 24 h (WRF real*4 and real*8 with 8 ranks)
_ACM7 = ["physics.bl_pbl_physics=7", "physics.sf_sfclay_physics=7"]
_MORR = ["physics.mp_physics=10"]
_SKV = ["physics.sst_skin=1", "physics.tmn_update=1"]
_R4E = R4 + [("r4m2", "r4", 2, None), ("r4m4", "r4", 4, None), ("r4m8", "r4", 8, None)]
_R4L = [("r4", "r4", 8, None), ("r8n", "r8", 8, None)]
for _n, _b, _s in (("P4", "run_FP2_P4_24h_r8", []), ("ACM7", "run_FP2_P4_24h_r8", _ACM7),
                   ("MORR", "run_FP2_P4_24h_r8", _MORR), ("SKV", "run_FP2_V1_24h_r8", _SKV)):
    CASES[_n + "S"] = (_b, 1, 10, _s, _R4E)
    if _n != "P4":
        CASES[_n + "L"] = (_b, 24, 60, _s, _R4L)
ORDER = ["N3L", "P4L", "ACM7L", "MORRL", "SKVL", "L", "N3", "P4", "P6", "P9", "V1", "T1", "SKV", "MORR", "ACM7",
         "SKVS", "ACM7S", "MORRS", "P4S"]
# CPPWRF runs per case: (tag, build); FP2's cases add GPU FP32 without FMA (cuda_nofma) and, for 1 h, CPU FP32 (cpu)
CPP = {c: [("fp32", "cuda"), ("fp64", "cuda_dp")] for c in CASES}
for _c in CASES:
    if _c.endswith("S"):
        CPP[_c] += [("nofma", "cuda_nofma"), ("cpu32", "cpu")]
    elif _c in ("ACM7L", "MORRL", "SKVL"):
        CPP[_c] += [("nofma", "cuda_nofma")]


def log(msg):
    os.makedirs(OUT, exist_ok=True)
    with open(OUT + "/fp1.log", "a") as f:
        f.write(time.strftime("%H:%M:%S ") + msg + "\n")


def setkey(nl, sec, key, val):
    m = re.search(rf"&{sec}\b(.*?)\n\s*/", nl, re.S | re.I)
    if not m:
        return nl.rstrip("\n") + f"\n\n&{sec}\n {key} = {val},\n/\n"
    body = m.group(1)
    pat = re.compile(rf"^\s*{key}\s*=.*$", re.M | re.I)
    line = f" {key} = {val},"
    body = pat.sub(line, body) if pat.search(body) else body + "\n" + line
    return nl[:m.start(1)] + body + nl[m.end(1):]


def geti(nl, key, d=0, default=0):
    m = re.search(rf"^\s*{key}\s*=\s*([^\n]*)$", nl, re.M | re.I)
    v = [x.strip() for x in m.group(1).split(",") if x.strip()] if m else []
    return int(float(v[d])) if len(v) > d else default


def base_dir(b):
    return b if b.startswith("/") else os.path.join(WR, b)


def make_dir(case, name, hours=None):
    """F/runs/CASE/NAME: BASE's inputs (links) and its namelist with the case's length and history interval."""
    b, h0, hist, sets, _ = CASES[case]
    hours = hours or h0
    base = base_dir(b)
    R = f"{F}/runs/{case}/{name}"
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
              "history_interval_s", "history_begin", "history_begin_h", "history_begin_m"):
        nl = re.sub(rf"^\s*{k}\s*=.*\n", "", nl, flags=re.M | re.I)
    per = lambda v: ", ".join([str(v)] * nd)  # noqa: E731
    for k, v in (("run_days", "0"), ("run_hours", str(hours)), ("run_minutes", "0"), ("run_seconds", "0"),
                 ("restart", ".false."), ("restart_interval", "100000"), ("frames_per_outfile", per(1000)),
                 ("history_interval", per(hist)), ("debug_level", "0")):
        nl = setkey(nl, "time_control", k, v)
    for s in sets:
        sk, v = s.split("=", 1)
        sec, k = sk.split(".", 1)
        nl = setkey(nl, sec, k, v)
    open(os.path.join(R, "namelist.input"), "w").write(nl)
    return R


def r4_tables(R):
    """WRF real*4 reads its own single-precision RRTMG / RRTM tables (BASE's are the *_DATA_DBL files of a real*8 run)
    and computes the Thompson tables itself."""
    run = f"{W}/WRF-r4/run"
    for f in os.listdir(run):
        if os.path.exists(f"{run}/{f}_DBL") and os.path.lexists(f"{R}/{f}"):
            os.remove(f"{R}/{f}")
            os.symlink(os.path.realpath(f"{run}/{f}"), f"{R}/{f}")
    for f in ("freezeH2O.dat", "qr_acr_qg_V4.dat", "qr_acr_qsV2.dat", "qr_acr_qg.dat", "qr_acr_qs.dat"):
        if os.path.lexists(f"{R}/{f}"):
            os.remove(f"{R}/{f}")


def prepare(cases):
    os.makedirs(OUT, exist_ok=True)
    plan, queue = [], []
    for c in cases:
        b = base_dir(CASES[c][0])
        if not os.path.exists(os.path.join(b, "namelist.input")):
            log(f"{c}: no inputs ({b})")
            continue
        for tag, prec, ranks, hours in CASES[c][4]:
            if prec == "r8" and NCPU < 32:  # the real*8 MPI runs only where the CPUs leave room for them
                continue
            R = make_dir(c, f"wrf_{tag}", hours)
            if prec == "r4":
                r4_tables(R)
            os.symlink(f"{W}/WRF-{prec}/main/wrf.exe", R + "/wrf.exe")
            plan.append(dict(case=c, tag=tag, prec=prec, ranks=min(ranks, CAP), dir=R,
                             hours=hours or CASES[c][1]))
        R = make_dir(c, "cpp")
        t = subprocess.run([sys.executable, SRC + "/tools/wrf2cppwrf.py", R + "/namelist.input", "--real-dir", R,
                            "--format", "wrfout", "-o", R + "/case.ini"], capture_output=True, text=True)
        open(R + "/translate.log", "w").write(t.stdout + t.stderr)
        if t.returncode:
            log(f"{c}: wrf2cppwrf.py FAILED: {(t.stdout + t.stderr).strip().splitlines()[-1:]}")
            continue
        for tag, build in CPP[c]:
            queue.append(f"{c} {tag} {build} {R} {R}/case.ini {R}/out_{tag}")
        log(f"{c}: prepared ({len(CASES[c][4])} WRF runs, CPPWRF FP32 and FP64)")
    # the real*4 runs of the long cases first (the session's critical path: P4L's 24 h in one process), then the
    # one-hour cases (the longest first), then the real*8 runs with MPI
    key = lambda p: (p["prec"] == "r8", -p["hours"] * (6 if p["ranks"] == 1 else 1), ORDER.index(p["case"]))  # noqa: E731
    json.dump(sorted(plan, key=key), open(F + "/wrf_plan.json", "w"), indent=1)
    queue.sort(key=lambda q: (ORDER.index(q.split()[0]), q.split()[1]))
    open(F + "/gpu_queue.txt", "w").write("\n".join(queue) + "\n")


def wrf_timings(R):
    p = R + "/rsl.out.0000"
    text = open(p, errors="replace").read() if os.path.exists(p) else ""
    main = [float(m.group(2)) for m in re.finditer(r"Timing for main: time \S+ on domain\s+(\d+):\s+([0-9.]+) elapsed", text)
            if m.group(1) == "1"]
    return "SUCCESS COMPLETE WRF" in text, len(main), sum(main)


def wrf():
    plan = json.load(open(F + "/wrf_plan.json"))
    pend, running, used = list(plan), [], 0
    env = dict(os.environ, OMP_NUM_THREADS="1")
    while pend or running:
        for p in list(pend):
            if used + p["ranks"] <= SLOTS or not running:
                cmd = f"ulimit -s unlimited; exec mpiexec -n {p['ranks']} ./wrf.exe > wrf.log 2>&1"
                p["t0"] = time.time()
                p["proc"] = subprocess.Popen(cmd, shell=True, cwd=p["dir"], env=env)
                running.append(p)
                pend.remove(p)
                used += p["ranks"]
                log(f"{p['case']} wrf_{p['tag']}: started ({p['ranks']} processes, {p['hours']} h; {used} of {SLOTS} slots)")
        time.sleep(10)
        for p in list(running):
            rc = p["proc"].poll()
            if rc is None:
                continue
            running.remove(p)
            used -= p["ranks"]
            ok, steps, step_s = wrf_timings(p["dir"])
            wall = time.time() - p["t0"]
            rec = dict(case=p["case"], tag=p["tag"], prec=p["prec"], ranks=p["ranks"], hours=p["hours"], exit=rc,
                       ok=ok and rc == 0, wall=round(wall, 1), steps=steps, main_s=round(step_s, 1))
            with open(OUT + "/runs.jsonl", "a") as f:
                f.write(json.dumps(rec) + "\n")
            os.makedirs(f"{OUT}/logs/{p['case']}", exist_ok=True)
            for x in ("rsl.out.0000", "rsl.error.0000", "wrf.log", "namelist.input"):
                if os.path.exists(f"{p['dir']}/{x}"):
                    os.system(f"tail -c 200000 {p['dir']}/{x} > {OUT}/logs/{p['case']}/wrf_{p['tag']}.{x}")
            log(f"{p['case']} wrf_{p['tag']}: {'ok' if rec['ok'] else 'FAILED'} ({wall:.0f} s)")
    open(F + "/wrf_done", "w").write("done\n")


def gpu(bdir):
    for line in open(F + "/gpu_queue.txt"):
        f = line.split()
        if len(f) != 6:
            continue
        c, tag, build, R, ini, out = f
        shutil.rmtree(out, ignore_errors=True)
        t0 = time.time()
        with open(f"{R}/{tag}.log", "w") as lf:
            if not os.path.exists(f"{bdir}/build-{build}/cppwrf"):
                lf.write(f"no build {build}\n")
                rc = 127
            else:
                rc = subprocess.run([f"{bdir}/build-{build}/cppwrf", ini, "--set", f"output.dir={out}"], cwd=R,
                                    stdout=lf, stderr=subprocess.STDOUT,
                                    env=dict(os.environ, OMP_NUM_THREADS="8" if build == "cpu" else "4")).returncode
        wall = time.time() - t0
        text = open(f"{R}/{tag}.log", errors="replace").read()
        done = re.findall(r"cppwrf: done, (\d+) steps, stepping ([0-9.]+) s", text)
        rec = dict(case=c, tag=tag, exit=rc, ok=rc == 0 and bool(done), wall=round(wall, 1),
                   stepping=float(done[0][1]) if done else None)
        with open(OUT + "/runs.jsonl", "a") as f:
            f.write(json.dumps(rec) + "\n")
        os.makedirs(f"{OUT}/logs/{c}", exist_ok=True)
        os.system(f"tail -c 100000 {R}/{tag}.log > {OUT}/logs/{c}/cpp_{tag}.log; cp {R}/case.ini {R}/translate.log "
                  f"{OUT}/logs/{c}/ 2> /dev/null")
        log(f"{c} {tag}: {'ok' if rec['ok'] else 'FAILED'} ({wall:.0f} s)")
    open(F + "/gpu_done", "w").write("done\n")


def domains(case):
    nl = open(os.path.join(base_dir(CASES[case][0]), "namelist.input")).read()
    return [f"d{i:02d}" for i in range(1, geti(nl, "max_dom", 0, 1) + 1)]


def stats():
    tool = SRC + "/tools/verify/fp32_accuracy.py"
    os.makedirs(OUT + "/stats", exist_ok=True)
    jobs = []
    for c in ORDER:
        R = f"{F}/runs/{c}"
        if not os.path.isdir(R):
            continue
        fp32, fp64 = f"{R}/cpp/out_fp32", f"{R}/cpp/out_fp64"
        pairs = [("fp32-r8", fp32, fp64)]
        for tag, _ in CPP[c][2:]:  # the other CPPWRF single-precision runs
            pairs += [(f"{tag}-r8", f"{R}/cpp/out_{tag}", fp64), (f"{tag}-r4", f"{R}/cpp/out_{tag}", f"{R}/wrf_r4")]
        for tag, prec, ranks, hours in CASES[c][4]:
            w = f"{R}/wrf_{tag}"
            if prec == "r4":
                n = "r4" if tag == "r4" else tag
                pairs += [(f"{n}-r8", w, fp64)] + ([("fp32-r4", fp32, w)] if tag == "r4" else [])
            else:
                pairs.append((f"{tag}-r8", w, fp64))
        base = base_dir(CASES[c][0])
        if glob.glob(base + "/wrfout_d01_*") and not c.endswith("L"):
            pairs.append(("fp64-ref", fp64, base))
        for d in domains(c):
            for name, a, b in pairs:
                jobs.append([sys.executable, tool, "compare", a, b, "--domain", d, "--case", c, "--pair", name,
                             "--out", f"{OUT}/stats/{c}_{d}_{name}.jsonl"])
    procs = []
    for j in jobs:  # 12 at a time
        while sum(p.poll() is None for p in procs) >= 12:
            time.sleep(1)
        procs.append(subprocess.Popen(j, stdout=subprocess.DEVNULL, stderr=open(OUT + "/stats.err", "a")))
    for p in procs:
        p.wait()
    files = sorted(glob.glob(OUT + "/stats/*.jsonl"))
    for m in ("rmse", "max_abs", "rel_l2"):
        subprocess.run([sys.executable, tool, "report", *files, "--metric", m, "--md", f"{OUT}/report_{m}.md"]
                       + (["--svg-dir", OUT + "/svg", "--svg-rel", "img"] if m == "rmse" else []))
    log(f"stats: {len(jobs)} comparisons, {len(files)} files")


def main():
    a = sys.argv[1:]
    if not a:
        print(__doc__)
        sys.exit(2)
    if a[0] == "prepare":
        prepare(a[1].split(",") if len(a) > 1 else ORDER)
    elif a[0] == "wrf":
        wrf()
    elif a[0] == "gpu":
        gpu(a[1])
    elif a[0] == "stats":
        stats()
    else:
        print(__doc__)
        sys.exit(2)


if __name__ == "__main__":
    main()
