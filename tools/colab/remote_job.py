#!/usr/bin/env python3
"""The Colab VM side of tools/colab/colab_run.sh (runs detached under nohup; /content/job/done marks the end).

Expects /content/src.tar.gz (git archive of the repository) and /content/data.tar (tools/colab/stage_data.py),
and the environment CPPWRF_JOB_CASES (wrf_suite.py case names, comma-separated), CPPWRF_JOB_FILES (cases of
tests/verify/run_wrf_files_test.sh, may be empty) and CPPWRF_JOB_BENCH (suite case names to time again without
output, may be empty; NAME@cuda or NAME@cuda_dp runs the case with that build). Builds the CUDA FP32 and FP64
binaries for the VM's GPU, runs the suite cases against the WRF references, the files test with the FP64 binary
(CPPWRF_FILES_NO_VALUES=1), times the bench cases, and leaves everything in /content/job (summary.txt, logs) and
/content/results.tar.gz.

Profiling (colab_run.sh --env / --nsys / --ncu):
- CPPWRF_JOB_ENV_KEYS: the variables of the environment that colab_run.sh --env set (e.g. CPPWRF_KTIME). The suite
  runs with them; each bench case runs once without them (the stepping time) and once more with them.
- CPPWRF_JOB_NSYS=1: each bench case also runs under nsys for its whole time and for half of it; the difference
  of the two gives the time loop per step (kernel launches and time, copies, synchronisations, allocations), the
  whole run the top kernels (profile/NAME.txt).
- CPPWRF_JOB_NCU=REGEX (demangled kernel names) or "top" (the nsys top kernels): ncu on those kernels' launches
  in the first 3 steps (achieved occupancy, local-memory traffic, DRAM bytes), joined with the registers and stack
  of every kernel (cuobjdump -res-usage).
The bench run without the extra variables also samples the GPU's used memory (nvidia-smi every 0.2 s): the peak
minus the idle value is the device memory the run reserves (cudaMemGetInfo's total - free of the only process).
The profiler runs leave out CPPWRF_KTIME, whose synchronisation after every launch would be in the trace.
"""
import csv
import glob
import os
import re
import shlex
import subprocess
import sys
import threading
import time

JOB = "/content/job"
SRC = "/content/cppwrf"
WR = "/content/data/wrf_runs"


def sh(cmd, log, cwd=None, env=None):
    with open(os.path.join(JOB, log), "a") as f:
        f.write(f"$ {cmd}\n")
        f.flush()
        t = time.time()
        r = subprocess.run(cmd, shell=True, cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT)
        f.write(f"[exit {r.returncode}, {time.time() - t:.0f} s]\n")
    return r.returncode


def out(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip()


def timings(log):
    """(steps, stepping seconds, ns/cell/step) of each domain's "cppwrf: done" line in the log."""
    return [(int(k), s, n) for k, s, n in
            re.findall(r"cppwrf: done, (\d+) steps, stepping ([0-9.]+) s, ([0-9.]+) ns/cell/step", open(log).read())]


class GpuMemory:
    """The GPU's peak used memory (nvidia-smi every 0.2 s) while the block runs, and the value before it."""

    def __enter__(self):
        q = ["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"]
        self.idle = int((out(shlex.join(q)).split() or ["0"])[0])
        self.vals = []
        self.p = subprocess.Popen(q + ["-lms", "200"], stdout=subprocess.PIPE, text=True)
        self.th = threading.Thread(target=self.read, daemon=True)
        self.th.start()
        return self

    def read(self):
        for line in self.p.stdout:
            if line.strip().isdigit():
                self.vals.append(int(line))

    def __exit__(self, *exc):
        self.p.terminate()
        self.th.join(5)

    def text(self):
        if not self.vals:
            return "? (no nvidia-smi samples)"
        return f"{max(self.vals) - self.idle} MiB (peak {max(self.vals)} MiB, idle {self.idle} MiB)"


def base(name):
    """A demangled kernel name without "void " and its parameter list (templates kept)."""
    name = name[5:] if name.startswith("void ") else name
    depth = 0
    for i, ch in enumerate(name):
        depth += (ch == "<") - (ch == ">")
        if ch == "(" and depth == 0 and not name.startswith("(anonymous", i):
            return name[:i]
    return name


def key(name):  # nsys, ncu and cu++filt space template arguments and print integer arguments differently
    k = re.sub(r"\s+", "", base(name).replace("anonymous namespace", "anonymous_namespace"))
    k = k.replace("(bool)1", "true").replace("(bool)0", "false")
    return re.sub(r"\((unsigned)?(int|long|short|char)\)(-?\d+)", r"\3", k)


def key_regex(k):
    """An ncu -k regex for key k: its launcher and functor, whatever the spelling of the template arguments."""
    m = re.match(r"(.*?<)([^<>,]+)", k)
    rx = re.escape(m.group(1)) + re.escape(m.group(2)) + r"\s*[<>,]" if m else re.escape(k)
    return rx.replace("anonymous_namespace", r"anonymous\snamespace")


def short(name):
    return base(name).replace("cppwrf::detail::", "").replace("cppwrf::", "")


def res_usage(binary):
    """key -> (registers, stack bytes, local bytes, shared bytes) of every kernel in the binary (cuobjdump)."""
    cur, raw = None, {}
    for line in out(f"cuobjdump -res-usage {binary}").splitlines():
        m = re.match(r"\s*Function (\S+):", line)
        if m:
            cur = m.group(1)
        v = dict(re.findall(r"\b(REG|STACK|SHARED|LOCAL):(\d+)", line))
        if cur and "REG" in v:
            old = raw.get(cur, (0, 0, 0, 0))
            raw[cur] = tuple(max(o, int(v.get(f, 0))) for o, f in zip(old, ("REG", "STACK", "LOCAL", "SHARED")))
            cur = None
    names = list(raw)
    try:
        dem = subprocess.run(["cu++filt"], input="\n".join(names) + "\n", capture_output=True,
                             text=True).stdout.splitlines()
    except OSError:
        dem = []
    if len(dem) != len(names):
        dem = names
    res = {key(d): raw[n] for n, d in zip(names, dem)}
    with open(os.path.join(JOB, f"res-usage-{os.path.basename(os.path.dirname(binary))}.txt"), "w") as f:
        f.write("REG STACK LOCAL SHARED kernel\n")
        for n, d in sorted(zip(names, dem), key=lambda x: -raw[x[0]][1]):
            f.write("%3d %6d %6d %6d %s\n" % (*raw[n], short(d)))
    return res


PROF = "/content/prof"
NSYS_REPORTS = ("cuda_gpu_kern_sum", "cuda_api_sum", "cuda_gpu_mem_size_sum")


def nsys_run(tag, part, cmd, run_time, env):
    """nsys of the case for run_time; (steps of the first domain, stepping seconds), {report: rows}."""
    rep, o = f"{PROF}/{tag}_{part}", os.path.join(JOB, "nsys", f"{tag}_{part}")
    os.makedirs(PROF, exist_ok=True)
    os.makedirs(os.path.dirname(o), exist_ok=True)
    log, d = f"nsys-{tag}-{part}.log", f"/content/bench/{tag}_nsys"
    c = cmd(d, run_time)
    sh(f"rm -rf {d} && nsys profile --trace=cuda --sample=none --cpuctxsw=none --force-overwrite=true -o {rep} "
       f"{shlex.join(c)}", log, cwd=SRC, env=env)
    sh(f"nsys stats -q -f csv -r {','.join(NSYS_REPORTS)} -o {o} {rep}.nsys-rep", log)
    for f in glob.glob(rep + ".*"):  # the reports are large: only the tables come back
        os.remove(f)
    rows = {}
    for r in NSYS_REPORTS:
        fs = glob.glob(f"{o}*{r}*.csv")
        rows[r] = list(csv.DictReader(open(fs[0]))) if fs else []
    done = timings(os.path.join(JOB, log))
    return (done[0][0], float(done[0][1])) if done else (0, 0.0), rows


def profile_nsys(tag, cmd, t, env, note):
    """The time loop per step from the whole and the half run; the top 15 kernels (keys) for ncu."""
    (sf, wf), f = nsys_run(tag, "full", cmd, t, env)
    (sh_, wh), h = nsys_run(tag, "half", cmd, t // 2, env)
    ds = sf - sh_
    if ds <= 0:
        note(f"nsys {tag}: no steps (see nsys-{tag}-*.log)")
        return []

    def table(rows, name, cols):
        d = {}
        for r in rows:
            v = d.setdefault(r.get(name, "?"), [0.0] * len(cols))
            for i, c in enumerate(cols):
                v[i] += float(r.get(c) or 0)
        return d

    kc = ("Total Time (ns)", "Instances")
    kf = table(f["cuda_gpu_kern_sum"], "Name", kc)
    kh = table(h["cuda_gpu_kern_sum"], "Name", kc)
    ac = ("Num Calls", "Total Time (ns)")
    af, ah = table(f["cuda_api_sum"], "Name", ac), table(h["cuda_api_sum"], "Name", ac)
    mc = ("Count", "Total (MB)")
    mf, mh = table(f["cuda_gpu_mem_size_sum"], "Operation", mc), table(h["cuda_gpu_mem_size_sum"], "Operation", mc)

    def per_step(a, b, n, i):
        return (a.get(n, [0, 0])[i] - b.get(n, [0, 0])[i]) / ds

    ks = {}  # kernels merged by key (one functor can be launched from several translation units)
    for n in kf:
        e = ks.setdefault(key(n), [n, 0.0, 0.0, 0.0, 0.0])
        e[1] += kf[n][0]
        e[2] += kf[n][1]
        e[3] += per_step(kf, kh, n, 0)
        e[4] += per_step(kf, kh, n, 1)
    launches = sum(e[4] for e in ks.values())
    kms = sum(e[3] for e in ks.values()) / 1e6
    top = sorted((x for x in ks.items() if x[1][4] > 0), key=lambda x: -x[1][3])[:15]
    api = sorted(((n, per_step(af, ah, n, 0), per_step(af, ah, n, 1) / 1e6) for n in af), key=lambda x: -x[1])
    mem = [(n, per_step(mf, mh, n, 0), per_step(mf, mh, n, 1)) for n in mf]
    os.makedirs(os.path.join(JOB, "profile"), exist_ok=True)
    with open(os.path.join(JOB, "profile", f"{tag}.txt"), "w") as o:
        o.write(f"{tag}: nsys runs of {sf} and {sh_} steps (first domain), stepping {wf:.3f} s and {wh:.3f} s under "
                f"nsys; per step of the time loop (difference / {ds} steps):\n")
        o.write(f"  wall {1e3 * (wf - wh) / ds:.3f} ms, kernels {kms:.3f} ms, {launches:.1f} kernel launches\n")
        o.write("  CUDA API calls per step (calls, ms):\n")
        for n, c, ms in api:
            if abs(c) >= 0.05:
                o.write(f"    {c:10.2f} {ms:9.3f}  {n}\n")
        o.write("  device memory operations per step (count, MB):\n")
        for n, c, mb in mem:
            o.write(f"    {c:10.2f} {mb:10.4f}  {n}\n")
        o.write("  top 15 kernels by time per step: ms/step, launches/step, us/launch, whole run ms, launches\n")
        for k, (n, tt, ni, dt, dn) in top:
            o.write(f"    {dt / 1e6:8.3f} {dn:7.1f} {dt / dn / 1e3 if dn else 0:9.1f} {tt / 1e6:10.1f} {ni:8.0f}  "
                    f"{short(n)}\n")
    calls = {n: c for n, c, _ in api}
    syncs = sum(c for n, c in calls.items() if "Synchronize" in n)
    allocs = sum(c for n, c in calls.items() if re.search(r"Malloc|Free|HostAlloc", n))
    copies = " ".join(f"{n.strip('[]').replace('CUDA ', '')} {c:.1f} ({mb:.3f} MB)" for n, c, mb in mem if abs(c) >= 0.05)
    note(f"nsys {tag}: per step {1e3 * (wf - wh) / ds:.2f} ms wall, {kms:.2f} ms kernels, {launches:.0f} launches, "
         f"{syncs:.1f} syncs, {allocs:.1f} alloc/free; copies: {copies or 'none'}")
    return [k for k, _ in top]


NCU_METRICS = ["gpu__time_duration.sum", "sm__warps_active.avg.pct_of_peak_sustained_active",
               "sm__maximum_warps_per_active_cycle_pct", "launch__registers_per_thread",
               "l1tex__t_bytes_pipe_lsu_mem_local_op_ld.sum", "l1tex__t_bytes_pipe_lsu_mem_local_op_st.sum",
               "dram__bytes.sum", "sm__throughput.avg.pct_of_peak_sustained_elapsed"]


def profile_ncu(tag, cmd, run_time, env, ncu, top, res, note):
    """ncu on the chosen kernels' launches of the first steps, aggregated per kernel (profile/NAME_ncu.txt)."""
    if ncu == "top":
        if not top:
            note(f"ncu {tag}: --ncu top needs the nsys top kernels")
            return
        rx = "|".join(key_regex(k) for k in top)
    else:
        rx = ncu
    rep, d = f"{PROF}/{tag}_ncu", f"/content/bench/{tag}_ncu"
    os.makedirs(PROF, exist_ok=True)
    for metrics in (NCU_METRICS, NCU_METRICS[:4]):  # the smaller set when a metric is unknown to this GPU
        for f in glob.glob(rep + ".*"):
            os.remove(f)
        sh(f"rm -rf {d} && ncu --kernel-name-base demangled -k {shlex.quote('regex:' + rx)} --launch-count 300 "
           f"--metrics {','.join(metrics)} -f -o {rep} {shlex.join(cmd(d, max(1, int(run_time))))}",
           f"ncu-{tag}.log", cwd=SRC, env=env)
        if os.path.exists(rep + ".ncu-rep"):
            break
    raw = os.path.join(JOB, "profile", f"{tag}_ncu.csv")
    os.makedirs(os.path.dirname(raw), exist_ok=True)
    sh(f"ncu --import {rep}.ncu-rep --csv --page raw --print-units base --print-kernel-base demangled > {raw}",
       f"ncu-{tag}.log")
    rows = list(csv.reader(open(raw))) if os.path.exists(raw) else []
    if len(rows) < 3 or "Kernel Name" not in rows[0]:
        note(f"ncu {tag}: no kernels profiled (see ncu-{tag}.log)")
        return
    hd = rows[0]
    agg = {}
    for r in rows[2:]:  # rows[1]: units
        d = dict(zip(hd, r))
        e = agg.setdefault(key(d["Kernel Name"]), {"name": d["Kernel Name"], "n": 0})
        e["n"] += 1
        for m in NCU_METRICS:
            try:
                e[m] = e.get(m, 0.0) + float(d.get(m, "").replace(",", ""))
            except ValueError:
                pass
    order = [k for k in top if k in agg] + sorted(k for k in agg if k not in top)
    with open(os.path.join(JOB, "profile", f"{tag}_ncu.txt"), "w") as o:
        o.write(f"{tag}: ncu, {run_time:.0f} s of run time, kernels {'nsys top 15' if ncu == 'top' else rx}; means per "
                "launch: launches, us, achieved / theoretical occupancy %, SM throughput %, registers (cuobjdump: "
                "stack bytes), local load / store KB, DRAM GB/s\n")
        for k in order:
            e, n = agg[k], agg[k]["n"]
            g = (lambda m: e.get(m, float("nan")) / n)
            r = res.get(k, (0, 0, 0, 0))
            us = g("gpu__time_duration.sum") / 1e3
            o.write(f"  {n:4d} {us:9.1f} {g(NCU_METRICS[1]):6.1f} {g(NCU_METRICS[2]):6.1f} {g(NCU_METRICS[7]):6.1f} "
                    f"{g(NCU_METRICS[3]):4.0f} ({r[1]:6d}) {g(NCU_METRICS[4]) / 1e3:10.1f} {g(NCU_METRICS[5]) / 1e3:10.1f} "
                    f"{g(NCU_METRICS[6]) / us / 1e3 if us else 0:8.1f}  {short(e['name'])}\n")
    note(f"ncu {tag}: {sum(e['n'] for e in agg.values())} launches of {len(agg)} kernels (profile/{tag}_ncu.txt)")


def main():
    os.makedirs(JOB, exist_ok=True)
    summary = open(os.path.join(JOB, "summary.txt"), "w")

    def note(s):
        summary.write(s + "\n")
        summary.flush()

    try:
        job(note)
    except Exception as e:  # noqa: BLE001 - colab_run.sh waits for the done file: always write it
        note(f"job error: {type(e).__name__}: {e}")
    summary.close()
    sh("tar czf /content/results.tar.gz -C /content job", "pack.log")
    open(os.path.join(JOB, "done"), "w").write("done\n")


def job(note):
    cases = os.environ.get("CPPWRF_JOB_CASES", "")
    files = [f for f in os.environ.get("CPPWRF_JOB_FILES", "").split(",") if f]
    bench = [b for b in os.environ.get("CPPWRF_JOB_BENCH", "").split(",") if b]
    gpu = out("nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader")
    cpu = out("grep -m1 'model name' /proc/cpuinfo")
    nvcc = out("nvcc --version | tail -2 | head -1")
    commit = open("/content/commit.txt").read().strip() if os.path.exists("/content/commit.txt") else "?"
    note(f"GPU: {gpu}")
    note(f"CPU: {out('nproc')} cores, {cpu}")
    note(f"nvcc: {nvcc}")
    note(f"commit: {commit}")

    sh("apt-get -qq update && apt-get -qq install -y libnetcdf-dev > /dev/null", "setup.log")
    sh("pip -q install netCDF4 numpy", "setup.log")
    sh(f"rm -rf {SRC} /content/data && mkdir -p {SRC} /content/data && tar xzf /content/src.tar.gz -C {SRC} "
       "&& tar xf /content/data.tar -C /content/data", "setup.log")
    nproc = int(out("nproc") or "4")
    bins = {}
    for name, dp in (("cuda", "OFF"), ("cuda_dp", "ON")):
        b = f"/content/build-{name}"
        t = time.time()
        rc = sh(f"cmake -S {SRC} -B {b} -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=ON -DCPPWRF_DOUBLE={dp} "
                f"&& cmake --build {b} -j{nproc} --target cppwrf", f"build-{name}.log")
        note(f"build {name}: exit {rc}, {time.time() - t:.0f} s")
        if rc == 0:
            bins[name] = f"{b}/cppwrf"

    env = dict(os.environ, CPPWRF_WRF_RUNS=WR, CPPWRF_WRF_REF="/content/data/wrf_ref", OMP_NUM_THREADS="4")
    if cases and bins:
        args = " ".join(f"--{k.replace('_', '-')} {v}" for k, v in bins.items())
        sh(f"python3 tests/verify/wrf_suite.py {args} --only {cases} --out /content/suite", "suite.log", cwd=SRC,
           env=env)
        for line in open(os.path.join(JOB, "suite.log")):
            if re.match(r"(PASS|FAIL|SKIP|ERROR|SUITE)\s", line):
                note("suite " + line.rstrip())

    if files and "cuda_dp" in bins:  # the GPU's wrfout / wrfrst files: variables, attributes and times as WRF's
        sh(f"bash tests/verify/run_wrf_files_test.sh /content/build-cuda_dp {' '.join(files)}", "files.log", cwd=SRC,
           env=dict(env, CPPWRF_FILES_NO_VALUES="1", CPPWRF_SCRATCH="/content/files"))
        for line in open(os.path.join(JOB, "files.log")):
            if re.match(r"(PASS|FAIL|SKIP|compared)\s", line):
                note("files " + line.rstrip())

    if bench:  # the suite cases again, without history output, for the stepping time
        sys.path.insert(0, os.path.join(SRC, "tests", "verify"))
        os.environ.update(CPPWRF_WRF_RUNS=WR)
        import wrf_suite
        keys = [k for k in os.environ.get("CPPWRF_JOB_ENV_KEYS", "").split(",") if k]
        plain = {k: v for k, v in env.items() if k not in keys}
        prof = {k: v for k, v in env.items() if k != "CPPWRF_KTIME"}
        nsys = os.environ.get("CPPWRF_JOB_NSYS") == "1"
        ncu = os.environ.get("CPPWRF_JOB_NCU", "")
        if nsys or ncu:
            res = {k: res_usage(b) for k, b in bins.items()}
        cases = {c[0]: c for c in wrf_suite.CASES}
        for spec in bench:
            name, _, over = spec.partition("@")
            if name not in cases:
                note(f"bench {spec}: no such suite case")
                continue
            _, key, cfg, sets, wo, t, tol = cases[name]
            key = over or key
            if key not in bins:
                continue
            tag = name + (f"@{over}" if over else "")

            def cmd(o, run_time):
                c = [bins[key], cfg, "--set", f"output.dir={o}"]
                for s in sets + [f"run.run_time={int(run_time)}", "output.out_interval=1000000",
                                 "output.diag_interval=1000000"]:
                    c += ["--set", s]
                return c

            o = f"/content/bench/{tag}"
            with GpuMemory() as mem:
                sh(f"rm -rf {o} && {shlex.join(cmd(o, t))}", f"bench-{tag}.log", cwd=SRC, env=plain)
            # one line per domain (a nest's parent first: its stepping time contains the nests' steps)
            done = timings(os.path.join(JOB, f"bench-{tag}.log"))
            note(f"bench {tag}: " + ("; ".join(f"{s} s stepping, {n} ns/cell/step, {k} steps" for k, s, n in done)
                                     if done else "no timing (see log)") + f"; device memory {mem.text()}")
            if keys:
                ev = " ".join(f"{k}={os.environ[k]}" for k in keys)
                sh(f"rm -rf {o} && {shlex.join(cmd(o, t))}", f"bench-env-{tag}.log", cwd=SRC, env=env)
                d2 = timings(os.path.join(JOB, f"bench-env-{tag}.log"))
                note(f"bench {tag} ({ev}): " + (f"{d2[0][1]} s stepping" if d2 else "no timing (see log)"))
            if not done:
                continue
            top = profile_nsys(tag, cmd, t, prof, note) if nsys else []
            if ncu:
                profile_ncu(tag, cmd, 3 * t / done[0][0], prof, ncu, top, res[key], note)


if __name__ == "__main__":
    main()
