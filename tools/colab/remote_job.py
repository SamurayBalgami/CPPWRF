#!/usr/bin/env python3
"""The Colab VM side of tools/colab/colab_run.sh (runs detached under nohup; /content/job/done marks the end).

Expects /content/src.tar.gz (git archive of the repository) and /content/data.tar (tools/colab/stage_data.py),
and the environment CPPWRF_JOB_CASES (wrf_suite.py case names, comma-separated), CPPWRF_JOB_FILES (cases of
tests/verify/run_wrf_files_test.sh, may be empty) and CPPWRF_JOB_BENCH (suite case names to time again without
output, may be empty). Builds the CUDA FP32 and FP64 binaries for the VM's GPU, runs the suite cases against the
WRF references, the files test with the FP64 binary (CPPWRF_FILES_NO_VALUES=1), times the bench cases, and leaves
everything in /content/job (summary.txt, logs) and /content/results.tar.gz.
"""
import os
import re
import subprocess
import sys
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
        for name, key, cfg, sets, wo, t, tol in wrf_suite.CASES:
            if name not in bench or key not in bins:
                continue
            o = f"/content/bench/{name}"
            cmd = [bins[key], cfg, "--set", f"output.dir={o}"]
            for s in sets + [f"run.run_time={int(t)}", "output.out_interval=1000000", "output.diag_interval=1000000"]:
                cmd += ["--set", s]
            sh("rm -rf " + o + " && " + " ".join(f"'{c}'" for c in cmd), f"bench-{name}.log", cwd=SRC, env=env)
            m = re.findall(r"stepping ([0-9.]+) s, ([0-9.]+) ns/cell/step", open(os.path.join(JOB, f"bench-{name}.log")).read())
            note(f"bench {name}: " + (f"{m[-1][0]} s stepping, {m[-1][1]} ns/cell/step" if m else "no timing (see log)"))


if __name__ == "__main__":
    main()
