#!/usr/bin/env python3
"""WRF's ideal cases as WRF output files: each CPU FP64 ideal case of wrf_suite.py runs with
output.format = wrfout and the namelist of its WRF reference run (history interval, run length, frames),
and every variable, dimension and attribute of the wrfout file must equal WRF's (tests/init_real/compare_real.py).

  ideal_files_test.py BUILD_DIR [CASE ...]

CPPWRF_IDEAL_TRANSLATED=1: the configuration is instead tools/wrf2cppwrf.py's translation of the reference
run's namelist.input (sounding: the input_sounding next to it), so the translator is checked as well.

References: ${CPPWRF_WRF_REF} (as wrf_suite.py). Outputs in ${CPPWRF_SCRATCH:-~/cppwrf_scratch}/ideal_files_test,
deleted for passing cases. Exit status 0 = all PASS.
"""
import os
import re
import shutil
import subprocess
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "tests", "verify"))
import wrf_suite  # noqa: E402


def config_value(cfg, section, key):
    sec = None
    for line in open(os.path.join(REPO, cfg)):
        line = line.split("#")[0].strip()
        m = re.match(r"\[(.+)\]", line)
        if m:
            sec = m.group(1)
        elif sec == section and "=" in line and line.split("=", 1)[0].strip() == key:
            return line.split("=", 1)[1].strip()
    return None


def main():
    exe = os.path.join(sys.argv[1], "cppwrf")
    want = sys.argv[2:]
    work = os.path.join(os.environ.get("CPPWRF_SCRATCH", os.path.expanduser("~/cppwrf_scratch")), "ideal_files_test")
    ok, n = True, 0
    for name, key, cfg, sets, wo, t, tol in wrf_suite.CASES:
        if not name.startswith("ideal_") or key != "cpu_dp" or (want and name not in want):
            continue
        ref = os.path.dirname(wo)
        if not os.path.exists(wo):
            print(f"SKIP {name}: {wo} not found")
            continue
        nl = open(os.path.join(ref, "namelist.input")).read()
        hs, hm = re.search(r"history_interval_s\s*=\s*(\d+)", nl), re.search(r"history_interval\s*=\s*(\d+)", nl)
        hist = int(hs.group(1)) if hs else int(hm.group(1)) * 60
        run = sum(int(m.group(1)) * f for f, k in ((3600, "hours"), (60, "minutes"), (1, "seconds"))
                  for m in [re.search(rf"run_{k}\s*=\s*(\d+)", nl)] if m)
        sets = [s for s in sets if not s.startswith(("run.run_time=", "output."))]
        dt = float(dict(s.split("=", 1) for s in sets).get("time.dt", config_value(cfg, "time", "dt")))
        out = os.path.join(work, name)
        shutil.rmtree(out, ignore_errors=True)
        if os.environ.get("CPPWRF_IDEAL_TRANSLATED") == "1":
            os.makedirs(out)
            ini = os.path.join(out, "translated.ini")
            t = subprocess.run([sys.executable, os.path.join(REPO, "tools", "wrf2cppwrf.py"),
                                os.path.join(ref, "namelist.input"), "-o", ini, "--format", "wrfout"],
                               capture_output=True, text=True)
            if t.returncode != 0:
                print(f"FAIL {name}: translation failed: {(t.stdout + t.stderr)[-300:]}")
                ok = False
                continue
            cmd = [exe, ini]
            sets = [f"output.wrf_namelist={ref}/namelist.input", f"output.dir={out}"]
        else:
            cmd = [exe, cfg]
            sets = sets + [f"run.run_time={run}", "output.format=wrfout", f"output.out_interval={round(hist / dt)}",
                           "output.diag_interval=100000", f"output.wrf_namelist={ref}/namelist.input",
                           f"output.dir={out}"]
        for s in sets:
            cmd += ["--set", s]
        r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"FAIL {name}: cppwrf failed: {(r.stdout + r.stderr)[-300:]}")
            ok = False
            continue
        c = subprocess.run([sys.executable, os.path.join(REPO, "tests/init_real/compare_real.py"),
                            os.path.join(out, os.path.basename(wo)), wo], capture_output=True, text=True)
        n += 1
        if c.returncode == 0:
            print(f"PASS {name}")
            shutil.rmtree(out, ignore_errors=True)
        else:
            print(f"FAIL {name}")
            print("\n".join((c.stdout + c.stderr).strip().splitlines()[:20]))
            ok = False
    print(f"compared {n} files")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
