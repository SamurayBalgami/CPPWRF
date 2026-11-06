#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Markdown tables of a tools/bench/colab_bench.sh run.

  bench_report.py RESULTS.txt [--rate CU_PER_HOUR] [--logs LOGS]

RESULTS.txt: the JSON lines of tools/bench/bench_run.py (runs, fastest WRF layouts, comparisons). Per case: WRF's
process layouts (history off), then every model at WRF's fastest layout with history off and on: wall time of the
process, stepping time (= seconds per simulated hour: every run is 1 h), time per step and ns per cell and step
(stepping time over the sum of cells x steps of all domains), the speed-up of the stepping and of the wall time
against WRF real*4 and real*8 at their fastest layouts, the time outside the stepping (wall - stepping: process
start, reading the inputs, initialisation, output), Colab compute units per simulated hour (wall time x the
runtime's rate) and the GPU's device memory (nvidia-smi's peak used memory minus the value before the run, of which
CPPWRF_MEMINFO's stack limit x resident threads is the local-memory reserve; the estimate for an RTX 2060, 6 GB with
30 SMs x 1024 resident threads, replaces that reserve by its own: an upper bound, as the RRTMG batch scratch that
CPPWRF sizes to the GPU's memory is smaller there). Then the comparisons of the outputs. The logs
(LOGS, default logs/ next to RESULTS.txt) give the CPPWRF_MEMINFO lines.
"""
import json
import os
import re
import sys

SM_2060 = 30 * 1024  # RTX 2060: 30 SMs x 1024 resident threads (the device-dependent local-memory reserve)


def stack_reserve(logs, r):
    """(stack limit bytes, resident threads) of CPPWRF_MEMINFO's line "at the end" in the run's log, or None."""
    f = f"{logs}/{r['case']}/cpp_{r['hist']}.{r['model']}_{r['hist']}.rep{r.get('rep', 1)}.log"
    if not logs or not os.path.exists(f):
        return None
    m = re.search(r"device memory at the end: \d+ MiB used of \d+ MiB; stack limit (\d+) bytes x (\d+) resident", open(f).read())
    return (int(m.group(1)), int(m.group(2))) if m else None

ORDER = ["wrf_r4", "wrf_r8", "gpu_fp32", "gpu_fp64", "cpu_fp64", "cpu_fp64_hwlog"]
NAMES = {"wrf_r4": "WRF real*4", "wrf_r8": "WRF real*8", "gpu_fp32": "CPPWRF GPU FP32", "gpu_fp64": "CPPWRF GPU FP64",
         "cpu_fp64": "CPPWRF CPU FP64", "cpu_fp64_hwlog": "CPPWRF CPU FP64 (native libmvec)"}
CASES = {"R1": "R1: dynamics only, 119 x 89 x 44, dt 60 s",
         "P4": "P4: full physics (KF, WSM6, YSU, Noah, RRTMG radt 1 min), 119 x 89 x 44, dt 60 s",
         "P4R15": "P4 with radt 15 min (production setting)",
         "N3": "N3: two-way nest, P4 physics on both domains",
         "L": "L: large domain, P4 physics with radt 15 min, 359 x 269 x 44 at 15 km (9.1 x P4), dt 60 s"}


def fmt(x, n=2):
    return "-" if x is None else f"{x:.{n}f}"


def report(path, rate=None, logs=None):
    if logs is None:
        logs = os.path.join(os.path.dirname(os.path.abspath(path)), "logs")
    recs = [json.loads(line) for line in open(path) if line.strip()]
    runs = [r for r in recs if r.get("kind") == "run"]
    best = {(r["case"], r["model"]): r["layout"] for r in recs if r.get("kind") == "best"}
    cases = []
    for r in recs:
        if r.get("case") and r["case"] not in cases and r.get("kind") in ("run", "best"):
            cases.append(r["case"])
    out = []
    w = out.append
    for c in cases:
        w(f"### {CASES.get(c, c)}\n")
        g = next((r["grid"] for r in runs if r["case"] == c and r.get("grid")), None)
        if g:
            w("Cells and steps per simulated hour: " + ", ".join(f"d{i + 1:02d} {a} x {b}" for i, (a, b) in enumerate(g))
              + "\n")
        scan = {}
        for r in runs:
            if r["case"] == c and r["model"].startswith("wrf") and r["hist"] == "off":
                scan.setdefault(r["layout"], {})[r["model"]] = r
        if scan:
            w("WRF process layouts (history off), stepping time in s (wall time in s):\n")
            w("| layout | real*4 | real*8 |")
            w("|---|---|---|")
            for lay, m in scan.items():
                cell = lambda r: (f"{r['stepping']:.2f} ({r['wall']:.1f})" if r.get("ok") else  # noqa: E731
                                  "failed") if r else "-"
                w(f"| {lay} | {cell(m.get('wrf_r4'))} | {cell(m.get('wrf_r8'))} |")
            w("")
        for r in runs:  # the serial run of bench_run.py --serial (a comparison, not a layout of the scan)
            if r["case"] == c and r["model"] == "wrf_r8" and r["layout"] == "dm1" and r.get("ok"):
                w(f"WRF real*8 in one process (dm1, history on): stepping {r['stepping']:.2f} s, wall {r['wall']:.1f} s.\n")
        sel = {}
        for r in runs:
            if r["case"] != c or not r.get("ok"):
                continue
            if r["model"].startswith("wrf") and r["layout"] != best.get((c, r["model"])):
                continue
            sel[(r["model"], r["hist"])] = r
        ref = {p: sel.get((f"wrf_{p}", "off")) for p in ("r4", "r8")}
        w("| model | layout | history | wall s | stepping s (per sim. hour) | s / step | ns / cell / step "
          "| speed-up vs WRF r4 (stepping / wall) | vs WRF r8 | wall - stepping s | CU / sim. hour | GPU memory MiB |")
        w("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for m in ORDER:
            for h in ("off", "on"):
                r = sel.get((m, h))
                if not r:
                    continue
                g = r.get("grid") or [[0, 1]]
                cs = sum(a * b for a, b in g)
                st = r["stepping"]
                rf = {p: sel.get((f"wrf_{p}", h)) or ref[p] for p in ("r4", "r8")}
                sp = {p: (f"{rf[p]['stepping'] / st:.1f} / {rf[p]['wall'] / r['wall']:.1f}" if rf[p] and st > 0 else "-")
                      for p in rf}
                mem = "-"
                if r.get("mem_peak") is not None and r.get("mem_idle") is not None:
                    used = r["mem_peak"] - r["mem_idle"]
                    sr = stack_reserve(logs, r)
                    if sr:  # the local-memory reserve is stack x resident threads: a smaller GPU reserves less
                        res = sr[0] * sr[1] / 2 ** 20
                        est = used - res + sr[0] * SM_2060 / 2 ** 20
                        mem = (f"{used} ({res:.0f} stack reserve); RTX 2060 at most ~{est:.0f}"
                               + (" (fits 6 GB)" if est < 5800 else ""))
                    else:
                        mem = f"{used}"
                cu = fmt(r["wall"] / 3600 * rate, 3) if rate else "-"
                w(f"| {NAMES.get(m, m)} | {r['layout']} | {h} | {r['wall']:.1f} | {st:.2f} | {st / g[0][1]:.3f} | "
                  f"{st / cs * 1e9:.1f} | {sp['r4']} | {sp['r8']} | {r['wall'] - st:.1f} | {cu} | {mem} |")
        fails = [r for r in runs if r["case"] == c and not r.get("ok")]
        if fails:
            w("\nFailed runs: " + "; ".join(f"{r['model']} {r.get('layout', '')} {r.get('hist', '')}: "
                                           f"{str(r.get('error'))[:200]}" for r in fails))
        w("\nOutputs compared:\n")
        for r in recs:
            if r.get("case") != c or r.get("kind") not in ("exact", "l2"):
                continue
            parts = []
            for d, v in r["domains"].items():
                if "error" in v:
                    parts.append(f"{d}: {v['error']}")
                elif r["kind"] == "exact":
                    parts.append(f"{d}: {v['frames']} frames, " + ("identical (tolerance 0)" if not v["differ"] else
                                 "differ: " + "; ".join(f"{k} {str(x)[:90]}" for k, x in list(v["differ"].items())[:8])
                                 + (f" ... ({len(v['differ'])} variables)" if len(v["differ"]) > 8 else "")))
                else:
                    parts.append(f"{d} at {v['time']}: " + ", ".join(f"{k} {x:.1e}" for k, x in v["l2"].items()
                                                                      if x is not None))
            w(f"- {r['what']}" + (" (relative L2)" if r["kind"] == "l2" else "") + ": " + " | ".join(parts))
        w("")
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    a = sys.argv[1:]
    rate = float(a[a.index("--rate") + 1]) if "--rate" in a else None
    sys.stdout.write(report(a[0], rate, a[a.index("--logs") + 1] if "--logs" in a else None))
