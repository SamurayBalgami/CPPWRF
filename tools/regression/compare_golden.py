#!/usr/bin/env python3
"""Compares a CPPWRF run output against a stored golden output.
Checks diag.csv (excluding timing columns) and the final snapshot.
Tolerances:
  FP64: rtol=1e-12, atol=1e-12
  FP32: rtol=1e-5,  atol=1e-5
  diag.csv values stored with fewer significant digits: at least half a unit in their last digit
"""
import csv
import glob
import json
import os
import sys
import numpy as np

def stored_rtol(text):
    """Relative rounding of a decimal written with the significant digits of `text` (0.5 ulp of the last digit)."""
    mant = text.strip().lower().split("e")[0].lstrip("+-").replace(".", "").lstrip("0")
    return 0.5 * 10.0 ** (1 - len(mant)) if mant else 0.0

def read_diag(d):
    path = os.path.join(d, "diag.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} not found")
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"{path} is empty")
    return rows

def find_last_step(d):
    files = glob.glob(os.path.join(d, "snap_*.json"))
    steps = []
    for f in files:
        base = os.path.basename(f)
        s = base.replace("snap_", "").replace(".json", "")
        if s.isdigit():
            steps.append(int(s))
    return max(steps) if steps else None

def read_snap(d, step):
    base = os.path.join(d, f"snap_{step:06d}")
    with open(base + ".json") as f:
        meta = json.load(f)
    out = {}
    with open(base + ".bin", "rb") as f:
        for v in meta["variables"]:
            dt = np.dtype({"float32": "<f4", "float64": "<f8"}[v["dtype"]])
            shape = tuple(v["shape"])
            f.seek(int(v["offset"]))
            count = int(np.prod(shape))
            a = np.frombuffer(f.read(count * dt.itemsize), dtype=dt)
            out[v["name"]] = a.reshape(shape)
    return out

def compare(run_dir, golden_dir, prec="float32"):
    if not os.path.exists(golden_dir):
        print(f"FAIL: Golden dir {golden_dir} does not exist.")
        return 1

    is_dp = prec in ("float64", "fp64") or "dp" in prec  # run.sh passes fp32 / fp64
    rtol = 1e-12 if is_dp else 1e-5
    atol = 1e-12 if is_dp else 1e-5

    # 1. Compare diag.csv
    try:
        run_diag = read_diag(run_dir)
        gold_diag = read_diag(golden_dir)
    except Exception as e:
        print(f"FAIL reading diag.csv: {e}")
        return 1

    if len(run_diag) != len(gold_diag):
        print(f"FAIL: diag.csv row count mismatch ({len(run_diag)} vs {len(gold_diag)})")
        return 1

    # Skip nondeterministic timing column 'wall_s'
    cols_to_check = [c for c in run_diag[0].keys() if c != "wall_s"]
    for row_idx, (r_row, g_row) in enumerate(zip(run_diag, gold_diag)):
        for c in cols_to_check:
            val_r = float(r_row[c])
            val_g = float(g_row[c])
            # a golden written with fewer digits (older diag.csv: 9 significant digits) is checked to the
            # precision it was stored with, not to the FP64 tolerance it cannot meet
            if not np.isclose(val_r, val_g, rtol=max(rtol, stored_rtol(g_row[c])), atol=atol):
                diff = abs(val_r - val_g)
                print(f"FAIL: diag.csv mismatch at row {row_idx}, col '{c}': run={val_r} gold={val_g} diff={diff:.3e}")
                return 1

    # 2. Compare final snapshot
    step_r = find_last_step(run_dir)
    step_g = find_last_step(golden_dir)
    if step_r is None or step_g is None:
        print(f"FAIL: snapshot not found (run step: {step_r}, gold step: {step_g})")
        return 1
    if step_r != step_g:
        print(f"FAIL: final snapshot step mismatch (run: {step_r}, gold: {step_g})")
        return 1

    snap_r = read_snap(run_dir, step_r)
    snap_g = read_snap(golden_dir, step_g)

    for vname in sorted(gold_snap := set(snap_g.keys())):
        if vname not in snap_r:
            print(f"FAIL: variable '{vname}' missing in run snapshot")
            return 1
        ar = snap_r[vname]
        ag = snap_g[vname]
        if ar.shape != ag.shape:
            print(f"FAIL: '{vname}' shape mismatch: run {ar.shape} vs gold {ag.shape}")
            return 1
        if not np.allclose(ar, ag, rtol=rtol, atol=atol):
            max_diff = np.max(np.abs(ar - ag))
            print(f"FAIL: '{vname}' value mismatch: max diff = {max_diff:.3e}")
            return 1

    print("PASS: Run output matches golden reference within tolerance.")
    return 0

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: compare_golden.py <run_dir> <golden_dir> [precision]")
        sys.exit(2)
    p = sys.argv[3] if len(sys.argv) > 3 else "float32"
    rc = compare(sys.argv[1], sys.argv[2], p)
    if rc:
        # a golden older than the code it checks fails for every change since; name its commit
        commit = os.path.join(os.path.dirname(os.path.abspath(sys.argv[2])), "COMMIT")
        if os.path.exists(commit):
            with open(commit) as f:
                print(f"       golden made at commit {f.read().strip()[:12]} (tools/regression/update_golden.sh "
                      "refreshes it once the change is verified)")
    sys.exit(rc)
