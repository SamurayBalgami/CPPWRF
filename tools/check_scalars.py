#!/usr/bin/env python3
"""Scalar checks on a CPPWRF output dir: tools/check_scalars.py <out_dir> [name ...]
For each scalar in the snapshots prints min, max, max|q-1| and the column-mass integral
sum(mu_full * q * |dnw|) relative to the initial snapshot (dry-mass weighted tracer mass)."""
import json, sys, glob, os
import numpy as np

def load(d, step_json):
    meta = json.load(open(step_json))
    raw = open(step_json[:-5] + ".bin", "rb").read()
    out = {}
    for v in meta["variables"]:
        dt = np.float32 if v["dtype"] == "float32" else np.float64
        n = int(np.prod(v["shape"]))
        out[v["name"]] = np.frombuffer(raw, dt, n, v["offset"]).reshape(v["shape"]).astype(np.float64)
    return out

def main():
    d = sys.argv[1]
    run = json.load(open(os.path.join(d, "run.json")))
    eta = np.array(run["grid"]["eta_levels"])
    dnw = -np.diff(eta)
    names = sys.argv[2:] or run["params"].get("scalars", [])
    snaps = sorted(glob.glob(os.path.join(d, "snap_*.json")))
    first = load(d, snaps[0])
    m0 = {n: np.sum(first["mu"] * first[n] * dnw[:, None, None]) for n in names}
    for s in snaps:
        f = load(d, s)
        row = [os.path.basename(s)[5:11]]
        for n in names:
            q = f[n]
            m = np.sum(f["mu"] * q * dnw[:, None, None])
            rel = (m - m0[n]) / m0[n] if m0[n] != 0 else m - m0[n]
            row.append(f"{n}: min {q.min():+.3e} max {q.max():.6e} max|q-1| {np.abs(q-1).max():.2e} dmass {rel:+.2e}")
        print("  ".join(row))

main()
