#!/usr/bin/env python3
"""Compare driver_noah_ice outputs (WRF Fortran vs C++): every output of both routines."""
import argparse
import sys
import numpy as np

GLACIAL = ["t1", "snowh", "sneqv", "sncovr", "snotime1", "ribb", "emissi", "albedo", "z0", "eta", "sheat",
           "eta_kinematic", "fdown", "esnow", "dew", "etp", "ssoil", "flx1", "flx2", "flx3", "snomlt", "runoff1",
           "q1", "ch", "alb", "z0brd", "embrd"]
SEAICE = ["t1", "snowh", "sneqv", "sncovr", "ribb", "emissi", "albedo", "z0", "z0brd", "eta", "sheat",
          "eta_kinematic", "fdown", "esnow", "dew", "etp", "ssoil", "flx1", "flx2", "flx3", "snomlt", "runoff1",
          "q1", "ch"]


def read(path):
    with open(path, "rb") as f:
        ncols, nsoil = np.frombuffer(f.read(8), dtype="=i4")
        data = np.frombuffer(f.read(), dtype="=f8")
    return int(ncols), int(nsoil), data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fortran")
    ap.add_argument("cpp")
    ap.add_argument("--tol", type=float, default=0.0)
    a = ap.parse_args()
    nf, sf, df = read(a.fortran)
    nc, sc, dc = read(a.cpp)
    if (nf, sf, df.size) != (nc, sc, dc.size):
        print(f"ERROR: shape mismatch: Fortran {nf}x{sf} ({df.size}), C++ {nc}x{sc} ({dc.size})")
        return 1
    names = [f"glacial.{n}" for n in GLACIAL + [f"stc{k + 1}" for k in range(sf)]] + \
            [f"seaice.{n}" for n in SEAICE + [f"stc{k + 1}" for k in range(sf)]]
    if len(names) * nf != df.size:
        print(f"ERROR: {df.size} values, expected {len(names)} x {nf}")
        return 1
    ok = True
    print(f"Noah ice outputs, {nf} columns, tol = {a.tol:e}")
    for n, name in enumerate(names):
        f = df[n * nf:(n + 1) * nf]
        c = dc[n * nf:(n + 1) * nf]
        if np.isnan(f).any() or np.isnan(c).any():
            bad = not np.array_equal(np.isnan(f), np.isnan(c))
        else:
            bad = False
        d = np.abs(np.nan_to_num(f) - np.nan_to_num(c))
        mx = float(d.max())
        bad = bad or mx > a.tol
        ok = ok and not bad
        print(f"{name:<24} max abs diff {mx:.3e}  range [{np.nanmin(f):.6g}, {np.nanmax(f):.6g}]  "
              f"{'FAIL (col ' + str(int(d.argmax()) + 1) + ')' if bad else 'OK'}")
    print("PASSED" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
