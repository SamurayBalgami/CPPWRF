#!/usr/bin/env python3
"""Compare the C++ P1 suite (suite_cpp.bin) with WRF v4.6.1 (suite_ref.bin).

Every record the C++ driver writes is compared with the WRF record of the same name
(P1_n = after first_rk_step_part1, P2_n = after first_rk_step_part2, OUTn = after
phy_prep_part2, n = 1..3). Pass: every element bit-identical or |c - f| <= tol*max(|c|,|f|)
(tol = 1e-12, FP64). Prints one table row per field (worst over the 3 calls) and, with
--markdown, the same table in Markdown.
"""
import argparse
import sys

import os

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from suite_io import read_records  # noqa: E402

GROUPS = [
    ("phy_prep", ["TH_PHY", "P_PHY", "PI_PHY", "T_PHY", "RHO", "U_PHY", "V_PHY", "Z", "Z_AT_W", "DZ8W", "P8W", "T8W",
                  "P_HYD", "P_HYD_W"]),
    ("sfclayrev", ["UST", "ZNT", "MOL", "RMOL", "PSIM", "PSIH", "GZ1OZ0", "FM", "FH", "ZOL", "BR", "WSPD", "REGIME",
                   "FLHC", "FLQC", "U10", "V10", "USTM", "CK", "CKA", "CD", "CDA", "PSFC"]),
    ("Noah / ice / sfcdiags", ["HFX", "QFX", "LH", "TSK", "QSFC", "SFCEXC", "T2", "TH2", "Q2", "GRDFLX", "TSLB",
                               "SMOIS", "SH2O", "SMCREL", "SNOW", "SNOWH", "SNOWC", "CANWAT", "ALBEDO", "ALBBCK",
                               "EMISS", "Z0", "LAI", "SNOTIME", "SMSTAV", "SMSTOT", "SFCRUNOFF", "UDRUNOFF", "ACSNOM",
                               "ACSNOW", "SNOPCX", "POTEVP", "NOAHRES", "SFCEVP", "ACHFX", "ACLHF", "ACGRDFLX",
                               "CHKLOWQ", "RAINBL"]),
    ("YSU", ["PBLH", "KPBL", "EXCH_H", "EXCH_M", "WSTAR_YSU", "DELTA_YSU"]),
    ("R*BLTEN (decoupled, OUT)", ["RUBLTEN", "RVBLTEN", "RTHBLTEN", "RQVBLTEN"]),
    ("part2 (coupled, *_tendf)", ["P2:RUBLTEN", "P2:RVBLTEN", "P2:RTHBLTEN", "P2:RQVBLTEN", "P2:RU_TENDF",
                                  "P2:RV_TENDF", "P2:T_TENDF", "P2:QV_TENDF"]),
]


def stats(c, f, tol):
    d = np.abs(c - f)
    scale = np.maximum(np.abs(c), np.abs(f))
    exact = int(np.sum(c == f))
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(d == 0, 0.0, d / np.where(scale > 0, scale, 1.0))
    bad = int(np.sum(d > tol * scale))
    return exact, c.size, float(d.max()) if d.size else 0.0, float(rel.max()) if rel.size else 0.0, bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("cpp")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--markdown", action="store_true")
    a = ap.parse_args()
    F = read_records(a.ref)
    C = read_records(a.cpp)
    rows = []
    ok = True
    seen = set()
    for group, names in GROUPS:
        for name in names:
            if name.startswith("P2:"):
                keys = [f"P2_{n}:{name[3:]}" for n in (1, 2, 3)]
            elif group == "R*BLTEN (decoupled, OUT)":
                keys = [f"OUT{n}:{name}" for n in (1, 2, 3)]
            else:
                keys = [f"P1_{n}:{name}" for n in (1, 2, 3)] + [f"OUT{n}:{name}" for n in (1, 2, 3)]
            keys = [k for k in keys if k in C and k in F]
            if not keys:
                print(f"missing: {name}", file=sys.stderr)
                ok = False
                continue
            tot = [0, 0, 0.0, 0.0, 0]
            for k in keys:
                seen.add(k)
                e, n, dmax, rmax, bad = stats(C[k], F[k], a.tol)
                tot[0] += e
                tot[1] += n
                tot[2] = max(tot[2], dmax)
                tot[3] = max(tot[3], rmax)
                tot[4] += bad
            passed = tot[4] == 0
            ok &= passed
            rows.append((group, name.replace("P2:", ""), tot[0], tot[1], tot[2], tot[3], "PASS" if passed else "FAIL"))
    # anything else written by both sides is compared as well
    extra_bad = []
    for k in sorted(set(C) & set(F) - seen):
        e, n, dmax, rmax, bad = stats(C[k], F[k], a.tol)
        if bad:
            extra_bad.append((k, rmax))
    for k, r in extra_bad:
        print(f"FAIL (extra record) {k}: max rel {r:.3e}")
        ok = False

    w = max(len(r[1]) for r in rows)
    print(f"{'group':26s} {'field':{w}s} {'bit-exact':>17s} {'max |diff|':>11s} {'max rel':>10s}")
    for g, n, e, tot, dmax, rmax, st in rows:
        print(f"{g:26s} {n:{w}s} {e:8d}/{tot:<8d} {dmax:11.3e} {rmax:10.3e} {st}")
    if a.markdown:
        print("\n| group | field | bit-identical elements | max abs diff | max rel diff | result |")
        print("|---|---|---|---|---|---|")
        for g, n, e, tot, dmax, rmax, st in rows:
            print(f"| {g} | {n} | {e}/{tot} | {dmax:.2e} | {rmax:.2e} | {st} |")
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
