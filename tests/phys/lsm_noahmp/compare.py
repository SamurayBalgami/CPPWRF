#!/usr/bin/env python3
"""
Differential validator for Noah-MP v5.1.2 (sf_surface_physics = 4).
Compares WRF v4.8.0 Fortran outputs (linked from libwrflib.a) vs C++ outputs
across 62 2D INOUT state variables, 77 2D OUT diagnostic variables, and 11 3D
snow/soil/albedo state arrays (150 variables total).
Target criterion: max_abs_diff == 0.0 in CPU FP64.
"""

import argparse
import struct
import sys
import numpy as np

INOUT_2D = [
    "tsk", "hfx", "qfx", "lh", "grdflx", "smstav", "smstot", "sfcrunoff", "udrunoff",
    "albedo", "snowc", "snow", "snowh", "canwat", "acsnom", "acsnow", "emiss", "qsfc",
    "z0", "znt", "rs", "isnowxy", "tvxy", "tgxy", "canicexy", "canliqxy", "eahxy",
    "tahxy", "cmxy", "chxy", "fwetxy", "sneqvoxy", "alboldxy", "qsnowxy", "qrainxy",
    "wslakexy", "zwtxy", "waxy", "wtxy", "lfmassxy", "rtmassxy", "stmassxy", "woodxy",
    "stblcpxy", "fastcpxy", "lai", "xsaixy", "taussxy", "smcwtdxy", "deeprechxy",
    "rechxy", "grainxy", "gddxy", "qtdrain", "acc_ssoilxy", "acc_qinsurxy",
    "acc_qsevaxy", "acc_dwaterxy", "acc_prcpxy", "acc_ecanxy", "acc_etranxy", "acc_edirxy",
]

OUT_2D = [
    "t2mvxy", "t2mbxy", "q2mvxy", "q2mbxy", "tradxy", "neexy", "gppxy", "nppxy",
    "fvegxy", "runsfxy", "runsbxy", "ecanxy", "edirxy", "etranxy", "fsaxy", "firaxy",
    "aparxy", "psnxy", "savxy", "sagxy", "rssunxy", "rsshaxy", "bgapxy", "wgapxy",
    "tgvxy", "tgbxy", "chvxy", "chbxy", "shgxy", "shcxy", "shbxy", "evgxy", "evbxy",
    "ghvxy", "ghbxy", "irgxy", "ircxy", "irbxy", "trxy", "evcxy", "chleafxy", "chucxy",
    "chv2xy", "chb2xy", "pahxy", "pahgxy", "pahvxy", "pahbxy", "qintsxy", "qintrxy",
    "qdripsxy", "qdriprxy", "qthrosxy", "qthrorxy", "qsnsubxy", "qsnfroxy",
    "qsubcxy", "qfrocxy", "qevacxy", "qdewcxy", "qfrzcxy", "qmeltcxy", "qsnbotxy",
    "qmeltxy", "pondingxy", "fpicexy", "rainlsm", "snowlsm", "forctlsm", "forcqlsm",
    "forcplsm", "forczlsm", "forcwlsm", "eflxbxy", "soilenergy", "snowenergy", "canhsxy",
]


def compare_noahmp_outputs(fortran_file, cpp_file, tol=0.0):
    with open(fortran_file, "rb") as ff, open(cpp_file, "rb") as cf:
        f_hdr = ff.read(8)
        c_hdr = cf.read(8)
        n_cols_f, nsoil_f = struct.unpack("=ii", f_hdr)
        n_cols_c, nsoil_c = struct.unpack("=ii", c_hdr)
        if (n_cols_f, nsoil_f) != (n_cols_c, nsoil_c):
            print(f"ERROR: Header mismatch: Fortran=({n_cols_f},{nsoil_f}), C++=({n_cols_c},{nsoil_c})")
            return False
        f_data = ff.read()
        c_data = cf.read()

    if len(f_data) != len(c_data):
        print(f"ERROR: Payload size mismatch: Fortran={len(f_data)} bytes, C++={len(c_data)} bytes")
        return False

    n_cols = n_cols_f
    nsoil = nsoil_f
    nsnow = 3

    vars_2d = INOUT_2D + OUT_2D
    vars_3d = [
        ("smois", nsoil),
        ("sh2o", nsoil),
        ("tslb", nsoil),
        ("smoiseq", nsoil),
        ("acc_etranixy", nsoil),
        ("tsnoxy", nsnow),
        ("snicexy", nsnow),
        ("snliqxy", nsnow),
        ("zsnsoxy", nsnow + nsoil),
        ("albsoildirxy", 2),
        ("albsoildifxy", 2),
    ]

    offset = 0
    all_pass = True
    global_max_abs = 0.0
    n_checked = 0

    print(f"Comparing Noah-MP v5.1.2 outputs ({n_cols} columns, nsoil={nsoil}, tol={tol:.1e}):")
    print(f"{'Variable':<16} | {'Shape':<10} | {'Active Pts':<10} | {'Max Abs Diff':<14} | {'Max Rel Diff':<14} | {'Status':<6}")
    print("-" * 82)

    for name in vars_2d:
        f_arr = np.frombuffer(f_data, dtype=np.float64, count=n_cols, offset=offset)
        c_arr = np.frombuffer(c_data, dtype=np.float64, count=n_cols, offset=offset)
        offset += n_cols * 8

        abs_diff = np.abs(f_arr - c_arr)
        denom = np.maximum(np.maximum(np.abs(f_arr), np.abs(c_arr)), 1e-12)
        rel_diff = abs_diff / denom

        max_abs = float(np.max(abs_diff))
        max_rel = float(np.max(rel_diff))
        global_max_abs = max(global_max_abs, max_abs)
        active_pts = int(np.count_nonzero((np.abs(f_arr) > 1e-15) & (f_arr != -9999.0)))
        n_checked += 1

        passed = (max_abs <= tol) or (max_rel <= tol)
        status = "PASS" if passed else "FAIL"
        if not passed:
            all_pass = False
            worst = int(np.argmax(abs_diff))
            print(f"{name:<16} | {f'({n_cols})':<10} | {active_pts:<10} | {max_abs:<14.6e} | {max_rel:<14.6e} | {status:<6}")
            print(f"  [DIFF] {name} col={worst}: wrf={f_arr[worst]:.17e} cpp={c_arr[worst]:.17e} diff={abs_diff[worst]:.3e}")
        elif name in ("tsk", "hfx", "qfx", "lh", "grdflx", "albedo", "snow", "snowh", "canwat", "qsfc", "z0", "tvxy", "tgxy", "t2mvxy", "t2mbxy", "q2mvxy", "q2mbxy", "tradxy", "fsaxy", "firaxy"):
            print(f"{name:<16} | {f'({n_cols})':<10} | {active_pts:<10} | {max_abs:<14.6e} | {max_rel:<14.6e} | {status:<6}")

    for name, nz in vars_3d:
        n3d = n_cols * nz
        f_arr = np.frombuffer(f_data, dtype=np.float64, count=n3d, offset=offset)
        c_arr = np.frombuffer(c_data, dtype=np.float64, count=n3d, offset=offset)
        offset += n3d * 8

        abs_diff = np.abs(f_arr - c_arr)
        denom = np.maximum(np.maximum(np.abs(f_arr), np.abs(c_arr)), 1e-12)
        rel_diff = abs_diff / denom

        max_abs = float(np.max(abs_diff))
        max_rel = float(np.max(rel_diff))
        global_max_abs = max(global_max_abs, max_abs)
        active_pts = int(np.count_nonzero((np.abs(f_arr) > 1e-15) & (f_arr != -9999.0)))
        n_checked += 1

        passed = (max_abs <= tol) or (max_rel <= tol)
        status = "PASS" if passed else "FAIL"
        shape_str = f"({n_cols},{nz})"
        print(f"{name:<16} | {shape_str:<10} | {active_pts:<10} | {max_abs:<14.6e} | {max_rel:<14.6e} | {status:<6}")
        if not passed:
            all_pass = False
            worst = int(np.argmax(abs_diff))
            k_w, c_w = divmod(worst, n_cols)
            print(f"  [DIFF] {name} col={c_w} k={k_w}: wrf={f_arr[worst]:.17e} cpp={c_arr[worst]:.17e} diff={abs_diff[worst]:.3e}")

    print("-" * 82)
    verdict = "PASS" if all_pass else "FAIL"
    print(f"Checked {n_checked} variables (62 2D INOUT + 77 2D OUT + 11 3D): Global Max Abs Diff = {global_max_abs:.6e} -> {verdict}")
    return all_pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fortran_file", type=str)
    parser.add_argument("cpp_file", type=str)
    parser.add_argument("--tol", type=float, default=0.0)
    args = parser.parse_args()
    ok = compare_noahmp_outputs(args.fortran_file, args.cpp_file, args.tol)
    sys.exit(0 if ok else 1)
