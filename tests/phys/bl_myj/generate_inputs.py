#!/usr/bin/env python3
"""
Generate realistic atmospheric columns from WRF run_P4_r8 / run_P1_r8 wrfout for bl_myj testing.
Produces real columns + synthetic edge cases (including nonzero terrain height HT and nonzero CT).
"""

import argparse
import os
import struct
import numpy as np
from netCDF4 import Dataset

def generate_bl_myj_inputs(wrfout_path, out_file, total_cols=2500):
    if not os.path.exists(wrfout_path):
        wrfout_path = os.path.join(os.environ.get("CPPWRF_WRF_RUNS", ""), "run_P1_r8/wrfout_d01_2026-10-12_00:00:00")
    ds = Dataset(wrfout_path)

    p_pert = ds.variables["P"][0]
    pb = ds.variables["PB"][0]
    pmid_grid = p_pert + pb

    theta_pert = ds.variables["T"][0]
    th_grid = theta_pert + 300.0

    qv_grid = np.maximum(ds.variables["QVAPOR"][0], 1e-12)
    qc_grid = np.maximum(ds.variables["QCLOUD"][0], 0.0) if "QCLOUD" in ds.variables else np.zeros_like(qv_grid)
    if "QICE" in ds.variables:
        qi_grid = np.maximum(ds.variables["QICE"][0], 0.0)
    else:
        qi_grid = np.zeros_like(qc_grid)

    ph = ds.variables["PH"][0]
    phb = ds.variables["PHB"][0]
    geopot = (ph + phb) / 9.81

    dz_grid = geopot[1:, :, :] - geopot[:-1, :, :]

    u_staggered = ds.variables["U"][0]
    v_staggered = ds.variables["V"][0]
    u_grid = 0.5 * (u_staggered[:, :, :-1] + u_staggered[:, :, 1:])
    v_grid = 0.5 * (v_staggered[:, :-1, :] + v_staggered[:, 1:, :])

    ht_grid = ds.variables["HGT"][0]
    tsk_grid = ds.variables["TSK"][0]
    xland_grid = ds.variables["XLAND"][0]
    ustar_grid = np.maximum(ds.variables["UST"][0], 1e-4)
    znt_grid = np.maximum(ds.variables["ZNT"][0], 1e-5)

    if "XICE" in ds.variables:
        sice_grid = ds.variables["XICE"][0]
    elif "SEAICE" in ds.variables:
        sice_grid = ds.variables["SEAICE"][0]
    else:
        sice_grid = np.zeros_like(tsk_grid)

    snow_grid = ds.variables["SNOW"][0] if "SNOW" in ds.variables else np.zeros_like(tsk_grid)
    elflx_grid = ds.variables["LH"][0] if "LH" in ds.variables else np.full_like(tsk_grid, 50.0)
    akhs_grid = np.maximum(ds.variables["AKHS"][0], 1e-4) if "AKHS" in ds.variables else np.full_like(tsk_grid, 0.01)
    akms_grid = np.maximum(ds.variables["AKMS"][0], 1e-4) if "AKMS" in ds.variables else np.full_like(tsk_grid, 0.02)
    ct_grid = ds.variables["CT"][0] if "CT" in ds.variables else np.zeros_like(tsk_grid)

    chklowq_grid = np.ones_like(tsk_grid)
    lowlyr_grid = np.ones_like(tsk_grid, dtype=np.int32)

    nz, ny, nx = dz_grid.shape

    r_d = 287.0
    cp = 1004.5
    capa = r_d / cp
    p1000mb = 100000.0

    exner_grid = (pmid_grid / p1000mb) ** capa
    t_grid = th_grid * exner_grid
    rho_grid = pmid_grid / (r_d * t_grid * (1.0 + 0.608 * qv_grid - qc_grid - qi_grid))

    pint_grid = np.zeros((nz + 1, ny, nx), dtype=np.float64)
    pint_grid[nz, :, :] = pmid_grid[nz - 1, :, :] - 0.5 * rho_grid[nz - 1, :, :] * 9.81 * dz_grid[nz - 1, :, :]
    pint_grid[0, :, :] = pmid_grid[0, :, :] + 0.5 * rho_grid[0, :, :] * 9.81 * dz_grid[0, :, :]
    for k in range(1, nz):
        pint_grid[k, :, :] = 0.5 * (pmid_grid[k - 1, :, :] + pmid_grid[k, :, :])

    tke_grid = np.zeros((nz, ny, nx), dtype=np.float64)
    for k in range(nz):
        tke_grid[k, :, :] = np.maximum(0.5 * (ustar_grid ** 2) * np.exp(-float(k) / 5.0), 0.1)

    exch_h_grid = np.zeros((nz, ny, nx), dtype=np.float64)

    psfc_grid = pint_grid[0, :, :]
    thz0_grid = tsk_grid * (p1000mb / psfc_grid) ** capa
    qz0_grid = qv_grid[0, :, :]
    qsfc_grid = qv_grid[0, :, :]
    uz0_grid = np.zeros_like(tsk_grid)
    vz0_grid = np.zeros_like(tsk_grid)

    n_real = min(total_cols - 202, ny * nx)
    n_synth = total_cols - n_real

    dz_flat = dz_grid.reshape(nz, -1)[:, :n_real]
    pmid_flat = pmid_grid.reshape(nz, -1)[:, :n_real]
    pint_flat = pint_grid.reshape(nz + 1, -1)[:, :n_real]
    th_flat = th_grid.reshape(nz, -1)[:, :n_real]
    t_flat = t_grid.reshape(nz, -1)[:, :n_real]
    exner_flat = exner_grid.reshape(nz, -1)[:, :n_real]
    qv_flat = qv_grid.reshape(nz, -1)[:, :n_real]
    qc_flat = qc_grid.reshape(nz, -1)[:, :n_real]
    qi_flat = qi_grid.reshape(nz, -1)[:, :n_real]
    u_flat = u_grid.reshape(nz, -1)[:, :n_real]
    v_flat = v_grid.reshape(nz, -1)[:, :n_real]
    rho_flat = rho_grid.reshape(nz, -1)[:, :n_real]
    tke_flat = tke_grid.reshape(nz, -1)[:, :n_real]
    exch_h_flat = exch_h_grid.reshape(nz, -1)[:, :n_real]

    ht_flat = ht_grid.reshape(-1)[:n_real]
    tsk_flat = tsk_grid.reshape(-1)[:n_real]
    xland_flat = xland_grid.reshape(-1)[:n_real]
    sice_flat = sice_grid.reshape(-1)[:n_real]
    snow_flat = snow_grid.reshape(-1)[:n_real]
    chklowq_flat = chklowq_grid.reshape(-1)[:n_real]
    lowlyr_flat = lowlyr_grid.reshape(-1)[:n_real]
    ustar_flat = ustar_grid.reshape(-1)[:n_real]
    znt_flat = znt_grid.reshape(-1)[:n_real]
    akhs_flat = akhs_grid.reshape(-1)[:n_real]
    akms_flat = akms_grid.reshape(-1)[:n_real]
    elflx_flat = elflx_grid.reshape(-1)[:n_real]
    qsfc_flat = qsfc_grid.reshape(-1)[:n_real]
    thz0_flat = thz0_grid.reshape(-1)[:n_real]
    qz0_flat = qz0_grid.reshape(-1)[:n_real]
    uz0_flat = uz0_grid.reshape(-1)[:n_real]
    vz0_flat = vz0_grid.reshape(-1)[:n_real]
    ct_flat = ct_grid.reshape(-1)[:n_real]

    rng = np.random.RandomState(42)
    s_dz = np.repeat(dz_flat[:, :1], n_synth, axis=1)
    s_pmid = np.repeat(pmid_flat[:, :1], n_synth, axis=1)
    s_pint = np.repeat(pint_flat[:, :1], n_synth, axis=1)
    s_th = np.repeat(th_flat[:, :1], n_synth, axis=1)
    s_t = np.repeat(t_flat[:, :1], n_synth, axis=1)
    s_exner = np.repeat(exner_flat[:, :1], n_synth, axis=1)
    s_qv = np.repeat(qv_flat[:, :1], n_synth, axis=1)
    s_qc = np.repeat(qc_flat[:, :1], n_synth, axis=1)
    s_qi = np.repeat(qi_flat[:, :1], n_synth, axis=1)
    s_u = np.repeat(u_flat[:, :1], n_synth, axis=1)
    s_v = np.repeat(v_flat[:, :1], n_synth, axis=1)
    s_rho = np.repeat(rho_flat[:, :1], n_synth, axis=1)
    s_tke = np.repeat(tke_flat[:, :1], n_synth, axis=1)
    s_exch_h = np.zeros((nz, n_synth), dtype=np.float64)

    s_ht = rng.uniform(0.0, 2500.0, n_synth)
    s_tsk = np.full(n_synth, 290.0)
    s_xland = np.where(rng.rand(n_synth) > 0.5, 2.0, 1.0)
    s_sice = np.where(rng.rand(n_synth) > 0.8, 1.0, 0.0)
    s_snow = np.where(rng.rand(n_synth) > 0.8, 50.0, 0.0)
    s_chklowq = np.where(rng.rand(n_synth) > 0.2, 1.0, 0.0)
    s_lowlyr = np.ones(n_synth, dtype=np.int32)
    s_ustar = rng.uniform(0.01, 1.5, n_synth)
    s_znt = rng.uniform(1e-4, 0.5, n_synth)
    s_akhs = rng.uniform(0.001, 0.05, n_synth)
    s_akms = rng.uniform(0.001, 0.08, n_synth)
    s_elflx = rng.uniform(-20.0, 400.0, n_synth)
    s_qsfc = rng.uniform(0.001, 0.02, n_synth)
    s_thz0 = s_tsk * 1.0
    s_qz0 = s_qsfc * 1.0
    s_uz0 = np.zeros(n_synth)
    s_vz0 = np.zeros(n_synth)
    s_ct = np.where(rng.rand(n_synth) > 0.7, 0.5, 0.0)

    s_u[:, :20] = 0.0
    s_v[:, :20] = 0.0
    s_u[:, 20:40] = 45.0
    s_v[:, 20:40] = 30.0
    for k in range(nz):
        s_th[k, 40:60] = 270.0 + k * 1.5
    for k in range(nz):
        s_th[k, 60:80] = 310.0 - k * 0.2
    for k in range(nz):
        s_qc[k, 80:100] = 2e-3 * np.exp(-float(k - 5) ** 2 / 10.0)
        s_qi[k, 80:100] = 1e-3 * np.exp(-float(k - 15) ** 2 / 12.0)
    s_tke[:, 100:120] = 10.0
    s_tke[:, 120:140] = 0.05

    all_dz = np.hstack([dz_flat, s_dz])
    all_pmid = np.hstack([pmid_flat, s_pmid])
    all_pint = np.hstack([pint_flat, s_pint])
    all_th = np.hstack([th_flat, s_th])
    all_t = np.hstack([t_flat, s_t])
    all_exner = np.hstack([exner_flat, s_exner])
    all_qv = np.hstack([qv_flat, s_qv])
    all_qc = np.hstack([qc_flat, s_qc])
    all_qi = np.hstack([qi_flat, s_qi])
    all_u = np.hstack([u_flat, s_u])
    all_v = np.hstack([v_flat, s_v])
    all_rho = np.hstack([rho_flat, s_rho])
    all_tke = np.hstack([tke_flat, s_tke])
    all_exch_h = np.hstack([exch_h_flat, s_exch_h])

    all_ht = np.concatenate([ht_flat, s_ht])
    all_tsk = np.concatenate([tsk_flat, s_tsk])
    all_xland = np.concatenate([xland_flat, s_xland])
    all_sice = np.concatenate([sice_flat, s_sice])
    all_snow = np.concatenate([snow_flat, s_snow])
    all_chklowq = np.concatenate([chklowq_flat, s_chklowq])
    all_lowlyr = np.concatenate([lowlyr_flat, s_lowlyr])
    all_ustar = np.concatenate([ustar_flat, s_ustar])
    all_znt = np.concatenate([znt_flat, s_znt])
    all_akhs = np.concatenate([akhs_flat, s_akhs])
    all_akms = np.concatenate([akms_flat, s_akms])
    all_elflx = np.concatenate([elflx_flat, s_elflx])
    all_qsfc = np.concatenate([qsfc_flat, s_qsfc])
    all_thz0 = np.concatenate([thz0_flat, s_thz0])
    all_qz0 = np.concatenate([qz0_flat, s_qz0])
    all_uz0 = np.concatenate([uz0_flat, s_uz0])
    all_vz0 = np.concatenate([vz0_flat, s_vz0])
    all_ct = np.concatenate([ct_flat, s_ct])

    n_cols = total_cols
    dt = 12.0
    stepbl = 1

    print(f"Writing {n_cols} columns (nz={nz}) to {out_file}...")
    with open(out_file, "wb") as f:
        f.write(struct.pack("iii", n_cols, nz, stepbl))
        f.write(struct.pack("d", dt))

        def write_arr(arr, dtype=np.float64):
            np.asarray(arr, dtype=dtype).tofile(f)

        write_arr(all_dz)
        write_arr(all_pmid)
        write_arr(all_pint)
        write_arr(all_th)
        write_arr(all_t)
        write_arr(all_exner)
        write_arr(all_qv)
        write_arr(all_qc)
        write_arr(all_qi)
        write_arr(all_u)
        write_arr(all_v)
        write_arr(all_rho)
        write_arr(all_tke)
        write_arr(all_exch_h)

        write_arr(all_ht)
        write_arr(all_tsk)
        write_arr(all_xland)
        write_arr(all_sice)
        write_arr(all_snow)
        write_arr(all_chklowq)
        write_arr(all_lowlyr, dtype=np.int32)
        write_arr(all_ustar)
        write_arr(all_znt)
        write_arr(all_akhs)
        write_arr(all_akms)
        write_arr(all_elflx)

        write_arr(all_qsfc)
        write_arr(all_thz0)
        write_arr(all_qz0)
        write_arr(all_uz0)
        write_arr(all_vz0)
        write_arr(all_ct)

    ds.close()
    print(f"Successfully generated {out_file} ({os.path.getsize(out_file)} bytes)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cols", type=int, default=2500)
    parser.add_argument("--out", type=str, default="bl_myj_input.bin")
    parser.add_argument("--wrfout", type=str, default=os.path.join(os.environ.get("CPPWRF_WRF_RUNS", ""), "run_P4_r8/wrfout_d01_2026-10-12_00:00:00"))
    args = parser.parse_args()
    generate_bl_myj_inputs(args.wrfout, args.out, args.cols)
