#!/usr/bin/env python3
"""
Generate input columns for Eta/MYJ surface layer (sf_myjsfc) differential test.
"""

import argparse
import os
import struct
import numpy as np

def generate_inputs(out_file, n_cols=2500, itimestep=2, wrfout_path=None):
    if wrfout_path is None:
        wrfout_path = os.path.join(os.environ.get("CPPWRF_WRF_RUNS", ""), "run_P4_r8/wrfout_d01_2026-10-12_00:00:00")
    if not os.path.exists(wrfout_path):
        wrfout_path = os.path.join(os.environ.get("CPPWRF_WRF_RUNS", ""), "run_P1_r8/wrfout_d01_2026-10-12_00:00:00")

    nz = 44
    dz_list, pmid_list, pint_list, th_list, t_list = [], [], [], [], []
    qv_list, qc_list, u_list, v_list, q2_list = [], [], [], [], []
    ht_list, tsk_list, xland_list, z0base_list, mavail_list = [], [], [], [], []
    lowlyr_list, ivgtyp_list = [], []
    qsfc_list, thz0_list, qz0_list, uz0_list, vz0_list = [], [], [], [], []
    ustar_list, znt_list, pblh_list, rmol_list, akhs_list, akms_list = [], [], [], [], [], []

    if not os.path.exists(wrfout_path):
        raise FileNotFoundError(f"{wrfout_path} (set CPPWRF_WRF_RUNS to the WRF reference runs)")
    try:
        import netCDF4 as nc
        if os.path.exists(wrfout_path):
            ds = nc.Dataset(wrfout_path)
            ny = ds.dimensions['south_north'].size
            nx = ds.dimensions['west_east'].size
            nz = ds.dimensions['bottom_top'].size

            ntimes = ds.dimensions['Time'].size
            times = [0, ntimes // 2, ntimes - 1]
            pts_per_time = (n_cols - 200) // len(times)

            np.random.seed(42)

            for tidx in times:
                U = ds.variables['U'][tidx]
                V = ds.variables['V'][tidx]
                T = ds.variables['T'][tidx]
                QVAPOR = ds.variables['QVAPOR'][tidx]
                QCLOUD = ds.variables['QCLOUD'][tidx] if 'QCLOUD' in ds.variables else np.zeros_like(QVAPOR)
                P = ds.variables['P'][tidx]
                PB = ds.variables['PB'][tidx]
                PH = ds.variables['PH'][tidx]
                PHB = ds.variables['PHB'][tidx]
                TSK = ds.variables['TSK'][tidx]
                PSFC = ds.variables['PSFC'][tidx]
                XLAND = ds.variables['XLAND'][tidx]
                ZNT = ds.variables['ZNT'][tidx]
                UST = ds.variables['UST'][tidx]
                PBLH = ds.variables['PBLH'][tidx]
                QSFC_in = ds.variables['QSFC'][tidx]

                p_full = P + PB
                th_full = T + 300.0
                ph_full = (PH + PHB) / 9.81
                t_abs = th_full * ((p_full / 100000.0) ** (287.0 / 1004.5))

                u_mass = 0.5 * (U[:, :, :-1] + U[:, :, 1:])
                v_mass = 0.5 * (V[:, :-1, :] + V[:, 1:, :])

                js = np.random.randint(0, ny, pts_per_time)
                is_ = np.random.randint(0, nx, pts_per_time)

                for j, i in zip(js, is_):
                    zh = ph_full[:, j, i]
                    dz = np.maximum(np.diff(zh), 1.0)

                    dz_list.append(dz)
                    pmid_list.append(p_full[:, j, i])

                    pint = np.zeros(nz + 1, dtype=np.float64)
                    pint[0] = PSFC[j, i]
                    for k in range(1, nz):
                        pint[k] = 0.5 * (p_full[k-1, j, i] + p_full[k, j, i])
                    pint[nz] = max(p_full[nz-1, j, i] - 0.5 * (p_full[nz-2, j, i] - p_full[nz-1, j, i]), 1000.0)
                    pint_list.append(pint)

                    th_list.append(th_full[:, j, i])
                    t_list.append(t_abs[:, j, i])
                    qv_list.append(np.maximum(QVAPOR[:, j, i], 1e-8))
                    qc_list.append(np.maximum(QCLOUD[:, j, i], 0.0))
                    u_list.append(u_mass[:, j, i])
                    v_list.append(v_mass[:, j, i])

                    tke = np.zeros(nz, dtype=np.float64)
                    for k in range(nz):
                        z_mid = zh[k] - zh[0]
                        pbl_val = max(PBLH[j, i], 100.0)
                        if z_mid < pbl_val:
                            tke[k] = max(UST[j, i] ** 2 * (1.0 - z_mid / pbl_val) ** 2 * 3.75, 0.1)
                        else:
                            tke[k] = 0.1
                    q2_list.append(tke)

                    ht_list.append(zh[0])
                    tsk_list.append(TSK[j, i])
                    xland_list.append(XLAND[j, i])
                    z0base_list.append(ZNT[j, i])
                    mavail_list.append(1.0 if XLAND[j, i] > 1.5 else 0.5)
                    lowlyr_list.append(1)
                    ivgtyp_list.append(16 if XLAND[j, i] > 1.5 else 7)

                    qsfc_list.append(QSFC_in[j, i])
                    thz0_list.append(th_full[0, j, i])
                    qz0_list.append(qv_list[-1][0])
                    uz0_list.append(0.0)
                    vz0_list.append(0.0)
                    ustar_list.append(max(UST[j, i], 0.05))
                    znt_list.append(ZNT[j, i])
                    pblh_list.append(max(PBLH[j, i], 50.0))
                    rmol_list.append(0.0)
                    akhs_list.append(0.0)
                    akms_list.append(0.0)

            ds.close()
    except Exception as e:
        raise RuntimeError(f"cannot read the real columns from {wrfout_path}") from e

    remaining = n_cols - len(ht_list)
    np.random.seed(12345)

    for c in range(remaining):
        scenario = c % 5
        z_interfaces = np.zeros(nz + 1)
        z_interfaces[0] = 0.0 if scenario != 4 else 3500.0
        cur_dz = 20.0
        for k in range(nz):
            z_interfaces[k+1] = z_interfaces[k] + cur_dz
            cur_dz = min(cur_dz * 1.1, 800.0)
        dz = np.diff(z_interfaces)
        ht = z_interfaces[0]

        psfc = 101325.0 if scenario != 4 else 65000.0
        tsk = 295.0
        xland = 2.0 if scenario in [0, 3] else 1.0
        z0base = 0.001 if xland > 1.5 else 0.1

        u_sfc = 5.0
        v_sfc = 2.0
        ustar = 0.2

        if scenario == 0:
            tsk, t_sfc, u_sfc, ustar = 280.0, 285.0, 3.0, 0.08
        elif scenario == 1:
            tsk, t_sfc, u_sfc, ustar, z0base = 315.0, 300.0, 4.0, 0.35, 0.15
        elif scenario == 2:
            tsk, t_sfc, u_sfc, v_sfc, ustar, z0base = 270.0, 275.0, 0.1, 0.05, 0.01, 0.05
        elif scenario == 3:
            tsk, t_sfc, u_sfc, v_sfc, ustar, z0base = 290.0, 289.0, 25.0, 15.0, 1.2, 0.005
        else:
            tsk, t_sfc, u_sfc, ustar = 282.0, 280.0, 8.0, 0.25

        t_arr = np.zeros(nz)
        th_arr = np.zeros(nz)
        p_arr = np.zeros(nz)
        pint = np.zeros(nz + 1)
        pint[0] = psfc
        qv_arr = np.zeros(nz)
        qc_arr = np.zeros(nz)
        u_arr = np.zeros(nz)
        v_arr = np.zeros(nz)
        q2_arr = np.zeros(nz)

        cur_p = psfc
        for k in range(nz):
            dp = 1.2 * 9.81 * dz[k]
            p_arr[k] = cur_p - 0.5 * dp
            pint[k+1] = max(cur_p - dp, 5000.0)
            cur_p = pint[k+1]

            t_arr[k] = t_sfc - 0.0065 * (z_interfaces[k+1] - z_interfaces[0])
            th_arr[k] = t_arr[k] * ((100000.0 / p_arr[k]) ** (287.0 / 1004.5))
            qv_arr[k] = max(0.012 * np.exp(-z_interfaces[k] / 2000.0), 1e-6)
            qc_arr[k] = 1e-4 if (k in [3, 4] and scenario == 0) else 0.0

            wind_factor = (z_interfaces[k+1] / 50.0) ** 0.15
            u_arr[k] = u_sfc * wind_factor
            v_arr[k] = v_sfc * wind_factor
            q2_arr[k] = max(2.0 * ustar ** 2 * np.exp(-z_interfaces[k] / 500.0), 0.2)

        dz_list.append(dz)
        pmid_list.append(p_arr)
        pint_list.append(pint)
        th_list.append(th_arr)
        t_list.append(t_arr)
        qv_list.append(qv_arr)
        qc_list.append(qc_arr)
        u_list.append(u_arr)
        v_list.append(v_arr)
        q2_list.append(q2_arr)

        ht_list.append(ht)
        tsk_list.append(tsk)
        xland_list.append(xland)
        z0base_list.append(z0base)
        mavail_list.append(1.0 if xland > 1.5 else 0.3)
        lowlyr_list.append(1)
        ivgtyp_list.append(16 if xland > 1.5 else 4)

        qsfc_list.append(qv_arr[0])
        thz0_list.append(th_arr[0])
        qz0_list.append(qv_arr[0])
        uz0_list.append(0.0)
        vz0_list.append(0.0)
        ustar_list.append(ustar)
        znt_list.append(z0base)
        pblh_list.append(500.0)
        rmol_list.append(0.0)
        akhs_list.append(0.0)
        akms_list.append(0.0)

    n_actual = len(ht_list)
    print(f"Writing {n_actual} columns (nz={nz}, itimestep={itimestep}) to {out_file}...")

    with open(out_file, "wb") as f:
        f.write(struct.pack("iiiii", n_actual, nz, itimestep, 1, 0))
        f.write(struct.pack("d", 100000.0))

        def write_3d(arr_list, levels=nz):
            data = np.zeros((levels, n_actual), dtype=np.float64)
            for c in range(n_actual):
                data[:, c] = arr_list[c][:levels]
            f.write(data.tobytes())

        def write_2d(val_list, dtype=np.float64):
            arr = np.array(val_list, dtype=dtype)
            f.write(arr.tobytes())

        write_3d(dz_list, nz)
        write_3d(pmid_list, nz)
        write_3d(pint_list, nz + 1)
        write_3d(th_list, nz)
        write_3d(t_list, nz)
        write_3d(qv_list, nz)
        write_3d(qc_list, nz)
        write_3d(u_list, nz)
        write_3d(v_list, nz)
        write_3d(q2_list, nz)

        write_2d(ht_list)
        write_2d(tsk_list)
        write_2d(xland_list)
        write_2d(z0base_list)
        write_2d(mavail_list)
        write_2d(lowlyr_list, np.int32)
        write_2d(ivgtyp_list, np.int32)

        write_2d(qsfc_list)
        write_2d(thz0_list)
        write_2d(qz0_list)
        write_2d(uz0_list)
        write_2d(vz0_list)
        write_2d(ustar_list)
        write_2d(znt_list)
        write_2d(pblh_list)
        write_2d(rmol_list)
        write_2d(akhs_list)
        write_2d(akms_list)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cols", type=int, default=2500)
    parser.add_argument("--itimestep", type=int, default=2)
    parser.add_argument("--out", type=str, default="sf_myjsfc_input.bin")
    parser.add_argument("--wrfout", type=str, default=None)
    args = parser.parse_args()
    generate_inputs(args.out, args.cols, args.itimestep, args.wrfout)
