#!/usr/bin/env python3
"""
Extract realistic coarse fields from WRF run_R0_r8 wrfout for the nesting tests.
"""
import os
import sys
import struct
import numpy as np
import netCDF4 as nc

def main():
    wrfout_path = os.path.expanduser("~/wrf_real/") + "run_R0_r8/wrfout_d01_2026-10-12_00:00:00"
    out_bin_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.environ.get("CPPWRF_SCRATCH", os.path.expanduser("~/cppwrf_scratch")), "nest_test/coarse_fields.bin")

    print(f"Reading {wrfout_path}...")
    ds = nc.Dataset(wrfout_path, 'r')

    # Dimensions
    nx = len(ds.dimensions['west_east'])        # 119
    ny = len(ds.dimensions['south_north'])      # 89
    nz = len(ds.dimensions['bottom_top'])       # 44

    # Variables at time index 0 as float64
    u = np.array(ds.variables['U'][0], dtype=np.float64)           # (44, 89, 120) -> (k, j, i)
    v = np.array(ds.variables['V'][0], dtype=np.float64)           # (44, 90, 119)
    w = np.array(ds.variables['W'][0], dtype=np.float64)           # (45, 89, 119)
    ph = np.array(ds.variables['PH'][0], dtype=np.float64)         # (45, 89, 119)
    t = np.array(ds.variables['T'][0], dtype=np.float64)           # (44, 89, 119)
    mu = np.array(ds.variables['MU'][0], dtype=np.float64)         # (89, 119)
    qvapor = np.array(ds.variables['QVAPOR'][0], dtype=np.float64) # (44, 89, 119)
    hgt = np.array(ds.variables['HGT'][0], dtype=np.float64)       # (89, 119)
    mub = np.array(ds.variables['MUB'][0], dtype=np.float64)       # (89, 119)
    phb = np.array(ds.variables['PHB'][0], dtype=np.float64)       # (45, 89, 119)

    msfuy = np.array(ds.variables['MAPFAC_UY'][0], dtype=np.float64) # (89, 120)
    msfvx = np.array(ds.variables['MAPFAC_VX'][0], dtype=np.float64) # (90, 119)
    msfty = np.array(ds.variables['MAPFAC_MY'][0], dtype=np.float64) # (89, 119)

    c1h = np.array(ds.variables['C1H'][0], dtype=np.float64)       # (44)
    c2h = np.array(ds.variables['C2H'][0], dtype=np.float64)       # (44)
    c3h = np.array(ds.variables['C3H'][0], dtype=np.float64)       # (44)
    c4h = np.array(ds.variables['C4H'][0], dtype=np.float64)       # (44)
    c1f = np.array(ds.variables['C1F'][0], dtype=np.float64)       # (45)
    c2f = np.array(ds.variables['C2F'][0], dtype=np.float64)       # (45)
    p_top = float(ds.variables['P_TOP'][0])

    ds.close()

    print(f"Writing binary output to {out_bin_path}...")
    with open(out_bin_path, 'wb') as f:
        # Header: nx, ny, nz (int32), p_top (float64) with no padding
        f.write(struct.pack('=iiid', nx, ny, nz, p_top))

        def write_3d(arr):
            arr_fort = np.transpose(arr, (2, 0, 1)) # -> (i, k, j)
            arr_fort.tofile(f)

        def write_2d(arr):
            arr_fort = np.transpose(arr, (1, 0)) # -> (i, j)
            arr_fort.tofile(f)

        def write_1d(arr):
            arr.tofile(f)

        write_3d(u)
        write_3d(v)
        write_3d(w)
        write_3d(ph)
        write_3d(t)
        write_2d(mu)
        write_3d(qvapor)
        write_2d(hgt)
        write_2d(mub)
        write_3d(phb)

        write_2d(msfuy)
        write_2d(msfvx)
        write_2d(msfty)

        write_1d(c1h)
        write_1d(c2h)
        write_1d(c3h)
        write_1d(c4h)
        write_1d(c1f)
        write_1d(c2f)

    print("Coarse fields extraction complete.")

if __name__ == '__main__':
    main()
