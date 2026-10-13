#!/usr/bin/env python3
"""Extract 2D slice reference fields from WRF em_hill2d_x wrfout NetCDF file.

Extracts W, U, T, PH at each history time (j=0 slice), plus HGT, ZNW, MUB, PHB at t=0
to a numpy .npz file.
"""

import argparse
import os
import sys
import netCDF4 as nc
import numpy as np


def extract_hill_ref(wrfout_path: str, out_npz_path: str):
    if not os.path.exists(wrfout_path):
        raise FileNotFoundError(f"wrfout file not found: {wrfout_path}")

    ds = nc.Dataset(wrfout_path, "r")

    # Time dimension
    nt = len(ds.dimensions["Time"])

    # Extract 3D/4D fields on j=0 slice
    # W: (Time, bottom_top_stag, south_north, west_east)
    w = ds.variables["W"][:, :, 0, :]
    # U: (Time, bottom_top, south_north, west_east_stag)
    u = ds.variables["U"][:, :, 0, :]
    # T: (Time, bottom_top, south_north, west_east)
    t = ds.variables["T"][:, :, 0, :]
    # PH: (Time, bottom_top_stag, south_north, west_east)
    ph = ds.variables["PH"][:, :, 0, :]

    # Static / initial fields at t=0
    # HGT: (Time, south_north, west_east)
    hgt = ds.variables["HGT"][0, 0, :]
    # ZNW: (Time, bottom_top_stag)
    znw = ds.variables["ZNW"][0, :]
    # MUB: (Time, south_north, west_east)
    mub = ds.variables["MUB"][0, 0, :]
    # PHB: (Time, bottom_top_stag, south_north, west_east)
    phb = ds.variables["PHB"][0, :, 0, :]

    # Times strings if present
    times = []
    if "Times" in ds.variables:
        raw_times = ds.variables["Times"][:]
        times = [b"".join(row).decode("utf-8").strip() for row in raw_times]

    os.makedirs(os.path.dirname(os.path.abspath(out_npz_path)), exist_ok=True)

    np.savez(
        out_npz_path,
        W=w,
        U=u,
        T=t,
        PH=ph,
        HGT=hgt,
        ZNW=znw,
        MUB=mub,
        PHB=phb,
        times=np.array(times),
    )

    print(f"Successfully extracted {nt} timesteps to {out_npz_path}")
    print(f"  W shape:   {w.shape}")
    print(f"  U shape:   {u.shape}")
    print(f"  T shape:   {t.shape}")
    print(f"  PH shape:  {ph.shape}")
    print(f"  HGT shape: {hgt.shape}")
    print(f"  ZNW shape: {znw.shape}")
    print(f"  MUB shape: {mub.shape}")
    print(f"  PHB shape: {phb.shape}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract 2D reference fields from em_hill2d_x wrfout."
    )
    parser.add_argument("wrfout", help="Path to wrfout file")
    parser.add_argument(
        "--out",
        default="/home/ostriquetrum/wrf_ref/run_hill/ref_fields.npz",
        help="Path to output .npz file",
    )
    args = parser.parse_args()
    extract_hill_ref(args.wrfout, args.out)


if __name__ == "__main__":
    main()
