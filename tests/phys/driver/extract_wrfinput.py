#!/usr/bin/env python3
"""
Extracts surface & soil initialization variables from wrfinput netCDF
for the C++ physics initialization smoke test.
"""
import argparse
import sys
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Extract wrfinput variables to binary for C++ test")
    parser.add_argument("--nc", required=True, help="Input wrfinput netCDF file")
    parser.add_argument("--out", required=True, help="Output binary file")
    args = parser.parse_args()

    try:
        import netCDF4 as nc
    except ImportError:
        print("netCDF4 not found", file=sys.stderr)
        sys.exit(1)

    ds = nc.Dataset(args.nc)

    nx = len(ds.dimensions['west_east'])
    ny = len(ds.dimensions['south_north'])
    num_soil_layers = len(ds.dimensions['soil_layers_stag']) if 'soil_layers_stag' in ds.dimensions else 4

    julday = int(getattr(ds, 'JULDAY', 285))
    cen_lat = float(getattr(ds, 'CEN_LAT', 39.0))
    mminlu = getattr(ds, 'MMINLU', 'MODIFIED_IGBP_MODIS_NOAH').strip()

    print(f"Extracting {args.nc}: nx={nx}, ny={ny}, nz_soil={num_soil_layers}, julday={julday}, cen_lat={cen_lat}, mminlu={mminlu}")

    # Read variables (time index 0)
    lu_index = np.array(ds['LU_INDEX'][0, :, :], dtype=np.float64)   # (ny, nx)
    ivgtyp   = np.array(ds['IVGTYP'][0, :, :], dtype=np.float64)     # (ny, nx)
    isltyp   = np.array(ds['ISLTYP'][0, :, :], dtype=np.float64)     # (ny, nx)
    xland    = np.array(ds['XLAND'][0, :, :], dtype=np.float64)      # (ny, nx)
    seaice   = np.array(ds['SEAICE'][0, :, :], dtype=np.float64)     # (ny, nx)
    tsk      = np.array(ds['TSK'][0, :, :], dtype=np.float64)        # (ny, nx)

    # 3D soil fields (time, soil_layers_stag, south_north, west_east)
    if 'TSLB' in ds.variables:
        tslb = np.array(ds['TSLB'][0, :, :, :], dtype=np.float64)    # (nz_soil, ny, nx)
    else:
        tslb = np.full((num_soil_layers, ny, nx), 280.0, dtype=np.float64)

    if 'SMOIS' in ds.variables:
        smois = np.array(ds['SMOIS'][0, :, :, :], dtype=np.float64)  # (nz_soil, ny, nx)
    else:
        smois = np.full((num_soil_layers, ny, nx), 0.25, dtype=np.float64)

    # Reference fields from wrfinput if present
    albbck = np.array(ds['ALBBCK'][0, :, :], dtype=np.float64) if 'ALBBCK' in ds.variables else np.zeros((ny, nx), dtype=np.float64)
    snoalb = np.array(ds['SNOALB'][0, :, :], dtype=np.float64) if 'SNOALB' in ds.variables else np.zeros((ny, nx), dtype=np.float64)

    with open(args.out, "wb") as f:
        # Header: nx, ny, num_soil_layers, julday (int32)
        np.array([nx, ny, num_soil_layers, julday], dtype=np.int32).tofile(f)
        # cen_lat (float64)
        np.array([cen_lat], dtype=np.float64).tofile(f)

        # 2D arrays: in C order row-major (ny, nx) -> element (x, y) = buffer[y * nx + x]
        lu_index.tofile(f)
        ivgtyp.tofile(f)
        isltyp.tofile(f)
        xland.tofile(f)
        seaice.tofile(f)
        tsk.tofile(f)

        # 3D arrays: (nz_soil, ny, nx) -> element (x, y, k) = buffer[x + y * nx + k * nx * ny]
        # In python (nz_soil, ny, nx).tofile() writes in C order (k, y, x) where x is fastest!
        tslb.tofile(f)
        smois.tofile(f)

        # Reference fields
        albbck.tofile(f)
        snoalb.tofile(f)

    print(f"Wrote {args.out} successfully.")

if __name__ == "__main__":
    main()
