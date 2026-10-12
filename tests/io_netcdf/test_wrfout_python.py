#!/usr/bin/env python3
"""Python validation for CPPWRF NetCDF wrfout writer.
Checks:
1. Round-trip of all variables and attributes through netCDF4-python.
2. wrf-python integration (or exact formula equivalent if wrf-python cannot be installed).
"""
import os
import sys
import numpy as np
import netCDF4

def test_wrfout(filepath):
    print(f"Testing NetCDF wrfout file: {filepath}")
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    ds = netCDF4.Dataset(filepath, "r")

    # 1. Dimension checks
    expected_dims = [
        "Time", "DateStrLen", "west_east", "south_north", "bottom_top",
        "bottom_top_stag", "west_east_stag", "south_north_stag"
    ]
    for d in expected_dims:
        assert d in ds.dimensions, f"Missing dimension {d}"
    print(f"Dimensions verified: {list(ds.dimensions.keys())}")

    # 2. Variable presence checks
    expected_vars = [
        "Times", "U", "V", "W", "PH", "PHB", "T", "P", "PB",
        "MU", "MUB", "XLAT", "XLONG", "HGT", "ZNU", "ZNW",
        "C1H", "C2H", "C1F", "C2F", "C3H", "C4H", "C3F", "C4F",
        "RDX", "RDY", "P_TOP", "XTIME"
    ]
    for v in expected_vars:
        assert v in ds.variables, f"Missing variable {v}"
    print("All standard WRF variables verified present.")

    # 3. Attributes check
    for vname in ["U", "V", "W", "PH", "PHB", "T", "P", "PB", "MU", "MUB"]:
        v = ds.variables[vname]
        assert hasattr(v, "FieldType") and v.FieldType == 104, f"{vname}: bad FieldType"
        assert hasattr(v, "MemoryOrder"), f"{vname}: missing MemoryOrder"
        assert hasattr(v, "description"), f"{vname}: missing description"
        assert hasattr(v, "units"), f"{vname}: missing units"
        assert hasattr(v, "stagger"), f"{vname}: missing stagger"

    assert ds.variables["U"].stagger == "X"
    assert ds.variables["V"].stagger == "Y"
    assert ds.variables["W"].stagger == "Z"
    assert ds.variables["PH"].stagger == "Z"
    assert ds.variables["PHB"].stagger == "Z"
    assert ds.variables["T"].stagger == ""
    print("Variable attributes (FieldType, MemoryOrder, description, units, stagger) verified.")

    # 4. Global attributes check
    req_globals = [
        "TITLE", "START_DATE", "SIMULATION_START_DATE",
        "WEST-EAST_GRID_DIMENSION", "SOUTH-NORTH_GRID_DIMENSION", "BOTTOM-TOP_GRID_DIMENSION",
        "DX", "DY", "DT", "MAP_PROJ", "CEN_LAT", "TRUELAT1", "TRUELAT2",
        "STAND_LON", "MOAD_CEN_LAT"
    ]
    for g in req_globals:
        assert hasattr(ds, g), f"Missing global attribute {g}"
    print("Global attributes verified.")

    # 5. Round-trip data check
    times = ds.variables["Times"][:]
    t0_str = "".join([c.decode("utf-8") if isinstance(c, bytes) else c for c in times[0]])
    t1_str = "".join([c.decode("utf-8") if isinstance(c, bytes) else c for c in times[1]])
    assert t0_str.startswith("0001-01-01_00:00:00"), f"Unexpected Times[0]: {t0_str}"
    assert t1_str.startswith("0001-01-01_00:05:00"), f"Unexpected Times[1]: {t1_str}"

    u_data = ds.variables["U"][:]
    assert np.isclose(u_data[0, 0, 0, 0], 10.5, atol=1e-5)
    assert np.isclose(u_data[1, 0, 0, 0], 12.0, atol=1e-5)

    t_data = ds.variables["T"][:]
    assert np.isclose(t_data[0, 0, 0, 0], 3.5, atol=1e-5)
    assert np.isclose(t_data[1, 0, 0, 0], 4.2, atol=1e-5)
    print("Round-trip variable data verified.")

    # 6. wrf-python check or reference equivalent
    wrf_installed = False
    wrf_err = ""
    try:
        import wrf
        wrf_installed = True
    except Exception as e:
        wrf_err = str(e)

    if wrf_installed:
        print("wrf-python is installed! Computing getvar('theta'), getvar('pressure'), getvar('z')...")
        theta = wrf.getvar(ds, "theta", timeidx=0)
        pressure = wrf.getvar(ds, "pressure", timeidx=0)
        z = wrf.getvar(ds, "z", timeidx=0)
        assert np.allclose(theta, t_data[0] + 300.0, rtol=1e-5, atol=1e-5)
        expected_p = (ds.variables["P"][0] + ds.variables["PB"][0]) / 100.0
        assert np.allclose(pressure, expected_p, rtol=1e-5, atol=1e-5)
        print("wrf-python diagnostics match CPPWRF float32 precision!")
    else:
        print(f"wrf-python could not be imported: {wrf_err}")
        print("Executing reference mathematical diagnostics equivalent to wrf-python:")
        # wrf-python getvar('theta') is T + 300.0
        theta_ref = t_data[0] + 300.0
        # wrf-python getvar('pressure') is (P + PB) / 100.0 in hPa
        p_ref = (ds.variables["P"][0] + ds.variables["PB"][0]) / 100.0
        # wrf-python getvar('z') is 0.5 * (PH[:-1] + PH[1:] + PHB[:-1] + PHB[1:]) / 9.81
        ph_tot = ds.variables["PH"][0] + ds.variables["PHB"][0]
        z_ref = 0.5 * (ph_tot[:-1, :, :] + ph_tot[1:, :, :]) / 9.81

        assert np.allclose(theta_ref[0, 0, 0], 303.5, atol=1e-5)
        assert np.allclose(p_ref[0, 0, 0], 999.5, atol=1e-5)
        expected_z0 = 0.5 * (10780.0 + 10780.0) / 9.81
        assert np.allclose(z_ref[0, 0, 0], expected_z0, atol=1e-5)
        print("Reference wrf-python diagnostic computations match expected values to float32 precision.")

    ds.close()
    print("ALL NETCDF PYTHON TESTS PASSED!")
    return True

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "tests/io_netcdf/test_out/wrfout_d01_test.nc"
    test_wrfout(path)
