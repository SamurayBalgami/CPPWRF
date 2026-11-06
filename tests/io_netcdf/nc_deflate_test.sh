#!/usr/bin/env bash
# NetCDF-4 files with the chunks compressed on host threads (src/wrf_io/nc_deflate.hpp) against HDF5's own deflate
# filter: tests/io_netcdf/test_nc_deflate.cpp writes the same file three times, with the precomputed chunks (default),
# with the filter compressing each chunk itself (CPPWRF_NC_THREADS=1) and with HDF5's own filter (CPPWRF_NC_THREADS=0);
# the files must be identical byte for byte.
#   usage: nc_deflate_test.sh TEST_EXE WORK_DIR
set -uo pipefail
EXE=$1 W=$2
rm -rf "$W" && mkdir -p "$W" || exit 1
"$EXE" "$W/pool.nc"; rc=$?
[ $rc = 0 ] || { echo "FAIL: precomputed chunks (exit $rc)"; exit $rc; }
CPPWRF_NC_THREADS=1 "$EXE" "$W/serial.nc" || { echo "FAIL: CPPWRF_NC_THREADS=1"; exit 1; }
CPPWRF_NC_THREADS=0 "$EXE" "$W/hdf5.nc" || { echo "FAIL: HDF5's own filter"; exit 1; }
bad=0
for f in pool serial; do
  if cmp "$W/$f.nc" "$W/hdf5.nc"; then echo "$f.nc: identical to HDF5's own filter ($(stat -c %s "$W/hdf5.nc") bytes)"
  else echo "FAIL: $f.nc differs from HDF5's own filter"; bad=1; fi
done
rm -rf "$W"
exit $bad
