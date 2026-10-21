#!/usr/bin/env bash
# Builds the reference for tests/phys/suite: UNMODIFIED WRF v4.8.0, em_real, serial gfortran,
# real*8 (configure -r8 = -fdefault-real-8) with configure.wrf's own flags (FCOPTIM -O2 -ftree-vectorize
# -funroll-loops, as wrf.exe): the C++ follows the loops gfortran vectorises with glibc's libmvec
# (_ZGVbN2vv_pow, _ZGVbN2v_exp, ...). The references of this repository are built with gfortran 14.2;
# another gfortran may vectorise other loops. No WRF source file is touched.
#
#   tests/phys/suite/build_wrf_r8.sh <dir>       -> <dir>/WRF (about 30 min on 4 cores)
# Needs gfortran, gcc, m4, csh, perl and NetCDF-Fortran (apt: gfortran libnetcdff-dev csh m4).
set -euo pipefail
DEST="${1:?usage: build_wrf_r8.sh <dir>}"
mkdir -p "$DEST"
cd "$DEST"
[ -d WRF ] || git clone -q --depth 1 --recurse-submodules --shallow-submodules --branch v4.8.0 https://github.com/wrf-model/WRF.git WRF
cd WRF
if [ -z "${NETCDF:-}" ]; then  # Debian/Ubuntu split NetCDF-C/Fortran: one prefix for configure
  NETCDF="$DEST/netcdf"
  mkdir -p "$NETCDF/include" "$NETCDF/lib"
  ln -sf /usr/include/netcdf* "$NETCDF/include/"
  ln -sf /usr/lib/x86_64-linux-gnu/libnetcdf* "$NETCDF/lib/"
fi
export NETCDF NETCDF_classic=1
printf '32\n1\n' | ./configure -r8 > configure.log 2>&1   # 32 = GNU serial, 1 = basic nesting
[ -f share/landread.c ] || cp share/landread.c.dist share/landread.c
J="-j ${JOBS:-4}" ./compile em_real > compile.log 2>&1
test -x main/wrf.exe
echo "WRF r8 reference built in $DEST/WRF"
