#!/usr/bin/env bash
# Builds the reference for tests/phys/suite: UNMODIFIED WRF v4.6.1, em_real, serial gfortran,
# real*8 (configure -r8 = -fdefault-real-8), the configure.wrf flags with one change:
# FCOPTIM -O2 -ftree-vectorize -funroll-loops -> -O2 -fno-tree-vectorize -funroll-loops -ffp-contract=off
# (gfortran -ftree-vectorize calls glibc's libmvec vector math, _ZGVbN2vv_pow / _ZGVbN2v_exp /
# _ZGVbN2v_log, whose results differ from scalar libm by ulps; all differential tests of this
# repository compile the Fortran with -fno-tree-vectorize). No WRF source file is touched.
#
#   tests/phys/suite/build_wrf_r8.sh <dir>       -> <dir>/WRF (about 30 min on 4 cores)
# Needs gfortran, gcc, m4, csh, perl and NetCDF-Fortran (apt: gfortran libnetcdff-dev csh m4).
set -euo pipefail
DEST="${1:?usage: build_wrf_r8.sh <dir>}"
mkdir -p "$DEST"
cd "$DEST"
[ -d WRF ] || git clone -q --depth 1 --branch v4.6.1 https://github.com/wrf-model/WRF.git WRF
cd WRF
if [ -z "${NETCDF:-}" ]; then  # Debian/Ubuntu split NetCDF-C/Fortran: one prefix for configure
  NETCDF="$DEST/netcdf"
  mkdir -p "$NETCDF/include" "$NETCDF/lib"
  ln -sf /usr/include/netcdf* "$NETCDF/include/"
  ln -sf /usr/lib/x86_64-linux-gnu/libnetcdf* "$NETCDF/lib/"
fi
export NETCDF NETCDF_classic=1
printf '32\n1\n' | ./configure -r8 > configure.log 2>&1   # 32 = GNU serial, 1 = basic nesting
sed -i 's/^FCOPTIM *=.*/FCOPTIM         =       -O2 -fno-tree-vectorize -funroll-loops -ffp-contract=off/' configure.wrf
[ -f share/landread.c ] || cp share/landread.c.dist share/landread.c
J="-j ${JOBS:-4}" ./compile em_real > compile.log 2>&1
test -x main/wrf.exe
if nm dyn_em/module_big_step_utilities_em.o | grep -q _ZGV; then
  echo "error: vector math calls remain in the WRF objects" >&2
  exit 1
fi
echo "WRF r8 reference built in $DEST/WRF"
