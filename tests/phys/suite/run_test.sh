#!/usr/bin/env bash
# P1 physics-suite differential test: ORIGINAL WRF v4.6.1 (real*8) vs the C++ PhysicsSuite.
#
#   WRF_DIR=<WRF built by build_wrf_r8.sh> tests/phys/suite/run_test.sh
#
# 1. generate_inputs.py writes the synthetic 24x20x44 tile (suite_input.bin)
# 2. fortran/suite_ref.F90 is linked against the WRF objects and runs WRF's start_domain and
#    3x (first_rk_step_part1, first_rk_step_part2, phy_prep_part2) -> suite_ref.bin
# 3. driver_suite.cpp runs the C++ suite on the same inputs (CPU FP64)  -> suite_cpp.bin
# 4. compare.py: every field bit-identical or <= 1e-12 relative
# The FP32 build of the same driver is compiled as well (it must compile; not compared).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/../../.." && pwd)"
WRF_DIR="${WRF_DIR:?set WRF_DIR to a WRF v4.6.1 -r8 build (tests/phys/suite/build_wrf_r8.sh)}"
B="$SCRIPT_DIR/build"
HEAVY="$REPO/tools/heavy.sh"
[ -x "$HEAVY" ] || HEAVY=""
mkdir -p "$B"
cd "$B"

echo "== 1. inputs"
python3 -I "$SCRIPT_DIR/generate_inputs.py" suite_input.bin
cp "$SCRIPT_DIR/namelist.input" .
for t in LANDUSE.TBL VEGPARM.TBL SOILPARM.TBL GENPARM.TBL; do cp "$WRF_DIR/run/$t" .; done

echo "== 2. WRF reference (fortran/suite_ref.F90 + unmodified WRF objects)"
W="$WRF_DIR"
INC="-I$W/main -I$W/frame -I$W/share -I$W/phys -I$W/dyn_em -I$W/inc -I$W/external/esmf_time_f90"
$HEAVY gfortran -c -O2 -ffree-form -ffree-line-length-none -fdefault-real-8 -fallow-argument-mismatch \
  $INC "$SCRIPT_DIR/fortran/suite_ref.F90" -o suite_ref.o
NCLIB="$(nf-config --flibs 2>/dev/null || echo "-L${NETCDF:-/usr}/lib -lnetcdff -lnetcdf")"
$HEAVY gfortran -o suite_ref suite_ref.o "$W/main/module_wrf_top.o" "$W/main/libwrflib.a" \
  "$W/external/fftpack/fftpack5/libfftpack.a" "$W/external/io_grib1/libio_grib1.a" \
  "$W/external/io_grib_share/libio_grib_share.a" "$W/external/io_int/libwrfio_int.a" \
  -L"$W/external/esmf_time_f90" -lesmf_time "$W/external/RSL_LITE/librsl_lite.a" \
  "$W/frame/module_internal_header_util.o" "$W/frame/pack_utils.o" \
  -L"$W/external/io_netcdf" -lwrfio_nf $NCLIB
./suite_ref > suite_ref.log 2>&1 || { tail -20 suite_ref.log; exit 1; }
tail -1 suite_ref.log

echo "== 3. C++ suite (CPU FP64, and FP32 compile check)"
SRC=("$SCRIPT_DIR/driver_suite.cpp" "$REPO/src/phys/driver/suite.cpp" "$REPO/src/phys/driver/phy_prep.cpp"
     "$REPO/src/phys/driver/tendency_coupling.cpp" "$REPO/src/phys/sf_sfclayrev/sfclayrev.cpp"
     "$REPO/src/phys/bl_ysu/ysu.cpp" "$REPO/src/phys/lsm_noah/noah.cpp" "$REPO/src/phys/lsm_noah/noah_tables.cpp"
     "$REPO/src/core/backend.cpp" "$REPO/src/core/halo.cpp" "$REPO/src/dyn/state.cpp")
$HEAVY g++ -std=c++20 -DCPPWRF_DOUBLE -O3 -march=native -ffp-contract=off -fopenmp -I "$REPO/src" "${SRC[@]}" -o driver_suite_fp64
$HEAVY g++ -std=c++20 -O3 -march=native -ffp-contract=off -fopenmp -I "$REPO/src" "${SRC[@]}" -o driver_suite_fp32
./driver_suite_fp64 suite_ref.bin suite_cpp.bin .

echo "== 4. compare"
python3 -I "$SCRIPT_DIR/compare.py" suite_ref.bin suite_cpp.bin "$@"
