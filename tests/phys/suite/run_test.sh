#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# P1 physics-suite differential test: ORIGINAL WRF v4.8.0 (real*8) vs the C++ PhysicsSuite.
#
#   WRF_DIR=<WRF built by build_wrf_r8.sh> tests/phys/suite/run_test.sh
#
# 1. generate_inputs.py writes the synthetic 24x20x44 tile (suite_input.bin)
# 2. fortran/suite_ref.F90 is linked against the WRF objects and runs WRF's start_domain and
#    3x (first_rk_step_part1, first_rk_step_part2, phy_prep_part2) -> suite_ref.bin
# 3. driver_suite.cpp runs the C++ suite on the same inputs (CPU FP64)  -> suite_cpp.bin
# 4. compare.py: every field bit-identical (--tol 0)
# The FP32 build of the same driver is compiled as well (it must compile; not compared).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/../../.." && pwd)"
. "$SCRIPT_DIR/../wrf_fortran.sh"
wrf_need r8   # tests/phys/suite/build_wrf_r8.sh builds one
B="${DIFF_TEST_BUILD_DIR:-$SCRIPT_DIR/build}/suite"
HEAVY="$REPO/tools/heavy.sh"
[ -x "$HEAVY" ] || HEAVY=""
mkdir -p "$B"
cd "$B"

echo "== 1. inputs"
python3 -I "$SCRIPT_DIR/generate_inputs.py" suite_input.bin
cp "$SCRIPT_DIR/namelist.input" .
for t in LANDUSE.TBL VEGPARM.TBL SOILPARM.TBL GENPARM.TBL; do cp -L "$W/run/$t" .; done

echo "== 2. WRF reference (fortran/suite_ref.F90 + unmodified WRF objects)"
wrf_fortran "$B/suite_ref" "$SCRIPT_DIR/fortran/suite_ref.F90" "$W/main/module_wrf_top.o"
./suite_ref > suite_ref.log 2>&1 || { tail -20 suite_ref.log; exit 1; }
tail -1 suite_ref.log

echo "== 3. C++ suite (CPU FP64, and FP32 compile check)"
# the physics sources of cppwrf_core (CMakeLists.txt) + the core pieces the driver uses
mapfile -t PHYS_SRC < <(find "$REPO/src/phys" -name '*.cpp' ! -name kessler.cpp \
                         ! -name wsm6.cpp | sort)
SRC=("$SCRIPT_DIR/driver_suite.cpp" "${PHYS_SRC[@]}" "$REPO/src/core/backend.cpp" "$REPO/src/core/column_work.cpp" "$REPO/src/core/halo.cpp"
     "$REPO/src/dyn/state.cpp")
$HEAVY g++ -std=c++20 -DCPPWRF_DOUBLE -O3 -march=native -ffp-contract=off -fopenmp -DCPPWRF_DATA_DIR=\"$REPO/data\" -I "$REPO/src" "${SRC[@]}" -lmvec -o driver_suite_fp64
$HEAVY g++ -std=c++20 -O3 -march=native -ffp-contract=off -fopenmp -DCPPWRF_DATA_DIR=\"$REPO/data\" -I "$REPO/src" "${SRC[@]}" -lmvec -o driver_suite_fp32
./driver_suite_fp64 suite_ref.bin suite_cpp.bin .

echo "== 4. compare"
python3 -I "$SCRIPT_DIR/compare.py" suite_ref.bin suite_cpp.bin --tol 0 "$@"
