#!/usr/bin/env bash
# Column differential test of the MYNN-EDMF PBL scheme (bl_pbl_physics = 5, WRF's default bl_mynn_* options)
# against WRF v4.8.0's own compiled mynnedmf_driver (${WRF_DIR}/main/libwrflib.a of a -r8 -fdefault-double-8
# build: its DOUBLE PRECISION locals are REAL*16 in a plain -r8 build). Pass: bit-identical (tolerance 0).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
[ -n "${WRF_DIR:-}" ] || { echo "SKIP: set WRF_DIR to a WRF v4.8.0 -r8 -fdefault-double-8 build"; exit 77; }
OUT_DIR="${1:-${PWD}/test_bl_mynn}"
mkdir -p "$OUT_DIR"; cd "$OUT_DIR"
gfortran -O2 -ffree-form -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 \
  -I"$WRF_DIR/main" -I"$WRF_DIR/phys" -I"$WRF_DIR/frame" -I"$WRF_DIR/share" \
  "$SCRIPT_DIR/driver_bl_mynn.f90" "$WRF_DIR/main/libwrflib.a" -lmvec -lm -o driver_bl_mynn_fortran
g++ -std=c++20 -DCPPWRF_DOUBLE -O3 -march=native -ffp-contract=off -I"$REPO_ROOT/src" \
  "$SCRIPT_DIR/driver_bl_mynn.cpp" "$REPO_ROOT/src/core/backend.cpp" -lmvec -lm -o driver_bl_mynn_cpu_fp64
fail=0
for cfg in "50 1 60 12000 1" "50 5 60 12000 2" "51 5 60 3000 3" "40 8 20 1000 4" "60 4 180 27000 5" "35 6 10 500 6"; do
  set -- $cfg
  python3 "$SCRIPT_DIR/generate_inputs.py" --cols 2000 --nz $1 --nsteps $2 --dt $3 --dx $4 --seed $5 in.bin
  ./driver_bl_mynn_fortran in.bin wrf.bin
  ./driver_bl_mynn_cpu_fp64 in.bin cpp.bin
  echo "== nz $1, $2 calls, dt $3, dx $4"
  python3 "$SCRIPT_DIR/compare.py" wrf.bin cpp.bin 2000 $1 || fail=1
done
[ $fail = 0 ] && echo "PASS" || { echo "FAIL"; exit 1; }
