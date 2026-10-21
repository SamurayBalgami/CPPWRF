#!/usr/bin/env bash
# Column differential test of the MYNN surface layer (sf_sfclay_physics = 5) against WRF v4.8.0's own
# compiled mynnsfc_driver (${WRF_DIR}/main/libwrflib.a, a -r8 build). Pass: bit-identical (tolerance 0).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
[ -n "${WRF_DIR:-}" ] || { echo "SKIP: set WRF_DIR to a WRF v4.8.0 -r8 build"; exit 77; }
OUT_DIR="${1:-${PWD}/test_sf_mynn}"
mkdir -p "$OUT_DIR"; cd "$OUT_DIR"
gfortran -O2 -ffree-form -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 \
  -I"$WRF_DIR/main" -I"$WRF_DIR/phys" -I"$WRF_DIR/frame" -I"$WRF_DIR/share" \
  "$SCRIPT_DIR/driver_sf_mynn.f90" "$WRF_DIR/main/libwrflib.a" -lmvec -lm -o driver_sf_mynn_fortran
g++ -std=c++20 -DCPPWRF_DOUBLE -O3 -march=native -ffp-contract=off -I"$REPO_ROOT/src" \
  "$SCRIPT_DIR/driver_sf_mynn.cpp" "$REPO_ROOT/src/core/backend.cpp" -lmvec -lm -o driver_sf_mynn_cpu_fp64
fail=0
for cfg in "1 1 1" "2 1 1" "2 5 1" "2 3 0"; do
  set -- $cfg
  python3 "$SCRIPT_DIR/generate_inputs.py" --cols 4000 --itimestep $1 --nsteps $2 --isfflx $3 --seed $1$2 in.bin
  ./driver_sf_mynn_fortran in.bin wrf.bin
  ./driver_sf_mynn_cpu_fp64 in.bin cpp.bin
  echo "== itimestep $1, $2 calls, isfflx $3"
  python3 "$SCRIPT_DIR/compare.py" wrf.bin cpp.bin 4000 || fail=1
done
[ $fail = 0 ] && echo "PASS" || { echo "FAIL"; exit 1; }
