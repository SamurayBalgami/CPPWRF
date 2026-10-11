#!/usr/bin/env bash
# Build (CPU FP32/FP64, CUDA FP32/FP64 if nvcc + a GPU are available) and run the M1 cases.
# Usage: tools/run_m1.sh [case ...]      (default: rest_state straka warm_bubble igw)
# Outputs: out/<case>_<backend>_<precision>/, verifier summary via tests/verify/verify.py.
set -euo pipefail
cd "$(dirname "$0")/.."
CASES=("$@")
[ ${#CASES[@]} -eq 0 ] && CASES=(rest_state straka warm_bubble igw)

build() {  # dir cuda double
  cmake -S . -B "$1" -G Ninja -DCPPWRF_ENABLE_CUDA="$2" -DCPPWRF_DOUBLE="$3" > /dev/null
  cmake --build "$1" > /dev/null
}
VARIANTS=("build-cpu:OFF:OFF:cpu_fp32" "build-cpu-dp:OFF:ON:cpu_fp64")
if command -v nvcc > /dev/null && nvidia-smi > /dev/null 2>&1; then
  VARIANTS+=("build-cuda:ON:OFF:cuda_fp32" "build-cuda-dp:ON:ON:cuda_fp64")
else
  echo "note: CUDA not usable (nvcc or driver missing); CPU only"
fi

status=0
for v in "${VARIANTS[@]}"; do
  IFS=: read -r dir cuda dbl tag <<< "$v"
  build "$dir" "$cuda" "$dbl"
  for c in "${CASES[@]}"; do
    out="out/${c}_${tag}"
    echo "== $c [$tag]"
    "./$dir/cppwrf" "configs/$c.ini" --set "output.dir=$out" | tail -1 || status=1
    python3 tests/verify/verify.py check "$c" "$out" || status=1
  done
done
exit $status
