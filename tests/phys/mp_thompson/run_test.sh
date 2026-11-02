#!/usr/bin/env bash
# Thompson (mp_physics = 8 and the aerosol-aware 28) tests: WRF v4.8.0's thompson_init tables and mp_gt_driver (module_mp_thompson.o of the
# -r8 -fdefault-double-8 build, WRF_R8D8_DIR) against the C++ tables (test_tables) and the C++ Thompson
# (driver_thompson_cpu of the CMake build), tolerance 0.
#   run_test.sh BUILD_DIR      (BUILD_DIR: the directory of driver_thompson_cpu; the work files go there)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="$(cd "${1:?usage: run_test.sh BUILD_DIR}" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
. "${SCRIPT_DIR}/../wrf_fortran.sh"
wrf_need r8d8
cd "${WORK_DIR}"

FORT_BIN="${WORK_DIR}/fortran_driver_thompson"
CPP_BIN="${WORK_DIR}/driver_thompson_cpu"
INPUT_BIN="${WORK_DIR}/input_thompson.bin"
FORT_OUT="${WORK_DIR}/output_fortran.bin"
CPP_OUT="${WORK_DIR}/output_cpp.bin"
wrf_fortran "${FORT_BIN}" "${SCRIPT_DIR}/fortran/driver_thompson.f90"
wrf_fortran "${WORK_DIR}/dump_tables_thompson" "${SCRIPT_DIR}/dump_tables_f.f90"

# 1. the tables: thompson_init computes them (no table files here, write_thompson_tables unset) and the C++
#    init_thompson_tables computes its own
echo "=== Thompson lookup tables (thompson_init) ==="
rm -f qr_acr_qg_V4.dat qr_acr_qsV2.dat freezeH2O.dat tables_fortran.bin
"${REPO_ROOT}/tools/heavy.sh" "${WORK_DIR}/dump_tables_thompson"  # writes tables_fortran.bin
"${REPO_ROOT}/tools/heavy.sh" "${WORK_DIR}/test_tables" "${WORK_DIR}/tables_fortran.bin"

# 2. the columns: the Fortran driver's thompson_init writes WRF's table files on its first run, which the C++
#    driver reads as the model reads a run's tables; 201 x 41: the odd tails of WRF's vectorised loops
for cfg in "2500 50" "201 41"; do
  set -- $cfg
  echo "=== Thompson Microphysics Differential Test ($1 cols x $2 levels) ==="
  python3 "${SCRIPT_DIR}/generate_inputs.py" --ncols "$1" --nlevels "$2"
  "${REPO_ROOT}/tools/heavy.sh" "${FORT_BIN}" "${INPUT_BIN}" "${FORT_OUT}"
  "${REPO_ROOT}/tools/heavy.sh" "${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}"
  python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPP_OUT}" --tol 0
  if [ "$1" = 2500 ] && [ -x "${WORK_DIR}/driver_thompson_cuda_dp" ]; then
    "${REPO_ROOT}/tools/heavy.sh" "${WORK_DIR}/driver_thompson_cuda_dp" "${INPUT_BIN}" "${WORK_DIR}/output_cpp_cuda_dp.bin"
    python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${WORK_DIR}/output_cpp_cuda_dp.bin" --tol 1e-12
  fi
done
# 3. mp_physics = 28 (aerosol aware, aer_init_opt 0, wif_input_opt 1): WRF's mp_gt_driver with nc, nwfa, nifa,
#    nwfa2d, nifa2d (fortran/driver_thompson_aero.f90) against driver_thompson_cpu --aero; CCN_ACTIVATE.BIN of W's run/
[ -f "$W/run/CCN_ACTIVATE.BIN" ] || { echo "no $W/run/CCN_ACTIVATE.BIN"; exit 1; }
ln -sf "$W/run/CCN_ACTIVATE.BIN" CCN_ACTIVATE.BIN
wrf_fortran "${WORK_DIR}/fortran_driver_thompson_aero" "${SCRIPT_DIR}/fortran/driver_thompson_aero.f90"
for cfg in "2500 50" "201 41"; do
  set -- $cfg
  echo "=== Aerosol-aware Thompson Differential Test ($1 cols x $2 levels) ==="
  python3 "${SCRIPT_DIR}/generate_inputs.py" --ncols "$1" --nlevels "$2" --aero
  "${REPO_ROOT}/tools/heavy.sh" "${WORK_DIR}/fortran_driver_thompson_aero" "${INPUT_BIN}" "${FORT_OUT}"
  "${REPO_ROOT}/tools/heavy.sh" "${CPP_BIN}" "${INPUT_BIN}" "${CPP_OUT}" --aero
  python3 "${SCRIPT_DIR}/compare.py" "${FORT_OUT}" "${CPP_OUT}" --tol 0 --aero
done
ls qr_acr_qg_V4.dat qr_acr_qsV2.dat freezeH2O.dat > /dev/null  # the C++ driver read WRF's table files
echo "=== Thompson Microphysics Differential Test PASSED ==="
