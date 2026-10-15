#!/usr/bin/env bash
# CPPWRF regression runner
# Usage: tools/regression/run.sh [--quick|--full] [--gpu]

set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

# Modest threading so that concurrent jobs on the machine are not starved
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
export OMP_WAIT_POLICY="${OMP_WAIT_POLICY:-passive}"

MODE="quick"
ENABLE_GPU=-1

for arg in "$@"; do
  case "$arg" in
    --quick) MODE="quick" ;;
    --full)  MODE="full" ;;
    --gpu)   ENABLE_GPU=1 ;;
    --no-gpu) ENABLE_GPU=0 ;;
    *) echo "Unknown option: $arg"; echo "Usage: $0 [--quick|--full] [--gpu|--no-gpu]"; exit 2 ;;
  esac
done

LOG_DIR="${SCRIPT_DIR}/logs"
RUNS_DIR="${SCRIPT_DIR}/runs"
GOLDEN_DIR="${SCRIPT_DIR}/golden"
mkdir -p "${LOG_DIR}" "${RUNS_DIR}" "${GOLDEN_DIR}"

DATE_STR="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
GIT_COMMIT="$(git rev-parse HEAD 2>/dev/null || echo "unknown")"
GIT_SHORT="$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")"

echo "=========================================================="
echo "CPPWRF Regression Suite"
echo "Commit: ${GIT_COMMIT}"
echo "Date:   ${DATE_STR}"
echo "Mode:   ${MODE}"
echo "OMP:    OMP_NUM_THREADS=${OMP_NUM_THREADS} OMP_WAIT_POLICY=${OMP_WAIT_POLICY}"
echo "Logs:   ${LOG_DIR}"
echo "Runs:   ${RUNS_DIR}"
echo "=========================================================="

OVERALL_PASS=1
declare -A BUILD_STATUS
declare -A CTEST_STATUS
declare -A CASE_STATUS
declare -A CASE_PERF
declare -A CASE_PERF_DELTA
declare -A GOLDEN_STATUS
declare -A VERIFY_STATUS
declare -A GPU_COMP_STATUS

CONFIGS=("rest_state" "straka" "warm_bubble" "igw" "squall2d")

# Check GPU availability
HAS_GPU=0
if [ ${ENABLE_GPU} -ne 0 ]; then
  if command -v nvcc >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
    HAS_GPU=1
  else
    if [ ${ENABLE_GPU} -eq 1 ]; then
      echo "Notice: GPU requested but nvcc or GPU driver not functional; GPU builds/runs will be skipped."
    fi
  fi
fi
echo "GPU Status: $([ ${HAS_GPU} -eq 1 ] && echo "Active (CUDA enabled)" || echo "Disabled/Not found")"

# Define build targets: name:enable_cuda:enable_double:backend:prec
TARGETS=(
  "build-reg-cpu:OFF:OFF:cpu:fp32"
  "build-reg-cpu-dp:OFF:ON:cpu:fp64"
)
if [ ${HAS_GPU} -eq 1 ]; then
  TARGETS+=(
    "build-reg-cuda:ON:OFF:cuda:fp32"
    "build-reg-cuda-dp:ON:ON:cuda:fp64"
  )
fi

echo ""
echo "--- 1. Building configurations & running CTest (wrapped in tools/heavy.sh) ---"
for tgt in "${TARGETS[@]}"; do
  IFS=: read -r bdir cuda dbl backend prec <<< "$tgt"
  echo ">> Building ${bdir} (CUDA=${cuda}, DOUBLE=${dbl})..."
  "${REPO_ROOT}/tools/heavy.sh" cmake -S . -B "${bdir}" -G Ninja -DCPPWRF_ENABLE_CUDA="${cuda}" -DCPPWRF_DOUBLE="${dbl}" > "${LOG_DIR}/${bdir}.cmake.log" 2>&1
  if "${REPO_ROOT}/tools/heavy.sh" cmake --build "${bdir}" -j4 > "${LOG_DIR}/${bdir}.build.log" 2>&1; then
    BUILD_STATUS["${bdir}"]="PASS"
    echo "   Build: PASS"
  else
    BUILD_STATUS["${bdir}"]="FAIL"
    OVERALL_PASS=0
    echo "   Build: FAIL (see ${LOG_DIR}/${bdir}.build.log)"
    continue
  fi

  echo ">> Running ctest in ${bdir}..."
  if "${REPO_ROOT}/tools/heavy.sh" ctest --test-dir "${bdir}" --timeout 120 --output-on-failure > "${LOG_DIR}/${bdir}.ctest.log" 2>&1; then
    CTEST_STATUS["${bdir}"]="PASS"
    echo "   CTest: PASS"
  else
    CTEST_STATUS["${bdir}"]="FAIL"
    OVERALL_PASS=0
    echo "   CTest: FAIL (see ${LOG_DIR}/${bdir}.ctest.log)"
  fi
done

echo ""
echo "--- 2. Running Simulation Cases (${MODE} mode, exclusive timing) ---"
PERF_HIST="${SCRIPT_DIR}/perf_history.csv"
if [ ! -f "${PERF_HIST}" ]; then
  echo "date,git_commit,case,backend,precision,ns_per_cell_step" > "${PERF_HIST}"
fi

for tgt in "${TARGETS[@]}"; do
  IFS=: read -r bdir cuda dbl backend prec <<< "$tgt"
  if [ "${BUILD_STATUS[${bdir}]:-FAIL}" != "PASS" ]; then
    echo "Skipping runs for ${bdir} due to build failure."
    continue
  fi

  bin_path="${bdir}/cppwrf"
  if [ ! -x "${bin_path}" ]; then
    echo "Executable ${bin_path} not found; skipping."
    continue
  fi

  for c in "${CONFIGS[@]}"; do
    case_key="${c}_${backend}_${prec}"
    out_dir="${RUNS_DIR}/${case_key}"
    rm -rf "${out_dir}"
    mkdir -p "${out_dir}"

    EXTRA_ARGS=()
    if [ "${MODE}" == "quick" ]; then
      case "$c" in
        rest_state)  EXTRA_ARGS=(--set run.run_time=20 --set output.out_interval=2 --set output.diag_interval=1) ;;
        straka)      EXTRA_ARGS=(--set run.run_time=2 --set output.out_interval=2 --set output.diag_interval=1) ;;
        warm_bubble) EXTRA_ARGS=(--set run.run_time=2 --set output.out_interval=2 --set output.diag_interval=1) ;;
        igw)         EXTRA_ARGS=(--set run.run_time=12 --set output.out_interval=2 --set output.diag_interval=1) ;;
        squall2d)    EXTRA_ARGS=(--set run.run_time=12 --set output.out_interval=2 --set output.diag_interval=1) ;;
      esac
    fi

    echo ">> Running ${c} [${backend} ${prec}] (tools/heavy.sh --exclusive)..."
    run_log="${LOG_DIR}/${case_key}.run.log"
    if "${REPO_ROOT}/tools/heavy.sh" --exclusive "${bin_path}" "configs/${c}.ini" --set "output.dir=${out_dir}" "${EXTRA_ARGS[@]}" > "${run_log}" 2>&1; then
      CASE_STATUS["${case_key}"]="PASS"
      echo "   Execution: PASS"
    else
      CASE_STATUS["${case_key}"]="FAIL (crash)"
      OVERALL_PASS=0
      echo "   Execution: FAIL (crashed, see ${run_log})"
      continue
    fi

    # Parse ns/cell/step
    perf_val=$(grep -oE "[0-9]+(\.[0-9]+)? ns/cell/step" "${run_log}" | head -1 | awk '{print $1}' || echo "")
    if [ -n "${perf_val}" ]; then
      CASE_PERF["${case_key}"]="${perf_val}"
      # Find previous entry in perf_history.csv
      prev_val=$(awk -F, -v c="${c}" -v b="${backend}" -v p="${prec}" '$3==c && $4==b && $5==p {last=$6} END {if(last) print last}' "${PERF_HIST}")
      if [ -n "${prev_val}" ]; then
        delta=$(python3 -c "new=float('${perf_val}'); old=float('${prev_val}'); diff=(new-old)/old*100; print(f'{diff:+.1f}%')")
        is_slow=$(python3 -c "new=float('${perf_val}'); old=float('${prev_val}'); print(1 if new > old*1.10 else 0)")
        if [ "${is_slow}" -eq 1 ]; then
          CASE_PERF_DELTA["${case_key}"]="${delta} (SLOWDOWN >10%)"
          echo "   Performance: ${perf_val} ns/cell/step (${delta} SLOWDOWN)"
        else
          CASE_PERF_DELTA["${case_key}"]="${delta}"
          echo "   Performance: ${perf_val} ns/cell/step (${delta})"
        fi
      else
        CASE_PERF_DELTA["${case_key}"]="baseline"
        echo "   Performance: ${perf_val} ns/cell/step (new baseline)"
      fi
      echo "${DATE_STR},${GIT_COMMIT},${c},${backend},${prec},${perf_val}" >> "${PERF_HIST}"
    else
      CASE_PERF["${case_key}"]="N/A"
      CASE_PERF_DELTA["${case_key}"]="N/A"
    fi

    # Run tests/verify/verify.py check
    if [ "${MODE}" == "quick" ] && [ "$c" == "straka" ]; then
      VERIFY_STATUS["${case_key}"]="SKIP (quick: 2s vs 900s)"
      echo "   Verifier: SKIP (quick mode: straka verification requires 900s)"
    else
      echo "   Running verify.py check..."
      if python3 tests/verify/verify.py check "${c}" "${out_dir}" > "${LOG_DIR}/${case_key}.verify.log" 2>&1; then
        VERIFY_STATUS["${case_key}"]="PASS"
        echo "   Verifier: PASS"
      else
        VERIFY_STATUS["${case_key}"]="FAIL"
        OVERALL_PASS=0
        echo "   Verifier: FAIL (see ${LOG_DIR}/${case_key}.verify.log)"
      fi
    fi

    # Golden comparison
    gold_dir="${GOLDEN_DIR}/${case_key}"
    echo "   Comparing against golden output (${gold_dir})..."
    if [ -d "${gold_dir}" ]; then
      if python3 tools/regression/compare_golden.py "${out_dir}" "${gold_dir}" "${prec}" > "${LOG_DIR}/${case_key}.golden.log" 2>&1; then
        GOLDEN_STATUS["${case_key}"]="PASS"
        echo "   Golden: PASS"
      else
        GOLDEN_STATUS["${case_key}"]="FAIL"
        OVERALL_PASS=0
        echo "   Golden: FAIL (see ${LOG_DIR}/${case_key}.golden.log)"
      fi
    else
      GOLDEN_STATUS["${case_key}"]="SKIP (no golden)"
      echo "   Golden: SKIP (golden dir not found, run update_golden.sh)"
    fi
  done
done

# CPU vs GPU Comparisons
echo ""
echo "--- 3. CPU vs GPU Numerical Comparisons ---"
if [ ${HAS_GPU} -eq 1 ]; then
  for prec in "fp32" "fp64"; do
    for c in "${CONFIGS[@]}"; do
      cpu_dir="${RUNS_DIR}/${c}_cpu_${prec}"
      gpu_dir="${RUNS_DIR}/${c}_cuda_${prec}"
      comp_key="${c}_${prec}"
      if [ -d "${cpu_dir}" ] && [ -d "${gpu_dir}" ]; then
        rtol="1e-3"
        if python3 tests/verify/verify.py compare "${cpu_dir}" "${gpu_dir}" --rtol "${rtol}" > "${LOG_DIR}/compare_${comp_key}.log" 2>&1; then
          GPU_COMP_STATUS["${comp_key}"]="PASS"
          echo ">> CPU vs GPU ${c} [${prec}]: PASS"
        else
          GPU_COMP_STATUS["${comp_key}"]="FAIL"
          OVERALL_PASS=0
          echo ">> CPU vs GPU ${c} [${prec}]: FAIL (see ${LOG_DIR}/compare_${comp_key}.log)"
        fi
      else
        GPU_COMP_STATUS["${comp_key}"]="SKIP"
        echo ">> CPU vs GPU ${c} [${prec}]: SKIP"
      fi
    done
  done
fi

# WRF Reference Fidelity (straka_wrf on CPU FP64)
echo ""
echo "--- 4. WRF Reference Fidelity (straka_wrf on CPU FP64) ---"
WRF_REF="${HOME}/wrf_ref/run_straka_d1/wrfout_d01_0001-01-01_00:00:00"
STRAKA_WRF_STATUS="SKIP"
STRAKA_WRF_DETAILS="reference file missing"
if [ -f "${WRF_REF}" ]; then
  bin_path_dp="build-reg-cpu-dp/cppwrf"
  if [ -x "${bin_path_dp}" ]; then
    out_straka_wrf="${RUNS_DIR}/straka_wrf_cpu_fp64"
    rm -rf "${out_straka_wrf}"
    mkdir -p "${out_straka_wrf}"
    echo ">> Running straka_wrf [cpu fp64] (tools/heavy.sh --exclusive)..."
    run_log="${LOG_DIR}/straka_wrf_cpu_fp64.run.log"
    if "${REPO_ROOT}/tools/heavy.sh" --exclusive "${bin_path_dp}" "configs/straka_wrf.ini" --set "output.dir=${out_straka_wrf}" > "${run_log}" 2>&1; then
      echo "   Execution: PASS"
      echo "   Running verify.py wrf..."
      if python3 tests/verify/verify.py wrf "${WRF_REF}" "${out_straka_wrf}" > "${LOG_DIR}/straka_wrf.verify.log" 2>&1; then
        STRAKA_WRF_STATUS="PASS"
        STRAKA_WRF_DETAILS="all steps match within RMS tol 0.05 K"
        echo ">> WRF Reference Comparison: PASS"
      else
        STRAKA_WRF_STATUS="FAIL"
        STRAKA_WRF_DETAILS="verify.py wrf reported differences exceeding tolerance"
        OVERALL_PASS=0
        echo ">> WRF Reference Comparison: FAIL (see ${LOG_DIR}/straka_wrf.verify.log)"
      fi
    else
      STRAKA_WRF_STATUS="FAIL (crash)"
      STRAKA_WRF_DETAILS="simulation run crashed"
      OVERALL_PASS=0
      echo ">> WRF Reference Run: FAIL (see ${run_log})"
    fi
  else
    STRAKA_WRF_STATUS="SKIP (no binary)"
    STRAKA_WRF_DETAILS="build-reg-cpu-dp/cppwrf not found"
    echo ">> WRF Reference Comparison: SKIP (build-reg-cpu-dp/cppwrf not found)"
  fi
else
  STRAKA_WRF_STATUS="SKIP (no ref)"
  STRAKA_WRF_DETAILS="${WRF_REF} not found"
  echo ">> WRF Reference Comparison: SKIP (${WRF_REF} not found)"
fi

echo ""
echo "--- 5. Writing Report (tools/regression/last_report.md) ---"
REPORT="${SCRIPT_DIR}/last_report.md"
{
  echo "# Regression Test Report"
  echo ""
  echo "- **Date**: ${DATE_STR}"
  echo "- **Git Commit**: \`${GIT_COMMIT}\` (\`${GIT_SHORT}\`)"
  echo "- **Mode**: \`${MODE}\`"
  echo "- **GPU Enabled**: $([ ${HAS_GPU} -eq 1 ] && echo "Yes" || echo "No")"
  echo "- **Overall Result**: **$([ ${OVERALL_PASS} -eq 1 ] && echo "PASS" || echo "FAIL")**"
  echo ""
  echo "## 1. Build and CTest Status"
  echo ""
  echo "| Configuration | Build | CTest |"
  echo "|---|---|---|"
  for tgt in "${TARGETS[@]}"; do
    IFS=: read -r bdir cuda dbl backend prec <<< "$tgt"
    echo "| \`${bdir}\` | ${BUILD_STATUS[${bdir}]:-N/A} | ${CTEST_STATUS[${bdir}]:-N/A} |"
  done
  if [ ${ENABLE_GPU} -eq 1 ] && [ ${HAS_GPU} -eq 0 ]; then
    echo "| \`build-reg-cuda*\` | SKIPPED | SKIPPED (GPU driver or nvcc missing) |"
  fi
  echo ""
  echo "## 2. Simulation Cases"
  echo ""
  echo "| Case | Backend | Precision | Run Status | Verifier | Golden Match | ns/cell/step | Perf Delta |"
  echo "|---|---|---|---|---|---|---|---|"
  for tgt in "${TARGETS[@]}"; do
    IFS=: read -r bdir cuda dbl backend prec <<< "$tgt"
    for c in "${CONFIGS[@]}"; do
      case_key="${c}_${backend}_${prec}"
      echo "| ${c} | ${backend} | ${prec} | ${CASE_STATUS[${case_key}]:-N/A} | ${VERIFY_STATUS[${case_key}]:-N/A} | ${GOLDEN_STATUS[${case_key}]:-N/A} | ${CASE_PERF[${case_key}]:-N/A} | ${CASE_PERF_DELTA[${case_key}]:-N/A} |"
    done
  done
  if [ ${ENABLE_GPU} -eq 1 ] && [ ${HAS_GPU} -eq 0 ]; then
    for c in "${CONFIGS[@]}"; do
      echo "| ${c} | cuda | fp32/fp64 | SKIPPED | SKIPPED | SKIPPED | N/A | N/A |"
    done
  fi

  if [ ${HAS_GPU} -eq 1 ]; then
    echo ""
    echo "## 3. CPU vs GPU Agreement"
    echo ""
    echo "| Case | Precision | Status |"
    echo "|---|---|---|"
    for prec in "fp32" "fp64"; do
      for c in "${CONFIGS[@]}"; do
        comp_key="${c}_${prec}"
        echo "| ${c} | ${prec} | ${GPU_COMP_STATUS[${comp_key}]:-N/A} |"
      done
    done
  fi

  echo ""
  echo "## 4. WRF Reference Fidelity (straka_wrf on CPU FP64)"
  echo ""
  echo "| Case | Configuration | Reference File | Status | Details |"
  echo "|---|---|---|---|---|"
  echo "| \`straka_wrf\` | \`cpu fp64\` | \`${WRF_REF}\` | ${STRAKA_WRF_STATUS} | ${STRAKA_WRF_DETAILS} |"

  echo ""
  echo "---"
  echo "*Report generated automatically by \`tools/regression/run.sh\`.*"
} > "${REPORT}"

echo "Report written to ${REPORT}."
if [ ${OVERALL_PASS} -eq 1 ]; then
  echo "SUCCESS: Regression suite PASSED."
  exit 0
else
  echo "FAILURE: One or more checks FAILED."
  exit 1
fi
