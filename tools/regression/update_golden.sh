#!/usr/bin/env bash
# CPPWRF Golden Reference Updater
# Usage: tools/regression/update_golden.sh [--quick|--full] [--gpu]
# Refreshes stored golden outputs under tools/regression/golden/ only on explicit request.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

MODE="quick"
ENABLE_GPU=0

for arg in "$@"; do
  case "$arg" in
    --quick) MODE="quick" ;;
    --full)  MODE="full" ;;
    --gpu)   ENABLE_GPU=1 ;;
    *) echo "Unknown option: $arg"; exit 2 ;;
  esac
done

echo "Refreshing golden references from HEAD..."
# Execute run.sh to produce fresh outputs
"${SCRIPT_DIR}/run.sh" "--${MODE}" $([ ${ENABLE_GPU} -eq 1 ] && echo "--gpu") || true

GOLDEN_DIR="${SCRIPT_DIR}/golden"
RUNS_DIR="${SCRIPT_DIR}/runs"
mkdir -p "${GOLDEN_DIR}"

GIT_COMMIT="$(git rev-parse HEAD 2>/dev/null || echo "unknown")"
DATE_STR="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

echo "${GIT_COMMIT}" > "${GOLDEN_DIR}/COMMIT"
echo "${DATE_STR}" > "${GOLDEN_DIR}/DATE"
echo "${MODE}" > "${GOLDEN_DIR}/MODE"

COUNT=0
for out in "${RUNS_DIR}"/*; do
  if [ -d "${out}" ]; then
    case_key="$(basename "${out}")"
    tgt_gold="${GOLDEN_DIR}/${case_key}"
    rm -rf "${tgt_gold}"
    mkdir -p "${tgt_gold}"

    if [ -f "${out}/diag.csv" ]; then
      cp "${out}/diag.csv" "${tgt_gold}/"
    fi

    # Find last snapshot
    last_snap_json=$(ls "${out}"/snap_*.json 2>/dev/null | sort | tail -1 || true)
    if [ -n "${last_snap_json}" ] && [ -f "${last_snap_json}" ]; then
      last_snap_bin="${last_snap_json%.json}.bin"
      cp "${last_snap_json}" "${tgt_gold}/"
      if [ -f "${last_snap_bin}" ]; then
        cp "${last_snap_bin}" "${tgt_gold}/"
      fi
    fi
    echo "Updated golden for ${case_key}"
    COUNT=$((COUNT + 1))
  fi
done

echo "Successfully refreshed ${COUNT} golden references at commit ${GIT_COMMIT}."
