#!/usr/bin/env bash
# tools/wrf2cppwrf.py unit tests (seconds). The checks against WRF reference runs skip without
# ${CPPWRF_WRF_RUNS} / ${CPPWRF_WRF_REF}; 77 (skipped) without python3.
command -v python3 > /dev/null || { echo "python3 not found"; exit 77; }
exec python3 "$(dirname "$0")/test_wrf2cppwrf.py" -v
