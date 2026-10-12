# Regression Test Report

- **Date**: 2026-10-12T23:11:31Z
- **Git Commit**: `9e5586e6dbd1d2011a6cfbe2b934fca9de2ddcb9` (`9e5586e`)
- **Mode**: `quick`
- **GPU Enabled**: No
- **Overall Result**: **FAIL**

## 1. Build and CTest Status

| Configuration | Build | CTest |
|---|---|---|
| `build-reg-cpu` | PASS | PASS |
| `build-reg-cpu-dp` | PASS | PASS |

## 2. Simulation Cases

| Case | Backend | Precision | Run Status | Verifier | Golden Match | ns/cell/step | Perf Delta |
|---|---|---|---|---|---|---|---|
| rest_state | cpu | fp32 | PASS | PASS | PASS | 1149.40 | -84.9% |
| straka | cpu | fp32 | PASS | FAIL | PASS | 739.05 | -73.1% |
| warm_bubble | cpu | fp32 | PASS | PASS | PASS | 773.71 | -83.4% |
| igw | cpu | fp32 | PASS | PASS | PASS | 1139.33 | -93.4% |
| rest_state | cpu | fp64 | PASS | PASS | PASS | 1560.61 | -75.9% |
| straka | cpu | fp64 | PASS | FAIL | PASS | 995.21 | -79.6% |
| warm_bubble | cpu | fp64 | PASS | PASS | PASS | 885.92 | -87.5% |
| igw | cpu | fp64 | PASS | PASS | PASS | 1055.49 | -71.2% |

---
*Report generated automatically by `tools/regression/run.sh`.*
