# RRTMG in single precision in CUDA FP32 builds

Measured 2026-11-16 on one laptop GPU (below). What changes: **GPU FP32 only.** GPU FP64, CPU FP64 and CPU FP32 give
the same bits as before (checks below).

## The idea

In the FP32 build (`-DCPPWRF_DOUBLE=OFF`) everything is `float` except RRTMG LW/SW: its device kernels (`KLwPrep`,
`KLwLayer`, `KLwGpoint`, `KSwPrep`, `KSwGpoint`, ...) compute in `double`. On consumer and workstation GPUs FP64 runs
at 1/64 of the FP32 rate, so in GPU FP32 runs with `radt = 1` RRTMG is most of the step.

WRF real*4 computes RRTMG in single precision (`rb = kind(1.0)`), so running it in `float` brings CPPWRF FP32 closer
to the build it is compared with in [FP32_ACCURACY.md](FP32_ACCURACY.md).

## The change

Mechanical, by the scripts in [tools/experiments/rrtmg_fp32](../tools/experiments/rrtmg_fp32): run in this order on a
clean checkout of c9c95c5 they give this change's `src/` exactly.

`rr` (`src/phys/ra_rrtmg_common/rrtmg_real.hpp`) is `cppwrf::real` in CUDA builds (float in FP32, double in FP64) and
`double` in CPU builds: the gain is a GPU one (FP64 at 1/32 to 1/64 of the FP32 rate), on the CPU the float RRTMG
was not faster in any measurement (below), and CPU FP32 results stay as they were.

1. `prep_exp1.py`: every `double` in `src/phys/ra_rrtmg_{lw,sw,common}` becomes `rr`. The binary tables stay
   `double` on disk: they are read into a `double` copy of the structs (`*_tables_f64.hpp`, namespace `f64`) and
   converted member by member.
2. `prep_exp2.py`: the radiation driver and Dudhia SW (they share RRTMG's strided views) too; every floating literal
   in these files becomes `rr(literal)` (an unsuffixed C++ literal is `double` and would pull the arithmetic back to
   FP64); `float` overloads of the Fortran-format readers (read as `double`, rounded once).
3. `prep_exp3.py` to `prep_exp5.py`: the call sites outside RRTMG. The old RRTM LW scheme (`ra_lw_physics = 1`,
   src/phys/ra_rrtm_lw) stays `double`; in CUDA FP32 builds the driver widens the columns it hands over and rounds
   the results once. `rrtmg_real.hpp` gets the CPU / CUDA choice above.

Where `rr` is `double`, `rr(x)` is `x` with its type, so the code is the same as before.

## Setup

RTX 4060 Laptop GPU (sm_89, 8 GB, driver 576.88), 16 host threads, WSL2 Ubuntu 26.04, CUDA 12.9 (nvcc 12.9.86) with
g++-14 as host compiler. glibc 2.43 needs a header workaround for CUDA 12.9 (`local_cuda.sh`, user-local, no change
to the system install). Case N3 (`configs/real_n3.ini`: P4 physics, two-way 5 km nest, `radt = 1`) from
`n3_data.tar`. Upstream = c9c95c5 unchanged; experiment = this branch. Both built with the same flags.

## Checks that the other builds are unchanged

| Check | Result |
|---|---|
| CPU FP64, N3 3 d01 steps, upstream vs this change | all snapshot files byte-identical |
| CPU FP32, N3 3 d01 steps | all snapshot files byte-identical |
| CPU FP32, N3 3 d01 steps, `ra_lw_physics = 1` | all snapshot files byte-identical |
| GPU FP64, N3 3 d01 steps, with RRTMG and with `ra_lw_physics = 1` | all snapshot files byte-identical |
| GPU FP64, N3 1 h, upstream vs this change (wrfout d01, d02) | max \|diff\| = 0 in every variable and frame |
| `examples/selftest/run.sh` cpu-fp64 and gpu-fp64 | PASS (the reference hash); gpu-fp32: PASS (finite) |
| GPU FP64 (upstream) vs WRF real*8 reference `run_N3h_r8` at 1 h | 0 in every variable of the reference file (RAIN excluded: the reference has RAINNC only) |
| GPU FP64 speed, N3 1 h | 43.66 s vs 43.86 s wall (same) |

## FP64 instructions left in RRTMG (cuobjdump -sass, FP32 build)

`sass_f64.sh`: SASS instructions DADD / DMUL / DFMA / DSETP / F2F.F64 / DMNMX / MUFU.*64H in the RRTMG batch kernels.

| | Upstream FP32 | Experiment FP32 |
|---|---|---|
| RRTMG kernels | 12956 | **0** |
| whole binary | 165687 | 151045 |

The rest of the FP32 binary still has FP64 instructions outside RRTMG (not examined here).

## Speed (GPU FP32, RTX 4060 Laptop)

Stepping time (`cppwrf: done ... stepping`, d01 includes its nest) and wall time of the run:

| Run | c9c95c5 | This change | Speed-up |
|---|---|---|---|
| N3 1 h, radt 1 (night), d01 stepping | 19.89 / 19.94 s | 11.61 / 11.60 s | **1.71x** |
| same, d02 stepping | 15.30 / 15.42 s | 7.89 / 7.91 s | **1.95x** |
| same, wall | 23.20 / 23.04 s | 14.71 / 14.69 s | **1.57x** |
| N3 1 h, radt 15, wall | 11.33 s | 10.67 s | 1.06x |
| N3 6 h, radt 1 (00-06 UTC, sunrise at about 04 UTC), d01 stepping | 209.1 s | 88.3 s | **2.37x** |
| same, d02 stepping | 188.2 s | 70.2 s | **2.68x** |
| same, wall | 212.5 s | 92.0 s | **2.31x** |
| N3 6 h, built with `--fmad=false`, wall | 265.1 s | 116.4 s | **2.28x** |

(two values: two runs, alternating.) The 6 h runs gain more than the 1 h run: after sunrise the SW g-point kernels
run too. For reference, GPU FP64 takes 322.3 s wall for the 6 h.

**CPU**, with `rr` = float also on the CPU (an earlier version of this change), CPU FP32 N3 with 16 threads: 1 h
391.7 s (c9c95c5) -> 459.7 s; 20 min, alternating, 206.5 / 184.3 s and 181.1 / 158.9 s (float first). This laptop's
CPU times are noisy (identical binaries differed by 20 % and more between consecutive runs), so these numbers only
say that float was not faster; CPU builds keep `double`.

## Accuracy (GPU FP32 against GPU FP64 = WRF real*8)

No WRF real*4 run was available here (the study's yardstick), so the comparison is c9c95c5 FP32 vs this change's
FP32, both against the GPU FP64 run of the same 6 h. Two FP32 realizations of each: the release build (CUDA FMA
contraction) and one built with `--fmad=false`, to see how much the RMSE moves between two equally valid FP32
answers. `tools/verify/fp32_accuracy.py compare`, `compare_6h.sh`.

Ratio R = mean RMSE of the two experiment runs / mean RMSE of the two upstream runs (R <= 1: the experiment is at
least as close to FP64):

| Domain, time | median R | range |
|---|---|---|
| d01 3 h | 0.91 | 0.81 (W) to 1.13 (QICE) |
| d01 6 h | 0.98 | 0.92 to 1.66 (QGRAUP; QSNOW 1.41, QRAIN 1.18) |
| d02 3 h | 1.08 | 0.79 to 1.44 (PBLH) |
| d02 6 h | 1.00 | 0.78 to 1.11 |

d01 at 6 h, RMSE against FP64 (upstream FMA, upstream no-FMA | experiment FMA, experiment no-FMA | R | upstream
FMA vs upstream no-FMA):

| Var | upstream | experiment | R | upstream vs upstream |
|---|---|---|---|---|
| T2 | 0.0359, 0.0318 | 0.0317, 0.034 | 0.97 | 0.0358 |
| U10 | 0.0733, 0.0788 | 0.0668, 0.0737 | 0.92 | 0.0664 |
| PSFC | 1.69, 1.83 | 1.73, 1.83 | 1.01 | 1.67 |
| TSK | 0.148, 0.143 | 0.144, 0.149 | 1.01 | 0.142 |
| SWDOWN | 9.51, 9.69 | 9.96, 9.92 | 1.03 | 9.2 |
| HFX | 2.72, 2.69 | 2.73, 2.8 | 1.02 | 2.7 |
| RAIN | 0.232, 0.219 | 0.189, 0.236 | 0.94 | 0.213 |
| T | 0.0423, 0.0411 | 0.0364, 0.0432 | 0.95 | 0.0393 |
| QVAPOR | 2.25e-05, 2.28e-05 | 2.12e-05, 2.28e-05 | 0.97 | 2.06e-05 |
| U | 0.0661, 0.067 | 0.0594, 0.0694 | 0.97 | 0.062 |
| W | 0.00865, 0.0087 | 0.00843, 0.00872 | 0.99 | 0.00827 |
| QCLOUD | 1.96e-06, 1.71e-06 | 1.55e-06, 1.81e-06 | 0.92 | 1.64e-06 |
| QRAIN | 6.47e-07, 7.21e-07 | 8.47e-07, 7.73e-07 | 1.18 | 6.53e-07 |
| QSNOW | 4.43e-07, 4.52e-07 | 6.26e-07, 6.34e-07 | 1.41 | 5.04e-07 |
| QGRAUP | 1.33e-07, 2.17e-07 | 2.97e-07, 2.84e-07 | 1.66 | 1.84e-07 |

Reading:
- The surface fields, the radiation (SWDOWN) and the dynamics are as close to FP64 as upstream FP32 is.
- d01 snow and graupel at 6 h are 1.4x / 1.7x in both realizations. Their bias is 1e-9 to 1e-10 with mixed signs
  (no systematic shift); the RMSE comes from a few points (max |diff| 2e-4); the two upstream FP32 builds differ
  from each other by 5e-7 (QSNOW), the experiment from FP64 by 6.3e-7. Not settled by two realizations of one 6 h
  run; d02 shows no such signal (R 0.95 / 1.00).

## Repository tests

`ctest -L unit` on the CPU FP64 build: 10 pass, 2 skipped, 3 fail (`python_package`, `post_output`, `verify_obs`)
in the Python test fixtures with numpy 2.3 (`'numpy.bytes_' object has no attribute 'encode'` in
`netCDF4.stringtochar`). The same three fail on c9c95c5 in the same environment: not this change.
`tools/gen_option_docs.py --check`: current. Not run here: the WRF identity suite (`tests/verify/wrf_suite.py`, no
reference data here) and the WRF differential tests (no WRF build).

## Not done

- The study's yardstick: WRF real*4 against real*8 for the same cases, the eleven 1 h cases and the 24 h / 72 h
  runs (tools/bench/colab_fp32.sh), and more realizations for the d01 snow / graupel signal.
- Other GPUs (the README's RTX PRO 6000 also runs FP64 at 1/64; datacenter GPUs with fast FP64 gain less).
- The table conversion is generated code (`convert_tables`); a check that keeps FP64 out of the RRTMG kernels
  (`sass_f64.sh`) could run in CI.
