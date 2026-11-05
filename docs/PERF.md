# CPPWRF memory-traffic budget

Model: each 3D field touched by a kernel is read (or written) once from DRAM per kernel,
i.e. stencil neighbours are assumed to hit cache; 1D base-state arrays and 2D fields are
negligible (O(1/nz)). Bytes per cell for FP32 (x2 for FP64). These are **counts from the
code, not measurements**; measured times are in docs/STATUS.md.

## Large step, per RK stage
| Kernel | reads (3D) | writes (3D) | B/cell FP32 |
|---|---|---|---|
| halo fills (8 prognostics/diagnostics + 6 derived) | halo only | halo only | ~0 (O(6/nx)) |
| KMu (mut, muu, muv) | 2D | 2D | ~0 |
| KPrep3 (ru, rv, rw, alt, php) | u v w al ph | ru rv rw alt php | 40 |
| KWW (calc_ww_cp, column) | ru rv (x2 passes) | ww | 12–20 |
| KTendU (advect_u + PGF [+tendf]) | u ru rv ww p ph alt al php [tendf] | ru_tend | 40–44 |
| KTendV (3D only) | v ru rv ww p ph alt al php [tendf] | rv_tend | 40–44 |
| KTendW (advect_w + pg_buoy_w) | w ru rv ww p [tendf] | rw_tend | 24–28 |
| KTendT (advect_scalar theta) | t ru rv ww [tendf] | t_tend | 20–24 |
| KTendPh (rhs_ph) | u v w ww ph | ph_tend | 24 |

Per stage: x-z (no KTendV): ~176 B/cell; 3D: ~220 B/cell.

## Once per time step
| Item | B/cell FP32 |
|---|---|
| copy *_2 -> *_1 (u v w t ph; mu is 2D) | 40 |
| diffusion, stage 1 only (KDiffU/V/W/T: field + alt -> tendf) | 36 (x-z) / 48 (3D) |

## Total large-step traffic per full time step
x-z: 3 x 176 + 40 + 36 ≈ **604 B/cell** (FP32); 3D: 3 x 220 + 40 + 48 ≈ **748 B/cell**.
At 336 GB/s this is a lower bound of ~1.8 ns/cell/step (x-z) for the large step alone on the
RTX 2060; the acoustic substeps (n_sound/2 + n_sound + 1 substeps) come on top.

## Known inefficiencies
- KTendU/V/W/T recompute both face fluxes per cell (2x flux arithmetic, no scratch arrays).
- Separate halo-fill launches per field (14 per stage x up to 2 launches); since F04 adjacent fills share one batch.
- `copy_dev` on CPU is a single-threaded memcpy.

## Profile (measured, Colab, commit 0577fb8)

Measured with `tools/colab/colab_run.sh --cases '' --bench ... --env CPPWRF_KTIME=1 --nsys --ncu top` in one session
per GPU: RTX PRO 6000 Blackwell Server ("G4", 188 SMs, sm_120, FP64 at 1/64 of FP32), L4 (58 SMs, sm_89, 1/64) and
A100-SXM4-40GB (108 SMs, sm_80, FP64 at 1/2 of FP32; P1 and P4 only: Colab assigned no H100 in 7 tries over an hour,
"Service Unavailable"; the A100 has the same half-rate FP64). CUDA 13.0, driver 580.82. Real-data domain
119 x 89 x 44 (N3: that domain and its two-way nest d02 with parent_time_step_ratio 3), dt 60 s, 1 h, no history
output. How the numbers are obtained:

- *stepping s/h*, *ms/step*: the bench run without any profiler (`cppwrf: done` line of the first domain; a nest
  parent's stepping time contains its nests' steps). *stepping with KTIME*: the same run with `CPPWRF_KTIME=1`
  (a device synchronisation after every launch: +0.15 to +1.1 s per hour).
- *per step of the time loop* (nsys): the difference of an nsys run of 1 h and one of 30 min, divided by the 30
  steps between them (initialisation cancels). *GPU busy* = kernel time / wall time under nsys.
- *host waits*: cudaStreamSynchronize (the WRF fatal-stop check `fatal_check`, a 4-byte copy from the device flag),
  cudaDeviceSynchronize (`dev_sync` after each step) and blocking cudaMemcpy.
- *device memory*: the GPU's peak used memory during the bench run (nvidia-smi every 0.2 s; idle 0 MiB), i.e.
  cudaMemGetInfo's total - free of the only process, including the local-memory reserve of the driver.
- kernel tables (G4): the 15 kernels with the most time per step (nsys), with ncu of their launches in the first
  2 steps (achieved / theoretical occupancy) and cuobjdump (registers, stack frame); *local ld+st* = local-memory
  bytes loaded + stored per launch (ncu `l1tex__t_bytes_pipe_lsu_mem_local_op_{ld,st}`). `KUpdate` merges two
  functors of the same demangled name: the scalar update with its advection (`dyn/scalars.cpp`, 50-320 us per
  launch) and the specified-zone update (`dyn/lateral_bc.cpp`, 3 us).
- The runs start at 00 UTC on a night-time domain: the RRTMG shortwave does (almost) no work (0.06 ms per step
  in P4); a daytime run adds it.

### Per case

| Case | GPU | stepping s/h | ms/step | kernels ms/step (nsys) | GPU busy | launches/step | host waits/step | D2H / H2D / D2D copies per step | malloc+free/step | device memory MiB | stepping with KTIME s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| P1 FP32 | G4 | 0.290 | 4.83 | 4.23 | 85 % | 396 | 4 | 3 (0.0 MB) / 0 / 8 (18.9 MB) | 0 | 1977 | 0.455 |
| P1 FP32 | L4 | 0.610 | 10.17 | 9.65 | 95 % | 396 | 4 | 3 (0.0 MB) / 0 / 8 (18.9 MB) | 0 | 1000 | 0.765 |
| P1 FP32 | A100 | 0.517 | 8.62 | 7.76 | 94 % | 396 | 4 | 3 (0.0 MB) / 0 / 8 (18.9 MB) | 0 | 1638 | 0.712 |
| P1 FP64 | G4 | 0.812 | 13.53 | 12.74 | 94 % | 396 | 4 | 3 (0.0 MB) / 0 / 8 (37.9 MB) | 0 | 3417 | 0.970 |
| P1 FP64 | L4 | 2.123 | 35.38 | 34.89 | 98 % | 396 | 4 | 3 (0.0 MB) / 0 / 8 (37.9 MB) | 0 | 1632 | 2.268 |
| P1 FP64 | A100 | 0.653 | 10.88 | 11.40 | 93 % | 396 | 4 | 3 (0.0 MB) / 0 / 8 (37.9 MB) | 0 | 2802 | 0.878 |
| P4 FP64 | G4 | 5.500 | 91.67 | 91.38 | 99 % | 633 | 5 | 4 (0.0 MB) / 0 / 15 (64.9 MB) | 0 | 20939 | 5.760 |
| P4 FP64 | L4 | 9.849 | 164.15 | 165.48 | 99 % | 633 | 5 | 4 (0.0 MB) / 0 / 15 (64.9 MB) | 0 | 8284 | 10.185 |
| P4 FP64 | A100 | 4.095 | 68.25 | 67.16 | 99 % | 633 | 5 | 4 (0.0 MB) / 0 / 15 (64.9 MB) | 0 | 16644 | 4.424 |
| C01 FP64 | G4 | 2.111 | 35.18 | 35.57 | 96 % | 641 | 4 | 3 (0.0 MB) / 0 / 21 (97.0 MB) | 0 | 19843 | 2.374 |
| C01 FP64 | L4 | 5.257 | 87.62 | 88.54 | 98 % | 641 | 4 | 3 (0.0 MB) / 0 / 21 (97.0 MB) | 0 | 7188 | 5.522 |
| D0 FP64 | G4 | 2.527 | 42.12 | 42.51 | 97 % | 635 | 4 | 3 (0.0 MB) / 0 / 17 (75.6 MB) | 0 | 21469 | 2.787 |
| D0 FP64 | L4 | 5.727 | 95.45 | 95.61 | 99 % | 635 | 4 | 3 (0.0 MB) / 0 / 17 (75.6 MB) | 0 | 8816 | 5.970 |
| N3 FP32 | G4 | 11.572 | 192.87 | 158.10 | 81 % | 2761 | 404 | 327 (70.5 MB) / 71 (12.9 MB) / 60 (130.8 MB) | 34 | 22835 | 12.671 |
| N3 FP32 | L4 | 21.511 | 358.52 | 267.29 | 74 % | 2761 | 404 | 327 (70.5 MB) / 71 (12.9 MB) / 60 (130.8 MB) | 34 | 10184 | 22.523 |
| N3 FP64 | G4 | 16.192 | 269.87 | 224.67 | 85 % | 2761 | 404 | 327 (140.5 MB) / 71 (25.7 MB) / 60 (261.5 MB) | 34 | 23555 | 17.232 |
| N3 FP64 | L4 | 35.154 | 585.90 | 475.46 | 80 % | 2761 | 404 | 327 (140.5 MB) / 71 (25.7 MB) / 60 (261.5 MB) | 34 | 10900 | 36.180 |

### Top 15 kernels per case

**P1 FP32** (G4; per step of the time loop; occupancy and local memory from ncu of the first 2 steps, registers and stack from cuobjdump; the last columns: the same kernel on the other GPUs)

| Kernel | ms/step | % of kernel time | launches/step | us/launch | achieved / theoretical occupancy % | registers | stack B | local ld+st MB/launch | L4 ms/step | A100 ms/step |
|---|---|---|---|---|---|---|---|---|---|---|
| KNoahRs | 1.021 | 24.1 | 1.0 | 1021 | 2 / 100 | 38 | 0 | 0.0 | 1.075 | 1.297 |
| ysu::YSURunner<48> | 0.496 | 11.7 | 1.0 | 496 | 8 / 33 | 128 | 4240 | 268.4 | 0.853 | 0.952 |
| KWSolve | 0.384 | 9.1 | 7.0 | 55 | 8 / 100 | 40 | 32 | 0.0 | 0.671 | 0.835 |
| KMuTCol | 0.251 | 5.9 | 7.0 | 36 | 7 / 83 | 48 | 0 | 0.0 | 0.495 | 0.492 |
| HaloMulti | 0.205 | 4.8 | 92.0 | 2 | 31 / 100 | 25 | 0 | 0.0 | 0.675 | 0.452 |
| KCoefW | 0.157 | 3.7 | 3.0 | 52 | 8 / 83 | 44 | 0 | 0.0 | 0.224 | 0.226 |
| KPhyPrep | 0.134 | 3.2 | 1.0 | 134 | 8 / 100 | 40 | 0 | 0.0 | 0.298 | 0.270 |
| KWW | 0.108 | 2.6 | 3.0 | 36 | 8 / 67 | 62 | 0 | 0.0 | (not top 15) | 0.144 |
| KUV | 0.103 | 2.4 | 7.0 | 15 | 73 / 83 | 46 | 0 | 0.0 | 0.567 | 0.267 |
| noah::NoahKernel | 0.091 | 2.2 | 1.0 | 91 | 7 / 25 | 133 | 576 | 13.9 | (not top 15) | (not top 15) |
| KWPRho | 0.078 | 1.8 | 7.0 | 11 | 62 / 100 | 26 | 0 | 0.0 | 0.245 | (not top 15) |
| KTendU | 0.071 | 1.7 | 3.0 | 24 | 70 / 83 | 48 | 0 | 0.0 | 0.266 | 0.168 |
| KUpdate | 0.070 | 1.7 | 37.0 | 2 | 10 / 100 | 17 | 0 | 0.0 | (not top 15) | 0.175 |
| KTendV | 0.067 | 1.6 | 3.0 | 22 | 72 / 83 | 44 | 0 | 0.0 | (not top 15) | 0.159 |
| KMuTPoint | 0.065 | 1.5 | 7.0 | 9 | 82 / 100 | 34 | 0 | 0.0 | 0.388 | 0.149 |

**P1 FP64** (G4; per step of the time loop; occupancy and local memory from ncu of the first 2 steps, registers and stack from cuobjdump; the last columns: the same kernel on the other GPUs)

| Kernel | ms/step | % of kernel time | launches/step | us/launch | achieved / theoretical occupancy % | registers | stack B | local ld+st MB/launch | L4 ms/step | A100 ms/step |
|---|---|---|---|---|---|---|---|---|---|---|
| ysu::YSURunner<48> | 1.293 | 10.1 | 1.0 | 1293 | 8 / 25 | 150 | 8464 | 540.1 | 2.385 | 1.409 |
| KNoahRs | 1.156 | 9.1 | 1.0 | 1156 | 2 / 100 | 40 | 0 | 0.0 | 1.096 | 1.615 |
| KWSolve | 0.875 | 6.9 | 7.0 | 125 | 8 / 75 | 56 | 40 | 0.0 | 1.488 | 1.087 |
| KTendW | 0.697 | 5.5 | 3.0 | 232 | 69 / 83 | 48 | 0 | 0.0 | 2.644 | 0.252 |
| KTendU | 0.692 | 5.4 | 3.0 | 231 | 63 / 75 | 52 | 0 | 0.0 | 2.555 | 0.279 |
| KTendV | 0.670 | 5.3 | 3.0 | 223 | 63 / 75 | 52 | 0 | 0.0 | 2.472 | 0.261 |
| noah::NoahKernel | 0.609 | 4.8 | 1.0 | 609 | 7 / 17 | 206 | 1152 | 27.8 | (not top 15) | (not top 15) |
| KMuTCol | 0.563 | 4.4 | 7.0 | 80 | 7 / 67 | 58 | 0 | 0.0 | (not top 15) | 0.648 |
| KTendT | 0.546 | 4.3 | 3.0 | 182 | 62 / 75 | 50 | 0 | 0.0 | 2.001 | (not top 15) |
| KUV | 0.519 | 4.1 | 7.0 | 74 | 80 / 100 | 40 | 0 | 0.0 | 1.990 | 0.474 |
| KCoefW | 0.476 | 3.7 | 3.0 | 159 | 8 / 67 | 60 | 0 | 0.0 | (not top 15) | 0.337 |
| KPRhoPhi | 0.457 | 3.6 | 4.0 | 114 | 80 / 100 | 36 | 0 | 0.0 | 1.707 | (not top 15) |
| KUpdate | 0.383 | 3.0 | 37.0 | 10 | 9 / 99 | 18 | 0 | 0.0 | 1.315 | 0.246 |
| KWW | 0.380 | 3.0 | 3.0 | 127 | 8 / 42 | 82 | 0 | 0.0 | (not top 15) | (not top 15) |
| KPhyPrep | 0.360 | 2.8 | 1.0 | 360 | 8 / 83 | 48 | 0 | 0.0 | (not top 15) | 0.378 |

**P4 FP64** (G4; per step of the time loop; occupancy and local memory from ncu of the first 2 steps, registers and stack from cuobjdump; the last columns: the same kernel on the other GPUs)

| Kernel | ms/step | % of kernel time | launches/step | us/launch | achieved / theoretical occupancy % | registers | stack B | local ld+st MB/launch | L4 ms/step | A100 ms/step |
|---|---|---|---|---|---|---|---|---|---|---|
| KLwPrep | 34.791 | 38.1 | 6.0 | 5799 | 4 / 17 | 192 | 1024 | 0.0 | 25.669 | 17.182 |
| KLwBand<60> | 20.891 | 22.9 | 6.0 | 3482 | 7 / 42 | 88 | 1920 | 484.3 | 49.188 | 21.932 |
| KLwLayer | 6.276 | 6.9 | 6.0 | 1046 | 14 / 17 | 250 | 64 | 0.2 | 18.624 | 7.498 |
| wsm6::Wsm6Kernel<double, 128> | 5.342 | 5.8 | 1.0 | 5342 | 7 / 25 | 148 | 49216 | 872.6 | 9.295 | 3.242 |
| cu_kf::KFRunner<double, 48> | 3.518 | 3.8 | 1.0 | 3518 | 5 / 25 | 168 | 21712 | 179.0 | 4.592 | 2.083 |
| KUpdate | 2.975 | 3.3 | 47.0 | 63 | 22 / 96 | 24 | 0 | 0.0 | 10.581 | 0.938 |
| ysu::YSURunner<48> | 1.438 | 1.6 | 1.0 | 1438 | 8 / 25 | 150 | 8464 | 655.5 | 2.517 | 1.413 |
| KNoahRs | 1.168 | 1.3 | 1.0 | 1168 | 2 / 100 | 40 | 0 | 0.0 | (not top 15) | 1.371 |
| KWSolve | 0.879 | 1.0 | 7.0 | 126 | 8 / 75 | 56 | 40 | 0.0 | 1.446 | 0.990 |
| KLwFinish | 0.837 | 0.9 | 6.0 | 139 | 4 / 67 | 64 | 0 | 0.0 | (not top 15) | 0.620 |
| KTendU | 0.689 | 0.8 | 3.0 | 230 | 63 / 75 | 52 | 0 | 0.0 | 2.463 | (not top 15) |
| KTendW | 0.684 | 0.7 | 3.0 | 228 | 69 / 83 | 48 | 0 | 0.0 | 2.513 | (not top 15) |
| KTendV | 0.669 | 0.7 | 3.0 | 223 | 63 / 75 | 52 | 0 | 0.0 | 2.385 | (not top 15) |
| KPRhoPhi | 0.567 | 0.6 | 5.0 | 113 | 80 / 100 | 36 | 0 | 0.0 | 2.078 | (not top 15) |
| KMuTCol | 0.566 | 0.6 | 7.0 | 81 | 7 / 67 | 58 | 0 | 0.0 | (not top 15) | 0.579 |

**C01 FP64** (G4; per step of the time loop; occupancy and local memory from ncu of the first 2 steps, registers and stack from cuobjdump; the last columns: the same kernel on the other GPUs)

| Kernel | ms/step | % of kernel time | launches/step | us/launch | achieved / theoretical occupancy % | registers | stack B | local ld+st MB/launch | L4 ms/step |
|---|---|---|---|---|---|---|---|---|---|
| mp_thompson::ThompsonKernel<double, 48> | 11.838 | 33.3 | 1.0 | 11838 | 8 / 17 | 172 | 63040 | 1273.7 | 17.501 |
| KUpdate | 6.045 | 17.0 | 59.0 | 102 | 31 / 93 | 29 | 0 | 0.0 | 21.381 |
| ysu::YSURunner<48> | 1.393 | 3.9 | 1.0 | 1393 | 8 / 25 | 150 | 8464 | 655.0 | 2.502 |
| KNoahRs | 1.163 | 3.3 | 1.0 | 1163 | 2 / 100 | 40 | 0 | 0.0 | (not top 15) |
| KWSolve | 0.881 | 2.5 | 7.0 | 126 | 8 / 75 | 56 | 40 | 0.0 | 1.462 |
| KFacePDX | 0.752 | 2.1 | 12.0 | 63 | 68 / 83 | 44 | 0 | 0.0 | 2.663 |
| KTendU | 0.689 | 1.9 | 3.0 | 230 | 63 / 75 | 52 | 0 | 0.0 | 2.466 |
| KTendW | 0.683 | 1.9 | 3.0 | 228 | 69 / 83 | 48 | 0 | 0.0 | 2.516 |
| KTendV | 0.669 | 1.9 | 3.0 | 223 | 63 / 75 | 52 | 0 | 0.0 | 2.386 |
| KFacePDY | 0.641 | 1.8 | 12.0 | 54 | 71 / 83 | 44 | 0 | 0.0 | 2.267 |
| noah::NoahKernel | 0.610 | 1.7 | 1.0 | 610 | 7 / 17 | 206 | 1152 | 27.8 | (not top 15) |
| KUpdatePD | 0.572 | 1.6 | 12.0 | 48 | 79 / 100 | 38 | 0 | 0.0 | 1.992 |
| KPRhoPhi | 0.569 | 1.6 | 5.0 | 114 | 80 / 100 | 36 | 0 | 0.0 | 2.090 |
| KMuTCol | 0.568 | 1.6 | 7.0 | 81 | 7 / 67 | 58 | 0 | 0.0 | (not top 15) |
| dyn::diffusion::KHorizDiffM1 | 0.549 | 1.5 | 12.0 | 46 | 68 / 83 | 44 | 0 | 0.0 | 1.977 |

**D0 FP64** (G4; per step of the time loop; occupancy and local memory from ncu of the first 2 steps, registers and stack from cuobjdump; the last columns: the same kernel on the other GPUs)

| Kernel | ms/step | % of kernel time | launches/step | us/launch | achieved / theoretical occupancy % | registers | stack B | local ld+st MB/launch | L4 ms/step |
|---|---|---|---|---|---|---|---|---|---|
| mp_thompson::ThompsonKernel<double, 48> | 10.550 | 24.8 | 1.0 | 10550 | 8 / 17 | 172 | 63040 | 1276.8 | 15.670 |
| KUpdate | 4.332 | 10.2 | 51.0 | 85 | 26 / 95 | 26 | 0 | 0.0 | 15.250 |
| dyn::diffusion::KHorizDiffS2 | 3.483 | 8.2 | 9.0 | 387 | 8 / 33 | 118 | 6144 | 93.2 | 6.062 |
| KTiedtke | 3.024 | 7.1 | 2.0 | 1512 | 7 / 42 | 96 | 49424 | 479.2 | 5.634 |
| KLwPrep | 2.286 | 5.4 | 0.4 | 5714 | 4 / 17 | 192 | 1024 | 0.0 | 1.686 |
| myjpbl::MYJPBLKernel | 1.353 | 3.2 | 1.0 | 1353 | 8 / 33 | 104 | 19632 | 670.2 | 2.291 |
| KLwBand<60> | 1.320 | 3.1 | 0.4 | 3300 | 7 / 42 | 88 | 1920 | 453.5 | 3.010 |
| KNoahRs | 1.169 | 2.8 | 1.0 | 1169 | 2 / 100 | 40 | 0 | 0.0 | (not top 15) |
| KWSolve | 0.878 | 2.1 | 7.0 | 125 | 8 / 75 | 56 | 40 | 0.0 | (not top 15) |
| KTendU | 0.691 | 1.6 | 3.0 | 230 | 63 / 75 | 52 | 0 | 0.0 | 2.473 |
| KTendV | 0.670 | 1.6 | 3.0 | 223 | 63 / 75 | 52 | 0 | 0.0 | 2.394 |
| KTendW | 0.647 | 1.5 | 3.0 | 216 | 68 / 83 | 48 | 0 | 0.0 | 2.393 |
| KMuTCol | 0.570 | 1.3 | 7.0 | 81 | 7 / 67 | 58 | 0 | 0.0 | (not top 15) |
| KPRhoPhi | 0.569 | 1.3 | 5.0 | 114 | 80 / 100 | 36 | 0 | 0.0 | 2.104 |
| KGwdo | 0.546 | 1.3 | 1.0 | 546 | 8 / 33 | 116 | 24032 | 363.2 | 1.591 |

**N3 FP32** (G4; per step of the time loop; occupancy and local memory from ncu of the first 2 steps, registers and stack from cuobjdump; the last columns: the same kernel on the other GPUs)

| Kernel | ms/step | % of kernel time | launches/step | us/launch | achieved / theoretical occupancy % | registers | stack B | local ld+st MB/launch | L4 ms/step |
|---|---|---|---|---|---|---|---|---|---|
| KLwPrep | 68.678 | 43.4 | 12.0 | 5723 | 4 / 17 | 192 | 1024 | 0.0 | 50.767 |
| KLwBand<60> | 38.786 | 24.5 | 12.0 | 3232 | 7 / 42 | 88 | 1920 | 457.2 | 96.511 |
| KLwLayer | 12.236 | 7.7 | 12.0 | 1020 | 15 / 17 | 250 | 64 | 0.0 | 37.165 |
| wsm6::Wsm6Kernel<float, 128> | 7.115 | 4.5 | 4.0 | 1779 | 8 / 25 | 136 | 24608 | 371.1 | 11.475 |
| KNoahRs | 3.824 | 2.4 | 4.0 | 956 | 2 / 100 | 38 | 0 | 0.0 | 3.848 |
| cu_kf::KFRunner<float, 48> | 2.485 | 1.6 | 4.0 | 621 | 7 / 33 | 128 | 10960 | 71.1 | 4.199 |
| ysu::YSURunner<48> | 1.910 | 1.2 | 4.0 | 478 | 8 / 33 | 128 | 4240 | 327.3 | 3.570 |
| nest::KBdySlabView3D | 1.892 | 1.2 | 44.0 | 43 | 9 / 42 | 96 | 0 | 0.0 | 4.332 |
| KLwFinish | 1.676 | 1.1 | 12.0 | 140 | 4 / 67 | 64 | 0 | 0.0 | (not top 15) |
| KWSolve | 1.647 | 1.0 | 28.0 | 59 | 8 / 100 | 40 | 32 | 0.0 | 2.667 |
| KUpdate | 1.173 | 0.7 | 209.0 | 6 | 23 / 100 | 20 | 0 | 0.0 | 3.848 |
| KMuTCol | 1.151 | 0.7 | 28.0 | 41 | 7 / 83 | 48 | 0 | 0.0 | 1.927 |
| HaloMulti | 0.954 | 0.6 | 406.0 | 2 | 34 / 100 | 25 | 0 | 0.0 | 3.382 |
| KRadColumn | 0.782 | 0.5 | 2.0 | 391 | 8 / 83 | 44 | 512 | 6.7 | (not top 15) |
| KWPRho | 0.701 | 0.4 | 28.0 | 25 | 61 / 100 | 26 | 0 | 0.0 | (not top 15) |

**N3 FP64** (G4; per step of the time loop; occupancy and local memory from ncu of the first 2 steps, registers and stack from cuobjdump; the last columns: the same kernel on the other GPUs)

| Kernel | ms/step | % of kernel time | launches/step | us/launch | achieved / theoretical occupancy % | registers | stack B | local ld+st MB/launch | L4 ms/step |
|---|---|---|---|---|---|---|---|---|---|
| KLwPrep | 66.463 | 29.6 | 12.0 | 5539 | 4 / 17 | 192 | 1024 | 0.0 | 50.850 |
| KLwBand<60> | 37.798 | 16.8 | 12.0 | 3150 | 7 / 42 | 88 | 1920 | 496.7 | 96.635 |
| wsm6::Wsm6Kernel<double, 128> | 20.284 | 9.0 | 4.0 | 5071 | 8 / 25 | 148 | 49216 | 863.6 | 38.434 |
| KLwLayer | 11.937 | 5.3 | 12.0 | 995 | 15 / 17 | 250 | 64 | 0.2 | 37.150 |
| KUpdate | 11.654 | 5.2 | 209.0 | 56 | 21 / 96 | 24 | 0 | 0.0 | 44.788 |
| cu_kf::KFRunner<double, 48> | 10.215 | 4.5 | 4.0 | 2554 | 6 / 25 | 168 | 21712 | 154.6 | 13.246 |
| ysu::YSURunner<48> | 5.400 | 2.4 | 4.0 | 1350 | 8 / 25 | 150 | 8464 | 657.3 | 9.830 |
| KNoahRs | 4.424 | 2.0 | 4.0 | 1106 | 2 / 100 | 40 | 0 | 0.0 | (not top 15) |
| KWSolve | 3.423 | 1.5 | 28.0 | 122 | 8 / 75 | 56 | 40 | 0.0 | 5.912 |
| KTendU | 2.648 | 1.2 | 12.0 | 221 | 63 / 75 | 52 | 0 | 0.0 | 10.204 |
| KTendW | 2.623 | 1.2 | 12.0 | 219 | 69 / 83 | 48 | 0 | 0.0 | 10.400 |
| KTendV | 2.576 | 1.1 | 12.0 | 215 | 63 / 75 | 52 | 0 | 0.0 | 9.931 |
| KMuTCol | 2.245 | 1.0 | 28.0 | 80 | 7 / 67 | 58 | 0 | 0.0 | (not top 15) |
| KPRhoPhi | 2.161 | 1.0 | 20.0 | 108 | 80 / 100 | 36 | 0 | 0.0 | 8.607 |
| KTendT | 2.084 | 0.9 | 12.0 | 174 | 63 / 75 | 50 | 0 | 0.0 | 7.978 |


### What the profile shows

- **The step is GPU-bound, not launch-bound.** Kernels keep the GPU busy 85-99 % of the step on a single domain
  (P1 FP32 on the G4: 4.23 of 5.00 ms; every physics case >= 96 %). 396-641 launches per step cost ~2.5 us of host
  time each but overlap with the running kernels. The idle part bounds what CUDA Graphs and launch fusion can win:
  at most 15 % for P1 FP32 on the G4, 5 % on the L4, 6 % / 2 % for P1 FP64, 1-4 % with physics.
- **RRTMG longwave** (radt = 1 in P4 and N3): 63 of 91 ms per step in P4 on the G4 (69 %), 118 of 225 ms in N3
  FP64. Its batch runs in chunks of 2048 columns (radiation.cpp:352): 6 chunks for 10 591 columns, each chunk's
  `KLwPrep` and `KLwFinish` a grid of 32 blocks x 64 threads, i.e. 32 of the G4's 188 SMs with 2 warps each
  (achieved occupancy 4 %), and `KLwBand` 32 x 16 blocks. `KLwPrep` is therefore *faster* on the L4 (25.7 ms per
  step) than on the G4 (34.8 ms): it is latency-bound serial per-column work and the bigger GPU stays idle. The
  `KLwBand` stack frame (1920 B) moves ~480 MB of local memory per launch.
- **`KNoahRs` is one GPU thread** walking the whole domain (surface_driver.cpp:109-126, launched on
  `Range2D{0, 1, 0, 1}`): 1.0-1.2 ms per step in every Noah case on both GPUs, 24 % of P1 FP32's kernel time on the
  G4 and 9 % of P1 FP64's. It only copies RC values (no arithmetic), so a parallel per-row carry (last active RC of
  each row, then a scan over rows) is bit-identical by construction.
- **Column physics run at 5-8 % achieved occupancy** on this domain (one thread per column: 90 blocks of 128 threads
  for 188 SMs) with large stack frames: Thompson (FP64, 48-level tier) 63 KB stack, 1.27 GB local-memory traffic per
  launch, 10.5-11.8 ms per step (25-33 % of C01 / D0); WSM6 49 KB, 0.87 GB, 5.3 ms; MYJ PBL 20 KB; Tiedtke 49 KB;
  YSU 8.5 KB and 0.55-0.66 GB. On bigger domains the occupancy rises with the column count; the local-memory
  traffic per column stays.
- **Device memory**: P1 (no WSM6 / Thompson) 2.0 GB FP32 and 3.4 GB FP64 on the G4; with WSM6 or Thompson 19.8-21.5
  GB (G4) and 7.2-8.8 GB (L4): the 64 KB stack limit those schemes set (physics_coupling.cpp:128, 145) reserves
  local memory for every resident thread of the GPU, ~16-18 GB on the G4's 188 SMs and ~5.6-7.2 GB on the L4.
  N3 FP64 needs 2.6 GB more than P4.
- **Two-way nest (N3)**: per parent step 327 device-to-host and 71 host-to-device `cudaMemcpy` (140 / 26 MB in FP64),
  17 `cudaMalloc` + 17 `cudaFree` + 17 `cudaMemset` (the smoother temporaries, nest_feedback.hpp:416, 466) and 404
  host waits; the GPU idles 15-19 % of the step on the G4 and 20-26 % on the L4. The single-domain cases have none
  of these in the time loop (3-4 4-byte fatal-flag reads and one `dev_sync` per step, 8-21 device-to-device copies
  from `copy_dev`).
- **Dycore**: the four tendency kernels (`KTendU/V/W/T`) take 0.22-0.23 ms per launch in FP64 on the G4 (63-69 %
  achieved occupancy, 82 % SM throughput: compute-bound on the 1/64-rate FP64 units) and 3.5x longer on the L4; the
  scalar transport (`KUpdate` of scalars.cpp, advection fluxes recomputed per cell, one launch per species and
  stage) 3-6 ms per step with moist species (10-17 % of C01 / D0). The acoustic column kernels (`KWSolve`,
  `KMuTCol`) run at 7-8 % occupancy, like the column physics.
- **Half-rate FP64 (A100)**: P1 FP64 runs 1.26x slower than P1 FP32 (G4: 2.8x, L4: 3.5x) and faster than on the G4
  (0.65 vs 0.81 s/h) although FP32 is slower (0.52 vs 0.29 s/h). The dycore tendency kernels drop out of the top
  kernels; what remains on top are the single-thread `KNoahRs` (1.6 ms per step, 14 % of P1 FP64), YSU and the
  latency-bound column kernels. P4: 4.1 s/h (G4 5.5), RRTMG LW still 47 of 67 ms per step (69 %). On such GPUs the
  FP64 arithmetic changes (rank 5 below) gain little; parallelism (ranks 1-3) is what counts.
- **The AR1 table (benchmarks) is confirmed** for the single-domain cases (P1 FP32 10.4 ns/cell/step on the G4,
  P4 197, D0 90, C01 76). For N3 it listed the d02 line of the bench output: the stepping time of the whole nested
  run is the parent's, 16.2 s/h FP64 and 11.6 s/h FP32 on the G4 (35.2 / 21.5 s/h on the L4), not 17.6 s/h (L4) and
  4.5 s/h (G4).

### AR1's top 10 (section 2.2), reordered by the measured gain

Gains are upper bounds derived from the measured time of the kernels or idle periods a change removes (G4 unless
stated); none is measured after a change yet.

| New rank | AR1 rank | Change | Measured basis | Gain (estimate from the profile) | AR1 estimate |
|---|---|---|---|---|---|
| 1 | 3 | RRTMG batch: all columns in one chunk (or a column x g-point grid for `KLwPrep` / `KLwFinish`), local memory of `KLwBand` | LW 63 / 91 ms (P4), 118 / 225 ms (N3 FP64); `KLwPrep` uses 32 of 188 SMs | P4 ~1.7x on the G4 (prep and finish once instead of 6 times: -30 ms, band -9 ms), ~1.15x on the L4 (58 SMs: the 32-block grids already fill more of it); operational radt (D0, radt 15): ~6 % | P4 1.5-2x, radt 10-25 % |
| 2 | (new) | `KNoahRs` in parallel (per-row carry + scan over rows) | 1.0-1.2 ms per step, one thread | P1 FP32 1.25x (G4), 1.12x (L4); P1 FP64 1.09x (G4), 1.03x (L4), 1.16x (A100); physics cases 1-3 % | - |
| 3 | 4 | Column-physics local memory: exact-nz tiers, k-major workspaces for Thompson / WSM6 / YSU / MYJ / Tiedtke / KF, lower stack limit | Thompson 25-33 % of C01 / D0, 1.27 GB local traffic per launch; 16-18 GB reserve (G4) | C01 / D0 10-25 % if the local traffic halves; device memory -16 GB (G4), -6 GB (L4) | Thompson 1.3-2x, C01 / D0 10-30 %, GBs |
| 4 | 1 | Nest feedback on the device, no allocations in the loop | N3: 327 + 71 blocking copies, 34 malloc / free, GPU idle 15-19 % (G4), 20-26 % (L4) | N3 15-19 % (G4), 20-26 % (L4) | N3 15-30 % |
| 5 | 6 | Face fluxes once (KTendU/V/W/T and the scalar advection in `KUpdate` / PD faces), `divv` once in KWW | KTend* 2.6 ms of 12.7 (P1 FP64); scalar `KUpdate` 3-6 ms (moist cases) | ~8 % (P1 FP64, G4 / L4), 10-15 % with many species (C01 / D0); little in FP32 and on half-rate FP64 GPUs (A100) | FP64 dycore 1.2-1.4x |
| 6 | 2 | CUDA Graphs per step kind, one fatal check per step | GPU idle 15 % P1 FP32 (G4), 5 % (L4); <= 4 % with physics | P1 FP32 <= 1.18x (G4), <= 1.05x (L4); FP64 / physics <= 6 % | P1 FP32 1.3-1.8x |
| 7 | 7 | Acoustic-step fusion | `KWSolve` + `KMuTCol` + `KCoefW` + `KUV` + `KWPRho` + `KMuTPoint` ~3 ms of 12.7 (P1 FP64), column kernels at 8 % occupancy | 5-10 % of the dycore (the column kernels are parallelism-bound, fusion saves traffic only) | dycore 10-25 % |
| 8 | 5 | Launch fusion (zero_tendf, AC* accumulations, couple / decouple, multi-species scalars, batched halos) | 92 `HaloMulti` (0.2 ms) and 37-209 small launches per step; part of the idle time of rank 6 | <= 5 % (P1 FP32), overlaps with rank 6 | P1 FP32 10-20 % |
| 9 | 8 | Allocate suite fields by package | P1 3.4 GB FP64 (7.3 KB per cell incl. context) | memory only | memory -20..-40 % |
| 10 | 9 | Fused asynchronous diagnostics | not in the bench runs (diag_interval off); 1 + n_q blocking reductions every 10 steps | <= 1-3 % | 1-3 % |
| 11 | 10 | Shared device table cache | start-up and nest memory | memory / start-up only | - |

### Nest feedback on the device (AR1 rank 1, F02; measured)

`med_nest_feedback`'s smoothers keep their temporaries per parent grid (`nest::FeedbackScratch`), and the physics-state
feedback of the history view (`WrfHistory::feedback_physics`) runs as device kernels on the parent's and the nest's
fields (pairs found once per nest, integer fields on the device); `feedback_from` asks once which fields the model
carries. N3 on the G4 (commit 7aa0970 vs the change, same session, each run twice, as the table above):

| Case | before s/h | after s/h | gain | D2H / H2D / memset per parent step (nsys) | malloc+free/step |
|---|---|---|---|---|---|
| N3 FP64 | 15.78, 15.86 | 14.46, 14.49 | 1.09x | 327 (140.5 MB) / 71 / 17 -> 28 (20.0 MB) / 0 / 0 | 34 -> 0 |
| N3 FP32 | 11.48, 11.49 | 10.02, 10.02 | 1.15x | (not profiled) | 34 -> 0 |

Kernel time per parent step is unchanged (227.0 / 227.2 ms, +71 small feedback launches): the gain is the host time
the GPU waited for. Every wrfout / wrfrst stream of N1, N2 (smooth_option 0, 1, 2) and N3 (1 h, wrfrst at 30 min) is
identical before and after at tolerance 0 in both precisions. The remaining 28 device-to-host copies per parent step
(apart from the fatal-flag reads) were not attributed here (HT_COARSE and d02's native history every 30 nest steps
are candidates).

## CUDA graphs for the time step (measured, Colab)

`dyn::RK3::step` runs inside a `StepScope` (core/backend.hpp): its launches and device copies are captured into
CUDA graphs on the one launch stream and replayed with `cudaGraphExecUpdate` (one instantiated graph per domain and
segment, up to 4 topologies each); an operation that waits for the device ends the segment. The fatal flags of the
step are copied asynchronously and read once when the host next waits; the host captures the next step while the
device runs this one (nests: up to the feedback's wait once per parent step). Same kernels, same order, same arguments: P1 has 1 segment of 400 nodes per step
(60 launches, 58 updates, 3 instantiations per hour), P4 1 of 652, N3 4 per parent step (d01 and 3 d02 steps).
`CPPWRF_GRAPHS=0` launches as before (same timings as origin/main), `CPPWRF_GRAPHS=stats` prints the counts.

Stepping s/h, 1 h, no history output, base (origin/main 7aa0970) / new / new with `CPPWRF_GRAPHS=0`, each twice,
interleaved in one session per GPU (L4: P1 FP32 with the run_P1_r8 input, two sessions):

| Case | GPU | base | new | speed-up | without graphs |
|---|---|---|---|---|---|
| P1 FP32 | G4 | 0.287, 0.287 | 0.255, 0.255 | 1.13x | 0.286, 0.287 |
| P1 FP64 | G4 | 0.846, 0.846 | 0.812, 0.812 | 1.04x | 0.845, 0.845 |
| P4 FP64 | G4 | 5.244, 5.244 | 5.189, 5.191 | 1.01x | 5.243, 5.242 |
| N3 FP64 | G4 | 15.739, 15.752 | 15.540, 15.557 | 1.01x | 15.752, 15.700 |
| N3 FP32 | G4 | 11.401, 11.406 | 11.266, 11.272 | 1.01x | 11.449, 11.454 |
| P1 FP32 | L4 | 0.604, 0.626, 0.632, 0.629 | 0.601, 0.609, 0.617, 0.619 | 1.02x | 0.609, 0.622, 0.628, 0.631 |
| P1 FP64 | L4 | 2.225, 2.313 | 2.218, 2.310 | 1.00x | 2.231, 2.323 |
| P4 FP64 | L4 | 9.939, 10.091 | 9.944, 10.095 | 1.00x | 9.965, 10.097 |
| N3 FP64 | L4 | 35.370, 35.879 | 35.485, 35.888 | 1.00x | 35.714, 35.958 |
| N3 FP32 | L4 | 21.451, 21.443 | 21.482, 21.382 | 1.00x | 21.511, 21.458 |

The L4 was already 95-99 % busy (profile above) and its host (Xeon 2.2 GHz) is slower: the graphs win little there.
Before the nests stopped waiting after each domain step, N3 on the L4 was 1 % slower with graphs (the capture was
not overlapped with the device's work).

nsys, P1 FP32 on the G4, per step: kernels 4.15 ms (base) / 4.13 ms (new), 396 kernels either way; the host makes 1
graph launch instead of 396 kernel launches, 1 event wait instead of 3 `cudaStreamSynchronize` (fatal flags) and 1
`cudaDeviceSynchronize`. The step now takes 4.25 ms for 4.13 ms of kernels (GPU busy ~97 %, before 87 %), close to
the bound of the profile above (rank 6). With physics the GPU was already 96-99 % busy, so the gain there is ~1 %;
the nests still wait once per parent step for the feedback.

## Launch fusion (AR1 rank 5, F04; measured, Colab G4)

Independent kernels of a step run as one launch (`LaunchBatch` / `parallel_fused`, core/backend.hpp): every part
keeps its range and parallel_for's thread mapping in x and y and gets its own blocks along z; the parts must not
touch each other's data. Batched: init_zero_tendency (all held tendencies), the lateral boundary kernels
(relax_bdy_dry, spec_bdy_dry, spec_bdyupdate of u / v and t / mu / muts, spec_bdyupdate_ph with zero_grad_bdy or
spec_bdyupdate(w), spec_bdy_final, the scalars' relax / spec / flow_dep_bdy), the scalar transport (each kernel of
rk_update_scalar / rk_update_scalar_pd for up to 8 species per launch; every species keeps its operations in WRF's
order, with its own work fields, allocated before the first step), calculate_phy_tend's couple and
phy_prep_part2's decouple, update_phy_ten (one launch per term of a tendency: each tendency gets its terms in
phy_ra_ten / phy_bl_ten / phy_cu_ten order), the 16 radiation accumulations, the small step's ru_m / rv_m edge sums,
and adjacent halo fills of single fields (`fill_halos` of `HaloField`: both ranks, explicit parities).

CUDA FP32 builds keep nvcc's default `--fmad=true`, and ptxas contracts `a*b + c` into an fma depending on the code
around it: `test_launch_batch` with a multiply-add kernel gave different FP32 values batched than alone. FP32 builds
therefore batch only kernels without multiply-add (`kNoMulAdd`: zeroing, spec_bdytend, flow_dep_bdy, zero_grad_bdy,
halo copies); the others launch as before. FP64 builds (`--fmad=false`) and the CPU (each part its own parallel_for)
batch everything.

Stepping s/h, 1 h, no history output, base origin/main 2d49050 / port/F04 564a5d7, each twice, interleaved in one
session; graph nodes per step = kernels + device copies (`CPPWRF_GRAPHS=stats`; N3: per parent step, 4 graph
launches); peak device memory from nvidia-smi:

| Case | base s/h | new s/h | gain | graph nodes per step | peak MiB |
|---|---|---|---|---|---|
| P1 FP32 | 0.251, 0.252 | 0.248, 0.247 | 1.02x | 400 -> 352 | 1981 -> 1979 |
| P1 FP64 | 0.806, 0.802 | 0.797, 0.794 | 1.01x | 400 -> 300 | 4059 -> 4061 |
| P4 FP64 | 5.137, 5.130 | 5.110, 5.104 | 1.005x | 652 -> 420 | 20953 -> 21193 |
| N3 FP64 | 14.115, 14.121 | 13.956, 13.944 | 1.01x | 2539 -> 1557 | 23599 -> 24075 |
| N3 FP32 | 9.689, 9.681 | 9.665, 9.663 | 1.00x | 2539 -> 2145 | 22857 -> 23177 |

With every kernel batched in FP32 too (commit e83ebd3, same session type) P1 FP32 had 300 nodes and 0.249 s/h and
N3 FP32 1557 nodes and 9.71 s/h (base 9.81-9.82 in that session); its FP32 outputs were byte-identical to base in
the suite cases, but the contraction cannot be proved for every code path, so FP32 keeps the restriction. As the
CUDA-graph measurements predicted, the step was already 97-99 % GPU-busy: fusion removes 25-40 % of the launches
and gains 0-2 %. The memory grows by the scalar work fields of up to 8 species (P4 FP64 +240 MiB, N3 +476 MiB).
Proof: every cuda_dp case against WRF at tolerance 0 (79 PASS), the 16 FP32 cases with base and new builds (120
snapshots and 3 restart files byte-identical), the CPU suite cases of the touched code at tolerance 0.

