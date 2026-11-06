# CPPWRF performance

## WRF v4.8.0 and CPPWRF on one Colab G4 runtime (measured)

Both models ran on the same Google Colab G4 runtime: an AMD EPYC 9B45 (24 cores, 48 hardware threads, 176 GB) and
an NVIDIA RTX PRO 6000 Blackwell Server Edition (96 GB, sm_120, driver 580.82, CUDA 13.0). `tools/bench/colab_bench.sh`
does the whole benchmark in one session (about 75 minutes, 12 compute units at 10.05 units per hour). The numbers
below are from that run, with CPPWRF at main b30ec5b.

**Builds.**
- WRF v4.8.0 is built on the runtime from the official release source (`v4.8.0.tar.gz` of github.com/wrf-model/WRF,
  md5 dbd299cf76f4c0b8d6523fd41cd77674). It uses WRF's own `configure` for GNU in the Debian 13 chroot of the suite
  (gfortran 14.2, MPICH 4.2.1, `FCOPTIM = -O2 -ftree-vectorize -funroll-loops`, no `-march`).
- WRF has four builds: dmpar (option 34) and dm+sm (option 35), each in real*4 (the usual production build) and
  real*8 (`configure -r8`, the suite's reference build).
- CPPWRF has three builds: CUDA FP32 and CUDA FP64 (the release builds of README.md), and the CPU FP64 build with
  OpenMP (gcc 14.2, `-O3 -march=native`). The CPU build runs with 48 threads.

**Cases.** Every run is 1 simulated hour. The same namelist goes to `wrf.exe` and to `tools/wrf2cppwrf.py`.
- R1: dynamics only.
- P4: KF, WSM6, YSU, Noah and RRTMG with radt 1 min.
- P4R15: P4 with radt 15 min, as production runs use.
- N3: a two-way nest with the P4 physics on both domains. d02 has 120 x 90 cells at 5 km with dt 20 s.
- L: a large domain, P4 physics with radt 15 min.

R1, P4 and N3 have 119 x 89 x 44 cells, 15 km, dt 60 s. L has 359 x 269 x 44 cells (9.1 times P4) at 15 km. It is
a Lambert domain centred at 39 N 35 E, made on the runtime:
- WPS 4.6.0's geogrid, ungrib and metgrid (the executables that made the suite's met_em files);
- WPS's low-resolution mandatory geography and the GFS 0.25 degree analysis and 6 h forecast of 2026-10-12 00 UTC
  (downloaded on the runtime);
- real.exe of the real*8 build.

The WPS steps took about 3 minutes with the downloads and real.exe 9 s, so a real-data domain was cheaper than a
grown idealized case.

All models read the same real*8 `wrfinput` / `wrfbdy` files. The real*4 builds of WRF and CPPWRF convert them on
reading, as WRF's NetCDF reader allows. The real*4 WRF runs read WRF's single-precision RRTMG tables.

**How the times are measured.**
- WRF's *stepping* time is the sum of its "Timing for main" lines for domain 1, minus the "Timing for Writing" lines
  inside them. A nest's steps are part of its parent's time.
- CPPWRF's stepping time is its first "cppwrf: done" line. That line is the parent domain's, so it contains the
  nests' steps. Its output is excluded, except for N3, where the parent's step contains d02's history output.
- *Wall* is the whole process (`mpiexec` or `cppwrf`). It includes start-up, reading the inputs, initialisation
  and output.
- *History off* means `history_interval = 1440`: only the frame at t = 0 is written. *History on* means a frame
  every 10 minutes (7 frames per domain).
- WRF ran with five process layouts, all with history off: 12, 24 and 48 MPI ranks, and the dm+sm build with
  24 x 2 and 12 x 4 ranks x threads. Ranks are bound to cores, or to hardware threads beyond 24 ranks. The fastest
  layout for each case and precision is "WRF best" below; it then ran again with history on.
- The GPU's history-off runs ran twice and the faster run counts (the first loads the module). Every other run ran
  once.
- An earlier session measured the same cases. Against it:
  - the GPU stepping times agree within 3.5 %;
  - WRF's times at the fastest layout agree within 3 % (P4 real*4: 9 %); single runs at other layouts vary up to
    15 %, and up to 90 % for the short 48-rank runs;
  - the CPU FP64 times agree within 8 % (R1: 23 %).

Stepping time, in seconds per simulated hour (history off):

| Case | WRF real*4 best | WRF real*8 best | CPPWRF GPU FP32 | CPPWRF GPU FP64 | CPPWRF CPU FP64 | GPU FP32 vs WRF r4 | GPU FP64 vs WRF r8 |
|---|---|---|---|---|---|---|---|
| R1 | 0.86 (12 x 4) | 2.15 (24) | 0.12 | 0.51 | 3.77 | 7.3x | 4.2x |
| P4 | 7.84 (24 x 2) | 9.72 (24 x 2) | 1.66 | 2.78 | 10.86 | 4.7x | 3.5x |
| P4R15 | 2.38 (12 x 4) | 4.34 (12 x 4) | 0.55 | 1.68 | 7.44 | 4.3x | 2.6x |
| N3 | 23.79 (24 x 2) | 34.02 (24) | 5.25 | 10.39 | 38.80 | 4.5x | 3.3x |
| L | 24.78 (24 x 2) | 45.06 (24) | 3.01 | 9.59 | 51.02 | 8.2x | 4.7x |

Time per step and per cell, in ns per cell per step (stepping time over the sum of cells x steps of all domains),
history off:

| Case | WRF r4 | WRF r8 | GPU FP32 | GPU FP64 | CPU FP64 |
|---|---|---|---|---|---|
| R1 | 30.6 | 76.8 | 4.2 | 18.2 | 134.7 |
| P4 | 280 | 348 | 59.3 | 99.4 | 389 |
| P4R15 | 85.1 | 155 | 19.7 | 60.2 | 266 |
| N3 | 210 | 300 | 46.3 | 91.5 | 342 |
| L | 97.2 | 177 | 11.8 | 37.6 | 200 |

Wall time of the whole run in seconds, history off / history on, and Colab compute units per simulated hour
(history off; the runtime costs the same whichever processor a run uses):

| Case | WRF real*4 best | WRF real*8 best | GPU FP32 | GPU FP64 | CPU FP64 | CU / h: WRF r4, GPU FP32 |
|---|---|---|---|---|---|---|
| R1 | 1.6 / 3.6 | 3.3 / 7.1 | 1.1 / 3.0 | 1.9 / 5.7 | 5.4 / 9.9 | 0.005, 0.003 |
| P4 | 9.0 / 12.8 | 11.4 / 21.0 | 3.1 / 6.7 | 4.8 / 12.3 | 12.5 / 21.5 | 0.025, 0.009 |
| P4R15 | 3.4 / 8.3 | 5.8 / 14.4 | 1.9 / 5.5 | 3.6 / 11.2 | 9.1 / 16.6 | 0.010, 0.005 |
| N3 | 25.1 / 30.6 | 36.0 / 46.3 | 7.1 / 11.4 | 12.8 / 21.4 | 41.0 / 47.6 | 0.070, 0.020 |
| L | 31.5 / 63.2 | 55.3 / 123.6 | 10.2 / 41.8 | 22.7 / 89.6 | 63.6 / 131.7 | 0.088, 0.029 |

Time outside the stepping (wall - stepping, in s), history off / on. It is start-up, reading the inputs,
initialisation and output. With history on, most of it is the 7 frames' output. For N3, CPPWRF's d02 output falls
inside the parent's stepping time.

| Case | WRF real*4 | WRF real*8 | GPU FP32 | GPU FP64 | CPU FP64 |
|---|---|---|---|---|---|
| R1 | 0.8 / 2.7 | 1.1 / 4.9 | 1.0 / 2.9 | 1.4 / 5.2 | 1.7 / 5.2 |
| P4 | 1.1 / 5.2 | 1.6 / 9.9 | 1.4 / 5.1 | 2.0 / 9.5 | 1.6 / 9.5 |
| P4R15 | 1.0 / 5.0 | 1.5 / 9.4 | 1.4 / 5.0 | 1.9 / 9.5 | 1.7 / 9.3 |
| N3 | 1.3 / 5.7 | 2.0 / 10.6 | 1.8 / 4.0 | 2.4 / 6.7 | 2.2 / 6.9 |
| L | 6.7 / 38.8 | 10.2 / 78.5 | 7.2 / 38.8 | 13.1 / 80.1 | 12.6 / 80.0 |

GPU memory: nvidia-smi's peak during the history-off run, in MiB. The *stack reserve* (CPPWRF_MEMINFO: stack limit
x 288768 resident threads on the G4) is the column physics' local memory, which CUDA reserves for every thread the
GPU can hold. An RTX 2060 (6 GB, 30 x 1024 resident threads) reserves 9.4 times less. The last column replaces the
G4's reserve by the RTX 2060's. It is an upper bound: CPPWRF sizes the RRTMG batch scratch to the GPU's memory, so
the scratch is smaller on a 6 GB GPU. No run used a 6 GB GPU.

| Case | FP32 peak (stack reserve) | FP32 on an RTX 2060, at most | FP64 peak (stack reserve) | FP64 on an RTX 2060, at most |
|---|---|---|---|---|
| R1 | 803 (282) | 551: fits | 987 (282) | 735: fits |
| P4, P4R15 | 8297 (3018) | 5600: fits | 13749 (7989) | 6610 |
| N3 | 9461 (3018) | 6764 | 15395 (7989) | 8256 |
| L | 20791 (3018) | 18094: does not fit | 29957 (7989) | 22818: does not fit |

**Agreement** (`tests/verify/exact.py`, every variable of every frame, history on):
- CPPWRF CPU FP64 against the suite's WRF real*8 references, which were made on Intel. P4 (7 frames) and N3 (d01 and
  d02, 7 frames each) are bit-identical. R1's reference covers 5 minutes only; the suite case `real_r1_fp64_cpu`
  checks it.
- CPPWRF GPU FP64 against CPPWRF CPU FP64 is bit-identical in all five cases, including N3's nest and the large
  domain L.
- The default CPU build computes glibc libmvec's vector log as libmvec does on Intel (`CPPWRF_PORTABLE_LIBMVEC`),
  so on this AMD CPU it does not match a WRF built here. The CPU FP64 build with this CPU's own libmvec log
  (`-DCPPWRF_PORTABLE_LIBMVEC=OFF`) is bit-identical over the hour to WRF real*8 built and run here in one process,
  for R1, P4 and L. L is a domain that the suite has never seen.
- WRF's own result depends on its process layout. The same real*8 executable run with 24 MPI ranks differs from
  its one-process run. Already in R1 (dynamics only) U differs by up to 1e-11 m/s after 10 minutes, and with the
  physics the difference grows to the size of the FP32 differences below. The most likely cause is the patch
  bounds of the vectorized loops. They decide which points take libmvec's vector functions and which take the
  scalar libm, and the two differ in the last bit. Exact comparisons are therefore made with one-process WRF,
  as the suite's references are.
- FP32: the relative L2 difference at 1 h (||a - b|| / ||b||) of GPU FP32 against WRF real*4 is U 2.0e-4 to
  7.8e-4, W 4.6e-3 to 2.0e-2, T 4.7e-6 to 5.4e-5 and QVAPOR 2.0e-5 to 6.8e-4 over the five cases. WRF real*4
  against WRF real*8 is of the same size (U 1.9e-4 to 7.7e-4, W 3.8e-3 to 1.9e-2). GPU FP32 is about as far from
  WRF real*8 as WRF real*4 is: closer in R1, P4, P4R15 and N3, and within 5 % of it on L.

**What the numbers show.**
- The GPU's advantage grows with the domain. GPU FP32 is 4.3x to 4.7x faster than WRF real*4 at its best layout on
  the 119 x 89 domain with physics, and 8.2x on L. GPU FP64 is 2.6x to 3.5x faster than WRF real*8 there, and 4.7x on
  L. Even FP64 on the GPU steps faster than WRF real*4 on all 48 hardware threads, in every case.
- The whole-run gain is smaller than the stepping gain. Start-up and output are serial in both models and take most
  of a short run.
  - With history off, GPU FP32 against WRF real*4 is 1.5x (R1) to 3.5x (N3) in wall time.
  - With history on, L's 7 frames cost CPPWRF FP32 about 32 s and FP64 about 67 s, as much as they cost WRF real*4
    (32 s) and real*8 (68 s). Output then decides the wall time. NetCDF output was the next bottleneck for large
    domains. Since IO1 the chunks are compressed on host threads, with the same bytes: L's hour with 7 frames and a
    restart file takes 13 s instead of 58 s in FP32 (section "NetCDF output: chunks compressed on host threads").
  - CPPWRF's start-up on L, with history off, takes 7.2 s in FP32 against 6.7 s for WRF real*4. It includes the
    reading of the 4.25 M-cell real*8 inputs.
- WRF gains little beyond one rank per core on this runtime: 48 ranks are slower than 24 in 9 of the 10 case and
  precision pairs (P4 real*4 aside). The dm+sm build was fastest in 7 of the 10, by up to 13 % over the best dmpar
  layout.
- WRF real*4 is only 1.2x to 2.5x faster than real*8, about 1.8x on L.
- CPPWRF's CPU FP64 build (OpenMP, 48 threads) runs at 0.6x to 0.9x the speed of WRF real*8 at its best MPI
  layout: 0.9x on P4, N3 and L, and 0.6x on R1 and P4R15, where the dynamics take a larger share.
- Per Colab compute unit the GPU runs are 1.5x to 3.5x cheaper than WRF on the same runtime, which costs the same
  whatever processor a run uses. WRF alone would cost less on a CPU-only runtime; that runtime was not measured.

The run's report (`report.md`: every layout, both history settings, the time outside the stepping, and every
comparison) and the logs are in the output directory of `colab_bench.sh`.

## Memory-traffic budget

Model: each 3D field touched by a kernel is read (or written) once from DRAM per kernel,
i.e. stencil neighbours are assumed to hit cache; 1D base-state arrays and 2D fields are
negligible (O(1/nz)). Bytes per cell for FP32 (x2 for FP64). These are **counts from the
code, not measurements**; the measured times are in the section above and the profile below.

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

### Column-physics local memory (AR1 rank 4, F06; measured, Colab G4)

The CUDA driver reserves a kernel's stack frame for every thread the GPU can hold at once (SMs x threads per SM:
188 x 1536 = 288768 on the G4, 30 x 1024 = 30720 on an RTX 2060), whatever the launch; the reserve grows at a
launch and never shrinks. The WSM6 and Thompson set-ups raised the stack limit to 64 KB, so every case with them
held 16-18 GB of device memory for stacks on the G4 before the first field.

- The stack limit is the largest frame launched: before a kernel's first launch its frame
  (`cudaFuncGetAttributes`) is reserved (`cudaDeviceSetLimit`, core/backend.cpp `reserve_local_memory`). When the
  device memory is short, the error names the kernel, its frame, the resident threads and the free memory; a failed
  allocation states its size and the free memory. `CPPWRF_MEMINFO=1` prints the device memory after the set-up and
  at the end (with the stack limit and its reserve).
- `core/column_work`: one device buffer shared by the column kernels (`column_for_work`), sized for the launch's
  columns (at most 1/8 of the device memory and 1/4 of what is free, else chunks of row pairs; reserved outside the
  time loop) and laid out as the driver lays out local memory: per warp, level k of an array of its 32 threads side
  by side. Thompson's column arrays (64 KB per column in FP64 at 48 levels) and MYNN's mass-flux plume arrays live
  there; the CPU keeps the local arrays.
- Grell-Freitas is compiled for a 48- and a 100-level tier, WSM6 for 48 / 64 / 128 levels (the tier of nz).
- Kain-Fritsch keeps its local arrays (29 KB FP64, 11 KB FP32 at 48 levels): in the work buffer its FP32 build
  (kf.cpp contracts multiply-adds) gave different P2 / P4 values. Also measured and dropped as slower: a k-major buffer
  (column stride) with Thompson's sedimentation arrays kept local (Thompson FP64 +2-4 %, P6 FP32 steps +7 %) and
  the warp layout with that split (Thompson FP32 +29 %).

Frames per thread (cuobjdump, 48-level tier, FP64 / FP32): Thompson 72128 / 34224 -> below YSU's 10768 / 4240 (YSU's
frame is now the stack limit of P6 and C01), MYNN 71312 / 35656 -> 38640 / 16104, Grell-Freitas 127536 / 63776 -> 68464 /
34240, WSM6 49216 / 24608 (128 levels) -> 18496 / 9248.

Base origin/main 2afa191 / port/F06 7fa5304, one G4 session, 1 h, no history output, each twice interleaved; peak
from nvidia-smi (200 ms), the stack limit and its reserve from `CPPWRF_MEMINFO`:

| Case | base s/h | new s/h | peak MiB base -> new | stack limit (reserve) | memory after set-up |
|---|---|---|---|---|---|
| P2 FP64 | 1.018, 1.018 | 1.029, 1.027 | 9085 -> 9085 | 29008 B, KF (7989 MiB) | 1308 |
| P5 FP64 | 1.446, 1.435 | 1.453, 1.446 | 19623 -> 6669 | 18496 B, WSM6 (5094) | 1542 |
| P6 FP64 | 1.703, 1.669 | 1.713, 1.709 | 21937 -> 5789 | 10768 B (2965) | 2690 |
| P8 FP64 | 1.247, 1.224 | 1.210, 1.208 | 36219 -> 19951 | 68464 B, GF (18854) | 1306 |
| C01 FP64 | 1.984, 1.936 | 1.991, 1.966 | 22009 -> 5861 | 10768 B (2965) | 2762 |
| MYNN FP64 | 2.550, 2.497 | 2.547, 2.507 | 23683 -> 15215 | 38640 B, MYNN (10641) | 4430 |
| P4 FP64 | 5.124, 5.120 | 5.139, 5.112 | 21193 -> 11135 | 29008 B (7989) | 3114 |
| N3 FP64 | 14.008, 13.808 | 14.025, 13.817 | 24075 -> 14017 | 29008 B (7989) | 5670 |
| P4 FP32 | 4.181, 4.176 | 4.189, 4.188 | 20743 -> 5715 | 10960 B (3018) | 2774 |
| P6 FP32 | 0.450, 0.444 | 0.477, 0.469 | 19449 -> 2955 | 4240 B (1168) | 1800 |
| N3 FP32 | 9.688, 9.697 | 9.727, 9.703 | 23177 -> 8149 | 10960 B (3018) | 4988 |

Kernel time per hour (`CPPWRF_KTIME=1`, base -> new): Thompson FP64 426 -> 412 ms (P6), 498 -> 509 (C01), 456 ->
457 (MYNN case); Thompson FP32 116 -> 134 ms (P6), 134 -> 158 (C01); MYNN 502 -> 498; Grell-Freitas 419 -> 369
(P8, the 48-level tier); KF and WSM6 unchanged within 1-2 %. The steps are 0-1 % slower in FP64 (P8 3 % faster),
P6 FP32 6 % slower (Thompson FP32 +16-18 %: its local arrays were cached better than the global buffer); the peak
memory falls by 9-17 GB on the G4 in every case with a Thompson, WSM6, MYNN or GF scheme.

Proof: every cuda_dp case against WRF at tolerance 0 (79 PASS, v05_files_links included); the 16 FP32 suite cases
with base and new builds (120 snapshots and 3 restart files byte-identical) and P5 / N3 FP32 with 10-minute
snapshots (21 byte-identical); the CPU suite cases of the touched schemes at tolerance 0.

Small GPUs: the reserve scales with the resident threads, so an RTX 2060 (6 GB, 30720 threads) reserves 321 MiB
for P4 / N3 in FP32 (KF's 10960 B; 64 KB before: 1.9 GiB). N3 FP32 needs about 5.0 GB after the set-up on the G4
(less on a smaller GPU's context) plus that reserve: it should fit 6 GB, with little margin (not run on such a GPU).
With F05's RRTMG batch merged, a 6 GB GPU's RRTMG scratch is 1/8 of the memory (0.75 GB, shared by the domains);
on the G4 F05 adds about 2.5 GB to the single-domain peaks above and 1.2 GB to N3. When memory is short the run
stops at the set-up or the first step with the message above instead of a failed cudaDeviceSetLimit.

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


## RRTMG batch restructuring (AR1 rank 3, F05; measured, Colab G4)

The device RRTMG batches (`RrtmgLwBatch`, `RrtmgSwBatch`) ran their stages in chunks of 2048 columns: 6 chunks for
the 10 591 columns of P4, each stage launch a grid of 32 blocks of 64 threads (32 of the G4's 188 SMs). Now:

- **All columns in one chunk** where the device memory allows (`rrtmg_chunk_columns`, phys/ra_rrtmg_common/
  batch_memory.cpp: the work arrays may take an eighth of the device memory and at most a quarter of what is free,
  at least 2048 columns; above that the fewest equal chunks). Per step 5 LW and 4 SW launches instead of 30 and 24.
  The per-column prep, which is latency-bound (the 2048-column launch and the 10 591-column launch take the same
  ~5 ms), runs once per step instead of 6 times.
- **Work arrays and tables shared**: one double and one int work buffer for every LW / SW batch of every domain
  (they run one after another on the one stream; the buffers grow to the largest request when a batch is
  constructed, outside the time loop), and one device copy per table set (LW, SW tables and constants) for all
  domains. The stages read the McICA cloud arrays only under the flags the same call set (as before, when a
  domain's last, smaller chunk reused the arrays of the full chunks with another stride).
- **LW setcoef per (column, layer)** (`KLwSetcoef`, the body of setcoef's layer loop as `rrtmg_lw_setcoef_layer`,
  laytrop counted per layer) after the per-column prep: KLwPrep 4.3 -> (chunking) 5.4 -> 4.3 ms per step in P4.

Same kernels and arithmetic per column (rrtmg_*_batch.cpp keep `--fmad=false` in both precisions); the CPU column
driver runs the same setcoef loop in order.

Stepping s/h, 1 h, no history output, base origin/main 2afa191 / port/F05 df400ae, built in the same session, each
twice, interleaved; graph nodes per graph launch (`CPPWRF_GRAPHS=stats`; N3: per parent step's graph launches);
peak device memory from nvidia-smi:

| Case | base s/h | new s/h | gain | graph nodes | peak MiB |
|---|---|---|---|---|---|
| P4 FP64 (radt 1, night) | 5.055, 5.074 | 2.963, 2.953 | 1.71x | 420 -> 376 | 21193 -> 23711 |
| D0 FP64 (radt 15) | 2.373, 2.370 | 2.256, 2.253 | 1.05x | 386.8 -> 383.9 | 23633 -> 26151 |
| N3 FP64 (radt 1, two domains) | 13.897, 13.870 | 9.764, 9.763 | 1.42x | 389.3 -> 367.3 | 24075 -> 25261 |
| N3 FP32 | 9.617, 9.620 | 4.956, 4.954 | 1.94x | 536.3 -> 514.3 | 23177 -> 24359 |
| V1 FP64 (radt 1, sunrise: RRTMG SW active) | 10.485, 10.478 | 5.603, 5.603 | 1.87x | 420 -> 376 | 21193 -> 23711 |

`CPPWRF_KTIME=1` (time per kernel type, 1 h = 60 steps; ms per step):

| Kernel | P4 base | P4 new | V1 base | V1 new |
|---|---|---|---|---|
| KLwPrep | 29.2 | 4.3 | 29.0 | 4.3 |
| KLwBand<60> | 19.9 | 10.7 | 18.9 | 10.4 |
| KLwLayer | 6.1 | 5.1 | 6.2 | 5.2 |
| KLwFinish | 0.86 | < 0.6 | 0.86 | < 0.6 |
| KSwPrep | - | - | 46.9 | 8.3 |
| KSwGpoint<48> | - | - | 42.2 | 35.0 |
| all kernels | 87.6 | 51.9 | 177.6 | 96.1 |

Memory: the shared work buffers hold all columns (P4: ~4.0 GB for the SW layout, the larger one) instead of the
two 2048-column buffers per domain (~1.5 GB): +2.5 GB for P4 and D0, +1.2 GB for N3 (whose two domains now share
them). What remains in the radiation is FP64 arithmetic at 1/64 rate (KLwBand, KLwLayer, KSwGpoint); the SW prep
still runs taumol_sw per column (8.3 ms per daytime step in V1), the next candidate for a per-layer split.

Proof: every cuda_dp case against WRF at tolerance 0 (SUITE PASS 79 FAIL 0 SKIP 0, v05_files_links.tar included),
the 16 FP32 cases with base and new builds (120 snapshots and 3 restart .nc files byte-identical), the CPU suite cases
of every radiation option and the smoke set at tolerance 0 (15 PASS).

## Single-GPU kernels (P01; measured, Colab G4)

Profiled on the G4 (nsys and ncu with `CPPWRF_GRAPHS=0`, `CPPWRF_KTIME=1`) on P4 (radt 1 and radt 15), P6
(Thompson), MYNN, N3 and L (a 360 x 270 domain made on the runtime with WPS 4.6.0 from GFS 0.25 deg, P4's physics,
radt 15), FP64 and FP32. The top kernels and what changed:

- **KNoahRs** (one thread for the whole domain: the RC of each land-ice point is the last RC that SFLX computed in
  WRF's i-j loop order; 10 ms per step in L): now per point (`KNoahRcLast`, one thread per row, finds the last SFLX
  point of each row; `KNoahRs` takes it from its own row or the nearest row before). Same values.
- **KUpdate / KTendU, V, W, T** (FP64 compute-bound at 81-87 % of the SM: every face flux of the scalar transport
  and of the rk_tendency advection was computed twice, once by each cell of the face): in FP64 builds the face
  fluxes are computed once per face into a scratch array (`KFaceX/Y/Z`, `KFace<KTend*, d>`; the same expressions)
  and both cells read them. FP32 builds keep the inline fluxes (those files build with `--fmad=true`, so the code
  stays textually unchanged there).
- **KLwBand** (one thread per (column, band), read-modify-write of the band sums; memory-bound: 35 % occupancy,
  ~2.9 GB of local-memory traffic per launch): now `KLwGpoint` (one thread per (column, g-point): the rtrnmc sweeps
  of one g-point as `rrtmg_lw_gpoint`, results in place in the work arrays) and `KLwLevel` (per (column, level): the
  band sums in g-point order, then the level fluxes). The CPU band loop calls the same g-point routine in the same
  order.
- **column_for blocks of one warp (32 x 1)** instead of 64 x 2 (`CPPWRF_COLUMN_BLOCK=XxY` overrides it; X a multiple
  of 32): the column physics (WSM6, KF, YSU, Thompson, MYNN) run at 5-8 % occupancy and are latency-bound, and
  smaller blocks spread the columns over more SMs (P4 radt 15: 1.74 -> 1.64 s/h FP64, 0.56 -> 0.49 FP32).

Tried and dropped: taumol per (layer, band) for the LW prep (KLwLayer 309 -> 324 ms per hour on P4: slower).

Stepping s/h, 1 h, no history output, base origin/main 1b4c654 / port/P01 98a2c03, built in the same session, each
twice, interleaved; peak device memory from nvidia-smi (`CPPWRF_MEMINFO=1`: stack limit 29008 B FP64, unchanged):

| Case | base s/h | new s/h | gain | peak MiB base -> new |
|---|---|---|---|---|
| P4 FP64 (radt 1) | 2.915, 2.897 | 2.568, 2.563 | 1.13x | 13653 -> 13671 |
| P4 FP32 | 1.741, 1.739 | 1.546, 1.546 | 1.13x | 8233 -> 8233 |
| P4 FP64, radt 15 | 1.737, 1.724 | 1.514, 1.499 | 1.15x | 13653 -> 13673 |
| P4 FP32, radt 15 | 0.558, 0.556 | 0.483, 0.482 | 1.15x | 8233 -> 8233 |
| P6 FP64 (Thompson) | 1.661, 1.649 | 1.417, 1.405 | 1.17x | 5789 -> 5799 |
| P6 FP32 | 0.461, 0.461 | 0.397, 0.396 | 1.16x | 2955 -> 2949 |
| MYNN FP64 | 2.289, 2.272 | 2.011, 1.994 | 1.14x | 17733 -> 17745 |
| MYNN FP32 | 0.725, 0.724 | 0.644, 0.642 | 1.13x | 10427 -> 10421 |
| N3 FP64 (two domains) | 9.575, 9.528 | 8.456, 8.447 | 1.13x | 15203 -> 15241 |
| N3 FP32 | 4.854, 4.843 | 4.393, 4.390 | 1.10x | 9331 -> 9333 |
| L FP64 (360 x 270 x 45, radt 15) | 9.974, 9.971 | 8.085, 8.082 | 1.23x | 29957 -> 30067 |
| L FP32 | 3.145, 3.142 | 2.642, 2.642 | 1.19x | 20791 -> 20791 |

`CPPWRF_KTIME=1`, FP64, ms per step (advection: KUpdate, KTend* and the new face kernels; LW: every KLw* kernel;
column physics: WSM6 / Thompson, KF, YSU):

| Kernels | P4 radt 15 base | new | P6 base | new | L base | new |
|---|---|---|---|---|---|---|
| advection | 8.00 | 5.95 | 9.52 | 7.04 | 63.2 | 42.5 |
| KNoahRs (+ KNoahRcLast) | 1.14 | 0.03 | 1.12 | 0.03 | 9.97 | 0.04 |
| column physics | 9.96 | 9.72 | 8.48 | 8.28 | 30.7 | 28.1 |
| all kernels | 32.2 | 28.6 | 30.9 | 27.0 | 173.1 | 137.3 |

On P4 (radt 1) the LW band stage went from KLwBand 10.7 + KLwLevel 0.26 ms to KLwGpoint 7.2 + KLwLevel 1.85 ms per
step. Memory:
the FP64 face scratch (three arrays of the domain with halos) adds ~20 MiB on P4 and ~110 MiB on L (0.4 %).

Proof (tolerance 0 throughout): history of P4, P4 radt 15, P6, MYNN, N3 and L every 20 min, base against new,
byte-identical in FP64 and FP32; every cuda_dp case against WRF; the 16 FP32 suite cases with base and new builds
byte-identical; the CPU suite cases of the changed code (dycore, radiation options, Noah, nests, ideal cases) against
WRF.

## Column physics in stages (P01b; measured, Colab G4)

Profile of main 97d1a48 (`CPPWRF_KTIME=1`, ms per hour of 60 steps): WSM6 1044 (L FP64), 300 (P4), 1153 (N3, two
domains); Thompson 417 FP64 / 131 FP32 (P6), 449 / 151 (MYNN case); MYNN PBL 490 / 129; YSU 326 and KF 321 (L
FP64). ncu showed the column kernels at 3-8 % achieved occupancy (one thread per column, 200+ registers).

What changed: the level loops of WSM6 and Thompson that read only their own level now run one thread per
(i, j, k); the vertical sweeps (sedimentation, the column reductions, the early return of a column without
microphysics) stay one thread per column, in WRF's order. Every value is the same statement on the same operands
(the loop bodies moved unchanged, `x[k]` on scalars or on the same column arrays).

- **WSM6, FP64 builds** (`wsm6_stages`, CPU and CUDA): Init (clamps, t, cpm, xl), A (qsat, rh, fall speeds), Sed
  (rain and snow + graupel sedimentation per column), B (melting, ice fall speed), Ice (ice sedimentation and the
  surface precipitation per column), C (freezing, warm and cold rain processes, saturation adjustment), per
  microphysics loop; the values carried between stages in a work array of 12 levels-sized fields (`Wsm6Work`, +390
  MiB on L, +40 MiB on P4). FP32 builds keep `Wsm6Kernel`: physics_coupling.cpp builds with `--fmad=true` there.
- **Thompson, FP32 and FP64** (thompson.cpp builds with `--fmad=false` in both): `ThompsonLevel0` (the first two
  level loops; no_micro and k_melting as level flags), `ThompsonStage1` (column: their reduction, twet, the early
  return), `ThompsonLevel` (the snow moments to the rain evaporation), `ThompsonStage3` (column: sedimentation,
  updates, totals, radar and effective radii) on the column arrays in the work buffer (`column_work_chunks`; a level
  kernel finds its column's arrays with `ColumnSlots(work, i, j, chunk, bx, by)`). The CPU runs the four parts per
  column.

Tried and dropped: the rain and the snow + graupel sedimentation of WSM6 on two threads per column (Sed 252 -> 297
ms per hour on L: slower).

Kernel time per hour (KTIME, base -> new): WSM6 L FP64 1044 -> 817 (StageC 291 and Sed 252 the largest: StageC is
FP64-bound at 86 % SM throughput, Sed is the sequential sweep), P4 300 -> 178; Thompson P6 FP64 417 -> 257, FP32
131 -> 78; MYNN case FP64 449 -> 288, FP32 151 -> 94.

Stepping s/h, 1 h, no history output, base origin/main 97d1a48 / port/P01b 4711d5e in one G4 session, each twice:

| Case | base s/h | new s/h | gain | peak MiB base -> new |
|---|---|---|---|---|
| P6 FP64 (Thompson) | 1.404, 1.403 | 1.260, 1.259 | 1.11x | 5799 -> 5801 |
| P6 FP32 | 0.394, 0.396 | 0.345, 0.344 | 1.15x | 2949 -> 2951 |
| MYNN FP64 | 1.992, 1.992 | 1.849, 1.849 | 1.08x | 17745 -> 17745 |
| MYNN FP32 | 0.640, 0.640 | 0.583, 0.584 | 1.10x | 10421 -> 10423 |
| P4 FP64 (WSM6) | 2.563, 2.563 | 2.449, 2.450 | 1.05x | 13671 -> 13713 |
| P4 FP32 | 1.547, 1.549 | 1.550, 1.547 | 1.00x | 8233 -> 8233 |
| N3 FP64 | 8.444, 8.459 | 8.026, 8.051 | 1.05x | 15241 -> 15327 |
| N3 FP32 | 4.388, 4.398 | 4.397, 4.399 | 1.00x | 9333 -> 9333 |
| L FP64 | 8.099, 8.090 | 7.874, 7.876 | 1.03x | 30067 -> 30455 |
| L FP32 | 2.640, 2.641 | 2.642, 2.640 | 1.00x | 20791 -> 20791 |

Memory grows at most 1.3 % (L FP64: the WSM6 work array). Proof (tolerance 0): history every 20 min of the ten
runs above, base against new, byte-identical; every cuda_dp case against WRF (79 PASS); the 16 FP32 suite cases
with base and new builds (120 snapshots, 3 restart files byte-identical); the CPU suite cases of WSM6, Thompson,
MYNN, Grell-Freitas (column_for_work) and nests against WRF and the WSM6, radar, Thompson, MYNN and suite
differential tests.

Open: MYNN PBL (472 ms per hour FP64, 124 FP32 on the MYNN case) is now the largest column kernel there; YSU and
KF (L FP64 ~320 ms each) build with `--fmad=true` in FP32, so a split there is FP64-only like WSM6's.

## NetCDF output: chunks compressed on host threads (IO1; measured, Colab G4)

WRF's NetCDF-4 history and restart files (external/io_netcdf `ext_ncd_write_field`) are chunked in halves of the two
horizontal dimensions and compressed with shuffle + deflate level 2. CPPWRF writes the same files with the same
NetCDF calls. HDF5's deflate filter compressed them one chunk after the other in the model's thread. On the large
domain L that was most of a run with history on (`CPPWRF_IO_TIMING=1` with HDF5's own filter, as in main d9c25e5;
G4, FP64):

| File (L, FP64) | values (device copies) | puts + sync (HDF5: shuffle, deflate, write) | total |
|---|---|---|---|
| wrfout, one frame | 0.33 s (0.22 s) | 7.8 s + 2.8 s | 11.1 s |
| wrfrst at 1 h | 0.82 s (0.54 s) | 22.1 s + 7.8 s | 31.4 s |

P4 FP64 takes 1.25 s per frame and 3.5 s for wrfrst; FP32 about half. The device-to-host copies, the
unstaggering and type conversion are 3 to 5 % of it, and the file system (Colab's local disk) is not the limit.

**Change** (src/wrf_io/nc_deflate.cpp, nc_out.cpp, stream_file.cpp):
- `StreamFile::write` gathers the values of all variables of a record first. For each deflated variable, `NcOut`
  queues every chunk on a pool of host threads. A pool thread cuts the chunk (an edge chunk is padded with the
  variable's fill value, as HDF5 pads it), shuffles it (H5Zshuffle.c) and calls `compress2` at the variable's level,
  the call of HDF5's filter (H5Zdeflate.c).
- HDF5's deflate filter is replaced through the library's `H5Z_register`, the function HDF5 registers its own filters
  with. The public `H5Zregister` refuses predefined filter ids. The new filter takes the precomputed output when the
  chunk HDF5 hands it has the same bytes (a sampled key, then the whole buffer compared). Otherwise it compresses the
  chunk itself, as HDF5's filter does. The class has the same id and name, so the filter pipeline messages keep their
  bytes too.
- The NetCDF / HDF5 calls and their order do not change, and only zlib runs on the pool threads. So every file has
  the same bytes, whether a chunk was precompressed or not. Without the symbol (another HDF5 build), HDF5's own
  filter stays.
- The host buffers of the field copies and of the staged records are kept from one field and record to the next.
  They are no longer new pages each time (L FP64: device copies 0.25 -> 0.08 s and staging 0.23 -> 0.12 s per frame).
- `CPPWRF_NC_THREADS` sets the pool size (default: the hardware threads, at most 64). 1 turns the precomputation off
  and 0 keeps HDF5's own filter. `CPPWRF_IO_TIMING=1` prints each record's times.

After the change (L FP64): 1.15 s per wrfout frame (values 0.22, staging 0.12, puts 0.63 of which 0.14 waiting for the
pool, sync 0.18) and 4.1 s for wrfrst. Every chunk (660 per frame, 1932 in wrfrst) came from the pool, and none was
compressed by the filter itself. 16 or 24 pool threads give the same times as 48.

Wall time of the whole process, 1 h with a history frame every 10 minutes (7 frames in one file) and wrfrst at 1 h,
base main d9c25e5 against port/IO1 in one G4 session, each run twice (CPU FP64: once, 48 OpenMP threads). "Off" is
history_interval = 1440: only the frame at t = 0 is written.

| Case | base wall (s) | new wall (s) | gain | history off: base / new (s) | output (on - off): base -> new (s) |
|---|---|---|---|---|---|
| L GPU FP64 | 123.0, 123.0 | 34.6, 31.7 | 3.7x | 22.0 / 16.5 | 101 -> 16.7 |
| L GPU FP32 | 58.1, 57.7 | 13.1, 13.6 | 4.3x | 9.8 / 7.0 | 48 -> 6.3 |
| P4 GPU FP64 | 15.9, 15.9 | 5.2, 5.4 | 3.0x | 4.4 / 3.8 | 11.5 -> 1.5 |
| P4 GPU FP32 | 8.7, 8.7 | 3.4, 3.4 | 2.5x | 2.9 / 2.6 | 5.8 -> 0.8 |
| L CPU FP64 | 169.3 | 80.7 | 2.1x | | |
| P4 CPU FP64 | 24.9 | 18.6 | 1.3x | | |

The stepping times are unchanged (L FP64 7.58 s, FP32 2.52 s). On L in FP32 the hour with 7 frames and a restart
file now takes 13 s instead of 58 s, against 10 s with history off in main. An earlier session (without the kept
buffers) gave L FP64 31.6 and 29.9 s, against 120.3 and 119.2 s in main. In the CPU build the device copies are
host copies with OpenMP, which run next to the pool's compression (P4: 0.5 s per frame); the CPU's stepping takes the
rest of the time.

**Proof.** Every output file is byte-identical to main's (`cmp`): wrfout, wrfrst and CPPWRF's restart_*.nc of L and
P4 in GPU FP64, GPU FP32 and CPU FP64 (history every 10 min, 1 h). `cppwrf_real`'s wrfinput / wrfbdy of L are
byte-identical too. `test_nc_deflate` (ctest) writes a file with float, double and int variables and odd dimensions
(edge chunks) three times: with the pool, with `CPPWRF_NC_THREADS=1` and with HDF5's own filter. The three files must
be byte-identical, and every chunk must come from the pool. It passes with HDF5 1.10.10 (Colab's Ubuntu) and 1.14.5
(Debian 13).

Open: the remaining output time (L FP64 1.15 s per frame) is HDF5's own work in the model's thread: the shuffle,
the copy of the precompressed chunk and the write. It could move to a background thread while the GPU steps, but
NetCDF is not thread-safe, and the model reads NetCDF files during the run (boundaries, wrfinput fields, sst_update).
