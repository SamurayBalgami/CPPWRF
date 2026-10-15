# CPPWRF: status and verification

CPPWRF reproduces WRF-ARW v4.6.1 statement by statement. This page lists what is implemented, how it
is verified against WRF, and the current results. All numbers come from `tests/verify/wrf_suite.py`
(41 cases, all passing).

## Implemented

| Component | WRF options |
|---|---|
| Dynamics | RK3 split-explicit, hybrid_opt 0/2, hypsometric_opt 1/2, map factors, Coriolis and curvature, terrain, specified lateral boundaries (spec_zone 1, relax_zone 4), w_damping, damp_opt 3, h_adv_order 5/6, v_adv_order 3/4, positive-definite moisture and scalar advection |
| Diffusion | diff_opt 1 with km_opt 1 (constant K) and 4 (2-D Smagorinsky), 6th-order diffusion |
| Microphysics | Kessler (1), WSM6 (6, with effective radii), Thompson (8) |
| Cumulus | Kain-Fritsch (1), Grell-Freitas (3, deep and shallow), New Tiedtke (16) |
| PBL, surface layer, land surface | YSU (1), revised MM5 surface layer (1), Noah LSM (2, including glacial and sea-ice points) |
| Radiation | RRTMG LW and SW (4), CAM ozone and greenhouse gases, cal_cldfra1, McICA |
| Nesting | static one-way and two-way nests (interp_fcn_sint, bdy_interp1, copy_fcn feedback, sm121/smdsm smoothers, blend_terrain, adjust_tempqv) |
| I/O | reads wrfinput_d0N / wrfbdy_d01, writes wrfout_d0N (NetCDF) and restart files; native binary output |
| Ideal cases | rest state, warm bubble, Straka density current, inertia-gravity wave, 2-D squall line, 2-D hill flow, quarter-circle supercell |

## Verification method

1. **Column and kernel tests** (`tests/phys`, `tests/dyn`, `tests/nest`): each scheme and kernel is
   driven with identical inputs in C++ and in the WRF Fortran (compiled with WRF's own flags) and the
   outputs are compared. Most of these tests require bit-for-bit agreement.
2. **Model runs against WRF** (`tests/verify/wrf_suite.py`): WRF v4.6.1 and CPPWRF are run from the same
   wrfinput/wrfbdy files (119 x 89 x 44, 15 km, dt 60 s; nests 120 x 90 x 44 at 5 km) and the wrfout
   files are compared field by field (relative L2 over the domain).
3. **Checkpoint dumps** (`tests/verify/wrf_dycore_dumps`): a WRF build that writes intermediate arrays
   inside solve_em and the physics drivers, used to locate the first differing operation.

Bit identity on the CPU requires reproducing WRF's arithmetic exactly: the order of every sum and
product, constant folding done by the Fortran compiler, the halo refresh schedule, and the vector math
functions (glibc libmvec) that gfortran uses in vectorised loops. CPPWRF calls the same libmvec
functions for the same loop iterations (`cppwrf::vm` in src/core/backend.hpp). The physics reference
runs use WRF built with `-fdefault-real-8 -fdefault-double-8`, so that DOUBLE PRECISION variables
(Thompson) are 8-byte as in CPPWRF.

## Results

CPU, FP64 against WRF real*8 (tolerance 0: every compared value identical):

| Case | Configuration | Duration |
|---|---|---|
| R1 | dynamics only | 5 steps |
| W1 | active w_damping (dt 450 s) | 10 steps |
| P1 | sfclayrev + Noah + YSU | 1 h |
| P2 | P1 + Kain-Fritsch | 1 h |
| P3 | P1 + RRTMG LW/SW (radt 1 min) | 3 steps |
| P4 | P2 + WSM6 + RRTMG (cloudy, effective radii) + diff_opt 1 / km_opt 4 | 1 h |
| P5 | P1 + WSM6 + diff_opt 1 / km_opt 4 | 1 h |
| P6 | P5 with Thompson | 1 h |
| P7 | P2 with New Tiedtke | 1 h |
| P8 | P2 with Grell-Freitas (deep and shallow) | 1 h |
| N1 | one-way nest (both domains) | 3 parent steps |
| N2 | two-way nest (both domains) | 1 h |

Restarts reproduce an uninterrupted run bit for bit (P1, P7, P8).

CUDA (relative L2 at 1 h unless noted):

| Case | Precision | U | W | THM | QV |
|---|---|---|---|---|---|
| R1 (5 steps) | FP64 | 2.6e-14 | 5.8e-12 | 3.7e-16 | 1.8e-15 |
| P1 | FP64 | 1.5e-9 | 2.0e-8 | 1.4e-10 | 8.9e-9 |
| P8 | FP64 | 3.7e-10 | 7.6e-9 | 9.8e-11 | 4.0e-9 |
| P7 | FP64 | 4.7e-7 | 1.1e-5 | 6.0e-8 | 1.2e-5 |
| P2 | FP64 | 5.2e-5 | 1.8e-3 | 9.4e-7 | 4.0e-5 |
| P4 | FP64 | 2.2e-4 | 3.6e-3 | 2.7e-6 | 2.1e-4 |
| N2 | FP64 | 4.8e-13 (d01), 2.3e-13 (d02) | | 9.8e-16 | |
| R0 (dynamics) | FP32 vs WRF FP32 | 1.6e-4 | 4.9e-3 | 1.0e-6 | 2.0e-5 |
| P1 | FP32 vs WRF FP32 | 1.9e-4 | 3.9e-3 | 7.5e-7 | 4.9e-5 |
| P4 / P6 / P8 | FP32 vs WRF FP32 | 2.2e-4 to 8.5e-4 | up to 4.0e-2 | up to 1.5e-5 | up to 6.3e-4 |

On the GPU the elementary functions (exp, log, pow) come from the CUDA math library, which rounds
differently from glibc, so CUDA FP64 agrees with WRF to roundoff rather than bit for bit. Schemes
with thresholds (the Kain-Fritsch trigger, WSM6 and Thompson saturation adjustments) amplify that
roundoff over an hour, which is why P2 and P4 show larger CUDA FP64 differences than P1. The CPU
build is the reference for bit identity.

Ideal cases: Straka density current against WRF em_grav2d_x, theta' RMS difference 7.4e-6 K at 900 s;
2-D hill flow against em_hill2d_x, w RMS difference 8e-5 m/s at 1 h.

## Performance

NVIDIA RTX 2060 (6 GB), real-data domain 119 x 89 x 44, 1 h:

| Run | CPPWRF | WRF v4.6.1 |
|---|---|---|
| Dry dynamics, FP32 | 1.26 s (GPU) | 13.2 s (4 OpenMP threads) |
| Dry dynamics, FP64 | 4.3 s (GPU, 1/32-rate FP64) | 37.9 s (serial, real*8) |

Dry dynamics run at about 32 ns per cell and step in FP32 on the GPU, about 52 ns with WSM6.
RRTMG runs on the GPU as a batched pipeline (docs/PHYSICS_SUITE.md section 10).

## Not yet implemented

MYJ PBL and the eta surface layer, sf_sfclay_physics 91, TKE closures (km_opt 2/3), diff_opt 2,
damp_opt 1/2, moving nests, a C++ real.exe (src/init_real is partial), and parallel (MPI) runs.
