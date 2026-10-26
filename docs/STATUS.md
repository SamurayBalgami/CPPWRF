# CPPWRF: status and verification

CPPWRF reproduces WRF-ARW v4.8.0 statement by statement. This page lists what is implemented, how it
is verified against WRF, and the current results. All numbers come from `tests/verify/wrf_suite.py`
(103 cases, all passing) and `tests/verify/run_wrf_files_test.sh` (WRF output files).

## Implemented

| Component | WRF options |
|---|---|
| Dynamics | RK3 split-explicit, hybrid_opt 0/2, hypsometric_opt 1/2, map factors, Coriolis and curvature, terrain, specified (spec_zone 1, relax_zone 4), open and periodic lateral boundaries, w_damping, damp_opt 1, 2 and 3, h_adv_order 5/6, v_adv_order 3/4, positive-definite moisture, scalar and TKE advection (h_adv_order 5 and 6) |
| Diffusion | diff_opt 1 with km_opt 1 (constant K), 3 (3-D Smagorinsky) and 4 (2-D Smagorinsky); diff_opt 2 with km_opt 2 (1.5-order TKE, advected TKE); 6th-order diffusion (diff_6th_opt 1, 2) |
| Microphysics | Kessler (1), WSM6 (6, with effective radii), Thompson (8) |
| Cumulus | Kain-Fritsch (1), Grell-Freitas-Li (3, deep and shallow), Tiedtke (6), New Tiedtke (16) |
| PBL, surface layer, land surface | YSU (1), MYJ (2), MYNN-EDMF (5, WRF's default bl_mynn_* options, icloud_bl 0 / 1), revised MM5 surface layer (1), Eta surface layer (2), MYNN surface layer (5), Noah LSM (2, including glacial and sea-ice points), Noah-MP v5.1.2 (4, default options), fractional sea ice (fractional_seaice = 1), SST update (sst_update = 1, wrflowinp), gravity wave drag (gwd_opt = 1, GWDO) |
| Radiation | RRTMG LW and SW (4), CAM ozone and greenhouse gases, cal_cldfra1, McICA |
| Nesting | static one-way and two-way nests (interp_fcn_sint, bdy_interp1, copy_fcn feedback, sm121/smdsm smoothers, blend_terrain, adjust_tempqv) |
| Namelist | physics_suite CONUS and tropical (the -1 options), in tools/wrf2cppwrf.py and in the WRF I/O schema |
| I/O | reads wrfinput_d0N / wrfbdy_d01; writes WRF v4.8.0 wrfout_d0N and wrfrst_d0N files (WRF I/O schema, iofields_filename) and its own restart files; native binary output; cppwrf_real (real.exe) |
| Ideal cases | WRF's ideal.exe for em_hill2d_x, em_quarter_ss, em_squall2d_x/y, em_grav2d_x; CPPWRF's own analytic tests (rest state, warm bubble, Straka, inertia-gravity wave, 2-D squall line) |

## Verification method

1. **Column and kernel tests** (`tests/phys`, `tests/dyn`, `tests/nest`): each scheme and
   kernel runs on the same inputs in C++ and in WRF v4.8.0's own Fortran (unmodified source or
   WRF's compiled objects, built with WRF's configure flags), and every output must be bit-identical
   in CPU FP64 (tolerance 0). This covers Kessler, WSM6, Thompson (with WRF's table files), YSU
   (ysu_topdown_pblmix 0 and 1), MYJ and the Eta surface layer, revised MM5 surface layer, Noah,
   Noah-MP, Kain-Fritsch, the diffusion and damping kernels (km_opt 1-4), the
   nest numerics, and the physics suite driver. Grell-Freitas-Li and New Tiedtke are checked by
   replaying the model's own column inputs through WRF's objects. Run them with
   `cmake -DCPPWRF_DOUBLE=ON -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_DIFF_TESTS=ON` and `ctest -L diff`,
   with `WRF_DIR` (a WRF -r8 build), `WRF_R8D8_DIR` (-r8 -double-8), `NETCDF` and
   `CPPWRF_WRF_RUNS` set.
2. **Model runs against WRF** (`tests/verify/wrf_suite.py`): WRF v4.8.0 and CPPWRF are run from the same
   wrfinput/wrfbdy files (119 x 89 x 44, 15 km, dt 60 s; nests 120 x 90 x 44 at 5 km) and the wrfout
   files are compared field by field (relative L2 over the domain).
3. **WRF output files** (`tests/verify/run_wrf_files_test.sh`, ctest `wrf_files_vs_wrf`): the
   wrfout and wrfrst files of the CPU FP64 runs, and the runs resumed from WRF's restart files,
   compared with WRF's files at tolerance 0 (values, variables, attributes, global attributes).
   With `CPPWRF_FILES_TRANSLATED=1` the CPPWRF configuration is generated from the WRF run's
   `namelist.input` by `tools/wrf2cppwrf.py`.
4. **Checkpoint dumps** (`tests/verify/wrf_dycore_dumps`): a WRF build that writes intermediate arrays
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
| P8 | P2 with Grell-Freitas-Li (deep and shallow) | 1 h |
| P9 | P1 with MYJ PBL and the Eta surface layer | 1 h |
| MYNN | P6 with the MYNN surface layer and the MYNN-EDMF PBL (icloud_bl 1) + RRTMG (radt 10 min) | 1 h |
| NMP1 | P1 with Noah-MP | 1 h |
| T1 | R1 with diff_opt 2 / km_opt 2 (TKE) | 5 steps |
| R1D6 | R1 with diff_opt 1 / km_opt 4, damp_opt 1, diff_6th_opt 2 | 10 min |
| N1 | one-way nest (both domains) | 3 parent steps |
| N2 | two-way nest (both domains) | 1 h |
| N2S1 / N2S2 | N2 with smooth_option 1 / 2 (smoothing of the fed-back fields, ht included) | 1 h |
| N3 | two-way nest with the P4 physics on both domains (physics feedback, nest ozone from the parent) | 1 h |
| V1 | tropics (Bay of Bengal), Mercator, 15 km, P4 physics | 6 h |
| V2 | Arctic (Svalbard), polar stereographic, 15 km, P4 physics (sea ice, glacier land ice) | 6 h |
| V3 | Alps, Lambert, 3 km, Thompson + MYJ + Eta + Noah-MP + RRTMG (radt 5), no cumulus | 2 h |
| T6 | P2 with Tiedtke (cu_physics 6) | 3 h |
| G2 | P2 with gravity wave drag (gwd_opt 1; ELVMAX = HGT + 3 VAR, WPS 4.6 writes none) | 3 h |
| F1 / F2 | V2 with fractional sea ice (real.exe with fractional_seaice 1), sfclayrev + YSU / Eta + MYJ | 3 h |
| D0 / D0V | WRF's default em_real namelist (physics_suite CONUS, fractional sea ice, gwd_opt 1, diff_opt 2 / km_opt 4, w_damping 0) on the P2 / V2 domains, configuration from tools/wrf2cppwrf.py | 6 h |
| SU1 / SU2 | D0V with sst_update 1 (wrflowinp from real.exe, read every 6 h); SU2: changed sea ice at the start (water/ice transitions) | 7 h / 20 min |
| HY1 | P1 with hypsometric_opt 1 (real.exe and wrf.exe) | 1 h |

V1-V3, T6, G2, F1, F2, D0, D0V, SU1 and SU2 are compared in every wrfout variable over the whole run; V1-V3 start from the
global GFS of 2026-10-12 00 UTC through WPS, so the domains, projections and surface types differ from the P
cases. The Thompson and MYNN references (P6, MYNN, mp8, V3, D0, D0V, SU1, SU2) come from WRF built with -r8 -double-8: with -r8 alone
gfortran makes Thompson's DOUBLE PRECISION variables 16-byte.

Bit identity holds for references made with gfortran 14.2 on the CPU family the references were made on
(Intel x86-64, glibc 2.41; another gfortran version vectorises other loops, e.g. in Thompson): WRF's
vectorised loops and CPPWRF call the same glibc libmvec functions, and glibc's `_ZGVbN2v_log` gives
different last bits on AMD EPYC (all scalar libm functions and the other libmvec functions tested agree).
WRF itself therefore differs between the two; on AMD, CPPWRF differs from Intel-made references by
about 1e-14 after one step. The suite runs on Intel machines (`tools/gcp/gcp_run.sh`).

Restarts reproduce an uninterrupted run bit for bit (P1, P7, P8, P9, NMP1, MYNN, T1; P1 also in FP32 on the GPU), and a run
resumed from WRF's own wrfrst files continues as WRF does (also the two-way physics nest N3, both domains).

Output files: the `wrfout_d0N` files of P1-P9, NMP1, MYNN, mp8, V1-V3, T6, G2, F1, F2, D0, D0V, SU1, HY1, R1, T1, N1, N2 and N3 (both
domains) and the `wrfrst_d0N` files of P1, P3, P4, P6, P7, P8, P9, NMP1, MYNN, mp8, N2, N3, T6, G2, F1, F2, D0, D0V and SU1 are
identical to WRF v4.8.0's. This covers
every variable, value, dimension, attribute and global attribute
(`tests/verify/run_wrf_files_test.sh`). Runs of these cases resumed from WRF's
`wrfrst` files, give the same files as WRF's resumed runs. The exceptions are values WRF leaves uninitialised:
RQCBLTEN / RQIBLTEN in runs without qc / qi, and QZ0 / USTM at fractional sea-ice points with the revised MM5
surface layer (restart fields). `cppwrf_real` output is identical to real.exe for 41 cases
([REAL.md](REAL.md)).

CUDA (relative L2 at 1 h unless noted):

| Case | Precision | U | W | THM | QV |
|---|---|---|---|---|---|
| R1 (5 steps) | FP64 | 2.6e-14 | 5.8e-12 | 3.7e-16 | 1.8e-15 |
| P1 | FP64 | 1.5e-9 | 2.0e-8 | 1.4e-10 | 8.9e-9 |
| P8 | FP64 | 3.7e-10 | 7.6e-9 | 9.8e-11 | 4.0e-9 |
| P7 | FP64 | 4.7e-7 | 1.1e-5 | 6.0e-8 | 1.2e-5 |
| P2 | FP64 | 5.2e-5 | 1.8e-3 | 9.4e-7 | 4.0e-5 |
| P4 | FP64 | 2.2e-4 | 3.6e-3 | 2.7e-6 | 2.1e-4 |
| P6 (Thompson) | FP64 | 7.6e-5 | 2.0e-3 | 5.0e-6 | 2.0e-5 |
| N2 | FP64 | 4.8e-13 (d01), 2.3e-13 (d02) | | 9.8e-16 | |
| R0 (dynamics) | FP32 vs WRF FP32 | 1.6e-4 | 4.9e-3 | 1.0e-6 | 2.0e-5 |
| P1 | FP32 vs WRF FP32 | 1.9e-4 | 3.9e-3 | 7.5e-7 | 4.9e-5 |
| P4 / P6 / P8 | FP32 vs WRF FP32 | 2.2e-4 to 8.5e-4 | up to 4.0e-2 | up to 1.5e-5 | up to 6.3e-4 |

On the GPU the elementary functions (exp, log, pow) come from the CUDA math library, which rounds
differently from glibc, so CUDA FP64 agrees with WRF to roundoff rather than bit for bit. Schemes
with thresholds (the Kain-Fritsch trigger, WSM6 and Thompson saturation adjustments) amplify that
roundoff over an hour, which is why P2 and P4 show larger CUDA FP64 differences than P1. The CPU
build is the reference for bit identity.

Ideal cases (CPPWRF's port of WRF's ideal.exe, `run.case = em_*`), CPU FP64 against WRF v4.8.0
real*8: em_hill2d_x (1 h), em_grav2d_x (Straka, 900 s), em_quarter_ss (Kessler supercell, 2 h, periodic and with open boundaries in x and y) and
em_squall2d_x / em_squall2d_y (open lateral boundaries in x / y, TKE diffusion, 10 min) are bit-identical
in every field; em_quarter_ss also with diff_opt 2 / km_opt 2 (and isfflx 0). Their wrfout files equal WRF's in every
variable and attribute (tests/verify/ideal_files_test.py). On the GPU em_squall2d_x is at roundoff in FP64 (U 1.3e-13 relative L2 at 10 min) and
within WRF's own single/double difference in FP32 (U 6.8e-5; WRF real*8 vs. WRF FP32: 8.2e-5).

Namelist variants (tools/gcp/gcp_matrix.sh, 3-minute runs from the P1, P4, P6, P8 and P9 inputs, every wrfout
variable at tolerance 0): kfeta_trigger 3, ysu_topdown_pblmix 0, icloud 0, cldovrlp 1-5 with idcor 0 and 1,
o3input 0, ghg_input 0, use_mp_re 0, rdmaxalb false, ifsnow 0, kim_tofd, isfflx 0, Kessler on real data,
radt / bldt / cudt > dt, Grell-Freitas clos_choice 1, 4, 10 and ishallow 0, hybrid_opt (input-file check),
diff_6th_opt 1 and 2, moist/scalar_adv_opt 0, w_damping 0 and w_crit_cfl, damp_opt 0, dampcoef / zdamp,
epssm, time_step_sound, smdiv / emdiv, c_s, km_opt 1 / 3 / 4 with diff_opt 1 (km_opt 3 also on em_quarter_ss with
diff_opt 2), diff_opt 2 / km_opt 4 (also mix_full_fields), mix_isotropic 1, PD advection with h_adv_order 6,
moist_mix2_off, MYNN (P6 inputs) with RRTMG, icloud_bl 0, bldt 2, WSM6, Kain-Fritsch and restart_interval 1 (the restart alarm in the first step), open-boundary 6th-order diffusion (em_squall2d_x, em_quarter_ss open), v_adv_order 4, and surface layer + LSM without a PBL scheme (bl_pbl_physics 0, with
diff_opt 1 or diff_opt 2 / km_opt 2) are bit-identical to WRF. Two-way nest N3 variants: smooth_option 0 / 1,
feedback 0, parent_time_step_ratio 2, per-domain cu_physics (cumulus off on d02), cudt, radt and bldt. New Tiedtke
(cu_physics 16, 30-minute runs) on the P1, P4 (also cudt 5), P6 and P9 inputs and on the N3s (also off on d02), N3sib
and N3gc nests, whose dx below 15 km takes the scheme's log scale factor, is bit-identical. (WRF leaves RQCBLTEN / RQIBLTEN of runs without qc / qi
and LH with isfflx 0 uninitialised; they are not compared.)

## Performance

NVIDIA RTX 2060 (6 GB), real-data domain 119 x 89 x 44, 1 h:

| Run | CPPWRF | WRF v4.6.1 |
|---|---|---|
| Dry dynamics, FP32 | 1.26 s (GPU) | 13.2 s (4 OpenMP threads) |
| Dry dynamics, FP64 | 4.3 s (GPU, 1/32-rate FP64) | 37.9 s (serial, real*8) |

Dry dynamics run at about 32 ns per cell and step in FP32 on the GPU, about 52 ns with WSM6.
RRTMG runs on the GPU as a batched pipeline (docs/PHYSICS_SUITE.md section 10).

## Not yet implemented

sf_sfclay_physics 91, km_opt 5, Noah-MP options other than the defaults (and crop/irrigation/
tile drainage), moving nests, and parallel (MPI) runs. Options without a port are rejected by
`tools/wrf2cppwrf.py` and by the model with an error.
