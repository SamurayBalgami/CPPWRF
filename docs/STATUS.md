# CPPWRF: status and verification

CPPWRF reproduces WRF-ARW v4.8.0 statement by statement. This page lists what is implemented, how it
is verified against WRF, and the current results. All numbers come from `tests/verify/wrf_suite.py`
(145 cases and 27 restart / Straka checks, all passing) and `tests/verify/run_wrf_files_test.sh` (WRF output files).

## Implemented

| Component | WRF options |
|---|---|
| Dynamics | RK3 split-explicit, adaptive time step (use_adaptive_time_step, single-domain real-data runs: target_cfl, target_hcfl, max_step_increase_pct, starting / max / min_time_step, step_to_output_time, restarts), hybrid_opt 0/2, hypsometric_opt 1/2, map factors, Coriolis and curvature, terrain, specified (spec_zone 1, relax_zone 4), open, periodic and symmetric lateral boundaries (symmetric: no wind across the walls), w_damping, damp_opt 1, 2 and 3, h_adv_order 5/6, v_adv_order 3/4, positive-definite moisture, scalar and TKE advection (h_adv_order 5 and 6) |
| Diffusion | diff_opt 1 with km_opt 1 (constant K), 3 (3-D Smagorinsky) and 4 (2-D Smagorinsky); diff_opt 2 with km_opt 2 (1.5-order TKE, advected TKE); 6th-order diffusion (diff_6th_opt 1, 2) |
| Microphysics | Kessler (1), WSM6 (6, with effective radii), Thompson (8) |
| Cumulus | Kain-Fritsch (1), Grell-Freitas-Li (3, deep and shallow), Tiedtke (6), New Tiedtke (16) |
| PBL, surface layer, land surface | YSU (1), MYJ (2), MYNN-EDMF (5, icloud_bl 0 / 1, with the surface layers 1, 2 and 5; bl_mynn_tkeadvect, bl_mynn_cloudpdf -2..2, bl_mynn_mixlength 0..2, bl_mynn_closure 2.5 / 2.6 / 3.0, bl_mynn_edmf 0 / 1, bl_mynn_edmf_mom, bl_mynn_edmf_tke, bl_mynn_cloudmix, bl_mynn_mixqt, bl_mynn_ess, bl_mynn_mixscalars, bl_mynn_output 0 / 1, tke_budget 0 / 1), MYNN surface layer options sf_mynn_sfcflux_land 0..3 and sf_mynn_sfcflux_water 0..4, fog deposition (grav_settling 1, 2), revised MM5 surface layer (1), Eta surface layer (2), MYNN surface layer (5), Noah LSM (2, including glacial and sea-ice points; usemonalb and rdlai2d: background albedo and LAI from wrfinput), Noah-MP v5.1.2 (4, default options), fractional sea ice (fractional_seaice = 1), SST update (sst_update = 1, wrflowinp), gravity wave drag (gwd_opt = 1, GWDO) |
| Radiation | RRTMG LW and SW (4), CAM ozone and greenhouse gases, cal_cldfra1, McICA |
| Nesting | static one-way and two-way nests (interp_fcn_sint, bdy_interp1, copy_fcn feedback, sm121/smdsm smoothers, blend_terrain, adjust_tempqv) |
| Namelist | physics_suite CONUS and tropical (the -1 options), in tools/wrf2cppwrf.py and in the WRF I/O schema |
| I/O | reads wrfinput_d0N / wrfbdy_d01; writes WRF v4.8.0 wrfout_d0N and wrfrst_d0N files (WRF I/O schema, iofields_filename, history_outname / rst_outname, frames_per_outfile, use_netcdf_classic) and its own restart files; native binary output; cppwrf_real (real.exe) |
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
   wrfinput/wrfbdy files (119 x 89 x 44, 15 km, dt 60 s; nests 120 x 90 x 44 at 5 km). The CPU FP64 cases
   compare CPPWRF's wrfout files with WRF's at tolerance 0 (tests/verify/exact.py: NaN equals only NaN, -0
   differs from +0) in every variable the reference holds, every reference frame up to the case's time, types and
   attributes, and check that the configuration sets the options of the reference run's namelist.input; the GPU
   and FP32 cases compare fields by relative L2 over the domain. The comparators' own test
   (tests/verify/selftest.py: one ulp, NaN, -0, a missing variable or frame must fail) runs before the VM suite.
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
| A01 / A01V | P1 / V2 with usemonalb and rdlai2d (wrfinput ALBBCK and LAI; V2: snow-covered points start with ALBEDO = SNOALB) | 1 h |
| P2 | P1 + Kain-Fritsch | 1 h |
| P3 | P1 + RRTMG LW/SW (radt 1 min) | 3 steps |
| P4 | P2 + WSM6 + RRTMG (cloudy, effective radii) + diff_opt 1 / km_opt 4 | 1 h |
| KESSKF | P4 with Kessler instead of WSM6 (Kain-Fritsch rain in Kessler's split sedimentation steps) | 1 h |
| P5 | P1 + WSM6 + diff_opt 1 / km_opt 4 | 1 h |
| P6 | P5 with Thompson | 1 h |
| P7 | P2 with New Tiedtke | 1 h |
| P8 | P2 with Grell-Freitas-Li (deep and shallow) | 1 h |
| P9 | P1 with MYJ PBL and the Eta surface layer | 1 h |
| MYNN | P6 with the MYNN surface layer and the MYNN-EDMF PBL (icloud_bl 1) + RRTMG (radt 10 min) | 1 h |
| MYNNSF1 / MYNNSF2 | MYNN with the revised MM5 / Eta surface layer instead of the MYNN one | 1 h |
| MYNNX1 / MYNNX2 | MYNN with non-default MYNN options (X1: tkeadvect, cloudpdf 1, mixlength 2, closure 3.0, edmf_tke 1, mixqt 1, output 1, tke_budget 1, sfcflux land 2 / water 3, grav_settling 2; X2: edmf 0, cloudmix 0, cloudpdf 0, mixlength 0, closure 2.5, ess 0, mixscalars 1, sfcflux land 3 / water 4, grav_settling 1) | 1 h |
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
cases. The Thompson and MYNN references (P6, MYNN, MYNNSF1, MYNNSF2, MYNNX1, MYNNX2, mp8, V3, D0, D0V, SU1, SU2) come from WRF built with -r8 -double-8: with -r8 alone
gfortran makes Thompson's DOUBLE PRECISION variables 16-byte.

Bit identity holds for references made with gfortran 14.2 on the CPU family the references were made on
(Intel x86-64, glibc 2.41; another gfortran version vectorises other loops, e.g. in Thompson): WRF's
vectorised loops and CPPWRF call the same glibc libmvec functions, and glibc's `_ZGVbN2v_log` gives
different last bits on AMD EPYC (all scalar libm functions and the other libmvec functions tested agree).
WRF itself therefore differs between the two; on AMD, CPPWRF differs from Intel-made references by
about 1e-14 after one step. The suite runs on Intel machines (`tools/gcp/gcp_run.sh`).
The CPU result does not depend on the OpenMP thread count or on the Intel machine: P4, P6, P8, P9, MYNN,
P1 with Noah-MP and the two-way nest N3s (1 h, both domains) and em_quarter_ss (30 min) with 1, 16 and 44
threads are bit-identical to serial WRF on c3 (Sapphire Rapids) and n2 machines (`run.threads` in
tools/gcp/gcp_matrix.sh). The required build flags are in [USER_GUIDE.md](USER_GUIDE.md#required-flags-and-determinism).

Restarts reproduce an uninterrupted run bit for bit (P1, A01, P7, P8, P9, NMP1, MYNN, MYNNSF1, MYNNSF2, MYNNX1, MYNNX2, T1; P1 also in FP32 on the GPU), and a run
resumed from WRF's own wrfrst files continues as WRF does (also the two-way physics nest N3, both domains).
The restart variants of tests/verify/matrix_restart_io.txt (tools/gcp/gcp_matrix.sh) write restart files in the
middle of a run and resume both models from WRF's: the wrfrst files of the straight and the resumed runs (data and
alarm attributes) and the resumed history are identical to WRF's for P1, P4 (radt 10), Kessler, Thompson, New
Tiedtke, Tiedtke, Grell-Freitas, MYJ, MYNN, Noah-MP, G2, V2, diff_opt 2 / km_opt 2, diff_6th_opt, A01, em_quarter_ss
(also diff_opt 2), SU1 across the 6 h wrflowinp and boundary update, and history intervals that are not a divisor
of the restart time (after a restart WRF's input_wrf takes the seconds to the next history time as the history
interval, which the resumed runs reproduce); the same holds for the two-way nest N3s (P4 physics, also with
Kessler, Thompson, New Tiedtke, Tiedtke, Grell-Freitas, MYJ, MYNN, diff_opt 2 / km_opt 2, A01, feedback 0), the
sibling nests N3sib and the nest in a nest N3gc. Where WRF's resumed run differs from its straight
run (sst_update, the history quirk), CPPWRF's differs in the same way.

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
| R1D6 (10 min) | FP64 | 3.1e-14 | 7.7e-12 | 3.9e-16 | 2.9e-15 |
| W1 (75 min, dt 450 s) | FP64 | 2.2e-12 | 1.3e-11 | 9.2e-15 | |
| P1 / A01 | FP64 | 1.5e-9 | 2.0e-8 | 1.4e-10 | 8.9e-9 |
| P1T2 | FP64 | 1.0e-9 | 2.4e-8 | 2.1e-10 | 5.9e-9 |
| HY1 | FP64 | 1.6e-11 | 3.0e-10 | 2.8e-12 | 1.8e-10 |
| P3 (RRTMG every step) | FP64 | 1.1e-10 | 1.9e-8 | 1.0e-11 | 4.4e-11 |
| P8 | FP64 | 3.6e-10 (L4), 8.1e-10 (G4) | 2.0e-8 / 2.2e-8 | 1.2e-10 / 7.8e-10 | 5.4e-9 / 1.5e-7 |
| P7 | FP64 | 4.7e-7 | 1.1e-5 | 6.0e-8 | 1.2e-5 |
| P2 | FP64 | 5.2e-5 | 1.8e-3 | 9.4e-7 | 4.0e-5 |
| P4 | FP64 | 2.2e-4 | 3.6e-3 | 2.7e-6 | 2.1e-4 |
| P5 (WSM6) | FP64 | 8.4e-5 | 1.1e-3 | 2.6e-6 | 5.1e-4 |
| P6 (Thompson) | FP64 | 3.2e-5 | 9.0e-4 | 4.0e-7 | 8.5e-6 |
| MYNN | FP64 | 2.6e-4 | 5.8e-3 | 5.5e-6 | 1.1e-4 |
| MYNNSF1 / MYNNSF2 | FP64 | 2.5e-4 / 2.6e-4 | 5.3e-3 / 5.4e-3 | 5.3e-6 / 5.6e-6 | 2.9e-4 / 1.1e-4 |
| MYNNX1 / MYNNX2 (largest of L4 and G4) | FP64 | 3.2e-4 / 2.1e-4 | 8.6e-3 / 7.3e-3 | 9.0e-6 / 4.8e-6 | 1.4e-4 / 6.7e-5 |
| KESSKF | FP64 | 6.9e-5 (L4), 9.2e-5 (G4) | 1.2e-3 | 1.1e-6 / 1.2e-6 | 4.4e-5 / 6.2e-5 |
| A01V (Arctic, V2 domain) | FP64 | 7.8e-5 | 3.2e-3 | 2.4e-6 | 1.8e-4 |
| V1 / V2 (2 h) | FP64 | 1.8e-3 / 1.1e-4 | 9.4e-2 / 3.8e-3 | 1.6e-5 / 3.6e-6 | 4.7e-4 / 3.5e-4 |
| V3 (30 min) | FP64 | 2.5e-7 | 5.4e-6 | 7.7e-8 | 1.7e-6 |
| F1 / F2 | FP64 | 1.3e-4 / 2.1e-4 | 3.4e-3 / 4.6e-3 | 4.8e-6 / 3.7e-6 | 2.3e-4 / 7.9e-4 |
| SU2 (20 min) | FP64 | 1.3e-5 | 4.1e-4 | 4.7e-7 | 6.1e-5 |
| N1 (3 min) | FP64 | 1.8e-13 (d01), 2.1e-13 (d02) | 6.3e-12 / 1.0e-11 | 2.3e-16 / 3.3e-16 | 4.7e-15 / 5.1e-15 |
| N2, N2S1, N2S2 | FP64 | 4.8e-13 (d01), 2.4e-13 (d02) | 6.8e-12 / 9.0e-12 | 9.8e-16 | 4.0e-14 / 2.4e-14 |
| N3 (largest of L4 and G4) | FP64 | 1.8e-4 (d01), 2.4e-4 (d02) | 4.8e-3 / 5.2e-3 | 2.9e-6 / 2.2e-6 | 9.4e-5 / 1.0e-4 |
| R0 (dynamics) | FP32 vs WRF FP32 | 1.6e-4 | 4.9e-3 | 1.0e-6 | 2.0e-5 |
| P1 | FP32 vs WRF FP32 | 1.9e-4 | 3.9e-3 | 7.5e-7 | 4.9e-5 |
| P3 | FP32 vs WRF FP32 | 1.9e-4 | 4.0e-3 | 7.5e-7 | 5.1e-5 |
| P4 / P6 / P8 | FP32 vs WRF FP32 | 2.2e-4 to 8.5e-4 | up to 4.0e-2 | up to 1.5e-5 | up to 6.3e-4 |
| N1 (1 h) | FP32 vs WRF FP32 | 1.9e-4 (d01), 2.8e-4 (d02) | 3.6e-3 / 4.7e-3 | 6.5e-7 / 9.7e-7 | 2.0e-5 / 3.1e-5 |

On the GPU the elementary functions (exp, log, pow) come from the CUDA math library, which rounds
differently from glibc, so CUDA FP64 agrees with WRF to roundoff rather than bit for bit. Schemes
with thresholds (the Kain-Fritsch trigger, WSM6 and Thompson saturation adjustments) amplify that
roundoff over an hour, which is why P2 and P4 show larger CUDA FP64 differences than P1. The CPU
build is the reference for bit identity. Every CPU FP64 case of the suite has a CUDA twin (FP64, or FP32 against
WRF's single-precision run), with tolerances of about three times the differences above (Colab L4 and RTX PRO 6000
Blackwell "G4"). On the GPU a run restarted after its first step reproduces the uninterrupted run bit for bit
(restart_<case>_fp64_cuda for P1, A01, P7, P8, P9, NMP1, MYNN, MYNNSF1, MYNNSF2, MYNNX1, MYNNX2, T1), and the
wrfout / wrfrst files of a GPU run (P1, P4, MYNN, N2, also resumed from WRF's wrfrst) have WRF's variables,
dimensions, attributes, global attributes and times (run_wrf_files_test.sh with CPPWRF_FILES_NO_VALUES=1).
The GPU results depend on the GPU architecture only through the compiler: with the same nvcc 13.0 and libdevice,
the G4 (sm_120) build differs from the L4 (sm_89) build from the first step in two surface fields (ZNT, SFROFF of
the surface layer and Noah, compiled with FMA contraction), while the compute_89 PTX run on the G4 and builds with
--fmad=false give the L4's result bit for bit on both. P8's Grell-Freitas thresholds amplify that difference to
QV 1.5e-7 on the G4 (5.4e-9 on the L4); both runs are deterministic.

Ideal cases (CPPWRF's port of WRF's ideal.exe, `run.case = em_*`), CPU FP64 against WRF v4.8.0
real*8: em_hill2d_x (1 h), em_grav2d_x (Straka, 900 s), em_quarter_ss (Kessler supercell, 2 h, periodic and with open boundaries in x and y) and
em_squall2d_x / em_squall2d_y (open lateral boundaries in x / y, TKE diffusion, 10 min) are bit-identical
in every field; em_quarter_ss also with diff_opt 2 / km_opt 2 (and isfflx 0). Their wrfout files equal WRF's in every
variable and attribute (tests/verify/ideal_files_test.py).
Symmetric lateral boundaries (WRF symmetric_xs/xe/ys/ye, rigid walls) with soundings without wind are bit-identical
too, with the storm and the cold pools against the walls: em_quarter_ss on all four sides (21 x 21 cells, 1 h), the
Straka density current and em_squall2d_x (TKE) in x (40 min); in the variant matrix also x or y alone, diff_opt 2 /
km_opt 2, km_opt 3 / 4, 6th-order diffusion, h_adv_order 6, moist_adv_opt 0, em_hill2d_x and winds along the walls. A
wind across a symmetric wall is rejected: WRF keeps the wall value of the normal velocity while the sign of muu flips
there, and its solution breaks down within a minute.
On the GPU em_squall2d_x is at roundoff in FP64 (U 1.3e-13 relative L2 at 10 min) and
within WRF's own single/double difference in FP32 (U 6.8e-5; WRF real*8 vs. WRF FP32: 8.2e-5); the other ideal
cases (also symmetric, open and diff_opt 2) agree with WRF in CUDA FP64 to 1e-16 .. 5e-11 (W).

Namelist variants (tools/gcp/gcp_matrix.sh, 3-minute runs from the P1, P4, P6, P8 and P9 inputs, every wrfout
variable at tolerance 0): kfeta_trigger 3, ysu_topdown_pblmix 0, icloud 0, cldovrlp 1-5 with idcor 0 and 1,
o3input 0, ghg_input 0, use_mp_re 0, rdmaxalb false, ifsnow 0, kim_tofd, isfflx 0, Kessler on real data (also with
Kain-Fritsch: cudt, kfeta_trigger 3, dt 30 s and 90 s, radiation off, 4.5 h with restart files),
radt / bldt / cudt > dt, Grell-Freitas clos_choice 1, 4, 10 and ishallow 0, hybrid_opt (input-file check),
diff_6th_opt 1 and 2, moist/scalar_adv_opt 0, w_damping 0 and w_crit_cfl, damp_opt 0, dampcoef / zdamp,
epssm, time_step_sound, smdiv / emdiv, c_s, km_opt 1 / 3 / 4 with diff_opt 1 (km_opt 3 also on em_quarter_ss with
diff_opt 2), diff_opt 2 / km_opt 4 (also mix_full_fields), mix_isotropic 1, PD advection with h_adv_order 6,
moist_mix2_off, MYNN (P6 inputs) with RRTMG, icloud_bl 0, bldt 2, WSM6, Kain-Fritsch and restart_interval 1 (the restart alarm in the first step), the MYNN PBL on the revised MM5 and Eta surface layers (P6 inputs, 3 min: without radiation, RRTMG with icloud_bl 0, bldt 2, WSM6, Kain-Fritsch, Noah-MP, isfflx 0, iz0tlnd 1, isftcflx 1, wrfrst with restart_interval 1; 1 h with RRTMG and wrfrst at 30 min; two-way nest N3s on d01 and d02), every non-default MYNN option value with a port (P6 inputs, MYNN surface layer, RRTMG: 3 min and 1 h each, wrfrst at 30 min with bl_mynn_tkeadvect, closure 3.0, grav_settling 2 and the combined X1 / X2 options, bldt 2, and the N3s two-way nest with tkeadvect, edmf 0 / closure 3.0 / mixlength 2, the MYNN surface-layer forms and grav_settling), grav_settling 1 and 2 with YSU, usemonalb and rdlai2d (together and alone, also false) on P1, P4, P9, P1 with Noah-MP, V2 (snow, sea ice, glacier) and the N3s nest (both domains), 3 min and 1 h, with wrfrst files at 30 min (P1) and 15 min (V2), and on SU1 over 6 h 20 min across the 6 h wrflowinp update of LAI and ALBBCK, open-boundary 6th-order diffusion (em_squall2d_x, em_quarter_ss open), v_adv_order 4, and surface layer + LSM without a PBL scheme (bl_pbl_physics 0, with
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

sf_sfclay_physics 91, km_opt 5, the MYNN PBL with fractional_seaice = 1, the MYNN options bl_mynn_edmf 2, bl_mynn_edmf_dd, bl_mynn_mixnumcon, bl_mynn_mixaerosols 0, bl_mynn_output 2, other bl_mynn_closure values and the GFS surface-layer forms (sf_mynn_sfcflux_land 4, sf_mynn_sfcflux_water 5), Noah-MP options other than the defaults (and crop/irrigation/
tile drainage), the adaptive time step with nests, sst_update or ideal cases, moving nests, ndown.exe (one-way nesting from a parent run's wrfout), I/O forms other than NetCDF
(io_form 2) and auxiliary history / input streams other than sst_update's wrflowinp (WRF writes none by default),
and parallel (MPI) runs. Options without a port are rejected by `tools/wrf2cppwrf.py` and by the model with an
error. Noah-MP with nests is rejected because wrf.exe 4.8.0 aborts with it (heap corruption in the nest run's
first history output).
