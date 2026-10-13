# CPPWRF status: M3 real-data dycore done (2026-10-13)

## Built
- M0: CMake/Ninja (CPU OpenMP + CUDA backends, FP32/FP64), single-source kernels, SoA fields,
  halo fill, WRF-named State/Tendencies/BaseState, WRF-style dry ideal init, I/O contract.
- M1: large step (rk_step_prep, 5th/3rd-order flux-form advection, PGF, pg_buoy_w, rhs_ph,
  constant-K diffusion), RK3 driver (WRF substep schedule, passes rk_step), diagnostics,
  acoustic small step. New WRF options: `dyn.h_adv_order` 5|6, `dyn.v_adv_order` 3|4.
- Run everything: `tools/run_m1.sh`. Tests: `ctest` 5/5 in build-cpu and build-cpu-dp
  (test_large_step adds diffusion-magnitude check: x 0.9992, z 0.9653 = 2nd-order FD factor).
- Determinism: two FP32 igw runs -> `cmp` of snap_000500.bin identical; diag.csv identical in columns 1-10
  (column 11, wall_s, is wall-clock time and differs by design).

## Results (tools/run_m1.sh, CPU = 12 threads, all verifier checks PASS)
| case | prec | backend | key metrics | stepping wall | ns/cell/step |
|---|---|---|---|---|---|
| rest_state 1 h (24x24x40) | FP32 | cpu | max w 8.24e-6 m/s, max th' 3.8e-6 K | 19.2 s | 2316 |
| rest_state | FP64 | cpu | max w 1.88e-13 m/s, max th' 7.1e-15 K | 17.6 s | 2120 |
| straka 100 m, 900 s | FP32 | cpu | min th' -9.456 K, front 15.346 km, sym 1.5e-5 | 39.7 s | 1347 |
| straka | FP64 | cpu | min th' -9.456 K, front 15.346 km, sym 7e-14 | 44.0 s | 1493 |
| warm_bubble 1000 s | FP32 | cpu | max w 15.018 m/s, sym 2.0e-4 | 28.5 s | 1427 |
| warm_bubble | FP64 | cpu | max w 15.018 m/s, sym 2.5e-12 | 25.3 s | 1265 |
| igw 3000 s | FP32 | cpu | th' [-0.00149, 0.00273] K, centroid err 769 m | 7.6 s | 5068* |
| igw | FP64 | cpu | th' [-0.00146, 0.00277] K, centroid err 1212 m | 2.0 s | 1312 |
(*) igw has 3000 cells: timing is overhead-dominated and noisy. No perf work done yet.
CUDA: builds with nvcc 12.4, not runnable (driver not loaded) -> no GPU numbers, no CPU-vs-GPU diff.

## M3: real-data dynamics identical to WRF (2026-10-13)
Same wrfinput/wrfbdy as WRF v4.6.1 (119x89x44, 15 km Lambert over Turkey, GFS
2026-10-12 00Z, dt 60 s, physics/diffusion/damping off, hybrid_opt 0, hypsometric_opt 2,
specified boundaries, qv advected with PD). tests/verify/compare_real.py, relative L2:

| run | precision | after | U | W | PH | THM | P |
|---|---|---|---|---|---|---|---|
| vs WRF real*8 (R1_r8) | CPU FP64 | 5 steps | 2.6e-14 | 5.8e-12 | 5.4e-16 | 3.7e-16 | 7.2e-15 |
| vs WRF real*8 (R0_r8) | CUDA FP64 | 1 h | 5.2e-14 | 8.1e-12 | 8.1e-16 | 1.1e-15 | 8.2e-15 |
| vs WRF FP32 (R0) | CUDA FP32 | 1 h | 2.8e-5 | 4.4e-3 | 4.5e-7 | 6.1e-7 | 4.5e-6 |

Ideal cases vs WRF: Straka (diff_opt=1) theta' RMS 7.4e-6 K at 900 s; em_hill2d_x w RMS
8e-5 m/s at 1 h. Wall time of the 1 h real run (WRF timed exclusive):
CPPWRF RTX 2060 FP32 1.26 s vs WRF FP32 4 OpenMP threads 13.2 s (~10x); FP64 4.3 s vs
WRF real*8 serial 37.9 s (~9x, consumer GPU at 1/32 FP64 rate). GPU dycore 3D FP32:
24.3 ns/cell/step (200x200x50).

Real data with WSM6 (mp_physics 6, run_mp6_r8): step 1 matches WRF real*8 to roundoff; from
step 2 ~600 ice cells differ (qi 2e-7, THM 5e-5). CPPWRF's own CPU FP64 vs CUDA FP64 diverge
the same way (qi 2.08e-7 in 604 cells, THM 5.0e-5 at step 2), and standalone WSM6 on WRF's
state is bit-exact with the Fortran: roundoff amplified by WSM6 thresholds near ice
saturation, not a coupling error. calc_cq is precomputed in WRF's summation order (KCq).

## M2 done (2026-10-13)
Items 1-2 were verified and committed (4f58fcc) and
wrote items 3-4. Commands below were run with OMP_NUM_THREADS=4 (the machine is shared).
- **1. Scalar transport** (scalars.*: advect_scalar / advect_scalar_pd with ru_m, rv_m, ww_m,
  rk_update_scalar): tests in test_scalars (constant tracer stays 1, PD positivity, mass).
- **2. Moist dynamics** (theta_m in t_2, calc_cq cqu/cqv/cqw in the PGF/advance_uv/advance_w/
  calc_coef_w, moist calc_p_rho_phi, moist ideal init): test_moist.
- **3. Physics coupling** (physics_coupling.*, docs/PHYSICS_INTERFACE.md): WRF
  moist_physics_prep_em -> microphysics -> moist_physics_finish_em -> calc_p_rho_phi after the
  RK loop; h_diabatic in rk_addtend_dry and removed in the last small_step_finish.
  Config: `phys.mp_physics` (0, 1 = Kessler), `phys.mp_tend_lim`. rainnc/rainncv in snapshots.
- **4. Moist case** `squall2d` (Weisman-Klemp 1982 sounding, 2 K bubble, Kessler):
  configs/squall2d.ini (120x40, dx 1 km, 1 h, no wind -> must stay symmetric).
- ctest 7/7 in CPU FP32, CPU FP64, CUDA FP32, CUDA FP64 (build-a0-*).

| squall2d check | result |
|---|---|
| deep convection | max w 31 m/s at 24 min, cold pool -10 K, upper warm anomaly +15 K |
| symmetry, FP64 CPU | 7e-13 (10 min) -> 1.5e-10 (1 h) |
| symmetry, FP32 CPU | 3e-4 (10 min) -> 0.08-0.15 (1 h): roundoff amplified by convection (FP64 shows the code is symmetric) |
| GPU vs CPU, FP64 | rel L2 <= 1.6e-11 at 1 h (w, th', qr, rainnc) |
| GPU vs CPU, FP32 | 1e-4 (10 min) -> 5e-2 (1 h): same chaotic growth |
| positivity | qc, qr >= 0 everywhere, every snapshot (PD advection) |
| water budget, no microphysics | total water constant to 1e-16 (transport is conservative) |
| water budget, Kessler | (water + rainnc) - W0 = -2.7e-4 W0 at 1 h, 0.75% of the rain. WRF Kessler's sedimentation uses the mass-level spacing (rdzk), not the layer thickness, so it is not exactly conservative; our Kessler is bit-exact with WRF's |
| timing (shared, loaded machine) | CUDA FP32 590 ns/cell/step vs CPU-4 FP32 ~9000: not a benchmark |

## REVIEW item 1: Straka front ~2% short (investigation, FP64 CPU, front / min th')
| variant | front km | min th' K |
|---|---|---|
| baseline (dt 1, n_sound 6, epssm .1, smdiv .1, emdiv .01, uniform_z) | 15.346 | -9.456 |
| dt 0.5 / n_sound 12 / epssm 0.02 / smdiv=emdiv=0 | 15.345–15.347 | -9.454..-9.508 |
| WRF linear eta levels | 15.316 | -9.461 |
| h_adv_order 6 / v_adv_order 4 / both (no upwind dissipation) | 15.345 / 15.341 / 15.343 | -9.507 / -9.596 / -9.608 |
| kvdif = 0 (khdif 75) | 15.618 | -12.85 |
Findings: (a) perturbation is -15 K in T (theta' = dT/exner(pb); initial min th' -16.62 K). ✔
(b) theta diffusion acts on theta - t_init (= theta - 300 here). ✔ (c) diffusion operator
magnitude verified against analytic -mu K k^2 (unit test); WRF's vertical form equals
mu K d2/dz2, i.e. Straka's constant-K Laplacian. ✔ (d) eta levels: -30 m.
Time step, acoustic parameters and advection dissipation do not move the front (<6 m); the
upwind dissipation explains ~0.15 K of the min th' gap. The front is controlled by vertical
diffusion (75 -> 0 moves it +272 m). I found no bug. Remaining candidates (not tested): front
measured at z = 50 m (lowest level) vs. Straka's surface value; WRF's first-stage-only
(forward-in-time) diffusion; mass-weighted diffusion along eta vs. along z. Verdict: within
the published model spread; no change made.

## Departures from WRF/TN
- Vertical diffusion BCs: zero flux at bottom/top (u, v, theta), w = 0 at both ends (WRF
  diff_opt=1: no vertical diffusion in the lowest layer). u diffusion has no u_base.
- No Coriolis/curvature (all M1 cases f = 0; igw without rotation). Eta default uniform_z.
- hypsometric_opt = 1 in base state/init/calc_p_rho_phi (WRF v4 default 2); keep until terrain.

## Open issues
- Root of the repo collects Fortran/test junk from the physics agents (*.mod, *.o, input_kf.bin,
  output_*.bin, test_kind*): they must build in their own directories.
- Machine oversubscription: every OpenMP process uses 12 threads by default. Set
  OMP_NUM_THREADS (2-4) and OMP_WAIT_POLICY=passive for every agent's runs.
- GPU untested (driver). CPU perf not optimised (M2, after GPU is measurable).
- Run script lives in tools/run_m1.sh (scripts/ is 's).

## Disagreements
- None.

## Next steps (after review and GPU driver)
- CUDA runs in FP32/FP64, CPU-vs-GPU differences, ns/cell/step per backend; then M2 perf
  (profile first: halo fills, small-step scratch arrays, flux recomputation).
