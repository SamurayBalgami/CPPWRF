# CPPWRF status — M1 (CPU), awaiting review

## Built
- M0: CMake/Ninja (CPU OpenMP + CUDA backends, FP32/FP64), single-source kernels, SoA fields,
  halo fill, WRF-named State/Tendencies/BaseState, WRF-style dry ideal init, I/O contract.
- M1: large step (rk_step_prep, 5th/3rd-order flux-form advection, PGF, pg_buoy_w, rhs_ph,
  constant-K diffusion), RK3 driver (WRF substep schedule, passes rk_step), diagnostics,
  acoustic small step. New WRF options: `dyn.h_adv_order` 5|6, `dyn.v_adv_order` 3|4.
- Run everything: `tools/run_m1.sh`. Tests: `ctest` 5/5 in build-cpu and build-cpu-dp
  (test_large_step adds diffusion-magnitude check: x 0.9992, z 0.9653 = 2nd-order FD factor).
- Determinism: two FP32 igw runs -> `cmp` of snap_000500.bin identical, diag.csv identical.

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
- GPU untested (driver). CPU perf not optimised (M2, after GPU is measurable).
- Run script lives in tools/run_m1.sh (scripts/ is 's).

## Disagreements
- None.

## Next steps (after review and GPU driver)
- CUDA runs in FP32/FP64, CPU-vs-GPU differences, ns/cell/step per backend; then M2 perf
  (profile first: halo fills, small-step scratch arrays, flux recomputation).
