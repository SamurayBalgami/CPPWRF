# CPPWRF status

## M0 — done (this commit)
- Build: CMake+Ninja; `CPPWRF_ENABLE_CUDA` (auto), `CPPWRF_DOUBLE`, sm_75. Build dirs used:
  `build-cpu` (FP32), `build-cpu-dp` (FP64), `build-cuda` (FP32, compiles; cannot run, see Env).
- Core: backend (`parallel_for`, `column_for`, deterministic `Reducer`), SoA fields (halo 3,
  pitch multiple of 32, j-stride 0 for ny=1), base-state fields stored 1D when flat
  (`BView3` with si=sj=0), halo fill (periodic/symmetric, parity-aware).
- `State`/`Tendencies`/`BaseState` with WRF Registry names (`src/dyn/state.hpp`).
- Cases: `rest_state`, `straka`, `warm_bubble`, `igw`; dry ideal init as WRF
  module_initialize_ideal (base state = sounding, mu'=0, p'=0, ph' hydrostatic).
- I/O contract via `cppwrf_io`: run.json, diag.csv, snap_*.bin/json.
- Large step (M1 groundwork, untested in time integration): rk_step_prep + rk_tendency
  (`src/dyn/large_step.cpp`), RK3 driver (`src/dyn/rk3.cpp`), diagnostics.
- Tests (`ctest --test-dir build-cpu`: 5/5 pass; FP64 likewise): `test_core` (layout,
  staggering, halo periodic/symmetric/corners, 2D y-periodicity, levels, Reducer),
  `test_large_step` (5th-order advection: FP64 order = 5.00; rest state -> all tendencies
  exactly 0; tendencies bitwise identical over repeated evaluations), io tests.
- `./build-cpu/cppwrf configs/rest_state.ini --set run.run_time=0` + verifier: all PASS.

## Interface notes for SmallStep::run_stage
- Driver copies `*_2 -> *_1` (u,v,w,t,ph,mu) at the start of each step (before stage 1).
- Before each stage the driver runs rk_step_prep: `mut/muu/muv` = mub + mu_2 (stage),
  `ru, rv, rw, ww, alt, php`; halos of u_2,v_2,w_2,t_2,ph_2,mu_2,p,al and of all derived
  fields are filled. SmallStep may overwrite mut/muu/muv (prep recomputes them).
- Tendencies as WRF: `ph_tend` includes the `+mut*g*w` term (rhs_ph, non_hydrostatic);
  `rw_tend` includes pg_buoy_w; `ru_tend/rv_tend` include horizontal_pressure_gradient;
  `*_tendf` (diffusion, stage 1 only) are already added into `*_tend`. `mu_tend` = 0.
- `top_lid = true` is used in the large-step PGF (dpn at the top extrapolated, cfn/cfn1).
- Stage schedule: dt/3 (1 step of dt/3), dt/2 (n/2 of dt/n), dt (n of dt/n);
  n = time_step_sound or the WRF auto rule. The driver needs `src/dyn/small_step.{hpp,cpp}`;
  CMake defines `CPPWRF_HAVE_SMALL_STEP` when small_step.cpp exists.

## Departures from WRF/TN (with reasons)
- Vertical diffusion BCs: zero flux at bottom/top for u, v, theta (Straka free-slip,
  insulating) and w = 0 at both ends. WRF diff_opt=1 instead sets the lowest-layer flux
  equal to the one above (no vertical diffusion in layer 1, meant for PBL schemes).
- Theta diffusion acts on theta - t_init (as WRF *_3dmp); u diffusion has no u_base.
- No Coriolis/curvature yet (all M1 cases use f = 0; igw is run without rotation).
- Eta levels default `uniform_z` (equal base-state dz) instead of WRF's linear eta, so that
  dx = dz = 100 m for Straka is honoured (max deviation 0.002 m).

## Open requests / notes
- `scripts/` belongs to , so the M1 run script is `tools/run_m1.sh` (not
  scripts/run_m1.sh as in the original prompt).
- M1 time integration is blocked on `src/dyn/small_step.*`.

## Environment
- g++ 14.2, cmake 3.31, ninja, numpy 2.2.4. nvcc 12.4 present (CUDA build compiles with
  no warnings), but `nvidia-smi` fails ("couldn't communicate with the NVIDIA driver"):
  driver not loaded (reboot likely needed). CUDA binaries exit with "no CUDA-capable device".

## Disagreements
- None so far.
