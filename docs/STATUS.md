# CPPWRF status

## Core types committed
- `src/dyn/state.hpp`: `State`, `Tendencies`, `BaseState` with WRF Registry names and the
  index conventions (C++ k = Fortran k - 1) documented at the top of the file.
- `src/core/halo.hpp`: `fill_halo(Field3|Field2, Grid[, Parity px, Parity py])`.
  Periodic and symmetric; U odd in x, V odd in y by default; 2D fields even.
- `src/core/field.hpp`: `View3/View2/View1` (state, tendencies), `BView3/BView2` (base state,
  si = sj = 0 when flat), `Array1` (1D vertical, host mirror + device copy).
- `src/core/backend.hpp`: `parallel_for(Range3D, k)`, `column_for(Range2D, k)`, `Reducer`.
- Interface notes for `SmallStep::run_stage`:
  * The RK3 driver copies `*_2 -> *_1` (u, v, w, t, ph, mu) at the start of each time step,
    before stage 1. SmallStep must not rely on doing that itself.
  * Before each stage the driver computes (rk_step_prep): `mut/muu/muv` from `mu_2`
    (stage values), `ru, rv, rw, ww`, `alt`, `php`, `p`, `al`; halos filled.
  * 2D runs (ny == 1): j-stride is 0, so y-differences vanish; driver skips v kernels.
  * `rdy = 0` for 2D runs.

## Environment
- g++ 14.2, cmake 3.31, ninja, numpy 2.2.4. nvcc 12.4 present; NVIDIA driver not loaded
  (`nvidia-smi`: "couldn't communicate with the NVIDIA driver"), so CUDA builds compile
  but cannot run yet.
