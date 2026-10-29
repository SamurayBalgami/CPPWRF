# CPPWRF User Guide

Welcome to **CPPWRF**, a modern, high-performance C++20 reimplementation of the Weather Research and Forecasting (WRF) Advanced Research WRF (ARW) dynamical core and physics suite.

This guide provides instructions for building CPPWRF on CPU and GPU architectures, running idealized and real-data cases, configuring model options via INI configuration files, translating WRF `namelist.input` files, interpreting output datasets, and executing the independent verification suite.

---

## Table of Contents
1. [Architecture & Design Principles](#1-architecture--design-principles)
2. [Building CPPWRF](#2-building-cppwrf)
   - [Prerequisites](#prerequisites)
   - [CPU Build (OpenMP, FP32 / FP64)](#cpu-build-openmp-fp32--fp64)
   - [CUDA GPU Build](#cuda-gpu-build)
   - [Build Flags & CMake Options](#build-flags--cmake-options)
3. [Running Simulations](#3-running-simulations)
   - [Command-Line Invocation](#command-line-invocation)
   - [Translating WRF Namelists (`wrf2cppwrf.py`)](#translating-wrf-namelists-wrf2cppwrfpy)
   - [Idealized Cases](#idealized-cases)
   - [Real-Data Simulations (WPS + real.exe)](#real-data-simulations-wps--realexe)
4. [Configuration Reference (INI Keys)](#4-configuration-reference-ini-keys)
   - [Summary Table](#configuration-keys-table)
5. [Outputs & Data Formats](#5-outputs--data-formats)
   - [Native Format (`run.json`, `diag.csv`, `snap_*.bin`)](#native-format)
   - [WRF NetCDF Output (`wrfout_d01_*`)](#wrf-netcdf-output)
   - [Restart Files & Resuming Simulations](#restart-files--resuming-simulations)
6. [Verification & Validation](#6-verification--validation)
   - [Independent Verification Harness (`verify.py`)](#independent-verification-harness-verifypy)
   - [Real-Data WRF Output Comparison (`compare_real.py`)](#real-data-wrf-output-comparison-compare_realpy)
   - [Unit & Differential Physics Tests](#unit--differential-physics-tests)
7. [Performance](#7-performance)

---

## 1. Architecture & Design Principles

CPPWRF is engineered to match WRF ARW (v4.8.0) numerics and equations discretely:

- **Equation Set:** Fully compressible, non-hydrostatic Euler equations in flux form with hybrid terrain-following vertical coordinates ($\eta$).
- **Time Integration:** 3rd-order Runge-Kutta (RK3) with split-explicit acoustic substeps.

- **Spatial Discretization:** Arakawa C-grid staggering; 5th- or 6th-order horizontal advection; 3rd- or 4th-order vertical advection.

- **Thermodynamic Consistency:** Diagnostic pressure and inverse density using moist potential temperature $\theta_m = \theta(1 + (R_v/R_d)q_v)$.

- **Single-Source Kernels:** Dynamics and physics kernels run identically on CPU (OpenMP parallel loops) and GPU (CUDA device kernels) without code duplication.

- **Strict IEEE Floating-Point:** Built with `-ffp-contract=off` (no fused multiply-add contraction) and without fast-math, so that the CPU FP64 results are bit-identical to WRF's and do not depend on the thread count ([Required Flags and Determinism](#required-flags-and-determinism)).

---

## 2. Building CPPWRF

### Prerequisites
- Modern C++20 compiler (`g++` >= 11 or `clang++` >= 14)
- CMake >= 3.24
- OpenMP runtime library
- NetCDF-C development headers (`libnetcdf-dev`)
- *(Optional)* NVIDIA CUDA Toolkit >= 12.0 for GPU offload

### CPU Build (OpenMP, FP32 / FP64)
To build standard single-precision (float32) for CPU:
```bash
cmake -B build-cpu -S . -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=OFF
cmake --build build-cpu -j4
```

To build double-precision (float64, `real = double`):
```bash
cmake -B build-cpu-dp -S . -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_DOUBLE=ON
cmake --build build-cpu-dp -j4
```

### CUDA GPU Build
To build the GPU-accelerated backend:
```bash
cmake -B build-cuda -S . -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=native
cmake --build build-cuda -j4
```

### Build Flags & CMake Options
| CMake Option | Default | Description |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | Compiler optimization level (`Release`, `Debug`, `RelWithDebInfo`) |
| `CPPWRF_DOUBLE` | `OFF` | Use double-precision floats (`real = double`) across dycore and physics |
| `CPPWRF_ENABLE_CUDA` | Auto | Enable CUDA compilation for NVIDIA GPUs |
| `CMAKE_CUDA_ARCHITECTURES` | `native` | Target GPU architecture (e.g. `80;90` for Ampere/Hopper) |
| `CPPWRF_NATIVE` | `ON` | Release builds: optimize for the host CPU (`-march=native`); `OFF` for portable binaries |
| `CPPWRF_DIFF_TESTS` | `OFF` | Register WRF Fortran differential test suites with `ctest` |

### Required Flags and Determinism
CPU FP64 builds give WRF v4.8.0 real*8's results bit for bit only with these flags; the CMake files set or
enforce them:

- **`-ffp-contract=off`** on every C++ source, in every build type and after `CMAKE_CXX_FLAGS` (and
  `-Xcompiler=-ffp-contract=off` for nvcc's host code). g++'s default is `-ffp-contract=fast`: with an FMA
  target (`-march=native`, `-march=skylake`) it fuses `a*b + c` into one rounding, which WRF's gfortran build
  (x86-64 without `-march`, so no FMA instructions) never does.
- **No fast-math:** `-Ofast`, `-ffast-math`, `-funsafe-math-optimizations`, `-fassociative-math`,
  `-freciprocal-math`, `-ffinite-math-only` and `-ffp-contract=fast|on` in `CMAKE_CXX_FLAGS` are a configure
  error. They reorder sums and products and let g++ call glibc's vector math functions in other loops than
  WRF does. In CUDA builds `--use_fast_math`, `-ftz=true`, `-prec-div=false` and `-prec-sqrt=false` in
  `CMAKE_CUDA_FLAGS` are an error too.
- **Optimisation level and instruction set** do not change the values (without fast-math g++ keeps every
  operation's order, and vectorised code rounds as scalar code): Release uses `-O3` and, with
  `CPPWRF_NATIVE=ON`, `-march=native`. The verification VMs build with `-DCPPWRF_NATIVE=OFF
  -DCMAKE_CXX_FLAGS="-march=skylake -mtune=skylake"`.
- **The WRF reference** must be WRF's own GNU serial build: `configure` option for gfortran (serial),
  `FCOPTIM = -O2 -ftree-vectorize -funroll-loops`, `-fdefault-real-8` (`-fdefault-real-8 -fdefault-double-8`
  for Thompson and MYNN), no `-march`, no OpenMP, with gfortran 14.2 and glibc 2.41 (Debian 13). Another
  gfortran version vectorises other loops (glibc libmvec), and `-march=native` lets gfortran contract FMAs.
- **Machine:** Intel x86-64 (the suite runs on Google Cloud c3 and n2 machines and on the reference
  workstation); on AMD EPYC glibc's vector `log` rounds differently, so WRF and CPPWRF differ from
  Intel-made references there.
- **Threads:** the CPU result does not depend on `OMP_NUM_THREADS`. Every kernel writes its own grid
  points, the radiation columns are independent (McICA seeds from each column's pressures), and the solution
  has no floating-point sums over the domain: the w_damping CFL maxima of the adaptive time step and the
  `diag.csv` sums use `Reducer`, which sums each level in a fixed order
  and then the levels in order, whatever the thread count.
- **CUDA:** nvcc's defaults (`-prec-div=true`, `-prec-sqrt=true`, `-ftz=false`); device code contracts
  FMAs (`--fmad=true`) except RRTMG's batch kernels, Thompson and the physics drivers (`--fmad=false`). GPU
  results agree with WRF to roundoff, not bit for bit (the CUDA math library rounds differently from glibc).

---

## 3. Running Simulations

### Command-Line Invocation
Simulations are launched by providing an INI configuration file to the executable:
```bash
./build-cpu/cppwrf configs/real_r0.ini
```

Individual configuration parameters can be overridden dynamically from the command line without modifying the `.ini` file:
```bash
# Override total simulation time to 2 minutes (120 s)
./build-cpu/cppwrf configs/real_r0.ini --set run.run_time=120

# Multiple overrides
./build-cpu/cppwrf configs/real_p4.ini --set output.format=both --set output.out_interval=5
```

### Translating WRF Namelists (`wrf2cppwrf.py`)
CPPWRF includes a standalone translator that converts WRF `namelist.input` files into valid CPPWRF `.ini` configurations:
```bash
# Basic conversion
python3 tools/wrf2cppwrf.py /path/to/wrf/namelist.input -o run.ini

# Real-data conversion with specified directory and NetCDF wrfout format
python3 tools/wrf2cppwrf.py /path/to/wrf/run/namelist.input -o run.ini --format wrfout
```

Key translator features:
- **f90nml-Free Parser:** Pure Python tokenizer supporting standard Fortran namelist syntax, repetition counts (`3*0`), exponents (`1.0d-4`), and inline comments.
- **Complete Option Coverage:**
  - **Time & I/O:** Translates `run_days/hours/minutes/seconds` (or `start_*`/`end_*` timestamps), `time_step` + `time_step_fract_num/den`, `time_step_sound`, `use_adaptive_time_step` with its &domains options, `history_interval_s` / `history_interval`, and `restart` / `restart_interval`.
  - **Dynamics (`[dyn]`):** Translates `diff_opt` (`0..2`), `km_opt` (`1..4`), `c_s`, `c_k`, `mix_isotropic`, `mix_upper_bound`, `mix_full_fields`, `tke_drag_coefficient`, `tke_heat_flux`, `tke_upper_bound`, `moist_mix2_off`, `scalar_mix2_off`, `tke_mix2_off`, `moist_mix6_off`, `scalar_mix6_off`, `tke_mix6_off`, `diff_6th_opt` (`0..2`), `diff_6th_factor`, `w_damping` (`0, 1`), `w_crit_cfl`, `damp_opt` (`0..3`), `zdamp`, `dampcoef`, `khdif`, `kvdif`, `moist_adv_opt` (`0, 1`), `scalar_adv_opt` (`0, 1`), `tke_adv_opt` (`0, 1`), `prandtl` ($\to$ `pr_inv`), `top_lid`, `h_mom_adv_order` ($\to$ `h_adv_order`), `v_mom_adv_order` ($\to$ `v_adv_order`), `epssm`, `smdiv`, `emdiv`, `etac`, `hybrid_opt` (`0..2`), and `hypsometric_opt` (`1, 2`).
  - **Physics (`[phys]` & `[moist]`):** Translates `mp_physics` (`0, 1, 6, 8`), `mp_tend_lim`, `hail_opt` (`0, 1`), `use_mp_re` (`0, 1`), `sf_sfclay_physics` (`0, 1, 2, 5`), `sf_surface_physics` (`0, 2, 4`), `num_soil_layers` (`4`), `opt_thcnd` (`1, 2`; Noah), `bl_pbl_physics` (`0, 1, 2, 5`), `cu_physics` (`0, 1, 3, 16`), `ishallow` (`0, 1`), `clos_choice` (`0..16`), `kfeta_trigger` (`1, 3`), `convtrans_avglen_m`, `ra_lw_physics` (`0, 4`), `ra_sw_physics` (`0, 4`), `radt`, `bldt`, `cudt`, `isfflx` (`0, 1`), `ysu_topdown_pblmix` (`0, 1`), `kim_tofd`, `tofd_factor`, `rdmaxalb`, `usemonalb` (every LSM: landuse_init keeps the wrfinput `ALBBCK`), `rdlai2d` (Noah; an urban point with `VEGFRA = 0` and `LAI <= 0` stops the run: WRF's result there is undefined), `icloud` (`0, 1`), `cldovrlp` (`1..5`), `idcor` (`0, 1`), `o3input` (`0, 2`), `ghg_input` (`0, 1`), `aer_opt` (`0`), and all 27 `&noah_mp` options when `sf_surface_physics = 4`.
  - **Nested Domains (`max_dom > 1`):** Translates multi-domain real-data namelists into `[domains]` (`max_dom`, `feedback`, `smooth_option`, `save_topo_from_real`, `blend_width`, `spec_bdy_width`, and `d02` placement keys) and per-child `[d02]`..`[d0N]` sections (`wrfinput`, `dt`, `out_interval`, `diag_interval`, `parent_id`, `i_parent_start`, `j_parent_start`, `parent_grid_ratio`, `parent_time_step_ratio`, and per-domain dynamics overrides).
- **Strict Validation & Rejection of Unsupported Options:**
  - **Rejected with Error (exit code 1):** Hydrostatic mode (`non_hydrostatic = .false.`), non-RK3 integration (`rk_ord != 3`), mismatched scalar vs. momentum advection orders (`h_sca_adv_order != h_mom_adv_order` or `v_sca_adv_order != v_mom_adv_order`), unsupported advection/diffusion/coordinate options (`moist_adv_opt`/`scalar_adv_opt`/`tke_adv_opt` outside `0..1`, `diff_6th_opt` outside `0..2`, `hybrid_opt` outside `0..2`, `hypsometric_opt` outside `1..2`), unpaired RRTMG radiation (`ra_lw_physics` vs. `ra_sw_physics`), `ishallow = 1` without `cu_physics = 3`, `clos_choice` outside `0..16`, `kfeta_trigger` outside `{1, 3}`, `num_soil_layers != 4` when an LSM is active, `cldovrlp` outside `1..5`, `idcor` outside `{0, 1}`, `o3input` outside `{0, 2}`, `ghg_input` outside `{0, 1}`, `aer_opt != 0`, unsupported `&noah_mp` option values, unimplemented physics switches (`iau`, `dyn_lightning_option`, `sf_urban_physics`, `cu_rad_feedback`, `shcu_physics`, `sf_lake_physics`, `topo_wind`, `slope_rad`, `topo_shading`, `swint_opt`, `seaice_thickness_opt`, `seaice_snowdepth_opt`, `ra_call_offset`, `shalwater_z0`, `scm_force_flux`), and invalid or unsupported nesting configurations (nested ideal cases, moving nests, vertical refinement, `feedback` outside `{0, 1}`, `smooth_option` outside `{0, 1, 2}`, invalid `parent_id` / `i_parent_start` / `j_parent_start` / `parent_grid_ratio` / `parent_time_step_ratio`, child dimensions not divisible by `parent_grid_ratio` or exceeding parent bounds, mismatched `e_vert` or `dx`/`dy`, `input_from_file = .false.`, or per-domain `[phys]` overrides in inline nested configs).
  - **Warnings & Notices:** Reports warnings on unimplemented scheme IDs (`--strict` promotes warnings to errors) and informational notices for ignored WPS/quilt metadata.

### Idealized Cases
WRF's own idealized cases are initialized by a port of WRF v4.8.0's `ideal.exe`
(`dyn_em/module_initialize_ideal.F`, then what `start_em` does with its output), from the case's
WRF `input_sounding`. In CPU FP64 they give WRF's (real*8) results bit for bit:
- **`em_hill2d_x`:** flow over a 2-D bell-shaped hill (WRF `ideal_case = 1`; `configs/wrf_hill2d_x.ini`).
- **`em_quarter_ss`:** 3-D supercell with Kessler microphysics (`ideal_case = 2`; `configs/wrf_quarter_ss.ini`).
- **`em_squall2d_x`, `em_squall2d_y`:** 2-D squall line with open lateral boundaries (`ideal_case = 4, 5`; `configs/wrf_squall2d_x.ini`).
- **`em_grav2d_x`:** Straka density current (`ideal_case = 6`; `configs/wrf_grav2d_x.ini`).

The namelist maps to `[domain] nx = e_we - 1`, `ny = e_sn - 1` (`ny = 1` for the x-z cases,
WRF `e_sn = 3`), `nz = e_vert - 1`, `dx`, `dy`, `ztop`, `bc_x`/`bc_y` (`&bdy_control`),
`[dyn] hybrid_opt`, `etac` and `[case] sounding_file`. As in WRF, ideal runs always use
`hypsometric_opt = 1`.

CPPWRF also has analytic test cases of its own (not WRF initializations):
- **`rest_state`:** Stably stratified atmosphere at rest testing hydrostatic balance.
- **`warm_bubble`:** Bryan & Fritsch (2002) dry rising thermal bubble in neutral atmosphere.
- **`straka`:** Straka et al. (1993) cold density current benchmark.
- **`igw`:** Skamarock & Klemp (1994) inertia-gravity wave.
- **`squall2d`:** Weisman & Klemp (1982) squall line with Kessler or WSM6 microphysics.

Run an idealized case:
```bash
./build-cpu/cppwrf configs/wrf_grav2d_x.ini
```

### Real-Data Simulations (WPS + real.exe)
Real-data simulations use initial conditions (`wrfinput_d01`) and lateral boundary updates (`wrfbdy_d01`) generated by the standard WRF Preprocessing System (WPS) and `real.exe`:

1. Run `geogrid.exe`, `ungrib.exe`, and `metgrid.exe` in WPS.
2. Run `real.exe` in WRF, or CPPWRF's `build/src/init_real/cppwrf_real` (identical output, see
   [REAL.md](REAL.md)), to generate `wrfinput_d01` and `wrfbdy_d01`.
3. In your `.ini` file:
   ```ini
   [run]
   case = real
   run_time = 3600

   [real]
   dir = /path/to/wrf_real/run_dir
   wrfinput = wrfinput_d01
   wrfbdy = wrfbdy_d01

   [time]
   dt = 60

   [output]
   format = wrfout
   out_interval = 10
   ```

4. Run `./build-cpu/cppwrf my_real_run.ini`.

---

## 4. Configuration Reference (INI Keys)

### Configuration Keys Table
Every option supported by the CPPWRF runtime is cataloged below, along with its type, default, WRF equivalent, and code source.

| Section | Key | Type | Default | Description | WRF Equivalent |
|---|---|---|---|---|---|
| `[run]` | `case` | string | *required* | Simulation case (`real`, `em_hill2d_x`, `em_quarter_ss`, `em_squall2d_x`, `em_squall2d_y`, `em_grav2d_x`, `straka`, `squall2d`, `warm_bubble`, `rest_state`, `igw`) | `ideal_case` |
| `[run]` | `run_time` | double | `0.0` | Total simulation duration in seconds | `run_days/hours/minutes/seconds` |
| `[run]` | `restart_file` | string | `""` | Path to NetCDF restart file (`restart_*.nc`) to resume execution | `restart` / restart filename |
| `[run]` | `start_date` | string | `"2000-01-01_00:00:00"` | UTC start timestamp string | `start_year/month/day/hour/etc` |
| `[real]` | `dir` | string | *required (real)* | Directory containing `wrfinput_d01` and `wrfbdy_d01` | WRF run directory |
| `[real]` | `wrfinput` | string | `"wrfinput_d01"` | Name of initial condition NetCDF file | `wrfinput_d01` |
| `[real]` | `wrfbdy` | string | `"wrfbdy_d01"` | Name of boundary condition NetCDF file | `wrfbdy_d01` |
| `[real]` | `bdy_interval` | double | `0.0` (auto) | Boundary condition update interval in seconds (0 = auto-detected from file) | `interval_seconds` |
| `[domain]` | `nx`, `ny`, `nz` | integer | *required (ideal)* | Grid dimensions in mass-point coordinates | `e_we - 1`, `e_sn - 1`, `e_vert - 1` |
| `[domain]` | `dx`, `dy` | double | *required (ideal)* | Horizontal grid spacing in meters | `dx`, `dy` |
| `[domain]` | `ztop` | double | *required (ideal)* | Model top elevation in meters | `ztop` |
| `[domain]` | `p_top` | double | Sounding $p(z_{top})$ | Pressure at top boundary in Pa | `p_top_requested` |
| `[domain]` | `bc_x`, `bc_y` | string | `"periodic"` | Lateral boundary condition (`periodic`, `symmetric`, `specified`, `open`) | `periodic_x`, `symmetric_xs`, `open_xs`/`open_xe`, etc. |
| `[domain]` | `eta` | string | `"uniform_z"` | Vertical coordinate distribution (`uniform_z`, `uniform_eta`, `list`) | `eta_levels` |
| `[domain]` | `eta_levels` | float list | `[]` | Explicit comma-separated list of full $\eta$ levels | `eta_levels` |
| `[domain]` | `map_factor` | double | `1.0` | Constant horizontal map projection scale factor | Map scale factor ($m$) |
| `[domain]` | `f0`, `e0` | double | `0.0` | Coriolis parameters $f$ and $e$ in $s^{-1}$ ($f$-plane) | Coriolis parameters |
| `[case]` | `sounding_file` | string | `"input_sounding"` | WRF `input_sounding` of an `em_*` case | `input_sounding` |
| `[case]` | `hm` | double | `100.0` | Hill height in m (`em_hill2d_x`; WRF hard-codes 100) | `hm` |
| `[case]` | `xa` | double | `5.0` | Hill half-width in grid points (`em_hill2d_x`; WRF hard-codes 5) | `xa` |
| `[time]` | `dt` | double | *required* | Large time step in seconds ($\Delta t$) | `time_step` |
| `[time]` | `time_step_sound` | integer | `0` (auto) | Acoustic small steps per RK stage (0 = auto-computed CFL rule) | `time_step_sound` |
| `[time]` | `adaptive` | bool | `false` | Adaptive time step (real-data, single domain; `dt` is then replaced by the starting step) | `use_adaptive_time_step` |
| `[time]` | `target_cfl`, `target_hcfl`, `max_step_increase_pct` | double, double, integer | `1.2`, `0.84`, `5` | Adaptive step: target vertical / horizontal CFL and the largest increase per step | same names (`&domains`) |
| `[time]` | `starting_time_step`, `max_time_step`, `min_time_step` (and `_den`) | integer | `-1` (4, 8, 3 × min(dx,dy)/1000 s) | Adaptive step limits in seconds (a fraction with `_den` > 0) | same names (`&domains`) |
| `[time]` | `step_to_output_time` | bool | `true` | Adaptive step: land on the history times | `step_to_output_time` |
| `[time]` | `history_interval_s`, `interval_seconds` | integer | `3600`, `21600` | Adaptive step: the history interval (seconds) and the boundary interval it steps to | `history_interval*`, `interval_seconds` |
| `[dyn]` | `khdif`, `kvdif` | double | `0.0` | Constant horizontal/vertical eddy diffusion in $m^2/s$ | `khdif`, `kvdif` |
| `[dyn]` | `diff_opt` | integer | `0` | Diffusion option (0 = off, 1 = coordinate-surface diffusion, 2 = physical-space diffusion) | `diff_opt` |
| `[dyn]` | `km_opt` | integer | `1` | Eddy viscosity closure (1 = constant, 2 = 1.5-order TKE, 3 = 3D Smagorinsky, 4 = 2D Smagorinsky) | `km_opt` |
| `[dyn]` | `c_s`, `c_k` | double | `0.25`, `0.15` | Smagorinsky and TKE closure constants | `c_s`, `c_k` |
| `[dyn]` | `mix_isotropic` | integer | `0` | Isotropic mixing length flag (0 = anisotropic, 1 = isotropic) | `mix_isotropic` |
| `[dyn]` | `mix_upper_bound` | double | `0.1` | Upper bound on eddy diffusivityCourant number | `mix_upper_bound` |
| `[dyn]` | `mix_full_fields` | boolean | `false` | Apply diffusion to full fields instead of perturbation fields | `mix_full_fields` |
| `[dyn]` | `tke_drag_coefficient`, `tke_heat_flux`, `tke_upper_bound` | double | `0.0`, `0.0`, `1000.0` | Surface drag, heat flux, and upper bound for prognostic TKE (`km_opt = 2`) | `tke_drag_coefficient`, `tke_heat_flux`, `tke_upper_bound` |
| `[dyn]` | `moist_mix2_off`, `scalar_mix2_off`, `tke_mix2_off` | boolean | `false` | Disable 2nd-order diffusion on moisture, scalars, or TKE | `moist_mix2_off`, `scalar_mix2_off`, `tke_mix2_off` |
| `[dyn]` | `moist_mix6_off`, `scalar_mix6_off`, `tke_mix6_off` | boolean | `false` | Disable 6th-order diffusion on moisture, scalars, or TKE | `moist_mix6_off`, `scalar_mix6_off`, `tke_mix6_off` |
| `[dyn]` | `diff_6th_opt` | integer | `0` | 6th-order spatial diffusion (0 = off, 1 = standard, 2 = monotonic) | `diff_6th_opt` |
| `[dyn]` | `diff_6th_factor` | double | `0.12` | 6th-order numerical diffusion coefficient | `diff_6th_factor` |
| `[dyn]` | `pr_inv` | double | `3.0` | Inverse turbulent Prandtl number ($Pr^{-1} = K_\theta / K_m$) | $1 / \text{prandtl}$ |
| `[dyn]` | `top_lid` | boolean | `false` | Rigid lid at top boundary (extrapolates $p'$) | `top_lid` |
| `[dyn]` | `h_adv_order` | integer | `5` | Horizontal advection order (5 = upwind, 6 = centered) | `h_mom_adv_order` |
| `[dyn]` | `v_adv_order` | integer | `3` | Vertical advection order (3 = upwind, 4 = centered) | `v_mom_adv_order` |
| `[dyn]` | `w_damping` | integer | `0` | Vertical velocity damping where $w$-CFL > 1 (0 = off, 1 = on) | `w_damping` |
| `[dyn]` | `w_crit_cfl` | double | `1.0` | Critical Courant number threshold for $w$-damping | `w_crit_cfl` |
| `[dyn]` | `damp_opt` | integer | `0` | Upper damping layer (0 = off, 1 = diffusive layer, 2 = Rayleigh relaxation, 3 = implicit Rayleigh damping of $w$) | `damp_opt` |
| `[dyn]` | `zdamp` | double | `5000.0` | Depth of upper Rayleigh damping sponge layer in meters | `zdamp` |
| `[dyn]` | `dampcoef` | double | `0.2` | Rayleigh damping coefficient ($s^{-1}$) | `dampcoef` |
| `[dyn]` | `moist_adv_opt` | integer | `1` | Moisture advection scheme (0 = standard, 1 = positive-definite) | `moist_adv_opt` |
| `[dyn]` | `scalar_adv_opt`, `tke_adv_opt` | integer | `1` | Scalar and TKE advection scheme (0 = standard, 1 = positive-definite) | `scalar_adv_opt`, `tke_adv_opt` |
| `[dyn]` | `hybrid_opt`, `etac`, `hypsometric_opt` | int/double | `0`/`2`, `0.2`, `1`/`2` | Hybrid vertical coordinate option, transition $\eta_c$, and hypsometric integration option | `hybrid_opt`, `etac`, `hypsometric_opt` |
| `[dyn]` | `epssm` | double | `0.1` | Acoustic small-step forward-centering factor | `epssm` |
| `[dyn]` | `smdiv`, `emdiv` | double | `0.1`, `0.01` | Acoustic divergence and external mode filters | `smdiv`, `emdiv` |
| `[moist]` | `species` | string list | `""` | Comma-separated list of active water species (`qv,qc,qr,qi,qs,qg`) | WRF moisture array |
| `[phys]` | `mp_physics` | integer | `0` | Microphysics scheme (0 = off, 1 = Kessler, 6 = WSM6, 8 = Thompson) | `mp_physics` |
| `[phys]` | `hail_opt`, `use_mp_re` | integer | `0`, `1` | WSM6 graupel/hail switch (0/1) and microphysics effective radii coupling to RRTMG (0/1) | `hail_opt`, `use_mp_re` |
| `[phys]` | `mp_tend_lim` | double | `10.0` | Heating rate limiter for microphysics ($K/s$) | `mp_tend_lim` |
| `[phys]` | `force_read_thompson` | bool | `false` | Thompson: `qr_acr_qg_V4.dat`, `qr_acr_qsV2.dat` and `freezeH2O.dat` must be in `table_dir` (else an error, as WRF stops) instead of being computed | `force_read_thompson` |
| `[phys]` | `sf_sfclay_physics` | integer | `1` | Surface layer formulation (0 = off, 1 = Revised MM5 / sfclayrev, 2 = Eta / MYJ, 5 = MYNN) | `sf_sfclay_physics` |
| `[phys]` | `isfflx`, `kim_tofd`, `tofd_factor` | int/bool/double | `1`, `false`, `0.003` | Surface heat/moisture flux switch and turbulent orographic form drag options | `isfflx`, `kim_tofd`, `tofd_factor` |
| `[phys]` | `sf_surface_physics`, `num_soil_layers`, `rdmaxalb`, `usemonalb`, `rdlai2d` | int/bool | `2`, `4`, `true`, `false`, `false` | Land surface model (0 = off, 2 = Unified Noah, 4 = Noah-MP v5.1.2, not with nests: wrf.exe 4.8.0 aborts), soil layers (4), snow albedo source, background albedo from wrfinput `ALBBCK` instead of LANDUSE.TBL / VEGPARM.TBL (snow-covered points start with `SNOALB`), Noah LAI from wrfinput `LAI` instead of VEGPARM.TBL | `sf_surface_physics`, `num_soil_layers`, `rdmaxalb`, `usemonalb`, `rdlai2d` |
| `[phys]` | `dveg`, `opt_*`, `soiltstep` | int/double | WRF defaults | Noah-MP sub-options when `sf_surface_physics = 4` (`dveg=4`, `opt_crs=1`, `opt_btr=1`, `opt_runsrf=3`, `opt_runsub=3`, `opt_sfc=1`, `opt_frz=1`, `opt_inf=1`, `opt_rad=3`, `opt_alb=2`, `opt_snf=1`, `opt_tbot=2`, `opt_stc=1`, `opt_gla=1`, `opt_rsf=1`, `opt_soil=1`, `opt_pedo=1`, `opt_crop=0`, `opt_irr=0`, `opt_irrm=0`, `opt_infdv=0`, `opt_tdrn=0`, `opt_tksno=1`, `opt_compact=2`, `opt_scf=1`, `opt_wetland=0`, `soiltstep=0.0`) | `&noah_mp` options |
| `[phys]` | `opt_thcnd` | integer | `1` | Noah soil thermal conductivity (1 = Peters-Lidard, 2 = McCumber-Pielke for soil types 3 and 4) | `opt_thcnd` |
| `[phys]` | `bl_pbl_physics`, `bldt`, `ysu_topdown_pblmix` | int/double | `1`, `0.0`, `1` | Planetary boundary layer scheme (0 = off, 1 = YSU with surface layer 1, 2 = MYJ with surface layer 2, 5 = MYNN-EDMF with surface layer 1, 2 or 5; MYNN not with `fractional_seaice = 1`), call interval (min), YSU top-down mixing | `bl_pbl_physics`, `bldt`, `ysu_topdown_pblmix` |
| `[phys]` | `bl_mynn_tkeadvect`, `bl_mynn_cloudpdf`, `bl_mynn_mixlength`, `bl_mynn_closure`, `bl_mynn_edmf`, `bl_mynn_edmf_mom`, `bl_mynn_edmf_tke`, `bl_mynn_cloudmix`, `bl_mynn_mixqt`, `bl_mynn_ess`, `bl_mynn_mixscalars`, `bl_mynn_output`, `tke_budget` | bool/int/double | WRF defaults (`false`, `2`, `1`, `2.6`, `1`, `1`, `0`, `1`, `0`, `1`, `0`, `0`, `0`) | MYNN-EDMF options (`bl_pbl_physics = 5`): advected QKE, subgrid cloud PDF (`-2..2`), mixing length (`0..2`), closure level (`2.5`, `2.6`, `3.0`), mass flux (`0, 1`), momentum / TKE transport by the plumes, cloud water / ice mixing, total-water mixing, mixing-length form, WRF's unconnected scalar mixing, EDMF_* output and TKE budget fields (`0, 1`); one value for all domains | same names |
| `[phys]` | `sf_mynn_sfcflux_land`, `sf_mynn_sfcflux_water` | integer | `0`, `1` | MYNN surface layer (`sf_sfclay_physics = 5`) roughness / thermal roughness forms over land (`0..3`) and water (`0..4`) | same names |
| `[phys]` | `grav_settling` | integer | `0` | Fog water deposition (1: settling velocity from qc, 2: Katata et al. by land-use class; USGS or MODIFIED_IGBP_MODIS_NOAH land use, needs qc, a PBL scheme and a surface layer) | `grav_settling` |
| `[phys]` | `cu_physics`, `cudt`, `ishallow`, `clos_choice`, `kfeta_trigger`, `convtrans_avglen_m` | int/double | `0`, `0.0`, `0`, `0`, `1`, `30.0` | Cumulus parameterization (0 = off, 1 = KF, 3 = Grell-Freitas, 16 = New Tiedtke), call interval, shallow convection, closure, trigger | `cu_physics`, `cudt`, `ishallow`, `clos_choice`, `kfeta_trigger`, `convtrans_avglen_m` |
| `[phys]` | `ra_lw_physics`, `ra_sw_physics` | integer | `0` | Longwave and shortwave radiation schemes (0 = off, 4 = RRTMG) | `ra_lw_physics`, `ra_sw_physics` |
| `[phys]` | `radt`, `icloud`, `cldovrlp`, `idcor`, `o3input`, `ghg_input`, `aer_opt` | int/double | `0.0`, `1`, `2`, `0`, `2`, `1`, `0` | Radiation call interval (min), cloud fraction, cloud overlap (`1..5`), decorrelation, ozone (`0, 2`), GHG (`0, 1`), aerosols (`0`) | `radt`, `icloud`, `cldovrlp`, `idcor`, `o3input`, `ghg_input`, `aer_opt` |
| `[phys]` | `table_dir` | string | `"."` | Path to directory containing physics tables (`LANDUSE.TBL`, `NoahmpTable.TBL`, etc.) | WRF run directory |
| `[output]` | `dir` | string | `"out/<case>"` | Directory for simulation output files | Output path |
| `[output]` | `format` | string | `"native"` | Output data format (`native`, `wrfout`, `both`) | `io_form_history` |
| `[output]` | `diag_interval` | integer | `10` | Frequency (in steps) for writing rows to `diag.csv` | Diagnostic interval |
| `[output]` | `out_interval` | integer | `nsteps` | Frequency (in steps) for writing field snapshots | `history_interval` |
| `[restart]` | `interval` | double | `0.0` | Restart snapshot write frequency in seconds (0 = disabled) | `restart_interval` |
| `[restart]` | `end` | double | `0.0` | No restart files after this time in seconds (0 = no end) | - |
| `[tracers]` | `names` | string list | `""` | Comma-separated names of passive scalar tracers | Tracer variables |
| `[domains]` | `max_dom`, `feedback`, `smooth_option`, `save_topo_from_real`, `blend_width`, `spec_bdy_width` | integer | `1`, `1`, `2`, `1`, `5`, `5` | Multi-domain nesting count, two-way feedback (`0, 1`), feedback smoother (`0..2`), terrain blending | `&domains` / `&bdy_control` |
| `[d02]`..`[d0N]` | `wrfinput`, `dt`, `out_interval`, `diag_interval`, `parent_id`, `i_parent_start`, `j_parent_start`, `parent_grid_ratio`, `parent_time_step_ratio` | string/int/double | *per-domain* | Child domain input file, timestep, output intervals, parent link, and per-domain dynamics overrides | `&domains` / `&time_control` / `&dynamics` |

---

## 5. Outputs & Data Formats

### Native Format
When `output.format = native` (or `both`), CPPWRF writes raw binary dumps:

1. **`run.json`:** Comprehensive simulation metadata containing grid parameters, vertical levels, physical constants, execution precision, and resolved configuration keys.
2. **`diag.csv`:** Tabular time-series diagnostics recorded every `output.diag_interval` steps. Columns include:
   - `step`: Time step index
   - `time`: Elapsed simulation seconds
   - `max_abs_u`, `max_abs_v`, `max_abs_w`: Peak velocity magnitudes ($m/s$)
   - `min_th_pert`, `max_th_pert`: Extrema of potential temperature perturbation $\theta'$ ($K$)
   - `dry_mass`, `theta_mass`: Total domain mass integrals
   - `total_energy`: Total dry energy
   - `wall_s`: Cumulative wall clock execution time

3. **`snap_<step>.bin` and `snap_<step>.json`:** Multi-field binary array dumps with companion JSON headers detailing variable offsets, shapes, and IEEE precision.

### WRF NetCDF Output
When `output.format = wrfout` (or `both`), CPPWRF writes `wrfout_d0N_<YYYY-MM-DD_HH:MM:SS>` files as
WRF v4.8.0 writes them: the history stream of the WRF I/O schema (`data/wrf_schema/wrf-4.8.0.tsv`) for
the run's options (package conditions, `iofields_filename`, `frames_per_outfile`), with WRF's
variables, dimensions, attributes and global attributes in WRF's order. The options come from the WRF
`namelist.input` next to the wrfinput file, or from `output.wrf_namelist`. In CPU FP64 builds the
files are identical to WRF real*8's (`tests/verify/run_wrf_files_test.sh`, ctest `wrf_files_vs_wrf`).
The one exception is RQCBLTEN / RQIBLTEN in runs without cloud water / ice, which WRF leaves
uninitialised.

File names follow WRF's `history_outname` and `rst_outname` (default `wrfout_d<domain>_<date>` and
`wrfrst_d<domain>_<date>`; `<domain>`, `<date>`, `<year>`, `<month>`, `<day>` and `<hour>` are replaced
as WRF's construct_filename2a does). The files are written to `output.dir`, so the translator rejects a
name with a directory part. With `frames_per_outfile`, a new file starts every that many frames, as in
WRF. CPPWRF writes NetCDF (`io_form_history = 2` and `io_form_restart = 2`); other io_forms (e.g.
11, parallel NetCDF) are errors; `use_netcdf_classic = .true.` writes NetCDF classic files, as in WRF. WRF writes no auxiliary history or input streams by default (every
`io_form_auxhistN` / `io_form_auxinputN` is 0). An auxhist stream with an interval and an io_form is
an error. Other stream options (auxinput streams other than sst_update's wrflowinp, auxinput4) are
errors when they change what wrf.exe reads or writes. ndown.exe (one-way
nesting from a parent run's wrfout) is not part of CPPWRF.

### Restart Files & Resuming Simulations
`restart.interval` (seconds) writes restart files. With `output.format = wrfout`, CPPWRF writes WRF's
`wrfrst_d0N_<date>` files, including the restart stream, the alarm global attributes and every
nest domain. In CPU FP64 builds these files are identical to WRF's. CPPWRF also writes its own
`restart_<step>.nc` files, which hold the complete model state at bitwise precision.
`restart.end` (seconds, default 0 = none) stops writing restart files after that time.

`run.restart_file` resumes a run from either kind of file. A WRF `wrfrst_d01_<date>` file can come
from WRF or from CPPWRF. CPPWRF reads the restart stream (prognostic state, base state, physics
state); fields that WRF does not restart start from WRF's initial value. A nested run reads
`wrfrst_d0N_<date>` from the same directory for each nest. Resumed from WRF's restart files, the CPU
FP64 results are identical to WRF's resumed run (`tests/verify/run_wrf_files_test.sh`, references
from `tests/verify/make_resumed_refs.sh`).

```ini
[restart]
interval = 3600   # restart files every 3600 seconds (1 hour)
```

To resume a simulation from a restart file:
```ini
[run]
case = real
run_time = 7200
restart_file = out/real/restart_000060.nc   # or .../wrfrst_d01_2026-10-12_01:00:00
```

When resuming, pre-restart rows in `diag.csv` are preserved automatically.

---

## 6. Verification & Validation

CPPWRF includes an independent, black-box verification suite located in `tests/verify/`.

### Independent Verification Harness (`verify.py`)
Run standard verification commands:
```bash
# Check idealized simulation physics & conservation laws
python3 tests/verify/verify.py check straka out/straka
python3 tests/verify/verify.py check warm_bubble out/warm_bubble
python3 tests/verify/verify.py check rest_state out/rest_state

# Field-by-field differential comparison against WRF reference run
python3 tests/verify/verify.py wrf /path/to/wrfout_d01 out/straka

# Compare two CPPWRF model runs (e.g. CPU vs GPU or FP32 vs FP64)
python3 tests/verify/verify.py compare out/cpu out/cuda --rtol 1e-4

# Grid convergence test (dx, dx/2, dx/4)
python3 tests/verify/verify.py convergence out/dx100 out/dx50 out/dx25
```

### Real-Data WRF Output Comparison (`compare_real.py`)
Compare CPPWRF real-data simulations directly against WRF reference runs:
```bash
python3 tests/verify/compare_real.py /path/to/wrfout_d01 out/real --interior 6
```

Reports max absolute difference and relative $L_2$ errors across all 3D fields (`U`, `V`, `W`, `PH`, `MU`, `P`, `THM`, `QVAPOR`) and 2D surface variables (`TSK`, `HFX`, `QFX`, `LH`, `UST`, `PBLH`, `RAINNC`).

### Unit & Differential Physics Tests
Run the automated CTest unit and regression suites:
```bash
# Core C++ unit tests
ctest --test-dir build-cpu --output-on-failure

# Python translation tool tests
python3 tests/tools/test_wrf2cppwrf.py -v
```

---

## 7. Performance

Measured on a workstation with an NVIDIA RTX 2060 (6 GB) and a 12-core x86-64 CPU, for the real-data
domain 119 x 89 x 44 (15 km, dt = 60 s, 1 h):

| Configuration | Time per cell and step |
|---|---|
| Dry dynamics, CUDA FP32 | ~32 ns |
| Dynamics + WSM6, CUDA FP32 | ~52 ns |

The 1 h dry run takes 1.26 s on the GPU in FP32 (WRF v4.6.1 with 4 OpenMP threads: 13.2 s) and 4.3 s
in FP64 (serial WRF real*8: 37.9 s). The memory-traffic budget of the large step is in docs/PERF.md.

GPU memory: the column physics keeps its per-column arrays in thread-local memory, which CUDA reserves for every
thread the GPU can hold at once, so the footprint grows with the GPU's size. On an NVIDIA L4 (24 GB) the same domain
in FP64 peaks at about 9.3 GB with Thompson + MYNN-EDMF + RRTMG and 11.3 GB for the two-way nest N3 with Thompson
(the domains share one set of Thompson tables).

`tools/heavy.sh` is a small helper that serialises heavy jobs (builds, runs, test suites) on a shared
workstation: `tools/heavy.sh <command>` runs in one of four job slots, `--light` runs at once with two
OpenMP threads, and `--exclusive` waits for an idle machine (use it for timings).
