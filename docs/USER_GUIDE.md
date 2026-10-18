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

- **Strict IEEE Floating-Point:** Built with `-ffp-contract=off` (no fused multiply-add contraction) to maintain bitwise determinism and exact verification against Fortran reference models.

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
| `CPPWRF_NATIVE` | `ON` | Optimize for host CPU architecture (`-march=native -ffp-contract=off`) |
| `CPPWRF_DIFF_TESTS` | `OFF` | Register WRF Fortran differential test suites with `ctest` |

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
- **Single-Domain Filtering:** Extracts domain 1 settings; issues descriptive warnings when nested domain arrays (`max_dom > 1`) are present.
- **Non-Silent Diagnostics:** Reports clear warnings on unsupported schemes (e.g. `diff_opt=3`, `mp_physics=2`) and informational notices for ignored WPS/quilt metadata.

### Idealized Cases
CPPWRF supports standard ARW idealized benchmarks:
- **`rest_state`:** Stably stratified atmosphere at rest testing hydrostatic balance.
- **`warm_bubble`:** Bryan & Fritsch (2002) dry rising thermal bubble in neutral atmosphere.
- **`straka`:** Straka et al. (1993) cold density current benchmark.
- **`igw`:** Skamarock & Klemp (1994) inertia-gravity wave.
- **`squall2d`:** Weisman & Klemp (1982) squall line with Kessler or WSM6 microphysics.
- **`hill2d`:** Flow over a 2D bell-shaped mountain from WRF `input_sounding`.
- **`quarter_ss`:** 3D supercell storm case from `input_sounding`.

Run an idealized case:
```bash
./build-cpu/cppwrf configs/straka_wrf.ini
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
| `[run]` | `case` | string | *required* | Simulation case (`real`, `hill2d`, `quarter_ss`, `straka`, `squall2d`, `warm_bubble`, `rest_state`, `igw`) | `ideal_case` |
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
| `[domain]` | `bc_x`, `bc_y` | string | `"periodic"` | Lateral boundary condition (`periodic`, `symmetric`, `specified`) | `periodic_x`, `symmetric_xs`, etc. |
| `[domain]` | `eta` | string | `"uniform_z"` | Vertical coordinate distribution (`uniform_z`, `uniform_eta`, `list`) | `eta_levels` |
| `[domain]` | `eta_levels` | float list | `[]` | Explicit comma-separated list of full $\eta$ levels | `eta_levels` |
| `[domain]` | `map_factor` | double | `1.0` | Constant horizontal map projection scale factor | Map scale factor ($m$) |
| `[domain]` | `f0`, `e0` | double | `0.0` | Coriolis parameters $f$ and $e$ in $s^{-1}$ ($f$-plane) | Coriolis parameters |
| `[case]` | `sounding_file` | string | `""` | Path to WRF format sounding file | `input_sounding` |
| `[case]` | `hm` | double | `100.0` | Mountain height in meters (for `hill2d`) | Hill height |
| `[case]` | `xa` | double | `5.0` | Mountain half-width in km (for `hill2d`) | Hill half-width |
| `[case]` | `u0` | double | `10.0` | Uniform inflow wind speed in m/s | Mean wind |
| `[case]` | `delt` | double | `3.0` | Thermal perturbation magnitude in K (for `quarter_ss`) | Bubble perturbation |
| `[time]` | `dt` | double | *required* | Large time step in seconds ($\Delta t$) | `time_step` |
| `[time]` | `time_step_sound` | integer | `0` (auto) | Acoustic small steps per RK stage (0 = auto-computed CFL rule) | `time_step_sound` |
| `[dyn]` | `khdif`, `kvdif` | double | `0.0` | Constant horizontal/vertical eddy diffusion in $m^2/s$ | `khdif`, `kvdif` |
| `[dyn]` | `diff_opt` | integer | `0` | Diffusion option (0 = off, 1 = coordinate-surface diffusion, 2 = physical-space diffusion) | `diff_opt` |
| `[dyn]` | `km_opt` | integer | `1` | Eddy viscosity closure (1 = constant, 2 = 1.5-order TKE, 4 = 2D Smagorinsky) | `km_opt` |
| `[dyn]` | `c_s` | double | `0.25` | Smagorinsky diffusion constant | `c_s` |
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
| `[dyn]` | `scalar_adv_opt` | integer | `1` | Passive tracer advection (0 = standard, 1 = positive-definite) | `scalar_adv_opt` |
| `[dyn]` | `epssm` | double | `0.1` | Acoustic small-step forward-centering factor | `epssm` |
| `[dyn]` | `smdiv`, `emdiv` | double | `0.1`, `0.01` | Acoustic divergence and external mode filters | `smdiv`, `emdiv` |
| `[moist]` | `species` | string list | `""` | Comma-separated list of active water species (`qv,qc,qr,qi,qs,qg`) | WRF moisture array |
| `[phys]` | `mp_physics` | integer | `0` | Microphysics scheme (0 = off, 1 = Kessler, 6 = WSM6, 8 = Thompson) | `mp_physics` |
| `[phys]` | `hail_opt` | integer | `0` | WSM6 6th class density (0 = graupel, 1 = hail) | `hail_opt` |
| `[phys]` | `mp_tend_lim` | double | `10.0` | Heating rate limiter for microphysics ($K/s$) | `mp_tend_lim` |
| `[phys]` | `sf_sfclay_physics` | integer | `1` | Surface layer formulation (0 = off, 1 = Revised MM5 / sfclayrev, 2 = Eta / MYJ) | `sf_sfclay_physics` |
| `[phys]` | `sf_surface_physics` | integer | `2` | Land surface model (0 = off, 2 = Unified Noah LSM 4-layer) | `sf_surface_physics` |
| `[phys]` | `bl_pbl_physics` | integer | `1` | Planetary boundary layer scheme (0 = off, 1 = YSU, 2 = MYJ) | `bl_pbl_physics` |
| `[phys]` | `cu_physics` | integer | `0` | Cumulus parameterization (0 = off, 1 = Kain-Fritsch, 3 = Grell-Freitas, 16 = New Tiedtke); GF also reads `ishallow` and `clos_choice` | `cu_physics` |
| `[phys]` | `ra_lw_physics` | integer | `0` | Longwave radiation scheme (0 = off, 4 = RRTMG LW) | `ra_lw_physics` |
| `[phys]` | `ra_sw_physics` | integer | `0` | Shortwave radiation scheme (0 = off, 4 = RRTMG SW) | `ra_sw_physics` |
| `[phys]` | `radt` | double | `0.0` | Radiation physics call interval (minutes; 0 = every step) | `radt` |
| `[phys]` | `table_dir` | string | `"."` | Path to directory containing physics tables (`LANDUSE.TBL`, etc.) | WRF run directory |
| `[output]` | `dir` | string | `"out/<case>"` | Directory for simulation output files | Output path |
| `[output]` | `format` | string | `"native"` | Output data format (`native`, `wrfout`, `both`) | `io_form_history` |
| `[output]` | `diag_interval` | integer | `10` | Frequency (in steps) for writing rows to `diag.csv` | Diagnostic interval |
| `[output]` | `out_interval` | integer | `nsteps` | Frequency (in steps) for writing field snapshots | `history_interval` |
| `[restart]` | `interval` | double | `0.0` | Restart snapshot write frequency in seconds (0 = disabled) | `restart_interval` |
| `[tracers]` | `names` | string list | `""` | Comma-separated names of passive scalar tracers | Tracer variables |

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

### Restart Files & Resuming Simulations
`restart.interval` (seconds) writes restart files. With `output.format = wrfout`, CPPWRF writes WRF's
`wrfrst_d0N_<date>` files, including the restart stream, the alarm global attributes and every
nest domain. In CPU FP64 builds these files are identical to WRF's. CPPWRF also writes its own
`restart_<step>.nc` files, which hold the complete model state at bitwise precision.

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

`tools/heavy.sh` is a small helper that serialises heavy jobs (builds, runs, test suites) on a shared
workstation: `tools/heavy.sh <command>` runs in one of four job slots, `--light` runs at once with two
OpenMP threads, and `--exclusive` waits for an idle machine (use it for timings).
