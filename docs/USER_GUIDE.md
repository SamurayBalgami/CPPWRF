# CPPWRF User Guide

CPPWRF is a C++20 / CUDA port of the WRF-ARW model v4.8.0: the dynamical core, the physics schemes listed below,
static nests, `real.exe` (`cppwrf_real`) and `ideal.exe`. It reads WRF's inputs (`namelist.input`, `wrfinput_d0N`,
`wrfbdy_d01`, `wrflowinp_d0N`, WRF's `run/` tables) and writes WRF's `wrfout_d0N`, `wrfrst_d0N` and `auxhist3`
files. This guide takes you from WRF's inputs to `wrfout` files: building, the end-to-end example, the translator
of `namelist.input`, every supported option with its WRF name, the options that are rejected and why, GPU builds,
determinism, restarts, nests and the output streams.

---

## Table of Contents
1. [What CPPWRF Is](#1-what-cppwrf-is)
2. [Quick Start](#2-quick-start)
3. [Building CPPWRF](#3-building-cppwrf)
   - [Prerequisites](#prerequisites)
   - [Builds](#builds)
   - [CMake Options](#cmake-options)
   - [Required Flags and Determinism](#required-flags-and-determinism)
4. [From WRF Inputs to wrfout](#4-from-wrf-inputs-to-wrfout)
   - [End-to-End Example (GFS, WPS, cppwrf_real, cppwrf)](#end-to-end-example-gfs-wps-cppwrf_real-cppwrf)
   - [Running from a WRF Run Directory](#running-from-a-wrf-run-directory)
   - [Translating WRF Namelists (`wrf2cppwrf.py`)](#translating-wrf-namelists-wrf2cppwrfpy)
   - [Command Line and Environment](#command-line-and-environment)
   - [Idealized Cases](#idealized-cases)
5. [Supported WRF Options](#5-supported-wrf-options)
   - [Rejected Options and Why](#rejected-options-and-why)
6. [Nested Domains](#6-nested-domains)
7. [GPU Runs](#7-gpu-runs)
8. [Configuration Reference (INI Keys)](#8-configuration-reference-ini-keys)
9. [Outputs & Data Formats](#9-outputs--data-formats)
   - [Native Format](#native-format)
   - [WRF NetCDF Output](#wrf-netcdf-output)
   - [Output Streams](#output-streams)
   - [Restart Files & Resuming Simulations](#restart-files--resuming-simulations)
10. [Verification & Validation](#10-verification--validation)
11. [Performance](#11-performance)

---

## 1. What CPPWRF Is

CPPWRF follows WRF v4.8.0 statement by statement: the same formulas, operation order, loop bounds, staggering, halo
updates, branches, constants and defaults. Every result is checked against WRF itself (docs/STATUS.md):

- **CPU, FP64** (`-DCPPWRF_DOUBLE=ON`, OpenMP): bit-identical to WRF v4.8.0 real*8 (`configure -r8`) run in one
  process, in every variable and output time of the wrfout and wrfrst files, for the 242 suite cases (real-data
  domains in the mid-latitudes, the tropics, the Arctic and the Alps, the physics suites, nests, restarts, 24 h and
  72 h forecasts, WRF's default namelist, WRF's idealized cases).
- **GPU, FP64** (CUDA): bit-identical to WRF and to the CPU build (79 GPU runs of the suite, on NVIDIA L4, RTX PRO
  6000 Blackwell and A100).
- **GPU, FP32**: the fast mode; its distance from WRF real*8 is that of stock WRF real*4 for most physics suites,
  and up to about twice it with ACM2 / Pleim-Xiu, Morrison or the skin SST (docs/FP32_ACCURACY.md). Use it for forecasts and
  ensembles on GPUs with slow FP64; use FP64 when the answer must be WRF real*8's (verification, comparisons with
  WRF runs, continuing a WRF run).

Every WRF option is either ported identically or rejected with an error; nothing is silently ignored. Not ported:
MPI runs (CPPWRF runs on one node: OpenMP threads or one GPU), moving nests, ndown.exe, I/O forms other than NetCDF,
and the schemes listed under [Rejected Options and Why](#rejected-options-and-why).

The design (one source for CPU and GPU kernels, the physics driver, the nests) is described in
docs/PHYSICS_SUITE.md, docs/PHYSICS_DRIVER.md, docs/NESTING.md and docs/NEST_NUMERICS.md.

---

## 2. Quick Start

```bash
git clone <repository> cppwrf && cd cppwrf
# CPU build, FP64 (Debian / Ubuntu: apt install build-essential cmake libnetcdf-dev python3-netcdf4 python3-numpy)
cmake -S . -B build-cpu -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_DOUBLE=ON
cmake --build build-cpu -j "$(nproc)"
ctest --test-dir build-cpu --output-on-failure

# an analytic case that needs no input: the Straka density current (900 s); output in out/straka
build-cpu/cppwrf configs/straka.ini

# WRF's em_quarter_ss supercell, initialised as WRF's ideal.exe does from WRF's own input_sounding
# (WRF_SRC: WRF v4.8.0's source tree; WRF does not need to be built)
build-cpu/cppwrf configs/wrf_quarter_ss.ini --set case.sounding_file=$WRF_SRC/test/em_quarter_ss/input_sounding
```

A real-data forecast from the GFS analysis to wrfout files is the
[end-to-end example](#end-to-end-example-gfs-wps-cppwrf_real-cppwrf):
`examples/quickstart/run.sh` (WPS, `cppwrf_real`, the translator and the model).

These commands, the GPU builds and the end-to-end example were run from a clean source tree on a Google Colab G4
runtime (AMD EPYC 9B45, NVIDIA RTX PRO 6000 Blackwell) with `tools/colab/colab_quickstart.sh`, which also checks that
the GPU FP64 output of the example equals the CPU FP64 output bit for bit.

---

## 3. Building CPPWRF

### Prerequisites
- A C++20 compiler: g++ 11 to 14 are tested (the suite uses Debian 13's g++ 14)
- CMake >= 3.24
- OpenMP
- NetCDF-C with its headers (`libnetcdf-dev`); NetCDF-4 / HDF5 support for compressed output
- Python 3 with `netCDF4` and `numpy` for the tools (`tools/wrf2cppwrf.py` itself needs only the standard library)
- *(GPU builds)* the NVIDIA CUDA Toolkit >= 12 and an NVIDIA GPU (tested: sm_80, sm_89, sm_120)

### Builds
One source tree gives four builds; build directories are independent.

```bash
# CPU, FP64: the configuration that is bit-identical to WRF real*8
cmake -S . -B build-cpu -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_DOUBLE=ON
cmake --build build-cpu -j "$(nproc)"

# GPU, FP32: the fast mode
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=ON -DCPPWRF_DOUBLE=OFF
cmake --build build -j "$(nproc)"

# GPU, FP64: bit-identical to WRF real*8 and to the CPU FP64 build
cmake -S . -B build-gpu-dp -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=ON -DCPPWRF_DOUBLE=ON
cmake --build build-gpu-dp -j "$(nproc)"

# CPU, FP32 (not part of the verification suite)
cmake -S . -B build-cpu-sp -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_DOUBLE=OFF
cmake --build build-cpu-sp -j "$(nproc)"
```

Each build makes `cppwrf` (the model), `src/init_real/cppwrf_real` (real.exe) and the unit tests. A build tree
reads `data/` (the RRTMG tables and the WRF I/O schema) from the source tree it was built from; `cmake --install`
copies the data next to the installed binaries, which find it from their own location (docs/INSTALL.md), and the
environment variable `CPPWRF_DATA_DIR` overrides both.

### CMake Options
| CMake Option | Default | Description |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `Release` (-O3), `RelWithDebInfo`, `Debug`; the values do not depend on it |
| `CPPWRF_DOUBLE` | `OFF` | `ON`: `real` is `double` in the dynamics and physics (FP64); `OFF`: FP32 |
| `CPPWRF_ENABLE_CUDA` | `ON` when a CUDA compiler is found | Build the CUDA backend (one GPU) instead of the OpenMP CPU backend |
| `CMAKE_CUDA_ARCHITECTURES` | `native` | The GPU of the build machine; e.g. `80;89;120` for other GPUs |
| `CPPWRF_NATIVE` | `ON` | Release builds: `-march=native`; `OFF` for portable binaries (same values either way) |
| `CPPWRF_PORTABLE_LIBMVEC` | `ON` | glibc's vector `log` / `log10` as on Intel CPUs, on every x86-64 CPU (the WRF references are Intel-made) |
| `CPPWRF_CUDA_GLIBC_MATH` | `ON` | CUDA FP64: glibc's math functions on the device (bit-identical to WRF); `OFF`: CUDA's math library |
| `CPPWRF_MPI` | `OFF` | `ON`: the program runs as MPI ranks under `mpirun` (needs an MPI library with its C API, e.g. MPICH or Open MPI; section 7, "Tiles and MPI ranks"); `OFF`: no MPI header or library is used |
| `CPPWRF_DIFF_TESTS` | `OFF` | Register the differential tests against WRF's compiled Fortran (need a WRF build, `WRF_DIR`) |

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
- **Machine:** any x86-64 CPU with SSE4.1, FMA and AVX2. glibc's vector `log` / `log10` round differently on AMD
  EPYC than on Intel (an approximate reciprocal, `rcpps`), so WRF differs between the two; the references are
  Intel-made, and CPPWRF computes these two functions as libmvec does on Intel on every CPU (CMake option
  `CPPWRF_PORTABLE_LIBMVEC`, default ON; OFF calls the machine's libmvec, which reproduces Intel-made references
  only on Intel).
- **Threads:** the CPU result does not depend on `OMP_NUM_THREADS`. Every kernel writes its own grid
  points, the radiation columns are independent (McICA seeds from each column's pressures), and the solution
  has no floating-point sums over the domain: the w_damping CFL maxima of the adaptive time step and the
  `diag.csv` sums use `Reducer`, which sums each level in a fixed order
  and then the levels in order, whatever the thread count.
- **CUDA:** nvcc's defaults (`-prec-div=true`, `-prec-sqrt=true`, `-ftz=false`). The FP64 build
  (`-DCPPWRF_DOUBLE=ON`) is bit-identical to WRF's real*8 build and to the CPU build: `--fmad=false` in every
  kernel and glibc's libm / libmvec algorithms on the device (`CPPWRF_CUDA_GLIBC_MATH`, default ON; OFF uses
  CUDA's math library, which agrees with WRF to roundoff). It costs 0 to 7 % of the FP64 stepping time. The FP32
  build contracts FMAs (`--fmad=true`) except RRTMG's batch kernels, Thompson and the physics drivers, and uses
  CUDA's math library: it agrees with WRF's single-precision run to roundoff.

---

## 4. From WRF Inputs to wrfout

A CPPWRF real-data run needs what wrf.exe needs: `namelist.input`, `wrfinput_d0N` and `wrfbdy_d01` (from WRF's
real.exe or from `cppwrf_real`), `wrflowinp_d0N` with `sst_update = 1`, `wrffdda_d0N` with `grid_fdda = 1`, and the table files of WRF's `run/`
directory. The RRTMG tables come with CPPWRF (`data/rrtmg`); everything else is read from WRF's `run/` directory, so
download WRF v4.8.0's source (WRF itself does not need to be built).

### End-to-End Example (GFS, WPS, cppwrf_real, cppwrf)

`examples/quickstart/run.sh` makes a 3-hour forecast from the GFS analysis with WRF's default domains and namelist:
WPS 4.6's default `namelist.wps` (a 150 x 130 domain at 15 km over the south-eastern United States with a 220 x 214
two-way nest at 5 km) and WRF v4.8.0's `test/em_real/namelist.input` (`physics_suite = 'CONUS'`: Thompson, Tiedtke,
MYJ, the Eta surface layer, Noah, RRTMG; `diff_opt 2`, `km_opt 4`, `gwd_opt 1`, `fractional_seaice 1`), with only the
dates, `run_hours = 3`, the low-resolution geography and the number of metgrid levels changed
(`examples/quickstart/namelist.wps`, `namelist.input`).

You need:
- a built WPS (4.6; its `geogrid.exe`, `ungrib.exe`, `metgrid.exe` and their `GEOGRID.TBL`, `METGRID.TBL`,
  `Vtable.GFS`);
- WPS's low-resolution static data (`geog_low_res_mandatory.tar.gz` from the WPS download page, unpacked:
  `WPS_GEOG_LOW_RES`);
- WRF v4.8.0's source tree (for the tables of `run/`);
- a CPPWRF build and Python 3 with `netCDF4`; internet access for two GFS files (about 500 MB each).

```bash
export WPS_DIR=/path/to/WPS WRF_DIR=/path/to/WRF-4.8.0 GEOG_DIR=/path/to/WPS_GEOG_LOW_RES
export CPPWRF_BUILD=build-cpu OMP_NUM_THREADS=$(nproc)
examples/quickstart/run.sh ~/qs
```

The script runs, in `~/qs`:

1. **WPS** (`~/qs/wps`): downloads the GFS 0.25-degree analysis and 3 h forecast of 2026-10-12 00 UTC (NOAA's open
   data on AWS; `GFS_DATE` and `GFS_CYCLE` choose another run, then change the dates of both namelists too), and runs
   `geogrid.exe`, `ungrib.exe` (Vtable.GFS) and `metgrid.exe`: `met_em.d01.*` and `met_em.d02.*`.
2. **cppwrf_real** (`~/qs/run`): sets `num_metgrid_levels` and `num_land_cat` from the met_em files (real.exe stops
   when they differ), then
   `cppwrf_real --namelist namelist.input --met-dir ~/qs/wps --out-dir .` writes `wrfinput_d01`, `wrfinput_d02` and
   `wrfbdy_d01`, as WRF's real.exe would (identical files, docs/REAL.md). WRF's own `real.exe` can be used instead.
3. **cppwrf**: links the files of WRF's `run/` directory next to the inputs (as for wrf.exe), translates the
   namelist (`tools/wrf2cppwrf.py namelist.input -o run.ini --format wrfout`) and runs
   `cppwrf run.ini --set output.dir=~/qs/out`: `out/d01/wrfout_d01_2026-10-12_00:00:00` ... `_03:00:00` and
   `out/d02/wrfout_d02_*` (a nested run writes each domain's files to `output.dir/d0N`; one file per hour, as
   `frames_per_outfile = 1`).

The same `run.ini` runs in every build: `CPPWRF_BUILD=build-gpu-dp` (or `build`) runs the model on the GPU. On the
Colab G4 runtime (AMD EPYC 9B45, RTX PRO 6000 Blackwell) the 3-hour run of both domains took 222 to 273 s of stepping
on 48 CPU threads in FP64, 42 s on the GPU in FP64 and 15 s in FP32 (the whole script 5 to 7 minutes, most of it the
CPU run and the GFS download), and the GPU FP64 output was bit-identical to the CPU FP64 output in every variable of
every wrfout file. WRF's default namelist sets `restart_interval = 7200` (minutes): no wrfrst file is due within the
run, but CPPWRF writes its own `restart_<step>.nc` at the last step of every domain (0.5 and 1.3 GB here;
`--set restart.interval=0` turns restart files off).

### Running from a WRF Run Directory

In a directory prepared for wrf.exe (`namelist.input`, `wrfinput_d0N`, `wrfbdy_d01` and the `run/` tables),
`cppwrf-wrf` is the drop-in wrf.exe: it translates `namelist.input`, runs the model and writes WRF's files where wrf.exe
writes them (`wrfout_d0N_<date>`, `wrfrst_d0N_<date>`, `auxhist3_d0N_<date>`, by `history_outname`, `rst_outname`,
`auxhist3_outname`; every domain in the run directory), with the logs `rsl.out.0000` / `rsl.error.0000` (WRF's
`Timing for main` and `Timing for Writing` lines, `d01 <date> wrf: SUCCESS COMPLETE WRF` at the end). `restart = .true.`
resumes from `wrfrst_d01_<start date>` (`rst_inname`) as wrf.exe does. An option without a port stops it at start-up
with WRF's `FATAL CALLED` block naming the option, its value and what is supported; any failure is a nonzero exit
status. `cppwrf-realexe` is the drop-in real.exe (`met_em.d0N.*` in the directory, `real_em: SUCCESS COMPLETE REAL_EM
INIT`) and `cppwrf-check namelist.input` lists, per option, whether it is ported, ignored or not ported, without
running anything (exit status 0 when CPPWRF runs the namelist):

```bash
$CPPWRF/tools/cppwrf-check namelist.input
$CPPWRF/tools/cppwrf-wrf --mode cpu-fp64 --threads $(nproc)   # or gpu-fp32, gpu-fp64
grep -q "SUCCESS COMPLETE WRF" rsl.out.0000 && ls wrfout_d0*
```

`--mode` picks the binary by what `cppwrf --version` reports (from `--exe`, `CPPWRF_EXE`, the installed `bin/`, the
source tree's `build*/` directories or `PATH`). CPPWRF's own files (the translated `run.ini`, `run.json`, `diag.csv`) go
to `cppwrf_work/`; `--set section.key=value` passes a configuration key on, `--echo` copies the log to the terminal.
The model runs as `cppwrf run.ini` does and the WRF files are byte for byte the same.

By hand, the same run is:

```bash
python3 $CPPWRF/tools/wrf2cppwrf.py namelist.input -o run.ini --format wrfout
OMP_NUM_THREADS=$(nproc) $CPPWRF/build-cpu/cppwrf run.ini --set output.dir=out
```

The translator takes the inputs (`[real] dir`) and the tables (`[phys] table_dir`) from the namelist's directory
(`--real-dir`, `--table-dir` choose others). The model reads the history and restart options it needs for the
wrfout / wrfrst files (iofields, file names, frames per file, NetCDF classic) from the `namelist.input` next to the
wrfinput file (or `output.wrf_namelist`). For results identical to WRF real*8 use inputs of real.exe's `-r8` build
or `cppwrf_real --precision double` (the default); inputs in single precision (WRF's default build) are read too.

### Translating WRF Namelists (`wrf2cppwrf.py`)

`tools/wrf2cppwrf.py` (Python 3, standard library only) turns a WRF `namelist.input` into a CPPWRF configuration
(`.ini`, section 8):

```bash
python3 tools/wrf2cppwrf.py /path/to/namelist.input -o run.ini --format wrfout
```

| Option | Meaning |
|---|---|
| `-o FILE` | the configuration (default: standard output) |
| `--format native\|wrfout\|both` | CPPWRF's native snapshots, WRF's NetCDF files, or both |
| `--real-dir DIR` | the directory of `wrfinput_d01` / `wrfbdy_d01` (default: the namelist's) |
| `--table-dir DIR` | WRF's tables (default: the namelist's directory when it has `LANDUSE.TBL`) |
| `--case NAME` | the case for ideal runs (`em_quarter_ss`, ...; default: real, or WRF's `ideal_case`) |
| `--explicit-files` | write `wrfinput` / `wrfbdy` names in `[real]` |
| `--strict` | warnings are errors |

The translator reads Fortran namelist syntax (repeat counts `3*0`, `1.0d-4`, comments) and applies WRF's own
processing: `physics_suite` and its `-1` options (CONUS, tropical), check_a_mundo's resets (for example
`grav_settling` with `mp_physics = 28`, the soil layers of each land surface model), the per-domain defaults of
the Registry. It exits with status 1 and an `[ERROR]` line for every problem:

- a key that is not a WRF 4.8.0 namelist option (wrf.exe stops reading the namelist too);
- an option value that CPPWRF does not port (section 5);
- a combination that wrf.exe rejects or that would make it fail (for example a longwave scheme without a shortwave
  scheme, `ishallow = 1` without Grell-Freitas, the slab land surface with `fractional_seaice = 1`).

Options that change nothing in a CPPWRF run (WPS and real.exe options such as `num_metgrid_levels`, `dzbot`, I/O
quilting) are listed as `[INFO]` in the configuration's header. `restart = .true.` is not translated: a restart is
`--set run.restart_file=<wrfrst file>` (section 9).

An option merged below the verified tier (`data/option_tiers.txt`, docs/STATUS.md) gives a `[WARNING]` naming its
tier; the model runs it only with `&cppwrf allow_unverified = .true. /` in the namelist, which the translator passes
on as `[cppwrf] allow_unverified = true` (wrf.exe skips the `&cppwrf` group; `cppwrf_real` reads it too).

### Command Line and Environment

```bash
cppwrf CONFIG.ini [--set section.key=value ...]
cppwrf_real --namelist namelist.input --met-dir DIR [--out-dir DIR] [--precision double|float] [--schema FILE]
```

`--set` overrides a key of the configuration (`--set run.run_time=7200 --set output.dir=out2`). A key that CPPWRF
does not know is an error at start-up (`unknown configuration key ...`), so a typing error cannot be ignored.

| Variable | Meaning |
|---|---|
| `OMP_NUM_THREADS` | CPU threads (CPU builds; GPU builds use host threads for RRTM LW and the NetCDF compression). The result does not depend on it |
| `CUDA_VISIBLE_DEVICES` | the GPU of a GPU build (one GPU per run) |
| `CPPWRF_NC_THREADS` | threads that compress the NetCDF-4 chunks (default: the hardware threads, at most 64; 1: no pool; 0: HDF5's own filter); the bytes are the same |
| `CPPWRF_IO_TIMING=1` | the time of each output record |
| `CPPWRF_MEMINFO=1` | GPU builds: the device memory in use at start-up and at the end |
| `CPPWRF_DEVICE_MEMORY_MIB=N` | GPU builds: size the work buffers (column physics, RRTMG) as on a device of N MiB |
| `CPPWRF_KTIME=1` | GPU builds: time every kernel (slower; a table at the end) |
| `CPPWRF_GRAPHS=0` | GPU builds: no CUDA graphs |
| `CPPWRF_TILES=PxQ` | run the domain as P x Q tiles (overrides `run.tiles`; section 7, "Tiles and MPI ranks") |
| `CPPWRF_RANKS=PxQ` | an MPI run's layout: P patches along x times Q along y, P Q = the rank count (overrides `run.ranks`; default: the most square layout with P <= Q) |
| `CPPWRF_RANK_DEVICES=d0,d1,...` | GPU builds under MPI: the device of rank r is `d[r % n]` (default: every rank on the current GPU) |
| `CPPWRF_TILE_DEVICES=d0,d1,...` | GPU builds: the device of tile t is `d[t % n]` (default: all tiles on the current GPU) |
| `CPPWRF_TILE_CHECK=1` | CPU builds: compare every launch of every time step of the tiles with a run without tiles (slow; the first point that differs names the launch, the field and the stale fields) |
| `CPPWRF_WRF_RUNS`, `CPPWRF_WRF_REF` | the WRF run and reference directories that `configs/real_*.ini` and the tests name |

### Idealized Cases
WRF's own idealized cases are initialized by a port of WRF v4.8.0's `ideal.exe`
(`dyn_em/module_initialize_ideal.F`, then what `start_em` does with its output), from the case's
WRF `input_sounding` (in WRF's `test/em_*` directories). In CPU FP64 they give WRF's (real*8) results bit for bit:
- **`em_hill2d_x`:** flow over a 2-D bell-shaped hill (WRF `ideal_case = 1`; `configs/wrf_hill2d_x.ini`).
- **`em_quarter_ss`:** 3-D supercell with Kessler microphysics (`ideal_case = 2`; `configs/wrf_quarter_ss.ini`;
  `configs/wrf_quarter_ss_sym.ini`: symmetric walls).
- **`em_squall2d_x`, `em_squall2d_y`:** 2-D squall line with open lateral boundaries (`ideal_case = 4, 5`;
  `configs/wrf_squall2d_x.ini`, `configs/wrf_squall2d_y.ini`).
- **`em_grav2d_x`:** Straka density current (`ideal_case = 6`; `configs/wrf_grav2d_x.ini`).

The configurations name the sounding as `${CPPWRF_WRF_REF}/run_*/input_sounding`; set `case.sounding_file` to the
file of WRF's test directory instead (section 2). `wrf2cppwrf.py` also translates the `namelist.input` of an ideal
case: the namelist maps to `[domain] nx = e_we - 1`, `ny = e_sn - 1` (`ny = 1` for the x-z cases, WRF `e_sn = 3`),
`nz = e_vert - 1`, `dx`, `dy`, `ztop`, `bc_x`/`bc_y` (`&bdy_control`), `[dyn] hybrid_opt`, `etac` and
`[case] sounding_file`. As in WRF, ideal runs always use `hypsometric_opt = 1`.

CPPWRF also has analytic test cases of its own (not WRF initializations), checked by `tests/verify/verify.py`:
- **`rest_state`:** stably stratified atmosphere at rest (hydrostatic balance).
- **`warm_bubble`:** Bryan & Fritsch (2002) dry rising thermal in a neutral atmosphere.
- **`straka`:** Straka et al. (1993) cold density current.
- **`igw`:** Skamarock & Klemp (1994) inertia-gravity wave.
- **`squall2d`:** Weisman & Klemp (1982) squall line with Kessler or WSM6 microphysics.

---

## 5. Supported WRF Options

Every option below is ported from WRF v4.8.0 and verified against WRF (docs/STATUS.md lists the cases); the INI
key of each is in section 8. Options of `namelist.input` that this table does not name keep WRF's default, and a
non-default value is an error unless it only concerns WPS / real.exe or I/O quilting.

| WRF namelist option | Supported values |
|---|---|
| **&time_control** | |
| `run_days/hours/minutes/seconds`, `start_*`, `end_*` | any (real-data runs from the start of wrfinput) |
| `interval_seconds` | taken from wrfbdy_d01 |
| `history_interval` (`_d`, `_h`, `_m`, `_s`), `frames_per_outfile`, `history_outname` | per domain |
| `restart_interval` (minutes; in a nested run a multiple of `time_step`), `rst_outname` | NetCDF wrfrst files of every domain |
| `io_form_history`, `io_form_restart`, `io_form_input`, `io_form_boundary` | 2 (NetCDF) |
| `use_netcdf_classic` | `.false.`, `.true.` |
| `iofields_filename` | the output streams as in WRF; lines that change an input stream (wrfinput, wrflowinp) are errors |
| `output_diagnostics`, `auxhist3_interval`, `io_form_auxhist3`, `frames_per_auxhist3`, `auxhist3_outname` | 0, 1 (with column physics, not with the adaptive time step) |
| `auxinput4_interval`, `io_form_auxinput4`, `auxinput4_inname` | with `sst_update = 1` (wrflowinp); each domain's interval a multiple of its time step (else WRF stops at the first read: "Could not find matching time in input file"; the translation and the model refuse it) |
| **&fdda** | |
| `grid_fdda` | 0, 1 (analysis nudging toward `wrffdda_d0N`, per domain; real-data runs with column physics, not with the adaptive time step) |
| `gfdda_inname`, `gfdda_interval_m`, `gfdda_end_h`, `io_form_gfdda` | the file name (`wrffdda_d<domain>` is real.exe's; WRF's default is `gfdda_d<domain>`), minutes between the analyses (> 0, a multiple of the domain's time step: as `auxinput4_interval`), the end of the nudging in hours; 2. WRF stops with `gfdda_begin_*`, `gfdda_end_y/d/m/s` or `gfdda_interval_y/d/h/s`, and so does the translation |
| `guv`, `gt`, `gq`, `fgdt`, `fgdtzero`, `if_no_pbl_nudging_uv/t/q`, `if_zfac_uv/t/q`, `k_zfac_uv/t/q`, `if_ramping`, `dtramp_min` | any (per domain but `if_ramping`, `dtramp_min`); `guv`, `gt`, `gq` >= 0 |
| **&domains** | |
| `time_step`, `time_step_fract_num/den` | any stable step |
| `max_dom`, `e_we`, `e_sn`, `e_vert`, `dx`, `dy`, `parent_id`, `i_parent_start`, `j_parent_start`, `parent_grid_ratio`, `parent_time_step_ratio` | static nests (section 6) |
| `feedback`, `smooth_option` | 0, 1; 0, 1, 2 |
| `use_adaptive_time_step`, `target_cfl`, `target_hcfl`, `max_step_increase_pct`, `starting_time_step`, `max_time_step`, `min_time_step` (and `_den`), `step_to_output_time` | single-domain real-data runs without `sst_update` |
| `p_top_requested`, `eta_levels`, `auto_levels_opt`, `dzbot`, `dzstretch_*`, `num_metgrid_*` | real.exe / cppwrf_real |
| **&dynamics** | |
| `rk_ord` | 3 |
| `non_hydrostatic`, `use_theta_m` | `.true.`, 1 |
| `diff_opt`, `km_opt` | 0, 1, 2; 1 (constant K), 2 (1.5-order TKE), 3 (3-D Smagorinsky), 4 (2-D Smagorinsky) |
| `khdif`, `kvdif`, `c_s`, `c_k`, `mix_isotropic`, `mix_upper_bound`, `mix_full_fields`, `tke_drag_coefficient`, `tke_heat_flux`, `tke_upper_bound`, `*_mix2_off`, `*_mix6_off` | any |
| `diff_6th_opt`, `diff_6th_factor` | 0, 1, 2 |
| `damp_opt`, `zdamp`, `dampcoef`, `w_damping`, `w_crit_cfl` | 0-3; any; 0, 1 |
| `h_mom_adv_order` = `h_sca_adv_order`, `v_mom_adv_order` = `v_sca_adv_order` | 5, 6; 3, 4 |
| `moist_adv_opt`, `scalar_adv_opt`, `tke_adv_opt` | 0, 1 |
| `hybrid_opt`, `etac` | 0, 2 (real.exe makes the coordinate) |
| `hypsometric_opt` | 1, 2 |
| `epssm`, `smdiv`, `emdiv`, `time_step_sound`, `top_lid`, `base_temp`, `iso_temp` | any |
| `gwd_opt` | 0, 1 (GWDO; the same on every domain where it is on) |
| **&bdy_control** | |
| `specified`, `spec_bdy_width`, `spec_zone`, `relax_zone` | real-data domain 1 (`spec_zone` 1, `relax_zone` = `spec_bdy_width` - 1) |
| `periodic_x/y`, `open_xs/xe/ys/ye`, `symmetric_xs/xe/ys/ye` | ideal cases (symmetric: no wind across the walls) |
| **&physics** | |
| `physics_suite` | `'CONUS'`, `'tropical'` and their `-1` entries |
| `mp_physics` | 0, 1 (Kessler), 6 (WSM6; `hail_opt` 0, 1), 8 (Thompson), 10 (Morrison two-moment: real data, `morr_rimed_ice` 0 / 1, `progn = 0`, `do_radar_ref = 0`), 28 (aerosol-aware Thompson: one domain, `aer_init_opt` 0, `wif_input_opt` 0 / 1, not with MYNN) |
| `use_mp_re`, `mp_tend_lim`, `do_radar_ref`, `force_read_thompson` | 0, 1; any; 0, 1; `.false.`, `.true.` |
| `cu_physics`, `cudt` | 0, 1 (Kain-Fritsch; `kfeta_trigger` 1, 3), 3 (Grell-Freitas-Li; `ishallow` 0, 1, `clos_choice` 0-16), 6 (Tiedtke), 16 (New Tiedtke); per domain |
| `cu_rad_feedback`, `shallowcu_forced_ra`, `convtrans_avglen_m` | per domain (KF and GF clouds in RRTMG) |
| `bl_pbl_physics`, `bldt` | 0, 1 (YSU; `ysu_topdown_pblmix` 0, 1), 2 (MYJ), 5 (MYNN-EDMF, below), 7 (ACM2, over the surface layer 1 or 7; a microphysics with qc and qi: `mp_physics` 6, 8, 28) |
| `bl_mynn_*`, `tke_budget`, `icloud_bl` | `tkeadvect`, `cloudpdf` -2..2, `mixlength` 0..2, `closure` 2.5 / 2.6 / 3.0, `edmf` 0 / 1, `edmf_mom`, `edmf_tke`, `cloudmix`, `mixqt`, `ess`, `mixscalars`, `output` 0 / 1, `tke_budget` 0 / 1, `icloud_bl` 0 / 1 |
| `sf_sfclay_physics` | 0, 1 (revised MM5; `isfflx` 0, 1, `isftcflx` 0-2, `iz0tlnd` 0-2, `kim_tofd`, `tofd_factor`), 2 (Eta), 5 (MYNN; `sf_mynn_sfcflux_land` 0-3, `sf_mynn_sfcflux_water` 0-4), 7 (Pleim-Xiu, with ACM2 and Noah, `fractional_seaice` 0, `isfflx` 1) |
| `grav_settling` | 0, 1, 2 |
| `sf_surface_physics`, `num_soil_layers` | 0, 1 (5-layer thermal diffusion, 5 layers; `ifsnow`), 2 (Noah, 4 layers; `opt_thcnd` 1, 2), 4 (Noah-MP v5.1.2, default `&noah_mp` options, one domain) |
| `rdmaxalb`, `usemonalb`, `rdlai2d`, `fractional_seaice`, `sst_update` | any; any; Noah; 0, 1; 0, 1 |
| `sst_skin`, `tmn_update`, `lagday` | 0, 1 (skin SST over water, module_sf_sstskin.F); 0, 1 (TMN from the daily mean TSK, module_sf_tmnupdate.F); with `tmn_update = 1` any number of days but 1 (WRF stops there); one value for all domains, Noah (`sf_surface_physics = 2`) only |
| `ra_lw_physics`, `ra_sw_physics`, `radt` | 0, 1 (RRTM LW / Dudhia SW, paired with each other or with RRTMG), 4 (RRTMG) |
| `icloud`, `cldovrlp`, `idcor`, `o3input`, `ghg_input`, `aer_opt`, `swrad_scat` | 0, 1; 0-5; 0, 1; 0, 2; 0, 1; 0; any |
| `bucket_mm`, `bucket_J`, `prec_acc_dt` | any (with column physics; `prec_acc_dt` per domain) |

### Rejected Options and Why

CPPWRF has no MPI: a run is one process on one node, with OpenMP threads or one GPU (WRF's own result changes with
its MPI layout; the references are one-process WRF runs). The options below stop the translator and the model with
an error, never silently:

| WRF option | Why |
|---|---|
| moving nests, vertical nest refinement, nests with their own start or end time, ndown.exe | not ported |
| nests in ideal cases; Noah-MP with nests | not ported; wrf.exe 4.8.0 itself aborts with Noah-MP on nests (heap corruption in the nest's first history output) |
| `io_form_*` other than 2, `nocolons`, auxiliary history / input streams other than auxhist3 and wrflowinp | CPPWRF reads and writes NetCDF; WRF writes no other stream by default |
| `restart = .true.` in the namelist | resume with `run.restart_file` (section 9) |
| `non_hydrostatic = .false.`, `rk_ord` 2, `use_theta_m = 0`, `h_sca_adv_order` different from `h_mom_adv_order`, `spec_exp` != 0 | not ported |
| `km_opt` 5 (SMS-3DTKE) | not ported |
| `mp_physics` other than 0, 1, 6, 8, 10, 28 (40: Morrison with aerosols) | not ported |
| `progn = 1`, `do_radar_ref = 1` with `mp_physics = 10`, `mp_physics = 10` in the ideal cases | not ported (prognostic droplet number needs module_mixactivate; Morrison's 10 cm reflectivity is DOUBLE PRECISION, quad in a -r8 WRF) |
| `cu_physics` other than 0, 1, 3, 6, 16 | not ported |
| `bl_pbl_physics` other than 0, 1, 2, 5, 7; `shcu_physics` | not ported |
| `sf_sfclay_physics` 91 and other values; `sf_urban_physics`, `sf_lake_physics` | not ported |
| `sf_sfclay_physics` 7 without ACM2 or Noah, with `fractional_seaice = 1` or `isfflx = 0` | not ported (the Pleim-Xiu surface layer is verified with ACM2 and Noah; its sea-ice wrapper is not ported) |
| `sf_surface_physics` 3, 5, 7, 8 (RUC, CLM, PX, SSiB); the `pxlsm_*` options and PX soil nudging | not ported (PX needs real.exe's 2-layer soil and its own state, LANDUSEF / SOILCTOP / SOILCBOT input and tables) |
| Noah-MP options other than the defaults (crop, irrigation, tile drainage, ...) | not ported |
| MYNN with `fractional_seaice = 1`; `bl_mynn_edmf` 2, `bl_mynn_edmf_dd`, `bl_mynn_mixnumcon`, `bl_mynn_mixaerosols` 0, `bl_mynn_output` 2, other `bl_mynn_closure` values, `sf_mynn_sfcflux_land` 4 / `_water` 5 | not ported |
| `mp_physics = 28` with `wif_input_opt = 2`, `use_aero_icbc`, `qna_update`, `dust_emis`, nests or MYNN | not ported |
| the Dudhia shortwave with Grell-Freitas and `cu_rad_feedback` | wrf.exe writes NaN cloud fractions there |
| `ra_lw_physics` / `ra_sw_physics` other than 0, 1, 4; a longwave scheme without a shortwave scheme or the reverse; `aer_opt` != 0 | not ported; wrf.exe stops |
| `nwp_diagnostics`, `do_avgflx_em`, `do_avgflx_cugd` and the other diagnostic packages | not ported (diagnostic_output_nwp, cal_helicity, the averaged fluxes) |
| `iau`, `dyn_lightning_option`, `topo_wind`, `slope_rad`, `topo_shading`, `swint_opt`, `seaice_thickness_opt`, `seaice_snowdepth_opt`, `ra_call_offset`, `shalwater_z0`, `scm_force_flux`, spectral nudging (`grid_fdda = 2`), surface analysis nudging (`grid_sfdda`: OBSGRID's `wrfsfdda_d0N`), observation nudging (`obs_nudge_opt`), stochastic physics, chemistry | not ported |
| the adaptive time step with nests, `sst_update`, `sst_skin`, `tmn_update`, `grid_fdda`, ideal cases or `output_diagnostics` | not ported |
| `sst_skin = 1` or `tmn_update = 1` with a land surface model other than Noah | not ported (no reference run); `tmn_update` resumes from wrfrst files only (NYEAR, NDAY) |
| a wind across a symmetric wall (ideal cases) | WRF's solution breaks down within a minute there (the wall keeps the normal velocity while the sign of `muu` flips) |
| real-data runs with `periodic_*`, `open_*` or `symmetric_*` sides | WRF's channel set-up is not ported; WRF stops with the others |

The translator prints the reason with each error. docs/STATUS.md ("Not yet implemented") keeps the list current.

---

## 6. Nested Domains

CPPWRF runs WRF's static nests: one-way (`feedback = 0`) and two-way (`feedback = 1`, `smooth_option` 0, 1, 2),
several children of one parent and nests inside nests, odd `parent_grid_ratio` and any `parent_time_step_ratio`, each
nest reading its own `wrfinput_d0N` (`input_from_file = .true.`, as real.exe or cppwrf_real makes them). The parent
forces the nest boundary as WRF's `med_nest_force` does (`bdy_interp1`, SINT), and two-way nests feed back with
`copy_fcn` and the sm121 / smdsm smoothers, `blend_terrain` and `adjust_tempqv`. `cu_physics`, `cudt`, `radt`, `bldt`,
`cu_rad_feedback` and `prec_acc_dt` may differ between domains; the other physics options are the same on every
domain, as the translator checks. Each domain's files go to `output.dir/d0N`. A nested configuration from the translator has a `[domains]` section and one
`[d0N]` section per nest (section 8). Every domain writes its own `wrfout_d0N` and `wrfrst_d0N` files, identical to
WRF's (cases N1, N2, N3, N3sib, N3gc of docs/STATUS.md). docs/NESTING.md and docs/NEST_NUMERICS.md describe the port.

---

## 7. GPU Runs

A CUDA build (`-DCPPWRF_ENABLE_CUDA=ON`) runs the whole model on one NVIDIA GPU: the dynamics, the column physics,
RRTMG and the nests share their kernels with the CPU build (`parallel_for` / `column_for` with `CPPWRF_HD` functors).
The RRTM longwave runs on host threads. The input, the configuration and the output files are the same as for the
CPU build.

- **FP64** (`-DCPPWRF_DOUBLE=ON`): bit-identical to WRF real*8 and to the CPU FP64 build. The CMake files compile
  every FP64 kernel with `--fmad=false` (no fused multiply-add, as gfortran on x86-64 without `-march`) and run
  glibc's libm / libmvec algorithms on the device (`CPPWRF_CUDA_GLIBC_MATH`); this costs 0 to 7 % of the stepping
  time. Verified on L4 (sm_89), RTX PRO 6000 Blackwell (sm_120) and A100 (sm_80).
- **FP32** (`-DCPPWRF_DOUBLE=OFF`): the fast mode, with CUDA's math library and FMA contraction (except in RRTMG's
  batch kernels, Thompson and the physics drivers). Its difference from WRF real*8 is of the size of stock WRF
  real*4's (docs/FP32_ACCURACY.md: the measurement, its method and tools/verify/fp32_accuracy.py, which compares any
  two sets of wrfout files).
- **Memory:** a kernel's local arrays are reserved by CUDA for every thread the GPU can hold at once, so that part
  of the footprint grows with the GPU (the largest frame is 10.8-12.4 KB per thread in FP64 and 4.2-11 KB in FP32:
  0.4 GB on an RTX 2060, 3 GB on an RTX PRO 6000). The large per-column arrays of the column physics live in a
  work buffer of at most 1/8 of the device memory instead. `CPPWRF_MEMINFO=1` prints the footprint. A 6 GB GPU
  runs every verified physics suite on a 120 x 90 x 45 domain in FP32 and FP64 ("GPU memory" below has the
  measured minimum per suite).
- **Determinism:** the GPU build has no atomics in floating-point sums and gives the same bits in every run.

### Tiles and MPI ranks (domain decomposition)

One domain, and its nests, can run as P x Q *patches* (docs/MPI.md has the design). A patch is a rectangle of the
horizontal grid with its own memory: its points plus a halo of 3 points, not the whole domain. A patch is run by a
*worker*, of one of two kinds that share all of the decomposition code (src/core/tiles.hpp):

- **Tiles**: `run.tiles = PxQ` (or `CPPWRF_TILES=PxQ`), host threads of one process, each with its own fields, its
  own CUDA stream and device (`CPPWRF_TILE_DEVICES`; by default every tile on the one GPU).
- **MPI ranks**: a build with `-DCPPWRF_MPI=ON` under `mpirun`, one patch per rank:

```bash
cmake -B build-mpi -DCMAKE_BUILD_TYPE=Release -DCPPWRF_DOUBLE=ON -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_MPI=ON
cmake --build build-mpi -j --target cppwrf
OMP_NUM_THREADS=4 mpirun -np 6 build-mpi/cppwrf configs/real_p4_r8.ini                 # 2 x 3 patches
CPPWRF_RANKS=3x2 OMP_NUM_THREADS=4 mpirun -np 6 build-mpi/cppwrf configs/real_p4_r8.ini
# a CUDA build with CPPWRF_MPI: two ranks on one GPU, or one GPU per rank
mpirun -np 2 build-cuda-mpi/cppwrf CONFIG.ini
CPPWRF_RANK_DEVICES=0,1 mpirun -np 2 build-cuda-mpi/cppwrf CONFIG.ini
```

  The layout is `run.ranks = PxQ` (or `CPPWRF_RANKS`), P patches along x times Q along y with P Q = the rank count;
  without it the most square layout with P <= Q, as WRF's. Each rank uses `OMP_NUM_THREADS` threads (CPU build), so
  ranks times threads should not exceed the cores. The same binary without `mpirun` (or with one rank) runs as one
  process. `run.tiles` in an MPI run is an error. At the end the run prints the ranks' peak host memory.

Inside a time step a worker computes only the points it owns (the edge patches also own the last faces and the
lateral halos). The workers meet in halo exchanges (every halo fill first receives, from the neighbouring workers,
the points they own within 3 points of the patch, and the sources of a periodic wrap), in the domain-wide
reductions (maxima; sums in the order of a run with one patch) and in the output: **one worker (patch 0) writes
every file**, from fields gathered one at a time, so the files are those of a run with one patch, in `output.dir`.
A point is always computed by one worker from the values a run with one patch has, so **the result is the same
bits for every worker count and layout**: CPU FP64 = one-process WRF real*8, CUDA FP64 the same; WRF's own MPI
runs do not have this property (docs/PERF.md). Nested runs decompose every domain into the same P x Q patches of
its own grid; the nest forcing and the feedback read, from the other domain's fields, the box over the worker's
patch.
Without `run.tiles` and outside `mpirun` nothing of this runs.

Limits: a worker's memory is its patch's fields; a real-data start reads only the patch's window of wrfinput's
3-D variables, but the 2-D variables, the wrfbdy records, a restart file's fields (one at a time) and the ideal
cases' initial arrays are whole on every worker while they are read, the lateral boundary slabs are whole edges,
and the writer holds the output frame it writes (memory that does not shrink with the worker count). CUDA graphs are off with several workers (the exchanges synchronise them from the host). A
patch must be at least 3 points wide; a 2-D (x-z) domain has one patch along y (`run.ranks = Nx1`), a y-z domain
one along x. `CPPWRF_TILE_CHECK` is for tiles without nests. Proved: the layouts and cases of docs/MPI.md, section
5, on one machine; runs on several GPUs or several nodes have not been made (docs/MPI.md, section 6).

The suite's cases `*_tPxQ` (tests/verify/cases/85_tiles.py) run other cases as tiles against the same WRF
references; `tools/colab/colab_tiles.sh` runs suite cases with several tile layouts and the check on a Colab
runtime, `tools/colab/colab_mpi.sh` runs suite cases under `mpirun` with several rank layouts (CPU, and
`--gpu-cases` with several ranks on the runtime's GPU), each against the WRF references and byte for byte against
the run with one rank; `tools/bench/colab_bench.sh --scaling` measures the strong scaling next to WRF dmpar.

---

## 8. Configuration Reference (INI Keys)

### Configuration Keys Table
Every option supported by the CPPWRF runtime is cataloged below, along with its type, default, WRF equivalent, and code source.

A key of the configuration file or of `--set` that CPPWRF does not read anywhere (a misspelled key, a key in the
wrong section, an option that is not ported) is an error at start-up (`unknown configuration key 'phys.mp_physic'`);
a known key the run does not read (an option of a scheme it does not use) is accepted with the warning
`unaccessed config key`, as WRF accepts namelist entries it does not need. The list of keys
comes from the sources at build time (`src/cases/config_keys.cmake`). `[domains]` keys are checked when the run is
nested; a `[d0N]` key reaches that domain's configuration and is checked there.

| Section | Key | Type | Default | Description | WRF Equivalent |
|---|---|---|---|---|---|
| `[run]` | `case` | string | *required* | Simulation case (`real`, `em_hill2d_x`, `em_quarter_ss`, `em_squall2d_x`, `em_squall2d_y`, `em_grav2d_x`, `straka`, `squall2d`, `warm_bubble`, `rest_state`, `igw`) | `ideal_case` |
| `[run]` | `run_time` | double | `0.0` | Total simulation duration in seconds | `run_days/hours/minutes/seconds` |
| `[run]` | `restart_file` | string | `""` | Path to NetCDF restart file (`restart_*.nc`) to resume execution | `restart` / restart filename |
| `[run]` | `start_date` | string | `"2000-01-01_00:00:00"` | UTC start timestamp string | `start_year/month/day/hour/etc` |
| `[run]` | `stop_on_nonfinite` | bool | `true` | Stop when a diagnostic step finds a NaN / Inf (false: go on, as WRF does) | none |
| `[run]` | `check_nonfinite_every_step` | bool | `false` | Look for NaN / Inf after every step (default: on the `output.diag_interval` steps); the message names the first non-finite field and point | none |
| `[run]` | `tiles` | string | `"1"` | `PxQ` (or `P`): the domain as P x Q tiles, host threads of one process (section 7, "Tiles and MPI ranks") | none |
| `[run]` | `ranks` | string | `""` | `PxQ`: the layout of an MPI run's ranks (P Q = the rank count; default the most square with P <= Q) | `nproc_x`, `nproc_y` |
| `[cppwrf]` | `allow_unverified` | bool | `false` | Run options below the verified tier (`data/option_tiers.txt`: column_verified, unverified); without it such an option stops the run, naming the option, value and tier. With it: a warning per option and the global attribute `CPPWRF_UNVERIFIED_OPTIONS` (`mp_physics=4:column_verified;...`) in every wrfout / wrfrst file (docs/STATUS.md) | `&cppwrf allow_unverified = .true.` (read by the translator and cppwrf_real; wrf.exe skips the group) |
| `[run]` | `restart_any_commit` | bool | `false` | Resume from a CPPWRF `restart_*.nc` file written by another commit | none |
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
| `[phys]` | `mp_physics` | integer | `0` | Microphysics scheme (0 = off, 1 = Kessler, 6 = WSM6, 8 = Thompson, 10 = Morrison two-moment: real data, `[moist] species` `qv,qc,qr,qi,qs,qg,qni,qns,qnr,qng`; 28 = aerosol-aware Thompson: real data, one domain, `[moist] species` with `qnc,qnwfa,qnifa,qnbca`, `CCN_ACTIVATE.BIN` in `table_dir`) | `mp_physics` |
| `[phys]` | `morr_rimed_ice` | integer | `1` | mp_physics = 10: 1 = hail, 0 = graupel as the rimed ice (WRF: any value other than 1 means graupel) | `morr_rimed_ice` |
| `[phys]` | `hail_opt`, `use_mp_re` | integer | `0`, `1` | WSM6 graupel/hail switch (0/1) and microphysics effective radii coupling to RRTMG (0/1) | `hail_opt`, `use_mp_re` |
| `[phys]` | `do_radar_ref` | integer | `0`, `1` | 10-cm radar reflectivity REFL_10CM (dBZ) in the history: WSM6 and Thompson (8, 28) at the steps that end on a history time (as WRF's diag_flag), 0 with other microphysics; not in restart files (0 after a restart until the next history time) | `do_radar_ref` |
| `[phys]` | `mp_tend_lim` | double | `10.0` | Heating rate limiter for microphysics ($K/s$) | `mp_tend_lim` |
| `[phys]` | `wif_input_opt` | integer | `0` | mp_physics = 28: 0 and 1 run alike (real.exe needs 1), 2 (black-carbon aerosols) is an error | `wif_input_opt` |
| `[phys]` | `force_read_thompson` | bool | `false` | Thompson: `qr_acr_qg_V4.dat`, `qr_acr_qsV2.dat` and `freezeH2O.dat` must be in `table_dir` (else an error, as WRF stops) instead of being computed | `force_read_thompson` |
| `[phys]` | `sf_sfclay_physics` | integer | `1` | Surface layer formulation (0 = off, 1 = Revised MM5 / sfclayrev, 2 = Eta / MYJ, 5 = MYNN, 7 = Pleim-Xiu with ACM2 and Noah) | `sf_sfclay_physics` |
| `[phys]` | `isfflx`, `kim_tofd`, `tofd_factor` | int/bool/double | `1`, `false`, `0.003` | Surface heat/moisture flux switch and turbulent orographic form drag options | `isfflx`, `kim_tofd`, `tofd_factor` |
| `[phys]` | `sf_surface_physics`, `num_soil_layers`, `rdmaxalb`, `usemonalb`, `rdlai2d` | int/bool | `2`, `4`, `true`, `false`, `false` | Land surface model (0 = off, 1 = 5-layer thermal diffusion (slab; not with `fractional_seaice = 1`), 2 = Unified Noah, 4 = Noah-MP v5.1.2, not with nests: wrf.exe 4.8.0 aborts), soil layers (4; 5 for the slab), snow albedo source, background albedo from wrfinput `ALBBCK` instead of LANDUSE.TBL / VEGPARM.TBL (snow-covered points start with `SNOALB`), Noah LAI from wrfinput `LAI` instead of VEGPARM.TBL | `sf_surface_physics`, `num_soil_layers`, `rdmaxalb`, `usemonalb`, `rdlai2d` |
| `[phys]` | `dveg`, `opt_*`, `soiltstep` | int/double | WRF defaults | Noah-MP sub-options when `sf_surface_physics = 4` (`dveg=4`, `opt_crs=1`, `opt_btr=1`, `opt_runsrf=3`, `opt_runsub=3`, `opt_sfc=1`, `opt_frz=1`, `opt_inf=1`, `opt_rad=3`, `opt_alb=2`, `opt_snf=1`, `opt_tbot=2`, `opt_stc=1`, `opt_gla=1`, `opt_rsf=1`, `opt_soil=1`, `opt_pedo=1`, `opt_crop=0`, `opt_irr=0`, `opt_irrm=0`, `opt_infdv=0`, `opt_tdrn=0`, `opt_tksno=1`, `opt_compact=2`, `opt_scf=1`, `opt_wetland=0`, `soiltstep=0.0`) | `&noah_mp` options |
| `[phys]` | `ifsnow` | integer | `1` | Slab (`sf_surface_physics = 1`): 1 limits the skin and top soil temperatures under snow to 273.16 K, 0 does not | `ifsnow` |
| `[phys]` | `opt_thcnd` | integer | `1` | Noah soil thermal conductivity (1 = Peters-Lidard, 2 = McCumber-Pielke for soil types 3 and 4) | `opt_thcnd` |
| `[phys]` | `bl_pbl_physics`, `bldt`, `ysu_topdown_pblmix` | int/double | `1`, `0.0`, `1` | Planetary boundary layer scheme (0 = off, 1 = YSU with surface layer 1, 2 = MYJ with surface layer 2, 5 = MYNN-EDMF with surface layer 1, 2 or 5; MYNN not with `fractional_seaice = 1`; 7 = ACM2 with surface layer 1 or 7, needs qc and qi), call interval (min), YSU top-down mixing | `bl_pbl_physics`, `bldt`, `ysu_topdown_pblmix` |
| `[phys]` | `bl_mynn_tkeadvect`, `bl_mynn_cloudpdf`, `bl_mynn_mixlength`, `bl_mynn_closure`, `bl_mynn_edmf`, `bl_mynn_edmf_mom`, `bl_mynn_edmf_tke`, `bl_mynn_cloudmix`, `bl_mynn_mixqt`, `bl_mynn_ess`, `bl_mynn_mixscalars`, `bl_mynn_output`, `tke_budget` | bool/int/double | WRF defaults (`false`, `2`, `1`, `2.6`, `1`, `1`, `0`, `1`, `0`, `1`, `0`, `0`, `0`) | MYNN-EDMF options (`bl_pbl_physics = 5`): advected QKE, subgrid cloud PDF (`-2..2`), mixing length (`0..2`), closure level (`2.5`, `2.6`, `3.0`), mass flux (`0, 1`), momentum / TKE transport by the plumes, cloud water / ice mixing, total-water mixing, mixing-length form, WRF's unconnected scalar mixing, EDMF_* output and TKE budget fields (`0, 1`); one value for all domains | same names |
| `[phys]` | `sf_mynn_sfcflux_land`, `sf_mynn_sfcflux_water` | integer | `0`, `1` | MYNN surface layer (`sf_sfclay_physics = 5`) roughness / thermal roughness forms over land (`0..3`) and water (`0..4`) | same names |
| `[phys]` | `grav_settling` | integer | `0` | Fog water deposition (1: settling velocity from qc, 2: Katata et al. by land-use class; USGS or MODIFIED_IGBP_MODIS_NOAH land use, needs qc, a PBL scheme and a surface layer) | `grav_settling` |
| `[phys]` | `cu_physics`, `cudt`, `ishallow`, `clos_choice`, `kfeta_trigger`, `convtrans_avglen_m` | int/double | `0`, `0.0`, `0`, `0`, `1`, `30.0` | Cumulus parameterization (0 = off, 1 = KF, 3 = Grell-Freitas, 16 = New Tiedtke), call interval, shallow convection, closure, trigger | `cu_physics`, `cudt`, `ishallow`, `clos_choice`, `kfeta_trigger`, `convtrans_avglen_m` |
| `[phys]` | `ra_lw_physics`, `ra_sw_physics` | integer | `0` | Longwave and shortwave radiation schemes (0 = off, 4 = RRTMG; `ra_lw_physics = 1`: RRTM longwave, `ra_sw_physics = 1`: Dudhia shortwave; RRTM and Dudhia pair with each other or with RRTMG) | `ra_lw_physics`, `ra_sw_physics` |
| `[phys]` | `swrad_scat` | double | `1.0` | Dudhia clear-sky scattering factor (CSSCA = swrad_scat * 1e-5) | `swrad_scat` |
| `[phys]` | `radt`, `icloud`, `cldovrlp`, `idcor`, `o3input`, `ghg_input`, `aer_opt` | int/double | `0.0`, `1`, `2`, `0`, `2`, `1`, `0` | Radiation call interval (min), cloud fraction, cloud overlap (`0..5`; 0 = clear sky only, the McICA clouds are left out), decorrelation, ozone (`0, 2`), GHG (`0, 1`), aerosols (`0`) | `radt`, `icloud`, `cldovrlp`, `idcor`, `o3input`, `ghg_input`, `aer_opt` |
| `[phys]` | `cu_rad_feedback` | bool | `false` | Cumulus clouds in the radiation (WRF's icloud_cu): Kain-Fritsch's `CLDFRA_DP` + `CLDFRA_SH` over the resolved cloud fraction and `QC_CU` / `QI_CU` times that fraction, Grell-Freitas' averaged `QC_CU` / `QI_CU` added to qc / qi before `cal_cldfra1`; no effect with Tiedtke / New Tiedtke; per domain (`[d0N] phys_cu_rad_feedback`) | `cu_rad_feedback` |
| `[phys]` | `bucket_mm`, `bucket_j`, `prec_acc_dt`, `prec_acc_opt` | double/bool | `-1`, `-1`, `0`, `prec_acc_dt > 0` | WRF's accumulation buckets and precipitation windows (diagnostic_output_calc; with the column physics): at each step that starts at a whole hour of the run, `RAINC` / `RAINNC` above `bucket_mm` lose `bucket_mm` and count `I_RAINC` / `I_RAINNC`; with `bucket_mm > 0` and `bucket_j > 0` the radiation accumulations `ACSWUPT` ... `ACLWDNBC` likewise (`I_ACSWUPT` ...); `PREC_ACC_C`, `PREC_ACC_NC`, `SNOW_ACC_NC` over `prec_acc_dt` minutes (per domain: `[d0N] phys_prec_acc_dt`); `prec_acc_opt`: the PREC_ACC fields on every domain (some domain's `prec_acc_dt > 0`) | `bucket_mm`, `bucket_J`, `prec_acc_dt` |
| `[phys]` | `output_diagnostics`, `auxhist3_interval` | int | `0`, `0` | WRF's climate diagnostics (clwrf_output_calc; with the column physics): over windows of `auxhist3_interval` minutes (per domain: `[d0N] phys_auxhist3_interval`) the minima, maxima, times of the extremes, means and standard deviations of T2, Q2, SKINTEMP, U10 / V10 (wind speed), RAINCV and RAINNCV (`T2MIN` ... `RAINNCVSTD`, restart fields) | `output_diagnostics`, `auxhist3_interval` |
| `[phys]` | `table_dir` | string | `"."` | Path to directory containing physics tables (`LANDUSE.TBL`, `NoahmpTable.TBL`, etc.) | WRF run directory |
| `[output]` | `dir` | string | `"out/<case>"` | Directory for simulation output files | Output path |
| `[output]` | `format` | string | `"native"` | Output data format (`native`, `wrfout`, `both`) | `io_form_history` |
| `[output]` | `diag_interval` | integer | `10` | Frequency (in steps) for writing rows to `diag.csv` | Diagnostic interval |
| `[output]` | `out_interval` | integer | `nsteps` | Frequency (in steps) for writing field snapshots | `history_interval` |
| `[output]` | `auxhist3` | bool | `false` | WRF's auxhist3 files (`auxhist3_outname`, default `auxhist3_d<domain>_<date>`; `frames_per_auxhist3`) with the `output_diagnostics` fields every `phys.auxhist3_interval` minutes; needs `format = wrfout` or `both` | `io_form_auxhist3 = 2` |
| `[restart]` | `interval` | double | `0.0` | Restart snapshot write frequency in seconds (0 = disabled) | `restart_interval` |
| `[restart]` | `end` | double | `0.0` | No restart files after this time in seconds (0 = no end) | - |
| `[restart]` | `native` | bool | `true` | Also write CPPWRF's own `restart_<step>.nc` files next to the WRF `wrfrst` files (`cppwrf-wrf` sets false) | - |
| `[output]` | `wrf_dir` | string | `output.dir` (per domain) | Directory of the WRF files (wrfout, wrfrst, auxhist3) of every domain; `cppwrf-wrf` sets the run directory | the run directory |
| `[output]` | `wrf_timing` | bool | `false` | Print WRF's `Timing for main` / `Timing for Writing` lines (the rsl logs of `cppwrf-wrf`) | - |
| `[tracers]` | `names` | string list | `""` | Comma-separated names of passive scalar tracers | Tracer variables |
| `[domains]` | `max_dom`, `feedback`, `smooth_option`, `save_topo_from_real`, `blend_width`, `spec_bdy_width` | integer | `1`, `1`, `2`, `1`, `5`, `5` | Multi-domain nesting count, two-way feedback (`0, 1`), feedback smoother (`0..2`), terrain blending | `&domains` / `&bdy_control` |
| `[d02]`..`[d0N]` | `wrfinput`, `dt`, `out_interval`, `diag_interval`, `parent_id`, `i_parent_start`, `j_parent_start`, `parent_grid_ratio`, `parent_time_step_ratio` | string/int/double | *per-domain* | Child domain input file, timestep, output intervals, parent link, and per-domain dynamics overrides | `&domains` / `&time_control` / `&dynamics` |

---

## 9. Outputs & Data Formats

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
`io_form_auxhistN` / `io_form_auxinputN` is 0). The auxhist3 stream of `output_diagnostics = 1` is written
(`io_form_auxhist3 = 2`, `auxhist3_outname`, `frames_per_auxhist3`; identical to WRF's in CPU FP64 builds);
another auxhist stream with an interval and an io_form is an error. Other stream options (auxinput streams other than sst_update's wrflowinp, auxinput4, and grid_fdda's wrffdda, auxinput10) are
errors when they change what wrf.exe reads or writes. ndown.exe (one-way
nesting from a parent run's wrfout) is not part of CPPWRF.

### Output Streams

| WRF stream | CPPWRF | Files |
|---|---|---|
| input (`wrfinput_d0N`, stream i:0) | read at the start of every domain | from real.exe or cppwrf_real |
| boundary (`wrfbdy_d01`) | read for domain 1 every `interval_seconds` | from real.exe or cppwrf_real |
| auxinput4 (`wrflowinp_d0N`) | read with `sst_update = 1` every `auxinput4_interval` | from real.exe or cppwrf_real |
| auxinput10 (`wrffdda_d0N`) | read with `grid_fdda = 1` at the start and every `gfdda_interval_m` up to `gfdda_end_h` | from real.exe or cppwrf_real (`gfdda_inname`) |
| history (h:0) | written every `history_interval` per domain | `wrfout_d0N_<date>` (`history_outname`, `frames_per_outfile`, `iofields_filename`) |
| restart (r) | written every `restart_interval` for every domain | `wrfrst_d0N_<date>` (`rst_outname`) and CPPWRF's `restart_<step>.nc` |
| auxhist3 (h:3) | written with `output_diagnostics = 1` every `auxhist3_interval` | `auxhist3_d0N_<date>` (`auxhist3_outname`, `frames_per_auxhist3`) |
| CPPWRF native | `output.format = native` or `both` | `run.json`, `diag.csv`, `snap_<step>.bin` / `.json` |

In a nested run the files of domain N are in `output.dir/d0N`. Any other auxiliary stream with an interval and an
io_form is an error (WRF writes none by default).

### Restart Files & Resuming Simulations
`restart.interval` (seconds) writes restart files. With `output.format = wrfout`, CPPWRF writes WRF's
`wrfrst_d0N_<date>` files, including the restart stream, the alarm global attributes and every
nest domain. In CPU FP64 builds these files are identical to WRF's. CPPWRF also writes its own
`restart_<step>.nc` files, which hold the complete model state at bitwise precision: at every restart time and, with
`restart.interval > 0`, at the last step of the run.
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

CPPWRF's own `restart_<step>.nc` files record a format version (`CPPWRF_RESTART_VERSION`), the commit of the sources
(`CPPWRF_COMMIT`) and the run's options (`CPPWRF_OPTIONS`: every configuration key except `[run]`, `[restart]`,
`[output]`, `[real]`, `time.history_interval_s`, `phys.auxhist3_interval` and `phys.table_dir`). A resumed run
refuses a file of another format version, of another commit (unless `run.restart_any_commit = true`) or of a run with
other options, and names the differing options. WRF's `wrfrst` files are the restart across versions and changed
options.

The restart alarm of a resumed run is WRF's: input_wrf sets its first ring to the restart time plus the file's
`WRF_ALARM_SECS_TIL_NEXT_RING_51`, it rings at the first step at or past that time and counts the interval from that
step (ESMF_ClockAdvance). With a `restart.interval` that is not a multiple of `time.dt` a resumed run therefore writes
its restart files at other times than the uninterrupted run (dt 50 s, interval 120 s, resumed at 150 s: 250, 400, 500,
650 s instead of 250, 400, 500, 600 s), and with `restart.interval = 0` it still writes one set at that first ring, as
WRF does. CPPWRF's own `restart_<step>.nc` (format version 2) carries the same attribute and resumes the same way.

When resuming, pre-restart rows in `diag.csv` are preserved automatically.

---

## 10. Verification & Validation

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
python3 -m unittest tests.tools.test_wrf2cppwrf
```

### The CPU Suite on GitHub or Colab
The model runs against WRF need the WRF reference runs (WRF's inputs, its wrfout / wrfrst files and its builds). They
are public: the tars of `tests/verify/vm_data.txt` are release assets of
[cppwrf-data](https://github.com/Xassemblianist/cppwrf-data) (release `data-v1`: parts of 512 MiB, `SHA256SUMS`,
`TARS.txt`, `INDEX.tsv.gz`; `diff_data.tar` is replaced by `diff_runs.tar`, its WRF run files without the reference
machine's WRF build trees, whose files name that machine's paths; `v05_files_data.tar.links` and `v05_suite_refs.tar.links`
restore four links between ideal reference directories that pointed into that machine's home, unpacked after a job's
tars), and the full CPU FP64 suite (wrf_suite.py's CPU cases, the output files tests, the
ideal files test, the differential and ctest tests, the translator tests and the comparators' self-test) runs on free
GitHub-hosted runners with the workflow `suite.yml` of that repository. The references are Intel-made; on an AMD CPU
the model gives the same results (`CPPWRF_PORTABLE_LIBMVEC`, see below) and so does WRF's side of the differential
tests (`tests/phys/libmvec_intel.cpp`), so every runner (AMD EPYC or Intel Xeon) is used.

```bash
# the full suite of HEAD (or of any local commit: tools/ci/ci_submit.sh COMMIT)
tools/ci/ci_submit.sh
# some cases (exact names), with the translator tests and the self-test; --full adds the other jobs
tools/ci/ci_submit.sh HEAD --cases real_r1_fp64_cpu,real_n3_r8_fp64_cpu
```

- `tools/ci/ci_submit.sh [COMMIT]` (default HEAD) starts the workflow and waits for it. A commit that is not on
  origin/main travels as a git bundle (only the objects origin/main lacks) attached to the `jobs` release of
  cppwrf-data and deleted after the run; nothing is pushed to this repository. The result files (logs, no data) go to
  `~/cppwrf_scratch/ci/<run id>/` and the script prints `tools/gcp/suite_status.sh`'s summary (`SUITE PASS n FAIL n
  SKIP n`, `files:`, `idealfiles:`, `diff:`, `ctest:`); exit 0 only for a run without FAIL or SKIP whose gate passed.
  Options: `--cases`, `--full`, `--data TAR,...` (instead of the cases' own tars), `--slow`, `--suite-jobs N` (default
  13), `--data-tag`, `--out`, `--hours`, `--no-wait`. It needs `gh` logged in with the right to start the workflow
  and upload release assets of cppwrf-data; no token is stored.
- Every push to main runs the full suite too (`.github/workflows/wrf_suite.yml` calls `suite.yml`).
- The workflow: `prep` takes the source and plans the jobs (`ci/plan.py`: the suite cases packed into jobs by their
  recorded times, `ci/times.txt`, and the download time of their tars; each job gets only the tars its cases read,
  `wrf_suite.py`'s `DATA`, or, for the files and translator tests, the tars holding the run directories they read),
  `build` compiles the CPU FP64 binary once as `tools/gcp/suite_body.sh` does (`-march=skylake`, Release, CUDA off),
  the `work` matrix runs the jobs in a Debian 13 container (gcc 14, gfortran 14, glibc 2.41, NetCDF 4.9.3, Python
  3.13: the packages of the machines that made the references) with OpenMP threads = the runner's vCPUs, and `gate`
  merges the logs and fails unless every job passed and the PASS lines are exactly the expected case list (a lost
  case is a FAIL, never a silent SKIP). At most 20 jobs run at once (the free plan's limit). The diff job takes
  WRF's build trees from the `wrf-build` release (built by `ci/wrf_build.sh` as below) instead of the reference
  machine's. The full suite (13 suite jobs, 2 files jobs, idealfiles, diff, checks) takes about 45 min from start to
  gate; the longest job about 40 min, of which up to 20 min is the data download.
- WRF references: the workflow `wrfref.yml` of cppwrf-data builds WRF v4.8.0 from the release source as the
  references were built (`configure` option 32, serial gfortran / gcc, basic nesting, `-r8`; and with
  `-fdefault-double-8`), keeps each build as an asset of its `wrf-build` release, runs `tools/gcp/wrf_ref_run.py`
  lines (`NAME BASE_RUN r8|r8d8 section.key=value ...`) on given data tars and returns the run directories as assets of
  its `wrfref` release; `rerun` runs WRF again in copies of existing reference directories and compares every output
  file's MD5 with the reference's (a different file is compared variable by variable, `ci/nc_same.py`). wrf.exe runs
  with CPPWRF's `tests/phys/libmvec_intel.cpp` preloaded (glibc's vector log / log10 with Intel's bits: on an AMD
  runner WRF's own libmvec calls change the results). The whole references are in `v05_suite_refs.tar`
  (`suite_data.tar` keeps only the compared variables). R1 1 h came back with the reference's MD5 on an AMD EPYC
  runner; P1 1 h equal in every variable but RQCBLTEN and RQIBLTEN (which the suite does not compare):
  `gh workflow run wrfref.yml -R Xassemblianist/cppwrf-data -f tag=NAME -f data=suite_data.tar,v05_suite_refs.tar -f rerun="wrf_runs/run_R1_r8:r8"`.
  A new reference is then staged into a new data tar (a later port task) and added to `vm_data.txt` and the release.
- Not on GitHub: the slow cases' data (`long_data.tar`, `long72_data.tar`, `long72n3_data.tar`) are not in the
  data release yet, and GPU runs need a GPU machine.

Colab is the other way (it needs a Colab subscription; the data come from the project's Drive):
`tools/colab/colab_suite.sh` (a Colab runtime, the same steps `tools/gcp/suite_body.sh` in a Debian 13 chroot; `--gpu G4 --any-cpu --jobs 12 --threads 4`
ran the full suite in about 42 min; GPU runs need it), with `--cases`, `--data`, `--slow`, `--diff-only [--faults DIR]`;
`tools/colab/colab_wrfref.sh` makes WRF references there,
`tools/colab/colab_mvec.sh` proves the portable vector log / log10 on a runtime of any CPU and
`tools/colab/colab_regression.sh` runs the legacy regression tool (`tools/regression`). The Google Cloud tools
(`tools/gcp/gcp_run.sh`, `gcp_wrfref.sh`) are kept but not used.

The suite's cases are defined in `tests/verify/cases/` (one file per group, `NN_<group>.py`, loaded in the order of
their names; `wrf_suite.py --list` prints them all). A group file defines `CASES`, `DATA` (the data tars each case
reads, in `vm_data.txt`'s order) and `STATUS` (its rows of the case table in docs/STATUS.md, which
`wrf_suite.py --status-table docs/STATUS.md` writes). `tests/verify/case_data.py LISTDIR [--write]` derives `DATA`
from the tars' file lists (`tools/colab/colab_tarlist.sh` reads them on a Colab runtime) and checks or writes it.

---

## 11. Performance

docs/PERF.md compares CPPWRF with WRF v4.8.0 on one Colab G4 runtime (AMD EPYC 9B45 with 48 hardware threads and an
RTX PRO 6000 Blackwell). For the real-data domain 119 x 89 x 44 with full physics (P4, 1 h), WRF real*4 at its
fastest MPI layout steps in 7.8 s. CPPWRF steps in 1.7 s on the GPU in FP32 and 2.8 s in FP64 (WRF real*8: 9.7 s).
On a 359 x 269 x 44 domain the GPU's lead grows to 8.2x in FP32 and 4.7x in FP64. `tools/bench/colab_bench.sh`
repeats the measurement. The memory-traffic budget of the large step is also in docs/PERF.md.

NetCDF-4 output (wrfout, wrfrst, auxhist3) is compressed as WRF compresses it (chunks of half the horizontal
dimensions, shuffle + deflate level 2), but the chunks of each record are compressed on a pool of host threads
before HDF5 writes them, with the same bytes as HDF5's own filter. `CPPWRF_NC_THREADS=N` sets the pool size (default:
the hardware threads, at most 64; 1: no pool; 0: HDF5's own filter) and `CPPWRF_IO_TIMING=1` prints the time of each
record (values, device copies, HDF5 calls). On the 359 x 269 x 44 domain a FP64 history frame takes about 1.2 s
instead of 11 s on the G4.

GPU memory: the fields, the physics tables, the work buffers (the column physics' per-column arrays and RRTMG's
batch scratch: together at most 1/4 of the device memory, less when little is free) and the kernels' local memory,
which CUDA reserves for every thread the GPU can hold at once (the largest frame times 30720 threads on an RTX 2060,
40960 on a T4, 89088 on an L4, 288768 on an RTX PRO 6000). The peak of a 1 h run on a 120 x 90 x 45 domain, measured
on a 16 GB T4 with 5800 MiB left free (what a 6 GB card offers) and `CPPWRF_DEVICE_MEMORY_MIB=5800`:

| Physics (mp / cu / pbl / lsm / radiation) | FP64 MiB | FP32 MiB |
|---|---|---|
| none / none / YSU / Noah / none (P1) | 1378 | 768 |
| none / none / MYJ / Noah / none (P9) | 1590 | 984 |
| none / none / YSU / Noah-MP / none (NMP1) | 1414 | 800 |
| none / Kain-Fritsch / YSU / Noah / none (G2) | 1738 | 1156 |
| none / Tiedtke / YSU / Noah / none (T6) | 1954 | 1058 |
| none / New Tiedtke / YSU / Noah / none (P7) | 2096 | 1182 |
| none / Grell-Freitas / YSU / Noah / none (P8) | 2102 | 1190 |
| WSM6 / Kain-Fritsch / YSU / Noah / RRTMG (P4) | 3098 | 2182 |
| aerosol-aware Thompson / none / YSU / Noah / none (C01) | 3150 | 1786 |
| WRF's CONUS suite: Thompson / Tiedtke / MYJ / Noah / RRTMG (D0) | 4052 | 2882 |
| Thompson / none / MYNN-EDMF / Noah / RRTMG (MYNN) | 4094 | 2688 |
| P4's physics on a two-way nest of two domains (N3) | 4646 | 3226 |

A larger domain adds its fields; a larger GPU adds local memory and work buffers (the same D0 FP64 run peaks at
8.7 GB on an RTX PRO 6000, which has 97 GB). On an RTX 2060 (6 GB) five steps of P4 and C01 in FP32 peak at 2.3 and
2.0 GB. When the memory is short the run stops at the set-up or the first step and names the kernel or the
allocation. `CPPWRF_DEVICE_MEMORY_MIB=N` sizes the work buffers for a device of N MiB (to reproduce a small GPU's
footprint on a large one, or to leave memory to another process).

`tools/heavy.sh` is a small helper that serialises heavy jobs (builds, runs, test suites) on a shared
workstation: `tools/heavy.sh <command>` runs in one of four job slots, `--light` runs at once with two
OpenMP threads, and `--exclusive` waits for an idle machine (use it for timings).
