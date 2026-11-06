# CPPWRF

CPPWRF is a C++20 / CUDA port of the WRF-ARW v4.8.0 model: the dynamical core, the physics, static nesting,
`real.exe` and `ideal.exe`. One source tree builds for the CPU (OpenMP) and for NVIDIA GPUs (CUDA), in FP64
or FP32. It reads WRF's `wrfinput` / `wrfbdy` files and writes WRF's `wrfout` and `wrfrst` files.

The port follows WRF statement by statement, and every result is checked against WRF itself:

- **CPU, FP64: bit-identical to WRF v4.8.0 real*8.** Every compared value of every variable and output time
  is equal at tolerance 0, for the real-data cases, the physics suites, nests, restarts and WRF's idealized cases.
- **GPU, FP64: bit-identical to WRF too.** All 79 CUDA FP64 runs of the suite match WRF and the CPU build at
  tolerance 0. They were checked on three GPU generations: L4 (sm_89), RTX PRO 6000 Blackwell (sm_120) and
  A100 (sm_80).
- **GPU, FP32:** the fast mode, compared with WRF by relative L2 error.

The details are in [docs/STATUS.md](docs/STATUS.md).

## What is ported

| Area | WRF options |
|---|---|
| Dynamics | RK3 split-explicit, adaptive time step, hybrid coordinate, map projections (Lambert, Mercator, polar), terrain; specified, open, periodic and symmetric boundaries; damping, 5th/6th-order advection, positive-definite moisture |
| Diffusion | diff_opt 1 and 2; km_opt 1-4 (constant K, Smagorinsky, 1.5-order TKE); 6th-order diffusion |
| Microphysics | Kessler, WSM6, Thompson, aerosol-aware Thompson (mp 28); radar reflectivity |
| Cumulus | Kain-Fritsch, Grell-Freitas-Li, Tiedtke, New Tiedtke |
| PBL and surface layer | YSU, MYJ, MYNN-EDMF with most of its options; revised MM5, Eta and MYNN surface layers; gravity wave drag |
| Land surface | 5-layer slab, Noah, Noah-MP (default options); fractional sea ice, SST update |
| Radiation | RRTMG LW and SW, RRTM LW, Dudhia SW, cumulus clouds in radiation, McICA |
| Nesting | static one-way and two-way nests with feedback and smoothing |
| Pre-processing | `cppwrf_real` (real.exe from WPS met_em files) and WRF's ideal.exe cases |
| Diagnostics and I/O | buckets, precipitation windows, output_diagnostics (auxhist3), WRF NetCDF wrfout / wrfrst, iofields |

WRF's `namelist.input` works as is: `tools/wrf2cppwrf.py` translates it, including `physics_suite`. An option
without a port is reported as an error and never silently ignored. Still missing: MPI runs, moving nests,
urban and lake models, and the other schemes listed in [docs/STATUS.md](docs/STATUS.md#not-yet-implemented).

## How it is verified

- **Column and kernel tests** (`tests/phys`, `tests/dyn`, `tests/nest`): each scheme runs on the same inputs in
  C++ and in WRF's own compiled Fortran objects, and the outputs must be bit-identical.
- **Model runs against WRF** (`tests/verify/wrf_suite.py`): 192 cases plus 33 restart checks. They cover
  real-data domains in the mid-latitudes, the tropics (Bay of Bengal), the Arctic (Svalbard) and the Alps
  at 3 km. They include 24 h and 72 h forecasts and WRF's default CONUS namelist. CPPWRF's wrfout files are
  compared with WRF's in every variable and output time.
- **WRF output files:** wrfout and wrfrst match WRF's files in values, attributes and metadata. Runs resumed
  from WRF's restart files match too.

To be bit-identical, the port must repeat WRF's arithmetic exactly. That means the order of every sum, the
constants gfortran folds, the halo schedule, and the glibc vector math (libmvec) that gfortran calls. On the
GPU, CPPWRF runs a device port of the same glibc functions and turns FMA contraction off. The CPU result is
the same on Intel and AMD processors (`CPPWRF_PORTABLE_LIBMVEC`).

## Performance

WRF v4.8.0 and CPPWRF ran on one Google Colab G4 runtime: an AMD EPYC 9B45 (24 cores, 48 hardware threads) and an
NVIDIA RTX PRO 6000 Blackwell (96 GB). WRF was built there from the release source with its own GNU `configure`
(MPI, and MPI + OpenMP), in real*4 and real*8; its fastest of five process layouts counts. Each run is 1 simulated
hour with the same namelist and inputs. The table gives the stepping time in seconds per simulated hour, without
start-up and history output ([docs/PERF.md](docs/PERF.md) has the wall times, other layouts and memory):

| Case | WRF real*4, best layout | WRF real*8, best layout | CPPWRF GPU FP32 | CPPWRF GPU FP64 |
|---|---|---|---|---|
| P4: 119 x 89 x 44, 15 km, full physics, radt 1 min | 7.84 | 9.72 | 1.66 (4.7x) | 2.78 (3.5x) |
| P4 with radt 15 min | 2.38 | 4.34 | 0.55 (4.3x) | 1.68 (2.6x) |
| N3: P4 with a two-way 5 km nest | 23.8 | 34.0 | 5.25 (4.5x) | 10.4 (3.3x) |
| L: 359 x 269 x 44 (9x P4), radt 15 min | 24.8 | 45.1 | 3.01 (8.2x) | 9.59 (4.7x) |

The speed-ups in brackets are against WRF real*4 for FP32 and against WRF real*8 for FP64. The GPU's advantage grows
with the domain size. Whole runs gain less: start-up and NetCDF output are serial in both models. Without history
output GPU FP32 runs 1.8x to 3.5x faster than WRF real*4 in wall time. CPPWRF compresses the NetCDF-4 chunks of
its history and restart files on host threads, and the files keep every byte: with a frame every
10 minutes and a restart file, L's hour takes 13 s in GPU FP32 (58 s with serial compression) and 33 s in FP64
(123 s).

In these runs, the GPU FP64 results are bit-identical to the CPU FP64 results. The CPU FP64 results are
bit-identical to WRF real*8 run in one process, which is how the references are made: WRF's own result changes
with its MPI layout. GPU FP32 differs from WRF real*4 by about as much as WRF real*4 differs from WRF real*8.

## Build

```bash
# GPU build, FP32
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=ON -DCPPWRF_DOUBLE=OFF
cmake --build build -j

# CPU build, FP64 (the configuration that is bit-identical to WRF real*8)
cmake -S . -B build-cpu -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_DOUBLE=ON
cmake --build build-cpu -j
```

You need a C++20 compiler, CMake >= 3.24, OpenMP and NetCDF-C. The CUDA toolkit >= 12 is optional.

## Run

From a WRF run directory:

```bash
python3 tools/wrf2cppwrf.py namelist.input -o run.ini --format wrfout
build/cppwrf run.ini
```

To make `wrfinput` and `wrfbdy` from WPS met_em files instead of with real.exe
([docs/REAL.md](docs/REAL.md)):

```bash
build/src/init_real/cppwrf_real --namelist namelist.input --met-dir /path/to/met_em --precision double
```

The example configurations in `configs/` include idealized cases (`build/cppwrf configs/straka.ini`) and the
suite's real-data cases. The real-data cases read the WRF run directories from `$CPPWRF_WRF_RUNS`.

## Documentation

- [docs/USER_GUIDE.md](docs/USER_GUIDE.md): building, running, configuration keys, output formats
- [docs/STATUS.md](docs/STATUS.md): what is ported, how it is verified, every case and result
- [docs/PERF.md](docs/PERF.md): the benchmark against WRF, the GPU profile and the performance roadmap
- [docs/PHYSICS_SUITE.md](docs/PHYSICS_SUITE.md), [docs/PHYSICS_DRIVER.md](docs/PHYSICS_DRIVER.md),
  [docs/PHYSICS_INTERFACE.md](docs/PHYSICS_INTERFACE.md): the physics port
- [docs/NESTING.md](docs/NESTING.md), [docs/NEST_NUMERICS.md](docs/NEST_NUMERICS.md): nesting
- [docs/REAL.md](docs/REAL.md): `cppwrf_real`

## License

GPL-3.0-or-later; see [LICENSE](LICENSE) and [NOTICE.md](NOTICE.md), which includes the WRF notice. WRF is
developed by NCAR/MMM and its partners. CPPWRF is an independent port and is not affiliated with or endorsed
by NCAR or UCAR.
