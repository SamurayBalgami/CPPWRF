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

The model steps entirely on the GPU with almost no host round trips, and each time step runs as
a CUDA graph. The table gives seconds of computation per simulated hour on a real-data domain
(119 x 89 x 44 cells, 15 km, dt 60 s), from [docs/PERF.md](docs/PERF.md):

| Case | RTX PRO 6000 Blackwell | A100 | L4 |
|---|---|---|---|
| P1 (YSU, Noah, sfclay), FP32 | 0.25 s | 0.52 s | 0.61 s |
| P1, FP64 | 0.81 s | 0.65 s | 2.12 s |
| P4 (adds KF, WSM6, RRTMG), FP64 | 5.5 s | 4.1 s | 9.8 s |

Most consumer and workstation GPUs run FP64 at 1/64 of the FP32 rate, so FP32 is the fast mode on them.
The profile ranks the next speed-ups: the RRTMG longwave takes 69 % of P4, and the column physics use the
GPU at low occupancy. Every speed-up must leave the FP64 results bit-identical.

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
- [docs/PERF.md](docs/PERF.md): GPU profile and the performance roadmap
- [docs/PHYSICS_SUITE.md](docs/PHYSICS_SUITE.md), [docs/PHYSICS_DRIVER.md](docs/PHYSICS_DRIVER.md),
  [docs/PHYSICS_INTERFACE.md](docs/PHYSICS_INTERFACE.md): the physics port
- [docs/NESTING.md](docs/NESTING.md), [docs/NEST_NUMERICS.md](docs/NEST_NUMERICS.md): nesting
- [docs/REAL.md](docs/REAL.md): `cppwrf_real`

## License

GPL-3.0-or-later; see [LICENSE](LICENSE) and [NOTICE.md](NOTICE.md), which includes the WRF notice. WRF is
developed by NCAR/MMM and its partners. CPPWRF is an independent port and is not affiliated with or endorsed
by NCAR or UCAR.
