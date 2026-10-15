# CPPWRF

CPPWRF is a C++20 / CUDA reimplementation of the WRF-ARW v4.6.1 dynamical core, physics and static
nesting, written from one source for CPU (OpenMP) and GPU (CUDA), in FP32 or FP64. It reads standard
`wrfinput_d0N` / `wrfbdy_d01` files and writes `wrfout_d0N` and restart files.

Every routine follows WRF statement by statement and is verified against WRF itself: on the CPU in FP64,
CPPWRF is bit-identical to WRF v4.6.1 real*8 for the real-data dynamics, the physics suites (surface
layer, Noah LSM, YSU, Kain-Fritsch, Grell-Freitas, New Tiedtke, WSM6, Thompson, RRTMG) and one- and
two-way nesting. The GPU build agrees to roundoff in FP64. See [docs/STATUS.md](docs/STATUS.md).

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=ON -DCPPWRF_DOUBLE=OFF
cmake --build build -j
build/cppwrf configs/straka.ini
```

Requirements: a C++20 compiler, CMake >= 3.24, OpenMP, NetCDF-C; optionally the CUDA toolkit >= 12.

## Running real-data cases

Configuration files are INI files (see [docs/USER_GUIDE.md](docs/USER_GUIDE.md)); values may refer to
environment variables as `${NAME}`. The real-data configs in `configs/` expect the WRF run directories
(wrfinput, wrfbdy and the reference wrfout files) under `$CPPWRF_WRF_RUNS`:

```bash
export CPPWRF_WRF_RUNS=/path/to/wrf_runs
build/cppwrf configs/real_p4_r8.ini
```

## Verification

`tests/verify/wrf_suite.py` runs every configuration that has a WRF reference and compares the output
field by field; `tests/phys`, `tests/dyn` and `tests/nest` hold column and kernel tests against the WRF
Fortran. Documentation: [docs/PHYSICS_SUITE.md](docs/PHYSICS_SUITE.md),
[docs/NEST_NUMERICS.md](docs/NEST_NUMERICS.md), [docs/NESTING.md](docs/NESTING.md).

## License

GPL-3.0-or-later, see [LICENSE](LICENSE) and [NOTICE.md](NOTICE.md) (includes the WRF notice).
