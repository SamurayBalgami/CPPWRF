# CPPWRF 1.0.0 release notes

CPPWRF 1.0.0 is a C++20 / CUDA port of the WRF-ARW model v4.8.0. It reads WRF's inputs and writes WRF's output
files, and its CPU FP64 and GPU FP64 builds give WRF v4.8.0 real*8's results bit for bit in every verified
configuration. CHANGELOG.md lists the contents; docs/USER_GUIDE.md explains the use; docs/STATUS.md lists every
verified case.

## Verified scope

What "verified" means here: WRF v4.8.0 (gfortran 14.2, glibc 2.41, Debian 13, `configure -r8`, serial build, one
process) and CPPWRF run from the same inputs, and every value of every variable of every output time is compared at
tolerance 0 (NaN equals only NaN, -0 differs from +0).

| Check | Scope |
|---|---|
| Model runs against WRF (`tests/verify/wrf_suite.py`) | 192 cases: 113 CPU FP64 and 64 GPU FP64 cases at tolerance 0, 15 GPU FP32 cases by relative L2; 8 of them are the 24 h and 72 h runs of `--slow` |
| Restarts | 33 checks: a run restarted from its own restart files equals the uninterrupted run (15 CPU and 15 GPU FP64 configurations, also one FP32), and the Straka checks |
| WRF output files | the wrfout and wrfrst files of the CPU FP64 runs equal WRF's in values, variables, dimensions, attributes and global attributes (31 runs), also when resumed from WRF's wrfrst files; the wrfout files of WRF's ideal cases (11 runs) |
| GPU FP64 | the 64 GPU FP64 cases and 15 GPU restart checks (79 runs) at tolerance 0 on NVIDIA L4 (sm_89), RTX PRO 6000 Blackwell (sm_120) and A100 (sm_80) |
| Schemes against WRF's own objects | the differential tests (`ctest -L diff`) link each scheme's driver to the objects of the WRF build that made the references: Kessler, WSM6, Thompson, YSU, MYJ, MYNN, the surface layers, Noah, Noah-MP, Kain-Fritsch, diffusion, nest numerics, the physics driver, Grell-Freitas and New Tiedtke column replays |
| real.exe | `cppwrf_real` output equals real.exe's for 43 cases (single domains, nests, `sst_update`, hybrid and hypsometric options) |
| Namelist variants | several hundred 3-minute and 1-hour variants of the namelist options against WRF (docs/STATUS.md, "Namelist variants") |

The configurations cover real-data domains in the mid-latitudes (15 km), the tropics (Bay of Bengal, Mercator), the
Arctic (Svalbard, polar stereographic, sea ice) and the Alps (3 km); 24 h and 72 h forecasts with restarts in the
middle; WRF's default `namelist.input` (physics_suite CONUS); one-way and two-way nests; WRF's idealized cases. The
CPU result does not depend on the OpenMP thread count, and it is the same on Intel and AMD CPUs
(`CPPWRF_PORTABLE_LIBMVEC`).

The quick start of README.md and docs/USER_GUIDE.md (the CPU FP64, GPU FP32 and GPU FP64 builds, `ctest`, the
translator's tests, the Straka and em_quarter_ss cases and the end-to-end example from the GFS analysis through WPS
and `cppwrf_real` to wrfout files) runs from a clean source tree: `tools/colab/colab_quickstart.sh` did it on a
Google Colab G4 runtime, where the example's GPU FP64 output equalled the CPU FP64 output bit for bit.

## Conditions of bit identity

- CPU FP64 builds with the flags the CMake files set (no FMA contraction, no fast-math; docs/USER_GUIDE.md,
  "Required Flags and Determinism"). Any x86-64 CPU with SSE4.1, FMA and AVX2.
- The comparison is with WRF real*8 (`configure -r8`; `-r8 -fdefault-double-8` for Thompson, MYNN and radar
  reflectivity, whose DOUBLE PRECISION variables are 16-byte with `-r8` alone) built with gfortran 14.2 and glibc
  2.41 and run in one process. WRF's own result changes with its MPI process layout, its compiler version
  (another gfortran vectorises other loops) and `-march` (FMA contraction).
- GPU FP64 builds with `--fmad=false` and `CPPWRF_CUDA_GLIBC_MATH` (the CMake defaults).
- GPU FP32 is not bit-identical to anything: it differs from WRF's single-precision build by about as much as WRF's
  single and double precision builds differ.

## Performance

On one Google Colab G4 runtime (AMD EPYC 9B45 with 48 hardware threads, NVIDIA RTX PRO 6000 Blackwell), stepping
time per simulated hour against WRF v4.8.0 at its fastest MPI layout on the same machine (docs/PERF.md):

| Case | WRF real*4 | WRF real*8 | CPPWRF GPU FP32 | CPPWRF GPU FP64 |
|---|---|---|---|---|
| P4: 119 x 89 x 44, 15 km, full physics, radt 1 min | 7.84 s | 9.72 s | 1.66 s | 2.78 s |
| N3: P4 with a two-way 5 km nest | 23.8 s | 34.0 s | 5.25 s | 10.4 s |
| L: 359 x 269 x 44, radt 15 min | 24.8 s | 45.1 s | 3.01 s | 9.59 s |

## Not included

MPI runs (a run is one process: OpenMP threads or one GPU), moving nests, vertical nest refinement, ndown.exe, I/O
forms other than NetCDF, and the options listed in docs/USER_GUIDE.md ("Rejected Options and Why"): urban and lake
models, other microphysics, cumulus, PBL and land surface schemes, Noah-MP's non-default options, FDDA, stochastic
physics, chemistry, the adaptive time step with nests. Each is rejected with an error by the translator and the model.

## Requirements

A C++20 compiler, CMake >= 3.24, OpenMP, NetCDF-C, Python 3 with `netCDF4` and `numpy`; the CUDA toolkit >= 12
and an NVIDIA GPU for the GPU builds. WRF's `run/` directory (WRF v4.8.0 source) for the physics tables, and WPS or
WRF's real.exe output for real-data inputs.

## Licence

GPL-3.0-or-later (LICENSE). NOTICE.md carries the notices of WRF (UCAR), RRTMG (AER), Noah-MP (UCAR) and the GNU C
Library, from which parts of CPPWRF are derived.
