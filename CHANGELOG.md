# Changelog

All notable changes of CPPWRF. Versions follow [Semantic Versioning](https://semver.org/): a change of the results of
a verified configuration, of a configuration key or of an output file is a new major version.

## 1.0.0

The first release: a C++20 / CUDA port of WRF-ARW v4.8.0 whose CPU FP64 and GPU FP64 output is bit-identical to WRF
v4.8.0 real*8 for every verified configuration. The release notes ([docs/RELEASE_NOTES.md](docs/RELEASE_NOTES.md))
give the exact verified scope.

### Model
- Dynamical core: RK3 split-explicit time integration with acoustic substeps, hybrid vertical coordinate
  (`hybrid_opt` 0 / 2), `hypsometric_opt` 1 / 2, map projections (Lambert, Mercator, polar stereographic), terrain,
  Coriolis and curvature terms; 5th / 6th-order horizontal and 3rd / 4th-order vertical advection, positive-definite
  moisture, scalar and TKE advection; specified, open, periodic and symmetric lateral boundaries; `w_damping`,
  `damp_opt` 1-3; the adaptive time step for single-domain real-data runs.
- Diffusion: `diff_opt` 1 and 2 with `km_opt` 1-4 (constant K, 1.5-order TKE, 3-D and 2-D Smagorinsky),
  6th-order diffusion (`diff_6th_opt` 1, 2).
- Microphysics: Kessler (1), WSM6 (6), Thompson (8), aerosol-aware Thompson (28); radar reflectivity
  (`do_radar_ref`).
- Cumulus: Kain-Fritsch (1), Grell-Freitas-Li (3), Tiedtke (6), New Tiedtke (16); cumulus clouds in the radiation
  (`cu_rad_feedback`).
- PBL and surface layer: YSU (1), MYJ (2), MYNN-EDMF (5) with its options; revised MM5 (1), Eta (2) and MYNN (5)
  surface layers; fog deposition (`grav_settling`); gravity wave drag (`gwd_opt` 1).
- Land surface: 5-layer thermal diffusion (1), Noah (2), Noah-MP v5.1.2 (4, default options); fractional sea ice,
  SST update (`sst_update`, wrflowinp).
- Radiation: RRTMG LW and SW (4), RRTM LW (1), Dudhia SW (1); McICA (`cldovrlp` 0-5), CAM ozone and greenhouse
  gases.
- Nests: static one-way and two-way nests (feedback, smoothing, terrain blending), several nests and nests in nests.
- Diagnostics: accumulation buckets (`bucket_mm`, `bucket_J`), precipitation windows (`prec_acc_dt`),
  `output_diagnostics` (auxhist3).
- Pre-processing: `cppwrf_real` (WRF's real.exe: wrfinput, wrfbdy, wrflowinp from WPS met_em files) and the port of
  WRF's ideal.exe for `em_hill2d_x`, `em_quarter_ss`, `em_squall2d_x/y`, `em_grav2d_x`.
- I/O: WRF's wrfout, wrfrst and auxhist3 NetCDF files (the WRF v4.8.0 I/O schema, `iofields_filename`, file names,
  frames per file, NetCDF classic), NetCDF-4 compression on host threads with the same bytes as HDF5's filter;
  restarts from WRF's or CPPWRF's restart files; CPPWRF's own restart files record the format version, the commit
  and the options.
- Backends: OpenMP CPU and CUDA GPU from one source, FP64 and FP32; on the GPU, CUDA graphs for the time step.

### Tools
- `tools/wrf2cppwrf.py`: WRF `namelist.input` to CPPWRF configuration, with `physics_suite`, WRF's check_a_mundo
  resets and an error for every option without a port and every unknown key.
- `examples/quickstart/run.sh`: GFS -> WPS -> cppwrf_real -> cppwrf -> wrfout with WRF's default domains and namelist.
- Verification: `tests/verify/wrf_suite.py` (model runs against WRF at tolerance 0), the WRF output files tests,
  the differential tests against WRF's compiled Fortran (`ctest -L diff`), unit tests, and the Colab / Google Cloud
  drivers in `tools/colab` and `tools/gcp`; `tools/bench` (the benchmark against WRF v4.8.0 in docs/PERF.md).

### Not included
MPI runs, moving nests, vertical nest refinement, ndown.exe, I/O forms other than NetCDF, and the schemes and options
listed in docs/USER_GUIDE.md ("Rejected Options and Why"); each is rejected with an error.
