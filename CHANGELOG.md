# Changelog

All notable changes of CPPWRF. Versions follow [Semantic Versioning](https://semver.org/): a change of the results of
a verified configuration, of a configuration key or of an output file is a new major version.

## Unreleased

### Added
- docs/FP32_ACCURACY.md "Root causes of the first-hour outliers": ACM7, MORR and SKV's FP32 excess at 1 h traced
  step by step to Kain-Fritsch deep / shallow and trigger decisions that are ties at the last bit (e.g. `TU95 >
  TV0` at 1 ulp in the first step), measured over five builds of the same code (GPU with and without FMA, CPU on
  two CPUs), perturbed starts and 24 h runs: no FP32 port defect, no output change. tools/bench/colab_fp32.sh
  --extra (GPU FP32 without FMA, CPU FP32) and --data; fp32_run.py's cases P4S / ACM7S / MORRS / SKVS, ACM7L / MORRL
  / SKVL and the perturbed-start E cases (data fp2_data.tar).
- BouLac PBL (bl_pbl_physics = 8), tier **verified**: WRF's module_bl_boulac.F (pbl_height, the length
  scales up / down of Bougeault and Lacarrere, the TKE equation and the implicit diffusion of u, v, theta, qv and qc
  with the Troen-Mahrt counter-gradient term) for one column, over the revised MM5 or the Eta surface layer, without
  the BEP urban terms; TKE_PBL, EL_PBL, EXCH_H / EXCH_M and the restart fluxes WU_TUR .. WQ_TUR as WRF writes them; no
  qi tendency (phy_bl_ten). Bit for bit with WRF's compiled routine in tests/phys/bl_boulac (ctest
  test_phys_bl_boulac) and suite cases BL1 / BL4 / BL9 against WRF real*8 (CPU and CUDA FP64 at tolerance 0).
  Shin-Hong (bl_pbl_physics = 11) stays an error.
- Verification tiers of options: `data/option_tiers.txt` lists option values merged before their whole-model WRF
  comparison as column_verified (bit-identical to WRF's Fortran in a column driver) or unverified (ported only).
  The model and cppwrf_real stop on such a value, naming option, value and tier, unless `[cppwrf] allow_unverified =
  true` (namelist.input: `&cppwrf allow_unverified = .true. /`, passed on by tools/wrf2cppwrf.py, which warns on
  them); with the switch each one is logged and the wrfout / wrfrst files get the global attribute
  `CPPWRF_UNVERIFIED_OPTIONS`. The table has no row: every option is verified, and no output byte changes. A new
  option's entry in this changelog names its tier.
- WSM3 and WSM5 microphysics, `mp_physics = 3` and `4` (phys/module_mp_wsm3.F, module_mp_wsm5.F): wsm3init /
  wsm5init, the semi-Lagrangian sedimentation, WRF's libmvec lanes of the vectorised saturation and diffusion
  loops, SNOWNC / SR, the effective radii for RRTMG (has_reqc = has_reqi = has_reqs = 1); real data. Tier:
  column_verified (tests/phys/mp_wsm35, bit for bit with the -r8 build's wsm3 / wsm5; suite cases W3 / W5 pending
  their WRF references). Not ported (errors): `do_radar_ref = 1` with them, the ideal cases, WSM7 (24), WSM6R (106).
- P3 microphysics with one ice category, `mp_physics = 50` (constant droplet number) and `51` (predicted droplet
  number QNCLOUD), phys/module_mp_p3.F: rime mass and volume QIR / QIB, TH_OLD / QV_OLD in the restart files, P3's
  REFL_10CM and effective radii for RRTMG; real data, one domain. Tier: verified (suite cases P3MP and P3NC against
  WRF real*8 of the -fdefault-double-8 build; column test tests/phys/mp_p3). Milbrandt-Yau (9), P3 52 / 53, P3 with
  nests, ideal cases or `do_radar_ref = 1` stay errors.
- Betts-Miller-Janjic cumulus (cu_physics = 2, WRF's module_cu_bmj.F): BMJINIT's look-up tables and SPLINE, BMJ's
  deep and shallow convection with the cloud efficiency CLDEFI (restart field), BMJDRV with its convective clouds
  (CCLDFRA, QCCONV, QICONV, CONVCLD; restart fields), phy_cu_ten (theta and qv), advance_ppt's HTOP / HBOT, and
  the adaptive time step's cudt = 0 (check_a_mundo); gfortran's libmvec lanes of the four vectorised loops.
  **Tier verified**: suite cases P9BMJ (3 steps, 1 h, restarted, cudt 5) and P4BMJ against WRF real*8 (CPU FP64
  and CUDA FP64 at tolerance 0); column driver tests/phys/cu_bmj (ctest test_phys_cu_bmj, label diff).
- bmj_rad_feedback (per domain, with BMJ): BMJ's clouds in the radiation (module_radiation_driver.F: QCCONV / QICONV
  added to qc / qi, CCLDFRA combined with CLDFRA); needs a microphysics scheme with cloud water and ice when the
  radiation is on. **Tier verified** (suite case P4BMJF, CPU and CUDA FP64 at tolerance 0).

### Fixed
- The restart alarm of a run resumed from a restart file (WRF's or CPPWRF's own) with a restart interval that is not
  a multiple of the time step: the first ring is at the restart time + WRF_ALARM_SECS_TIL_NEXT_RING_51 and the later
  ones are counted from that ring's step (share/input_wrf.F, ESMF_ClockAdvance), so the wrfrst files and their alarm
  attributes are WRF's; with restart_interval = 0 the resumed run writes the one set WRF writes at that ring. CPPWRF's
  restart_<step>.nc files are format version 2 (they carry WRF_ALARM_SECS_TIL_NEXT_RING_51). New suite cases RSA
  (P4, dt 50 s, restart every 2 min, the wrfrst files compared too) and its resumed runs, CPU and CUDA FP64 at
  tolerance 0 (data port_af2_data.tar); ctest restart_alarm.
- auxinput4_interval (sst_update) and gfdda_interval_m (grid_fdda) that are not multiples of the domain's time step:
  WRF stops at the first read ("Could not find matching time in input file"); tools/wrf2cppwrf.py and the model now
  refuse them (before, the model ran on with the first record). ctest input_interval.
- WRF's clock values on the fixed time step: the physics' xtime (radiation's coszen and hour angle, grid nudging),
  julian and the year (radiation, tmn_update, Noah-MP), solve_em's curr_secs2 / curr_mins2 (buckets, prec_acc_dt,
  grid nudging), clwrf's xtime and the history files' XTIME now come from the exact clock, n times the domain's
  rational time step (a nest's: the parent's / parent_time_step_ratio), as WRF's ESMF clock gives them and as the
  adaptive time step already did. Before, they were computed from (itimestep - 1) * dt, which differs in the last
  bit past one day (XTIME = days*1440 + seconds/60; e.g. dt 16 s with radt 10 min, from the radiation call of
  27.7 h on) and from the first steps on with a nest step that is a fraction of a second (a parent dt not divisible
  by parent_time_step_ratio). New suite cases P4D16 (30 h, dt 16 s, radt 10 min) and N3D50 (N3 with d01 dt 50 s,
  nest dt 50/3 s), CPU and CUDA FP64 at tolerance 0 (data port_af1_data.tar). Results of the suite's earlier cases
  are unchanged.

### Performance
- GPU memory: the column arrays of Tiedtke, New Tiedtke, Grell-Freitas-Li and GWDO (FP32 and FP64) and of
  Kain-Fritsch and MYJ (FP64) live in the column work buffer instead of the kernels' stack frames, which CUDA
  reserves for every thread the GPU can hold. The frames of 20-80 KB per thread are gone; the largest left is
  YSU's / Noah's (11-12 KB in FP64). A 1 h run on 120 x 90 x 45 peaks 4.7-15 GB lower on an RTX PRO 6000 with these
  schemes (WRF's CONUS suite in FP64: 20.7 -> 8.7 GB) and 1.1-3.1 GB lower on an L4; every verified suite fits a
  6 GB GPU in FP32 and FP64 (docs/USER_GUIDE.md "GPU memory" has the minimum per suite, docs/PERF.md the tables).
  The moved kernels are up to 11 % (FP64) and 30 % (FP32) slower on the RTX PRO 6000; the steps of Grell-Freitas
  and Tiedtke runs 3-7 %. Results are unchanged: CUDA FP64 at tolerance 0, FP32 byte-identical.
- `CPPWRF_DEVICE_MEMORY_MIB=N`: the work buffers (column physics, RRTMG) are sized as on a device of N MiB.
- tools/colab/colab_run.sh `--base REV` (both builds in one session: times, memory, kernel times, outputs byte for
  byte) and `--device-mib N` (the runs see N MiB of free device memory).

### Added
- Domain decomposition with divided memory, and MPI (docs/MPI.md, USER_GUIDE section 7): a domain and its nests run
  as P x Q patches, each worker holding its patch and a halo of 3 points under the domain's indices. The workers
  are tiles (`run.tiles`, as before, now also with nests) or MPI ranks: a build with `-DCPPWRF_MPI=ON` (off by
  default) under `mpirun -np N`, layout `run.ranks = PxQ` / `CPPWRF_RANKS` (default the most square with P <= Q),
  `CPPWRF_RANK_DEVICES` for the GPUs of a CUDA build. One worker writes every output file. The result has the same
  bits for every worker count and layout (CPU FP64 = one-process WRF real*8): proved on one machine for 1 to 8
  ranks in 13 layouts on eight cases, 53 suite cases in three layouts, and up to 48 ranks against the one-process
  run; not run on several GPUs or nodes. A real-data start reads only the worker's window of wrfinput's 3-D
  variables. Results without tiles and outside `mpirun` are unchanged.
- tools/colab/colab_mpi.sh (suite cases under `mpirun` with several layouts on a Colab runtime, CPU and GPU, each
  against the WRF references and byte for byte against one rank) and tools/bench/colab_bench.sh `--scaling` (strong
  scaling next to WRF real*8 dmpar).

### Changed
- Tiles no longer write `output.dir.tileN` copies: tile 0 writes the run's files from gathered fields.
- An MPI run prints the peak host memory of the writer, of the largest other rank and of all ranks.
### Verification
- New suite cases against WRF v4.8.0 real*8 (tests/verify/cases/69_coverage.py, data port_cov1_data.tar and
  port_cov1b_data.tar; 3 steps and 1 h on the CPU, 1 h on the GPU, every variable at tolerance 0) for option values
  that had no suite case: RRTMG `cldovrlp` 1, 3, 4 (with `idcor` 1) and 5 (OVL1, OVL3, OVL4I, OVL5), `icloud` 0
  (ICL0), `diff_6th_opt` 1 (D6O1), `km_opt` 3 with `diff_opt` 1 (KM3H), `tke_adv_opt` 0 with `km_opt` 2 (TKA0),
  `moist_adv_opt` / `scalar_adv_opt` 0 (ADV0) and `bl_mynn_cloudpdf` -1 / -2 (MPDFM1, MPDFM2). `km_opt` 3 with
  `diff_opt` 2 (KM3) is compared over two steps only: on this 15 km domain WRF itself blows up in the third step
  (NaN in every field, which wrf.exe does not stop on); CPPWRF is identical up to that point and then stops with its
  non-finite check. All verified tier.
- Verification coverage (docs/STATUS.md): the audit's order check (serial reversed and i-outer launch order against
  the OpenMP run) and the pattern-initialised, bounds-checked build on the configurations whose inputs were not on
  the workstation, and an NVIDIA compute-sanitizer pass (memcheck, racecheck, initcheck, synccheck; FP64 and FP32)
  over P4, MYNN, Thompson, Morrison, RRTMG, the two-way nest, a 2x2 tiled run and an ideal case.

### Tools and documentation
- The CPU FP64 suite on free GitHub-hosted runners (docs/USER_GUIDE.md, "The CPU Suite on GitHub"): the reference
  data (the tars of tests/verify/vm_data.txt) are public release assets of github.com/Xassemblianist/cppwrf-data,
  whose workflow suite.yml runs wrf_suite.py's CPU cases, the output files tests, the ideal files test, the
  differential and ctest tests, the translator tests and the comparators' self-test in a Debian 13 container, a
  gate failing unless every expected case PASSed. tools/ci/ci_submit.sh runs it for any local commit (a git bundle
  of the commits origin/main lacks) and prints suite_status.sh's summary; .github/workflows/wrf_suite.yml runs it on
  every push to main. Its workflow wrfref.yml builds WRF v4.8.0 as the references were built (wrf.exe with the Intel
  vector log / log10 of tests/phys/libmvec_intel.cpp preloaded) and makes WRF references or runs existing ones again
  (R1 1 h reproduced byte for byte on an AMD runner). Colab (tools/colab) stays the other way, and the one for GPU
  runs.
- `cppwrf-post` (docs/POST.md): forecaster fields from wrfout files with wrf-python's methods (ported to numpy,
  identical to wrf-python on single precision files), written as CF-1.8 NetCDF and / or GRIB2 (Lambert, polar
  stereographic, Mercator and latitude-longitude grids; WMO parameters; simple packing, JPEG 2000 / CCSDS through
  ecCodes); a post.yaml / post.ini configuration.
- `cppwrf-verify-obs` (docs/VERIFY_OBS.md): a run against surface observations (IEM METARs, NCEI ISD or a CSV
  file), with continuous, wind and precipitation contingency scores per lead time and station, SVG plots, a
  Markdown summary, and the comparison of two runs.
- docs/FP32_ACCURACY.md: how close GPU FP32 is to WRF real*8 compared with stock WRF real*4 (R = FP32's difference /
  WRF real*4's, measured on eleven domains at 1 h and over 24 h / 72 h), its method and when to use FP32 or FP64.
- tools/verify/fp32_accuracy.py: max |d|, RMSE, relative L2 and bias of any two sets of wrfout files per frame and
  variable (surface fields, accumulated and hourly rain with category scores, 3-D state, hydrometeors), and the report
  of the ratio between pairs with growth figures (SVG).
- tools/bench/colab_fp32.sh (with fp32_body.sh, fp32_gpu.sh, fp32_run.py): the study on one Colab runtime: WRF
  v4.8.0 real*4 and real*8 built from the release source, CPPWRF CUDA FP32 and FP64, twelve configurations
  (1 h; 24 h and 72 h from the suite's long reference inputs), statistics on the runtime.
- The FP32 cases real_acm7_fp32_1h, real_morr_fp32_1h and real_skv_fp32_1h have measured tolerances (about 3x the
  values of the study); real_skv_fp32_1h no longer lists TSK, which tests/verify/compare_real.py does not compare.

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
- Microphysics: Kessler (1), WSM6 (6), Thompson (8), Morrison two-moment (10, `morr_rimed_ice` 0 / 1),
  aerosol-aware Thompson (28); radar reflectivity (`do_radar_ref`).
- Cumulus: Kain-Fritsch (1), Grell-Freitas-Li (3), Tiedtke (6), New Tiedtke (16); cumulus clouds in the radiation
  (`cu_rad_feedback`).
- PBL and surface layer: YSU (1), MYJ (2), MYNN-EDMF (5) with its options, ACM2 (7); revised MM5 (1), Eta (2),
  MYNN (5) and Pleim-Xiu (7) surface layers; fog deposition (`grav_settling`); gravity wave drag (`gwd_opt` 1).
- Land surface: 5-layer thermal diffusion (1), Noah (2), Noah-MP v5.1.2 (4, default options); fractional sea ice,
  SST update (`sst_update`, wrflowinp).
- Radiation: RRTMG LW and SW (4), RRTM LW (1), Dudhia SW (1); McICA (`cldovrlp` 0-5), CAM ozone and greenhouse
  gases.
- Nests: static one-way and two-way nests (feedback, smoothing, terrain blending), several nests and nests in nests.
- Four-dimensional data assimilation: grid analysis nudging (`grid_fdda = 1`, wrffdda_d0N).
- Tiles: one domain as P x Q tiles (`run.tiles`) with the same bits as without tiles.
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
