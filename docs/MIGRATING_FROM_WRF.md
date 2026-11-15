# Migrating from WRF to CPPWRF

CPPWRF is a port of WRF-ARW v4.8.0, not a new model. It reads your `namelist.input` and WRF's input files and writes
WRF's output files. The options it ports give WRF's numbers. This page says what stays the same, which options are
ported, what differs, and how to check CPPWRF against your own WRF runs.

## What is identical

- **Inputs:**
  - `namelist.input`, unchanged;
  - `wrfinput_d0N` and `wrfbdy_d01` from WRF's real.exe or from `cppwrf-realexe`;
  - `wrflowinp_d0N` (`sst_update = 1`) and `wrffdda_d0N` (`grid_fdda = 1`);
  - WRF's restart files (`wrfrst_d0N_<date>`);
  - the tables of WRF's `run/` directory;
  - the `input_sounding` of WRF's idealized cases.
- **Outputs:** `wrfout_d0N_<date>`, `wrfrst_d0N_<date>` and `auxhist3_d0N_<date>` (`output_diagnostics`), with
  WRF's variables, dimensions, attributes and global attributes in WRF's order. `history_outname`, `rst_outname`,
  `frames_per_outfile`, `iofields_filename` and `use_netcdf_classic` work as in WRF. Post-processing that reads
  wrfout files (wrf-python, NCL, UPP, ncview) reads CPPWRF's.
- **The numbers:**
  - `cpu-fp64` and `gpu-fp64` give every value of every variable at every output time exactly as WRF v4.8.0 real*8
    does, at tolerance 0. That holds for every configuration of the verification suite (docs/RELEASE_NOTES.md).
  - The comparison is with WRF built with gfortran 14.2 on Debian 13 and run in one process. WRF's own result
    changes with its MPI layout, its compiler and its flags. So WRF on your cluster, with 64 MPI ranks, does not
    match itself on one rank bit for bit either.
  - `gpu-fp32` is the fast mode. It differs from WRF real*8 by about as much as WRF's own real*4 build does.
- **The run directory:** `cppwrf-wrf` and `cppwrf-realexe` are drop-in replacements for `wrf.exe` and `real.exe`
  (below). They keep WRF's file names and places, the `rsl.out.0000` / `rsl.error.0000` logs with `Timing for main`
  lines, `SUCCESS COMPLETE WRF`, and a nonzero exit status on failure.

## From wrf.exe to cppwrf-wrf

| WRF | CPPWRF |
|---|---|
| `./real.exe` (or `mpirun -np N ./real.exe`) | `cppwrf-realexe` |
| `mpirun -np N ./wrf.exe` | `cppwrf-wrf --mode cpu-fp64 --threads N`, or `cppwrf-wrf --mode gpu-fp32` / `gpu-fp64` |
| `grep SUCCESS rsl.out.0000` | the same |
| `restart = .true.` with `wrfrst_d01_<date>` in the directory | the same |
| a namelist option WRF has | `cppwrf-check namelist.input`: ported, ignored or not ported, per option |

In a batch script, the change is one line:

```bash
#SBATCH --gres=gpu:1
cd $RUNDIR
cppwrf-wrf --mode gpu-fp32 || exit 1      # was: srun ./wrf.exe
```

CPPWRF's own files go to `cppwrf_work/` in the run directory: the translated configuration `run.ini`, `run.json` and
`diag.csv`. They are not needed afterwards. `cppwrf-wrf --dry-run` shows the translation and the command without
running.

## What differs

- **Parallelism.** A CPPWRF run is one process: OpenMP threads on the CPU, or one GPU. There is no MPI, so a domain
  must fit in one machine's memory, or in one GPU's (`cppwrf-fit` estimates it). The CPU result does not depend on
  the thread count.
- **Options without a port** stop the run at start-up with the option, the value and what is supported. They are:
  - MPI runs and moving nests;
  - vertical nest refinement, ndown.exe;
  - the urban and lake models;
  - spectral and observation nudging;
  - stochastic physics and chemistry;
  - I/O forms other than NetCDF;
  - the other schemes not listed below.

  docs/STATUS.md keeps the list current.
- **Restart files.** CPPWRF reads and writes WRF's `wrfrst` files. It can also write its own `restart_<step>.nc`
  (`cppwrf-wrf --native-restart`), which keep every bit of its state.
- **Before the drop-ins** (CPPWRF up to 1.0.0), a run needed `tools/wrf2cppwrf.py namelist.input -o run.ini` and
  wrote into `out/`, with a nested run's files in `out/d0N`. That way still works (docs/USER_GUIDE.md, "Translating
  WRF Namelists"), and the files are the same.
- **real.exe:** `cppwrf-realexe` writes real*8 files by default, as WRF's `-r8` build does. That is what gives WRF
  real*8's results. `--precision float` writes WRF's default single-precision files.

## Choosing the mode

| You want | Mode |
|---|---|
| The fastest run, for forecasts and ensembles | `gpu-fp32`. Its difference from WRF real*8 is the size of WRF real*4's. |
| WRF real*8's numbers, fast | `gpu-fp64`. It is bit-identical to WRF real*8 and about 2 to 3 times slower than FP32. |
| WRF real*8's numbers without a GPU | `cpu-fp64`. On one node it runs at about 0.6 to 0.9 times the speed of WRF real*8 on all its cores. |
| To check a GPU result | `cpu-fp64` on the same case, then `cppwrf-compare`. FP64 must be identical. |

docs/PERF.md has the measured speeds, and `cppwrf-bench` measures your machine.

## Comparing your WRF run with CPPWRF

Run CPPWRF in a copy of your WRF run directory, then compare the two sets of wrfout files:

```bash
cp -as $WRFRUN $CPPRUN          # links to the inputs and tables
rm $CPPRUN/wrfout_* $CPPRUN/wrfrst_* $CPPRUN/rsl.*
(cd $CPPRUN && cppwrf-wrf --mode cpu-fp64)
cppwrf-compare $WRFRUN $CPPRUN                       # every variable of every common frame, every domain
cppwrf-compare $WRFRUN/wrfout_d01_2026-10-12_00:00:00 $CPPRUN/wrfout_d01_2026-10-12_00:00:00 -v T,U,QVAPOR
```

For each variable, `cppwrf-compare` prints:

- the frames that differ;
- the largest absolute difference, the RMSE and the largest distance in units in the last place (ulp);
- the first differing point (frame, k, j, i).

`cppwrf-compare --exact` exits with status 1 at the first difference, and `--json` writes everything. It uses the
same comparison code as the verification suite (`tests/verify/exact.py`).

What to expect when you compare with your own WRF:

- **Your WRF is a real*8 build, run in one process** (`-np 1`, no OpenMP), built with gfortran on x86-64 without
  `-march`. Then `cpu-fp64` and `gpu-fp64` should be identical to it. If WRF was built with another gfortran than
  14.2, some vectorised loops may differ in the last bit, and those differences grow over a forecast.
- **Your WRF is the usual real*4 build, or runs with MPI.** Then expect differences of the size of WRF real*4
  against real*8: a relative L2 error of about 1e-4 to 1e-3 for U after an hour, more for W. Compare `gpu-fp32` with
  WRF real*4 by these norms, not by bits.
- **Some variables are left out of exact comparisons.** RQCBLTEN and RQIBLTEN in runs without cloud water or ice are
  uninitialised memory in WRF.

## The options CPPWRF ports

These tables are generated from the translator's own tables (`tools/gen_option_docs.py`; ctest
`docs_options_current` fails when they are stale). `cppwrf-check namelist.input` applies them to a namelist.

<!-- BEGIN generated by tools/gen_option_docs.py: do not edit by hand -->

### Physics and dynamics schemes

Each option accepts these values; any other value stops `cppwrf-wrf` (and the translator) at start-up.

| Option | Values CPPWRF runs |
|---|---|
| `mp_physics` | 0: None / off; 1: Kessler warm-rain; 3: WSM3 3-class simple ice (column-verified); 4: WSM5 5-class single-moment (column-verified); 6: WSM6 6-class single-moment; 8: Thompson; 10: Morrison two-moment (morr_rimed_ice 0 / 1, constant droplet number: progn = 0); 28: Thompson aerosol-aware (aer_init_opt = 0, wif_input_opt 0 / 1); 50: P3, one ice category, constant droplet number; 51: P3, one ice category, predicted droplet number |
| `sf_sfclay_physics` | 0: None / off; 1: Revised MM5 surface layer (sfclayrev); 2: Eta similarity surface layer (MYJ); 5: MYNN surface layer; 7: Pleim-Xiu surface layer (with the ACM2 PBL and Noah) |
| `sf_surface_physics` | 0: None / off; 1: 5-layer thermal diffusion land surface (slab); 2: Unified Noah land surface model (4 soil layers); 4: Noah-MP land surface model (v5.1.2, 4 soil layers) |
| `bl_pbl_physics` | 0: None / off; 1: YSU (Yonsei University) planetary boundary layer; 2: MYJ (Mellor-Yamada-Janjic) planetary boundary layer; 5: MYNN-EDMF planetary boundary layer; 7: ACM2 (Asymmetric Convective Model 2) planetary boundary layer; 8: BouLac (Bougeault-Lacarrere) TKE planetary boundary layer (column-verified, without BEP urban) |
| `cu_physics` | 0: None / off; 1: Kain-Fritsch cumulus parameterization; 2: Betts-Miller-Janjic cumulus parameterization; 3: Grell-Freitas cumulus parameterization; 6: Tiedtke cumulus parameterization; 16: New Tiedtke cumulus parameterization |
| `ra_lw_physics` | 0: None / off; 1: RRTM longwave radiation; 4: RRTMG longwave radiation |
| `ra_sw_physics` | 0: None / off; 1: Dudhia shortwave radiation; 4: RRTMG shortwave radiation |
| `diff_opt` | 0: None / off; 1: Cartesian or coordinate-surface diffusion; 2: Physical-space (full metric) diffusion |
| `km_opt` | 1: Constant eddy viscosity (khdif / kvdif); 2: 1.5-order prognostic TKE; 3: 3D Smagorinsky; 4: 2D horizontal Smagorinsky |
| `damp_opt` | 0: None / off; 1: Upper diffusive damping layer; 2: Rayleigh relaxation layer; 3: Implicit Rayleigh damping of w in upper layer |
| `h_mom_adv_order`, `h_sca_adv_order` | 5: 5th-order upwind; 6: 6th-order centered |
| `v_mom_adv_order`, `v_sca_adv_order` | 3: 3rd-order upwind; 4: 4th-order centered |

Options a nest may set differently from its parent: `bldt`, `bmj_rad_feedback`, `cu_physics`, `cu_rad_feedback`, `cudt`, `prec_acc_dt`, `radt`.

### physics_suite

| Suite | Options it sets (the `-1` entries) |
|---|---|
| `CONUS` | mp_physics = 8, cu_physics = 6, ra_lw_physics = 4, ra_sw_physics = 4, bl_pbl_physics = 2, sf_sfclay_physics = 2, sf_surface_physics = 2 |
| `TROPICAL` | mp_physics = 6, cu_physics = 16, ra_lw_physics = 4, ra_sw_physics = 4, bl_pbl_physics = 1, sf_sfclay_physics = 1, sf_surface_physics = 2 |

### Noah-MP options (`sf_surface_physics = 4`)

| Option | Values CPPWRF runs |
|---|---|
| `dveg` | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 |
| `opt_crs` | 1, 2 |
| `opt_btr` | 1, 2, 3 |
| `opt_runsrf` | 1, 2, 3, 4 |
| `opt_runsub` | 1, 2, 3, 4 |
| `opt_sfc` | 1, 2 |
| `opt_frz` | 1, 2 |
| `opt_inf` | 1, 2 |
| `opt_rad` | 1, 2, 3 |
| `opt_alb` | 1, 2 |
| `opt_snf` | 1, 2, 3, 4 |
| `opt_tbot` | 1, 2 |
| `opt_stc` | 1, 2, 3 |
| `opt_gla` | 1, 2 |
| `opt_rsf` | 1, 2, 3, 4 |
| `opt_soil` | 1 |
| `opt_pedo` | 1 |
| `opt_crop` | 0 |
| `opt_irr` | 0 |
| `opt_irrm` | 0 |
| `opt_infdv` | 0 |
| `opt_tdrn` | 0 |
| `opt_tksno` | 1, 2, 3, 4, 5 |
| `opt_compact` | 1, 2 |
| `opt_scf` | 1, 2, 3, 4, 5 |
| `opt_wetland` | 0, 1, 2 |

### Options that must keep WRF's default

CPPWRF runs these features only at WRF's default value (another value is an error).

| Option | WRF default | Feature without a port |
|---|---|---|
| `&dynamics.momentum_adv_opt` | `1` | WENO momentum advection |
| `&dynamics.phi_adv_z` | `1` | geopotential advection without double-staggered omega |
| `&physics.no_mp_heating` | `0` | microphysics without latent heating |
| `&physics.mp_zero_out` | `0` | microphysics zero-out of small moist values |
| `&bdy_control.polar` | `.false.` | polar filters and boundaries |
| `&bdy_control.have_bcs_moist` | `.false.` | moist species from the boundary file |
| `&bdy_control.have_bcs_scalar` | `.false.` | scalar species from the boundary file |
| `&dynamics.use_baseparam_fr_nml` | `.false.` | base-state parameters from the namelist |
| `&dynamics.use_input_w` | `.false.` | vertical velocity from the input file |
| `&domains.rebalance` | `0` | rebalancing of the input state |
| `&time_control.cycling` | `.false.` | cycling of the input state |
| `&dfi_control.dfi_opt` | `0` | digital filter initialisation |
| `&domains.allowed` | `.true.` | domains that are not integrated (grid_allowed) |
| `&dynamics.use_q_diabatic` | `0` | qv / qc diabatic terms in the scalar tendencies (q_diabatic_add) |
| `&dynamics.rad_nudge` | `0` | theta relaxation to the initial profile (theta_relaxation) |
| `&dynamics.do_avgflx_em` | `0` | time-averaged mass fluxes in the history output (AVGFLX_*) |
| `&dynamics.do_avgflx_cugd` | `0` | time-averaged convective mass fluxes in the history output (AVGFLX_*) |
| `&physics.acc_phy_tend` | `0` | accumulated physics tendencies in the history output (ATHCUTEN, ...) |
| `&physics.ua_phys` | `.false.` | UA Noah snow physics |
| `&physics.sf_surface_mosaic` | `0` | Noah mosaic land use |
| `&physics.sf_ocean_physics` | `0` | ocean mixed layer / 3-D ocean model |
| `&physics.omlcall` | `0` | ocean mixed layer |
| `&physics.aercu_opt` | `0` | aerosol-aware cumulus |
| `&physics.windfarm_opt` | `0` | wind farm parameterisation |
| `&physics.sf_surf_irr_scheme` | `0` | irrigation |
| `&physics.distributed_ahe_opt` | `0` | anthropogenic heat (AHE of the urban canopy model) |
| `&physics.kf_edrates` | `0` | KF entrainment / detrainment rate outputs (UDR_KF, DDR_KF, UER_KF, DER_KF, TIMEC_KF) |
| `&physics.ra_sw_eclipse` | `0` | solar eclipse in the shortwave radiation (ECOBSC, ECMASK) |
| `&diags.solar_diagnostics` | `0` | solar diagnostics (radiation_driver's Solar_step every step, SZA, CLRNIDX ...) |
| `&fdda.grid_sfdda` | `0` | surface grid nudging (needs wrfsfdda_d0N from OBSGRID) |
| `&fdda.obs_nudge_opt` | `0` | observation nudging |
| `&fire.ifire` | `0` | wildland fire model |
| `&stoch.multi_perturb` | `0` | stochastic multi-physics perturbations |
| `&stoch.sppt` | `0` | stochastically perturbed physics tendencies |
| `&stoch.skebs` | `0` | stochastic kinetic-energy backscatter |
| `&stoch.rand_perturb` | `0` | random perturbations |
| `&stoch.spp` | `0` | stochastically perturbed parameters |
| `&stoch.spp_pbl` | `0` | stochastically perturbed PBL parameters |
| `&stoch.spp_conv` | `0` | stochastically perturbed convection parameters |
| `&stoch.spp_lsm` | `0` | stochastically perturbed land-surface parameters |
| `&stoch.perturb_bdy` | `0` | perturbed lateral boundary tendencies (spec_bdy_dry_perturb) |
| `&stoch.perturb_chem_bdy` | `0` | perturbed chemistry boundaries |
| `&bdy_control.constant_bc` | `.false.` | constant lateral boundaries (zero_bdytend) |
| `&bdy_control.multi_bdy_files` | `.false.` | one boundary file per interval |

### Switches without a port (must be off)

| Option | Feature |
|---|---|
| `shcu_physics` | Separate shallow cumulus scheme |
| `sf_lake_physics` | Lake physics model |
| `topo_wind` | Topographic wind correction |
| `slope_rad` | Slope-dependent radiation |
| `topo_shading` | Topographic shading |
| `swint_opt` | Shortwave radiation interpolation between radt calls |
| `seaice_thickness_opt` | Variable sea-ice thickness option |
| `seaice_snowdepth_opt` | Sea-ice snow depth option |
| `seaice_albedo_opt` | Sea-ice albedo option |
| `ra_call_offset` | Radiation call offset |
| `shalwater_z0` | Shallow-water roughness formulation |
| `scm_force_flux` | Single-column model prescribed surface flux |
| `diff_6th_slopeopt` | 6th-order diffusion slope limiter |
| `sfs_opt` | nonlinear backscatter / anisotropy (NBA) subfilter stresses |
| `m_opt` | subfilter stress output |
| `tracer_mix2_off` | tracer 2nd-order diffusion switch |
| `tracer_mix6_off` | tracer 6th-order diffusion switch |
| `chem_mix2_off` | chemistry 2nd-order diffusion switch |
| `chem_mix6_off` | chemistry 6th-order diffusion switch |
| `zadvect_implicit` | implicit-explicit vertical advection (IEVA) |
| `pert_coriolis` | perturbation Coriolis |
| `tracer_opt` | passive tracers |

### Options that change nothing in a CPPWRF run

WRF reads them but they do not change the run here (WPS / real.exe options, infrastructure, or options of features that are off); `cppwrf-check` lists them as ignored or ported with the reason.

| Option | Why |
|---|---|
| `adaptation_domain` | adaptive time step only |
| `adjust_heights` | real.exe only |
| `aggregate_lu` | real.exe only |
| `auto_levels_opt` | real.exe only |
| `auxinput4_inname` | wrflowinp file name; used with sst_update = 1 |
| `auxinput4_interval` | wrflowinp interval; used with sst_update = 1 |
| `base_temp` | Base potential temperature (read from sounding or wrfinput) |
| `coriolis2d` | not used by wrf.exe |
| `dampcoef` | damp_opt = 0 |
| `debug_level` | WRF internal logging verbosity |
| `diff_6th_factor` | Sixth-order diffusion factor (default 0.12) |
| `dk_zfac_ph` | spectral nudging (grid_fdda = 2) only |
| `dk_zfac_q` | spectral nudging (grid_fdda = 2) only |
| `dk_zfac_t` | spectral nudging (grid_fdda = 2) only |
| `dk_zfac_uv` | spectral nudging (grid_fdda = 2) only |
| `do_coriolis` | not used by wrf.exe |
| `do_curvature` | not used by wrf.exe |
| `do_gradp` | not used by wrf.exe |
| `dtramp_min` | grid_fdda = 1 only |
| `dx` | Grid spacing (real case: the namelist value, checked against wrfinput as input_wrf does) |
| `dy` | Grid spacing (real case: the namelist value, checked against wrfinput as input_wrf does) |
| `dzbot` | real.exe vertical level generation (eta levels read from wrfinput_d01) |
| `dzstretch_s` | real.exe vertical level generation (eta levels read from wrfinput_d01) |
| `dzstretch_u` | real.exe vertical level generation (eta levels read from wrfinput_d01) |
| `e_sn` | Grid dimension (real case: checked against wrfinput as input_wrf does) |
| `e_vert` | Grid dimension (real case: checked against wrfinput as input_wrf does) |
| `e_we` | Grid dimension (real case: checked against wrfinput as input_wrf does) |
| `emdiv` | External mode filter coefficient (default 0.01) |
| `end_day` | End date component |
| `end_hour` | End date component |
| `end_minute` | End date component |
| `end_month` | End date component |
| `end_second` | End date component |
| `end_year` | End date component |
| `epssm` | Acoustic divergence damping parameter (default 0.1) |
| `eta_levels` | real.exe only |
| `extrap_type` | real.exe only |
| `feedback` | Two-way nesting feedback |
| `fgdt` | grid_fdda = 1 only |
| `fgdtzero` | grid_fdda = 1 only |
| `force_sfc_in_vinterp` | real.exe only |
| `frames_per_outfile` | Frame count per NetCDF file handled by output format |
| `gfdda_begin_d` | grid_fdda = 1 only |
| `gfdda_begin_h` | grid_fdda = 1 only |
| `gfdda_begin_m` | grid_fdda = 1 only |
| `gfdda_begin_s` | grid_fdda = 1 only |
| `gfdda_begin_y` | grid_fdda = 1 only |
| `gfdda_end_d` | grid_fdda = 1 only |
| `gfdda_end_h` | grid_fdda = 1 only |
| `gfdda_end_m` | grid_fdda = 1 only |
| `gfdda_end_s` | grid_fdda = 1 only |
| `gfdda_end_y` | grid_fdda = 1 only |
| `gfdda_inname` | grid_fdda = 1 only |
| `gfdda_interval` | grid_fdda = 1 only |
| `gfdda_interval_d` | grid_fdda = 1 only |
| `gfdda_interval_h` | grid_fdda = 1 only |
| `gfdda_interval_m` | grid_fdda = 1 only |
| `gfdda_interval_s` | grid_fdda = 1 only |
| `gfdda_interval_y` | grid_fdda = 1 only |
| `gph` | spectral nudging (grid_fdda = 2) only |
| `gq` | grid_fdda = 1 only |
| `gq_sfc` | surface analysis nudging (grid_sfdda) only |
| `grid_id` | Domain grid hierarchy identifier |
| `gt` | grid_fdda = 1 only |
| `gt_sfc` | surface analysis nudging (grid_sfdda) only |
| `guv` | grid_fdda = 1 only |
| `guv_sfc` | surface analysis nudging (grid_sfdda) only |
| `h_sca_adv_order` | Scalar horizontal advection order (coupled with h_mom_adv_order in CPPWRF) |
| `hybrid_opt` | Real case hybrid coordinate flag read directly from wrfinput_d01 global attribute |
| `i_parent_start` | Domain nesting starting index |
| `if_no_pbl_nudging_ph` | spectral nudging (grid_fdda = 2) only |
| `if_no_pbl_nudging_q` | grid_fdda = 1 only |
| `if_no_pbl_nudging_t` | grid_fdda = 1 only |
| `if_no_pbl_nudging_uv` | grid_fdda = 1 only |
| `if_ramping` | grid_fdda = 1 only |
| `if_zfac_ph` | spectral nudging (grid_fdda = 2) only |
| `if_zfac_q` | grid_fdda = 1 only |
| `if_zfac_t` | grid_fdda = 1 only |
| `if_zfac_uv` | grid_fdda = 1 only |
| `ifsnow` | Snow effects flag (used by the slab LSM only) |
| `ignore_iofields_warning` | read by the history output from namelist.input |
| `input_from_file` | Real case initialization always reads wrfinput |
| `interp_theta` | real.exe only |
| `interp_type` | real.exe only |
| `interval_seconds` | Lateral boundary interval (read directly from wrfbdy_d01 NetCDF in real case) |
| `io_form_auxinput4` | wrflowinp format (NetCDF); used with sst_update = 1 |
| `io_form_gfdda` | grid_fdda = 1 only |
| `io_form_restart` | No restart files (restart_interval = 0) |
| `io_form_sgfdda` | surface analysis nudging (grid_sfdda) only |
| `iofields_filename` | Custom NetCDF field selection (handled by WrfHistory via output.wrf_namelist) |
| `j_parent_start` | Domain nesting starting index |
| `k_zfac_ph` | spectral nudging (grid_fdda = 2) only |
| `k_zfac_q` | grid_fdda = 1 only |
| `k_zfac_t` | grid_fdda = 1 only |
| `k_zfac_uv` | grid_fdda = 1 only |
| `ktrop` | spectral nudging (grid_fdda = 2) only |
| `lagrange_order` | real.exe only |
| `linear_interp` | real.exe only |
| `lowest_lev_from_sfc` | real.exe only |
| `max_dz` | real.exe only |
| `max_step_increase_pct` | adaptive time step only |
| `max_time_step` | adaptive time step only |
| `max_time_step_den` | adaptive time step only |
| `max_ts_level` | time series (tslist) only |
| `max_ts_locs` | time series (tslist) only |
| `maxw_above_this_level` | real.exe only |
| `maxw_horiz_pres_diff` | real.exe only |
| `min_time_step` | adaptive time step only |
| `min_time_step_den` | adaptive time step only |
| `nested` | Child domain nested lateral boundary flag (handled automatically for [d02]..[d0N]) |
| `nio_groups` | WRF quilt I/O server |
| `nio_tasks_per_group` | WRF quilt I/O server |
| `nssl_2moment_on` | NSSL microphysics only (check_a_mundo resets it) |
| `nssl_3moment` | NSSL microphysics only (check_a_mundo resets it) |
| `nssl_ccn_on` | NSSL microphysics only (check_a_mundo resets it) |
| `nssl_density_on` | NSSL microphysics only (check_a_mundo resets it) |
| `nssl_hail_on` | NSSL microphysics only (check_a_mundo resets it) |
| `num_land_cat` | Land-use category count (real case: checked against wrfinput as input_wrf does) |
| `num_metgrid_levels` | Metgrid vertical levels used by real.exe |
| `num_metgrid_soil_levels` | Metgrid soil levels used by real.exe |
| `num_soil_layers` | 4 soil layers (Noah, Noah-MP), 5 (slab) |
| `override_restart_timers` | restart runs only |
| `p_top_requested` | Real case p_top read directly from wrfinput_d01 |
| `parent_grid_ratio` | Nesting grid refinement ratio |
| `parent_id` | Domain nesting parent identifier |
| `parent_time_step_ratio` | Nesting timestep ratio |
| `qv_max_flag` | real.exe only |
| `qv_max_p_safe` | real.exe only |
| `qv_max_value` | real.exe only |
| `qv_min_flag` | real.exe only |
| `qv_min_p_safe` | real.exe only |
| `qv_min_value` | real.exe only |
| `rh2qv_method` | real.exe only |
| `rh2qv_wrt_liquid` | real.exe only |
| `rinblw` | surface analysis nudging (grid_sfdda) only |
| `s_sn` | Domain starting index in Y (always 1) |
| `s_vert` | Domain starting index in Z (always 1) |
| `s_we` | Domain starting index in X (always 1) |
| `sfcp_to_sfcp` | real.exe only |
| `sgfdda_begin_d` | surface analysis nudging (grid_sfdda) only |
| `sgfdda_begin_h` | surface analysis nudging (grid_sfdda) only |
| `sgfdda_begin_m` | surface analysis nudging (grid_sfdda) only |
| `sgfdda_begin_s` | surface analysis nudging (grid_sfdda) only |
| `sgfdda_begin_y` | surface analysis nudging (grid_sfdda) only |
| `sgfdda_end_d` | surface analysis nudging (grid_sfdda) only |
| `sgfdda_end_h` | surface analysis nudging (grid_sfdda) only |
| `sgfdda_end_m` | surface analysis nudging (grid_sfdda) only |
| `sgfdda_end_s` | surface analysis nudging (grid_sfdda) only |
| `sgfdda_end_y` | surface analysis nudging (grid_sfdda) only |
| `sgfdda_inname` | surface analysis nudging (grid_sfdda) only |
| `sgfdda_interval` | surface analysis nudging (grid_sfdda) only |
| `sgfdda_interval_d` | surface analysis nudging (grid_sfdda) only |
| `sgfdda_interval_h` | surface analysis nudging (grid_sfdda) only |
| `sgfdda_interval_m` | surface analysis nudging (grid_sfdda) only |
| `sgfdda_interval_s` | surface analysis nudging (grid_sfdda) only |
| `sgfdda_interval_y` | surface analysis nudging (grid_sfdda) only |
| `shallowcu_forced_ra` | KF-CuP (cu_physics = 10) only |
| `smdiv` | Acoustic divergence damping parameter (default 0.1) |
| `smooth_cg_topo` | real.exe only |
| `smooth_option` | Nesting boundary smoothing |
| `snicar_aerosol_readtable` | SNICAR (opt_alb = 3) only |
| `snicar_bandnumber_opt` | SNICAR (opt_alb = 3) only |
| `snicar_dustoptics_opt` | SNICAR (opt_alb = 3) only |
| `snicar_rtsolver_opt` | SNICAR (opt_alb = 3) only |
| `snicar_snowbc_intmix` | SNICAR (opt_alb = 3) only |
| `snicar_snowdust_intmix` | SNICAR (opt_alb = 3) only |
| `snicar_snowoptics_opt` | SNICAR (opt_alb = 3) only |
| `snicar_snowshape_opt` | SNICAR (opt_alb = 3) only |
| `snicar_solarspec_opt` | SNICAR (opt_alb = 3) only |
| `snicar_use_aerosol` | SNICAR (opt_alb = 3) only |
| `snicar_use_oc` | SNICAR (opt_alb = 3) only |
| `spec_bdy_width` | Real case lateral boundary width (must be 5; read from wrfbdy_d01) |
| `start_day` | Start date component |
| `start_hour` | Start date component |
| `start_minute` | Start date component |
| `start_month` | Start date component |
| `start_second` | Start date component |
| `start_year` | Start date component |
| `starting_time_step` | adaptive time step only |
| `starting_time_step_den` | adaptive time step only |
| `step_to_output_time` | adaptive time step only |
| `swrad_scat` | read by the Dudhia shortwave scheme only |
| `t_extrap_type` | real.exe only |
| `target_cfl` | adaptive time step only |
| `target_hcfl` | adaptive time step only |
| `tc` | tc.exe (bogus vortex) only |
| `time_step_fract_den` | Combined into time.dt if non-zero |
| `time_step_fract_num` | Combined into time.dt if non-zero |
| `top_radiation` | not used by wrf.exe |
| `tracer_pblmix` | tracer_opt = 0 |
| `trop_horiz_pres_diff` | real.exe only |
| `ts_buf_size` | time series (tslist) only |
| `tslist_unstagger_winds` | time series (tslist) only |
| `use_levels_below_ground` | real.exe only |
| `use_maxw_level` | real.exe only |
| `use_netcdf_classic` | NetCDF classic format of the wrfout / wrfrst files (read by the model from the WRF namelist) |
| `use_sh_qv` | real.exe only |
| `use_surface` | real.exe only |
| `use_tavg_for_tsk` | real.exe only |
| `use_trop_level` | real.exe only |
| `v_sca_adv_order` | Scalar vertical advection order (coupled with v_mom_adv_order in CPPWRF) |
| `wif_input_opt` | aerosol-aware microphysics (mp_physics = 28) only |
| `write_hist_at_0h_rst` | restart runs only |
| `xwavenum` | spectral nudging (grid_fdda = 2) only |
| `ywavenum` | spectral nudging (grid_fdda = 2) only |
| `zap_close_levels` | real.exe only |
| `zdamp` | damp_opt = 0 |

Every other option of WRF 4.8.0's namelist (2231 options in all) is an error when it is away from its Registry default and CPPWRF has no port of it; an option that is not in WRF 4.8.0 is an error as in wrf.exe.

<!-- END generated -->
