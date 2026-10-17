# Physics driver and dycore coupling

This page describes where the column physics is called inside a CPPWRF time step, which diagnostic
fields it receives, and how its tendencies return to the dynamical core. Every step follows the
corresponding WRF routine; the WRF file is named in each section. The C++ code is in
`src/phys/driver/` (suite, radiation, phy_prep, tendency coupling, initialisation) and
`src/dyn/rk3.cpp`, `src/dyn/physics_coupling.cpp` (call sites in the time step).

## 1. Order of calls in one time step

WRF `solve_em` (dyn_em/solve_em.F) and CPPWRF `RK3::step` (src/dyn/rk3.cpp) run the same sequence:

| # | WRF routine | CPPWRF | When |
|---|---|---|---|
| 1 | `first_rk_step_part1`: `init_zero_tendency` | `LargeStep::zero_tendf` | RK stage 1 |
| 2 | `phy_prep` | `run_phy_prep` | RK stage 1 |
| 3 | `radiation_driver` (RRTMG LW/SW, every `stepra` steps) | `RadiationDriver::run` | RK stage 1 |
| 4 | `surface_driver`: surface layer, then the land surface model | `PhysicsSuite::surface_driver` | RK stage 1 |
| 5 | `pbl_driver` (every `stepbl` steps) | `PhysicsSuite::pbl_driver` | RK stage 1 |
| 6 | `cumulus_driver` (every `stepcu` steps) | `PhysicsSuite::cumulus_driver` | RK stage 1 |
| 7 | `first_rk_step_part2`: diffusion metrics and eddy viscosities, `calculate_phy_tend`, `update_phy_ten`, `conv_t_tendf_to_moist` (use_theta_m = 1) | `DiffusionDriver::prepare`, `PhysicsSuite::part2` | RK stage 1 |
| 8 | `rk_tendency`, acoustic steps, scalar transport | `LargeStep`, `SmallStep`, `ScalarTransport` | RK stages 1-3 |
| 9 | `advance_ppt`, `phy_prep_part2` (advective forcing for New Tiedtke / Grell-Freitas) | `PhysicsSuite::after_rk` | after RK3 |
| 10 | `moist_physics_prep_em`, `microphysics_driver`, `moist_physics_finish_em` | `PhysicsCoupling::step` | after RK3 |
| 11 | `calc_p_rho_phi` (tile without the specified zone), `spec_bdy_final`, `set_w_surface` | `SmallStep`, `LateralBC::final` | end of step |

The physics in rows 2-7 therefore sees the state at the start of the step. Microphysics (row 10)
adjusts the state after the RK3 integration and does not produce tendencies.

## 2. Diagnostics from `phy_prep`

`phy_prep` (dyn_em/module_big_step_utilities_em.F) is ported in `src/phys/driver/phy_prep.cpp`
(`KPhyPrep`). Fields at mass points unless noted:

| Field | Definition |
|---|---|
| `th_phy` | dry potential temperature: `(t + t0) / (1 + rvovrd*qv)` with use_theta_m = 1, else `t + t0` |
| `p_phy` | `p + pb` |
| `pi_phy` | `(p_phy / p1000mb)**rcp`, evaluated with glibc libmvec lanes where gfortran vectorises the loop |
| `t_phy` | `th_phy * pi_phy` |
| `rho` | `1 / alt * (1 + qv)` |
| `u_phy`, `v_phy` | averages of the two neighbouring C-grid velocities |
| `z_at_w`, `dz8w`, `z` | heights of the full levels, layer thicknesses, mass-level heights from `ph + phb` |
| `p8w`, `t8w` (W levels) | interface pressure and temperature: `fnm`/`fnp` weights inside, extrapolation with `cf1..cf3` at the surface and `cfn`, `cfn1` at the top |
| `p_hyd_w`, `p_hyd` | hydrostatic pressure at full and half levels, integrated downward from `p_top` with the hybrid coefficients and the total moisture |

## 3. Inputs and outputs of each scheme

Scheme selection follows the WRF namelist options (`sf_sfclay_physics`, `sf_surface_physics`,
`bl_pbl_physics`, `cu_physics`, `ra_lw_physics`, `ra_sw_physics`, `mp_physics`).

### 3.1 Surface layer (module_surface_driver.F)

Before the surface layer, `psfc = p_hyd_w(k=1)` and `u_phytmp`, `v_phytmp` (= `u_phy`, `v_phy`, no
moving frame) are set.

- **Revised MM5 (`sf_sfclay_physics = 1`, src/phys/sf_sfclayrev):** lowest level of `u_phytmp`,
  `v_phytmp`, `t_phy`, `qv`, `p_hyd`, `dz8w`, plus `psfc`, `tsk`, `xland`, `mavail`, `znt`, `pblh`.
  Outputs `ust`, `znt`, `mol`, `rmol`, `zol`, `regime`, `psim`, `psih`, `fm`, `fh`, `br`, `wspd`,
  `chs`, `chs2`, `cqs2`, `flhc`, `flqc`, `qgh`, `qsfc`, `hfx`, `qfx`, `lh`, `u10`, `v10`, `t2`,
  `th2`, `q2`; then `ch = chs`.
- **Eta / MYJ (`sf_sfclay_physics = 2`, src/phys/sf_myjsfc):** full columns of `p_hyd`, `p_hyd_w`,
  `th_phy`, `t_phy`, `qv`, `qc`, `u_phy`, `v_phy` and the MYJ `tke_pbl`; outputs as above plus the
  roughness-level values `thz0`, `qz0`, `uz0`, `vz0`, `akhs`, `akms` and the shelter diagnostics.

### 3.2 Land surface: Noah (`sf_surface_physics = 2`, src/phys/lsm_noah)

The Noah driver (phys/module_sf_noahdrv.F) computes the forcing per point (`sfcprs` from `p8w`,
`sfctmp` from `t_phy`, `q2k`, `zlvl` from `dz8w`, downward radiation `glw`/`swdown` with `albedo`,
precipitation `rainbl` with the frozen fraction `sr` from the microphysics of the previous step),
calls SFLX (SFLX_GLACIAL for land ice, SFLX_SEAICE for sea ice) and writes back `tsk`, `hfx`, `qfx`,
`lh`, `grdflx`, `tslb`, `smois`, `sh2o`, `snow`, `snowh`, `snowc`, `canwat`, `albedo`, `emiss` and
the accumulated budget fields.

### 3.3 Planetary boundary layer (module_pbl_driver.F)

- **YSU (`bl_pbl_physics = 1`, src/phys/bl_ysu):** `u_phytmp`, `v_phytmp`, `t_phy`, `qv`, `qc`,
  `qi`, `p_hyd` (P3D), `p_hyd_w` (P3DI), `pi_phy`, `dz8w`, `rthraten`, and the surface fields `psfc`,
  `znt`, `ust`, `fm`, `fh`, `xland`, `hfx`, `qfx`, `wspd`, `br`. Outputs `rublten`, `rvblten`,
  `rthblten`, `rqvblten`, `rqcblten`, `rqiblten`, `exch_h`, `exch_m`, `pblh`, `kpbl`, `wstar`,
  `delta`.
- **MYJ (`bl_pbl_physics = 2`, src/phys/bl_myj):** prognostic `tke_pbl` and `el_pbl`; same
  tendency outputs, plus `exch_h`, `mixht`, `pblh`, `kpbl`.

### 3.4 Cumulus (module_cumulus_driver.F)

Kain-Fritsch (`cu_physics = 1`), Grell-Freitas (`3`, deep and shallow) and New Tiedtke (`16`) read
`t_phy`, `qv`, the condensates, `p_hyd`, `p_hyd_w`, `pi_phy`, `dz8w`, `rho`, `u_phy`, `v_phy` and
`w`; New Tiedtke and Grell-Freitas also use the advective forcing `rthften`, `rqvften` from
`phy_prep_part2` of the previous step. Outputs: `rthcuten`, `rqvcuten`, `rqccuten`, `rqicuten`
(KF also `rqrcuten`, `rqscuten`; New Tiedtke also the momentum tendencies `rucuten`, `rvcuten`),
`raincv`, `pratec`, `nca` (KF), `cutop`, `cubot`.

### 3.5 Radiation (module_radiation_driver.F)

RRTMG LW/SW (`ra_lw_physics = ra_sw_physics = 4`) with CAM ozone and greenhouse gases, `cal_cldfra1`
and, with WSM6, the effective radii from `mp_wsm6_effectRad`. Inputs `t_phy`, `p_hyd`, `p_hyd_w`,
`pi_phy`, `rho`, `dz8w`, `t8w`, moisture, `tsk`, `albedo`, `emiss`, `xland`, `xice`, `snow`.
Outputs `rthraten`, `glw`, `gsw`, `swdown`, `coszen`, `cldfra`. Between radiation calls the stored
`rthraten` is reused, as in WRF.

### 3.6 Microphysics (module_microphysics_driver.F)

Kessler (`mp_physics = 1`), WSM6 (`6`) and Thompson (`8`) run after the RK3 integration on the
updated state. `moist_physics_prep_em` builds `th_phy`, `pii`, `rho`, `dz8w`, `z`; the scheme
updates theta and the moisture species directly and produces `rainnc`, `rainncv`, `snownc`,
`graupelnc`, `sr`; `moist_physics_finish_em` converts back to `theta_m - t0` and stores the heating
rate `h_diabatic`.

## 4. Return of the tendencies to the dynamics

`module_physics_addtendc.F` is ported in `src/phys/driver/tendency_coupling.cpp`:

1. `calculate_phy_tend` couples each tendency with the column mass of the hybrid coordinate,
   `ten = (c1(k)*mut + c2(k)) * ten` (mass points; U and V points use the averaged column mass).
2. `update_phy_ten` adds the coupled radiation, PBL and cumulus tendencies to `t_tendf`, `ru_tendf`,
   `rv_tendf` and `q_tendf` (`add_a2a`, and `add_a2c_u` / `add_a2c_v` for momentum, with WRF's loop
   bounds for specified and periodic boundaries).
3. With `use_theta_m = 1`, `conv_t_tendf_to_moist` converts the dry theta tendency to moist theta
   using the moisture tendency, with WRF's one-step lag of `t_1`.

`rk_tendency` then combines these physics tendencies with the dynamical tendencies in every RK stage.

## 5. Initialisation for real data

`src/phys/driver/physics_init.cpp` and `src/cases/real_physics.cpp` follow `phy_init`
(phys/module_physics_init.F): `landuse_init` with `LANDUSE.TBL` (season from the date and latitude,
sea-ice overrides), the Noah tables (`SOILPARM.TBL`, `VEGPARM.TBL`, `GENPARM.TBL`) and `lsminit`
(liquid soil water from `FRH2O`), the surface layer tables (`sf_sfclayrev_init`, `myjsfcinit`),
`myjpblinit`, the KF lookup tables (`kf_lutab`), the WSM6 and Thompson tables, and the RRTMG
absorption data (`data/rrtmg`, WRF 4.6.1 images).
