# P1 physics suite: WRF-identical driver glue (sfclayrev + Noah + YSU)

`src/phys/driver/suite.{hpp,cpp}` (class `PhysicsSuite`, a `dyn::ColumnPhysics`) is a
statement-by-statement port of the path WRF v4.6.1 executes for

```
sf_sfclay_physics = 1 (sfclayrev)   sf_surface_physics = 2 (Noah, num_soil_layers = 4)
bl_pbl_physics = 1 (YSU)            bldt = 0, isfflx = 1, fractional_seaice = 0
mp_physics = cu_physics = ra_lw_physics = ra_sw_physics = 0
real data (specified boundaries), use_theta_m = 1, sf_urban_physics = 0, no lake, sst_update = 0
```

Differential test against the **original** WRF code: `tests/phys/suite/` (below). Result: all
94 compared fields, 3 consecutive physics calls with the state carried over, **bit-identical**
in CPU FP64 (table at the end).

## 1. What WRF executes for these options

Per time step (`solve_em`): `grid%itimestep = grid%itimestep + 1` (solve_em.F:371), then in RK
stage 1:

**first_rk_step_part1** (module_first_rk_step_part1.F)
1. `init_zero_tendency` (:181) - done by the C++ dycore before `ColumnPhysics::part1`.
2. `phy_prep` (:207, module_big_step_utilities_em.F:4730): th_phy (moist theta -> dry theta with
   qv), p_phy = p+pb, pi_phy, t_phy, rho = 1/alt*(1+qv), u_phy, v_phy, z_at_w, dz8w, z, p8w, t8w
   (with the bottom/top extrapolations), p_hyd_w (hydrostatic, from mut and all moist
   species), p_hyd.
3. radiation_driver: nothing for ra = 0 (GLW, SWDOWN, GSW, RTHRATEN keep their values).
4. `surface_driver` (:593). Arguments to note: **`P_PHY=grid%p_hyd`, `P8W=grid%p_hyd_w`**,
   `PI_PHY=pi_phy`, `T_PHY=t_phy`, `SR=grid%sr` present (so `FRPCPN = .true.`, FFROZP = SR).
   Inside (module_surface_driver.F):
   - every call: QGH, CHS, CPM, CHS2 = 0 and u/v_phytmp = 0 (locals of surface_driver), RAINBL +=
     RAINCV + RAINNCV + RAINSHV, max(.,0) (:1539-1570);
   - run decision (:1855-1895): itimestep == 1, or bldt == 0 / stepbl == 1, or mod(itimestep,
     stepbl) == 0; DTBL = dt (bldt = 0) else dt*stepbl;
   - PSFC = p8w(kts) (= p_hyd_w), u_phytmp = u_phy + u_frame minus surface currents (:1982);
   - `SFCLAYREV` wrapper (:2116, module_sf_sfclayrev.F) on level kts; then `ch = chs`;
   - `lsm` (module_sf_noahdrv.F:38): itimestep == 1 water/sea-ice soil init; per point Q2SAT,
     DQSDT2 (with the snow branch), TH2 (Exner), LWDN, SOLNET = SWDOWN*(1-ALBEDO), PRCP =
     RAINBL/DTBL, soil type 14 -> 7, SNOWH fix, urban -> ISURBAN; sea ice (XICE >= 0.5) is
     skipped (SH2O = 1, LAI = 0.01); land ice (IVGTYP == ISICE) -> `SFLX_GLACIAL`; other land ->
     `SFLX`; post-processing (CANWAT, SNOW, ALBBCK = ALB, Z0 = Z0BRD, ZNT = Z0K, TSK, HFX, LH, QFX,
     GRDFLX, POTEVP, CHS2 = CQS2, QSFC = Q1/(1-Q1), NOAHRES, SMSTAV, SMSTOT, SMCREL, runoffs,
     ACSNOW, ACSNOM, SNOPCX);
   - `seaice_noah` (module_sf_noah_seaice_drv.F:14, `SFLX_SEAICE`, DT = model dt);
   - SFCEVP, SFCEXC = CHS, ACHFX, ACLHF, ACGRDFLX; `SFCDIAGS` (:2994; T2, TH2, Q2 over all points);
   - RAINBL = 0; Q2 = min(Q2, 1.05*qv(kts)) over land (:4436-4453).
5. `pbl_driver` (:1112, module_pbl_driver.F): same run decision; u/v_phytmp = u/v_phy + frame,
   PSFC = p8w(kms), RTHBLTEN/RUBLTEN/RVBLTEN/RQVBLTEN/RQCBLTEN (= RQIBLTEN if qi) = 0; `ysu`
   with `P3D=p_hyd, P3DI=p_hyd_w, PI3D=pi_phy, T3D=t_phy, PSIM=grid%fm, PSIH=grid%fh, DT=dtbl,
   HPBL=pblh, WSTAR=grid%wstar_ysu, DELTA=grid%delta_ysu`; the module_bl_ysu.F wrapper divides
   the temperature tendency of bl_ysu_run by pi again (`NEED_B4B_DURING_CCPP_TESTING`).
6. cumulus_driver: nothing for cu = 0.

**first_rk_step_part2** (module_first_rk_step_part2.F)
7. `calculate_phy_tend` (:392, module_em.F:2079): R*BLTEN *= (c1h*mut + c2h).
8. `update_phy_ten` (:793, module_physics_addtendc.F): `phy_bl_ten` YSU: add_a2a(t_tendf,
   RTHBLTEN), add_a2c_u(ru_tendf, RUBLTEN), add_a2c_v(rv_tendf, RVBLTEN), add_a2a(moist_tend qv),
   (qc, qi if present); specified boundaries: the outer row is excluded.
9. `conv_t_tendf_to_moist` (:997): theta tendency -> moist-theta tendency with t_1, moist_old
   (qv) and the qv tendency. (This belongs to the physics coupling: it was missing in CPPWRF.)

**after the RK loop** (solve_em.F:3580)
10. `phy_prep_part2`: R*BLTEN /= (c1h*muts + c2h).

**State carried between steps** (WRF grid fields): TSK, TSLB, SMOIS, SH2O, SMCREL, SNOW,
SNOWH, SNOWC, CANWAT, ALBEDO, ALBBCK, EMISS, Z0, ZNT, LAI, SNOTIME, UST, MOL, RMOL, PBLH, KPBL,
HFX, QFX, LH, QSFC, FLHC, FLQC, BR, WSPD, FM, FH, ZOL, REGIME, U10, V10, T2, TH2, Q2, PSFC,
GRDFLX, RAINBL, the accumulations (ACSNOM, ACSNOW, SNOPCX, POTEVP, SFCEVP, SFCRUNOFF, UDRUNOFF,
ACHFX, ACLHF, ACGRDFLX), R*BLTEN (decoupled), EXCH_H/M. QGH, CHS, CHS2, CQS2, CPM are locals of
surface_driver and psim, psih, gz1oz0, chklowq, th_phy, pi_phy, p_phy, dz8w, p8w, t8w locals
of solve_em: they are recomputed each step (scratch members in C++).

## 2. Call map (WRF file:line -> C++ file:line)

| WRF v4.6.1 | C++ |
|---|---|
| solve_em.F:371 `itimestep + 1` | src/phys/driver/suite.cpp:605 (`PhysicsSuite::part1`) |
| module_big_step_utilities_em.F:4730 `phy_prep` | src/phys/driver/phy_prep.cpp:16 (`KPhyPrep`), :88 |
| module_surface_driver.F:1539-1570 zeroed locals, RAINBL | suite.cpp:46 (`KSfcAlways`) |
| module_surface_driver.F:1855-1895, 1955-1975 run decision, DTBL | suite.cpp:630 |
| module_surface_driver.F:1980-1997 PSFC, u/v_phytmp | suite.cpp:68 (`KSfcPrep`) |
| module_surface_driver.F:2116 `SFCLAYREV` (module_sf_sfclayrev.F:16) | suite.cpp:641 -> src/phys/sf_sfclayrev |
| module_surface_driver.F:2133 `ch = chs` | suite.cpp:652 |
| module_sf_noahdrv.F:747-1100 lsm before SFLX | suite.cpp:101 (`KLsmPre`) |
| module_sf_noahdrv.F:1105 `SFLX` | suite.cpp:716 -> src/phys/lsm_noah/noah.hpp (`NoahKernel`) |
| module_sf_noahdrv.F:1140-1170 ICE == -1, `SFLX_GLACIAL` | suite.cpp:255 (`KLsmGlacial`), src/phys/driver/noah_ice.hpp:284 |
| module_sf_noahdrv.F:1215-1600 lsm after SFLX | suite.cpp:309 (`KLsmPost`) |
| module_sf_noah_seaice_drv.F:14 `seaice_noah`, module_sf_noah_seaice.F `SFLX_SEAICE` | suite.cpp:366 (`KSeaiceNoah`), noah_ice.hpp:627 |
| module_surface_driver.F:2984-2994 accumulations, `SFCDIAGS` (module_sf_sfcdiags.F:7) | suite.cpp:461 (`KLsmDiag`) |
| module_surface_driver.F:4436-4453 RAINBL = 0, Q2 limit | suite.cpp:491 (`KSfcEnd`) |
| module_pbl_driver.F:930-1045 run decision, phytmp, PSFC, zeroed R*BLTEN | suite.cpp:791, :505 (`KPblPrep`) |
| module_pbl_driver.F:1217 `ysu` (module_bl_ysu.F:398 `bl_ysu_run`) | suite.cpp:802 -> src/phys/bl_ysu/ysu.hpp |
| module_em.F:2079 `calculate_phy_tend` | src/phys/driver/tendency_coupling.cpp:81 |
| module_physics_addtendc.F:28/210 `update_phy_ten`/`phy_bl_ten` (add_a2a :2285, add_a2c_u :2381, add_a2c_v :2435) | tendency_coupling.cpp:96 |
| module_big_step_utilities_em.F:6674 `conv_t_tendf_to_moist` | tendency_coupling.cpp:111 |
| module_big_step_utilities_em.F:4980 `phy_prep_part2` | tendency_coupling.cpp:116 (from suite.cpp:618 `after_rk`) |

Constants are `module_model_constants` evaluated in the working precision
(src/phys/driver/wrf_constants.hpp). All fields are allocated once in the constructor; the Noah
tables are read once (VEGPARM/SOILPARM/GENPARM, incl. the LCZ categories) and the sfclayrev psi
tables are built once and uploaded to the device (`d_sfclay_tables_`). Kernels are trivially
copyable structs (`column_for`), no host pointers, literals as `real(...)`; the FP32 build
compiles (`-DCPPWRF_DOUBLE=OFF`).

**Not ported** (the constructor rejects these options): fractional_seaice = 1,
seaice_thickness_opt / seaice_snowdepth_opt = 1 (need ICEDEPTH/SNOWSI), urban physics, MYJ-type
QZ0 test (`myj`), the obs-nudging diagnostics URATX/VRATX/TRATX (not used by the physics).
phy_init is not part of the suite: the caller fills the initial surface state (what wrfinput
and `phy_init` provide; the test takes it from WRF after `start_domain`).

## 3. Corrections found by the test (all now bit-identical to WRF)

- **`P_PHY = grid%p_hyd`.** WRF 4.6.1 passes the hydrostatic pressure (and `P8W = p_hyd_w`) to
  surface_driver and pbl_driver; sfclayrev and YSU see p_hyd, not p+pb. The draft's "SFCLAYREV
  gets p_hyd instead of p_phy" is what WRF does; the C++ suite keeps it.
- **conv_t_tendf_to_moist** (use_theta_m = 1) was missing.
- **YSU** (src/phys/bl_ysu/ysu.hpp): the `ysu_topdown_pblmix` PBL-height pass runs for every
  column in WRF (it can switch pblflg on), the port only ran it where pblflg was already true
  -> e.g. PBLH 61 m vs 97 m and 518 m vs 69 m in columns with zero or negative buoyancy flux
  (the likely source of the 10-50 % PBLH/RTHBLTEN differences reported for the draft); RTHBLTEN is `(ttend*pi)/pi` as in the WRF wrapper;
  `-3.*x**2./hpbl**2.` and `ust**2` keep WRF's grouping. The YSU wind/temperature outputs go to
  grid%wstar_ysu / delta_ysu.
- **Noah SFLX** (src/phys/lsm_noah/noah.hpp): `T14 = T1*T1; T14 = T14*T14` grouping in SNOPAC;
  rdlai2d / usemonalb are namelist options (WRF default .false.: LAI and background albedo from
  the table, ALB is INOUT and becomes ALBBCK); the standalone "skip vegtyp 16 / soil 14" test
  must not apply in the driver (MODIS category 16 is barren land): the driver passes an
  `active` mask instead.
- CHS, CHS2, CQS2 are locals of surface_driver in WRF; CHS is compared through SFCEXC (= CHS),
  CHS2/CQS2 through T2/Q2 (SFCDIAGS divides by them) and Q2SAT/QGH through the Noah outputs.

## 4. Differential test (tests/phys/suite)

Reference: **unmodified** WRF v4.6.1 compiled by `tests/phys/suite/build_wrf_r8.sh`
(`configure -r8`, option 32 serial gfortran, the configure flags except
`-ftree-vectorize` -> `-fno-tree-vectorize -ffp-contract=off`: with -ftree-vectorize gfortran
calls glibc's libmvec vector math (`_ZGVbN2vv_pow`, `_ZGVbN2v_exp`, `_ZGVbN2v_log`), which
differs from scalar libm by ulps; the other differential tests of this repository use
-fno-tree-vectorize as well).

`fortran/suite_ref.F90` links against the WRF objects and calls WRF's own code only:
the framework setup of `wrf_init` (namelist.input), `start_domain` (phy_init, table reads),
and 3x `first_rk_step_part1`, `first_rk_step_part2`, `phy_prep_part2` exactly as solve_em calls
them. **No stub, no changed WRF statement.** The driver replaces only (a) `med_initialdata_input`:
the fields of a wrfinput are copied into `grid%` from `suite_input.bin` (plus the global
attributes MMINLU = MODIFIED_IGBP_MODIS_NOAH, ISWATER 17, ISICE 15, ISURBAN 13, map_proj 1) and
(b) the dynamics between the calls: a deterministic change of t_2, u_2, v_2, qv, mu_2, and
RAINNCV per call (precipitation, since mp = 0). Every input of each call is dumped (`INn:*`) and
read by the C++ driver; the surface state is taken from WRF only at call 1 and then carried by
the C++ suite itself.

Tile (`generate_inputs.py`): 24x20x44, p_top 50 hPa, 15 km, dt 60 s, ridge up to 1.9 km, 96 open
water, 36 sea-ice and 21 land-ice (glacier) points, 140 snow points, 17 MODIS land categories
(urban 13, barren 16) and 17 soil types (incl. a land point with soil 14), stable (snow/ice
inversions) and unstable (warm land, sea) columns, frozen and unfrozen soil, rain on snow
and freezing rain (T1 <= 273.15).

```bash
# once (about 30 min): unmodified WRF v4.6.1, real*8, scalar libm
tests/phys/suite/build_wrf_r8.sh /path/to/wrfref
# the test (needs gfortran, g++, python3 + numpy); CPU FP64 compared, FP32 compiled
WRF_DIR=/path/to/wrfref/WRF tests/phys/suite/run_test.sh [--markdown]
# or with CMake: -DCPPWRF_DIFF_TESTS=ON, WRF_DIR exported, ctest -L diff (diff_phys_suite)
```

Pass criterion: every element bit-identical or `|c - f| <= 1e-12 max(|c|,|f|)`.
Result (CPU FP64, 3 calls, P1_n = after part1, P2_n = after part2, OUTn = after phy_prep_part2):
94/94 fields bit-identical (`compare.py` prints the per-field table).
