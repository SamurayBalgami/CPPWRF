# Physics: ports, driver glue and verification

This document describes how the WRF v4.6.1 physics is reproduced in CPPWRF, scheme by scheme, and
how each part is verified. Sections 1-5 cover the column-physics driver (surface layer, land surface,
PBL); sections 6-13 add microphysics, cumulus and radiation; section 14 lists the WRF behaviours that
had to be reproduced for bit-identical results. Case names (P1-P8) refer to the reference runs listed
in docs/STATUS.md.

## 0. The P1 driver (sfclayrev + Noah + YSU)

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
seaice_thickness_opt / seaice_snowdepth_opt = 1 (need ICEDEPTH/SNOWSI), urban physics, the obs-nudging diagnostics URATX/VRATX/TRATX (not used by the physics).
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

## 5. Model wiring and real-data acceptance

`phys.sf_sfclay_physics / sf_surface_physics / bl_pbl_physics` in the config (real data only)
make `Model::build` call `init_real_physics` (src/cases/real_physics.cpp): the suite is created,
its initial surface state filled as WRF start_domain does for a non-restart run (the wrfinput
fields; phy_init: UST = 1e-4, IVGTYP = nint(LU_INDEX), DX2D = dx, landuse_init from LANDUSE.TBL
for the season of JULDAY; bl_init/LSMINIT: SH2O with FRH2O for frozen layers, CANWAT = 0;
CTOPO = CTOPO2 = 1), and attached to the RK3 integrator (`set_column_physics`). The surface
fields and R*BLTEN, EXCH_H, TSLB, SMOIS, SH2O are written to the native snapshots.

**Correction found by the real-data run:** `conv_t_tendf_to_moist` gets `grid%t_1`, which
small_step_prep refreshes only after first_rk_step_part2: in step n it is theta_m at the start of
step n-1 (the initial state in step 1). The suite keeps that lagged copy (`wrf_t1_lag`); with the
dycore's current t_1 step 2 already differed by 1e-9 (U) from WRF. The differential test passes
WRF's t_1 itself and sets `wrf_t1_lag = false`.

Acceptance against WRF v4.6.1 real*8, run P1_r8 (`configs/real_p1_r8.ini`, tests/verify/wrf_suite.py):
CPU FP64 is bit-identical to WRF for 1 h (every compared field, including the surface state). CUDA FP64
agrees to roundoff (U 1.5e-9, THM 1.4e-10 relative L2 at 1 h); the land surface amplifies the
difference between the CUDA and glibc elementary functions over the hour. CUDA FP32 against
single-precision WRF: U 1.9e-4, THM 7.5e-7.

## 6. Microphysics with the column physics (P5 = P1 + WSM6 + diff_opt 1, km_opt 4)

Reference: the WRF v4.6.1 real*8 run P5_r8 (the P4 inputs with ra_lw/ra_sw/cu_physics = 0; history
every step; `configs/real_p5_r8.ini`). CPU FP64 is bit-identical to WRF for 1 h, including all water
species. WSM6 amplifies roundoff strongly once ice is present (ice processes near saturation
thresholds), so the CUDA builds differ from WRF by more than in the dry case. The column replay
(tests/phys/mp_wsm6) runs CPPWRF's WSM6 inputs through WRF's compiled mp_wsm6: identical output.

## 7. Kain-Fritsch cumulus (P2, `cu_physics = 1`)

### What WRF executes for `cu_physics = 1` (`KFETASCHEME`)

1. **`cu_init` / `kf_eta_init`** (`phys/module_physics_init.F:1945`, `phys/module_cu_kfeta.F:3084`):
   - `stepcu = max(nint(cudt * 60. / dt), 1)`
   - `NCA = -100.`, `W0AVG = 0.`, `RTHCUTEN = RQVCUTEN = RQCCUTEN = RQRCUTEN = RQICUTEN = RQSCUTEN = 0.`
   - `kf_lutab(SVP1, SVP2, SVP3, SVPT0)`: builds lookup tables `TTAB(250,220)`, `QSTAB(250,220)`, `THE0K(220)`, `ALU(200)`, `RDPR`, `RDTHK`, `PLUTOP` (`src/phys/cu_kf/kf_tables.cpp`, uploaded once to `d_kf_tables_`).
2. **`first_rk_step_part1` -> `cumulus_driver`** (`dyn_em/module_first_rk_step_part1.F:1376`, `phys/module_cumulus_driver.F:1015`):
   - Run decision (`module_cumulus_driver.F:830-864`): `itimestep == 1 || cudt == 0 || stepcu == 1 || itimestep % stepcu == 0`.
   - Passes `P = grid%p_hyd`, `PI = pi_phy`, `T = t_phy`, `TH = th_phy`, `U = u_phy`, `V = v_phy`, `W = grid%w_2`, `RHO = grid%rho`, `DZ8W = dz8w`, `HTOP = grid%cutop`, `HBOT = grid%cubot`, `WARM_RAIN = grid%warm_rain` (`mp_physics == 1`).
   - Inside `KF_eta_CPS` (`phys/module_cu_kfeta.F:29`):
     - Updates running-mean vertical velocity every cumulus step: `W0 = 0.5*(w(k) + w(k+1))`, `W0AVG = (W0AVG * (2*STEPCU - 1) + W0) / (2*STEPCU)`.
     - Checks `NCA(i,j) >= 0.5*DT`: if active from a prior deep-convection event, returns early keeping existing `R*CUTEN`, `RAINCV`, `PRATEC`.
     - Otherwise zeroes `R*CUTEN`, `RAINCV`, `PRATEC`, `CLDFRA_DP/SH`, `QC_CU/QI_CU`, sets `CUTOP = 1`, `CUBOT = KTE + 1`, and runs `KF_eta_PARA` (`src/phys/cu_kf/kf.hpp`).
     - Note on `mp_physics = 0`: `RQCCUTEN` is passed unconditionally in `cumulus_driver`, so `PRESENT(rqccuten)` is true in `KF_eta_CPS` and `RQCCUTEN = DQRDT + DQCDT` (since `f_qr = .false.`), while `calculate_phy_tend`, `update_phy_ten`, and `phy_prep_part2` guard `RQCCUTEN` with `P_QC >= PARAM_FIRST_SCALAR` (`s.P_QC >= 0`).
3. **`first_rk_step_part2`** (`src/phys/driver/tendency_coupling.cpp`):
   - `calculate_phy_tend` (`dyn_em/module_em.F:2208-2261`): multiplies `RUCUTEN, RVCUTEN, RTHCUTEN, RQVCUTEN` (and `RQCCUTEN, RQRCUTEN, RQICUTEN, RQSCUTEN` when the corresponding `P_Q* >= 0`) by `mu = c1h(k)*mut + c2h(k)`.
   - `update_phy_ten` / `phy_cu_ten` (`phys/module_physics_addtendc.F:1185-1225`): adds `RTHCUTEN` to `t_tendf`, `RQVCUTEN` to `moist_tend(P_QV)` (and `RQCCUTEN, RQRCUTEN, RQICUTEN, RQSCUTEN` when `P_Q* >= 0`) before `conv_t_tendf_to_moist`.
4. **After the RK loop (`solve_em.F:3560`, before `phy_prep_part2`)**:
   - `advance_ppt` (`phys/module_physics_addtendc.F:2059-2235`, `PhysicsSuite::advance_ppt` / `KAdvancePpt` in `src/phys/driver/suite.cpp`):
     - `RAINCV = DT * PRATEC`, `RAINSHV = DT * PRATESH`
     - `RAINC += DT * PRATEC`, `RAINSH += DT * PRATESH`, `CUPPT += DT * PRATEC / 1000`
     - Updates `HTOP = max(CUTOP, HTOP)` and `HBOT = min(CUBOT, HBOT)` where `RAINCV > 0`
     - For `KFETASCHEME`: if `NCA(i,j) > 0`, checks `nint(NCA(i,j)/DT) <= 1` (zeroing `RTHCUTEN, RQVCUTEN, RQCCUTEN, RQRCUTEN, RQICUTEN, RQSCUTEN` when convection expires) and decrements `NCA(i,j) -= DT`.
   - `phy_prep_part2` (`dyn_em/module_big_step_utilities_em.F:5105-5167`): decouples `R*CUTEN /= (c1h(k)*muts + c2h(k))`.
   - On the following step, `surface_driver` (`KSfcAlways`) adds the previous step's `RAINCV` into `RAINBL` before `cumulus_driver` runs.

### Real-data acceptance against WRF v4.6.1 (`run_P2_r8` and `run_P2`)

CPU FP64 is bit-identical to WRF for 1 h (state, RAINC, RAINCV, NCA, W0AVG and the cumulus tendencies).
CUDA FP64 at 1 h: U 5.2e-5, THM 9.4e-7 relative L2; CUDA FP32 against single-precision WRF: U 8.4e-4.

**Why the GPU results drift (a 1-ulp trigger):**
In `KF_eta_PARA` (`module_cu_kfeta.F:929-1022`), whenever a 2-layer updraft source layer (`KPBL = LC + 1`) is saturated (`q0 = qes`, so `TLCL = min(TLCL, TMIX) = TMIX`) and has `WKL < 1e-4` (`DTLCL = 0`, `DTRH = 0`), the parcel LCL height is `ZLCL = ZMIX = (dp_0 z_0 + dp_1 z_1)/(dp_0 + dp_1)`, `DLP = (ZLCL - z_0)/(z_1 - z_0) = dp_1/(dp_0 + dp_1)`, and `TENV = T_0 + (T_1 - T_0)*DLP = TMIX = TLCL`. Thus `TLCL + DTLCL + DTRH - TENV == 0` **identically in exact arithmetic**, and the buoyancy test `IF (TLCL + DTLCL + DTRH .LT. TENV) CYCLE usl` tests `0 < 0` up to 1 ULP (`5.68e-14 K`) of floating-point rounding. At step 6 (`t = 300 -> 360 s`), column `(j=77, i=8)` at `lc = 8` satisfies this exact degeneracy (`0.0` vs `-5.68e-14 K`), so that any rounding difference in the inputs (as between the CUDA and glibc math libraries) can switch shallow convection on or off in that column. The CPU build reproduces WRF's arithmetic and makes the same decision.

## 8. Radiation (P3 = P1 + RRTMG LW/SW)

`src/phys/driver/radiation.{hpp,cpp}` ports radiation_driver for ra_lw_physics = ra_sw_physics = 4
with o3input = 2 (CAM ozone: `src/phys/ra_rrtmg_common/cam_ozone.*` = oznini, ozn_time_int,
ozn_p_int), ghg_input = 1 (`cam_gases.*` = read_CAMgases), aer_opt = 0, swint_opt = 0: radconst,
calc_coszen at xtime + radt/2, RTHRATEN = LW (+ SW in daytime), GSW, GLW, SWDOWN = GSW/(1-ALBEDO);
P_PHY = p_hyd, P8W = p_hyd_w as in first_rk_step_part1. The RRTMG ports were checked
column by column against dumps of a debug build of WRF (the
RRTMG_LWRAD/RRTMG_SWRAD arguments and the rrtmg_lw/rrtmg_sw core arguments, 10591 columns, the SW
at 09 UTC of a 9 h run): **LW and SW are bit-identical to WRF** in both the cores and the wrappers
after these corrections:
- LW band 16 upper atmosphere: WRF's nspb(16) = 0, so absb(1..2, ig) is used for every layer;
- gfortran evaluates `p**4` as `(p*p)*(p*p)` (36 places in taumol); the SW band-24 sum order
  (ozone before the continuum); inatm `amttl = amttl + coldry + wkl` order;
- NLAYERS = kme + nint(p_top/400 Pa) - 1 with kme = nz + 1 (one layer was missing);
- o3input = 2 ozone and ghg_input = 1 gases (the port had the climatology / fixed gases only);
- SW core as WRF's rrtmg_sw: zepzen 1e-10, band albedos, adjflux = adjes*scon/rrsw_scon, the
  heating rate in the core with heatfac/pdp and the top layer set to 0, RTHRATENSW =
  (swhr/86400)/pi3d.
CAM ozone is bit-identical on all columns, and so is COSZEN (calc_coszen's sin/cos are libmvec lanes,
section 14). Tables: a directory with rrtmg_lw/sw_tables.bin (dumped from WRF's
rrtmg_lw/swinit with the _DBL data, ozone*.formatted, CAMtr_volume_mixing_ratio; config key
phys.rrtmg_tables). Execution: section 10.
Clouds: section 9. Not yet: aer_opt, swint_opt.
Acceptance (wrf_suite.py): real_p3_r8_fp64_cpu bit-identical (steps 1-3); real_p3_fp32_1h U 1.9e-4,
THM 7.4e-7.

## 9. Clouds in the radiation (P4 = P2 + WSM6 + RRTMG + diff_opt 1 / km_opt 4)

Reference: the WRF v4.6.1 real*8 run P4_r8 (kf + WSM6 + RRTMG radt = 1,
icloud = 1, use_mp_re = 1 (default); steps/ = 3 minutes with history every step and RQRCUTEN,
RQICUTEN, RQSCUTEN, RE_CLOUD, RE_ICE, RE_SNOW in the history; `configs/real_p4_r8.ini`).

What WRF executes (and the port):
- **phy_init**: with use_mp_re = 1, RRTMG LW + SW and mp_physics = WSM6 (also WSM3/5/7, WDM*,
  Thompson, ...) `has_reqc = has_reqi = has_reqs = 1`; RE_CLOUD/RE_ICE/RE_SNOW start at
  RE_QC_BG/RE_QI_BG/RE_QS_BG = 2.49/4.99/9.99 um (`PhysicsParams::microp_re`, set in model.cpp;
  `PhysicsCoupling::re_cloud/re_ice/re_snow`).
- **module_mp_wsm6.F after every WSM6 call**: `mp_wsm6_effectRad_run` with t = th*pii of the new
  column, the background values and RE_Q*_MAX = 50/125/999 um (`Wsm6EffRadKernel`,
  src/phys/mp_wsm6/wsm6.hpp, run on the microphysics columns only). The per-row has_qc/qi/qs
  flags only skip work.
- **radiation_driver**: CLDFRA = 0, then with icloud = 1 and F_QC or F_QI `cal_cldfra1`
  (Randall / Hong-Pan, t = t_phy, p = p_hyd; `driver::cal_cldfra1`); icloud_cu = 0, no KF cloud
  fraction. CLDFRA is a suite field (wrfout, restart).
- **RRTMG_LWRAD / RRTMG_SWRAD** (`rrtmg_lw_column` / `rrtmg_sw_column` take re_cloud/re_ice/re_snow,
  nullptr when has_req* = 0): qc/qi/qs clipped at 0; inflg = 3 with re_cloud (recloud =
  max(2.5, re*1e6), 10.5 over ocean / 7.5 over land when <= 2.5 in a cloudy layer), inflg = iceflg
  = 4 with re_ice (max(5, re*1e6), retab with index <= 75 when <= 5 in a cloudy layer), inflg =
  iceflg = 5 with re_snow (max(10, re*1e6)); the P3 special case (no re_snow) is ported too. With
  iceflg >= 4 the ice path excludes the snow; with iceflg = 5 the snow path uses snow_mass_factor
  (resnow > 130 um: (130/resnow)^2 and resnow = 130; WRF's `gicewp = gicewp + 1% of the snow`
  happens after cicewp is set and has no effect). reliq/reice from relcalc/reicalc are replaced
  by the microphysics sizes (inflg >= 3, iceflg >= 4); the Fu 1.0315 factor only for iceflg = 3.
  The cores take inflag/iceflag/liqflag (they were fixed at 2/3/1).
- **SW cldprmc_sw**: the asymmetry factor uses WRF's grouping
  `scatice*((gice - forwice)/(1 - forwice))` (and the same for snow); the port had
  `(scatice*(gice - forwice))/(1 - forwice)`.

Verification:
- P4 (wrf_suite real_p4_r8_fp64_cpu and real_p4_r8_fp64_cpu_1h): CPU FP64 bit-identical to WRF for
  1 h, including every physics tendency, the water species, CLDFRA and the effective radii.
- P4 1 h CUDA FP64 vs WRF (real_p4_r8_fp64_cuda_1h): U 2.2e-4, THM 2.7e-6, QV 2.1e-4 (rel L2), the
  size of P2's (the cumulus trigger, section 7).
- Restart (P4, CPU FP64): step 3 straight vs restarted after step 1 bit-identical (CLDFRA,
  RE_* carried).
- Cloudy RRTMG column by column: debug-WRF dumps of RRTMG_LWRAD / RRTMG_SWRAD at 09 UTC of a 9 h
  P4 run (5431 cloudy columns, WSM6 radii, snow): LW and SW wrappers
  **bit-identical to WRF on 10591/10591 columns**, and cal_cldfra1 equals WRF's CLDFRA exactly.
  GPU batch path on the same dumps: LW 8884, SW 10490 columns bit-identical, the rest at roundoff
  (GLW 1e-16, RTHRATENLW 6e-15, GSW 3e-16, RTHRATENSW 3e-14).

## 10. RRTMG on the GPU

RRTMG_LWRAD / RRTMG_SWRAD are written as stages over strided per-column views
(`src/phys/ra_rrtmg_common/strided.hpp`: x[i] = p[i*s]; one column in local arrays with s = 1, or a
batch of columns in device memory with the column index fastest, s = ncol):

| stage | LW (`rrtmg_lw_stages.hpp`) | SW (`rrtmg_sw_stages.hpp`) | GPU kernel over |
|---|---|---|---|
| prep | RRTMG_LWRAD profiles, buffer layers, inirad/CAM ozone, cloud properties, McICA subcolumns, inatm, setcoef, secdiff | RRTMG_SWRAD profiles, ozone, cloud properties, McICA, inatm_sw, setcoef_sw, taumol_sw | column |
| layer | cldprmc, taumol (one layer, all bands), cloudy-layer flag | - | column x layer |
| band / g-point | rtrnmc's g-point loop of one band (band sums of every level) | spcvmc_sw for one g-point (its cldprmc_sw, delta scaling, reftra, vrtqdr), contributions to every level | column x band (16) / column x g-point (112) |
| level | band sums in band order (fluxfac, fnet) | g-point sums in WRF's order (also UV/visible, near-IR) | column x level |
| finish | heating rates, GLW, RTHRATENLW | heating rates, GSW, RTHRATENSW | column |

Every value keeps WRF's arithmetic and summation order (the McICA random numbers are drawn
subcolumn by subcolumn in WRF's order; the subcolumn water paths are the layer paths where the
subcolumn is cloudy; WRF's clear-sky assignments clrdrad = drad stay inside the band loop).
CPU builds run the stages column by column (`rrtmg_lw_column` / `rrtmg_sw_column`, OpenMP over
columns): **bit-identical to WRF** (P3 night and 09 UTC dumps, 10591/10591 columns, cores and
wrappers). CUDA builds (`RrtmgLwBatch`, `RrtmgSwBatch`, `RadiationDriver` device path: inputs
staged once, coszen, CAM ozone and cal_cldfra1 per column on the device, no host transfers) run in
chunks of 2048 columns (~0.7 GB of work arrays) and are compiled with `--fmad=false` (the table
indices int(tblint*x + 0.5) jump with contracted products): on the same dumps 10484/10591 (LW) and
10488/10591 (SW) columns bit-identical to WRF, the rest at roundoff (GLW 3e-16, RTHRATENLW 7e-15,
GSW 2e-15, RTHRATENSW 2e-13; device exp/log). The batch path also runs on the CPU
(tests/phys/ra_rrtmg_{lw,sw}/batch_check.cpp: bit-identical, but 2.4x slower than the column path
there because of the column-fastest layout).

Speed (RTX 2060 mobile, FP64, 119 x 89 x 44, one call for all 10591 columns): LW 0.17 s
(band 40 %, layer 33 %, prep 25 %), SW 0.37 s in daytime (g-point 83 %); the host path took
3.4 s (LW) single-threaded. Model, radt = 1 (radiation every step): P3 CUDA FP64 3200 ->
510 ns/cell/step, P4 1 h 3234 -> 650 ns/cell/step.

## 11. Thompson microphysics (P6, mp_physics = 8)

Port: `src/phys/mp_thompson/` (column scheme `thompson.hpp`, tables `thompson_tables.*`), coupled in
`src/dyn/physics_coupling.cpp` with the scalars qni, qnr (`[moist] species = ...,qni,qnr`); column
differential test against WRF's Fortran: tests/phys/mp_thompson (<= 1e-12).

**WRF reference build.** The stock real*8 WRF (`configure -r8`: PROMOTION = -fdefault-real-8 only)
promotes DOUBLE PRECISION to 16-byte quad in gfortran. Thompson declares its tables and many locals
DOUBLE PRECISION, so that build computes them in software quad: qr_acr_qg alone did not finish in
8 h. The Thompson references are therefore made with a WRF build from the same source
with `-fdefault-real-8 -fdefault-double-8` (DOUBLE PRECISION stays 8 bytes, as the author intends
and as in the single-precision WRF; qr_acr_qg takes 280 s). The schemes verified before (WSM6, KF,
YSU, Noah, sfclayrev, RRTMG, Kessler) have no DOUBLE PRECISION in their used code paths.

**Tables.** WRF writes qr_acr_qg_V4.dat, qr_acr_qsV2.dat and freezeH2O.dat on its first run and reads
them afterwards. CPPWRF reads them from `phys.table_dir` (default real.dir; either byte order) and
then uses WRF's values bit for bit; without them it computes the tables (OpenMP, ~2.4 min on 3
threads), which agree with WRF's to 1e-15 (a few ulp, 0 for the sacr1 tables).

**Acceptance** (P6 = P5 with mp_physics = 8: P1 column physics + Thompson + diff_opt 1 / km_opt 4 on
the P4 inputs; reference run P6_r8, FP32: P6): CPU FP64 bit-identical to WRF for 1 h, including
QRAIN, QICE, QNICE and QNRAIN. CUDA FP32 against single-precision WRF, 1 h: U 2.2e-4, THM 2.2e-6.
The column replay (CPPWRF's Thompson inputs through WRF's compiled mp_gt_driver) agrees bit for bit.
wrf_suite: real_p6_r8_fp64_cpu, real_p6_r8_fp64_cpu_1h, real_p6_fp32_1h.

## 12. New Tiedtke cumulus (P7, cu_physics = 16)

Port: `src/phys/cu_ntiedtke/ntiedtke.hpp` = WRF 4.6.1 `phys/physics_mmm/cu_ntiedtke.F90` for one column
(cumastrn and the routines it calls, statement for statement; integer powers as products, cutypen's
INTEGER `zcbase`, cuflxn's `ktopm2 = 2`). The flags WRF computes over a row of columns (llo3, the
early exits) do not depend on the other columns: a whole row (lq = 119) through WRF's object gives
the same result as one column at a time. Driver glue (`src/phys/driver/suite.cpp`, `KNTiedtke`) =
`module_cu_ntiedtke.F` cu_ntiedtke_driver: pre_run (zi from dz8w, omega = -0.5 g rho (w(k)+w(k+1)),
arrays flipped top to bottom, the forcing zeroed at itimestep 1), cu_ntiedtke_run, post_run
(R*CUTEN = (new - old) / (stepcu dt), RTHCUTEN divided by pi, RAINCV = rn / stepcu, PRATEC).

**Advective forcing (RTHFTEN, RQVFTEN).** WRF hands the scheme the large-scale tendencies:
- `rk_tendency`, `set_tend(RTHFTEN, t_tend, msfty)`: the theta_m advection only, every RK stage, so the
  last stage counts (`LargeStep` `KTendT`, hook `ColumnPhysics::th_forcing`);
- `rk_scalar_tend`, `set_tend(RQVFTEN, advect_tend, msfty)` for qv (`ScalarTransport`, hook
  `qv_forcing`; the positive-definite tendency in the last stage with moist_adv_opt = 1);
- on a specified boundary the advection routines still form a tendency in the spec zone: no x flux
  divergence in the outer columns, no y divergence in the outer rows, the vertical term everywhere.
  The C++ forms exactly that for the forcing (the dycore's own spec-zone values are untouched); the
  outer column feeds RUCUTEN at i = 1 through add_a2c_u;
- `phy_prep_part2`: both divided by (c1 muts + c2); with use_theta_m = 1 RTHFTEN becomes a dry-theta
  tendency, th/(t+t0) (RTHFTEN - rv/rd th RQVFTEN);
- `cumulus_driver`: RTHFTEN = (RTHFTEN + RTHRATEN + RTHBLTEN) pi, RQVFTEN += RQVBLTEN.
Coupling: `phy_cu_ten` adds RTHCUTEN, RQVCUTEN, RQCCUTEN, RQICUTEN (add_a2a) and RUCUTEN / RVCUTEN
(add_a2c_u / add_a2c_v); `advance_ppt` only accumulates RAINC (CASE DEFAULT). RTHFTEN and RQVFTEN are
restart fields (restart_p7_fp64: bit-identical).

**Column check** (tests/phys/cu_ntiedtke): the scheme's inputs and outputs dumped from CPPWRF at steps 2
and 4 of P7, replayed through WRF's compiled cu_ntiedtke_run (r8d8): 10591/10591 columns agree.

**Acceptance** (P7 = P2 with cu_physics = 16, reference run P7_r8, FP32: P7): CPU FP64 bit-identical to
WRF for 1 h. CUDA FP64 at 1 h: U 4.7e-7, THM 6.0e-8; CUDA FP32 against single-precision WRF: U 4.9e-4.
wrf_suite: real_p7_r8_fp64_cpu, real_p7_r8_fp64_cpu_1h, real_p7_r8_fp64_cuda_1h, real_p7_fp32_1h,
restart_p7_fp64.

## 13. Grell-Freitas-Li cumulus (P8, cu_physics = 3, GFLSCHEME)

Port: `src/phys/cu_gfl/` = WRF 4.8.0 `phys/GFL/module_cu_gfl_deep.F90` (`cu_gfl_deep_run`),
`phys/GFL/module_cu_gfl_sh.F90` (`cu_gfl_sh_run`, `ishallow = 1`, `ichoice = 3`) and
`phys/GFL/WRF/module_cu_gfl_wrf.F90` (`cu_grell_freitas_li`) with WRF's settings (`imid = 0`,
`sub3d = 0`, `dicycle = 0`, `irainevap = 1`, `icoldpool = 1`, `autoconv = 1`, `aeroevap = 1`,
`ccn = 1500.`, `ccnclean = 250.`, no stochastic perturbations, `clos_choice` = `phys.clos_choice`).
Driver glue `KGrellFreitas` in `src/phys/driver/suite.cpp`; it uses the same `RTHFTEN`/`RQVFTEN`
forcing as New Tiedtke (`cu_grell_freitas_li` combines it with `RTHRATEN`/`RTHBLTEN` itself). All
columns in the domain are active (`its..itf`, `jts..jtf`), `RUCUTEN`/`RVCUTEN` are computed (including
cold-pool gust-front momentum impact and ECMWF-style kinetic energy dissipation heating) but WRF's
`phy_cu_ten` (`CASE (GDSCHEME, G3SCHEME, GFLSCHEME)`) does not add them to the model momentum
tendencies, and `RQCCUTEN`/`RQICUTEN` split at `T = 258 K`.

Exactness details in `src/phys/cu_gfl/`:
- **Vector pow and operation order.** `gfortran` vectorises the `zu` profile loop (`get_zu_zd_pdf_fim`)
  and calls glibc's `libmvec` `_ZGVbN2vv_pow` for pairs of iterations (scalar `pow` for an odd remainder
  iteration), evaluating `(((x**(a - 1.) * b) / beta) * (1. - x)**(b - 1.))`; it can differ from scalar
  `pow` or left-to-right multiplication in the last bit. The CPU path calls the same `libmvec` function
  (`cppwrf_core` links `libmvec`) and matches Fortran's evaluation order (`rates_up_pdf` `smooth_Heaviside`
  and `irainevap`'s `edt * evfact * (sig_init * sig_init)`); the GPU uses its own `pow`.
- **Column-local inversion layers.** Unlike WRF 4.6.1 GF, `get_inversion_layers` in GFL v4.8.0
  initialises `sec_deriv(:) = 0.` at entry, starts at `ilev = max(4, kstart(i)+1)`, and resets
  `sec_deriv(:) = 1.e9` inside the column loop, so there is no row-carried state or uninitialised read
  and the entire scheme runs in a single parallel pass over all columns.
- **Constants.** `cu_grell_freitas_li` receives WRF's model constants (`xlv = 2.5e6`, `cp = 1004.5`,
  `g = 9.81`, `r_v = 461.6`) as `real(kind=kind_phys)` for `dhdt`, `omeg`, and `mconv`, while
  `module_cu_gfl_deep` and `module_cu_gfl_sh` use their own module parameters (`g = 9.81`, `cp = 1004.`,
  `xlv = 2.5e6`, `r_v = 461.`, `tcrit = 258.`, `fluxtune = 1.5`).

**Column replay** (`tests/phys/cu_gfl`): `cu_grell_freitas_li`'s inputs dumped from CPPWRF at P8 steps
1, 5 and 10 run through WRF v4.8.0's `r8d8` objects: every column bit-identical (deep only and
shallow + deep).

**Acceptance** (P8 = P2 with `cu_physics = 3`, `ishallow = 1`; reference run `P8_r8`, FP32: `P8`): CPU
FP64 bit-identical to WRF for 10 min (`steps/`) and 1 h. CUDA FP64 at 1 h: U 3.6e-10, W 2.0e-8,
THM 1.2e-10, QV 5.4e-9; CUDA FP32 against single-precision WRF: U 4.3e-4, W 2.2e-2, THM 9.8e-6,
QV 2.9e-4. `wrf_suite`: `real_p8_r8_fp64_cpu_10min`, `real_p8_r8_fp64_cpu_1h`,
`real_p8_r8_fp64_cuda_1h`, `real_p8_fp32_1h`, `restart_p8_fp64`.

## 14. WRF behaviours reproduced for bit identity (P4-P6, nesting)

Found with WRF checkpoint dumps (tests/verify/wrf_dycore_dumps) and column replays through WRF's own
r8d8 objects (WSM6: mp_wsm6; Thompson: mp_gt_driver; Noah: SFLX with SOIL_VEG_GEN_PARM for MODIS).

- Diffusion (diff_opt 1, km_opt 4): WRF's left-to-right products in cal_deform_and_div (defor12
  `dhat*0.25*tmpz*(sum rdzw)` with each term's own rdzw order, defor13/23 `(u-u_base-u(k-1)+u_base(k-1))
  *0.5*(rdz+rdz)`), smag2d mlen_h `((dx/msftx)*dy)/msfty`, horizontal_diffusion_3dmp
  `f(i+1)-f(i)-b(i+1)+b(i)`.
- Microphysics coupling: moist_physics_prep_em and the post-microphysics calc_p_rho_phi run on the tile
  without the spec zone, so the libmvec pow lanes of pii start at i = 1.
- The LSM gets the microphysics' RAINNCV and SR (they were never passed before); SNOPAC's SNCOVR = 0 is
  written back; WDFCND `(500*SICEMAX)**3.` is a libm pow.
- Kain-Fritsch: no clipping of PPTFLX, PRATEC, RAINCV at 0 (WRF has none; a -2e-18 RAINCV reaches the
  LSM).
- Thompson: libmvec pairs in the rain lambda loops (reversed k) and the SLP-paired
  `(lamr+fv_r)**(-cre)`; table bins `xDx(2)` folded at compile time (correctly rounded) and `xDx(3:)`
  through vector exp; `10.**nn` and `Dc**nu_c` with integer exponents are libgcc's __powidf2
  (cppwrf::powi); the rain-graupel tables (dimNRHG = 1 when not hail-aware) are indexed with
  idx_bg = 5 as WRF does, i.e. 4 steps along idx_r1 (an out-of-bounds read of WRF's, clamped at the end);
  operation orders of the autoconversion, tcond and the number-tendency limits.
- Radiation: calc_coszen's sin/cos are libmvec lanes of the 0..nx-1 loop.
- Nesting: med_initialdata_input runs start_domain(nest) after reading wrfinput, so adjust_tempqv sees
  the recomputed p; lbc_fcx_gcx `0.1/dt*(n)/(relax_zone-1)` left to right; couple_or_uncouple_em reads
  mu's halo at ids/jds; no halo refresh after feedback (update_after_feedback_em is a DM exchange only);
  history of a nest is written in med_before_solve_io, i.e. after the parent's next med_nest_force.
