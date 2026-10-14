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

| run | U | W | THM | QV | surface (HFX, UST, PBLH, TSK, T2, TSLB, R*BLTEN, EXCH_H) |
|---|---|---|---|---|---|
| CPU FP64, 3 steps (history every step) | 2.2e-14 | 4.6e-12 | 1.8e-16 | 1.0e-15 | <= 8e-12 (max rel) |
| CUDA FP64, 1 h | 6.0e-10 | 1.8e-08 | 1.2e-10 | 3.1e-09 | TSK 6e-5 max rel at 1 h |
| CUDA FP32 vs WRF FP32 (P1), 1 h | 1.9e-04 | 3.9e-03 | 7.5e-07 | 4.9e-05 | |

(relative L2 over the domain unless noted.) Over an hour the land surface amplifies roundoff:
CPPWRF CPU FP64 vs CUDA FP64 differ by the same amount as either does from WRF (TSK 7.5e-5,
TSLB 1e-4 max rel at 1 h, from 5e-13 at 10 min), so this is the sensitivity of the problem,
not a port difference.

## 6. Microphysics with the column physics (P5 = P1 + WSM6 + diff_opt 1, km_opt 4)

Reference: WRF v4.6.1 real*8 run ~/wrf_real/run_P5_r8 (run_P4_r8 inputs with
ra_lw/ra_sw/cu_physics = 0; history every step; `configs/real_p5_r8.ini`). Step 1 is at
roundoff (U 7e-15, QCLOUD 3e-12, QICE 5e-13; wrf_suite.py real_p5_r8_fp64_cpu). From step 2,
QICE differs by ~4e-5 (rel L2). Proven not to be a port difference: a debug build of WRF dumped
its WSM6 inputs at step 2; they equal CPPWRF's to 5e-12 (th 2e-15, qi 2.5e-12, rho 6e-14), the
scalar Fortran driver (tests/phys/mp_wsm6) reproduces WRF's QICE bit for bit from WRF's inputs
and CPPWRF's WSM6 bit for bit from CPPWRF's inputs. WSM6 amplifies the input roundoff by ~1e8
(ice processes at thresholds once ice exists): CPPWRF CPU FP64 vs CUDA FP64 differ by the same
7e-5 at step 2. Microphysics runs can therefore only be checked at roundoff for the first step;
later steps and FP32 are statistical comparisons.


## 7. Stage 2 (P2): Kain-Fritsch cumulus (`cu_physics = 1`)

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

| run | U | W | THM | QV | surface & cumulus (`RAINC`, `RAINCV`, `NCA`, `W0AVG`, `R*CUTEN`, `HFX`, `UST`, `PBLH`, `TSK`) |
|---|---|---|---|---|---|
| CPU FP64, 3 steps (`t = 60, 120, 180 s`) | 2.20e-14 | 4.60e-12 | 1.77e-16 | 1.01e-15 | `NCA` bit-exact (`0.0`); `RAINC` 1.2e-13, `RAINCV` 3.5e-13, `W0AVG` 1.5e-12, `RTHCUTEN` 8.5e-13, `RQVCUTEN` 7.4e-13, `RQCCUTEN` 2.5e-12 (max rel) |
| CUDA FP64, 3 steps (`t = 60, 120, 180 s`) | 1.79e-13 | 6.29e-12 | 2.30e-16 | 4.75e-15 | matches CPU FP64 to `1e-13` |
| CUDA FP64, 10 min (`t = 600 s`, outside `(j=77,i=8)` box) | 2.87e-13 | 6.81e-12 | 4.79e-16 | 1.43e-14 | `NCA` bit-exact across all 10,591 columns (`193` deep, `413` shallow); `RAINC` max diff `7.9e-13` |
| CUDA FP64, 1 h (`t = 3600 s`, whole domain) | 8.29e-05 | 1.74e-03 | 1.35e-06 | 5.16e-05 | `RAINC` relL2 `1.0e-4` (max rel `2.5e-4`), `TSK` relL2 `8.8e-6`, `HFX` relL2 `1.6e-4` |
| CUDA FP32 vs WRF FP32 (`run_P2`), 1 h | 7.86e-04 | 4.19e-02 | 1.09e-05 | 3.65e-04 | `TSK` relL2 `2.4e-5`, `HFX` relL2 `8.0e-4`, `UST` relL2 `8.6e-4` |

**Note on the 1-ULP shallow-convection trigger threshold at `t = 360 s` (step 6):**
In `KF_eta_PARA` (`module_cu_kfeta.F:929-1022`), whenever a 2-layer updraft source layer (`KPBL = LC + 1`) is saturated (`q0 = qes`, so `TLCL = min(TLCL, TMIX) = TMIX`) and has `WKL < 1e-4` (`DTLCL = 0`, `DTRH = 0`), the parcel LCL height is `ZLCL = ZMIX = (dp_0 z_0 + dp_1 z_1)/(dp_0 + dp_1)`, `DLP = (ZLCL - z_0)/(z_1 - z_0) = dp_1/(dp_0 + dp_1)`, and `TENV = T_0 + (T_1 - T_0)*DLP = TMIX = TLCL`. Thus `TLCL + DTLCL + DTRH - TENV == 0` **identically in exact arithmetic**, and the buoyancy test `IF (TLCL + DTLCL + DTRH .LT. TENV) CYCLE usl` tests `0 < 0` up to 1 ULP (`5.68e-14 K`) of floating-point rounding. At step 6 (`t = 300 -> 360 s`), column `(j=77, i=8)` at `lc = 8` satisfies this exact degeneracy (`0.0` vs `-5.68e-14 K`), flipping shallow convection at `lc = 8` for that single step while CPU FP64 and CUDA FP64 agree to `1e-13` and the rest of the domain matches WRF to `<= 2.9e-13` at `t = 600 s`.


## 8. Radiation (P3 = P1 + RRTMG LW/SW)

`src/phys/driver/radiation.{hpp,cpp}` ports radiation_driver for ra_lw_physics = ra_sw_physics = 4
with o3input = 2 (CAM ozone: `src/phys/ra_rrtmg_common/cam_ozone.*` = oznini, ozn_time_int,
ozn_p_int), ghg_input = 1 (`cam_gases.*` = read_CAMgases), aer_opt = 0, swint_opt = 0: radconst,
calc_coszen at xtime + radt/2, RTHRATEN = LW (+ SW in daytime), GSW, GLW, SWDOWN = GSW/(1-ALBEDO);
P_PHY = p_hyd, P8W = p_hyd_w as in first_rk_step_part1. The RRTMG ports were checked
column by column against dumps of a debug build of WRF (`~/cppwrf_scratch/wrf_dbg_src`: the
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
CAM ozone is bit-identical on all columns; COSZEN differs by 1 ulp (the WRF binary's vectorised
sin/cos, libmvec). Tables: `~/wrf_real/rrtmg_tables` (rrtmg_lw/sw_tables.bin dumped from WRF's
rrtmg_lw/swinit with the _DBL data, ozone*.formatted, CAMtr_volume_mixing_ratio; config key
phys.rrtmg_tables). The columns run on the host with OpenMP for now (GPU port pending).
Clouds: section 9. Not yet: aer_opt, swint_opt.
Acceptance (wrf_suite.py): real_p3_r8_fp64_cpu steps 1-3 U 2.2e-14, THM 1.8e-16; P3 1 h CPU FP64
interior U 7e-12; real_p3_fp32_1h U 1.9e-4, THM 7.3e-7.


## 9. Clouds in the radiation (P4 = P2 + WSM6 + RRTMG + diff_opt 1 / km_opt 4)

Reference: WRF v4.6.1 real*8 run ~/wrf_real/run_P4_r8 (kf + WSM6 + RRTMG radt = 1,
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
- WSM6 on real data: `CPPWRF_WSM6_DUMP` debug dumps (not in main) of the first and second WSM6
  calls of P4/P5 through WRF's Fortran WSM6 (tests/phys/fortran, WRF and scalar flags): C++ =
  Fortran bit for bit, including rain/ice/snow sedimentation. RE_ICE/RE_SNOW after step 1 equal
  WRF's to 2e-12 / exactly.
- P4 step 1 (wrf_suite real_p4_r8_fp64_cpu): every physics tendency (RTHCUTEN, RQ*CUTEN, R*BLTEN,
  RTHRATEN, GLW) and U at roundoff. Two columns differ in QRAIN/QCLOUD (1e-7): a debug-WRF dump of
  WRF's WSM6 inputs equals CPPWRF's to roundoff, but at the cloud top of KF column (j=76, i=12) the
  advected qr is a cancellation residue, +3.3e-25 in WRF and -1.1e-25 in CPPWRF, and WSM6 turns
  that sign into a 12% QRAIN difference one level below (WRF's sign in CPPWRF's input reproduces
  WRF's QRAIN bit for bit). Step 2 radiation: of 5341 cloudy columns only the 157 whose CLDFRA
  differs (from that roundoff) differ in GLW (McICA subcolumns flip with cldf); the others agree
  to < 1e-9.
- P4 1 h CUDA FP64 vs WRF (real_p4_r8_fp64_cuda_1h): U 2.7e-4, THM 3.7e-6, QV 2.5e-4 (rel L2), the
  same size as P2's.
- Restart (P4, CPU FP64): step 3 straight vs restarted after step 1 bit-identical (CLDFRA,
  RE_* carried).
