# Physics interface (how src/phys/ schemes plug into the dycore)

The dycore side is src/dyn/physics_coupling.{hpp,cpp} (microphysics) and src/dyn/column_physics.hpp
(the interface of the column physics driver, src/phys/driver). Both follow WRF solve_em.

## Call order per time step
1. RK3 loop (3 stages): large step, acoustic small steps, scalar transport, calc_p_rho_phi.
2. `PhysicsCoupling::step(s, base, dt)` (once per step, full dt), as WRF does after the RK loop:
   1. **moist_physics_prep_em** (`KPrep`): computes, at mass points,
      `th_phy = (t_2 + t0) / (1 + Rv/Rd qv)` (dry theta, K), `rho = 1/(al + alb)` (dry density),
      `pii = ((p + pb)/p0)^(Rd/cp)`, `z` (mass-level height, m), `dz8w` (layer thickness, m),
      and saves `h_diabatic = th_phy` and `qv_diabatic = qv` (values before microphysics).
   2. **microphysics** (`column_for`, one thread per column): updates `th_phy` and the moist
      species `s.q_2[P_QV], [P_QC], [P_QR], ...` in place; accumulates `rainnc`, sets `rainncv`.
   3. **moist_physics_finish_em** (`KFinish`, use_theta_m = 1): new `t_2` (theta_m - t0) from
      the theta and qv changes; `h_diabatic` = the heating rate (K/s, theta_m units), with
      WRF's `mp_tend_lim` limit.
3. calc_p_rho_phi again (p and al from the new theta_m and moisture).

During the next step, `h_diabatic` is added to the theta tendency in every RK stage
(rk_addtend_dry: `t_tend += (c1h mut + c2h) h_diabatic`) and removed in the last
small_step_finish, so the acoustic/gravity-wave dynamics feel the latent heating but theta is
not changed twice.

## What a microphysics scheme must provide
- A header-only kernel struct (as `phys::kessler::KesslerKernel`): trivially copyable,
  `CPPWRF_HD void operator()(int i, int j) const`, Views as members, no allocation.
- Arguments with WRF names and units: `th` (dry theta), moist species (mixing ratios, kg/kg,
  relative to dry air), `rho`, `pii`, `p` (Pa, if needed: pass `p + pb`), `dz8w`, `z`, `dt`,
  `rainnc`, `rainncv` (mm), plus scheme-specific fields.
- To add a scheme: a new `mp_physics` value in `PhysicsParams`, the launch in
  `PhysicsCoupling::step`, and the species it needs in `moist.species` (e.g. WSM6: qv, qc, qr,
  qi, qs, qg). A new scheme is wired in only after its column differential test against WRF passes.

## Tendency-type physics (PBL, radiation, cumulus)
These run in RK stage 1 through `dyn::ColumnPhysics` (`part1` = first_rk_step_part1, `part2` =
first_rk_step_part2, `after_rk` = phy_prep_part2), implemented by `phys::driver::PhysicsSuite`. As in
WRF they produce uncoupled tendencies (RTHBLTEN, RUCUTEN, RTHRATEN, ...), which calculate_phy_tend
multiplies by the column mass and update_phy_ten adds to `Tendencies::*_tendf` and `q_tendf`; those are
held fixed for RK stages 2 and 3. The full call map is in docs/PHYSICS_SUITE.md.
