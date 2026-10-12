# Physics interface (how src/phys/ schemes plug into the dycore)

(dycore side: src/dyn/physics_coupling.{hpp,cpp}). Physics ports
target this interface. It follows WRF solve_em exactly.

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
  qi, qs, qg).

## Tendency-type physics (PBL, radiation, cumulus)
Not wired yet. As in WRF, they will run in RK stage 1 and add mass-coupled tendencies
(`RTHBLTEN`-style, multiplied by `c1h*mut + c2h`) into `Tendencies::*_tendf` and `q_tendf`,
which are held fixed for stages 2 and 3. Their kernels should return **uncoupled** tendencies
(K/s, kg/kg/s, m/s^2) in their own arrays; the coupling does the mass multiplication.
