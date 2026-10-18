# Dycore checkpoint dumps: CPPWRF vs a debug WRF, statement by statement

`solve_em_dumps.patch` (against WRF 4.8.0 `dyn_em/solve_em.F`, r8 build, serial) writes, for one
time step, the arrays at seven checkpoints of every RK stage:

| file | where in solve_em | arrays |
|---|---|---|
| vm_tend_rk*N* | before small_step_prep | ru/rv/rw/ph/t/mu_tend, ru, rv, php, alt, ww, muu, mut, mu_2, mub, u_2, v_2 |
| vm_s0_rk*N* | after the small-step set_physical_bc calls | p, al, t_save, u_save, v_save, ph_2, php, alt, cqu, pb, muu, mu_2 |
| vm_s1_rk*N* | after advance_uv + spec_bdyupdate (1st acoustic step) | u_2, v_2 |
| vm_s2_rk*N* | before advance_w | t_2, ww, mu_2, muts |
| vm_s3_rk*N* | after sumflux/spec_bdyupdate_ph/calc_p_rho | w_2, ph_2, p, al |
| vm_state_rk*N* | before calc_p_rho_phi | u_2, v_2, w_2, t_2, ph_2, mu_2, muts, qv, ww, ru_m, rv_m, ww_m |
| vm_p_rk*N* | after calc_p_rho_phi | p, al |

(The patch also carries older moist/WSM6/radiation debug writes.) The step is set in the patch
(`grid%itimestep == 3`). CPPWRF writes the same checkpoints with

    CPPWRF_VM_DUMP=dir CPPWRF_VM_DUMP_STEP=3 build/cppwrf configs/real_r1.ini --set run.run_time=180

and `python3 vmcmp.py <wrf dir> <cppwrf dir>` prints, per checkpoint and array, the largest
difference and the fraction of points that differ (WRF files are big-endian with memory halos).

Found this way (2026-10-15), all now as WRF:
- start_em's p_surf uses `(...)**0.5`, i.e. pow, not sqrt (1 ulp in a few columns);
- glibc libmvec in vectorised loops: calc_p_rho_phi (vpow, LOG), with the scalar libm for a trailing
  odd iteration and in the interior-only call after spec_bdy_final (`cppwrf::vm`, core/backend.hpp);
- summation orders: calc_mu_uv `0.5*(MU(i)+MU(i-1)+MUB(i)+MUB(i-1))` (edges double the column),
  calc_ww_cp's own `0.5*(MUP(i)+MUB(i)+MUP(i-1)+MUB(i-1))` with the halo column, advance_uv
  `((ph(k+1)-ph_m(k+1))+(ph(k)-ph_m(k)))` and dx = 1./rdx, small-step al = `-1./m * (...)`,
  spec_bdyupdate_ph `field*m_old/m_new`;
- small_step_prep at rk_step 1 copies MUUS = MUU, MUVS = MUV;
- advect_scalar_pd has its own flux5/flux3 statement functions (folded constants 37./60. ...) and,
  at i = ids+1 only, cr = (ru/mu)*dt/dx;
- halos: no refresh after spec_bdy_final, and mu_2's halo is not refreshed before rk_step_prep.
