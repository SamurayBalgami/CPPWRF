! C binding of one RK stage of WRF v4.8.0's scalar transport for one moist species, called as
! dyn_em/solve_em.F calls it (one tile = the whole domain), for tests/dyn/scalars/test_scalars.cpp:
!   last stage with positive-definite advection: rk_update_scalar_pd, then set_physical_bc3d on moist_old
!   set_physical_bc3d on moist (the halo the previous step left), rk_scalar_tend, rk_update_scalar
! Diffusion, IEVA, cumulus forcing and the lateral boundary forcing are off (tested elsewhere).
! Arrays are WRF memory arrays: (ims:ime, kms:kme, jms:jme) with a halo of H = 4 points (bdyzone) around the
! staggered domain (ide = nx + 1, jde = ny + 1, kde = nz + 1), 2-D fields (ims:ime, jms:jme).
MODULE scalar_bridge
  USE ISO_C_BINDING
  USE module_configure, ONLY: grid_config_rec_type
  USE module_bc, ONLY: set_physical_bc3d
  USE module_em, ONLY: rk_scalar_tend, rk_update_scalar, rk_update_scalar_pd
  IMPLICIT NONE
  INTEGER, PARAMETER :: H = 4, SZ = 1  ! H = bdyzone: the memory halo of a serial WRF
  INTEGER :: ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte
  TYPE(grid_config_rec_type) :: cf

CONTAINS

  ! bc: 0 periodic, 1 specified, 2 open, 3 symmetric (in x and y alike)
  SUBROUTINE setup(nx, ny, nz, bc, h_order, v_order)
    INTEGER, INTENT(IN) :: nx, ny, nz, bc, h_order, v_order
    ids = 1; ide = nx + 1; jds = 1; jde = ny + 1; kds = 1; kde = nz + 1
    ims = ids - H; ime = ide + H; jms = jds - H; jme = jde + H; kms = kds; kme = kde
    its = ids; ite = ide; jts = jds; jte = jde; kts = kds; kte = kde
    cf%specified = bc == 1
    cf%nested = .FALSE.
    cf%periodic_x = bc == 0
    cf%periodic_y = bc == 0
    cf%open_xs = bc == 2
    cf%open_xe = bc == 2
    cf%open_ys = bc == 2
    cf%open_ye = bc == 2
    cf%symmetric_xs = bc == 3
    cf%symmetric_xe = bc == 3
    cf%symmetric_ys = bc == 3
    cf%symmetric_ye = bc == 3
    cf%polar = .FALSE.
    cf%h_sca_adv_order = h_order
    cf%v_sca_adv_order = v_order
    cf%rk_ord = 3
    cf%diff_opt = 0
    cf%cu_physics = 0
    cf%zadvect_implicit = 0
    cf%spec_bdy_width = 5
  END SUBROUTINE setup

  SUBROUTINE c_scalar_stage(nx, ny, nz, bc, h_order, v_order, adv_opt, rk_step, dt_rk, &
                            q1, q2, q_tend, ru, rv, ww, muts, mub, mu1, mu2, &
                            c1h, c2h, c1f, c2f, fnm, fnp, rdnw, &
                            msfux, msfuy, msfvx, msfvy, msftx, msfty, rdx, rdy) BIND(C, NAME="scalar_stage_c")
    INTEGER(C_INT), VALUE :: nx, ny, nz, bc, h_order, v_order, adv_opt, rk_step
    REAL(C_DOUBLE), VALUE :: dt_rk, rdx, rdy
    REAL(C_DOUBLE), DIMENSION(1-H:nx+1+H, 1:nz+1, 1-H:ny+1+H), INTENT(INOUT) :: q1, q2, q_tend
    REAL(C_DOUBLE), DIMENSION(1-H:nx+1+H, 1:nz+1, 1-H:ny+1+H), INTENT(IN) :: ru, rv, ww
    REAL(C_DOUBLE), DIMENSION(1-H:nx+1+H, 1-H:ny+1+H), INTENT(IN) :: muts, mub, mu1, mu2, &
                                                                  msfux, msfuy, msfvx, msfvy, msftx, msfty
    REAL(C_DOUBLE), DIMENSION(1:nz+1), INTENT(IN) :: c1h, c2h, c1f, c2f, fnm, fnp, rdnw
    REAL, DIMENSION(1-H:nx+1+H, 1:nz+1, 1-H:ny+1+H) :: advect_tend, h_tendency, z_tendency, rqvften, &
                                                       zero3, wwE, wwI
    REAL, DIMENSION(1-H:nx+1+H, 1-H:ny+1+H) :: msfvx_inv
    REAL, DIMENSION(1:nz+1) :: zero1

    CALL setup(nx, ny, nz, bc, h_order, v_order)
    zero3 = 0.; zero1 = 0.; advect_tend = 0.
    msfvx_inv = 1. / msfvx

    IF (adv_opt == 1 .AND. rk_step == cf%rk_ord) THEN  ! solve_em: physics tendency first, then the bc
      CALL rk_update_scalar_pd(1, 1, q1(ims,kms,jms), q_tend(ims,kms,jms), c1h, c2h, mu1, mu1, mub, &
                               rk_step, dt_rk, SZ, cf, &
                               ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                               its, ite, jts, jte, kts, kte)
      CALL set_physical_bc3d(q1, 'p', cf, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                             ids, ide, jds, jde, kds, kde, its, ite, jts, jte, kts, kte)
    END IF
    CALL set_physical_bc3d(q2, 'p', cf, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                           ids, ide, jds, jde, kds, kde, its, ite, jts, jte, kts, kte)

    CALL rk_scalar_tend(1, 1, cf, .FALSE., rk_step, dt_rk, ru, rv, ww, wwE, wwI, zero3, zero3, &
                        muts, mub, mu1, c1h, c2h, c1f, c2f, zero3, &
                        q1(ims,kms,jms), q2(ims,kms,jms), q_tend(ims,kms,jms), &
                        advect_tend, h_tendency, z_tendency, rqvften, zero1, .TRUE., fnm, fnp, &
                        msfux, msfuy, msfvx, msfvx_inv, msfvy, msftx, msfty, &
                        rdx, rdy, zero1, rdnw, 0., 0., zero3, 0, 0., adv_opt, zero3, zero3, .TRUE., .TRUE., &
                        ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                        its, ite, jts, jte, kts, kte)

    CALL rk_update_scalar(scs=1, sce=1, scalar_1=q1(ims,kms,jms), scalar_2=q2(ims,kms,jms), &
                          sc_tend=q_tend(ims,kms,jms), advect_tend=advect_tend, &
                          h_tendency=h_tendency, z_tendency=z_tendency, msftx=msftx, msfty=msfty, &
                          c1=c1h, c2=c2h, mu_old=mu1, mu_new=mu2, mu_base=mub, &
                          rk_step=rk_step, dt=dt_rk, spec_zone=SZ, config_flags=cf, tenddec=.FALSE., &
                          ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
                          ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
                          its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
  END SUBROUTINE c_scalar_stage

END MODULE scalar_bridge
