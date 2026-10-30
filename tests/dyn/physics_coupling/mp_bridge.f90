! C binding of WRF v4.8.0's microphysics coupling after the RK3 loop, as dyn_em/solve_em.F runs it
! with mp_physics = 1 (one tile = the whole domain): moist_physics_prep_em and moist_physics_finish_em
! on the tile without the spec zone (sz = 1 with specified boundaries, in x only unless periodic_x),
! kessler between them on the same columns, as microphysics_driver calls it; then
! set_physical_bc3d(h_diabatic, 'p'). For tests/dyn/physics_coupling/test_physics_coupling.cpp.
! Arrays are WRF memory arrays (ims:ime, kms:kme, jms:jme) with a halo of H = 4 (bdyzone) around
! the staggered domain (ide = nx + 1, jde = ny + 1, kde = nz + 1); 2-D (ims:ime, jms:jme).
MODULE mp_bridge
  USE ISO_C_BINDING
  USE module_configure, ONLY: grid_config_rec_type
  USE module_state_description, ONLY: P_QV, P_QC
  USE module_model_constants, ONLY: xlv, cp, ep_2, svp1, svp2, svp3, svpt0, rhowater, p1000mb
  USE module_big_step_utilities_em, ONLY: moist_physics_prep_em, moist_physics_finish_em
  USE module_mp_kessler, ONLY: kessler
  USE module_bc, ONLY: set_physical_bc3d
  IMPLICIT NONE
  INTEGER, PARAMETER :: H = 4
  TYPE(grid_config_rec_type) :: cf

CONTAINS

  ! bc: 0 periodic, 1 specified, 2 open, 3 symmetric (in x and y alike)
  SUBROUTINE c_mp_step(nx, ny, nz, bc, dt, mp_tend_lim, t2, t1, qv, qc, qr, h_diabatic, qv_diabatic, qc_diabatic, &
                       al, alb, p, pb, ph, phb, fnm, fnp, rainnc, rainncv, &
                       th_phy, pii, rho, z, dz8w) BIND(C, NAME="mp_step_c")
    INTEGER(C_INT), VALUE :: nx, ny, nz, bc
    REAL(C_DOUBLE), VALUE :: dt, mp_tend_lim
    REAL(C_DOUBLE), DIMENSION(1-H:nx+1+H, 1:nz+1, 1-H:ny+1+H), INTENT(INOUT) :: t2, t1, qv, qc, qr, &
        h_diabatic, qv_diabatic, qc_diabatic, th_phy, pii, rho, z, dz8w
    REAL(C_DOUBLE), DIMENSION(1-H:nx+1+H, 1:nz+1, 1-H:ny+1+H), INTENT(IN) :: al, alb, p, pb, ph, phb
    REAL(C_DOUBLE), DIMENSION(1:nz+1), INTENT(IN) :: fnm, fnp
    REAL(C_DOUBLE), DIMENSION(1-H:nx+1+H, 1-H:ny+1+H), INTENT(INOUT) :: rainnc, rainncv
    REAL, DIMENSION(1-H:nx+1+H, 1:nz+1, 1-H:ny+1+H) :: p8w, pf, z_at_w
    INTEGER :: ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, sz

    ids = 1; ide = nx + 1; jds = 1; jde = ny + 1; kds = 1; kde = nz + 1
    ims = ids - H; ime = ide + H; jms = jds - H; jme = jde + H; kms = kds; kme = kde
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
    cf%use_theta_m = 1
    cf%no_mp_heating = 0
    cf%mp_tend_lim = mp_tend_lim
    P_QV = 2; P_QC = 3  ! moist(:,:,:,P_QV) present (set_scalar_indices_from_config)
    p8w = 0.; pf = 0.; z_at_w = 0.

    sz = 0
    IF (cf%specified) sz = 1  ! grid%spec_zone
    ! solve_em scalar_tile_loop_1a / 2
    IF (cf%periodic_x) THEN
      its = ids; ite = ide - 1
    ELSE
      its = ids + sz; ite = ide - 1 - sz
    END IF
    jts = jds + sz; jte = jde - 1 - sz

    CALL moist_physics_prep_em(t2, t1, 300., rho, al, alb, p, p8w, p1000mb, pb, ph, phb, th_phy, pii, pf, &
                               z, z_at_w, dz8w, dt, h_diabatic, qv, qv_diabatic, qc, qc_diabatic, cf, fnm, fnp, &
                               ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                               its, ite, jts, jte, kds, kde)
    CALL kessler(th_phy, qv, qc, qr, rho, pii, dt, z, xlv, cp, ep_2, svp1, svp2, svp3, svpt0, rhowater, dz8w, &
                 rainnc, rainncv, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                 its, ite, jts, jte, kds, kde)
    CALL moist_physics_finish_em(t2, t1, 300., rho, th_phy, h_diabatic, dt, qv, qv_diabatic, qc, qc_diabatic, &
                                 pf, cf, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kds, kde)
    CALL set_physical_bc3d(h_diabatic, 'p', cf, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                           ids, ide, jds, jde, kds, kde, ids, ide, jds, jde, kds, kde)
    ! solve_em tile_bc_loop_2 (before spec_bdy_final): set_phys_bc_dry_2 on t_2, then the moist species
    CALL set_physical_bc3d(t2, 't', cf, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                           ids, ide, jds, jde, kds, kde, ids, ide, jds, jde, kds, kde)
    CALL set_physical_bc3d(qv, 'p', cf, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                           ids, ide, jds, jde, kds, kde, ids, ide, jds, jde, kds, kde)
    CALL set_physical_bc3d(qc, 'p', cf, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                           ids, ide, jds, jde, kds, kde, ids, ide, jds, jde, kds, kde)
    CALL set_physical_bc3d(qr, 'p', cf, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                           ids, ide, jds, jde, kds, kde, ids, ide, jds, jde, kds, kde)
  END SUBROUTINE c_mp_step

END MODULE mp_bridge
