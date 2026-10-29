! C bindings of WRF v4.8.0's lateral boundary routines, called as dyn_em/solve_em.F calls them
! (specified = .true., or nested = .true. after lbc_set_nested_c(1); one tile = the whole domain), for
! tests/dyn/lateral_bc/test_lateral_bc.cpp.
! Arrays are WRF memory arrays: (ims:ime, kms:kme, jms:jme) with a halo of H points around the
! staggered domain (ide = nx + 1, jde = ny + 1, kde = nz + 1), as WRF allocates them; 2-D fields
! (ims:ime, jms:jme). Boundary data of a variable come packed as xs, xe, ys, ye, then the
! tendencies btxs, btxe, btys, btye: (jms:jme, kds:kde, W) in x, (ims:ime, kds:kde, W) in y
! (k extent 1 for mu).
MODULE lbc_bridge
  USE ISO_C_BINDING
  USE module_configure, ONLY: grid_config_rec_type
  USE module_bc, ONLY: spec_bdyupdate, zero_grad_bdy, flow_dep_bdy, spec_bdy_final
  USE module_bc_em, ONLY: lbc_fcx_gcx, relax_bdy_dry, spec_bdy_dry, relax_bdy_scalar, spec_bdy_scalar, &
                          spec_bdyupdate_ph
  IMPLICIT NONE
  INTEGER, PARAMETER :: H = 3, SBW = 5, SZ = 1, RZ = 4
  INTEGER :: ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte
  TYPE(grid_config_rec_type) :: cf
  LOGICAL :: nested = .FALSE.  ! a nest (specified = .false., nested = .true.) instead of a real-data domain

CONTAINS

  SUBROUTINE c_set_nested(n) BIND(C, NAME="lbc_set_nested_c")
    INTEGER(C_INT), VALUE :: n
    nested = n /= 0
  END SUBROUTINE c_set_nested

  SUBROUTINE setup(nx, ny, nz)
    INTEGER, INTENT(IN) :: nx, ny, nz
    ids = 1; ide = nx + 1; jds = 1; jde = ny + 1; kds = 1; kde = nz + 1
    ims = ids - H; ime = ide + H; jms = jds - H; jme = jde + H; kms = kds; kme = kde
    its = ids; ite = ide; jts = jds; jte = jde; kts = kds; kte = kde
    cf%specified = .NOT. nested
    cf%nested = nested
    cf%periodic_x = .FALSE.
    cf%periodic_y = .FALSE.
    cf%symmetric_xs = .FALSE.
    cf%symmetric_xe = .FALSE.
    cf%symmetric_ys = .FALSE.
    cf%symmetric_ye = .FALSE.
    cf%open_xs = .FALSE.
    cf%open_xe = .FALSE.
    cf%open_ys = .FALSE.
    cf%open_ye = .FALSE.
    cf%polar = .FALSE.
    cf%spec_exp = 0.0
    cf%spec_bdy_width = SBW
  END SUBROUTINE setup

  ! offsets of the 8 slabs of one variable's packed boundary data (k extent nk)
  SUBROUTINE slab_offsets(nk, o)
    INTEGER, INTENT(IN) :: nk
    INTEGER, INTENT(OUT) :: o(8)
    INTEGER :: nxs, nys, s
    nxs = (jme - jms + 1) * nk * SBW
    nys = (ime - ims + 1) * nk * SBW
    o(1) = 1
    DO s = 2, 8
      IF (MOD(s - 2, 4) < 2) THEN
        o(s) = o(s - 1) + nxs
      ELSE
        o(s) = o(s - 1) + nys
      END IF
    END DO
  END SUBROUTINE slab_offsets

  SUBROUTINE c_lbc_fcx_gcx(nx, ny, nz, dt, fcx, gcx) BIND(C, NAME="lbc_fcx_gcx_c")
    INTEGER(C_INT), VALUE :: nx, ny, nz
    REAL(C_DOUBLE), VALUE :: dt
    REAL(C_DOUBLE), INTENT(OUT) :: fcx(SBW), gcx(SBW)
    CALL setup(nx, ny, nz)
    CALL lbc_fcx_gcx(fcx, gcx, SBW, SZ, RZ, dt, cf%spec_exp, cf%specified, cf%nested)
  END SUBROUTINE c_lbc_fcx_gcx

  ! solve_em 944: relax_bdy_dry into u_save, v_save, ph_save, t_save, w_save and mu_tend
  SUBROUTINE c_relax_bdy_dry(nx, ny, nz, u_save, v_save, ph_save, t_save, w_save, mu_tend, &
                             c1h, c2h, c1f, c2f, ru, rv, ph, t, w, mu, mut, &
                             ub, vb, phb, tb, wb, mub, fcx, gcx, dtbc) BIND(C, NAME="relax_bdy_dry_c")
    INTEGER(C_INT), VALUE :: nx, ny, nz
    REAL(C_DOUBLE), VALUE :: dtbc
    REAL(C_DOUBLE) :: u_save(*), v_save(*), ph_save(*), t_save(*), w_save(*), mu_tend(*)
    REAL(C_DOUBLE) :: c1h(*), c2h(*), c1f(*), c2f(*), ru(*), rv(*), ph(*), t(*), w(*), mu(*), mut(*)
    REAL(C_DOUBLE) :: ub(*), vb(*), phb(*), tb(*), wb(*), mub(*), fcx(SBW), gcx(SBW)
    INTEGER :: o(8), om(8)
    CALL setup(nx, ny, nz)
    CALL slab_offsets(kde - kds + 1, o)
    CALL slab_offsets(1, om)
    CALL relax_bdy_dry(cf, u_save, v_save, ph_save, t_save, w_save, mu_tend, c1h, c2h, c1f, c2f, &
                       ru, rv, ph, t, w, mu, mut, &
                       ub(o(1)), ub(o(2)), ub(o(3)), ub(o(4)), vb(o(1)), vb(o(2)), vb(o(3)), vb(o(4)), &
                       phb(o(1)), phb(o(2)), phb(o(3)), phb(o(4)), tb(o(1)), tb(o(2)), tb(o(3)), tb(o(4)), &
                       wb(o(1)), wb(o(2)), wb(o(3)), wb(o(4)), mub(om(1)), mub(om(2)), mub(om(3)), mub(om(4)), &
                       ub(o(5)), ub(o(6)), ub(o(7)), ub(o(8)), vb(o(5)), vb(o(6)), vb(o(7)), vb(o(8)), &
                       phb(o(5)), phb(o(6)), phb(o(7)), phb(o(8)), tb(o(5)), tb(o(6)), tb(o(7)), tb(o(8)), &
                       wb(o(5)), wb(o(6)), wb(o(7)), wb(o(8)), mub(om(5)), mub(om(6)), mub(om(7)), mub(om(8)), &
                       SBW, SZ, RZ, dtbc, fcx, gcx, &
                       ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                       ids, ide, jds, jde, kds, kde, its, ite, jts, jte, kts, kte)
  END SUBROUTINE c_relax_bdy_dry

  ! solve_em 988: spec_bdy_dry into ru_tend, rv_tend, ph_tend, t_tend, rw_tend, mu_tend
  SUBROUTINE c_spec_bdy_dry(nx, ny, nz, ru_tend, rv_tend, ph_tend, t_tend, rw_tend, mu_tend, &
                            ub, vb, phb, tb, wb, mub) BIND(C, NAME="spec_bdy_dry_c")
    INTEGER(C_INT), VALUE :: nx, ny, nz
    REAL(C_DOUBLE) :: ru_tend(*), rv_tend(*), ph_tend(*), t_tend(*), rw_tend(*), mu_tend(*)
    REAL(C_DOUBLE) :: ub(*), vb(*), phb(*), tb(*), wb(*), mub(*)
    INTEGER :: o(8), om(8)
    CALL setup(nx, ny, nz)
    CALL slab_offsets(kde - kds + 1, o)
    CALL slab_offsets(1, om)
    CALL spec_bdy_dry(cf, ru_tend, rv_tend, ph_tend, t_tend, rw_tend, mu_tend, &
                      ub(o(1)), ub(o(2)), ub(o(3)), ub(o(4)), vb(o(1)), vb(o(2)), vb(o(3)), vb(o(4)), &
                      phb(o(1)), phb(o(2)), phb(o(3)), phb(o(4)), tb(o(1)), tb(o(2)), tb(o(3)), tb(o(4)), &
                      wb(o(1)), wb(o(2)), wb(o(3)), wb(o(4)), mub(om(1)), mub(om(2)), mub(om(3)), mub(om(4)), &
                      ub(o(5)), ub(o(6)), ub(o(7)), ub(o(8)), vb(o(5)), vb(o(6)), vb(o(7)), vb(o(8)), &
                      phb(o(5)), phb(o(6)), phb(o(7)), phb(o(8)), tb(o(5)), tb(o(6)), tb(o(7)), tb(o(8)), &
                      wb(o(5)), wb(o(6)), wb(o(7)), wb(o(8)), mub(om(5)), mub(om(6)), mub(om(7)), mub(om(8)), &
                      SBW, SZ, &
                      ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                      ids, ide, jds, jde, kds, kde, its, ite, jts, jte, kts, kte)
  END SUBROUTINE c_spec_bdy_dry

  ! solve_em 1350-1486, 1613: spec_bdyupdate (two_d: 2-D field, k extent 1, as the mu calls)
  SUBROUTINE c_spec_bdyupdate(nx, ny, nz, field, tend, dts, var, two_d) BIND(C, NAME="spec_bdyupdate_c")
    INTEGER(C_INT), VALUE :: nx, ny, nz, two_d
    REAL(C_DOUBLE), VALUE :: dts
    CHARACTER(C_CHAR), VALUE :: var
    REAL(C_DOUBLE) :: field(*), tend(*)
    CALL setup(nx, ny, nz)
    IF (two_d /= 0) THEN
      CALL spec_bdyupdate(field, tend, dts, var, cf, SZ, ids, ide, jds, jde, 1, 1, ims, ime, jms, jme, 1, 1, &
                          ids, ide, jds, jde, 1, 1, its, ite, jts, jte, 1, 1)
    ELSE
      CALL spec_bdyupdate(field, tend, dts, var, cf, SZ, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                          ids, ide, jds, jde, kds, kde, its, ite, jts, jte, kts, kte)
    END IF
  END SUBROUTINE c_spec_bdyupdate

  ! solve_em 1591, then 1603 (zero_grad_bdy on w) or, in a nest, 1613 (spec_bdyupdate of w with rw_tend)
  SUBROUTINE c_update_ph_w(nx, ny, nz, ph_save, ph, ph_tend, mu_tend, muts, c1f, c2f, dts, w, rw_tend) &
      BIND(C, NAME="update_ph_w_c")
    INTEGER(C_INT), VALUE :: nx, ny, nz
    REAL(C_DOUBLE), VALUE :: dts
    REAL(C_DOUBLE) :: ph_save(*), ph(*), ph_tend(*), mu_tend(*), muts(*), c1f(*), c2f(*), w(*), rw_tend(*)
    CALL setup(nx, ny, nz)
    CALL spec_bdyupdate_ph(ph_save, ph, ph_tend, mu_tend, muts, c1f, c2f, dts, 'h', cf, SZ, &
                           ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                           ids, ide, jds, jde, kds, kde, its, ite, jts, jte, kts, kte)
    IF (nested) THEN
      CALL spec_bdyupdate(w, rw_tend, dts, 'h', cf, SZ, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                          ids, ide, jds, jde, kds, kde, its, ite, jts, jte, kts, kte)
    ELSE
      CALL zero_grad_bdy(w, 'w', cf, SZ, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                         ids, ide, jds, jde, kds, kde, its, ite, jts, jte, kts, kte)
    END IF
  END SUBROUTINE c_update_ph_w

  ! solve_em 2314 and 2331: relax_bdy_scalar, then spec_bdy_scalar into moist_tend
  SUBROUTINE c_relax_spec_scalar(nx, ny, nz, q_tend, q, mut, c1h, c2h, qb, fcx, gcx, dtbc, relax, spec) &
      BIND(C, NAME="relax_spec_scalar_c")
    INTEGER(C_INT), VALUE :: nx, ny, nz, relax, spec
    REAL(C_DOUBLE), VALUE :: dtbc
    REAL(C_DOUBLE) :: q_tend(*), q(*), mut(*), c1h(*), c2h(*), qb(*), fcx(SBW), gcx(SBW)
    INTEGER :: o(8)
    CALL setup(nx, ny, nz)
    CALL slab_offsets(kde - kds + 1, o)
    IF (relax /= 0) &
      CALL relax_bdy_scalar(q_tend, q, mut, c1h, c2h, qb(o(1)), qb(o(2)), qb(o(3)), qb(o(4)), &
                            qb(o(5)), qb(o(6)), qb(o(7)), qb(o(8)), SBW, SZ, RZ, dtbc, fcx, gcx, cf, &
                            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                            ids, ide, jds, jde, kds, kde, its, ite, jts, jte, kts, kte)
    IF (spec /= 0) &
      CALL spec_bdy_scalar(q_tend, qb(o(1)), qb(o(2)), qb(o(3)), qb(o(4)), qb(o(5)), qb(o(6)), qb(o(7)), qb(o(8)), &
                           SBW, SZ, cf, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                           ids, ide, jds, jde, kds, kde, its, ite, jts, jte, kts, kte)
  END SUBROUTINE c_relax_spec_scalar

  ! solve_em 2500: flow_dep_bdy of a moist species without boundary data
  SUBROUTINE c_flow_dep_bdy(nx, ny, nz, q, ru_m, rv_m) BIND(C, NAME="flow_dep_bdy_c")
    INTEGER(C_INT), VALUE :: nx, ny, nz
    REAL(C_DOUBLE) :: q(*), ru_m(*), rv_m(*)
    CALL setup(nx, ny, nz)
    CALL flow_dep_bdy(q, ru_m, rv_m, cf, SZ, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                      ids, ide, jds, jde, kds, kde, its, ite, jts, jte, kts, kte)
  END SUBROUTINE c_flow_dep_bdy

  ! solve_em 4711-4795: spec_bdy_final (two_d: mu, k extent 1)
  SUBROUTINE c_spec_bdy_final(nx, ny, nz, field, mu, c1, c2, msf, fb, var, dtbc, two_d) &
      BIND(C, NAME="spec_bdy_final_c")
    INTEGER(C_INT), VALUE :: nx, ny, nz, two_d
    REAL(C_DOUBLE), VALUE :: dtbc
    CHARACTER(C_CHAR), VALUE :: var
    REAL(C_DOUBLE) :: field(*), mu(*), c1(*), c2(*), msf(*), fb(*)
    INTEGER :: o(8)
    CALL setup(nx, ny, nz)
    IF (two_d /= 0) THEN
      CALL slab_offsets(1, o)
      CALL spec_bdy_final(field, mu, c1, c2, msf, fb(o(1)), fb(o(2)), fb(o(3)), fb(o(4)), &
                          fb(o(5)), fb(o(6)), fb(o(7)), fb(o(8)), var, cf, SBW, SZ, dtbc, &
                          ids, ide, jds, jde, 1, 1, ims, ime, jms, jme, 1, 1, &
                          ids, ide, jds, jde, 1, 1, its, ite, jts, jte, 1, 1)
    ELSE
      CALL slab_offsets(kde - kds + 1, o)
      CALL spec_bdy_final(field, mu, c1, c2, msf, fb(o(1)), fb(o(2)), fb(o(3)), fb(o(4)), &
                          fb(o(5)), fb(o(6)), fb(o(7)), fb(o(8)), var, cf, SBW, SZ, dtbc, &
                          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                          ids, ide, jds, jde, kds, kde, its, ite, jts, jte, kts, kte)
    END IF
  END SUBROUTINE c_spec_bdy_final

END MODULE lbc_bridge
