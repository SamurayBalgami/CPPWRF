MODULE bdy_fortran_bridge
  USE ISO_C_BINDING
  USE module_configure, ONLY: grid_config_rec_type
  USE module_bc, ONLY: relax_bdytend_core, spec_bdytend, spec_bdyupdate, zero_grad_bdy, &
                       flow_dep_bdy, spec_bdy_final, set_physical_bc3d, set_physical_bc2d
  USE module_bc_em, ONLY: lbc_fcx_gcx, spec_bdyupdate_ph, mass_weight
  IMPLICIT NONE

CONTAINS

  SUBROUTINE init_config_flags(config_flags)
    TYPE(grid_config_rec_type), INTENT(OUT) :: config_flags
    config_flags%specified = .TRUE.
    config_flags%nested = .FALSE.
    config_flags%periodic_x = .FALSE.
    config_flags%periodic_y = .FALSE.
    config_flags%symmetric_xs = .FALSE.
    config_flags%symmetric_xe = .FALSE.
    config_flags%symmetric_ys = .FALSE.
    config_flags%symmetric_ye = .FALSE.
    config_flags%open_xs = .FALSE.
    config_flags%open_xe = .FALSE.
    config_flags%open_ys = .FALSE.
    config_flags%open_ye = .FALSE.
    config_flags%polar = .FALSE.
    config_flags%spec_exp = 0.0_8
  END SUBROUTINE init_config_flags

  SUBROUTINE bridge_lbc_fcx_gcx(fcx, gcx, spec_bdy_width, spec_zone, relax_zone, dt, spec_exp) BIND(C, NAME="bridge_lbc_fcx_gcx")
    INTEGER(C_INT), VALUE :: spec_bdy_width, spec_zone, relax_zone
    REAL(C_DOUBLE), VALUE :: dt, spec_exp
    REAL(C_DOUBLE), INTENT(OUT) :: fcx(spec_bdy_width), gcx(spec_bdy_width)
    CALL lbc_fcx_gcx(fcx, gcx, spec_bdy_width, spec_zone, relax_zone, dt, spec_exp, .TRUE., .FALSE.)
  END SUBROUTINE bridge_lbc_fcx_gcx

  SUBROUTINE bridge_relax_bdytend_core_3d(field, field_tend, &
       bxs, bxe, bys, bye, btxs, btxe, btys, btye, fcx, gcx, dtbc, &
       var_code, nx, ny, nz, spec_bdy_width, spec_zone, relax_zone) BIND(C, NAME="bridge_relax_bdytend_core_3d")
    INTEGER(C_INT), VALUE :: nx, ny, nz, spec_bdy_width, spec_zone, relax_zone
    CHARACTER(C_CHAR), VALUE :: var_code
    REAL(C_DOUBLE), VALUE :: dtbc
    REAL(C_DOUBLE), INTENT(IN) :: fcx(spec_bdy_width), gcx(spec_bdy_width)
    REAL(C_DOUBLE), INTENT(IN) :: field(*), bxs(*), bxe(*), bys(*), bye(*)
    REAL(C_DOUBLE), INTENT(IN) :: btxs(*), btxe(*), btys(*), btye(*)
    REAL(C_DOUBLE), INTENT(INOUT) :: field_tend(*)

    INTEGER :: ids, ide, jds, jde, kds, kde
    INTEGER :: ims, ime, jms, jme, kms, kme
    INTEGER :: its, ite, jts, jte, kts, kte
    INTEGER :: ni, nj, nk

    TYPE(grid_config_rec_type) :: config_flags
    CALL init_config_flags(config_flags)

    ids = 1; ide = nx + 1
    jds = 1; jde = ny + 1
    kds = 1; kde = nz + 1

    ni = nx; IF (var_code == 'u' .OR. var_code == 'U') ni = nx + 1
    nj = ny; IF (var_code == 'v' .OR. var_code == 'V') nj = ny + 1
    nk = nz; IF (var_code == 'h' .OR. var_code == 'H' .OR. var_code == 'w' .OR. var_code == 'W') nk = nz + 1

    ims = 1; ime = ni
    jms = 1; jme = nj
    kms = 1; kme = nk

    its = 1; ite = ni
    jts = 1; jte = nj
    kts = 1; kte = nk

    CALL relax_bdytend_core(field, field_tend, &
         bxs, bxe, bys, bye, btxs, btxe, btys, btye, &
         var_code, config_flags, spec_bdy_width, spec_zone, relax_zone, &
         dtbc, fcx, gcx, &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte, &
         its, ite, jts, jte, kts, kte, &
         ims, ime, jms, jme, kms, kme)
  END SUBROUTINE bridge_relax_bdytend_core_3d

  SUBROUTINE bridge_relax_bdytend_core_2d(field, field_tend, &
       bxs, bxe, bys, bye, btxs, btxe, btys, btye, fcx, gcx, dtbc, &
       nx, ny, spec_bdy_width, spec_zone, relax_zone) BIND(C, NAME="bridge_relax_bdytend_core_2d")
    INTEGER(C_INT), VALUE :: nx, ny, spec_bdy_width, spec_zone, relax_zone
    REAL(C_DOUBLE), VALUE :: dtbc
    REAL(C_DOUBLE), INTENT(IN) :: fcx(spec_bdy_width), gcx(spec_bdy_width)
    REAL(C_DOUBLE), INTENT(IN) :: field(*), bxs(*), bxe(*), bys(*), bye(*)
    REAL(C_DOUBLE), INTENT(IN) :: btxs(*), btxe(*), btys(*), btye(*)
    REAL(C_DOUBLE), INTENT(INOUT) :: field_tend(*)

    INTEGER :: ids, ide, jds, jde
    INTEGER :: ims, ime, jms, jme
    INTEGER :: its, ite, jts, jte

    TYPE(grid_config_rec_type) :: config_flags
    CALL init_config_flags(config_flags)

    ids = 1; ide = nx + 1
    jds = 1; jde = ny + 1
    ims = 1; ime = nx
    jms = 1; jme = ny
    its = 1; ite = nx
    jts = 1; jte = ny

    CALL relax_bdytend_core(field, field_tend, &
         bxs, bxe, bys, bye, btxs, btxe, btys, btye, &
         'm', config_flags, spec_bdy_width, spec_zone, relax_zone, &
         dtbc, fcx, gcx, &
         ids, ide, jds, jde, 1, 1, &
         ims, ime, jms, jme, 1, 1, &
         its, ite, jts, jte, 1, 1, &
         its, ite, jts, jte, 1, 1, &
         ims, ime, jms, jme, 1, 1)
  END SUBROUTINE bridge_relax_bdytend_core_2d

  SUBROUTINE bridge_spec_bdytend_3d(field_tend, &
       bxs, bxe, bys, bye, btxs, btxe, btys, btye, &
       var_code, nx, ny, nz, spec_bdy_width, spec_zone) BIND(C, NAME="bridge_spec_bdytend_3d")
    INTEGER(C_INT), VALUE :: nx, ny, nz, spec_bdy_width, spec_zone
    CHARACTER(C_CHAR), VALUE :: var_code
    REAL(C_DOUBLE), INTENT(IN) :: bxs(*), bxe(*), bys(*), bye(*)
    REAL(C_DOUBLE), INTENT(IN) :: btxs(*), btxe(*), btys(*), btye(*)
    REAL(C_DOUBLE), INTENT(OUT) :: field_tend(*)

    INTEGER :: ids, ide, jds, jde, kds, kde
    INTEGER :: ims, ime, jms, jme, kms, kme
    INTEGER :: its, ite, jts, jte, kts, kte
    INTEGER :: ni, nj, nk

    TYPE(grid_config_rec_type) :: config_flags
    CALL init_config_flags(config_flags)

    ids = 1; ide = nx + 1
    jds = 1; jde = ny + 1
    kds = 1; kde = nz + 1

    ni = nx; IF (var_code == 'u' .OR. var_code == 'U') ni = nx + 1
    nj = ny; IF (var_code == 'v' .OR. var_code == 'V') nj = ny + 1
    nk = nz; IF (var_code == 'h' .OR. var_code == 'H' .OR. var_code == 'w' .OR. var_code == 'W') nk = nz + 1

    ims = 1; ime = ni
    jms = 1; jme = nj
    kms = 1; kme = nk

    its = 1; ite = ni
    jts = 1; jte = nj
    kts = 1; kte = nk

    CALL spec_bdytend(field_tend, &
         bxs, bxe, bys, bye, btxs, btxe, btys, btye, &
         var_code, config_flags, spec_bdy_width, spec_zone, &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte, &
         its, ite, jts, jte, kts, kte)
  END SUBROUTINE bridge_spec_bdytend_3d

  SUBROUTINE bridge_spec_bdytend_2d(field_tend, &
       bxs, bxe, bys, bye, btxs, btxe, btys, btye, &
       nx, ny, spec_bdy_width, spec_zone) BIND(C, NAME="bridge_spec_bdytend_2d")
    INTEGER(C_INT), VALUE :: nx, ny, spec_bdy_width, spec_zone
    REAL(C_DOUBLE), INTENT(IN) :: bxs(*), bxe(*), bys(*), bye(*)
    REAL(C_DOUBLE), INTENT(IN) :: btxs(*), btxe(*), btys(*), btye(*)
    REAL(C_DOUBLE), INTENT(OUT) :: field_tend(*)

    INTEGER :: ids, ide, jds, jde
    INTEGER :: ims, ime, jms, jme
    INTEGER :: its, ite, jts, jte

    TYPE(grid_config_rec_type) :: config_flags
    CALL init_config_flags(config_flags)

    ids = 1; ide = nx + 1
    jds = 1; jde = ny + 1
    ims = 1; ime = nx
    jms = 1; jme = ny
    its = 1; ite = nx
    jts = 1; jte = ny

    CALL spec_bdytend(field_tend, &
         bxs, bxe, bys, bye, btxs, btxe, btys, btye, &
         'm', config_flags, spec_bdy_width, spec_zone, &
         ids, ide, jds, jde, 1, 1, &
         ims, ime, jms, jme, 1, 1, &
         its, ite, jts, jte, 1, 1, &
         its, ite, jts, jte, 1, 1)
  END SUBROUTINE bridge_spec_bdytend_2d

  SUBROUTINE bridge_spec_bdyupdate_3d(field, field_tend, dt, var_code, &
       nx, ny, nz, spec_zone) BIND(C, NAME="bridge_spec_bdyupdate_3d")
    INTEGER(C_INT), VALUE :: nx, ny, nz, spec_zone
    REAL(C_DOUBLE), VALUE :: dt
    CHARACTER(C_CHAR), VALUE :: var_code
    REAL(C_DOUBLE), INTENT(INOUT) :: field(*)
    REAL(C_DOUBLE), INTENT(IN) :: field_tend(*)

    INTEGER :: ids, ide, jds, jde, kds, kde
    INTEGER :: ims, ime, jms, jme, kms, kme
    INTEGER :: its, ite, jts, jte, kts, kte
    INTEGER :: ni, nj, nk

    TYPE(grid_config_rec_type) :: config_flags
    CALL init_config_flags(config_flags)

    ids = 1; ide = nx + 1
    jds = 1; jde = ny + 1
    kds = 1; kde = nz + 1

    ni = nx; IF (var_code == 'u' .OR. var_code == 'U') ni = nx + 1
    nj = ny; IF (var_code == 'v' .OR. var_code == 'V') nj = ny + 1
    nk = nz; IF (var_code == 'h' .OR. var_code == 'H' .OR. var_code == 'w' .OR. var_code == 'W') nk = nz + 1

    ims = 1; ime = ni
    jms = 1; jme = nj
    kms = 1; kme = nk

    its = 1; ite = ni
    jts = 1; jte = nj
    kts = 1; kte = nk

    CALL spec_bdyupdate(field, field_tend, dt, var_code, config_flags, spec_zone, &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte, &
         its, ite, jts, jte, kts, kte)
  END SUBROUTINE bridge_spec_bdyupdate_3d

  SUBROUTINE bridge_spec_bdyupdate_2d(field, field_tend, dt, &
       nx, ny, spec_zone) BIND(C, NAME="bridge_spec_bdyupdate_2d")
    INTEGER(C_INT), VALUE :: nx, ny, spec_zone
    REAL(C_DOUBLE), VALUE :: dt
    REAL(C_DOUBLE), INTENT(INOUT) :: field(*)
    REAL(C_DOUBLE), INTENT(IN) :: field_tend(*)

    INTEGER :: ids, ide, jds, jde
    INTEGER :: ims, ime, jms, jme
    INTEGER :: its, ite, jts, jte

    TYPE(grid_config_rec_type) :: config_flags
    CALL init_config_flags(config_flags)

    ids = 1; ide = nx + 1
    jds = 1; jde = ny + 1
    ims = 1; ime = nx
    jms = 1; jme = ny
    its = 1; ite = nx
    jts = 1; jte = ny

    CALL spec_bdyupdate(field, field_tend, dt, 'm', config_flags, spec_zone, &
         ids, ide, jds, jde, 1, 1, &
         ims, ime, jms, jme, 1, 1, &
         its, ite, jts, jte, 1, 1, &
         its, ite, jts, jte, 1, 1)
  END SUBROUTINE bridge_spec_bdyupdate_2d

  SUBROUTINE bridge_spec_bdyupdate_ph(ph_save, field, field_tend, mu_tend, muts, &
       c1, c2, dt, nx, ny, nz, spec_zone) BIND(C, NAME="bridge_spec_bdyupdate_ph")
    INTEGER(C_INT), VALUE :: nx, ny, nz, spec_zone
    REAL(C_DOUBLE), VALUE :: dt
    REAL(C_DOUBLE), INTENT(IN) :: c1(nz + 1), c2(nz + 1)
    REAL(C_DOUBLE), INTENT(IN) :: ph_save(*)
    REAL(C_DOUBLE), INTENT(INOUT) :: field(*)
    REAL(C_DOUBLE), INTENT(IN) :: field_tend(*)
    REAL(C_DOUBLE), INTENT(IN) :: mu_tend(*)
    REAL(C_DOUBLE), INTENT(IN) :: muts(*)

    INTEGER :: ids, ide, jds, jde, kds, kde
    INTEGER :: ims, ime, jms, jme, kms, kme
    INTEGER :: its, ite, jts, jte, kts, kte

    TYPE(grid_config_rec_type) :: config_flags
    CALL init_config_flags(config_flags)

    ids = 1; ide = nx + 1
    jds = 1; jde = ny + 1
    kds = 1; kde = nz + 1
    ims = 1; ime = nx
    jms = 1; jme = ny
    kms = 1; kme = nz + 1
    its = 1; ite = nx
    jts = 1; jte = ny
    kts = 1; kte = nz + 1

    CALL spec_bdyupdate_ph(ph_save, field, field_tend, mu_tend, muts, &
         c1, c2, dt, 'h', config_flags, spec_zone, &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte, &
         its, ite, jts, jte, kts, kte)
  END SUBROUTINE bridge_spec_bdyupdate_ph

  SUBROUTINE bridge_zero_grad_bdy(field, nx, ny, nz, spec_zone) BIND(C, NAME="bridge_zero_grad_bdy")
    INTEGER(C_INT), VALUE :: nx, ny, nz, spec_zone
    REAL(C_DOUBLE), INTENT(INOUT) :: field(*)

    INTEGER :: ids, ide, jds, jde, kds, kde
    INTEGER :: ims, ime, jms, jme, kms, kme
    INTEGER :: its, ite, jts, jte, kts, kte

    TYPE(grid_config_rec_type) :: config_flags
    CALL init_config_flags(config_flags)

    ids = 1; ide = nx + 1
    jds = 1; jde = ny + 1
    kds = 1; kde = nz + 1
    ims = 1; ime = nx
    jms = 1; jme = ny
    kms = 1; kme = nz + 1
    its = 1; ite = nx
    jts = 1; jte = ny
    kts = 1; kte = nz + 1

    CALL zero_grad_bdy(field, 'w', config_flags, spec_zone, &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte, &
         its, ite, jts, jte, kts, kte)
  END SUBROUTINE bridge_zero_grad_bdy

  SUBROUTINE bridge_flow_dep_bdy(field, u, v, nx, ny, nz, spec_zone) BIND(C, NAME="bridge_flow_dep_bdy")
    INTEGER(C_INT), VALUE :: nx, ny, nz, spec_zone
    REAL(C_DOUBLE), INTENT(INOUT) :: field(*)
    REAL(C_DOUBLE), INTENT(IN) :: u(*)
    REAL(C_DOUBLE), INTENT(IN) :: v(*)

    INTEGER :: ids, ide, jds, jde, kds, kde
    INTEGER :: ims, ime, jms, jme, kms, kme
    INTEGER :: its, ite, jts, jte, kts, kte

    TYPE(grid_config_rec_type) :: config_flags
    CALL init_config_flags(config_flags)

    ids = 1; ide = nx + 1
    jds = 1; jde = ny + 1
    kds = 1; kde = nz + 1
    ims = 1; ime = nx + 1
    jms = 1; jme = ny + 1
    kms = 1; kme = nz
    its = 1; ite = nx
    jts = 1; jte = ny
    kts = 1; kte = nz

    CALL flow_dep_bdy(field, u, v, config_flags, spec_zone, &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte, &
         its, ite, jts, jte, kts, kte)
  END SUBROUTINE bridge_flow_dep_bdy

  SUBROUTINE bridge_spec_bdy_final_3d(field, mu, c1, c2, msf, &
       bxs, bxe, bys, bye, btxs, btxe, btys, btye, &
       var_code, dtbc, nx, ny, nz, spec_bdy_width, spec_zone) BIND(C, NAME="bridge_spec_bdy_final_3d")
    INTEGER(C_INT), VALUE :: nx, ny, nz, spec_bdy_width, spec_zone
    CHARACTER(C_CHAR), VALUE :: var_code
    REAL(C_DOUBLE), VALUE :: dtbc
    REAL(C_DOUBLE), INTENT(INOUT) :: field(*)
    REAL(C_DOUBLE), INTENT(IN) :: mu(*), c1(*), c2(*), msf(*)
    REAL(C_DOUBLE), INTENT(IN) :: bxs(*), bxe(*), bys(*), bye(*)
    REAL(C_DOUBLE), INTENT(IN) :: btxs(*), btxe(*), btys(*), btye(*)

    INTEGER :: ids, ide, jds, jde, kds, kde
    INTEGER :: ims, ime, jms, jme, kms, kme
    INTEGER :: its, ite, jts, jte, kts, kte
    INTEGER :: ni, nj, nk

    TYPE(grid_config_rec_type) :: config_flags
    CALL init_config_flags(config_flags)

    ids = 1; ide = nx + 1
    jds = 1; jde = ny + 1
    kds = 1; kde = nz + 1

    ni = nx; IF (var_code == 'u' .OR. var_code == 'U') ni = nx + 1
    nj = ny; IF (var_code == 'v' .OR. var_code == 'V') nj = ny + 1
    nk = nz; IF (var_code == 'h' .OR. var_code == 'H' .OR. var_code == 'w' .OR. var_code == 'W') nk = nz + 1

    ims = 1; ime = ni
    jms = 1; jme = nj
    kms = 1; kme = nk

    its = 1; ite = ni
    jts = 1; jte = nj
    kts = 1; kte = nk

    CALL spec_bdy_final(field, mu, c1, c2, msf, &
         bxs, bxe, bys, bye, btxs, btxe, btys, btye, &
         var_code, config_flags, spec_bdy_width, spec_zone, dtbc, &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte, &
         its, ite, jts, jte, kts, kte)
  END SUBROUTINE bridge_spec_bdy_final_3d

  SUBROUTINE bridge_spec_bdy_final_2d(field, &
       bxs, bxe, bys, bye, btxs, btxe, btys, btye, &
       dtbc, nx, ny, spec_bdy_width, spec_zone) BIND(C, NAME="bridge_spec_bdy_final_2d")
    INTEGER(C_INT), VALUE :: nx, ny, spec_bdy_width, spec_zone
    REAL(C_DOUBLE), VALUE :: dtbc
    REAL(C_DOUBLE), INTENT(INOUT) :: field(*)
    REAL(C_DOUBLE), INTENT(IN) :: bxs(*), bxe(*), bys(*), bye(*)
    REAL(C_DOUBLE), INTENT(IN) :: btxs(*), btxe(*), btys(*), btye(*)

    INTEGER :: ids, ide, jds, jde
    INTEGER :: ims, ime, jms, jme
    INTEGER :: its, ite, jts, jte
    REAL(C_DOUBLE) :: dummy_mu(nx, ny), dummy_msf(nx, ny)
    REAL(C_DOUBLE) :: dummy_c1(1), dummy_c2(1)

    TYPE(grid_config_rec_type) :: config_flags
    CALL init_config_flags(config_flags)
    dummy_mu = 1.0_8
    dummy_msf = 1.0_8
    dummy_c1 = 1.0_8
    dummy_c2 = 0.0_8

    ids = 1; ide = nx + 1
    jds = 1; jde = ny + 1
    ims = 1; ime = nx
    jms = 1; jme = ny
    its = 1; ite = nx
    jts = 1; jte = ny

    CALL spec_bdy_final(field, dummy_mu, dummy_c1, dummy_c2, dummy_msf, &
         bxs, bxe, bys, bye, btxs, btxe, btys, btye, &
         'm', config_flags, spec_bdy_width, spec_zone, dtbc, &
         ids, ide, jds, jde, 1, 1, &
         ims, ime, jms, jme, 1, 1, &
         its, ite, jts, jte, 1, 1, &
         its, ite, jts, jte, 1, 1)
  END SUBROUTINE bridge_spec_bdy_final_2d

  SUBROUTINE bridge_set_physical_bc3d(dat, var_code, nx, ny, nz, halo) BIND(C, NAME="bridge_set_physical_bc3d")
    INTEGER(C_INT), VALUE :: nx, ny, nz, halo
    CHARACTER(C_CHAR), VALUE :: var_code
    REAL(C_DOUBLE), INTENT(INOUT) :: dat(*)

    INTEGER :: ids, ide, jds, jde, kds, kde
    INTEGER :: ims, ime, jms, jme, kms, kme
    INTEGER :: its, ite, jts, jte, kts, kte
    INTEGER :: ni, nj, nk

    TYPE(grid_config_rec_type) :: config_flags
    CALL init_config_flags(config_flags)

    ids = 1; ide = nx + 1
    jds = 1; jde = ny + 1
    kds = 1; kde = nz + 1

    ni = nx; IF (var_code == 'u' .OR. var_code == 'U') ni = nx + 1
    nj = ny; IF (var_code == 'v' .OR. var_code == 'V') nj = ny + 1
    nk = nz; IF (var_code == 'h' .OR. var_code == 'H' .OR. var_code == 'w' .OR. var_code == 'W') nk = nz + 1

    ims = 1 - halo; ime = ni + halo
    jms = 1 - halo; jme = nj + halo
    kms = 1; kme = nk

    its = 1; ite = ide
    jts = 1; jte = jde
    kts = 1; kte = nk

    CALL set_physical_bc3d(dat, var_code, config_flags, &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte, &
         its, ite, jts, jte, kts, kte)
  END SUBROUTINE bridge_set_physical_bc3d

  SUBROUTINE bridge_set_physical_bc2d(dat, var_code, nx, ny, halo) BIND(C, NAME="bridge_set_physical_bc2d")
    INTEGER(C_INT), VALUE :: nx, ny, halo
    CHARACTER(C_CHAR), VALUE :: var_code
    REAL(C_DOUBLE), INTENT(INOUT) :: dat(*)

    INTEGER :: ids, ide, jds, jde
    INTEGER :: ims, ime, jms, jme
    INTEGER :: its, ite, jts, jte
    INTEGER :: ni, nj

    TYPE(grid_config_rec_type) :: config_flags
    CALL init_config_flags(config_flags)

    ids = 1; ide = nx + 1
    jds = 1; jde = ny + 1

    ni = nx; IF (var_code == 'u' .OR. var_code == 'U') ni = nx + 1
    nj = ny; IF (var_code == 'v' .OR. var_code == 'V') nj = ny + 1

    ims = 1 - halo; ime = ni + halo
    jms = 1 - halo; jme = nj + halo

    its = 1; ite = ide
    jts = 1; jte = jde

    CALL set_physical_bc2d(dat, var_code, config_flags, &
         ids, ide, jds, jde, &
         ims, ime, jms, jme, &
         its, ite, jts, jte, &
         its, ite, jts, jte)
  END SUBROUTINE bridge_set_physical_bc2d

END MODULE bdy_fortran_bridge
