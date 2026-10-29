! C binding of WRF v4.8.0's set_physical_bc3d / set_physical_bc2d (share/module_bc.F) for
! tests/dyn/halo/test_halo.cpp, one tile = the whole domain as in a serial run.
! Arrays are WRF memory arrays: (ims:ime, kms:kme, jms:jme) or (ims:ime, jms:jme) with a halo of
! H = 4 points (bdyzone) around the staggered domain (ide = nx + 1, jde = ny + 1, kde = nz + 1).
MODULE halo_bridge
  USE ISO_C_BINDING
  USE module_configure, ONLY: grid_config_rec_type
  USE module_bc, ONLY: set_physical_bc3d, set_physical_bc2d
  IMPLICIT NONE
  INTEGER, PARAMETER :: H = 4
  TYPE(grid_config_rec_type) :: cf

CONTAINS

  ! bc: 0 periodic, 1 specified, 2 open, 3 symmetric (in x and y alike); var: WRF's variable letter;
  ! rank 3 or 2
  SUBROUTINE c_phys_bc(nx, ny, nz, bc, var, rank, dat) BIND(C, NAME="phys_bc_c")
    INTEGER(C_INT), VALUE :: nx, ny, nz, bc, rank
    CHARACTER(KIND=C_CHAR), VALUE :: var
    REAL(C_DOUBLE), INTENT(INOUT) :: dat(*)
    INTEGER :: ids, ide, jds, jde, kds, kde, ims, ime, jms, jme
    ids = 1; ide = nx + 1; jds = 1; jde = ny + 1; kds = 1; kde = nz + 1
    ims = ids - H; ime = ide + H; jms = jds - H; jme = jde + H
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
    IF (rank == 3) THEN
      CALL bc3(dat)
    ELSE
      CALL bc2(dat)
    END IF
  CONTAINS
    SUBROUTINE bc3(d)
      REAL, INTENT(INOUT) :: d(ims:ime, kds:kde, jms:jme)
      CALL set_physical_bc3d(d, var, cf, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kds, kde, &
                             ids, ide, jds, jde, kds, kde, ids, ide, jds, jde, kds, kde)
    END SUBROUTINE bc3
    SUBROUTINE bc2(d)
      REAL, INTENT(INOUT) :: d(ims:ime, jms:jme)
      CALL set_physical_bc2d(d, var, cf, ids, ide, jds, jde, ims, ime, jms, jme, &
                             ids, ide, jds, jde, ids, ide, jds, jde)
    END SUBROUTINE bc2
  END SUBROUTINE c_phys_bc

END MODULE halo_bridge
