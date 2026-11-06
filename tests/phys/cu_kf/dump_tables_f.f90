! SPDX-License-Identifier: GPL-3.0-or-later
! WRF's Kain-Fritsch lookup tables: kf_lutab of module_cu_kfeta.o with the saturation constants of
! module_model_constants (as kfinit gets them); dump_tables_c.c writes the module's private tables.
PROGRAM dump_tables_f
  USE module_cu_kfeta, ONLY: kf_lutab
  USE module_model_constants, ONLY: svp1, svp2, svp3, svpt0
  IMPLICIT NONE
  INTERFACE
    SUBROUTINE dump_kf_tables()
    END SUBROUTINE dump_kf_tables
  END INTERFACE

  CALL kf_lutab(svp1, svp2, svp3, svpt0)
  CALL dump_kf_tables()
  PRINT *, "WRF kf_lutab tables written to tables_fortran.bin"
END PROGRAM dump_tables_f
