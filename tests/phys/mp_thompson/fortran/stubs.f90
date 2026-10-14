MODULE module_wrf_error
  IMPLICIT NONE
  CHARACTER(LEN=512) :: wrf_err_message
CONTAINS
  SUBROUTINE wrf_error_fatal(msg)
    CHARACTER(LEN=*), INTENT(IN) :: msg
    WRITE(*,*) 'FATAL ERROR in Thompson: ', TRIM(msg)
    STOP 1
  END SUBROUTINE wrf_error_fatal

  SUBROUTINE wrf_message(msg)
    CHARACTER(LEN=*), INTENT(IN) :: msg
    WRITE(*,*) TRIM(msg)
  END SUBROUTINE wrf_message

  SUBROUTINE wrf_debug(level, msg)
    INTEGER, INTENT(IN) :: level
    CHARACTER(LEN=*), INTENT(IN) :: msg
    ! silent
  END SUBROUTINE wrf_debug
END MODULE module_wrf_error

MODULE module_model_constants
  IMPLICIT NONE
  REAL, PARAMETER :: RE_QC_BG = 2.49E-6
  REAL, PARAMETER :: RE_QI_BG = 4.99E-6
  REAL, PARAMETER :: RE_QS_BG = 9.99E-6
END MODULE module_model_constants

MODULE module_timing
  IMPLICIT NONE
CONTAINS
  SUBROUTINE start_timing
  END SUBROUTINE start_timing

  SUBROUTINE end_timing(str)
    CHARACTER(LEN=*), INTENT(IN) :: str
  END SUBROUTINE end_timing
END MODULE module_timing

MODULE module_domain
  IMPLICIT NONE
END MODULE module_domain

MODULE module_dm
  IMPLICIT NONE
END MODULE module_dm

SUBROUTINE nl_get_force_read_thompson(id, flag)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: id
  LOGICAL, INTENT(OUT) :: flag
  flag = .FALSE.
END SUBROUTINE nl_get_force_read_thompson

SUBROUTINE nl_get_write_thompson_tables(id, flag)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: id
  LOGICAL, INTENT(OUT) :: flag
  flag = .FALSE.
END SUBROUTINE nl_get_write_thompson_tables

SUBROUTINE nl_get_write_thompson_mp38table(id, flag)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: id
  LOGICAL, INTENT(OUT) :: flag
  flag = .FALSE.
END SUBROUTINE nl_get_write_thompson_mp38table

LOGICAL FUNCTION wrf_dm_on_monitor()
  IMPLICIT NONE
  wrf_dm_on_monitor = .TRUE.
END FUNCTION wrf_dm_on_monitor
