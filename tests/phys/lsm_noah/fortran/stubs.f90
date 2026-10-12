MODULE module_model_constants
  IMPLICIT NONE
  REAL, PARAMETER :: CP = 1004.5
  REAL, PARAMETER :: R_D = 287.0
  REAL, PARAMETER :: XLF = 3.50E5
  REAL, PARAMETER :: XLV = 2.5E6
  REAL, PARAMETER :: RHOWATER = 1000.0
  REAL, PARAMETER :: STBOLT = 5.67051E-8
  REAL, PARAMETER :: KARMAN = 0.4
END MODULE module_model_constants

MODULE module_wrf_error
  IMPLICIT NONE
CONTAINS
  SUBROUTINE wrf_error_fatal(msg)
    CHARACTER(LEN=*), INTENT(IN) :: msg
    WRITE(*,*) 'FATAL ERROR: ', TRIM(msg)
    STOP 1
  END SUBROUTINE wrf_error_fatal

  SUBROUTINE wrf_message(msg)
    CHARACTER(LEN=*), INTENT(IN) :: msg
    ! silent in harness or write(*,*) TRIM(msg)
  END SUBROUTINE wrf_message
END MODULE module_wrf_error
