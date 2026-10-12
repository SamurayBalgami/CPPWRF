MODULE module_wrf_error
  IMPLICIT NONE
CONTAINS
  SUBROUTINE wrf_error_fatal(msg)
    CHARACTER(LEN=*), INTENT(IN) :: msg
    WRITE(*,*) 'FATAL ERROR in KF: ', TRIM(msg)
    STOP 1
  END SUBROUTINE wrf_error_fatal

  SUBROUTINE wrf_message(msg)
    CHARACTER(LEN=*), INTENT(IN) :: msg
    ! silent
  END SUBROUTINE wrf_message
END MODULE module_wrf_error
