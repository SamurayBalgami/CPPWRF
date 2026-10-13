! Standalone WRF support modules for RRTMG differential testing
MODULE module_model_constants
  IMPLICIT NONE
  REAL, PARAMETER :: r_d = 287.0
  REAL, PARAMETER :: cp  = 7.0 * r_d / 2.0  ! 1004.5 J/(kg K)
END MODULE module_model_constants

MODULE module_state_description
  IMPLICIT NONE
  INTEGER, PARAMETER :: FER_MP_HIRES        = 5
  INTEGER, PARAMETER :: FER_MP_HIRES_ADVECT = 15
  INTEGER, PARAMETER :: ETAMP_HWRF          = 8
END MODULE module_state_description

MODULE module_wrf_error
  IMPLICIT NONE
CONTAINS
  SUBROUTINE wrf_error_fatal(msg)
    CHARACTER(LEN=*), INTENT(IN) :: msg
    WRITE(*,*) 'FATAL: ', TRIM(msg)
    STOP 1
  END SUBROUTINE wrf_error_fatal

  SUBROUTINE wrf_message(msg)
    CHARACTER(LEN=*), INTENT(IN) :: msg
  END SUBROUTINE wrf_message

  SUBROUTINE wrf_debug(level, msg)
    INTEGER, INTENT(IN) :: level
    CHARACTER(LEN=*), INTENT(IN) :: msg
  END SUBROUTINE wrf_debug
END MODULE module_wrf_error

MODULE module_ra_clWRF_support
  IMPLICIT NONE
CONTAINS
  SUBROUTINE read_CAMgases(yr, julian, READtrFILE, model, co2vmr, n2ovmr, ch4vmr, cfc11vmr, cfc12vmr)
    INTEGER, INTENT(IN)              :: yr
    REAL, INTENT(IN)                 :: julian
    CHARACTER(LEN=*), INTENT(IN)     :: model
    REAL(8), INTENT(OUT)            :: co2vmr, n2ovmr, ch4vmr, cfc11vmr, cfc12vmr
    LOGICAL                          :: READtrFILE
    co2vmr   = 379.e-6_8
    n2ovmr   = 319.e-9_8
    ch4vmr   = 1774.e-9_8
    cfc11vmr = 251.e-12_8
    cfc12vmr = 538.e-12_8
  END SUBROUTINE read_CAMgases
END MODULE module_ra_clWRF_support

LOGICAL FUNCTION wrf_dm_on_monitor()
  IMPLICIT NONE
  wrf_dm_on_monitor = .TRUE.
END FUNCTION wrf_dm_on_monitor

SUBROUTINE wrf_dm_bcast_bytes(buf, nbytes)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nbytes
  INTEGER, INTENT(INOUT) :: buf(*)
END SUBROUTINE wrf_dm_bcast_bytes
