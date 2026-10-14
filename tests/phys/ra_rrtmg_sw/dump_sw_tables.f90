PROGRAM dump_sw_tables
  USE module_ra_rrtmg_sw, only: rrtmg_swinit
  USE rrsw_con
  USE rrsw_wvn
  USE rrsw_tbl
  USE rrsw_cld
  USE rrsw_ref
  USE rrsw_aer
  IMPLICIT NONE

  INTEGER :: uout, ios

  ! Initialize tables by reading RRTMG_SW_DATA
  CALL rrtmg_swinit(.TRUE., 1, 1, 1, 1, 1, 50, 1, 1, 1, 1, 1, 50, 1, 1, 1, 1, 1, 50)

  OPEN(NEWUNIT=uout, FILE="rrtmg_sw_tables.bin", ACCESS='STREAM', FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN', STATUS='REPLACE', IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'Error opening rrtmg_sw_tables.bin'
    STOP 1
  END IF

  ! 1. Lookup tables from rrsw_tbl: exp_tbl(0:10000)
  WRITE(uout) exp_tbl

  ! 2. Spectral info from rrsw_wvn
  WRITE(uout) wavenum1
  WRITE(uout) wavenum2
  WRITE(uout) delwave
  WRITE(uout) ngc
  WRITE(uout) ngs
  WRITE(uout) ngn
  WRITE(uout) ngb
  WRITE(uout) ngm
  WRITE(uout) wt
  WRITE(uout) rwgt

  ! 3. Cloud tables from rrsw_cld
  WRITE(uout) extliq1
  WRITE(uout) ssaliq1
  WRITE(uout) asyliq1
  WRITE(uout) extice2
  WRITE(uout) ssaice2
  WRITE(uout) asyice2
  WRITE(uout) extice3
  WRITE(uout) ssaice3
  WRITE(uout) asyice3
  WRITE(uout) fdlice3
  WRITE(uout) abari
  WRITE(uout) bbari
  WRITE(uout) cbari
  WRITE(uout) dbari
  WRITE(uout) ebari
  WRITE(uout) fbari

  ! 4. Reference profiles from rrsw_ref: preflog(59), tref(59)
  WRITE(uout) preflog
  WRITE(uout) tref

  ! 5. Aerosol properties from rrsw_aer
  WRITE(uout) rsrtaua
  WRITE(uout) rsrpiza
  WRITE(uout) rsrasya

  ! 6. Band absorption tables
  CALL dump_b16(uout)
  CALL dump_b17(uout)
  CALL dump_b18(uout)
  CALL dump_b19(uout)
  CALL dump_b20(uout)
  CALL dump_b21(uout)
  CALL dump_b22(uout)
  CALL dump_b23(uout)
  CALL dump_b24(uout)
  CALL dump_b25(uout)
  CALL dump_b26(uout)
  CALL dump_b27(uout)
  CALL dump_b28(uout)
  CALL dump_b29(uout)

  CLOSE(uout)
  WRITE(*,*) 'Successfully dumped rrtmg_sw_tables.bin'

CONTAINS

  SUBROUTINE dump_b16(u)
    USE rrsw_kg16
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, selfref, forref, sfluxref, rayl, strrat1, layreffr
  END SUBROUTINE dump_b16

  SUBROUTINE dump_b17(u)
    USE rrsw_kg17
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, selfref, forref, sfluxref, rayl, strrat, layreffr
  END SUBROUTINE dump_b17

  SUBROUTINE dump_b18(u)
    USE rrsw_kg18
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, selfref, forref, sfluxref, rayl, strrat, layreffr
  END SUBROUTINE dump_b18

  SUBROUTINE dump_b19(u)
    USE rrsw_kg19
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, selfref, forref, sfluxref, rayl, strrat, layreffr
  END SUBROUTINE dump_b19

  SUBROUTINE dump_b20(u)
    USE rrsw_kg20
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, selfref, forref, sfluxref, absch4, rayl, layreffr
  END SUBROUTINE dump_b20

  SUBROUTINE dump_b21(u)
    USE rrsw_kg21
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, selfref, forref, sfluxref, rayl, strrat, layreffr
  END SUBROUTINE dump_b21

  SUBROUTINE dump_b22(u)
    USE rrsw_kg22
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, selfref, forref, sfluxref, rayl, strrat, layreffr
  END SUBROUTINE dump_b22

  SUBROUTINE dump_b23(u)
    USE rrsw_kg23
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, selfref, forref, sfluxref, rayl, givfac, layreffr
  END SUBROUTINE dump_b23

  SUBROUTINE dump_b24(u)
    USE rrsw_kg24
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, selfref, forref, sfluxref, abso3a, abso3b, rayla, raylb, strrat, layreffr
  END SUBROUTINE dump_b24

  SUBROUTINE dump_b25(u)
    USE rrsw_kg25
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, sfluxref, abso3a, abso3b, rayl, layreffr
  END SUBROUTINE dump_b25

  SUBROUTINE dump_b26(u)
    USE rrsw_kg26
    INTEGER, INTENT(IN) :: u
    WRITE(u) sfluxref, rayl
  END SUBROUTINE dump_b26

  SUBROUTINE dump_b27(u)
    USE rrsw_kg27
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, sfluxref, rayl, scalekur, layreffr
  END SUBROUTINE dump_b27

  SUBROUTINE dump_b28(u)
    USE rrsw_kg28
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, sfluxref, rayl, strrat, layreffr
  END SUBROUTINE dump_b28

  SUBROUTINE dump_b29(u)
    USE rrsw_kg29
    INTEGER, INTENT(IN) :: u
    WRITE(u) ka, kb, selfref, forref, sfluxref, absh2o, absco2, rayl, layreffr
  END SUBROUTINE dump_b29

END PROGRAM dump_sw_tables
