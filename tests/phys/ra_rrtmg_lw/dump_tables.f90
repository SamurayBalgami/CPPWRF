PROGRAM dump_tables
  USE module_ra_rrtmg_lw, only: rrtmg_lwinit
  USE rrlw_con
  USE rrlw_wvn
  USE rrlw_tbl
  USE rrlw_cld
  USE rrlw_ref
  IMPLICIT NONE

  INTEGER :: uout, ios

  ! Initialize tables by reading RRTMG_LW_DATA
  CALL rrtmg_lwinit(5000.0_8, .TRUE., 1, 1, 1, 1, 1, 50, 1, 1, 1, 1, 1, 50, 1, 1, 1, 1, 1, 50)

  OPEN(NEWUNIT=uout, FILE="rrtmg_lw_tables.bin", ACCESS='STREAM', FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN', STATUS='REPLACE', IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'Error opening rrtmg_lw_tables.bin'
    STOP 1
  END IF

  ! Lookup tables from rrlw_tbl: tau_tbl(0:10000), exp_tbl(0:10000), tfn_tbl(0:10000)
  WRITE(uout) tau_tbl
  WRITE(uout) exp_tbl
  WRITE(uout) tfn_tbl

  ! Spectral info from rrlw_wvn: totplnk(181, 16), totplk16(181), delwave(16), wavenum1(16), wavenum2(16)
  WRITE(uout) totplnk
  WRITE(uout) totplk16
  WRITE(uout) delwave
  WRITE(uout) wavenum1
  WRITE(uout) wavenum2
  WRITE(uout) ngc
  WRITE(uout) ngs
  WRITE(uout) ngn
  WRITE(uout) ngb
  WRITE(uout) ngm
  WRITE(uout) wt
  WRITE(uout) rwgt

  ! Cloud tables from rrlw_cld: abscld1, absice0, absice1, absice2, absice3, absliq0, absliq1
  WRITE(uout) abscld1
  WRITE(uout) absice0
  WRITE(uout) absice1
  WRITE(uout) absice2
  WRITE(uout) absice3
  WRITE(uout) absliq0
  WRITE(uout) absliq1

  ! Reference profiles from rrlw_ref: preflog(59), tref(59), chi_mls(7, 59)
  WRITE(uout) preflog
  WRITE(uout) tref
  WRITE(uout) chi_mls

  CALL dump_b01(uout)
  CALL dump_b02(uout)
  CALL dump_b03(uout)
  CALL dump_b04(uout)
  CALL dump_b05(uout)
  CALL dump_b06(uout)
  CALL dump_b07(uout)
  CALL dump_b08(uout)
  CALL dump_b09(uout)
  CALL dump_b10(uout)
  CALL dump_b11(uout)
  CALL dump_b12(uout)
  CALL dump_b13(uout)
  CALL dump_b14(uout)
  CALL dump_b15(uout)
  CALL dump_b16(uout)

  CLOSE(uout)
  WRITE(*,*) 'Successfully dumped rrtmg_lw_tables.bin'

CONTAINS

  SUBROUTINE dump_b01(u)
    USE rrlw_kg01
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, ka_mn2, kb_mn2, selfref, forref
  END SUBROUTINE dump_b01

  SUBROUTINE dump_b02(u)
    USE rrlw_kg02
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, selfref, forref
  END SUBROUTINE dump_b02

  SUBROUTINE dump_b03(u)
    USE rrlw_kg03
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, ka_mn2o, kb_mn2o, selfref, forref
  END SUBROUTINE dump_b03

  SUBROUTINE dump_b04(u)
    USE rrlw_kg04
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, selfref, forref
  END SUBROUTINE dump_b04

  SUBROUTINE dump_b05(u)
    USE rrlw_kg05
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, ka_mo3, ccl4, selfref, forref
  END SUBROUTINE dump_b05

  SUBROUTINE dump_b06(u)
    USE rrlw_kg06
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, ka, ka_mco2, cfc11adj, cfc12, selfref, forref
  END SUBROUTINE dump_b06

  SUBROUTINE dump_b07(u)
    USE rrlw_kg07
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, ka_mco2, kb_mco2, selfref, forref
  END SUBROUTINE dump_b07

  SUBROUTINE dump_b08(u)
    USE rrlw_kg08
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, ka_mco2, kb_mco2, ka_mn2o, kb_mn2o, ka_mo3, cfc12, cfc22adj, selfref, forref
  END SUBROUTINE dump_b08

  SUBROUTINE dump_b09(u)
    USE rrlw_kg09
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, ka_mn2o, kb_mn2o, selfref, forref
  END SUBROUTINE dump_b09

  SUBROUTINE dump_b10(u)
    USE rrlw_kg10
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, selfref, forref
  END SUBROUTINE dump_b10

  SUBROUTINE dump_b11(u)
    USE rrlw_kg11
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, ka_mo2, kb_mo2, selfref, forref
  END SUBROUTINE dump_b11

  SUBROUTINE dump_b12(u)
    USE rrlw_kg12
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, ka, selfref, forref
  END SUBROUTINE dump_b12

  SUBROUTINE dump_b13(u)
    USE rrlw_kg13
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, ka_mco2, ka_mco, kb_mo3, selfref, forref
  END SUBROUTINE dump_b13

  SUBROUTINE dump_b14(u)
    USE rrlw_kg14
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, selfref, forref
  END SUBROUTINE dump_b14

  SUBROUTINE dump_b15(u)
    USE rrlw_kg15
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, ka, ka_mn2, selfref, forref
  END SUBROUTINE dump_b15

  SUBROUTINE dump_b16(u)
    USE rrlw_kg16
    INTEGER, INTENT(IN) :: u
    WRITE(u) fracrefa, fracrefb, ka, kb, selfref, forref
  END SUBROUTINE dump_b16

END PROGRAM dump_tables
