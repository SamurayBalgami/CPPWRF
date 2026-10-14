PROGRAM dump_tables_f
  USE module_mp_thompson
  IMPLICIT NONE

  REAL, DIMENSION(2,2,2) :: hgt
  INTEGER :: u

  hgt = 100.0

  PRINT *, "Starting thompson_init..."
  CALL thompson_init(hgt=hgt, &
                     ids=1, ide=2, jds=1, jde=2, kds=1, kde=2, &
                     ims=1, ime=2, jms=1, jme=2, kms=1, kme=2, &
                     its=1, ite=2, jts=1, jte=2, kts=1, kte=2)
  PRINT *, "thompson_init finished successfully!"

  OPEN(NEWUNIT=u, FILE="tables_fortran.bin", FORM="unformatted", ACCESS="stream", STATUS="replace")

  ! 1. Graupel-rain collection (5 tables, 37x37x1x37x37)
  WRITE(u) tcg_racg
  WRITE(u) tmr_racg
  WRITE(u) tcr_gacr
  WRITE(u) tnr_racg
  WRITE(u) tnr_gacr

  ! 2. Snow-rain collection (12 tables, 37x9x37x37)
  WRITE(u) tcs_racs1
  WRITE(u) tmr_racs1
  WRITE(u) tcs_racs2
  WRITE(u) tmr_racs2
  WRITE(u) tcr_sacr1
  WRITE(u) tms_sacr1
  WRITE(u) tcr_sacr2
  WRITE(u) tms_sacr2
  WRITE(u) tnr_racs1
  WRITE(u) tnr_racs2
  WRITE(u) tnr_sacr1
  WRITE(u) tnr_sacr2

  ! 3. Rain freezing (4 tables, 37x37x45x55)
  WRITE(u) tpi_qrfz
  WRITE(u) tpg_qrfz
  WRITE(u) tni_qrfz
  WRITE(u) tnr_qrfz

  ! 4. Cloud water freezing (2 tables, 37x100x45x55)
  WRITE(u) tpi_qcfz
  WRITE(u) tni_qcfz

  ! 5. Ice to snow autoconversion (3 tables, 64x55)
  WRITE(u) tps_iaus
  WRITE(u) tni_iaus
  WRITE(u) tpi_ide

  ! 6. Collision efficiencies (2 tables, 100x100)
  WRITE(u) t_Efrw
  WRITE(u) t_Efsw

  ! 7. Droplet evaporation (2 tables, 100x37x100)
  WRITE(u) tpc_wev
  WRITE(u) tnc_wev

  CLOSE(u)

  PRINT *, "Dumped all Fortran tables to tables_fortran.bin successfully."
  PRINT *, "Sample check:"
  PRINT *, "tcg_racg(10,10,1,10,10) = ", tcg_racg(10,10,1,10,10)
  PRINT *, "tcs_racs1(10,5,10,10)   = ", tcs_racs1(10,5,10,10)
  PRINT *, "tpi_qrfz(10,10,20,20)   = ", tpi_qrfz(10,10,20,20)
  PRINT *, "tpi_qcfz(10,50,20,20)   = ", tpi_qcfz(10,50,20,20)
  PRINT *, "tps_iaus(20,20)         = ", tps_iaus(20,20)
  PRINT *, "t_Efrw(50,50)           = ", t_Efrw(50,50)
  PRINT *, "tpc_wev(50,20,50)       = ", tpc_wev(50,20,50)

END PROGRAM dump_tables_f
