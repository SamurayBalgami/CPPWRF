PROGRAM dump_tables_f
  USE module_cu_kfeta, ONLY: kf_lutab, ttab, qstab, the0k, alu, rdpr, rdthk, plutop
  IMPLICIT NONE

  REAL :: svp1, svp2, svp3, svpt0
  INTEGER :: u

  svp1 = 0.6112
  svp2 = 17.67
  svp3 = 29.65
  svpt0 = 273.15

  CALL kf_lutab(svp1, svp2, svp3, svpt0)

  OPEN(NEWUNIT=u, FILE="tables_fortran.bin", FORM="unformatted", ACCESS="stream", STATUS="replace")
  WRITE(u) rdpr, rdthk, plutop
  WRITE(u) the0k
  WRITE(u) alu
  WRITE(u) ttab
  WRITE(u) qstab
  CLOSE(u)

  PRINT *, "Fortran tables dumped successfully."
  PRINT *, "rdpr =", rdpr, "rdthk =", rdthk, "plutop =", plutop
  PRINT *, "the0k(1) =", the0k(1), "the0k(220) =", the0k(220)
  PRINT *, "ttab(1,1) =", ttab(1,1), "ttab(250,220) =", ttab(250,220)
  PRINT *, "qstab(1,1) =", qstab(1,1), "qstab(250,220) =", qstab(250,220)
END PROGRAM dump_tables_f
