! The aerosol functions of WRF's module_mp_thompson (activ_ncloud, Eff_aero, iceDeMott, iceKoop) on the samples of
! aero_funcs.py: thompson_init (its table files in the working directory) for av_g / bv_g, then table_ccnAct
! (CCN_ACTIVATE.BIN). Built with WRF's BYTESWAPIO; the sample files are opened CONVERT='native'.
PROGRAM aero_funcs
  USE module_mp_thompson, ONLY: thompson_init, table_ccnAct, activ_ncloud, Eff_aero, iceDeMott, iceKoop
  IMPLICIT NONE
  INTEGER :: u, n, i
  REAL*8, ALLOCATABLE :: tt(:), ww(:), nccn(:), d(:), da(:), visc(:), rhoa(:), tk(:), sp(:)
  REAL*8, ALLOCATABLE :: tc(:), rho(:), nifa(:), tkp(:), qv(:), qvs(:), na(:), dt(:), out(:,:)
  REAL*8 :: hgt(1, 2, 1)
  CHARACTER(LEN=1) :: species(0:2) = (/'r', 's', 'g'/)
  CHARACTER(LEN=256) :: infile, outfile
  CALL GET_COMMAND_ARGUMENT(1, infile)
  CALL GET_COMMAND_ARGUMENT(2, outfile)
  OPEN(NEWUNIT=u, FILE=TRIM(infile), FORM="unformatted", ACCESS="stream", STATUS="old", CONVERT="native")
  READ(u) n
  ALLOCATE(tt(n), ww(n), nccn(n), d(n), da(n), visc(n), rhoa(n), tk(n), sp(n))
  ALLOCATE(tc(n), rho(n), nifa(n), tkp(n), qv(n), qvs(n), na(n), dt(n), out(n, 4))
  READ(u) tt, ww, nccn, d, da, visc, rhoa, tk, sp, tc, rho, nifa, tkp, qv, qvs, na, dt
  CLOSE(u)
  hgt = 100.0
  CALL nl_set_write_thompson_tables(1, .TRUE.)
  CALL thompson_init(hgt=hgt, ids=1, ide=2, jds=1, jde=2, kds=1, kde=2, ims=1, ime=1, jms=1, jme=1, kms=1, kme=2, &
                     its=1, ite=1, jts=1, jte=1, kts=1, kte=2)
  CALL table_ccnAct
  DO i = 1, n
    out(i, 1) = activ_ncloud(tt(i), ww(i), nccn(i))
    out(i, 2) = Eff_aero(d(i), da(i), visc(i), rhoa(i), tk(i), species(NINT(sp(i))))
    out(i, 3) = iceDeMott(tc(i), qvs(i), qvs(i), qvs(i), rho(i), nifa(i))
    out(i, 4) = iceKoop(tkp(i), qv(i), qvs(i), na(i), dt(i))
  END DO
  OPEN(NEWUNIT=u, FILE=TRIM(outfile), FORM="unformatted", ACCESS="stream", STATUS="replace", CONVERT="native")
  WRITE(u) out
  CLOSE(u)
END PROGRAM aero_funcs
