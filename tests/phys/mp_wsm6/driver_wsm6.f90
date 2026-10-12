program driver_wsm6
  use module_mp_wsm6, only: wsm6
  use mp_wsm6, only: mp_wsm6_init
  implicit none

  character(len=256) :: infile, outfile
  character(len=256) :: errmsg
  integer :: errflg
  integer :: n_cols, nz, hail_opt, nargs
  real*8 :: delt, g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin, xls, xlv0, xlf0
  real*8 :: den0, denr, dens, cliq, cice, psat

  real*8, allocatable, dimension(:,:,:) :: th, q, qc, qr, qi, qs, qg
  real*8, allocatable, dimension(:,:,:) :: den, pii, p, delz
  real*8, allocatable, dimension(:,:,:) :: re_cloud, re_ice, re_snow
  real*8, allocatable, dimension(:,:) :: rain, rainncv, snow, snowncv, graupel, graupelncv, sr
  real*8 :: re_qc_bg, re_qi_bg, re_qs_bg, re_qc_max, re_qi_max, re_qs_max
  integer :: has_reqc, has_reqi, has_reqs, do_radar_ref
  logical :: diagflag

  nargs = command_argument_count()
  if (nargs < 2) then
    infile = 'wsm6_input.bin'
    outfile = 'wsm6_output_fortran.bin'
  else
    call get_command_argument(1, infile)
    call get_command_argument(2, outfile)
  end if

  open(unit=10, file=trim(infile), status='old', access='stream', form='unformatted')
  read(10) n_cols, nz, delt, g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin, xls, xlv0, xlf0, &
           den0, denr, dens, cliq, cice, psat, hail_opt

  allocate(th(n_cols, nz, 1), q(n_cols, nz, 1), qc(n_cols, nz, 1), qr(n_cols, nz, 1), &
           qi(n_cols, nz, 1), qs(n_cols, nz, 1), qg(n_cols, nz, 1), &
           den(n_cols, nz, 1), pii(n_cols, nz, 1), p(n_cols, nz, 1), delz(n_cols, nz, 1), &
           re_cloud(n_cols, nz, 1), re_ice(n_cols, nz, 1), re_snow(n_cols, nz, 1), &
           rain(n_cols, 1), rainncv(n_cols, 1), snow(n_cols, 1), snowncv(n_cols, 1), &
           graupel(n_cols, 1), graupelncv(n_cols, 1), sr(n_cols, 1))

  read(10) th
  read(10) q
  read(10) qc
  read(10) qr
  read(10) qi
  read(10) qs
  read(10) qg
  read(10) den
  read(10) pii
  read(10) p
  read(10) delz
  read(10) rain
  read(10) rainncv
  read(10) snow
  read(10) snowncv
  read(10) graupel
  read(10) graupelncv
  read(10) sr
  close(10)

  ! Initialize WSM6
  call mp_wsm6_init(den0, denr, dens, cliq, cpv, hail_opt, errmsg, errflg)
  if (errflg /= 0) then
    print *, "Error in mp_wsm6_init: ", trim(errmsg)
    stop 1
  end if

  has_reqc = 0
  has_reqi = 0
  has_reqs = 0
  re_qc_bg = 0.d0
  re_qi_bg = 0.d0
  re_qs_bg = 0.d0
  re_qc_max = 0.d0
  re_qi_max = 0.d0
  re_qs_max = 0.d0
  diagflag = .false.
  do_radar_ref = 0

  call wsm6(th, q, qc, qr, qi, qs, qg, den, pii, p, delz, &
            delt, g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin, &
            xls, xlv0, xlf0, den0, denr, cliq, cice, psat, &
            rain, rainncv, snow, snowncv, graupel, graupelncv, sr, &
            diagflag=diagflag, do_radar_ref=do_radar_ref, &
            has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs, &
            re_qc_bg=re_qc_bg, re_qi_bg=re_qi_bg, re_qs_bg=re_qs_bg, &
            re_qc_max=re_qc_max, re_qi_max=re_qi_max, re_qs_max=re_qs_max, &
            re_cloud=re_cloud, re_ice=re_ice, re_snow=re_snow, &
            ids=1, ide=n_cols, jds=1, jde=1, kds=1, kde=nz, &
            ims=1, ime=n_cols, jms=1, jme=1, kms=1, kme=nz, &
            its=1, ite=n_cols, jts=1, jte=1, kts=1, kte=nz, &
            errmsg=errmsg, errflg=errflg)

  if (errflg /= 0) then
    print *, "Error in wsm6: ", trim(errmsg)
    stop 1
  end if

  open(unit=20, file=trim(outfile), status='replace', access='stream', form='unformatted')
  write(20) n_cols, nz
  write(20) th
  write(20) q
  write(20) qc
  write(20) qr
  write(20) qi
  write(20) qs
  write(20) qg
  write(20) rain
  write(20) rainncv
  write(20) snow
  write(20) snowncv
  write(20) graupel
  write(20) graupelncv
  write(20) sr
  close(20)

  deallocate(th, q, qc, qr, qi, qs, qg, den, pii, p, delz, &
             re_cloud, re_ice, re_snow, &
             rain, rainncv, snow, snowncv, graupel, graupelncv, sr)

end program driver_wsm6
