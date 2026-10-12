program driver_kessler
  use module_mp_kessler, only: kessler
  implicit none

  character(len=256) :: infile, outfile
  integer :: n_cols, nz
  real*8 :: dt, xlv, cp, ep2, svp1, svp2, svp3, svpt0, rhowater

  real*8, allocatable, dimension(:,:,:) :: t, qv, qc, qr, rho, pii, dz8w, z
  real*8, allocatable, dimension(:,:) :: rainnc, rainncv
  integer :: nargs

  nargs = command_argument_count()
  if (nargs < 2) then
    infile = 'kessler_input.bin'
    outfile = 'kessler_output_fortran.bin'
  else
    call get_command_argument(1, infile)
    call get_command_argument(2, outfile)
  end if

  open(unit=10, file=trim(infile), status='old', access='stream', form='unformatted')
  read(10) n_cols, nz, dt, xlv, cp, ep2, svp1, svp2, svp3, svpt0, rhowater

  allocate(t(n_cols, nz, 1))
  allocate(qv(n_cols, nz, 1))
  allocate(qc(n_cols, nz, 1))
  allocate(qr(n_cols, nz, 1))
  allocate(rho(n_cols, nz, 1))
  allocate(pii(n_cols, nz, 1))
  allocate(dz8w(n_cols, nz, 1))
  allocate(z(n_cols, nz, 1))
  allocate(rainnc(n_cols, 1))
  allocate(rainncv(n_cols, 1))

  read(10) t
  read(10) qv
  read(10) qc
  read(10) qr
  read(10) rho
  read(10) pii
  read(10) dz8w
  read(10) z
  read(10) rainnc
  read(10) rainncv
  close(10)

  call kessler( t, qv, qc, qr, rho, pii, dt, z, xlv, cp, &
                ep2, svp1, svp2, svp3, svpt0, rhowater, dz8w, &
                rainnc, rainncv, &
                1, n_cols, 1, 1, 1, nz, &
                1, n_cols, 1, 1, 1, nz, &
                1, n_cols, 1, 1, 1, nz )

  open(unit=20, file=trim(outfile), status='replace', access='stream', form='unformatted')
  write(20) n_cols, nz
  write(20) t
  write(20) qv
  write(20) qc
  write(20) qr
  write(20) rainnc
  write(20) rainncv
  close(20)

  deallocate(t, qv, qc, qr, rho, pii, dz8w, z, rainnc, rainncv)

end program driver_kessler
