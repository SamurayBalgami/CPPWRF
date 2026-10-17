program driver_sf_myjsfc_fortran
  use module_sf_myjsfc, only: myjsfc, myjsfcinit
  implicit none

  character(len=256) :: infile, outfile
  integer :: n_cols, nz, itimestep, isurban, iz0tlnd
  real(8) :: p1000mb

  ! 3D atmospheric arrays: (n_cols, nz, 1) or (n_cols, nz+1, 1)
  real(8), allocatable :: dz(:,:,:), pmid(:,:,:), pint(:,:,:), th(:,:,:)
  real(8), allocatable :: t(:,:,:), qv(:,:,:), qc(:,:,:), u(:,:,:), v(:,:,:), q2(:,:,:)

  ! 2D arrays: (n_cols, 1)
  real(8), allocatable :: ht(:,:), tsk(:,:), xland(:,:), z0base(:,:), mavail(:,:)
  integer, allocatable :: lowlyr(:,:), ivgtyp(:,:)
  real(8), allocatable :: seamask(:,:), xice(:,:)

  ! Inout arrays: (n_cols, 1)
  real(8), allocatable :: qsfc(:,:), thz0(:,:), qz0(:,:), uz0(:,:), vz0(:,:)
  real(8), allocatable :: ustar(:,:), znt(:,:), pblh(:,:), rmol(:,:), akhs(:,:), akms(:,:)

  ! Output arrays: (n_cols, 1)
  real(8), allocatable :: rib(:,:), chs(:,:), chs2(:,:), cqs2(:,:), hfx(:,:), qfx(:,:)
  real(8), allocatable :: flx_lh(:,:), flhc(:,:), flqc(:,:), qgh(:,:), cpm(:,:), ct(:,:)
  real(8), allocatable :: u10(:,:), v10(:,:), t02(:,:), th02(:,:), tshltr(:,:)
  real(8), allocatable :: th10(:,:), q02(:,:), qshltr(:,:), q10(:,:), pshltr(:,:)
  real(8), allocatable :: u10e(:,:), v10e(:,:)

  integer :: i, k, iunit, ounit

  if (command_argument_count() >= 2) then
    call get_command_argument(1, infile)
    call get_command_argument(2, outfile)
  else
    infile = "sf_myjsfc_input.bin"
    outfile = "sf_myjsfc_output_fortran.bin"
  end if

  open(newunit=iunit, file=trim(infile), status='old', access='stream', form='unformatted')
  read(iunit) n_cols, nz, itimestep, isurban, iz0tlnd
  read(iunit) p1000mb

  allocate(dz(n_cols, nz, 1), pmid(n_cols, nz, 1), pint(n_cols, nz+1, 1), th(n_cols, nz, 1))
  allocate(t(n_cols, nz, 1), qv(n_cols, nz, 1), qc(n_cols, nz, 1))
  allocate(u(n_cols, nz, 1), v(n_cols, nz, 1), q2(n_cols, nz, 1))

  allocate(ht(n_cols, 1), tsk(n_cols, 1), xland(n_cols, 1), z0base(n_cols, 1), mavail(n_cols, 1))
  allocate(lowlyr(n_cols, 1), ivgtyp(n_cols, 1))
  allocate(seamask(n_cols, 1), xice(n_cols, 1))

  allocate(qsfc(n_cols, 1), thz0(n_cols, 1), qz0(n_cols, 1), uz0(n_cols, 1), vz0(n_cols, 1))
  allocate(ustar(n_cols, 1), znt(n_cols, 1), pblh(n_cols, 1), rmol(n_cols, 1), akhs(n_cols, 1), akms(n_cols, 1))

  allocate(rib(n_cols, 1), chs(n_cols, 1), chs2(n_cols, 1), cqs2(n_cols, 1), hfx(n_cols, 1), qfx(n_cols, 1))
  allocate(flx_lh(n_cols, 1), flhc(n_cols, 1), flqc(n_cols, 1), qgh(n_cols, 1), cpm(n_cols, 1), ct(n_cols, 1))
  allocate(u10(n_cols, 1), v10(n_cols, 1), t02(n_cols, 1), th02(n_cols, 1), tshltr(n_cols, 1))
  allocate(th10(n_cols, 1), q02(n_cols, 1), qshltr(n_cols, 1), q10(n_cols, 1), pshltr(n_cols, 1))
  allocate(u10e(n_cols, 1), v10e(n_cols, 1))

  read(iunit) ((dz(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((pmid(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((pint(i, k, 1), i=1,n_cols), k=1,nz+1)
  read(iunit) ((th(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((t(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((qv(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((qc(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((u(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((v(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((q2(i, k, 1), i=1,n_cols), k=1,nz)

  read(iunit) (ht(i, 1), i=1,n_cols)
  read(iunit) (tsk(i, 1), i=1,n_cols)
  read(iunit) (xland(i, 1), i=1,n_cols)
  read(iunit) (z0base(i, 1), i=1,n_cols)
  read(iunit) (mavail(i, 1), i=1,n_cols)
  read(iunit) (lowlyr(i, 1), i=1,n_cols)
  read(iunit) (ivgtyp(i, 1), i=1,n_cols)

  read(iunit) (qsfc(i, 1), i=1,n_cols)
  read(iunit) (thz0(i, 1), i=1,n_cols)
  read(iunit) (qz0(i, 1), i=1,n_cols)
  read(iunit) (uz0(i, 1), i=1,n_cols)
  read(iunit) (vz0(i, 1), i=1,n_cols)
  read(iunit) (ustar(i, 1), i=1,n_cols)
  read(iunit) (znt(i, 1), i=1,n_cols)
  read(iunit) (pblh(i, 1), i=1,n_cols)
  read(iunit) (rmol(i, 1), i=1,n_cols)
  read(iunit) (akhs(i, 1), i=1,n_cols)
  read(iunit) (akms(i, 1), i=1,n_cols)

  close(iunit)

  do i = 1, n_cols
    seamask(i, 1) = xland(i, 1) - 1.0_8
    xice(i, 1) = 0.0_8
  end do

  ! Initialize tables via WRF's compiled libwrflib.a
  call myjsfcinit(lowlyr, ustar, znt, seamask, xice, ivgtyp, .true., .false., &
                  1, n_cols, 1, 1, 1, nz, &
                  1, n_cols, 1, 1, 1, nz, &
                  1, n_cols, 1, 1, 1, nz)

  ! Run MYJSFC via WRF's compiled libwrflib.a
  call myjsfc(itimestep, ht, dz, &
              pmid, pint, th, t, qv, qc, u, v, q2, &
              tsk, qsfc, thz0, qz0, uz0, vz0, &
              lowlyr, xland, ivgtyp, isurban, iz0tlnd, &
              ustar, znt, z0base, pblh, mavail, rmol, &
              akhs, akms, &
              rib, &
              chs, chs2, cqs2, hfx, qfx, flx_lh, flhc, flqc, &
              qgh, cpm, ct, &
              u10, v10, t02, th02, tshltr, th10, q02, qshltr, q10, pshltr, &
              p1000mb, u10e, v10e, &
              1, n_cols, 1, 1, 1, nz, &
              1, n_cols, 1, 1, 1, nz, &
              1, n_cols, 1, 1, 1, nz)

  print *, "Fortran Col 1: chs2 =", chs2(1,1), "chs =", chs(1,1), "ustar =", ustar(1,1)

  open(newunit=ounit, file=trim(outfile), status='replace', access='stream', form='unformatted')
  ! 11 inout variables
  write(ounit) (qsfc(i, 1), i=1,n_cols)
  write(ounit) (thz0(i, 1), i=1,n_cols)
  write(ounit) (qz0(i, 1), i=1,n_cols)
  write(ounit) (uz0(i, 1), i=1,n_cols)
  write(ounit) (vz0(i, 1), i=1,n_cols)
  write(ounit) (ustar(i, 1), i=1,n_cols)
  write(ounit) (znt(i, 1), i=1,n_cols)
  write(ounit) (pblh(i, 1), i=1,n_cols)
  write(ounit) (rmol(i, 1), i=1,n_cols)
  write(ounit) (akhs(i, 1), i=1,n_cols)
  write(ounit) (akms(i, 1), i=1,n_cols)

  ! 24 output variables
  write(ounit) (rib(i, 1), i=1,n_cols)
  write(ounit) (chs(i, 1), i=1,n_cols)
  write(ounit) (chs2(i, 1), i=1,n_cols)
  write(ounit) (cqs2(i, 1), i=1,n_cols)
  write(ounit) (hfx(i, 1), i=1,n_cols)
  write(ounit) (qfx(i, 1), i=1,n_cols)
  write(ounit) (flx_lh(i, 1), i=1,n_cols)
  write(ounit) (flhc(i, 1), i=1,n_cols)
  write(ounit) (flqc(i, 1), i=1,n_cols)
  write(ounit) (qgh(i, 1), i=1,n_cols)
  write(ounit) (cpm(i, 1), i=1,n_cols)
  write(ounit) (ct(i, 1), i=1,n_cols)
  write(ounit) (u10(i, 1), i=1,n_cols)
  write(ounit) (v10(i, 1), i=1,n_cols)
  write(ounit) (t02(i, 1), i=1,n_cols)
  write(ounit) (th02(i, 1), i=1,n_cols)
  write(ounit) (tshltr(i, 1), i=1,n_cols)
  write(ounit) (th10(i, 1), i=1,n_cols)
  write(ounit) (q02(i, 1), i=1,n_cols)
  write(ounit) (qshltr(i, 1), i=1,n_cols)
  write(ounit) (q10(i, 1), i=1,n_cols)
  write(ounit) (pshltr(i, 1), i=1,n_cols)
  write(ounit) (u10e(i, 1), i=1,n_cols)
  write(ounit) (v10e(i, 1), i=1,n_cols)
  close(ounit)

end program driver_sf_myjsfc_fortran
