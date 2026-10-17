program driver_bl_myj_fortran
  use module_bl_myjpbl, only: myjpbl
  implicit none

  character(len=256) :: infile, outfile
  integer :: n_cols, nz, stepbl
  real(8) :: dt

  ! 3D atmospheric arrays: (n_cols, nz, 1) or (n_cols, nz+1, 1)
  real(8), allocatable :: dz(:,:,:), pmid(:,:,:), pint(:,:,:), th(:,:,:)
  real(8), allocatable :: t(:,:,:), exner(:,:,:), qv(:,:,:), qc(:,:,:), qi(:,:,:)
  real(8), allocatable :: u(:,:,:), v(:,:,:), rho(:,:,:), tke_myj(:,:,:), exch_h(:,:,:)

  ! 2D arrays: (n_cols, 1)
  real(8), allocatable :: ht(:,:), tsk(:,:), xland(:,:), sice(:,:), snow(:,:)
  real(8), allocatable :: chklowq(:,:), ustar(:,:), znt(:,:), akhs(:,:), akms(:,:), elflx(:,:)
  integer, allocatable :: lowlyr(:,:)

  ! Inouts & Outputs
  real(8), allocatable :: qsfc(:,:), thz0(:,:), qz0(:,:), uz0(:,:), vz0(:,:), ct(:,:)
  real(8), allocatable :: pblh(:,:), mixht(:,:)
  integer, allocatable :: kpbl(:,:)

  ! 3D Outputs: (n_cols, nz, 1)
  real(8), allocatable :: el_myj(:,:,:), rublten(:,:,:), rvblten(:,:,:), rthblten(:,:,:)
  real(8), allocatable :: rqvblten(:,:,:), rqcblten(:,:,:), rqiblten(:,:,:)

  ! Dummy optional hydrometeors (not used by MYJ level 2.5)
  real(8), allocatable :: qs_dummy(:,:,:), qr_dummy(:,:,:), qg_dummy(:,:,:)
  real(8), allocatable :: rqs_dummy(:,:,:), rqr_dummy(:,:,:), rqg_dummy(:,:,:)

  integer :: i, k, iunit, ounit

  if (command_argument_count() >= 2) then
    call get_command_argument(1, infile)
    call get_command_argument(2, outfile)
  else
    infile = "bl_myj_input.bin"
    outfile = "bl_myj_output_fortran.bin"
  end if

  open(newunit=iunit, file=trim(infile), status='old', access='stream', form='unformatted')
  read(iunit) n_cols, nz, stepbl
  read(iunit) dt

  ! Allocate arrays
  allocate(dz(n_cols, nz, 1), pmid(n_cols, nz, 1), pint(n_cols, nz+1, 1), th(n_cols, nz, 1))
  allocate(t(n_cols, nz, 1), exner(n_cols, nz, 1), qv(n_cols, nz, 1), qc(n_cols, nz, 1), qi(n_cols, nz, 1))
  allocate(u(n_cols, nz, 1), v(n_cols, nz, 1), rho(n_cols, nz, 1), tke_myj(n_cols, nz, 1), exch_h(n_cols, nz, 1))

  allocate(ht(n_cols, 1), tsk(n_cols, 1), xland(n_cols, 1), sice(n_cols, 1), snow(n_cols, 1))
  allocate(chklowq(n_cols, 1), lowlyr(n_cols, 1), ustar(n_cols, 1), znt(n_cols, 1))
  allocate(akhs(n_cols, 1), akms(n_cols, 1), elflx(n_cols, 1))

  allocate(qsfc(n_cols, 1), thz0(n_cols, 1), qz0(n_cols, 1), uz0(n_cols, 1), vz0(n_cols, 1), ct(n_cols, 1))
  allocate(pblh(n_cols, 1), kpbl(n_cols, 1), mixht(n_cols, 1))

  allocate(el_myj(n_cols, nz, 1), rublten(n_cols, nz, 1), rvblten(n_cols, nz, 1), rthblten(n_cols, nz, 1))
  allocate(rqvblten(n_cols, nz, 1), rqcblten(n_cols, nz, 1), rqiblten(n_cols, nz, 1))

  ! Read 3D arrays
  read(iunit) ((dz(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((pmid(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((pint(i, k, 1), i=1,n_cols), k=1,nz+1)
  read(iunit) ((th(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((t(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((exner(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((qv(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((qc(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((qi(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((u(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((v(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((rho(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((tke_myj(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((exch_h(i, k, 1), i=1,n_cols), k=1,nz)

  ! Read 2D arrays
  read(iunit) (ht(i, 1), i=1,n_cols)
  read(iunit) (tsk(i, 1), i=1,n_cols)
  read(iunit) (xland(i, 1), i=1,n_cols)
  read(iunit) (sice(i, 1), i=1,n_cols)
  read(iunit) (snow(i, 1), i=1,n_cols)
  read(iunit) (chklowq(i, 1), i=1,n_cols)
  read(iunit) (lowlyr(i, 1), i=1,n_cols)
  read(iunit) (ustar(i, 1), i=1,n_cols)
  read(iunit) (znt(i, 1), i=1,n_cols)
  read(iunit) (akhs(i, 1), i=1,n_cols)
  read(iunit) (akms(i, 1), i=1,n_cols)
  read(iunit) (elflx(i, 1), i=1,n_cols)

  ! Read Inouts
  read(iunit) (qsfc(i, 1), i=1,n_cols)
  read(iunit) (thz0(i, 1), i=1,n_cols)
  read(iunit) (qz0(i, 1), i=1,n_cols)
  read(iunit) (uz0(i, 1), i=1,n_cols)
  read(iunit) (vz0(i, 1), i=1,n_cols)
  read(iunit) (ct(i, 1), i=1,n_cols)

  close(iunit)

  ! Run MYJPBL from WRF's compiled libwrflib.a
  call myjpbl(dt, stepbl, ht, dz, &
              pmid, pint, th, t, exner, qv, qc, qi, qs_dummy, qr_dummy, qg_dummy, &
              u, v, rho, tsk, qsfc, chklowq, thz0, qz0, uz0, vz0, &
              lowlyr, xland, sice, snow, &
              tke_myj, exch_h, ustar, znt, el_myj, pblh, kpbl, ct, &
              akhs, akms, elflx, mixht, &
              rublten, rvblten, rthblten, rqvblten, rqcblten, &
              rqiblten, rqs_dummy, rqr_dummy, rqg_dummy, &
              1, n_cols, 1, 1, 1, nz, &
              1, n_cols, 1, 1, 1, nz, &
              1, n_cols, 1, 1, 1, nz)

  print *, "Fortran Col 1: pblh =", pblh(1,1), "kpbl =", kpbl(1,1), "tke_sfc =", tke_myj(1,1,1)

  ! Write output binary
  open(newunit=ounit, file=trim(outfile), status='replace', access='stream', form='unformatted')
  ! 2D variables (7 variables, each n_cols)
  write(ounit) (qsfc(i, 1), i=1,n_cols)
  write(ounit) (thz0(i, 1), i=1,n_cols)
  write(ounit) (qz0(i, 1), i=1,n_cols)
  write(ounit) (ct(i, 1), i=1,n_cols)
  write(ounit) (pblh(i, 1), i=1,n_cols)
  write(ounit) (real(kpbl(i, 1), 8), i=1,n_cols)
  write(ounit) (mixht(i, 1), i=1,n_cols)

  ! 3D variables (9 variables, each n_cols * nz, written level by level)
  write(ounit) ((tke_myj(i, k, 1), i=1,n_cols), k=1,nz)
  write(ounit) ((exch_h(i, k, 1), i=1,n_cols), k=1,nz)
  write(ounit) ((el_myj(i, k, 1), i=1,n_cols), k=1,nz)
  write(ounit) ((rublten(i, k, 1), i=1,n_cols), k=1,nz)
  write(ounit) ((rvblten(i, k, 1), i=1,n_cols), k=1,nz)
  write(ounit) ((rthblten(i, k, 1), i=1,n_cols), k=1,nz)
  write(ounit) ((rqvblten(i, k, 1), i=1,n_cols), k=1,nz)
  write(ounit) ((rqcblten(i, k, 1), i=1,n_cols), k=1,nz)
  write(ounit) ((rqiblten(i, k, 1), i=1,n_cols), k=1,nz)
  close(ounit)

end program driver_bl_myj_fortran
