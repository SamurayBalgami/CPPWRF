program driver_ysu_fortran
  use module_bl_ysu, only: ysu
  use bl_ysu, only: bl_ysu_init
  implicit none

  character(len=256) :: infile, outfile
  integer :: n_cols, nz
  integer :: ysu_topdown_pblmix, flag_qc_i, flag_qi_i
  logical :: flag_qc, flag_qi
  real(8) :: dt
  real(8) :: cp, g, rovcp, rd, rovg, ep1, ep2, karman, xlv, rv

  ! 3D inputs
  real(8), allocatable :: u3d(:,:,:), v3d(:,:,:), t3d(:,:,:), qv3d(:,:,:)
  real(8), allocatable :: qc3d(:,:,:), qi3d(:,:,:), p3d(:,:,:), p3di(:,:,:), pi3d(:,:,:)
  real(8), allocatable :: dz8w(:,:,:), rthraten(:,:,:)

  ! 2D inputs
  real(8), allocatable :: psfc(:,:), znt(:,:), ust(:,:), psim(:,:), psih(:,:)
  real(8), allocatable :: xland(:,:), hfx(:,:), qfx(:,:), wspd(:,:), br(:,:)
  real(8), allocatable :: u10(:,:), v10(:,:), uoce(:,:), voce(:,:)
  real(8), allocatable :: ctopo(:,:), ctopo2(:,:)

  ! 3D outputs
  real(8), allocatable :: rublten(:,:,:), rvblten(:,:,:), rthblten(:,:,:)
  real(8), allocatable :: rqvblten(:,:,:), rqcblten(:,:,:), rqiblten(:,:,:)
  real(8), allocatable :: exch_h(:,:,:), exch_m(:,:,:)

  ! 2D outputs
  real(8), allocatable :: hpbl(:,:), wstar(:,:), delta(:,:)
  integer, allocatable :: kpbl2d(:,:)

  character(len=256) :: errmsg
  integer :: errflg
  integer :: i, k, iunit

  if (command_argument_count() >= 2) then
    call get_command_argument(1, infile)
    call get_command_argument(2, outfile)
  else
    infile = "ysu_input.bin"
    outfile = "ysu_output_fortran.bin"
  end if

  open(newunit=iunit, file=trim(infile), status='old', access='stream', form='unformatted')
  read(iunit) n_cols, nz
  read(iunit) dt
  read(iunit) ysu_topdown_pblmix, flag_qc_i, flag_qi_i
  read(iunit) cp, g, rovcp, rd, rovg, ep1, ep2, karman, xlv, rv

  flag_qc = (flag_qc_i == 1)
  flag_qi = (flag_qi_i == 1)

  allocate(u3d(n_cols, nz, 1), v3d(n_cols, nz, 1), t3d(n_cols, nz, 1), qv3d(n_cols, nz, 1))
  allocate(qc3d(n_cols, nz, 1), qi3d(n_cols, nz, 1), p3d(n_cols, nz, 1), p3di(n_cols, nz+1, 1), pi3d(n_cols, nz, 1))
  allocate(dz8w(n_cols, nz, 1), rthraten(n_cols, nz, 1))

  allocate(psfc(n_cols, 1), znt(n_cols, 1), ust(n_cols, 1), psim(n_cols, 1), psih(n_cols, 1))
  allocate(xland(n_cols, 1), hfx(n_cols, 1), qfx(n_cols, 1), wspd(n_cols, 1), br(n_cols, 1))
  allocate(u10(n_cols, 1), v10(n_cols, 1), uoce(n_cols, 1), voce(n_cols, 1))
  allocate(ctopo(n_cols, 1), ctopo2(n_cols, 1))
  ctopo = 1.0_8
  ctopo2 = 1.0_8

  allocate(rublten(n_cols, nz, 1), rvblten(n_cols, nz, 1), rthblten(n_cols, nz, 1))
  allocate(rqvblten(n_cols, nz, 1), rqcblten(n_cols, nz, 1), rqiblten(n_cols, nz, 1))
  allocate(exch_h(n_cols, nz, 1), exch_m(n_cols, nz, 1))
  allocate(hpbl(n_cols, 1), wstar(n_cols, 1), delta(n_cols, 1), kpbl2d(n_cols, 1))

  ! Read 3D arrays: element (i, k) is at i + (k-1)*n_cols in binary stream
  read(iunit) ((u3d(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((v3d(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((t3d(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((qv3d(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((qc3d(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((qi3d(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((p3d(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((p3di(i, k, 1), i=1,n_cols), k=1,nz+1)
  read(iunit) ((pi3d(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((dz8w(i, k, 1), i=1,n_cols), k=1,nz)
  read(iunit) ((rthraten(i, k, 1), i=1,n_cols), k=1,nz)

  ! Read 2D arrays
  read(iunit) (psfc(i, 1), i=1,n_cols)
  read(iunit) (znt(i, 1), i=1,n_cols)
  read(iunit) (ust(i, 1), i=1,n_cols)
  read(iunit) (psim(i, 1), i=1,n_cols)
  read(iunit) (psih(i, 1), i=1,n_cols)
  read(iunit) (xland(i, 1), i=1,n_cols)
  read(iunit) (hfx(i, 1), i=1,n_cols)
  read(iunit) (qfx(i, 1), i=1,n_cols)
  read(iunit) (wspd(i, 1), i=1,n_cols)
  read(iunit) (br(i, 1), i=1,n_cols)
  read(iunit) (u10(i, 1), i=1,n_cols)
  read(iunit) (v10(i, 1), i=1,n_cols)
  read(iunit) (uoce(i, 1), i=1,n_cols)
  read(iunit) (voce(i, 1), i=1,n_cols)

  close(iunit)

  call bl_ysu_init(errmsg, errflg)

  call ysu(u3d=u3d, v3d=v3d, t3d=t3d, qv3d=qv3d, qc3d=qc3d, qi3d=qi3d, &
           p3d=p3d, p3di=p3di, pi3d=pi3d, &
           rublten=rublten, rvblten=rvblten, rthblten=rthblten, &
           rqvblten=rqvblten, rqcblten=rqcblten, rqiblten=rqiblten, &
           flag_qc=flag_qc, flag_qi=flag_qi, &
           cp=cp, g=g, rovcp=rovcp, rd=rd, rovg=rovg, ep1=ep1, ep2=ep2, &
           karman=karman, xlv=xlv, rv=rv, &
           dz8w=dz8w, psfc=psfc, &
           znt=znt, ust=ust, hpbl=hpbl, psim=psim, psih=psih, &
           xland=xland, hfx=hfx, qfx=qfx, wspd=wspd, br=br, &
           dt=dt, kpbl2d=kpbl2d, &
           exch_h=exch_h, exch_m=exch_m, &
           wstar=wstar, delta=delta, &
           u10=u10, v10=v10, &
           uoce=uoce, voce=voce, &
           ctopo=ctopo, ctopo2=ctopo2, &
           rthraten=rthraten, ysu_topdown_pblmix=ysu_topdown_pblmix, &
           idiff=0, flag_bep=.false., &
           ids=1, ide=n_cols, jds=1, jde=1, kds=1, kde=nz, &
           ims=1, ime=n_cols, jms=1, jme=1, kms=1, kme=nz+1, &
           its=1, ite=n_cols, jts=1, jte=1, kts=1, kte=nz, &
           errmsg=errmsg, errflg=errflg)

  ! Write output binary
  open(newunit=iunit, file=trim(outfile), status='replace', access='stream', form='unformatted')
  write(iunit) ((rublten(i, k, 1), i=1,n_cols), k=1,nz)
  write(iunit) ((rvblten(i, k, 1), i=1,n_cols), k=1,nz)
  write(iunit) ((rthblten(i, k, 1), i=1,n_cols), k=1,nz)
  write(iunit) ((rqvblten(i, k, 1), i=1,n_cols), k=1,nz)
  write(iunit) ((rqcblten(i, k, 1), i=1,n_cols), k=1,nz)
  write(iunit) ((rqiblten(i, k, 1), i=1,n_cols), k=1,nz)
  write(iunit) ((exch_h(i, k, 1), i=1,n_cols), k=1,nz)
  write(iunit) ((exch_m(i, k, 1), i=1,n_cols), k=1,nz)
  write(iunit) (hpbl(i, 1), i=1,n_cols)
  write(iunit) (wstar(i, 1), i=1,n_cols)
  write(iunit) (delta(i, 1), i=1,n_cols)
  write(iunit) (real(kpbl2d(i, 1), 8), i=1,n_cols)
  close(iunit)

end program driver_ysu_fortran
