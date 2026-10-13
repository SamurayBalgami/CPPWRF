program driver_sfclayrev_fortran
  use module_sf_sfclayrev, only: sfclayrev
  use sf_sfclayrev, only: sf_sfclayrev_init
  implicit none

  character(len=256) :: infile, outfile
  integer :: n_cols, nz
  integer :: isfflx, shalwater_z0, isftcflx, iz0tlnd, scm_force_flux
  real(8) :: cp, g, rovcp, r, xlv, svp1, svp2, svp3, svpt0, ep1, ep2, karman, p1000mb

  ! 3D atmospheric arrays at k=1 (or 3D with nz levels)
  real(8), allocatable :: u3d(:,:,:), v3d(:,:,:), t3d(:,:,:), qv3d(:,:,:), p3d(:,:,:), dz8w(:,:,:)
  ! 2D surface arrays
  real(8), allocatable :: psfc(:,:), tsk(:,:), xland(:,:), dx(:,:), mavail(:,:), pblh(:,:)
  real(8), allocatable :: lakemask(:,:), water_depth(:,:)
  ! Inout arrays
  real(8), allocatable :: regime(:,:), hfx(:,:), qfx(:,:), qsfc(:,:), mol(:,:), rmol(:,:)
  real(8), allocatable :: gz1oz0(:,:), wspd(:,:), br(:,:), psim(:,:), psih(:,:), fm(:,:), fh(:,:)
  real(8), allocatable :: znt(:,:), zol(:,:), ust(:,:), cpm(:,:), chs2(:,:), cqs2(:,:), chs(:,:)
  real(8), allocatable :: flhc(:,:), flqc(:,:), qgh(:,:), ustm(:,:)
  ! Output arrays
  real(8), allocatable :: lh(:,:), u10(:,:), v10(:,:), th2(:,:), t2(:,:), q2(:,:)
  real(8), allocatable :: ck(:,:), cka(:,:), cd(:,:), cda(:,:)

  character(len=256) :: errmsg
  integer :: errflg
  integer :: i, iunit

  if (command_argument_count() >= 2) then
    call get_command_argument(1, infile)
    call get_command_argument(2, outfile)
  else
    infile = "sfclayrev_input.bin"
    outfile = "sfclayrev_output_fortran.bin"
  end if

  open(newunit=iunit, file=trim(infile), status='old', access='stream', form='unformatted')
  read(iunit) n_cols, nz
  read(iunit) isfflx, shalwater_z0, isftcflx, iz0tlnd, scm_force_flux
  read(iunit) cp, g, rovcp, r, xlv, svp1, svp2, svp3, svpt0, ep1, ep2, karman, p1000mb

  ! Allocate arrays: dim1 is 1:n_cols, dim2 is 1:nz for 3D, dim2 is 1:1 for 2D (WRF j-loop)
  allocate(u3d(n_cols, nz, 1), v3d(n_cols, nz, 1), t3d(n_cols, nz, 1))
  allocate(qv3d(n_cols, nz, 1), p3d(n_cols, nz, 1), dz8w(n_cols, nz, 1))
  allocate(psfc(n_cols, 1), tsk(n_cols, 1), xland(n_cols, 1), dx(n_cols, 1))
  allocate(mavail(n_cols, 1), pblh(n_cols, 1), lakemask(n_cols, 1), water_depth(n_cols, 1))

  allocate(regime(n_cols, 1), hfx(n_cols, 1), qfx(n_cols, 1), qsfc(n_cols, 1))
  allocate(mol(n_cols, 1), rmol(n_cols, 1), gz1oz0(n_cols, 1), wspd(n_cols, 1), br(n_cols, 1))
  allocate(psim(n_cols, 1), psih(n_cols, 1), fm(n_cols, 1), fh(n_cols, 1))
  allocate(znt(n_cols, 1), zol(n_cols, 1), ust(n_cols, 1), cpm(n_cols, 1))
  allocate(chs2(n_cols, 1), cqs2(n_cols, 1), chs(n_cols, 1), flhc(n_cols, 1))
  allocate(flqc(n_cols, 1), qgh(n_cols, 1), ustm(n_cols, 1))

  allocate(lh(n_cols, 1), u10(n_cols, 1), v10(n_cols, 1), th2(n_cols, 1), t2(n_cols, 1), q2(n_cols, 1))
  allocate(ck(n_cols, 1), cka(n_cols, 1), cd(n_cols, 1), cda(n_cols, 1))

  ! Read input arrays (level 1 / k=1 slice stored in binary)
  read(iunit) (u3d(i, 1, 1), i=1,n_cols)
  read(iunit) (v3d(i, 1, 1), i=1,n_cols)
  read(iunit) (t3d(i, 1, 1), i=1,n_cols)
  read(iunit) (qv3d(i, 1, 1), i=1,n_cols)
  read(iunit) (p3d(i, 1, 1), i=1,n_cols)
  read(iunit) (dz8w(i, 1, 1), i=1,n_cols)
  read(iunit) (psfc(i, 1), i=1,n_cols)
  read(iunit) (tsk(i, 1), i=1,n_cols)
  read(iunit) (xland(i, 1), i=1,n_cols)
  read(iunit) (dx(i, 1), i=1,n_cols)
  read(iunit) (mavail(i, 1), i=1,n_cols)
  read(iunit) (pblh(i, 1), i=1,n_cols)
  read(iunit) (lakemask(i, 1), i=1,n_cols)
  read(iunit) (water_depth(i, 1), i=1,n_cols)

  ! Read initial inout arrays
  read(iunit) (regime(i, 1), i=1,n_cols)
  read(iunit) (hfx(i, 1), i=1,n_cols)
  read(iunit) (qfx(i, 1), i=1,n_cols)
  read(iunit) (qsfc(i, 1), i=1,n_cols)
  read(iunit) (mol(i, 1), i=1,n_cols)
  read(iunit) (rmol(i, 1), i=1,n_cols)
  read(iunit) (gz1oz0(i, 1), i=1,n_cols)
  read(iunit) (wspd(i, 1), i=1,n_cols)
  read(iunit) (br(i, 1), i=1,n_cols)
  read(iunit) (psim(i, 1), i=1,n_cols)
  read(iunit) (psih(i, 1), i=1,n_cols)
  read(iunit) (fm(i, 1), i=1,n_cols)
  read(iunit) (fh(i, 1), i=1,n_cols)
  read(iunit) (znt(i, 1), i=1,n_cols)
  read(iunit) (zol(i, 1), i=1,n_cols)
  read(iunit) (ust(i, 1), i=1,n_cols)
  read(iunit) (cpm(i, 1), i=1,n_cols)
  read(iunit) (chs2(i, 1), i=1,n_cols)
  read(iunit) (cqs2(i, 1), i=1,n_cols)
  read(iunit) (chs(i, 1), i=1,n_cols)
  read(iunit) (flhc(i, 1), i=1,n_cols)
  read(iunit) (flqc(i, 1), i=1,n_cols)
  read(iunit) (qgh(i, 1), i=1,n_cols)
  read(iunit) (ustm(i, 1), i=1,n_cols)

  close(iunit)

  ! Initialize lookup tables
  call sf_sfclayrev_init(errmsg, errflg)

  ! Call WRF sfclayrev subroutine
  call sfclayrev(u3d=u3d, v3d=v3d, t3d=t3d, qv3d=qv3d, p3d=p3d, dz8w=dz8w, &
                 cp=cp, g=g, rovcp=rovcp, r=r, xlv=xlv, psfc=psfc, &
                 chs=chs, chs2=chs2, cqs2=cqs2, cpm=cpm, &
                 znt=znt, ust=ust, pblh=pblh, mavail=mavail, zol=zol, mol=mol, regime=regime, &
                 psim=psim, psih=psih, fm=fm, fh=fh, &
                 xland=xland, hfx=hfx, qfx=qfx, lh=lh, tsk=tsk, flhc=flhc, flqc=flqc, &
                 qgh=qgh, qsfc=qsfc, rmol=rmol, &
                 u10=u10, v10=v10, th2=th2, t2=t2, q2=q2, &
                 gz1oz0=gz1oz0, wspd=wspd, br=br, isfflx=isfflx, dx=dx, &
                 svp1=svp1, svp2=svp2, svp3=svp3, svpt0=svpt0, ep1=ep1, ep2=ep2, &
                 karman=karman, p1000mb=p1000mb, lakemask=lakemask, &
                 ids=1, ide=n_cols, jds=1, jde=1, kds=1, kde=nz, &
                 ims=1, ime=n_cols, jms=1, jme=1, kms=1, kme=nz, &
                 its=1, ite=n_cols, jts=1, jte=1, kts=1, kte=nz, &
                 ustm=ustm, ck=ck, cka=cka, cd=cd, cda=cda, &
                 isftcflx=isftcflx, iz0tlnd=iz0tlnd, &
                 shalwater_z0=shalwater_z0, water_depth=water_depth, &
                 scm_force_flux=scm_force_flux, errmsg=errmsg, errflg=errflg)

  ! Write output binary
  open(newunit=iunit, file=trim(outfile), status='replace', access='stream', form='unformatted')
  write(iunit) (regime(i, 1), i=1,n_cols)
  write(iunit) (hfx(i, 1), i=1,n_cols)
  write(iunit) (qfx(i, 1), i=1,n_cols)
  write(iunit) (lh(i, 1), i=1,n_cols)
  write(iunit) (qsfc(i, 1), i=1,n_cols)
  write(iunit) (mol(i, 1), i=1,n_cols)
  write(iunit) (rmol(i, 1), i=1,n_cols)
  write(iunit) (gz1oz0(i, 1), i=1,n_cols)
  write(iunit) (wspd(i, 1), i=1,n_cols)
  write(iunit) (br(i, 1), i=1,n_cols)
  write(iunit) (psim(i, 1), i=1,n_cols)
  write(iunit) (psih(i, 1), i=1,n_cols)
  write(iunit) (fm(i, 1), i=1,n_cols)
  write(iunit) (fh(i, 1), i=1,n_cols)
  write(iunit) (znt(i, 1), i=1,n_cols)
  write(iunit) (zol(i, 1), i=1,n_cols)
  write(iunit) (ust(i, 1), i=1,n_cols)
  write(iunit) (cpm(i, 1), i=1,n_cols)
  write(iunit) (chs2(i, 1), i=1,n_cols)
  write(iunit) (cqs2(i, 1), i=1,n_cols)
  write(iunit) (chs(i, 1), i=1,n_cols)
  write(iunit) (flhc(i, 1), i=1,n_cols)
  write(iunit) (flqc(i, 1), i=1,n_cols)
  write(iunit) (qgh(i, 1), i=1,n_cols)
  write(iunit) (ustm(i, 1), i=1,n_cols)
  write(iunit) (u10(i, 1), i=1,n_cols)
  write(iunit) (v10(i, 1), i=1,n_cols)
  write(iunit) (th2(i, 1), i=1,n_cols)
  write(iunit) (t2(i, 1), i=1,n_cols)
  write(iunit) (q2(i, 1), i=1,n_cols)
  write(iunit) (ck(i, 1), i=1,n_cols)
  write(iunit) (cka(i, 1), i=1,n_cols)
  write(iunit) (cd(i, 1), i=1,n_cols)
  write(iunit) (cda(i, 1), i=1,n_cols)
  close(iunit)

end program driver_sfclayrev_fortran
