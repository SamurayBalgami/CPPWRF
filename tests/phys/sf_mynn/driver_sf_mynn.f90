! WRF side of the MYNN surface-layer differential test: WRF v4.8.0's compiled mynnsfc_init / mynnsfc_driver
! (main/libwrflib.a of a -r8 build) on ni columns (j = 1), with the arguments module_surface_driver.F passes.
! usage: driver_sf_mynn_fortran IN.bin OUT.bin    (record order: fields.txt)
program driver_sf_mynn
  use module_sf_mynnsfc_driver, only: mynnsfc_init, mynnsfc_driver
  implicit none
  integer :: ni, itimestep, isfflx, nsteps, n, s
  real(8), allocatable :: in3(:,:,:,:), in2(:,:,:), io2(:,:,:)
  character(len=512) :: fin, fout, errmsg
  integer :: errflg
  integer, parameter :: n3 = 8, n2 = 9, nio = 34
  logical :: restart
  real :: xice_threshold
  call get_command_argument(1, fin)
  call get_command_argument(2, fout)
  open(10, file=fin, access='stream', form='unformatted', status='old')
  read(10) ni, itimestep, isfflx, nsteps
  allocate(in3(ni, 2, 1, n3), in2(ni, 1, n2), io2(ni, 1, nio))
  do n = 1, n3
    read(10) in3(:, 1, 1, n)
    read(10) in3(:, 2, 1, n)
  end do
  do n = 1, n2
    read(10) in2(:, 1, n)
  end do
  do n = 1, nio
    read(10) io2(:, 1, n)
  end do
  close(10)
  call mynnsfc_init(.true., errmsg, errflg)
  restart = .false.
  xice_threshold = 0.5
  do s = 1, nsteps
    call mynnsfc_driver(u3d=in3(:,:,:,1), v3d=in3(:,:,:,2), t3d=in3(:,:,:,3), qv3d=in3(:,:,:,4), &
         p3d=in3(:,:,:,5), dz8w=in3(:,:,:,6), th3d=in3(:,:,:,7), rho3d=in3(:,:,:,8), &
         psfcpa=in2(:,:,1), mavail=in2(:,:,2), xland=in2(:,:,3), tsk=in2(:,:,4), snowh=in2(:,:,5), &
         dx=in2(:,:,6), xice=in2(:,:,7), lakemask=in2(:,:,8), water_depth=in2(:,:,9), &
         pblh=io2(:,:,1), chs=io2(:,:,2), chs2=io2(:,:,3), cqs=io2(:,:,4), cqs2=io2(:,:,5), cpm=io2(:,:,6), &
         znt=io2(:,:,7), ust=io2(:,:,8), ustm=io2(:,:,9), zol=io2(:,:,10), mol=io2(:,:,11), rmol=io2(:,:,12), &
         psim=io2(:,:,13), psih=io2(:,:,14), qgh=io2(:,:,15), hfx=io2(:,:,16), qfx=io2(:,:,17), lh=io2(:,:,18), &
         flhc=io2(:,:,19), flqc=io2(:,:,20), qsfc=io2(:,:,21), u10=io2(:,:,22), v10=io2(:,:,23), &
         th2=io2(:,:,24), t2=io2(:,:,25), q2=io2(:,:,26), gz1oz0=io2(:,:,27), wspd=io2(:,:,28), br=io2(:,:,29), &
         ch=io2(:,:,30), ck=io2(:,:,31), cka=io2(:,:,32), cd=io2(:,:,33), cda=io2(:,:,34), &
         itimestep=itimestep + s - 1, spp_pbl=0, &
         ids=1, ide=ni+1, jds=1, jde=2, kds=1, kde=3, ims=1, ime=ni, jms=1, jme=1, kms=1, kme=2, &
         its=1, ite=ni, jts=1, jte=1, kts=1, kte=2, &
         xice_threshold=xice_threshold, sf_mynn_sfcflux_land=0, sf_mynn_sfcflux_water=1, &
         isfflx=isfflx, shalwater_z0=0, flag_lsm=2, restart=restart, errmsg=errmsg, errflg=errflg)
  end do
  open(11, file=fout, access='stream', form='unformatted', status='replace')
  do n = 1, nio
    write(11) io2(:, 1, n)
  end do
  close(11)
end program driver_sf_mynn
