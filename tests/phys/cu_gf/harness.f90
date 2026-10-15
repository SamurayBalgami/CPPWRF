! Runs WRF's GFDRV (module_cu_gf_wrfdrv) on the inputs CPPWRF dumped (CPPWRF_GF_DUMP) and compares
! RTHCUTEN, RQVCUTEN, RUCUTEN, RAINCV with CPPWRF's.  usage: harness dump.bin [ishallow] [ichoice]
program harness
  use module_cu_gf_wrfdrv, only: gfdrv
  implicit none
  integer :: nx, ny, nz, i, j, k, ios, ishallow, ichoice, nbad, ncol
  real(8), allocatable, dimension(:,:,:) :: u, v, t, q, p, pi, rho, dz8w, w, p8w, rthften, rthraten, rthblten, &
       rqvften, rqvblten, rthcuten, rqvcuten, rqccuten, rqicuten, rucuten, rvcuten, c_rth, c_rqv, c_ru, c_rv
  real(8), allocatable, dimension(:,:) :: ht, hfx, qfx, xland, rkpbl, raincv, pratec, htop, hbot, c_rain, c_prat
  integer, allocatable :: kpbl(:,:), ktop_deep(:,:), k22s(:,:), kbcs(:,:), ktops(:,:)
  real(8), allocatable :: xmbs(:,:), gdc(:,:,:), gdc2(:,:,:)
  real(8) :: dt, dx, e, m, emax
  character(len=256) :: fn, a
  call get_command_argument(1, fn)
  ishallow = 1; ichoice = 0
  if (command_argument_count() >= 2) then; call get_command_argument(2, a); read(a,*) ishallow; endif
  if (command_argument_count() >= 3) then; call get_command_argument(3, a); read(a,*) ichoice; endif
  open(10, file=trim(fn), access='stream', form='unformatted', status='old')
  read(10) nx, ny, nz
  allocate(u(nx+1,nz+1,ny+1)); u = 0
  v = u; t = u; q = u; p = u; pi = u; rho = u; dz8w = u; w = u; p8w = u; rthften = u; rthraten = u
  rthblten = u; rqvften = u; rqvblten = u; rthcuten = u; rqvcuten = u; rqccuten = u; rqicuten = u; rucuten = u
  rvcuten = u; c_rth = u; c_rqv = u; c_ru = u; c_rv = u
  allocate(ht(nx+1,ny+1)); ht = 0
  hfx = ht; qfx = ht; xland = ht; rkpbl = ht; raincv = ht; pratec = ht; htop = ht; hbot = ht; c_rain = ht; c_prat = ht
  allocate(kpbl(nx+1,ny+1)); kpbl = 1
  ktop_deep = kpbl; k22s = kpbl; kbcs = kpbl; ktops = kpbl; xmbs = ht; gdc = u; gdc2 = u
  call r3(u, nz); call r3(v, nz); call r3(t, nz); call r3(q, nz); call r3(p, nz); call r3(pi, nz); call r3(rho, nz)
  call r3(dz8w, nz); call r3(w, nz+1); call r3(p8w, nz+1)
  call r3(rthften, nz); call r3(rthraten, nz); call r3(rthblten, nz); call r3(rqvften, nz); call r3(rqvblten, nz)
  call r2(ht); call r2(hfx); call r2(qfx); call r2(xland); call r2(rkpbl)
  call r3(c_rth, nz); call r3(c_rqv, nz); call r3(c_ru, nz); call r3(c_rv, nz); call r2(c_rain); call r2(c_prat)
  kpbl = nint(rkpbl)
  dt = 60.d0; dx = 15000.d0
  if (command_argument_count() >= 4) then; call get_command_argument(4, a); read(a,*) dt; endif
  if (command_argument_count() >= 5) then; call get_command_argument(5, a); read(a,*) dx; endif
  call gfdrv(spp_conv=0, dt=dt, dx=dx, rho=rho, raincv=raincv, pratec=pratec, u=u, v=v, t=t, w=w, q=q, p=p, pi=pi, &
       dz8w=dz8w, p8w=p8w, htop=htop, hbot=hbot, ktop_deep=ktop_deep, gdc=gdc, gdc2=gdc2, k22_shallow=k22s, kbcon_shallow=kbcs, ktop_shallow=ktops, xmb_shallow=xmbs, &
       ht=ht, hfx=hfx, qfx=qfx, xland=xland, kpbl=kpbl, &
       ichoice=ichoice, ishallow_g3=ishallow, &
       ids=1, ide=nx+1, jds=1, jde=ny+1, kds=1, kde=nz+1, ims=1, ime=nx+1, jms=1, jme=ny+1, kms=1, kme=nz+1, &
       its=1, ite=nx, jts=1, jte=ny, kts=1, kte=nz, periodic_x=.false., periodic_y=.false., &
       rqvcuten=rqvcuten, rqccuten=rqccuten, rqicuten=rqicuten, rqvften=rqvften, rthften=rthften, &
       rthcuten=rthcuten, rthraten=rthraten, rqvblten=rqvblten, rthblten=rthblten, dudt_phy=rucuten, dvdt_phy=rvcuten)
  nbad = 0; ncol = 0; emax = 0
  do j = 1, ny
    do i = 1, nx
      m = maxval(abs(rthcuten(i,1:nz,j)))
      if (m > 0 .or. maxval(abs(c_rth(i,1:nz,j))) > 0) ncol = ncol + 1
      e = maxval(abs(rthcuten(i,1:nz,j) - c_rth(i,1:nz,j))) / max(m, 1d-30)
      e = max(e, maxval(abs(rqvcuten(i,1:nz,j) - c_rqv(i,1:nz,j))) / max(maxval(abs(rqvcuten(i,1:nz,j))), 1d-30))
      e = max(e, maxval(abs(rucuten(i,1:nz,j) - c_ru(i,1:nz,j))) / max(maxval(abs(rucuten(i,1:nz,j))), 1d-30))
      e = max(e, abs(raincv(i,j) - c_rain(i,j)) / max(abs(raincv(i,j)), 1d-30))
      emax = max(emax, e)
      if (e > 1d-12) then
        nbad = nbad + 1
        if (nbad <= 15) write(*,'(a,2i4,a,es10.3,a,2es14.6)') 'i,j(0-based)', i-1, j-1, ' rel', e, &
             ' rain wrf/cpp', raincv(i,j), c_rain(i,j)
      end if
    end do
  end do
  write(*,*) 'convecting columns', ncol, ' differing >1e-12:', nbad, ' max rel', emax
contains
  subroutine r3(x, nk)
    real(8) :: x(:,:,:)
    integer :: nk, ii, jj, kk
    do jj = 1, ny; do kk = 1, nk; do ii = 1, nx
      read(10) x(ii,kk,jj)
    end do; end do; end do
  end subroutine
  subroutine r2(x)
    real(8) :: x(:,:)
    integer :: ii, jj
    do jj = 1, ny; do ii = 1, nx
      read(10) x(ii,jj)
    end do; end do
  end subroutine
end program
