! Runs WRF's cu_grell_freitas_li (module_cu_gfl) on the inputs CPPWRF dumped
! (CPPWRF_GFL_DUMP) and compares RTHCUTEN, RQVCUTEN, RQCCUTEN, RQICUTEN,
! RUCUTEN, RVCUTEN, RAINCV, PRATEC, HTOP, HBOT with CPPWRF's.
! usage: harness dump.bin [ishallow=1] [ichoice=0] [dt=60] [dx=15000] [itimestep=1]
program harness
  use module_cu_gfl, only: cu_grell_freitas_li
  implicit none
  integer :: nx, ny, nz, i, j, ishallow, ichoice, itimestep, nbad, ncol
  real(8), allocatable, dimension(:,:,:) :: u, v, t, q, p, pi, rho, dz8w, w, p8w, rthften, rthraten, rthblten, &
       rqvften, rqvblten, rthcuten, rqvcuten, rqccuten, rqicuten, rucuten, rvcuten, &
       c_rth, c_rqv, c_rqc, c_rqi, c_ru, c_rv, gdc, gdc2
  real(8), allocatable, dimension(:,:) :: ht, hfx, qfx, xland, gsw, rkpbl, raincv, pratec, htop, hbot, xmbs, &
       c_rain, c_prat, c_htop, c_hbot
  integer, allocatable, dimension(:,:) :: kpbl, k22s, kbcs, ktops
  real(8) :: dt, dx, e, aerr, m, emax, amax
  character(len=256) :: fn, a

  call get_command_argument(1, fn)
  ishallow = 1; ichoice = 0; dt = 60.d0; dx = 15000.d0; itimestep = 1
  if (command_argument_count() >= 2) then; call get_command_argument(2, a); read(a,*) ishallow; endif
  if (command_argument_count() >= 3) then; call get_command_argument(3, a); read(a,*) ichoice; endif
  if (command_argument_count() >= 4) then; call get_command_argument(4, a); read(a,*) dt; endif
  if (command_argument_count() >= 5) then; call get_command_argument(5, a); read(a,*) dx; endif
  if (command_argument_count() >= 6) then; call get_command_argument(6, a); read(a,*) itimestep; endif

  open(10, file=trim(fn), access='stream', form='unformatted', status='old')
  read(10) nx, ny, nz
  allocate(u(nx+1,nz+1,ny+1)); u = 0
  v = u; t = u; q = u; p = u; pi = u; rho = u; dz8w = u; w = u; p8w = u; rthften = u; rthraten = u
  rthblten = u; rqvften = u; rqvblten = u; rthcuten = u; rqvcuten = u; rqccuten = u; rqicuten = u
  rucuten = u; rvcuten = u; c_rth = u; c_rqv = u; c_rqc = u; c_rqi = u; c_ru = u; c_rv = u
  gdc = u; gdc2 = u
  allocate(ht(nx+1,ny+1)); ht = 0
  hfx = ht; qfx = ht; xland = ht; gsw = ht; rkpbl = ht; raincv = ht; pratec = ht; htop = ht; hbot = ht; xmbs = ht
  c_rain = ht; c_prat = ht; c_htop = ht; c_hbot = ht
  allocate(kpbl(nx+1,ny+1)); kpbl = 1
  k22s = kpbl; kbcs = kpbl; ktops = kpbl

  call r3(u, nz); call r3(v, nz); call r3(t, nz); call r3(q, nz); call r3(p, nz); call r3(pi, nz); call r3(rho, nz)
  call r3(dz8w, nz); call r3(w, nz+1); call r3(p8w, nz+1)
  call r3(rthften, nz); call r3(rthraten, nz); call r3(rthblten, nz); call r3(rqvften, nz); call r3(rqvblten, nz)
  call r2(ht); call r2(hfx); call r2(qfx); call r2(xland); call r2(rkpbl)
  call r3(c_rth, nz); call r3(c_rqv, nz); call r3(c_rqc, nz); call r3(c_rqi, nz); call r3(c_ru, nz); call r3(c_rv, nz)
  call r2(c_rain); call r2(c_prat); call r2(c_htop); call r2(c_hbot)
  close(10)
  kpbl = nint(rkpbl)

  call cu_grell_freitas_li( &
       itimestep=itimestep, dt=dt, dxcell=dx, &
       u=u, v=v, w=w, t=t, q=q, rho=rho, p=p, pi=pi, p8w=p8w, dz8w=dz8w, &
       ht=ht, xland=xland, hfx=hfx, qfx=qfx, gsw=gsw, rqvften=rqvften, rthften=rthften, &
       rthblten=rthblten, rqvblten=rqvblten, rthraten=rthraten, kpbl=kpbl, &
       xlv=2.5d6, cp=1004.5d0, g=9.81d0, r_v=461.6d0, &
       ichoice_deep=ichoice, ishallow_g3=ishallow, &
       htop=htop, hbot=hbot, k22_shallow=k22s, kbcon_shallow=kbcs, ktop_shallow=ktops, &
       xmb_shallow=xmbs, raincv=raincv, pratec=pratec, gdc=gdc, gdc2=gdc2, &
       rthcuten=rthcuten, rqvcuten=rqvcuten, rqccuten=rqccuten, rqicuten=rqicuten, &
       rucuten=rucuten, rvcuten=rvcuten, &
       ims=1, ime=nx+1, jms=1, jme=ny+1, kms=1, kme=nz+1, &
       ids=1, ide=nx+1, jds=1, jde=ny+1, kds=1, kde=nz+1, &
       its=1, ite=nx, jts=1, jte=ny, kts=1, kte=nz)

  nbad = 0; ncol = 0; emax = 0.d0; amax = 0.d0
  do j = 1, ny
    do i = 1, nx
      m = maxval(abs(rthcuten(i,1:nz,j)))
      if (m > 0.d0 .or. maxval(abs(c_rth(i,1:nz,j))) > 0.d0) ncol = ncol + 1
      aerr = max( &
           maxval(abs(rthcuten(i,1:nz,j) - c_rth(i,1:nz,j))), &
           maxval(abs(rqvcuten(i,1:nz,j) - c_rqv(i,1:nz,j))), &
           maxval(abs(rqccuten(i,1:nz,j) - c_rqc(i,1:nz,j))), &
           maxval(abs(rqicuten(i,1:nz,j) - c_rqi(i,1:nz,j))), &
           maxval(abs(rucuten(i,1:nz,j)  - c_ru(i,1:nz,j))), &
           maxval(abs(rvcuten(i,1:nz,j)  - c_rv(i,1:nz,j))), &
           abs(raincv(i,j) - c_rain(i,j)), &
           abs(pratec(i,j) - c_prat(i,j)), &
           abs(htop(i,j)   - c_htop(i,j)), &
           abs(hbot(i,j)   - c_hbot(i,j)))
      e = maxval(abs(rthcuten(i,1:nz,j) - c_rth(i,1:nz,j))) / max(m, 1d-30)
      e = max(e, maxval(abs(rqvcuten(i,1:nz,j) - c_rqv(i,1:nz,j))) / max(maxval(abs(rqvcuten(i,1:nz,j))), 1d-30))
      e = max(e, maxval(abs(rqccuten(i,1:nz,j) - c_rqc(i,1:nz,j))) / max(maxval(abs(rqccuten(i,1:nz,j))), 1d-30))
      e = max(e, maxval(abs(rqicuten(i,1:nz,j) - c_rqi(i,1:nz,j))) / max(maxval(abs(rqicuten(i,1:nz,j))), 1d-30))
      e = max(e, maxval(abs(rucuten(i,1:nz,j)  - c_ru(i,1:nz,j)))  / max(maxval(abs(rucuten(i,1:nz,j))),  1d-30))
      e = max(e, maxval(abs(rvcuten(i,1:nz,j)  - c_rv(i,1:nz,j)))  / max(maxval(abs(rvcuten(i,1:nz,j))),  1d-30))
      e = max(e, abs(raincv(i,j) - c_rain(i,j)) / max(abs(raincv(i,j)), 1d-30))
      e = max(e, abs(pratec(i,j) - c_prat(i,j)) / max(abs(pratec(i,j)), 1d-30))
      e = max(e, abs(htop(i,j)   - c_htop(i,j)) / max(abs(htop(i,j)),   1d-30))
      e = max(e, abs(hbot(i,j)   - c_hbot(i,j)) / max(abs(hbot(i,j)),   1d-30))
      amax = max(amax, aerr)
      emax = max(emax, e)
      if (aerr > 0.d0) then
        nbad = nbad + 1
        if (nbad <= 15) write(*,'(a,2i4,a,es10.3,a,es10.3,a,2es14.6)') 'i,j(0-based)', i-1, j-1, &
             ' abs', aerr, ' rel', e, ' rain wrf/cpp', raincv(i,j), c_rain(i,j)
      end if
    end do
  end do
  write(*,*) 'convecting columns', ncol, ' differing >0:', nbad, ' max abs', amax, ' max rel', emax
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
