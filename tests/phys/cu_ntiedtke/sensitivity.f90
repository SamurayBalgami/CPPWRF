program sens3
  use cu_ntiedtke, only: cu_ntiedtke_run, cu_ntiedtke_init
  implicit none
  integer, parameter :: km = 44
  real(8) :: sc(6), a(km, 11), hi(km+1, 2), outc(km, 6), rnc, a0(km,11), hi0(km+1,2), r1(km,11), r2(km+1,2)
  real(8), dimension(1,km) :: u, v, t, qv, qc, qi, qvf, tf, oz, om, ap
  real(8), dimension(1,km+1) :: zz, aph
  real(8) :: evap(1), hfx(1), dx(1), rn(1), pert, d, dmax
  integer :: lndj(1), col, errflg, trial, target, kmx
  character(len=256) :: fn, errmsg, ts
  real(8), parameter :: rd = 287., cp = 7.*rd/2.
  call get_command_argument(1, fn); call get_command_argument(2, ts); read(ts,*) target
  call get_command_argument(3, ts); read(ts,*) pert
  call cu_ntiedtke_init(cp, rd, 461.6d0, 2.5d6, 2.85d6, 3.5d5, 9.81d0, errmsg, errflg)
  open(10, file=trim(fn), access='stream', form='unformatted', status='old')
  do col = 1, target
    read(10) sc, a0, hi0, outc, rnc
  end do
  dmax = 0
  do trial = 1, 500
    call random_number(r1); call random_number(r2)
    a = a0 * (1d0 + pert*(2*r1-1)); hi = hi0 * (1d0 + pert*(2*r2-1))
    u(1,:) = a(:,1); v(1,:) = a(:,2); t(1,:) = a(:,3); qv(1,:) = a(:,4); qc(1,:) = a(:,5); qi(1,:) = a(:,6)
    qvf(1,:) = a(:,7); tf(1,:) = a(:,8); oz(1,:) = a(:,9); om(1,:) = a(:,10); ap(1,:) = a(:,11)
    zz(1,:) = hi(:,1); aph(1,:) = hi(:,2)
    evap(1) = sc(1); hfx(1) = sc(2); lndj(1) = nint(sc(3)); dx(1) = sc(5); rn = 0
    call cu_ntiedtke_run(u, v, t, qv, qc, qi, qvf, tf, oz, zz, om, ap, aph, evap, hfx, rn, lndj, 1, km, km+1, &
                         sc(4), dx, errmsg, errflg)
    d = maxval(abs((t(1,:)-a(:,3)) - (outc(:,3)-a0(:,3))))
    if (d > dmax) then; dmax = d; kmx = maxloc(abs((t(1,:)-a(:,3)) - (outc(:,3)-a0(:,3))),1); end if
  end do
  write(*,*) 'max |dT_conv - base| over trials', dmax, ' top-index', kmx
end program
