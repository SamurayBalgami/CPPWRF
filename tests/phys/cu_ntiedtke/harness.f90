program harn
  use ccpp_kind_types, only: kind_phys
  use cu_ntiedtke, only: cu_ntiedtke_run, cu_ntiedtke_init
  implicit none
  integer, parameter :: km = 44
  real(8) :: sc(6), a(km, 11), hi(km+1, 2), outc(km, 6), rnc
  real(8), dimension(1,km) :: u, v, t, qv, qc, qi, qvf, tf, oz, om, ap
  real(8), dimension(1,km+1) :: zz, aph
  real(8) :: evap(1), hfx(1), dx(1), rn(1), d, dmax(7), rel
  integer :: lndj(1), ios, col, nbad, errflg, k
  character(len=256) :: fn, errmsg
  real(8), parameter :: rd = 287., cp = 7.*rd/2.
  call get_command_argument(1, fn)
  call cu_ntiedtke_init(cp, rd, 461.6d0, 2.5d6, 2.85d6, 3.5d5, 9.81d0, errmsg, errflg)
  open(10, file=trim(fn), access='stream', form='unformatted', status='old')
  col = 0; nbad = 0; dmax = 0
  do
    read(10, iostat=ios) sc
    if (ios /= 0) exit
    read(10) a, hi, outc, rnc
    col = col + 1
    u(1,:) = a(:,1); v(1,:) = a(:,2); t(1,:) = a(:,3); qv(1,:) = a(:,4); qc(1,:) = a(:,5); qi(1,:) = a(:,6)
    qvf(1,:) = a(:,7); tf(1,:) = a(:,8); oz(1,:) = a(:,9); om(1,:) = a(:,10); ap(1,:) = a(:,11)
    zz(1,:) = hi(:,1); aph(1,:) = hi(:,2)
    evap(1) = sc(1); hfx(1) = sc(2); lndj(1) = nint(sc(3)); dx(1) = sc(5)
    rn = 0
    call cu_ntiedtke_run(u, v, t, qv, qc, qi, qvf, tf, oz, zz, om, ap, aph, evap, hfx, rn, lndj, 1, km, km+1, &
                         sc(4), dx, errmsg, errflg)
    rel = 0
    rel = max(rel, maxval(abs(u(1,:) - outc(:,1))) / max(1d-30, maxval(abs(u(1,:)))))
    rel = max(rel, maxval(abs(v(1,:) - outc(:,2))) / max(1d-30, maxval(abs(v(1,:)))))
    rel = max(rel, maxval(abs(t(1,:) - outc(:,3))) / max(1d-30, maxval(abs(t(1,:)))))
    rel = max(rel, maxval(abs(qv(1,:) - outc(:,4))) / max(1d-30, maxval(abs(qv(1,:)))))
    rel = max(rel, abs(rn(1) - rnc) / max(1d-30, abs(rn(1))))
    if (rel > 0d0) then
      nbad = nbad + 1
      if (nbad <= 10) write(*,'(a,i6,a,es10.3,a,2es14.6,a,2es14.6)') 'col', col, ' rel', rel, ' rn f/c', rn(1), rnc, &
        ' tbot f/c', t(1,km), outc(km,3)
    end if
  end do
  write(*,*) 'columns', col, ' mismatching', nbad
end program
