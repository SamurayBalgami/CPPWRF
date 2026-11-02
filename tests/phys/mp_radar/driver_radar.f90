! WSM6's radar reflectivity (do_radar_ref = 1): mp_wsm6_init's radar_init tables (physics_mmm/mp_radar.F90) and
! refl10cm_wsm6 per column as module_mp_wsm6.F calls it after the scheme (t1d = th*pii, refl_10cm = MAX(-35., dBZ)),
! WRF's objects of the -r8 -fdefault-double-8 build. Input: the WSM6 test's wsm6_input.bin (mp_wsm6/generate_inputs.py).
! Output (real*8, stream): xxds, xdts, xxdg, xdtg (50 each), simpson(51), k_w, pi5, lamda4, m_w_0, m_i_0
! (re, im), xcre, xcse, xcge, xcrg, xcsg, xcgg (4 each), xam_r, xbm_r, xmu_r, xobmr, xam_s, xbm_s, xmu_s, xoams,
! xobms, xocms, xam_g, xbm_g, xmu_g, xoamg, xobmg, xocmg, xorg2, xosg2, xogg2, then refl_10cm(n_cols, nz).
program driver_radar
  use mp_wsm6, only: mp_wsm6_init, refl10cm_wsm6
  use mp_radar
  implicit none

  character(len=256) :: infile, outfile, errmsg
  integer :: errflg, n_cols, nz, hail_opt, i, k
  real*8 :: c(19)
  real*8, allocatable, dimension(:,:) :: th, q, qc, qr, qi, qs, qg, den, pii, p, refl
  real*8, allocatable, dimension(:) :: t1d, p1d, qv1d, qr1d, qs1d, qg1d, dbz

  call get_command_argument(1, infile)
  call get_command_argument(2, outfile)
  open(unit=10, file=trim(infile), status='old', access='stream', form='unformatted')
  read(10) n_cols, nz, c, hail_opt
  allocate(th(n_cols, nz), q(n_cols, nz), qc(n_cols, nz), qr(n_cols, nz), qi(n_cols, nz), qs(n_cols, nz), &
           qg(n_cols, nz), den(n_cols, nz), pii(n_cols, nz), p(n_cols, nz), refl(n_cols, nz))
  allocate(t1d(nz), p1d(nz), qv1d(nz), qr1d(nz), qs1d(nz), qg1d(nz), dbz(nz))
  read(10) th, q, qc, qr, qi, qs, qg, den, pii, p
  close(10)

  ! den0, denr, dens, cliq, cpv (c: delt, g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin, xls, xlv0, xlf0, den0, denr,
  ! dens, cliq, cice, psat)
  call mp_wsm6_init(c(14), c(15), c(16), c(17), c(4), hail_opt, errmsg, errflg)
  if (errflg /= 0) then
    print *, "Error in mp_wsm6_init: ", trim(errmsg)
    stop 1
  end if

  refl = 0.
  do i = 1, n_cols
    do k = 1, nz
      t1d(k) = th(i, k)*pii(i, k)
      p1d(k) = p(i, k)
      qv1d(k) = q(i, k)
      qr1d(k) = qr(i, k)
      qs1d(k) = qs(i, k)
      qg1d(k) = qg(i, k)
    end do
    call refl10cm_wsm6(qv1d, qr1d, qs1d, qg1d, t1d, p1d, dbz, 1, nz)
    do k = 1, nz
      refl(i, k) = MAX(-35., dbz(k))
    end do
  end do

  open(unit=11, file=trim(outfile), status='replace', access='stream', form='unformatted')
  write(11) xxds, xdts, xxdg, xdtg, simpson, k_w, pi5, lamda4
  write(11) real(m_w_0, 8), aimag(m_w_0), real(m_i_0, 8), aimag(m_i_0)
  write(11) real(xcre, 8), real(xcse, 8), real(xcge, 8), real(xcrg, 8), real(xcsg, 8), real(xcgg, 8)
  write(11) real((/ xam_r, xbm_r, xmu_r, xobmr, xam_s, xbm_s, xmu_s, xoams, xobms, xocms, &
                    xam_g, xbm_g, xmu_g, xoamg, xobmg, xocmg, xorg2, xosg2, xogg2 /), 8)
  write(11) refl
  close(11)
end program driver_radar
