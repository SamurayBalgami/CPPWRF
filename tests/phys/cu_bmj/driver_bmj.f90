! SPDX-License-Identifier: GPL-3.0-or-later
! BMJ (cu_physics = 2) column test, the WRF side: BMJINIT and BMJDRV of WRF v4.8.0's own module_cu_bmj.o (the -r8
! build, WRF_DIR), every column as one i of a row (ims:ime = 1:ncols, j = 1), the constants the cumulus driver passes.
!   driver_bmj_fortran IN OUT   (IN: generate_inputs.py; OUT: see write_out below)
program driver_bmj
  use module_cu_bmj, only : bmjinit, bmjdrv, qs0_exp, sqs_exp, ptbl_exp
  use module_model_constants, only : cp, r_d, xlv, xls, g, svpt0, ep_1
  implicit none
  integer :: ncols, nz, stepcu, ifb, k, i, u
  real :: dt
  real, allocatable, dimension(:,:,:) :: t, qv, p, pi, rho, dz8w, p8w, th, rthcuten, rqvcuten, rqccuten, rqrcuten, &
                                         ccldfra, qcconv, qiconv
  real, allocatable, dimension(:,:) :: xland, cldefi, raincv, pratec, cutop, cubot, convcld, cldefi0
  real, allocatable :: buf(:)
  integer, allocatable :: kpbl(:,:), lowlyr(:,:)
  logical, allocatable :: cu_act_flag(:,:)
  logical :: fb
  character(len=512) :: fin, fout
  integer :: ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte

  call get_command_argument(1, fin)
  call get_command_argument(2, fout)
  open(newunit=u, file=trim(fin), access='stream', form='unformatted', status='old', convert='little_endian')
  read(u) ncols, nz, stepcu, ifb
  read(u) dt
  fb = ifb /= 0
  ids = 1; ide = ncols + 1; jds = 1; jde = 2; kds = 1; kde = nz + 1
  ims = 1; ime = ncols; jms = 1; jme = 1; kms = 1; kme = nz + 1
  its = 1; ite = ncols; jts = 1; jte = 1; kts = 1; kte = nz
  allocate(t(ims:ime,kms:kme,jms:jme))
  allocate(qv, p, pi, rho, dz8w, p8w, th, rthcuten, rqvcuten, rqccuten, rqrcuten, ccldfra, qcconv, qiconv, mold=t)
  allocate(xland(ims:ime,jms:jme))
  allocate(cldefi, raincv, pratec, cutop, cubot, convcld, cldefi0, mold=xland)
  allocate(kpbl(ims:ime,jms:jme), lowlyr(ims:ime,jms:jme), cu_act_flag(ims:ime,jms:jme))
  t = 0; qv = 0; p = 0; pi = 1; rho = 0; dz8w = 0; p8w = 0; th = 0
  ccldfra = 0; qcconv = 0; qiconv = 0; convcld = 0; raincv = 0; pratec = 0; kpbl = 1
  allocate(buf(nz * ncols))
  call rd3(t, nz); call rd3(qv, nz); call rd3(p, nz); call rd3(pi, nz); call rd3(rho, nz); call rd3(dz8w, nz)
  deallocate(buf); allocate(buf((nz + 1) * ncols))
  call rd3(p8w, nz + 1)
  deallocate(buf); allocate(buf(ncols))
  read(u) buf; xland(:, 1) = buf
  read(u) buf; cldefi0(:, 1) = buf
  close(u)
  th = t / pi
  ! bmjinit (restart = .false.: CLDEFI = AVGEFI, LOWLYR = 1); the test's own CLDEFI afterwards (a later step)
  call bmjinit(rthcuten, rqvcuten, rqccuten, rqrcuten, cldefi, lowlyr, cp, r_d, .false., .true., &
               ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  cldefi = cldefi0
  call bmjdrv(ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, &
              dt, 2, stepcu, ccldfra, convcld, raincv, pratec, cutop, cubot, kpbl, &
              th, t, qv, qcconv, qiconv, fb, p8w, p, pi, rho, dz8w, &
              cp, r_d, xlv, xls, g, svpt0, ep_1, cldefi, lowlyr, xland, cu_act_flag, &
              rthcuten, rqvcuten)
  open(newunit=u, file=trim(fout), access='stream', form='unformatted', status='replace', convert='little_endian')
  write(u) qs0_exp, sqs_exp, ptbl_exp
  call wr3(rthcuten); call wr3(rqvcuten); call wr3(ccldfra); call wr3(qcconv); call wr3(qiconv)
  write(u) raincv(:, 1), pratec(:, 1), cutop(:, 1), cubot(:, 1), cldefi(:, 1), convcld(:, 1)
  close(u)
contains
  subroutine rd3(a, n)
    real, intent(inout) :: a(ims:ime, kms:kme, jms:jme)
    integer, intent(in) :: n
    read(u) buf
    do i = 1, ncols
      do k = 1, n
        a(i, k, 1) = buf((i - 1) * n + k)
      end do
    end do
  end subroutine rd3
  subroutine wr3(a)
    real, intent(in) :: a(ims:ime, kms:kme, jms:jme)
    write(u) ((a(i, k, 1), k = 1, nz), i = 1, ncols)
  end subroutine wr3
end program driver_bmj
