! SPDX-License-Identifier: GPL-3.0-or-later
! WRF v4.8.0's BOULAC (module_bl_boulac.o of the -r8 build) for the columns of generate_inputs.py, as the PBL driver
! calls it without BEP (flag_bep = .false., idiff = 0; the BEP arrays are passed and not read).
program driver_boulac_fortran
  use module_bl_boulac, only: boulac
  implicit none
  character(len=256) :: infile, outfile
  integer :: n, nz, i, k, iunit
  real(8) :: dt, cp, g
  real(8), allocatable, dimension(:,:,:) :: u, v, th, rho, qv, qc, dz8w, tke, dlk, wu, wv, wt, wq, exch_h, exch_m, &
       rublten, rvblten, rthblten, rqvblten, rqcblten, bep
  real(8), allocatable, dimension(:,:) :: ust, hfx, qfx, pblh, frc
  call get_command_argument(1, infile)
  call get_command_argument(2, outfile)
  open(newunit=iunit, file=trim(infile), status='old', access='stream', form='unformatted', convert='little_endian')
  read(iunit) n, nz
  read(iunit) dt, cp, g
  allocate(u(n,nz+1,1), v(n,nz+1,1), th(n,nz+1,1), rho(n,nz+1,1), qv(n,nz+1,1), qc(n,nz+1,1), dz8w(n,nz+1,1), &
           tke(n,nz+1,1), dlk(n,nz+1,1), wu(n,nz+1,1), wv(n,nz+1,1), wt(n,nz+1,1), wq(n,nz+1,1), exch_h(n,nz+1,1), &
           exch_m(n,nz+1,1), rublten(n,nz+1,1), rvblten(n,nz+1,1), rthblten(n,nz+1,1), rqvblten(n,nz+1,1), &
           rqcblten(n,nz+1,1), bep(n,nz+1,1), ust(n,1), hfx(n,1), qfx(n,1), pblh(n,1), frc(n,1))
  u = 0; v = 0; th = 0; rho = 0; qv = 0; qc = 0; dz8w = 0; tke = 0; dlk = 0; wu = 0; wv = 0; wt = 0; wq = 0
  exch_h = 0; exch_m = 0; rublten = 0; rvblten = 0; rthblten = 0; rqvblten = 0; rqcblten = 0; bep = 0
  frc = 0; pblh = 0
  read(iunit) ((u(i,k,1), i=1,n), k=1,nz)
  read(iunit) ((v(i,k,1), i=1,n), k=1,nz)
  read(iunit) ((th(i,k,1), i=1,n), k=1,nz)
  read(iunit) ((rho(i,k,1), i=1,n), k=1,nz)
  read(iunit) ((qv(i,k,1), i=1,n), k=1,nz)
  read(iunit) ((qc(i,k,1), i=1,n), k=1,nz)
  read(iunit) ((dz8w(i,k,1), i=1,n), k=1,nz)
  read(iunit) ((tke(i,k,1), i=1,n), k=1,nz)
  read(iunit) (ust(i,1), i=1,n)
  read(iunit) (hfx(i,1), i=1,n)
  read(iunit) (qfx(i,1), i=1,n)
  close(iunit)
  call boulac(frc_urb2d=frc, idiff=0, flag_bep=.false., dz8w=dz8w, dt=dt, u_phy=u, v_phy=v, th_phy=th, rho=rho, &
       qv_curr=qv, qc_curr=qc, hfx=hfx, qfx=qfx, ustar=ust, cp=cp, g=g, rublten=rublten, rvblten=rvblten, &
       rthblten=rthblten, rqvblten=rqvblten, rqcblten=rqcblten, tke=tke, dlk=dlk, wu=wu, wv=wv, wt=wt, wq=wq, &
       exch_h=exch_h, exch_m=exch_m, pblh=pblh, a_u_bep=bep, a_v_bep=bep, a_t_bep=bep, a_q_bep=bep, a_e_bep=bep, &
       b_u_bep=bep, b_v_bep=bep, b_t_bep=bep, b_q_bep=bep, b_e_bep=bep, dlg_bep=bep, dl_u_bep=bep, sf_bep=bep, &
       vl_bep=bep, ids=1, ide=n+1, jds=1, jde=2, kds=1, kde=nz+1, ims=1, ime=n, jms=1, jme=1, kms=1, kme=nz+1, &
       its=1, ite=n, jts=1, jte=1, kts=1, kte=nz)
  open(newunit=iunit, file=trim(outfile), status='replace', access='stream', form='unformatted', &
       convert='little_endian')
  write(iunit) ((tke(i,k,1), i=1,n), k=1,nz)
  write(iunit) ((dlk(i,k,1), i=1,n), k=1,nz)
  write(iunit) ((wu(i,k,1), i=1,n), k=1,nz)
  write(iunit) ((wv(i,k,1), i=1,n), k=1,nz)
  write(iunit) ((wt(i,k,1), i=1,n), k=1,nz)
  write(iunit) ((wq(i,k,1), i=1,n), k=1,nz)
  write(iunit) ((exch_h(i,k,1), i=1,n), k=1,nz)
  write(iunit) ((exch_m(i,k,1), i=1,n), k=1,nz)
  write(iunit) ((rublten(i,k,1), i=1,n), k=1,nz)
  write(iunit) ((rvblten(i,k,1), i=1,n), k=1,nz)
  write(iunit) ((rthblten(i,k,1), i=1,n), k=1,nz)
  write(iunit) ((rqvblten(i,k,1), i=1,n), k=1,nz)
  write(iunit) ((rqcblten(i,k,1), i=1,n), k=1,nz)
  write(iunit) (pblh(i,1), i=1,n)
  close(iunit)
end program driver_boulac_fortran
