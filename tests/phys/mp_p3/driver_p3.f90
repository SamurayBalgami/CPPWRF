! SPDX-License-Identifier: GPL-3.0-or-later
! WRF v4.8.0's P3 (module_mp_p3.o of the -r8 -fdefault-double-8 build): p3_init('.', 1, .false.) and
! mp_p3_wrapper_wrf on one (i, k) slab (j = 1) for nsteps calls, as module_microphysics_driver.F calls it for
! mp_physics = 50 (no nc_3d) and 51 (nc_3d = qnc). Writes the state and diagnostics after every call.
program driver_p3
  use microphy_p3, only: p3_init, mp_p3_wrapper_wrf
  implicit none
  character(len=256) :: infile, outfile
  integer :: n, nz, nsteps, predict, start_it, istat, s
  real*8 :: dt
  real*8, allocatable, dimension(:,:,:) :: th, qv, qc, qr, qnr, qi, qni, qir, qib, qnc, tho, qvo, p, dz, w, pii, &
                                           refl, effc, effi, vmi, di, rhopo
  real*8, allocatable, dimension(:,:) :: rainnc, rainncv, sr, snownc, snowncv

  call get_command_argument(1, infile)
  call get_command_argument(2, outfile)
  open(10, file=trim(infile), access='stream', form='unformatted', status='old', convert='little_endian')
  read(10) n, nz, nsteps, predict, dt
  read(10) start_it
  allocate(th(n,nz,1), qv(n,nz,1), qc(n,nz,1), qr(n,nz,1), qnr(n,nz,1), qi(n,nz,1), qni(n,nz,1), qir(n,nz,1), &
           qib(n,nz,1), qnc(n,nz,1), tho(n,nz,1), qvo(n,nz,1), p(n,nz,1), dz(n,nz,1), w(n,nz,1), pii(n,nz,1), &
           refl(n,nz,1), effc(n,nz,1), effi(n,nz,1), vmi(n,nz,1), di(n,nz,1), rhopo(n,nz,1))
  allocate(rainnc(n,1), rainncv(n,1), sr(n,1), snownc(n,1), snowncv(n,1))
  read(10) th; read(10) qv; read(10) qc; read(10) qr; read(10) qnr; read(10) qi; read(10) qni; read(10) qir
  read(10) qib; read(10) qnc; read(10) tho; read(10) qvo; read(10) p; read(10) dz
  close(10)
  w = 0.; pii = 0.; refl = 0.; effc = 0.; effi = 0.; vmi = 0.; di = 0.; rhopo = 0.
  rainnc = 0.; rainncv = 0.; sr = 0.; snownc = 0.; snowncv = 0.

  call p3_init('.', 1, .false., 'WRF', istat, .false.)

  open(11, file=trim(outfile), access='stream', form='unformatted', status='replace', convert='little_endian')
  do s = 1, nsteps
    if (predict == 1) then
      call mp_p3_wrapper_wrf(ITIMESTEP=start_it+s-1, TH_3d=th, QV_3d=qv, QC_3d=qc, QR_3d=qr, QNR_3d=qnr,       &
             QI1_3d=qi, QIR1_3d=qir, QNI1_3d=qni, QIB1_3d=qib, th_old_3d=tho, qv_old_3d=qvo, nc_3d=qnc,        &
             PII=pii, P=p, DT=dt, DZ=dz, W=w, RAINNC=rainnc, RAINNCV=rainncv, SR=sr, SNOWNC=snownc,             &
             SNOWNCV=snowncv, N_ICECAT=1,                                                                      &
             IDS=1, IDE=n+1, JDS=1, JDE=2, KDS=1, KDE=nz+1, IMS=1, IME=n, JMS=1, JME=1, KMS=1, KME=nz,          &
             ITS=1, ITE=n, JTS=1, JTE=1, KTS=1, KTE=nz, diag_zdbz_3d=refl, diag_effc_3d=effc, diag_effi_3d=effi, &
             diag_vmi_3d=vmi, diag_di_3d=di, diag_rhopo_3d=rhopo)
    else
      call mp_p3_wrapper_wrf(ITIMESTEP=start_it+s-1, TH_3d=th, QV_3d=qv, QC_3d=qc, QR_3d=qr, QNR_3d=qnr,       &
             QI1_3d=qi, QIR1_3d=qir, QNI1_3d=qni, QIB1_3d=qib, th_old_3d=tho, qv_old_3d=qvo,                   &
             PII=pii, P=p, DT=dt, DZ=dz, W=w, RAINNC=rainnc, RAINNCV=rainncv, SR=sr, SNOWNC=snownc,             &
             SNOWNCV=snowncv, N_ICECAT=1,                                                                      &
             IDS=1, IDE=n+1, JDS=1, JDE=2, KDS=1, KDE=nz+1, IMS=1, IME=n, JMS=1, JME=1, KMS=1, KME=nz,          &
             ITS=1, ITE=n, JTS=1, JTE=1, KTS=1, KTE=nz, diag_zdbz_3d=refl, diag_effc_3d=effc, diag_effi_3d=effi, &
             diag_vmi_3d=vmi, diag_di_3d=di, diag_rhopo_3d=rhopo)
    endif
    write(11) th, qv, qc, qr, qnr, qi, qni, qir, qib, qnc, tho, qvo, refl, effc, effi, vmi, di, rhopo
    write(11) rainnc, rainncv, sr, snownc, snowncv
  enddo
  close(11)
end program driver_p3
