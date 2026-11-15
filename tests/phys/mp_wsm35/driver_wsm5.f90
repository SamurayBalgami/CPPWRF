! SPDX-License-Identifier: GPL-3.0-or-later
! The WRF side of the WSM5 / WSM3 column test: wsm5init + wsm5 (or wsm3init + wsm3) of the -r8 build's
! module_mp_wsm5.o / module_mp_wsm3.o on an ni x nj slab (its..ite = 1..ni: the lanes of the vectorised i loops),
! nsteps calls in a row, with the microphysics_driver's arguments (module_model_constants) and the effective radii
! (has_reqc = has_reqi = has_reqs = 1, as with RRTMG). Output: wsm5init's module variables, then per step the state,
! RE_CLOUD / RE_ICE / RE_SNOW and the precipitation (layout of generate_inputs.py / compare.py).
program driver_wsm
  use module_mp_wsm5, only: wsm5, wsm5init, qc0, qck1, pidnc, g1pbr, g3pbr, g4pbr, g5pbro2, pvtr, pacrr, precr1, &
      precr2, roqimax, g1pbs, g3pbs, g4pbs, g5pbso2, pvts, pacrs, precs1, precs2, pidn0r, pidn0s, xlv1, pacrc, pi, &
      rslopermax, rslopesmax, rsloperbmax, rslopesbmax, rsloper2max, rslopes2max, rsloper3max, rslopes3max
  use module_mp_wsm3, only: wsm3, wsm3init
  use module_mp_wsm3, only: w3qc0 => qc0, w3qck1 => qck1, w3pidnc => pidnc, w3g1pbr => g1pbr, w3g3pbr => g3pbr, &
      w3g4pbr => g4pbr, w3g5pbro2 => g5pbro2, w3pvtr => pvtr, w3pacrr => pacrr, w3precr1 => precr1, &
      w3precr2 => precr2, w3roqimax => roqimax, w3g1pbs => g1pbs, w3g3pbs => g3pbs, w3g4pbs => g4pbs, &
      w3g5pbso2 => g5pbso2, w3pvts => pvts, w3pacrs => pacrs, w3precs1 => precs1, w3precs2 => precs2, &
      w3pidn0r => pidn0r, w3pidn0s => pidn0s, w3xlv1 => xlv1, w3pi => pi, w3rslopermax => rslopermax, &
      w3rslopesmax => rslopesmax, w3rsloperbmax => rsloperbmax, w3rslopesbmax => rslopesbmax, &
      w3rsloper2max => rsloper2max, w3rslopes2max => rslopes2max, w3rsloper3max => rsloper3max, &
      w3rslopes3max => rslopes3max
  use module_model_constants, only: g, cp, cpv, r_d, r_v, svpt0, ep_1, ep_2, epsilon, xls, xlv, xlf, rhoair0, &
      rhowater, rhosnow, cliq, cice, psat
  implicit none
  character(len=256) :: infile, outfile
  integer :: hdr(5), ni, nj, nz, nsteps, scheme, s
  real*8 :: dt
  real*8, allocatable, dimension(:,:,:) :: th, q, qc, qr, qi, qs, den, pii, p, delz, w, re_c, re_i, re_s, refl
  real*8, allocatable, dimension(:,:) :: rain, rainncv, snow, snowncv, sr
  real*8 :: c(33)
  call get_command_argument(1, infile)
  call get_command_argument(2, outfile)
  open(10, file=trim(infile), status='old', access='stream', form='unformatted')
  read(10) hdr, dt
  ni = hdr(1); nj = hdr(2); nz = hdr(3); nsteps = hdr(4); scheme = hdr(5)
  allocate(th(ni,nz,nj), q(ni,nz,nj), qc(ni,nz,nj), qr(ni,nz,nj), qi(ni,nz,nj), qs(ni,nz,nj), den(ni,nz,nj), &
           pii(ni,nz,nj), p(ni,nz,nj), delz(ni,nz,nj), w(ni,nz,nj), re_c(ni,nz,nj), re_i(ni,nz,nj), &
           re_s(ni,nz,nj), refl(ni,nz,nj), rain(ni,nj), rainncv(ni,nj), snow(ni,nj), snowncv(ni,nj), sr(ni,nj))
  read(10) th, q, qc, qr, qi, qs, den, pii, p, delz, w, rain, snow
  close(10)
  rainncv = 0.d0; snowncv = 0.d0; sr = 0.d0; re_c = 0.d0; re_i = 0.d0; re_s = 0.d0; refl = 0.d0
  if (scheme == 5) then
    call wsm5init(rhoair0, rhowater, rhosnow, cliq, cpv, .true.)
    c = (/ qc0, qck1, pidnc, g1pbr, g3pbr, g4pbr, g5pbro2, pvtr, pacrr, precr1, precr2, roqimax, g1pbs, g3pbs, &
           g4pbs, g5pbso2, pvts, pacrs, precs1, precs2, pidn0r, pidn0s, xlv1, pacrc, pi, rslopermax, rslopesmax, &
           rsloperbmax, rslopesbmax, rsloper2max, rslopes2max, rsloper3max, rslopes3max /)
  else
    call wsm3init(rhoair0, rhowater, rhosnow, cliq, cpv, .true.)
    c = (/ w3qc0, w3qck1, w3pidnc, w3g1pbr, w3g3pbr, w3g4pbr, w3g5pbro2, w3pvtr, w3pacrr, w3precr1, w3precr2, &
           w3roqimax, w3g1pbs, w3g3pbs, w3g4pbs, w3g5pbso2, w3pvts, w3pacrs, w3precs1, w3precs2, w3pidn0r, &
           w3pidn0s, w3xlv1, 0.0, w3pi, w3rslopermax, w3rslopesmax, w3rsloperbmax, w3rslopesbmax, w3rsloper2max, &
           w3rslopes2max, w3rsloper3max, w3rslopes3max /)
  end if
  open(20, file=trim(outfile), status='replace', access='stream', form='unformatted')
  write(20) c
  do s = 1, nsteps
    if (scheme == 5) then
      call wsm5(TH=th, Q=q, QC=qc, QR=qr, QI=qi, QS=qs, DEN=den, PII=pii, P=p, DELZ=delz, DELT=dt, G=g, CPD=cp, &
                CPV=cpv, RD=r_d, RV=r_v, T0C=svpt0, EP1=ep_1, EP2=ep_2, QMIN=epsilon, XLS=xls, XLV0=xlv, XLF0=xlf, &
                DEN0=rhoair0, DENR=rhowater, CLIQ=cliq, CICE=cice, PSAT=psat, RAIN=rain, RAINNCV=rainncv, &
                SNOW=snow, SNOWNCV=snowncv, SR=sr, REFL_10CM=refl, diagflag=.false., do_radar_ref=0, &
                has_reqc=1, has_reqi=1, has_reqs=1, re_cloud=re_c, re_ice=re_i, re_snow=re_s, &
                IDS=1, IDE=ni+1, JDS=1, JDE=nj+1, KDS=1, KDE=nz+1, IMS=1, IME=ni, JMS=1, JME=nj, KMS=1, KME=nz, &
                ITS=1, ITE=ni, JTS=1, JTE=nj, KTS=1, KTE=nz)
    else
      call wsm3(TH=th, Q=q, QCI=qc, QRS=qr, W=w, DEN=den, PII=pii, P=p, DELZ=delz, DELT=dt, G=g, CPD=cp, CPV=cpv, &
                RD=r_d, RV=r_v, T0C=svpt0, EP1=ep_1, EP2=ep_2, QMIN=epsilon, XLS=xls, XLV0=xlv, XLF0=xlf, &
                DEN0=rhoair0, DENR=rhowater, CLIQ=cliq, CICE=cice, PSAT=psat, RAIN=rain, RAINNCV=rainncv, &
                SNOW=snow, SNOWNCV=snowncv, SR=sr, has_reqc=1, has_reqi=1, has_reqs=1, re_cloud=re_c, &
                re_ice=re_i, re_snow=re_s, &
                IDS=1, IDE=ni+1, JDS=1, JDE=nj+1, KDS=1, KDE=nz+1, IMS=1, IME=ni, JMS=1, JME=nj, KMS=1, KME=nz, &
                ITS=1, ITE=ni, JTS=1, JTE=nj, KTS=1, KTE=nz)
    end if
    write(20) th, q, qc, qr, qi, qs, re_c, re_i, re_s, rain, rainncv, snow, snowncv, sr
  end do
  close(20)
end program driver_wsm
