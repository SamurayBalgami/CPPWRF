! WRF side of the MYNN-EDMF differential test: WRF v4.8.0's compiled mynnedmf_driver (main/libwrflib.a of a
! -r8 -fdefault-double-8 build) on ni columns (j = 1) with WRF's default bl_mynn_* options, for nsteps calls
! (initflag = 1 at the first); between calls u, v, theta, qv, qc, qi take the tendencies (x = x + dx/dt * dt).
! usage: driver_bl_mynn_fortran IN.bin OUT.bin    (record order: generate_inputs.py, compare.py)
program driver_bl_mynn
  use module_bl_mynnedmf_driver, only: mynnedmf_driver
  implicit none
  integer :: ni, nz, nsteps, s, n, initflag
  real(8) :: delt, dx
  integer, parameter :: n3 = 14, n2 = 12, no3 = 25, no2 = 7
  real(8), allocatable :: in3(:,:,:,:), in2(:,:,:), o3(:,:,:,:), o2(:,:,:), zero3(:,:,:)
  real(8), allocatable :: rqsblten(:,:,:), rqncblten(:,:,:), rqniblten(:,:,:), rqnwfablten(:,:,:), &
       rqnifablten(:,:,:), rqnbcablten(:,:,:), qnc(:,:,:), qni(:,:,:), qnwfa(:,:,:), qnifa(:,:,:), qnbca(:,:,:), &
       sub_thl3d(:,:,:), sub_sqv3d(:,:,:), det_thl3d(:,:,:), det_sqv3d(:,:,:), dqke(:,:,:), qwt(:,:,:), &
       qshear(:,:,:), qbuoy(:,:,:), qdiss(:,:,:), maxwidth_dd(:,:), maxmf_dd(:,:), maxtkeprod(:,:), &
       cldtop_cooling(:,:), ent_eff(:,:)
  integer, allocatable :: kpbl(:,:)
  character(len=512) :: fin, fout
  call get_command_argument(1, fin)
  call get_command_argument(2, fout)
  open(10, file=fin, access='stream', form='unformatted', status='old')
  read(10) ni, nz, nsteps
  read(10) delt, dx
  allocate(in3(ni, nz, 1, n3), in2(ni, 1, n2), o3(ni, nz, 1, no3), o2(ni, 1, no2), zero3(ni, nz, 1))
  do n = 1, n3
    read(10) in3(:, :, 1, n)
  end do
  do n = 1, n2
    read(10) in2(:, 1, n)
  end do
  close(10)
  o3 = 0.d0
  o2 = 0.d0
  zero3 = 0.d0
  allocate(rqsblten, rqncblten, rqniblten, rqnwfablten, rqnifablten, rqnbcablten, qnc, qni, qnwfa, qnifa, qnbca, &
           sub_thl3d, sub_sqv3d, det_thl3d, det_sqv3d, dqke, qwt, qshear, qbuoy, qdiss, mold=zero3)
  allocate(maxwidth_dd(ni, 1), maxmf_dd(ni, 1), maxtkeprod(ni, 1), cldtop_cooling(ni, 1), ent_eff(ni, 1), kpbl(ni, 1))
  qnc = 0.d0; qni = 0.d0; qnwfa = 0.d0; qnifa = 0.d0; qnbca = 0.d0
  kpbl = 0
  do s = 1, nsteps
    initflag = merge(1, 0, s == 1)
    call mynnedmf_driver(ids=1, ide=ni+1, jds=1, jde=2, kds=1, kde=nz+1, ims=1, ime=ni, jms=1, jme=1, kms=1, kme=nz, &
         its=1, ite=ni, jts=1, jte=1, kts=1, kte=nz, &
         flag_qc=.true., flag_qi=.true., flag_qs=.true., flag_qnc=.false., flag_qni=.false., &
         flag_qnifa=.false., flag_qnwfa=.false., flag_qnbca=.false., initflag=initflag, &
         restart=.false., cycling=.false., delt=delt, dxc=dx, xland=in2(:,:,1), ps=in2(:,:,4), ts=in2(:,:,2), &
         qsfc=in2(:,:,3), ust=in2(:,:,5), ch=in2(:,:,6), hfx=in2(:,:,7), qfx=in2(:,:,8), wspd=in2(:,:,9), &
         znt=in2(:,:,10), uoce=in2(:,:,11), voce=in2(:,:,12), &
         dz=in3(:,:,:,1), u=in3(:,:,:,2), v=in3(:,:,:,3), w=in3(:,:,:,4), th=in3(:,:,:,5), t3d=in3(:,:,:,6), &
         p=in3(:,:,:,7), exner=in3(:,:,:,8), rho=in3(:,:,:,9), qv=in3(:,:,:,10), qc=in3(:,:,:,11), &
         qi=in3(:,:,:,12), qs=in3(:,:,:,13), qnc=qnc, qni=qni, qnifa=qnifa, qnwfa=qnwfa, qnbca=qnbca, &
         rthraten=in3(:,:,:,14), pblh=o2(:,:,1), kpbl=kpbl, maxwidth_dd=maxwidth_dd, &
         cldfra_bl=o3(:,:,:,19), qc_bl=o3(:,:,:,17), qi_bl=o3(:,:,:,18), maxwidth=o2(:,:,2), maxmf=o2(:,:,3), &
         ztop_plume=o2(:,:,4), excess_h=o2(:,:,5), excess_q=o2(:,:,6), maxmf_dd=maxmf_dd, &
         maxtkeprod=maxtkeprod, cldtop_cooling=cldtop_cooling, ent_eff=ent_eff, &
         qke=o3(:,:,:,1), qke_adv=o3(:,:,:,10), tsq=o3(:,:,:,7), qsq=o3(:,:,:,8), cov=o3(:,:,:,9), &
         el_pbl=o3(:,:,:,2), rublten=o3(:,:,:,11), rvblten=o3(:,:,:,12), rthblten=o3(:,:,:,13), &
         rqvblten=o3(:,:,:,14), rqcblten=o3(:,:,:,15), rqiblten=o3(:,:,:,16), rqsblten=rqsblten, &
         rqncblten=rqncblten, rqniblten=rqniblten, rqnifablten=rqnifablten, rqnwfablten=rqnwfablten, &
         rqnbcablten=rqnbcablten, edmf_a=o3(:,:,:,20), edmf_w=o3(:,:,:,21), edmf_qt=o3(:,:,:,22), &
         edmf_thl=o3(:,:,:,23), edmf_ent=o3(:,:,:,24), edmf_qc=o3(:,:,:,25), &
         sub_thl3d=sub_thl3d, sub_sqv3d=sub_sqv3d, det_thl3d=det_thl3d, det_sqv3d=det_sqv3d, &
         exch_h=o3(:,:,:,5), exch_m=o3(:,:,:,6), dqke=dqke, qwt=qwt, qshear=qshear, qbuoy=qbuoy, qdiss=qdiss, &
         sh3d=o3(:,:,:,3), sm3d=o3(:,:,:,4), spp_pbl=0, pattern_spp_pbl=zero3, &
         bl_mynn_tkeadvect=.false., tke_budget=0, bl_mynn_cloudpdf=2, bl_mynn_mixlength=1, &
         bl_mynn_closure=2.6d0, bl_mynn_edmf=1, bl_mynn_edmf_mom=1, bl_mynn_edmf_tke=0, bl_mynn_output=1, &
         bl_mynn_mixscalars=0, bl_mynn_mixaerosols=1, bl_mynn_mixnumcon=0, bl_mynn_cloudmix=1, &
         bl_mynn_mixqt=0, bl_mynn_edmf_dd=0, bl_mynn_ess=1)
    o2(:,:,7) = real(kpbl, 8)
    if (s < nsteps) then  ! advance the state with the tendencies
      in3(:,:,:,2) = in3(:,:,:,2) + o3(:,:,:,11) * delt
      in3(:,:,:,3) = in3(:,:,:,3) + o3(:,:,:,12) * delt
      in3(:,:,:,5) = in3(:,:,:,5) + o3(:,:,:,13) * delt
      in3(:,:,:,10) = max(in3(:,:,:,10) + o3(:,:,:,14) * delt, 0.d0)
      in3(:,:,:,11) = max(in3(:,:,:,11) + o3(:,:,:,15) * delt, 0.d0)
      in3(:,:,:,12) = max(in3(:,:,:,12) + o3(:,:,:,16) * delt, 0.d0)
      in3(:,:,:,6) = in3(:,:,:,5) * in3(:,:,:,8)
    end if
  end do
  open(11, file=fout, access='stream', form='unformatted', status='replace')
  do n = 1, no3
    write(11) o3(:, :, 1, n)
  end do
  do n = 1, no2
    write(11) o2(:, 1, n)
  end do
  close(11)
end program driver_bl_mynn
