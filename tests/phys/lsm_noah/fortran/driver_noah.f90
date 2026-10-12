program driver_noah
  use module_sf_noahlsm
  implicit none

  character(len=256) :: infile, outfile, table_dir
  integer :: n_cols, nsoil, isurban_in
  real*8 :: dt_in, zlvl_in
  real*8, dimension(8) :: sldpth_in

  integer :: nargs, c, k, lc, iindex, num_slope
  character(len=4) :: sltype_in, lutype_in
  character(len=256) :: llanduse_str = 'USGS', lsoil_str = 'STAS'

  ! Inputs per column
  integer, allocatable, dimension(:) :: vegtyp, soiltyp, slopetyp
  real*8, allocatable, dimension(:) :: shdfac, shdmin, shdmax, alb, snoalb, tbot, z0brd, embrd
  real*8, allocatable, dimension(:) :: lwdn, soldn, solnet, sfcprs, prcp, sfctmp, th2, q2, sfcspd, q2sat, dqsdt2, ffrozp

  ! State variables (INOUT)
  real*8, allocatable, dimension(:) :: cmc, t1, snowh, sneqv, albedo, ch, cm, emissi, z0, xlai, snotime1, ribb
  real*8, allocatable, dimension(:,:) :: stc, smc, sh2o

  ! Outputs
  real*8, allocatable, dimension(:) :: eta, sheat, eta_kinematic, fdown, ec, edir, ett, esnow, drip, dew, beta, etp, ssoil
  real*8, allocatable, dimension(:) :: flx1, flx2, flx3, snomlt, sncovr, runoff1, runoff2, runoff3, rc, pc, rsmin, soilw, soilm, q1
  real*8, allocatable, dimension(:,:) :: et, smav

  ! Temporary column arrays and scalar outputs for SFLX
  real*8, dimension(8) :: stc_col, smc_col, sh2o_col, et_col, smav_col
  real*8 :: eta_val, sheat_val, eta_kin_val, fdown_val, ec_val, edir_val, ett_val, esnow_val
  real*8 :: drip_val, dew_val, beta_val, etp_val, ssoil_val, flx1_val, flx2_val, flx3_val
  real*8 :: snomlt_val, sncovr_val, r1_val, r2_val, r3_val, rc_val, pc_val, rsmin_val, soilw_val, soilm_val, q1_val
  real*8 :: cur_cmc, cur_t1, cur_snowh, cur_sneqv, cur_albedo, cur_ch, cur_cm, cur_emissi, cur_z0, cur_xlai, cur_snotime1, cur_ribb
  real*8 :: cur_shdfac, cur_alb, cur_embrd
  real*8 :: dummy_val, flx4_val, fvb_val, fbur_val, fgsn_val, rcs_val, rct_val, rcq_val, rcsoil_val
  real*8 :: smcwlt_val, smcdry_val, smcref_val, smcmax_val
  real*8 :: sfhead1rt_val, infxs1rt_val, etpnd1_val, aoasis_val
  real*8 :: xsda_qfx_val, hfx_phy_val, qfx_phy_val, xqnorm_val, hcpct_fasdas_val
  integer :: nroot_val, opt_thcnd_val, fasdas_val

  nargs = command_argument_count()
  if (nargs < 1) then
    infile = 'noah_input.bin'
  else
    call get_command_argument(1, infile)
  end if

  if (nargs < 2) then
    outfile = 'noah_output_fortran.bin'
  else
    call get_command_argument(2, outfile)
  end if

  if (nargs < 3) then
    table_dir = '.'
  else
    call get_command_argument(3, table_dir)
  end if

  ! Read GENPARM.TBL
  open(19, file=trim(table_dir)//'/GENPARM.TBL', status='old', action='read')
  read(19,*)
  read(19,*)
  read(19,*) num_slope
  slpcats = num_slope
  do lc = 1, slpcats
     read(19,*) slope_data(lc)
  end do
  read(19,*)
  read(19,*) sbeta_data
  read(19,*)
  read(19,*) fxexp_data
  read(19,*)
  read(19,*) csoil_data
  read(19,*)
  read(19,*) salp_data
  read(19,*)
  read(19,*) refdk_data
  read(19,*)
  read(19,*) refkdt_data
  read(19,*)
  read(19,*) frzk_data
  read(19,*)
  read(19,*) zbot_data
  read(19,*)
  read(19,*) czil_data
  read(19,*)
  read(19,*) smlow_data
  read(19,*)
  read(19,*) smhigh_data
  read(19,*)
  read(19,*) lvcoef_data
  close(19)

  ! Read SOILPARM.TBL
  open(19, file=trim(table_dir)//'/SOILPARM.TBL', status='old', action='read')
  read(19,*)
  read(19,'(A4)') sltype_in
  read(19,*) slcats, iindex
  do lc = 1, slcats
     read(19,*) iindex, bb(lc), drysmc(lc), f11(lc), maxsmc(lc), &
                refsmc(lc), satpsi(lc), satdk(lc), satdw(lc),    &
                wltsmc(lc), qtz(lc)
  end do
  close(19)

  ! Read VEGPARM.TBL
  open(19, file=trim(table_dir)//'/VEGPARM.TBL', status='old', action='read')
  read(19,*)
  read(19,*) lutype_in
  read(19,*) lucats, iindex
  do lc = 1, lucats
     read(19,*) iindex, shdtbl(lc), nrotbl(lc), rstbl(lc), rgltbl(lc), hstbl(lc), &
                snuptbl(lc), maxalb(lc), laimintbl(lc), laimaxtbl(lc), emissmintbl(lc), &
                emissmaxtbl(lc), albedomintbl(lc), albedomaxtbl(lc), z0mintbl(lc), z0maxtbl(lc), &
                ztopvtbl(lc), zbotvtbl(lc)
  end do
  read(19,*)
  read(19,*) topt_data
  read(19,*)
  read(19,*) cmcmax_data
  read(19,*)
  read(19,*) cfactr_data
  read(19,*)
  read(19,*) rsmax_data
  read(19,*)
  read(19,*) bare
  read(19,*)
  read(19,*) natural
  close(19)

  ! Open input file
  open(unit=10, file=trim(infile), status='old', access='stream', form='unformatted')
  read(10) n_cols, nsoil, dt_in, zlvl_in, isurban_in
  read(10) sldpth_in(1:nsoil)

  allocate(vegtyp(n_cols), soiltyp(n_cols), slopetyp(n_cols))
  allocate(shdfac(n_cols), shdmin(n_cols), shdmax(n_cols), alb(n_cols), snoalb(n_cols))
  allocate(tbot(n_cols), z0brd(n_cols), embrd(n_cols))
  allocate(lwdn(n_cols), soldn(n_cols), solnet(n_cols), sfcprs(n_cols), prcp(n_cols))
  allocate(sfctmp(n_cols), th2(n_cols), q2(n_cols), sfcspd(n_cols), q2sat(n_cols), dqsdt2(n_cols), ffrozp(n_cols))

  allocate(cmc(n_cols), t1(n_cols), snowh(n_cols), sneqv(n_cols), albedo(n_cols), ch(n_cols), cm(n_cols))
  allocate(emissi(n_cols), z0(n_cols), xlai(n_cols), snotime1(n_cols), ribb(n_cols))
  allocate(stc(n_cols, nsoil), smc(n_cols, nsoil), sh2o(n_cols, nsoil))

  allocate(eta(n_cols), sheat(n_cols), eta_kinematic(n_cols), fdown(n_cols), ec(n_cols), edir(n_cols))
  allocate(ett(n_cols), esnow(n_cols), drip(n_cols), dew(n_cols), beta(n_cols), etp(n_cols), ssoil(n_cols))
  allocate(flx1(n_cols), flx2(n_cols), flx3(n_cols), snomlt(n_cols), sncovr(n_cols))
  allocate(runoff1(n_cols), runoff2(n_cols), runoff3(n_cols), rc(n_cols), pc(n_cols), rsmin(n_cols))
  allocate(soilw(n_cols), soilm(n_cols), q1(n_cols))
  allocate(et(n_cols, nsoil), smav(n_cols, nsoil))

  read(10) vegtyp
  read(10) soiltyp
  read(10) slopetyp
  read(10) shdfac
  read(10) shdmin
  read(10) shdmax
  read(10) alb
  read(10) snoalb
  read(10) tbot
  read(10) z0brd
  read(10) embrd
  read(10) lwdn
  read(10) soldn
  read(10) solnet
  read(10) sfcprs
  read(10) prcp
  read(10) sfctmp
  read(10) th2
  read(10) q2
  read(10) sfcspd
  read(10) q2sat
  read(10) dqsdt2
  read(10) ffrozp

  read(10) cmc
  read(10) t1
  read(10) snowh
  read(10) sneqv
  read(10) albedo
  read(10) ch
  read(10) cm
  read(10) emissi
  read(10) z0
  read(10) xlai
  read(10) snotime1
  read(10) ribb
  read(10) stc
  read(10) smc
  read(10) sh2o
  close(10)

  ! Run SFLX on each column
  opt_thcnd_val = 1
  aoasis_val = 1.0
  fasdas_val = 0
  xsda_qfx_val = 0.0
  hfx_phy_val = 0.0
  qfx_phy_val = 0.0
  xqnorm_val = 0.0
  sfhead1rt_val = 0.0
  infxs1rt_val = 0.0
  etpnd1_val = 0.0

  do c = 1, n_cols
     ! Skip water points
     if (vegtyp(c) == 16 .or. soiltyp(c) == 14) then
        eta(c) = 0.0
        sheat(c) = 0.0
        eta_kinematic(c) = 0.0
        fdown(c) = 0.0
        ec(c) = 0.0
        edir(c) = 0.0
        ett(c) = 0.0
        esnow(c) = 0.0
        drip(c) = 0.0
        dew(c) = 0.0
        beta(c) = 0.0
        etp(c) = 0.0
        ssoil(c) = 0.0
        flx1(c) = 0.0
        flx2(c) = 0.0
        flx3(c) = 0.0
        snomlt(c) = 0.0
        sncovr(c) = 0.0
        runoff1(c) = 0.0
        runoff2(c) = 0.0
        runoff3(c) = 0.0
        rc(c) = 0.0
        pc(c) = 0.0
        rsmin(c) = 0.0
        soilw(c) = 0.0
        soilm(c) = 0.0
        q1(c) = 0.0
        do k = 1, nsoil
           et(c, k) = 0.0
           smav(c, k) = 0.0
        end do
        cycle
     end if

     cur_cmc = cmc(c)
     cur_t1 = t1(c)
     cur_snowh = snowh(c)
     cur_sneqv = sneqv(c)
     cur_albedo = albedo(c)
     cur_ch = ch(c)
     cur_cm = cm(c)
     cur_emissi = emissi(c)
     cur_z0 = z0(c)
     cur_xlai = xlai(c)
     cur_snotime1 = snotime1(c)
     cur_ribb = ribb(c)
     cur_shdfac = shdfac(c)
     cur_alb = alb(c)
     cur_embrd = embrd(c)

     do k = 1, nsoil
        stc_col(k) = stc(c, k)
        smc_col(k) = smc(c, k)
        sh2o_col(k) = sh2o(c, k)
     end do

     dummy_val = 0.0
     pc_val = 0.0
     rc_val = 0.0

     call SFLX (c, 1, ffrozp(c), isurban_in, dt_in, zlvl_in, nsoil, sldpth_in(1:nsoil), &
                .FALSE., &
                llanduse_str, lsoil_str, &
                lwdn(c), soldn(c), solnet(c), sfcprs(c), prcp(c), sfctmp(c), q2(c), sfcspd(c), &
                dummy_val, dummy_val, dummy_val, &
                th2(c), q2sat(c), dqsdt2(c), &
                vegtyp(c), soiltyp(c), slopetyp(c), cur_shdfac, shdmin(c), shdmax(c), &
                cur_alb, snoalb(c), tbot(c), z0brd(c), cur_z0, cur_emissi, cur_embrd, &
                cur_cmc, cur_t1, stc_col(1:nsoil), smc_col(1:nsoil), sh2o_col(1:nsoil), &
                cur_snowh, cur_sneqv, cur_albedo, cur_ch, cur_cm, &
                eta_val, sheat_val, eta_kin_val, fdown_val, &
                ec_val, edir_val, et_col(1:nsoil), ett_val, esnow_val, drip_val, dew_val, &
                beta_val, etp_val, ssoil_val, &
                flx1_val, flx2_val, flx3_val, &
                flx4_val, fvb_val, fbur_val, fgsn_val, .FALSE., &
                snomlt_val, sncovr_val, &
                r1_val, r2_val, r3_val, &
                rc_val, pc_val, rsmin_val, cur_xlai, rcs_val, rct_val, rcq_val, rcsoil_val, &
                soilw_val, soilm_val, q1_val, smav_col(1:nsoil), &
                .TRUE., .TRUE., &
                cur_snotime1, &
                cur_ribb, &
                smcwlt_val, smcdry_val, smcref_val, smcmax_val, nroot_val, &
                sfhead1rt_val, &
                infxs1rt_val, etpnd1_val, opt_thcnd_val, aoasis_val, &
                xsda_qfx_val, hfx_phy_val, qfx_phy_val, xqnorm_val, &
                fasdas_val, hcpct_fasdas_val)

     ! Store back outputs
     eta(c) = eta_val
     sheat(c) = sheat_val
     eta_kinematic(c) = eta_kin_val
     fdown(c) = fdown_val
     ec(c) = ec_val
     edir(c) = edir_val
     ett(c) = ett_val
     esnow(c) = esnow_val
     drip(c) = drip_val
     dew(c) = dew_val
     beta(c) = beta_val
     etp(c) = etp_val
     ssoil(c) = ssoil_val
     flx1(c) = flx1_val
     flx2(c) = flx2_val
     flx3(c) = flx3_val
     snomlt(c) = snomlt_val
     sncovr(c) = sncovr_val
     runoff1(c) = r1_val
     runoff2(c) = r2_val
     runoff3(c) = r3_val
     rc(c) = rc_val
     pc(c) = pc_val
     rsmin(c) = rsmin_val
     soilw(c) = soilw_val
     soilm(c) = soilm_val
     q1(c) = q1_val

     cmc(c) = cur_cmc
     t1(c) = cur_t1
     snowh(c) = cur_snowh
     sneqv(c) = cur_sneqv
     albedo(c) = cur_albedo
     emissi(c) = cur_emissi
     z0(c) = cur_z0
     xlai(c) = cur_xlai
     snotime1(c) = cur_snotime1
     ribb(c) = cur_ribb

     do k = 1, nsoil
        stc(c, k) = stc_col(k)
        smc(c, k) = smc_col(k)
        sh2o(c, k) = sh2o_col(k)
        et(c, k) = et_col(k)
        smav(c, k) = smav_col(k)
     end do
  end do

  ! Write output binary file
  open(unit=20, file=trim(outfile), status='replace', access='stream', form='unformatted')
  write(20) n_cols, nsoil

  ! 2D outputs
  write(20) eta
  write(20) sheat
  write(20) eta_kinematic
  write(20) fdown
  write(20) ec
  write(20) edir
  write(20) ett
  write(20) esnow
  write(20) drip
  write(20) dew
  write(20) beta
  write(20) etp
  write(20) ssoil
  write(20) flx1
  write(20) flx2
  write(20) flx3
  write(20) snomlt
  write(20) sncovr
  write(20) runoff1
  write(20) runoff2
  write(20) runoff3
  write(20) rc
  write(20) pc
  write(20) rsmin
  write(20) soilw
  write(20) soilm
  write(20) q1

  ! Inout updated states
  write(20) cmc
  write(20) t1
  write(20) snowh
  write(20) sneqv
  write(20) albedo
  write(20) emissi
  write(20) z0
  write(20) xlai
  write(20) snotime1
  write(20) ribb

  ! 3D outputs (nsoil x n_cols)
  write(20) stc
  write(20) smc
  write(20) sh2o
  write(20) et
  write(20) smav
  close(20)

  print *, 'Fortran Noah run complete: ', n_cols, 'columns processed'
end program driver_noah
