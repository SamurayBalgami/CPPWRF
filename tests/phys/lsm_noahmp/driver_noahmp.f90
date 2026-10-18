program driver_noahmp_fortran
  use NoahmpIOVarType, only: NoahmpIO_type
  use NoahmpWRFinitMod, only: NoahmpWRFinit
  use NoahmpWRFmainMod, only: NoahmpWRFmain
  implicit none

  character(len=256) :: infile, outfile, mminlu
  type(NoahmpIO_type) :: NoahmpIO

  integer :: n_cols, nsoil, n_steps, yr, step, i, k, iunit, ounit
  real(8) :: dt, dx, julian_start, julian, xice_threshold, wtddt, soiltstep
  real(8), allocatable :: dzs(:)

  ! Options (WRF v4.8.0 default options)
  integer :: idveg, iopt_crs, iopt_btr, iopt_runsrf, iopt_runsub
  integer :: iopt_sfc, iopt_frz, iopt_inf, iopt_rad, iopt_alb
  integer :: iopt_snf, iopt_tbot, iopt_stc, iopt_gla, iopt_rsf
  integer :: iopt_soil, iopt_pedo, iopt_crop, iopt_irr, iopt_irrm
  integer :: iopt_infdv, iopt_tdrn, iopt_tksno, iopt_compact, iopt_scf, iopt_wetland
  integer :: iz0tlnd, sf_urban_physics, stepwtd
  integer :: snicar_bandnumber_opt, snicar_solarspec_opt, snicar_snowoptics_opt, snicar_dustoptics_opt
  integer :: snicar_rtsolver_opt, snicar_snowshape_opt
  logical :: snicar_use_aerosol, snicar_snowbc_intmix, snicar_snowdust_intmix
  logical :: snicar_use_oc, snicar_aerosol_readtable
  logical :: restart, allowed_to_read, fndsoilw, fndsnowh

  ! 2D Integer arrays (n_cols, 1)
  integer, allocatable :: ivgtyp(:,:), isltyp(:,:), isnowxy(:,:), cropcat(:,:), pgsxy(:,:)
  integer, allocatable :: irnumsi(:,:), irnummi(:,:), irnumfi(:,:)

  ! 2D Real input arrays (n_cols, 1)
  real(8), allocatable :: xlat(:,:), xlong(:,:), coszen(:,:), vegfra(:,:), shdmax(:,:)
  real(8), allocatable :: tmn(:,:), xland(:,:), xice(:,:)
  real(8), allocatable :: swdown(:,:), swddir(:,:), swddif(:,:), glw(:,:)
  real(8), allocatable :: rainbl(:,:), sr(:,:)
  real(8), allocatable :: mp_rainc(:,:), mp_rainnc(:,:), mp_shcv(:,:)
  real(8), allocatable :: mp_snow(:,:), mp_graup(:,:), mp_hail(:,:)
  real(8), allocatable :: dz8w_1(:,:), t_phy_1(:,:), qv_curr_1(:,:), u_phy_1(:,:), v_phy_1(:,:)
  real(8), allocatable :: p8w_1(:,:), p8w_2(:,:)
  real(8), allocatable :: planting(:,:), harvest(:,:), season_gdd(:,:), td_fraction(:,:)
  real(8), allocatable :: irfract(:,:), sifract(:,:), mifract(:,:), fifract(:,:)
  real(8), allocatable :: soilcl1(:,:), soilcl2(:,:), soilcl3(:,:), soilcl4(:,:)
  real(8), allocatable :: fdepthxy(:,:), ht(:,:), msftx(:,:), msfty(:,:), rechclim(:,:)
  real(8), allocatable :: areaxy(:,:), riverbedxy(:,:), eqzwt(:,:), rivercondxy(:,:), pexpxy(:,:)
  real(8), allocatable :: qrfsxy(:,:), qspringsxy(:,:), qslatxy(:,:), chstarxy(:,:)
  real(8), allocatable :: fsatxy(:,:), wsurfxy(:,:)

  ! 3D Real input arrays
  real(8), allocatable :: dz8w(:,:,:), t_phy(:,:,:), qv_curr(:,:,:), u_phy(:,:,:), v_phy(:,:,:), p8w(:,:,:)
  real(8), allocatable :: croptype(:,:,:), soilcomp(:,:,:)

  ! 2D Real INOUT state arrays (n_cols, 1)
  real(8), allocatable :: tsk(:,:), hfx(:,:), qfx(:,:), lh(:,:), grdflx(:,:)
  real(8), allocatable :: smstav(:,:), smstot(:,:), sfcrunoff(:,:), udrunoff(:,:)
  real(8), allocatable :: albedo(:,:), snowc(:,:), snow(:,:), snowh(:,:), canwat(:,:)
  real(8), allocatable :: acsnom(:,:), acsnow(:,:), emiss(:,:), qsfc(:,:), z0(:,:), znt(:,:), rs(:,:)
  real(8), allocatable :: tvxy(:,:), tgxy(:,:), canicexy(:,:), canliqxy(:,:), eahxy(:,:), tahxy(:,:)
  real(8), allocatable :: cmxy(:,:), chxy(:,:), fwetxy(:,:), sneqvoxy(:,:), alboldxy(:,:)
  real(8), allocatable :: qsnowxy(:,:), qrainxy(:,:), wslakexy(:,:), zwtxy(:,:), waxy(:,:), wtxy(:,:)
  real(8), allocatable :: lfmassxy(:,:), rtmassxy(:,:), stmassxy(:,:), woodxy(:,:)
  real(8), allocatable :: stblcpxy(:,:), fastcpxy(:,:), lai(:,:), xsaixy(:,:), taussxy(:,:)
  real(8), allocatable :: smcwtdxy(:,:), deeprechxy(:,:), rechxy(:,:), grainxy(:,:), gddxy(:,:), qtdrain(:,:)
  real(8), allocatable :: acc_ssoilxy(:,:), acc_qinsurxy(:,:), acc_qsevaxy(:,:), acc_dwaterxy(:,:)
  real(8), allocatable :: acc_prcpxy(:,:), acc_ecanxy(:,:), acc_etranxy(:,:), acc_edirxy(:,:)
  real(8), allocatable :: irwatsi(:,:), irwatmi(:,:), irwatfi(:,:), ireloss(:,:)
  real(8), allocatable :: irsivol(:,:), irmivol(:,:), irfivol(:,:), irrsplh(:,:)

  ! 2D Real OUT diagnostic arrays (n_cols, 1)
  real(8), allocatable :: t2mvxy(:,:), t2mbxy(:,:), q2mvxy(:,:), q2mbxy(:,:)
  real(8), allocatable :: tradxy(:,:), neexy(:,:), gppxy(:,:), nppxy(:,:), fvegxy(:,:)
  real(8), allocatable :: runsfxy(:,:), runsbxy(:,:), ecanxy(:,:), edirxy(:,:), etranxy(:,:)
  real(8), allocatable :: fsaxy(:,:), firaxy(:,:), aparxy(:,:), psnxy(:,:), savxy(:,:), sagxy(:,:)
  real(8), allocatable :: rssunxy(:,:), rsshaxy(:,:), bgapxy(:,:), wgapxy(:,:), tgvxy(:,:), tgbxy(:,:)
  real(8), allocatable :: chvxy(:,:), chbxy(:,:), shgxy(:,:), shcxy(:,:), shbxy(:,:)
  real(8), allocatable :: evgxy(:,:), evbxy(:,:), ghvxy(:,:), ghbxy(:,:)
  real(8), allocatable :: irgxy(:,:), ircxy(:,:), irbxy(:,:), trxy(:,:), evcxy(:,:)
  real(8), allocatable :: chleafxy(:,:), chucxy(:,:), chv2xy(:,:), chb2xy(:,:)
  real(8), allocatable :: pahxy(:,:), pahgxy(:,:), pahvxy(:,:), pahbxy(:,:)
  real(8), allocatable :: qintsxy(:,:), qintrxy(:,:), qdripsxy(:,:), qdriprxy(:,:)
  real(8), allocatable :: qthrosxy(:,:), qthrorxy(:,:), qsnsubxy(:,:), qsnfroxy(:,:)
  real(8), allocatable :: qsubcxy(:,:), qfrocxy(:,:), qevacxy(:,:), qdewcxy(:,:)
  real(8), allocatable :: qfrzcxy(:,:), qmeltcxy(:,:), qsnbotxy(:,:), qmeltxy(:,:), pondingxy(:,:)
  real(8), allocatable :: fpicexy(:,:), rainlsm(:,:), snowlsm(:,:)
  real(8), allocatable :: forctlsm(:,:), forcqlsm(:,:), forcplsm(:,:), forczlsm(:,:), forcwlsm(:,:)
  real(8), allocatable :: eflxbxy(:,:), soilenergy(:,:), snowenergy(:,:), canhsxy(:,:)

  ! 3D Real INOUT state arrays
  real(8), allocatable :: smois(:,:,:), sh2o(:,:,:), tslb(:,:,:), smoiseq(:,:,:), acc_etranixy(:,:,:)
  real(8), allocatable :: tsnoxy(:,:,:), snicexy(:,:,:), snliqxy(:,:,:), zsnsoxy(:,:,:)
  real(8), allocatable :: albsoildirxy(:,:,:), albsoildifxy(:,:,:)
  real(8), allocatable :: snrdsxy(:,:,:), snfrxy(:,:,:), bcphixy(:,:,:), bcphoxy(:,:,:)
  real(8), allocatable :: ocphixy(:,:,:), ocphoxy(:,:,:)
  real(8), allocatable :: dust1xy(:,:,:), dust2xy(:,:,:), dust3xy(:,:,:), dust4xy(:,:,:), dust5xy(:,:,:)
  real(8), allocatable :: massconcbcphixy(:,:,:), massconcbcphoxy(:,:,:)
  real(8), allocatable :: massconcocphixy(:,:,:), massconcocphoxy(:,:,:)
  real(8), allocatable :: massconcdust1xy(:,:,:), massconcdust2xy(:,:,:), massconcdust3xy(:,:,:)
  real(8), allocatable :: massconcdust4xy(:,:,:), massconcdust5xy(:,:,:)

  if (command_argument_count() >= 2) then
    call get_command_argument(1, infile)
    call get_command_argument(2, outfile)
  else
    infile = "noahmp_in.bin"
    outfile = "noahmp_wrf.bin"
  end if

  open(newunit=iunit, file=trim(infile), status='old', access='stream', form='unformatted')
  read(iunit) n_cols, nsoil, n_steps, yr
  read(iunit) dt, dx, julian_start, xice_threshold
  allocate(dzs(nsoil))
  read(iunit) (dzs(k), k=1,nsoil)

  ! Allocate 2D integer arrays
  allocate(ivgtyp(n_cols, 1), isltyp(n_cols, 1), isnowxy(n_cols, 1), cropcat(n_cols, 1), pgsxy(n_cols, 1))
  allocate(irnumsi(n_cols, 1), irnummi(n_cols, 1), irnumfi(n_cols, 1))
  isnowxy = 0; cropcat = 0; pgsxy = 0; irnumsi = 0; irnummi = 0; irnumfi = 0

  ! Allocate 2D real inputs
  allocate(xlat(n_cols, 1), xlong(n_cols, 1), coszen(n_cols, 1), vegfra(n_cols, 1), shdmax(n_cols, 1))
  allocate(tmn(n_cols, 1), xland(n_cols, 1), xice(n_cols, 1))
  allocate(swdown(n_cols, 1), swddir(n_cols, 1), swddif(n_cols, 1), glw(n_cols, 1))
  allocate(rainbl(n_cols, 1), sr(n_cols, 1))
  allocate(mp_rainc(n_cols, 1), mp_rainnc(n_cols, 1), mp_shcv(n_cols, 1))
  allocate(mp_snow(n_cols, 1), mp_graup(n_cols, 1), mp_hail(n_cols, 1))
  allocate(dz8w_1(n_cols, 1), t_phy_1(n_cols, 1), qv_curr_1(n_cols, 1), u_phy_1(n_cols, 1), v_phy_1(n_cols, 1))
  allocate(p8w_1(n_cols, 1), p8w_2(n_cols, 1))
  allocate(planting(n_cols, 1), harvest(n_cols, 1), season_gdd(n_cols, 1), td_fraction(n_cols, 1))
  allocate(irfract(n_cols, 1), sifract(n_cols, 1), mifract(n_cols, 1), fifract(n_cols, 1))
  allocate(soilcl1(n_cols, 1), soilcl2(n_cols, 1), soilcl3(n_cols, 1), soilcl4(n_cols, 1))
  allocate(fdepthxy(n_cols, 1), ht(n_cols, 1), msftx(n_cols, 1), msfty(n_cols, 1), rechclim(n_cols, 1))
  allocate(areaxy(n_cols, 1), riverbedxy(n_cols, 1), eqzwt(n_cols, 1), rivercondxy(n_cols, 1), pexpxy(n_cols, 1))
  allocate(qrfsxy(n_cols, 1), qspringsxy(n_cols, 1), qslatxy(n_cols, 1), chstarxy(n_cols, 1))
  allocate(fsatxy(n_cols, 1), wsurfxy(n_cols, 1))
  swddir = 0.0d0; swddif = 0.0d0
  planting = 0.0d0; harvest = 0.0d0; season_gdd = 0.0d0; td_fraction = 0.0d0
  irfract = 0.0d0; sifract = 0.0d0; mifract = 0.0d0; fifract = 0.0d0
  soilcl1 = 0.0d0; soilcl2 = 0.0d0; soilcl3 = 0.0d0; soilcl4 = 0.0d0
  fdepthxy = 0.0d0; ht = 0.0d0; msftx = 1.0d0; msfty = 1.0d0; rechclim = 0.0d0
  areaxy = 0.0d0; riverbedxy = 0.0d0; eqzwt = 0.0d0; rivercondxy = 0.0d0; pexpxy = 0.0d0
  qrfsxy = 0.0d0; qspringsxy = 0.0d0; qslatxy = 0.0d0; chstarxy = 0.0d0
  fsatxy = 0.0d0; wsurfxy = 0.0d0

  ! Allocate 3D real inputs
  allocate(dz8w(n_cols, 2, 1), t_phy(n_cols, 2, 1), qv_curr(n_cols, 2, 1))
  allocate(u_phy(n_cols, 2, 1), v_phy(n_cols, 2, 1), p8w(n_cols, 2, 1))
  allocate(croptype(n_cols, 5, 1), soilcomp(n_cols, nsoil*2, 1))
  dz8w = 0.0d0; t_phy = 0.0d0; qv_curr = 0.0d0; u_phy = 0.0d0; v_phy = 0.0d0; p8w = 0.0d0
  croptype = 0.0d0; soilcomp = 0.0d0

  ! Allocate 2D INOUT state arrays
  allocate(tsk(n_cols, 1), hfx(n_cols, 1), qfx(n_cols, 1), lh(n_cols, 1), grdflx(n_cols, 1))
  allocate(smstav(n_cols, 1), smstot(n_cols, 1), sfcrunoff(n_cols, 1), udrunoff(n_cols, 1))
  allocate(albedo(n_cols, 1), snowc(n_cols, 1), snow(n_cols, 1), snowh(n_cols, 1), canwat(n_cols, 1))
  allocate(acsnom(n_cols, 1), acsnow(n_cols, 1), emiss(n_cols, 1), qsfc(n_cols, 1), z0(n_cols, 1), znt(n_cols, 1), rs(n_cols, 1))
  allocate(tvxy(n_cols, 1), tgxy(n_cols, 1), canicexy(n_cols, 1), canliqxy(n_cols, 1), eahxy(n_cols, 1), tahxy(n_cols, 1))
  allocate(cmxy(n_cols, 1), chxy(n_cols, 1), fwetxy(n_cols, 1), sneqvoxy(n_cols, 1), alboldxy(n_cols, 1))
  allocate(qsnowxy(n_cols, 1), qrainxy(n_cols, 1), wslakexy(n_cols, 1), zwtxy(n_cols, 1), waxy(n_cols, 1), wtxy(n_cols, 1))
  allocate(lfmassxy(n_cols, 1), rtmassxy(n_cols, 1), stmassxy(n_cols, 1), woodxy(n_cols, 1))
  allocate(stblcpxy(n_cols, 1), fastcpxy(n_cols, 1), lai(n_cols, 1), xsaixy(n_cols, 1), taussxy(n_cols, 1))
  allocate(smcwtdxy(n_cols, 1), deeprechxy(n_cols, 1), rechxy(n_cols, 1), grainxy(n_cols, 1), gddxy(n_cols, 1), qtdrain(n_cols, 1))
  allocate(acc_ssoilxy(n_cols, 1), acc_qinsurxy(n_cols, 1), acc_qsevaxy(n_cols, 1), acc_dwaterxy(n_cols, 1))
  allocate(acc_prcpxy(n_cols, 1), acc_ecanxy(n_cols, 1), acc_etranxy(n_cols, 1), acc_edirxy(n_cols, 1))
  allocate(irwatsi(n_cols, 1), irwatmi(n_cols, 1), irwatfi(n_cols, 1), ireloss(n_cols, 1))
  allocate(irsivol(n_cols, 1), irmivol(n_cols, 1), irfivol(n_cols, 1), irrsplh(n_cols, 1))
  hfx = 0.0d0; qfx = 0.0d0; lh = 0.0d0; grdflx = 0.0d0
  smstav = 0.0d0; smstot = 0.0d0; sfcrunoff = 0.0d0; udrunoff = 0.0d0
  albedo = 0.18d0; snowc = 0.0d0; acsnom = 0.0d0; acsnow = 0.0d0
  emiss = 0.98d0; qsfc = 0.005d0; z0 = 0.05d0; znt = 0.05d0; rs = 0.0d0
  tvxy = 0.0d0; tgxy = 0.0d0; canicexy = 0.0d0; canliqxy = 0.0d0; eahxy = 0.0d0; tahxy = 0.0d0
  cmxy = 0.0d0; chxy = 0.0d0; fwetxy = 0.0d0; sneqvoxy = 0.0d0; alboldxy = 0.0d0
  qsnowxy = 0.0d0; qrainxy = 0.0d0; wslakexy = 0.0d0; zwtxy = 0.0d0; waxy = 0.0d0; wtxy = 0.0d0
  lfmassxy = 0.0d0; rtmassxy = 0.0d0; stmassxy = 0.0d0; woodxy = 0.0d0
  stblcpxy = 0.0d0; fastcpxy = 0.0d0; xsaixy = 0.0d0; taussxy = 0.0d0
  smcwtdxy = 0.0d0; deeprechxy = 0.0d0; rechxy = 0.0d0; grainxy = 0.0d0; gddxy = 0.0d0; qtdrain = 0.0d0
  acc_ssoilxy = 0.0d0; acc_qinsurxy = 0.0d0; acc_qsevaxy = 0.0d0; acc_dwaterxy = 0.0d0
  acc_prcpxy = 0.0d0; acc_ecanxy = 0.0d0; acc_etranxy = 0.0d0; acc_edirxy = 0.0d0
  irwatsi = 0.0d0; irwatmi = 0.0d0; irwatfi = 0.0d0; ireloss = 0.0d0
  irsivol = 0.0d0; irmivol = 0.0d0; irfivol = 0.0d0; irrsplh = 0.0d0

  ! Allocate 2D OUT diagnostic arrays
  allocate(t2mvxy(n_cols, 1), t2mbxy(n_cols, 1), q2mvxy(n_cols, 1), q2mbxy(n_cols, 1))
  allocate(tradxy(n_cols, 1), neexy(n_cols, 1), gppxy(n_cols, 1), nppxy(n_cols, 1), fvegxy(n_cols, 1))
  allocate(runsfxy(n_cols, 1), runsbxy(n_cols, 1), ecanxy(n_cols, 1), edirxy(n_cols, 1), etranxy(n_cols, 1))
  allocate(fsaxy(n_cols, 1), firaxy(n_cols, 1), aparxy(n_cols, 1), psnxy(n_cols, 1), savxy(n_cols, 1), sagxy(n_cols, 1))
  allocate(rssunxy(n_cols, 1), rsshaxy(n_cols, 1), bgapxy(n_cols, 1), wgapxy(n_cols, 1), tgvxy(n_cols, 1), tgbxy(n_cols, 1))
  allocate(chvxy(n_cols, 1), chbxy(n_cols, 1), shgxy(n_cols, 1), shcxy(n_cols, 1), shbxy(n_cols, 1))
  allocate(evgxy(n_cols, 1), evbxy(n_cols, 1), ghvxy(n_cols, 1), ghbxy(n_cols, 1))
  allocate(irgxy(n_cols, 1), ircxy(n_cols, 1), irbxy(n_cols, 1), trxy(n_cols, 1), evcxy(n_cols, 1))
  allocate(chleafxy(n_cols, 1), chucxy(n_cols, 1), chv2xy(n_cols, 1), chb2xy(n_cols, 1))
  allocate(pahxy(n_cols, 1), pahgxy(n_cols, 1), pahvxy(n_cols, 1), pahbxy(n_cols, 1))
  allocate(qintsxy(n_cols, 1), qintrxy(n_cols, 1), qdripsxy(n_cols, 1), qdriprxy(n_cols, 1))
  allocate(qthrosxy(n_cols, 1), qthrorxy(n_cols, 1), qsnsubxy(n_cols, 1), qsnfroxy(n_cols, 1))
  allocate(qsubcxy(n_cols, 1), qfrocxy(n_cols, 1), qevacxy(n_cols, 1), qdewcxy(n_cols, 1))
  allocate(qfrzcxy(n_cols, 1), qmeltcxy(n_cols, 1), qsnbotxy(n_cols, 1), qmeltxy(n_cols, 1), pondingxy(n_cols, 1))
  allocate(fpicexy(n_cols, 1), rainlsm(n_cols, 1), snowlsm(n_cols, 1))
  allocate(forctlsm(n_cols, 1), forcqlsm(n_cols, 1), forcplsm(n_cols, 1), forczlsm(n_cols, 1), forcwlsm(n_cols, 1))
  allocate(eflxbxy(n_cols, 1), soilenergy(n_cols, 1), snowenergy(n_cols, 1), canhsxy(n_cols, 1))
  t2mvxy = 0.0d0; t2mbxy = 0.0d0; q2mvxy = 0.0d0; q2mbxy = 0.0d0
  tradxy = 0.0d0; neexy = 0.0d0; gppxy = 0.0d0; nppxy = 0.0d0; fvegxy = 0.0d0
  runsfxy = 0.0d0; runsbxy = 0.0d0; ecanxy = 0.0d0; edirxy = 0.0d0; etranxy = 0.0d0
  fsaxy = 0.0d0; firaxy = 0.0d0; aparxy = 0.0d0; psnxy = 0.0d0; savxy = 0.0d0; sagxy = 0.0d0
  rssunxy = 0.0d0; rsshaxy = 0.0d0; bgapxy = 0.0d0; wgapxy = 0.0d0; tgvxy = 0.0d0; tgbxy = 0.0d0
  chvxy = 0.0d0; chbxy = 0.0d0; shgxy = 0.0d0; shcxy = 0.0d0; shbxy = 0.0d0
  evgxy = 0.0d0; evbxy = 0.0d0; ghvxy = 0.0d0; ghbxy = 0.0d0
  irgxy = 0.0d0; ircxy = 0.0d0; irbxy = 0.0d0; trxy = 0.0d0; evcxy = 0.0d0
  chleafxy = 0.0d0; chucxy = 0.0d0; chv2xy = 0.0d0; chb2xy = 0.0d0
  pahxy = 0.0d0; pahgxy = 0.0d0; pahvxy = 0.0d0; pahbxy = 0.0d0
  qintsxy = 0.0d0; qintrxy = 0.0d0; qdripsxy = 0.0d0; qdriprxy = 0.0d0
  qthrosxy = 0.0d0; qthrorxy = 0.0d0; qsnsubxy = 0.0d0; qsnfroxy = 0.0d0
  qsubcxy = 0.0d0; qfrocxy = 0.0d0; qevacxy = 0.0d0; qdewcxy = 0.0d0
  qfrzcxy = 0.0d0; qmeltcxy = 0.0d0; qsnbotxy = 0.0d0; qmeltxy = 0.0d0; pondingxy = 0.0d0
  fpicexy = 0.0d0; rainlsm = 0.0d0; snowlsm = 0.0d0
  forctlsm = 0.0d0; forcqlsm = 0.0d0; forcplsm = 0.0d0; forczlsm = 0.0d0; forcwlsm = 0.0d0
  eflxbxy = 0.0d0; soilenergy = 0.0d0; snowenergy = 0.0d0; canhsxy = 0.0d0

  ! Allocate 3D INOUT arrays
  allocate(smois(n_cols, nsoil, 1), sh2o(n_cols, nsoil, 1), tslb(n_cols, nsoil, 1))
  allocate(smoiseq(n_cols, nsoil, 1), acc_etranixy(n_cols, nsoil, 1))
  allocate(tsnoxy(n_cols, -2:0, 1), snicexy(n_cols, -2:0, 1), snliqxy(n_cols, -2:0, 1))
  allocate(zsnsoxy(n_cols, -2:nsoil, 1))
  allocate(albsoildirxy(n_cols, 1:2, 1), albsoildifxy(n_cols, 1:2, 1))
  allocate(snrdsxy(n_cols, -2:0, 1), snfrxy(n_cols, -2:0, 1), bcphixy(n_cols, -2:0, 1), bcphoxy(n_cols, -2:0, 1))
  allocate(ocphixy(n_cols, -2:0, 1), ocphoxy(n_cols, -2:0, 1))
  allocate(dust1xy(n_cols, -2:0, 1), dust2xy(n_cols, -2:0, 1), dust3xy(n_cols, -2:0, 1))
  allocate(dust4xy(n_cols, -2:0, 1), dust5xy(n_cols, -2:0, 1))
  allocate(massconcbcphixy(n_cols, -2:0, 1), massconcbcphoxy(n_cols, -2:0, 1))
  allocate(massconcocphixy(n_cols, -2:0, 1), massconcocphoxy(n_cols, -2:0, 1))
  allocate(massconcdust1xy(n_cols, -2:0, 1), massconcdust2xy(n_cols, -2:0, 1), massconcdust3xy(n_cols, -2:0, 1))
  allocate(massconcdust4xy(n_cols, -2:0, 1), massconcdust5xy(n_cols, -2:0, 1))
  smoiseq = 0.0d0; acc_etranixy = 0.0d0
  tsnoxy = 0.0d0; snicexy = 0.0d0; snliqxy = 0.0d0; zsnsoxy = 0.0d0
  albsoildirxy = 0.0d0; albsoildifxy = 0.0d0
  snrdsxy = 0.0d0; snfrxy = 0.0d0; bcphixy = 0.0d0; bcphoxy = 0.0d0; ocphixy = 0.0d0; ocphoxy = 0.0d0
  dust1xy = 0.0d0; dust2xy = 0.0d0; dust3xy = 0.0d0; dust4xy = 0.0d0; dust5xy = 0.0d0
  massconcbcphixy = 0.0d0; massconcbcphoxy = 0.0d0; massconcocphixy = 0.0d0; massconcocphoxy = 0.0d0
  massconcdust1xy = 0.0d0; massconcdust2xy = 0.0d0; massconcdust3xy = 0.0d0
  massconcdust4xy = 0.0d0; massconcdust5xy = 0.0d0

  ! Read 2D integer inputs
  read(iunit) (ivgtyp(i, 1), i=1,n_cols)
  read(iunit) (isltyp(i, 1), i=1,n_cols)

  ! Read 2D real inputs
  read(iunit) (xlat(i, 1), i=1,n_cols)
  read(iunit) (xlong(i, 1), i=1,n_cols)
  read(iunit) (coszen(i, 1), i=1,n_cols)
  read(iunit) (vegfra(i, 1), i=1,n_cols)
  read(iunit) (shdmax(i, 1), i=1,n_cols)
  read(iunit) (tmn(i, 1), i=1,n_cols)
  read(iunit) (xland(i, 1), i=1,n_cols)
  read(iunit) (xice(i, 1), i=1,n_cols)
  read(iunit) (swdown(i, 1), i=1,n_cols)
  read(iunit) (glw(i, 1), i=1,n_cols)
  read(iunit) (rainbl(i, 1), i=1,n_cols)
  read(iunit) (sr(i, 1), i=1,n_cols)
  read(iunit) (mp_rainc(i, 1), i=1,n_cols)
  read(iunit) (mp_rainnc(i, 1), i=1,n_cols)
  read(iunit) (mp_shcv(i, 1), i=1,n_cols)
  read(iunit) (mp_snow(i, 1), i=1,n_cols)
  read(iunit) (mp_graup(i, 1), i=1,n_cols)
  read(iunit) (mp_hail(i, 1), i=1,n_cols)
  read(iunit) (dz8w_1(i, 1), i=1,n_cols)
  read(iunit) (t_phy_1(i, 1), i=1,n_cols)
  read(iunit) (qv_curr_1(i, 1), i=1,n_cols)
  read(iunit) (u_phy_1(i, 1), i=1,n_cols)
  read(iunit) (v_phy_1(i, 1), i=1,n_cols)
  read(iunit) (p8w_1(i, 1), i=1,n_cols)
  read(iunit) (p8w_2(i, 1), i=1,n_cols)
  read(iunit) (tsk(i, 1), i=1,n_cols)
  read(iunit) (snow(i, 1), i=1,n_cols)
  read(iunit) (snowh(i, 1), i=1,n_cols)
  read(iunit) (canwat(i, 1), i=1,n_cols)
  read(iunit) (lai(i, 1), i=1,n_cols)

  ! Read 3D initial soil state (k=1..nsoil, i=1..n_cols)
  read(iunit) ((tslb(i, k, 1), i=1,n_cols), k=1,nsoil)
  read(iunit) ((smois(i, k, 1), i=1,n_cols), k=1,nsoil)
  read(iunit) ((sh2o(i, k, 1), i=1,n_cols), k=1,nsoil)
  close(iunit)

  do i = 1, n_cols
    dz8w(i, 1, 1) = dz8w_1(i, 1)
    t_phy(i, 1, 1) = t_phy_1(i, 1)
    qv_curr(i, 1, 1) = qv_curr_1(i, 1)
    u_phy(i, 1, 1) = u_phy_1(i, 1)
    v_phy(i, 1, 1) = v_phy_1(i, 1)
    p8w(i, 1, 1) = p8w_1(i, 1)
    p8w(i, 2, 1) = p8w_2(i, 1)
  end do

  ! Default WRF v4.8.0 Noah-MP options
  mminlu = "MODIFIED_IGBP_MODIS_NOAH"
  idveg = 4
  iopt_crs = 1
  iopt_btr = 1
  iopt_runsrf = 3
  iopt_runsub = 3
  iopt_sfc = 1
  iopt_frz = 1
  iopt_inf = 1
  iopt_rad = 3
  iopt_alb = 2
  iopt_snf = 1
  iopt_tbot = 2
  iopt_stc = 1
  iopt_gla = 1
  iopt_rsf = 1
  iopt_soil = 1
  iopt_pedo = 1
  iopt_crop = 0
  iopt_irr = 0
  iopt_irrm = 0
  iopt_infdv = 0
  iopt_tdrn = 0
  iopt_tksno = 1
  iopt_compact = 2
  iopt_scf = 1
  iopt_wetland = 0
  iz0tlnd = 0
  sf_urban_physics = 0
  soiltstep = 0.0d0
  wtddt = 30.0d0
  restart = .false.
  allowed_to_read = .true.
  fndsoilw = .true.
  fndsnowh = .true.
  snicar_bandnumber_opt = 1
  snicar_solarspec_opt = 1
  snicar_snowoptics_opt = 3
  snicar_dustoptics_opt = 1
  snicar_rtsolver_opt = 2
  snicar_snowshape_opt = 3
  snicar_use_aerosol = .true.
  snicar_snowbc_intmix = .true.
  snicar_snowdust_intmix = .false.
  snicar_use_oc = .false.
  snicar_aerosol_readtable = .true.

  ! 1. Call NoahmpWRFinit (note ide=n_cols+1, jde=2 because NoahmpWRFinit uses min0(ite, ide-1))
  call NoahmpWRFinit(NoahmpIO, trim(mminlu), snow, snowh, canwat, isltyp, ivgtyp, xlat, &
    tslb, smois, sh2o, dzs, fndsoilw, fndsnowh, &
    tsk, isnowxy, tvxy, tgxy, canicexy, tmn, xice, &
    canliqxy, eahxy, tahxy, cmxy, chxy, &
    fwetxy, sneqvoxy, alboldxy, qsnowxy, qrainxy, wslakexy, zwtxy, waxy, &
    wtxy, tsnoxy, zsnsoxy, snicexy, snliqxy, lfmassxy, rtmassxy, &
    stmassxy, woodxy, stblcpxy, fastcpxy, xsaixy, lai, &
    grainxy, gddxy, &
    croptype, cropcat, &
    irnumsi, irnummi, irnumfi, irwatsi, &
    irwatmi, irwatfi, ireloss, irsivol, &
    irmivol, irfivol, irrsplh, &
    t2mvxy, t2mbxy, chstarxy, fsatxy, wsurfxy, &
    snrdsxy, snfrxy, bcphixy, bcphoxy, ocphixy, ocphoxy, dust1xy, &
    dust2xy, dust3xy, dust4xy, dust5xy, massconcbcphixy, massconcbcphoxy, &
    massconcocphixy, massconcocphoxy, massconcdust1xy, massconcdust2xy, &
    massconcdust3xy, massconcdust4xy, massconcdust5xy, &
    albsoildirxy, albsoildifxy, &
    nsoil, restart, &
    allowed_to_read, iopt_runsub, iopt_crop, iopt_irr, iopt_irrm, &
    sf_urban_physics, iopt_soil, iopt_alb, iopt_wetland, &
    snicar_snowoptics_opt, snicar_dustoptics_opt, snicar_solarspec_opt, &
    snicar_bandnumber_opt, &
    1, n_cols+1, 1, 2, 1, 2, &
    1, n_cols,   1, 1, 1, 2, &
    1, n_cols,   1, 1, 1, 2, &
    smoiseq, smcwtdxy, rechxy, deeprechxy, qtdrain, areaxy, dx, dx, msftx, msfty, &
    wtddt, stepwtd, dt, qrfsxy, qspringsxy, qslatxy, &
    fdepthxy, ht, riverbedxy, eqzwt, rivercondxy, pexpxy, rechclim)

  ! 2. Call NoahmpWRFmain for n_steps
  do step = 1, n_steps
    julian = julian_start + dble(step - 1) * dt / 86400.0d0
    call NoahmpWRFmain(NoahmpIO, step, yr, julian, coszen, xlat, xlong, &
      dz8w, dt, dzs, nsoil, dx, &
      ivgtyp, isltyp, vegfra, shdmax, tmn, &
      xland, xice, xice_threshold, cropcat, &
      planting, harvest, season_gdd, &
      idveg, iopt_crs, iopt_btr, iopt_runsub, iopt_sfc, iopt_frz, &
      iopt_inf, iopt_rad, iopt_alb, iopt_snf, iopt_tbot, iopt_stc, &
      iopt_gla, iopt_rsf, iopt_soil, iopt_pedo, iopt_crop, iopt_irr, &
      iopt_irrm, iopt_infdv, iopt_tdrn, soiltstep, &
      iopt_runsrf, iopt_tksno, iopt_compact, iopt_scf, iopt_wetland, &
      iz0tlnd, sf_urban_physics, &
      soilcomp, soilcl1, soilcl2, soilcl3, soilcl4, &
      t_phy, qv_curr, u_phy, v_phy, swdown, swddir, &
      swddif, glw, &
      p8w, rainbl, sr, &
      irfract, sifract, mifract, fifract, &
      tsk, hfx, qfx, lh, grdflx, smstav, &
      smstot, sfcrunoff, udrunoff, albedo, snowc, smois, &
      sh2o, tslb, snow, snowh, canwat, acsnom, &
      acsnow, emiss, qsfc, &
      z0, znt, &
      irnumsi, irnummi, irnumfi, irwatsi, irwatmi, irwatfi, &
      ireloss, irsivol, irmivol, irfivol, irrsplh, mminlu, &
      isnowxy, tvxy, tgxy, canicexy, canliqxy, eahxy, &
      tahxy, cmxy, chxy, fwetxy, sneqvoxy, alboldxy, &
      qsnowxy, qrainxy, wslakexy, zwtxy, waxy, wtxy, tsnoxy, &
      zsnsoxy, snicexy, snliqxy, lfmassxy, rtmassxy, stmassxy, &
      woodxy, stblcpxy, fastcpxy, lai, xsaixy, taussxy, &
      smoiseq, smcwtdxy, deeprechxy, rechxy, grainxy, gddxy, pgsxy, &
      qtdrain, td_fraction, &
      t2mvxy, t2mbxy, q2mvxy, q2mbxy, &
      tradxy, neexy, gppxy, nppxy, fvegxy, runsfxy, &
      runsbxy, ecanxy, edirxy, etranxy, fsaxy, firaxy, &
      aparxy, psnxy, savxy, sagxy, rssunxy, rsshaxy, &
      bgapxy, wgapxy, tgvxy, tgbxy, chvxy, chbxy, &
      shgxy, shcxy, shbxy, evgxy, evbxy, ghvxy, &
      ghbxy, irgxy, ircxy, irbxy, trxy, evcxy, &
      chleafxy, chucxy, chv2xy, chb2xy, rs, &
      qintsxy, qintrxy, qdripsxy, &
      qdriprxy, qthrosxy, qthrorxy, &
      qsnsubxy, qsnfroxy, qsubcxy, &
      qfrocxy, qevacxy, qdewcxy, qfrzcxy, qmeltcxy, &
      qsnbotxy, qmeltxy, pondingxy, pahxy, pahgxy, pahvxy, pahbxy, &
      fpicexy, rainlsm, snowlsm, forctlsm, forcqlsm, forcplsm, forczlsm, &
      forcwlsm, acc_ssoilxy, acc_qinsurxy, acc_qsevaxy, acc_etranixy, &
      eflxbxy, soilenergy, snowenergy, canhsxy, &
      acc_dwaterxy, acc_prcpxy, acc_ecanxy, acc_etranxy, acc_edirxy, &
      fsatxy, wsurfxy, &
      snicar_bandnumber_opt, snicar_solarspec_opt, &
      snicar_snowoptics_opt, snicar_dustoptics_opt, &
      snicar_rtsolver_opt, snicar_snowshape_opt, &
      snicar_use_aerosol, snicar_snowbc_intmix, &
      snicar_snowdust_intmix, snicar_use_oc, &
      snicar_aerosol_readtable, snrdsxy, snfrxy, bcphixy, bcphoxy, &
      ocphixy, ocphoxy, dust1xy, dust2xy, dust3xy, dust4xy, dust5xy, &
      massconcbcphixy, massconcbcphoxy, massconcocphixy, &
      massconcocphoxy, massconcdust1xy, massconcdust2xy, &
      massconcdust3xy, massconcdust4xy, massconcdust5xy, &
      albsoildirxy, albsoildifxy, &
      1, n_cols+1, 1, 2, 1, 2, &
      1, n_cols,   1, 1, 1, 2, &
      1, n_cols,   1, 1, 1, 1, &
      mp_rainc, mp_rainnc, mp_shcv, mp_snow, mp_graup, mp_hail)
  end do

  print *, "Fortran Col 1: tsk =", tsk(1,1), "hfx =", hfx(1,1), "lh =", lh(1,1)

  open(newunit=ounit, file=trim(outfile), status='replace', access='stream', form='unformatted')
  write(ounit) n_cols, nsoil

  ! 62 2D INOUT state arrays
  write(ounit) (tsk(i, 1), i=1,n_cols)
  write(ounit) (hfx(i, 1), i=1,n_cols)
  write(ounit) (qfx(i, 1), i=1,n_cols)
  write(ounit) (lh(i, 1), i=1,n_cols)
  write(ounit) (grdflx(i, 1), i=1,n_cols)
  write(ounit) (smstav(i, 1), i=1,n_cols)
  write(ounit) (smstot(i, 1), i=1,n_cols)
  write(ounit) (sfcrunoff(i, 1), i=1,n_cols)
  write(ounit) (udrunoff(i, 1), i=1,n_cols)
  write(ounit) (albedo(i, 1), i=1,n_cols)
  write(ounit) (snowc(i, 1), i=1,n_cols)
  write(ounit) (snow(i, 1), i=1,n_cols)
  write(ounit) (snowh(i, 1), i=1,n_cols)
  write(ounit) (canwat(i, 1), i=1,n_cols)
  write(ounit) (acsnom(i, 1), i=1,n_cols)
  write(ounit) (acsnow(i, 1), i=1,n_cols)
  write(ounit) (emiss(i, 1), i=1,n_cols)
  write(ounit) (qsfc(i, 1), i=1,n_cols)
  write(ounit) (z0(i, 1), i=1,n_cols)
  write(ounit) (znt(i, 1), i=1,n_cols)
  write(ounit) (rs(i, 1), i=1,n_cols)
  write(ounit) (dble(isnowxy(i, 1)), i=1,n_cols)
  write(ounit) (tvxy(i, 1), i=1,n_cols)
  write(ounit) (tgxy(i, 1), i=1,n_cols)
  write(ounit) (canicexy(i, 1), i=1,n_cols)
  write(ounit) (canliqxy(i, 1), i=1,n_cols)
  write(ounit) (eahxy(i, 1), i=1,n_cols)
  write(ounit) (tahxy(i, 1), i=1,n_cols)
  write(ounit) (cmxy(i, 1), i=1,n_cols)
  write(ounit) (chxy(i, 1), i=1,n_cols)
  write(ounit) (fwetxy(i, 1), i=1,n_cols)
  write(ounit) (sneqvoxy(i, 1), i=1,n_cols)
  write(ounit) (alboldxy(i, 1), i=1,n_cols)
  write(ounit) (qsnowxy(i, 1), i=1,n_cols)
  write(ounit) (qrainxy(i, 1), i=1,n_cols)
  write(ounit) (wslakexy(i, 1), i=1,n_cols)
  write(ounit) (zwtxy(i, 1), i=1,n_cols)
  write(ounit) (waxy(i, 1), i=1,n_cols)
  write(ounit) (wtxy(i, 1), i=1,n_cols)
  write(ounit) (lfmassxy(i, 1), i=1,n_cols)
  write(ounit) (rtmassxy(i, 1), i=1,n_cols)
  write(ounit) (stmassxy(i, 1), i=1,n_cols)
  write(ounit) (woodxy(i, 1), i=1,n_cols)
  write(ounit) (stblcpxy(i, 1), i=1,n_cols)
  write(ounit) (fastcpxy(i, 1), i=1,n_cols)
  write(ounit) (lai(i, 1), i=1,n_cols)
  write(ounit) (xsaixy(i, 1), i=1,n_cols)
  write(ounit) (taussxy(i, 1), i=1,n_cols)
  write(ounit) (smcwtdxy(i, 1), i=1,n_cols)
  write(ounit) (deeprechxy(i, 1), i=1,n_cols)
  write(ounit) (rechxy(i, 1), i=1,n_cols)
  write(ounit) (grainxy(i, 1), i=1,n_cols)
  write(ounit) (gddxy(i, 1), i=1,n_cols)
  write(ounit) (qtdrain(i, 1), i=1,n_cols)
  write(ounit) (acc_ssoilxy(i, 1), i=1,n_cols)
  write(ounit) (acc_qinsurxy(i, 1), i=1,n_cols)
  write(ounit) (acc_qsevaxy(i, 1), i=1,n_cols)
  write(ounit) (acc_dwaterxy(i, 1), i=1,n_cols)
  write(ounit) (acc_prcpxy(i, 1), i=1,n_cols)
  write(ounit) (acc_ecanxy(i, 1), i=1,n_cols)
  write(ounit) (acc_etranxy(i, 1), i=1,n_cols)
  write(ounit) (acc_edirxy(i, 1), i=1,n_cols)

  ! 77 2D OUT diagnostic arrays
  write(ounit) (t2mvxy(i, 1), i=1,n_cols)
  write(ounit) (t2mbxy(i, 1), i=1,n_cols)
  write(ounit) (q2mvxy(i, 1), i=1,n_cols)
  write(ounit) (q2mbxy(i, 1), i=1,n_cols)
  write(ounit) (tradxy(i, 1), i=1,n_cols)
  write(ounit) (neexy(i, 1), i=1,n_cols)
  write(ounit) (gppxy(i, 1), i=1,n_cols)
  write(ounit) (nppxy(i, 1), i=1,n_cols)
  write(ounit) (fvegxy(i, 1), i=1,n_cols)
  write(ounit) (runsfxy(i, 1), i=1,n_cols)
  write(ounit) (runsbxy(i, 1), i=1,n_cols)
  write(ounit) (ecanxy(i, 1), i=1,n_cols)
  write(ounit) (edirxy(i, 1), i=1,n_cols)
  write(ounit) (etranxy(i, 1), i=1,n_cols)
  write(ounit) (fsaxy(i, 1), i=1,n_cols)
  write(ounit) (firaxy(i, 1), i=1,n_cols)
  write(ounit) (aparxy(i, 1), i=1,n_cols)
  write(ounit) (psnxy(i, 1), i=1,n_cols)
  write(ounit) (savxy(i, 1), i=1,n_cols)
  write(ounit) (sagxy(i, 1), i=1,n_cols)
  write(ounit) (rssunxy(i, 1), i=1,n_cols)
  write(ounit) (rsshaxy(i, 1), i=1,n_cols)
  write(ounit) (bgapxy(i, 1), i=1,n_cols)
  write(ounit) (wgapxy(i, 1), i=1,n_cols)
  write(ounit) (tgvxy(i, 1), i=1,n_cols)
  write(ounit) (tgbxy(i, 1), i=1,n_cols)
  write(ounit) (chvxy(i, 1), i=1,n_cols)
  write(ounit) (chbxy(i, 1), i=1,n_cols)
  write(ounit) (shgxy(i, 1), i=1,n_cols)
  write(ounit) (shcxy(i, 1), i=1,n_cols)
  write(ounit) (shbxy(i, 1), i=1,n_cols)
  write(ounit) (evgxy(i, 1), i=1,n_cols)
  write(ounit) (evbxy(i, 1), i=1,n_cols)
  write(ounit) (ghvxy(i, 1), i=1,n_cols)
  write(ounit) (ghbxy(i, 1), i=1,n_cols)
  write(ounit) (irgxy(i, 1), i=1,n_cols)
  write(ounit) (ircxy(i, 1), i=1,n_cols)
  write(ounit) (irbxy(i, 1), i=1,n_cols)
  write(ounit) (trxy(i, 1), i=1,n_cols)
  write(ounit) (evcxy(i, 1), i=1,n_cols)
  write(ounit) (chleafxy(i, 1), i=1,n_cols)
  write(ounit) (chucxy(i, 1), i=1,n_cols)
  write(ounit) (chv2xy(i, 1), i=1,n_cols)
  write(ounit) (chb2xy(i, 1), i=1,n_cols)
  write(ounit) (pahxy(i, 1), i=1,n_cols)
  write(ounit) (pahgxy(i, 1), i=1,n_cols)
  write(ounit) (pahvxy(i, 1), i=1,n_cols)
  write(ounit) (pahbxy(i, 1), i=1,n_cols)
  write(ounit) (qintsxy(i, 1), i=1,n_cols)
  write(ounit) (qintrxy(i, 1), i=1,n_cols)
  write(ounit) (qdripsxy(i, 1), i=1,n_cols)
  write(ounit) (qdriprxy(i, 1), i=1,n_cols)
  write(ounit) (qthrosxy(i, 1), i=1,n_cols)
  write(ounit) (qthrorxy(i, 1), i=1,n_cols)
  write(ounit) (qsnsubxy(i, 1), i=1,n_cols)
  write(ounit) (qsnfroxy(i, 1), i=1,n_cols)
  write(ounit) (qsubcxy(i, 1), i=1,n_cols)
  write(ounit) (qfrocxy(i, 1), i=1,n_cols)
  write(ounit) (qevacxy(i, 1), i=1,n_cols)
  write(ounit) (qdewcxy(i, 1), i=1,n_cols)
  write(ounit) (qfrzcxy(i, 1), i=1,n_cols)
  write(ounit) (qmeltcxy(i, 1), i=1,n_cols)
  write(ounit) (qsnbotxy(i, 1), i=1,n_cols)
  write(ounit) (qmeltxy(i, 1), i=1,n_cols)
  write(ounit) (pondingxy(i, 1), i=1,n_cols)
  write(ounit) (fpicexy(i, 1), i=1,n_cols)
  write(ounit) (rainlsm(i, 1), i=1,n_cols)
  write(ounit) (snowlsm(i, 1), i=1,n_cols)
  write(ounit) (forctlsm(i, 1), i=1,n_cols)
  write(ounit) (forcqlsm(i, 1), i=1,n_cols)
  write(ounit) (forcplsm(i, 1), i=1,n_cols)
  write(ounit) (forczlsm(i, 1), i=1,n_cols)
  write(ounit) (forcwlsm(i, 1), i=1,n_cols)
  write(ounit) (eflxbxy(i, 1), i=1,n_cols)
  write(ounit) (soilenergy(i, 1), i=1,n_cols)
  write(ounit) (snowenergy(i, 1), i=1,n_cols)
  write(ounit) (canhsxy(i, 1), i=1,n_cols)

  ! 11 3D arrays
  write(ounit) ((smois(i, k, 1), i=1,n_cols), k=1,nsoil)
  write(ounit) ((sh2o(i, k, 1), i=1,n_cols), k=1,nsoil)
  write(ounit) ((tslb(i, k, 1), i=1,n_cols), k=1,nsoil)
  write(ounit) ((smoiseq(i, k, 1), i=1,n_cols), k=1,nsoil)
  write(ounit) ((acc_etranixy(i, k, 1), i=1,n_cols), k=1,nsoil)
  write(ounit) ((tsnoxy(i, k, 1), i=1,n_cols), k=-2,0)
  write(ounit) ((snicexy(i, k, 1), i=1,n_cols), k=-2,0)
  write(ounit) ((snliqxy(i, k, 1), i=1,n_cols), k=-2,0)
  write(ounit) ((zsnsoxy(i, k, 1), i=1,n_cols), k=-2,nsoil)
  write(ounit) ((albsoildirxy(i, k, 1), i=1,n_cols), k=1,2)
  write(ounit) ((albsoildifxy(i, k, 1), i=1,n_cols), k=1,2)

  close(ounit)
end program driver_noahmp_fortran
