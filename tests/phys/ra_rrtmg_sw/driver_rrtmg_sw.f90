PROGRAM driver_rrtmg_sw
  USE module_ra_rrtmg_sw
  IMPLICIT NONE

  INTEGER :: num_cols, nz
  REAL(8), PARAMETER :: g = 9.81_8
  REAL(8), PARAMETER :: r_d = 287.0_8
  REAL(8), PARAMETER :: solcon = 1367.0_8

  ! 2D Solar and Surface arrays: (num_cols, 1)
  REAL(8), ALLOCATABLE :: coszen(:,:), albedo(:,:), xland(:,:), xice(:,:), snow(:,:), xlat(:,:)
  REAL(8), ALLOCATABLE :: gsw(:,:), swcf(:,:)
  REAL(8), ALLOCATABLE :: swupt(:,:), swuptc(:,:), swuptcln(:,:)
  REAL(8), ALLOCATABLE :: swdnt(:,:), swdntc(:,:), swdntcln(:,:)
  REAL(8), ALLOCATABLE :: swupb(:,:), swupbc(:,:), swupbcln(:,:)
  REAL(8), ALLOCATABLE :: swdnb(:,:), swdnbc(:,:), swdnbcln(:,:)
  REAL(8), ALLOCATABLE :: swvisdir(:,:), swvisdif(:,:)
  REAL(8), ALLOCATABLE :: swnirdir(:,:), swnirdif(:,:)
  REAL(8), ALLOCATABLE :: coszr(:,:), xlong(:,:), tsk(:,:), obscur(:,:)
  REAL(8), ALLOCATABLE :: alswvisdir(:,:), alswvisdif(:,:), alswnirdir(:,:), alswnirdif(:,:)
  REAL(8), ALLOCATABLE :: swddir(:,:), swddni(:,:), swddif(:,:), swdownc(:,:), swddnic(:,:), swddirc(:,:)

  ! 3D Layer arrays: (num_cols, nz, 1)
  REAL(8), ALLOCATABLE :: p3d(:,:,:), pi3d(:,:,:), dz8w(:,:,:), t3d(:,:,:), rho3d(:,:,:)
  REAL(8), ALLOCATABLE :: qv3d(:,:,:), qc3d(:,:,:), qr3d(:,:,:), qi3d(:,:,:), qs3d(:,:,:), qg3d(:,:,:)
  REAL(8), ALLOCATABLE :: cldfra3d(:,:,:)
  REAL(8), ALLOCATABLE :: rthratensw(:,:,:), rthratenswc(:,:,:)
  REAL(8), ALLOCATABLE :: o33d(:,:,:), lradius(:,:,:), iradius(:,:,:)
  REAL(8), ALLOCATABLE :: re_cloud(:,:,:), re_ice(:,:,:), re_snow(:,:,:)
  REAL(8), ALLOCATABLE :: f_ice_phy(:,:,:), f_rain_phy(:,:,:), qndrop3d(:,:,:)
  REAL(8), ALLOCATABLE :: aerod(:,:,:,:)
  REAL(8), ALLOCATABLE :: tauaer300(:,:,:), tauaer400(:,:,:), tauaer600(:,:,:), tauaer999(:,:,:)
  REAL(8), ALLOCATABLE :: gaer300(:,:,:), gaer400(:,:,:), gaer600(:,:,:), gaer999(:,:,:)
  REAL(8), ALLOCATABLE :: waer300(:,:,:), waer400(:,:,:), waer600(:,:,:), waer999(:,:,:)
  REAL(8), POINTER :: tauaer3d_sw(:,:,:,:) => NULL()
  REAL(8), POINTER :: ssaaer3d_sw(:,:,:,:) => NULL()
  REAL(8), POINTER :: asyaer3d_sw(:,:,:,:) => NULL()

  ! 3D Interface arrays: (num_cols, nz+1, 1)
  REAL(8), ALLOCATABLE :: p8w(:,:,:), t8w(:,:,:)

  INTEGER :: i, k, ios, uout, uin
  CHARACTER(LEN=256) :: infile, outfile

  infile = "input_rrtmg_sw.bin"
  outfile = "output_fortran.bin"

  OPEN(NEWUNIT=uin, FILE=TRIM(infile), ACCESS='STREAM', FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN', STATUS='OLD', IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'Error opening input file: ', TRIM(infile)
    STOP 1
  END IF

  READ(uin, IOSTAT=ios) num_cols, nz
  IF (ios /= 0) THEN
    WRITE(*,*) 'Error reading header, iostat =', ios
    STOP 1
  END IF
  WRITE(*,*) 'Header: num_cols =', num_cols, ' nz =', nz

  ! Allocate arrays
  ALLOCATE(coszen(num_cols, 1), albedo(num_cols, 1), xland(num_cols, 1), &
           xice(num_cols, 1), snow(num_cols, 1), xlat(num_cols, 1))
  ALLOCATE(gsw(num_cols, 1), swcf(num_cols, 1))
  ALLOCATE(swupt(num_cols, 1), swuptc(num_cols, 1), swuptcln(num_cols, 1))
  ALLOCATE(swdnt(num_cols, 1), swdntc(num_cols, 1), swdntcln(num_cols, 1))
  ALLOCATE(swupb(num_cols, 1), swupbc(num_cols, 1), swupbcln(num_cols, 1))
  ALLOCATE(swdnb(num_cols, 1), swdnbc(num_cols, 1), swdnbcln(num_cols, 1))
  ALLOCATE(swvisdir(num_cols, 1), swvisdif(num_cols, 1))
  ALLOCATE(swnirdir(num_cols, 1), swnirdif(num_cols, 1))
  ALLOCATE(coszr(num_cols, 1), xlong(num_cols, 1), tsk(num_cols, 1), obscur(num_cols, 1))
  ALLOCATE(alswvisdir(num_cols, 1), alswvisdif(num_cols, 1), alswnirdir(num_cols, 1), alswnirdif(num_cols, 1))
  ALLOCATE(swddir(num_cols, 1), swddni(num_cols, 1), swddif(num_cols, 1), swdownc(num_cols, 1), swddnic(num_cols, 1), swddirc(num_cols, 1))

  ALLOCATE(p3d(num_cols, nz, 1), pi3d(num_cols, nz, 1), dz8w(num_cols, nz, 1), &
           t3d(num_cols, nz, 1), rho3d(num_cols, nz, 1))
  ALLOCATE(qv3d(num_cols, nz, 1), qc3d(num_cols, nz, 1), qr3d(num_cols, nz, 1), &
           qi3d(num_cols, nz, 1), qs3d(num_cols, nz, 1), qg3d(num_cols, nz, 1))
  ALLOCATE(cldfra3d(num_cols, nz, 1))
  ALLOCATE(rthratensw(num_cols, nz, 1), rthratenswc(num_cols, nz, 1))
  ALLOCATE(o33d(num_cols, nz, 1), lradius(num_cols, nz, 1), iradius(num_cols, nz, 1))
  ALLOCATE(re_cloud(num_cols, nz, 1), re_ice(num_cols, nz, 1), re_snow(num_cols, nz, 1))
  ALLOCATE(f_ice_phy(num_cols, nz, 1), f_rain_phy(num_cols, nz, 1), qndrop3d(num_cols, nz, 1))
  ALLOCATE(aerod(num_cols, nz, 1, 6))
  ALLOCATE(tauaer300(num_cols, nz, 1), tauaer400(num_cols, nz, 1), tauaer600(num_cols, nz, 1), tauaer999(num_cols, nz, 1))
  ALLOCATE(gaer300(num_cols, nz, 1), gaer400(num_cols, nz, 1), gaer600(num_cols, nz, 1), gaer999(num_cols, nz, 1))
  ALLOCATE(waer300(num_cols, nz, 1), waer400(num_cols, nz, 1), waer600(num_cols, nz, 1), waer999(num_cols, nz, 1))

  ALLOCATE(p8w(num_cols, nz+1, 1), t8w(num_cols, nz+1, 1))

  ! Read 2D arrays
  READ(uin) (coszen(i, 1), i=1,num_cols)
  READ(uin) (albedo(i, 1), i=1,num_cols)
  READ(uin) (xland(i, 1), i=1,num_cols)
  READ(uin) (xice(i, 1), i=1,num_cols)
  READ(uin) (snow(i, 1), i=1,num_cols)
  READ(uin) (xlat(i, 1), i=1,num_cols)

  ! Read 3D layer arrays: (nz, num_cols) in stream order
  DO k = 1, nz
    DO i = 1, num_cols
      READ(uin) p3d(i, k, 1)
    END DO
  END DO
  DO k = 1, nz
    DO i = 1, num_cols
      READ(uin) pi3d(i, k, 1)
    END DO
  END DO
  DO k = 1, nz
    DO i = 1, num_cols
      READ(uin) dz8w(i, k, 1)
    END DO
  END DO
  DO k = 1, nz
    DO i = 1, num_cols
      READ(uin) t3d(i, k, 1)
    END DO
  END DO
  DO k = 1, nz
    DO i = 1, num_cols
      READ(uin) rho3d(i, k, 1)
    END DO
  END DO
  DO k = 1, nz
    DO i = 1, num_cols
      READ(uin) qv3d(i, k, 1)
    END DO
  END DO
  DO k = 1, nz
    DO i = 1, num_cols
      READ(uin) qc3d(i, k, 1)
    END DO
  END DO
  DO k = 1, nz
    DO i = 1, num_cols
      READ(uin) qr3d(i, k, 1)
    END DO
  END DO
  DO k = 1, nz
    DO i = 1, num_cols
      READ(uin) qi3d(i, k, 1)
    END DO
  END DO
  DO k = 1, nz
    DO i = 1, num_cols
      READ(uin) qs3d(i, k, 1)
    END DO
  END DO
  DO k = 1, nz
    DO i = 1, num_cols
      READ(uin) qg3d(i, k, 1)
    END DO
  END DO
  DO k = 1, nz
    DO i = 1, num_cols
      READ(uin) cldfra3d(i, k, 1)
    END DO
  END DO

  ! Read 3D interface arrays: (nz+1, num_cols) in stream order
  DO k = 1, nz + 1
    DO i = 1, num_cols
      READ(uin) p8w(i, k, 1)
    END DO
  END DO
  DO k = 1, nz + 1
    DO i = 1, num_cols
      READ(uin) t8w(i, k, 1)
    END DO
  END DO

  CLOSE(uin)

  ! Initialize other input fields
  xlong(:, :) = 0.0_8
  tsk(:, :) = t8w(:, 1, :)
  obscur(:, :) = 0.0_8
  alswvisdir(:, :) = albedo(:, :)
  alswvisdif(:, :) = albedo(:, :)
  alswnirdir(:, :) = albedo(:, :)
  alswnirdif(:, :) = albedo(:, :)
  coszr(:, :) = coszen(:, :)

  o33d(:, :, :) = 0.0_8
  lradius(:, :, :) = 0.0_8
  iradius(:, :, :) = 0.0_8
  re_cloud(:, :, :) = 0.0_8
  re_ice(:, :, :) = 0.0_8
  re_snow(:, :, :) = 0.0_8
  f_ice_phy(:, :, :) = 0.0_8
  f_rain_phy(:, :, :) = 0.0_8
  qndrop3d(:, :, :) = 0.0_8
  aerod(:, :, :, :) = 0.0_8
  tauaer300(:, :, :) = 0.0_8
  tauaer400(:, :, :) = 0.0_8
  tauaer600(:, :, :) = 0.0_8
  tauaer999(:, :, :) = 0.0_8
  gaer300(:, :, :) = 0.0_8
  gaer400(:, :, :) = 0.0_8
  gaer600(:, :, :) = 0.0_8
  gaer999(:, :, :) = 0.0_8
  waer300(:, :, :) = 0.0_8
  waer400(:, :, :) = 0.0_8
  waer600(:, :, :) = 0.0_8
  waer999(:, :, :) = 0.0_8

  ! Initialize RRTMG-SW tables
  CALL rrtmg_swinit(.TRUE., 1, num_cols, 1, 1, 1, nz, 1, num_cols, 1, 1, 1, nz, 1, num_cols, 1, 1, 1, nz)

  ! Run RRTMG_SWRAD
  CALL rrtmg_swrad( &
      rthratensw=rthratensw, rthratenswc=rthratenswc, &
      swupt=swupt, swuptc=swuptc, swuptcln=swuptcln, &
      swdnt=swdnt, swdntc=swdntc, swdntcln=swdntcln, &
      swupb=swupb, swupbc=swupbc, swupbcln=swupbcln, &
      swdnb=swdnb, swdnbc=swdnbc, swdnbcln=swdnbcln, &
      swcf=swcf, gsw=gsw, &
      xtime=0.0_8, gmt=0.0_8, xlat=xlat, xlong=xlong, &
      radt=0.0_8, degrad=3.141592653589793_8/180.0_8, declin=0.0_8, &
      coszr=coszr, julday=1, solcon=solcon, &
      albedo=albedo, t3d=t3d, t8w=t8w, tsk=tsk, &
      p3d=p3d, p8w=p8w, pi3d=pi3d, rho3d=rho3d, &
      dz8w=dz8w, cldfra3d=cldfra3d, lradius=lradius, iradius=iradius, &
      is_cammgmp_used=.FALSE., r=r_d, g=g, &
      re_cloud=re_cloud, re_ice=re_ice, re_snow=re_snow, &
      has_reqc=0, has_reqi=0, has_reqs=0, &
      icloud=1, warm_rain=.FALSE., &
      cldovrlp=2, idcor=0, &
      f_ice_phy=f_ice_phy, f_rain_phy=f_rain_phy, &
      xland=xland, xice=xice, snow=snow, &
      qv3d=qv3d, qc3d=qc3d, qr3d=qr3d, &
      qi3d=qi3d, qs3d=qs3d, qg3d=qg3d, &
      o3input=0, o33d=o33d, &
      aer_opt=0, aerod=aerod, no_src=0, &
      alswvisdir=alswvisdir, alswvisdif=alswvisdif, &
      alswnirdir=alswnirdir, alswnirdif=alswnirdif, &
      swvisdir=swvisdir, swvisdif=swvisdif, &
      swnirdir=swnirdir, swnirdif=swnirdif, &
      sf_surface_physics=0, &
      f_qv=.TRUE., f_qc=.TRUE., f_qr=.TRUE., f_qi=.TRUE., f_qs=.TRUE., f_qg=.TRUE., &
      tauaer300=tauaer300, tauaer400=tauaer400, tauaer600=tauaer600, tauaer999=tauaer999, &
      gaer300=gaer300, gaer400=gaer400, gaer600=gaer600, gaer999=gaer999, &
      waer300=waer300, waer400=waer400, waer600=waer600, waer999=waer999, &
      aer_ra_feedback=0, progn=0, calc_clean_atm_diag=0, &
      qndrop3d=qndrop3d, f_qndrop=.FALSE., mp_physics=6, &
      ids=1, ide=num_cols, jds=1, jde=1, kds=1, kde=nz, &
      ims=1, ime=num_cols, jms=1, jme=1, kms=1, kme=nz, &
      its=1, ite=num_cols, jts=1, jte=1, kts=1, kte=nz, &
      tauaer3d_sw=tauaer3d_sw, ssaaer3d_sw=ssaaer3d_sw, asyaer3d_sw=asyaer3d_sw, &
      swddir=swddir, swddni=swddni, swddif=swddif, &
      swdownc=swdownc, swddnic=swddnic, swddirc=swddirc, &
      xcoszen=coszen, yr=2000, julian=1.0_8, ghg_input=0, &
      obscur=obscur, proceed_cmaq_sw=.FALSE. &
  )

  ! Write Fortran outputs
  OPEN(NEWUNIT=uout, FILE=TRIM(outfile), ACCESS='STREAM', FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN', STATUS='REPLACE', IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'Error creating output file: ', TRIM(outfile)
    STOP 1
  END IF

  WRITE(uout) (gsw(i, 1), i=1,num_cols)
  WRITE(uout) (swcf(i, 1), i=1,num_cols)
  WRITE(uout) (swupt(i, 1), i=1,num_cols)
  WRITE(uout) (swuptc(i, 1), i=1,num_cols)
  WRITE(uout) (swdnt(i, 1), i=1,num_cols)
  WRITE(uout) (swdntc(i, 1), i=1,num_cols)
  WRITE(uout) (swupb(i, 1), i=1,num_cols)
  WRITE(uout) (swupbc(i, 1), i=1,num_cols)
  WRITE(uout) (swdnb(i, 1), i=1,num_cols)
  WRITE(uout) (swdnbc(i, 1), i=1,num_cols)
  WRITE(uout) (swvisdir(i, 1), i=1,num_cols)
  WRITE(uout) (swvisdif(i, 1), i=1,num_cols)
  WRITE(uout) (swnirdir(i, 1), i=1,num_cols)
  WRITE(uout) (swnirdif(i, 1), i=1,num_cols)

  DO k = 1, nz
    DO i = 1, num_cols
      WRITE(uout) rthratensw(i, k, 1)
    END DO
  END DO
  DO k = 1, nz
    DO i = 1, num_cols
      WRITE(uout) rthratenswc(i, k, 1)
    END DO
  END DO

  CLOSE(uout)
  WRITE(*,*) 'Successfully ran Fortran RRTMG-SW and wrote ', TRIM(outfile)

END PROGRAM driver_rrtmg_sw
