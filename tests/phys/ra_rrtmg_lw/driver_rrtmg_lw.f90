PROGRAM driver_rrtmg_lw
  USE module_ra_rrtmg_lw
  IMPLICIT NONE

  INTEGER :: num_cols, nz, magic
  REAL(8) :: p_top
  INTEGER :: icloud, cldovrlp, idcor, mp_physics
  REAL(8), PARAMETER :: g = 9.81_8
  REAL(8), PARAMETER :: r_d = 287.0_8

  ! 2D Surface arrays: (num_cols, 1)
  REAL(8), ALLOCATABLE :: tsk(:,:), emiss(:,:), xland(:,:), xice(:,:), snow(:,:), xlat(:,:)
  REAL(8), ALLOCATABLE :: glw(:,:), olr(:,:), lwcf(:,:)
  REAL(8), ALLOCATABLE :: lwupt(:,:), lwuptc(:,:), lwuptcln(:,:)
  REAL(8), ALLOCATABLE :: lwdnt(:,:), lwdntc(:,:), lwdntcln(:,:)
  REAL(8), ALLOCATABLE :: lwupb(:,:), lwupbc(:,:), lwupbcln(:,:)
  REAL(8), ALLOCATABLE :: lwdnb(:,:), lwdnbc(:,:), lwdnbcln(:,:)

  ! 3D Layer arrays: (num_cols, nz, 1)
  REAL(8), ALLOCATABLE :: p3d(:,:,:), pi3d(:,:,:), dz8w(:,:,:), t3d(:,:,:), rho3d(:,:,:)
  REAL(8), ALLOCATABLE :: qv3d(:,:,:), qc3d(:,:,:), qr3d(:,:,:), qi3d(:,:,:), qs3d(:,:,:), qg3d(:,:,:)
  REAL(8), ALLOCATABLE :: cldfra3d(:,:,:)
  REAL(8), ALLOCATABLE :: rthratenlw(:,:,:), rthratenlwc(:,:,:)
  REAL(8), ALLOCATABLE :: o33d(:,:,:)
  REAL(8), ALLOCATABLE :: re_cloud(:,:,:), re_ice(:,:,:), re_snow(:,:,:)
  REAL(8), ALLOCATABLE :: f_ice_phy(:,:,:), f_rain_phy(:,:,:)

  ! 3D Interface arrays: (num_cols, nz+1, 1)
  REAL(8), ALLOCATABLE :: p8w(:,:,:), t8w(:,:,:)

  INTEGER :: i, k, ios, uout, uin
  LOGICAL :: warm_rain, is_cammgmp_used
  LOGICAL :: f_qv, f_qc, f_qr, f_qi, f_qs, f_qg
  CHARACTER(LEN=256) :: infile, outfile

  infile = "input_rrtmg_lw.bin"
  outfile = "output_fortran.bin"

  OPEN(NEWUNIT=uin, FILE=TRIM(infile), ACCESS='STREAM', FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN', STATUS='OLD', IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'Error opening input file: ', TRIM(infile)
    STOP 1
  END IF

  ! Read header
  READ(uin, IOSTAT=ios) magic, num_cols, nz, p_top, icloud, cldovrlp, idcor, mp_physics
  IF (ios /= 0) THEN
    WRITE(*,*) 'Error reading header, iostat =', ios
    STOP 1
  END IF
  WRITE(*,*) 'Header: num_cols =', num_cols, ' nz =', nz, ' p_top =', p_top, ' icloud =', icloud

  ! Allocate arrays
  ALLOCATE(tsk(num_cols, 1), emiss(num_cols, 1), xland(num_cols, 1), &
           xice(num_cols, 1), snow(num_cols, 1), xlat(num_cols, 1))
  ALLOCATE(glw(num_cols, 1), olr(num_cols, 1), lwcf(num_cols, 1))
  ALLOCATE(lwupt(num_cols, 1), lwuptc(num_cols, 1), lwuptcln(num_cols, 1))
  ALLOCATE(lwdnt(num_cols, 1), lwdntc(num_cols, 1), lwdntcln(num_cols, 1))
  ALLOCATE(lwupb(num_cols, 1), lwupbc(num_cols, 1), lwupbcln(num_cols, 1))
  ALLOCATE(lwdnb(num_cols, 1), lwdnbc(num_cols, 1), lwdnbcln(num_cols, 1))

  ALLOCATE(p3d(num_cols, nz, 1), pi3d(num_cols, nz, 1), dz8w(num_cols, nz, 1), &
           t3d(num_cols, nz, 1), rho3d(num_cols, nz, 1))
  ALLOCATE(qv3d(num_cols, nz, 1), qc3d(num_cols, nz, 1), qr3d(num_cols, nz, 1), &
           qi3d(num_cols, nz, 1), qs3d(num_cols, nz, 1), qg3d(num_cols, nz, 1))
  ALLOCATE(cldfra3d(num_cols, nz, 1))
  ALLOCATE(rthratenlw(num_cols, nz, 1), rthratenlwc(num_cols, nz, 1))
  ALLOCATE(o33d(num_cols, nz, 1))
  ALLOCATE(re_cloud(num_cols, nz, 1), re_ice(num_cols, nz, 1), re_snow(num_cols, nz, 1))
  ALLOCATE(f_ice_phy(num_cols, nz, 1), f_rain_phy(num_cols, nz, 1))

  ALLOCATE(p8w(num_cols, nz+1, 1), t8w(num_cols, nz+1, 1))

  ! Read 2D arrays
  READ(uin) (tsk(i, 1), i=1,num_cols)
  READ(uin) (emiss(i, 1), i=1,num_cols)
  READ(uin) (xland(i, 1), i=1,num_cols)
  READ(uin) (xice(i, 1), i=1,num_cols)
  READ(uin) (snow(i, 1), i=1,num_cols)
  READ(uin) (xlat(i, 1), i=1,num_cols)

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

  CLOSE(uin)

  ! Initialize output arrays to zero
  rthratenlw = 0.0_8
  rthratenlwc = 0.0_8
  glw = 0.0_8
  olr = 0.0_8
  lwcf = 0.0_8
  lwupt = 0.0_8
  lwuptc = 0.0_8
  lwuptcln = 0.0_8
  lwdnt = 0.0_8
  lwdntc = 0.0_8
  lwdntcln = 0.0_8
  lwupb = 0.0_8
  lwupbc = 0.0_8
  lwupbcln = 0.0_8
  lwdnb = 0.0_8
  lwdnbc = 0.0_8
  lwdnbcln = 0.0_8

  o33d = 0.0_8
  re_cloud = 10.0e-6_8
  re_ice = 10.0e-6_8
  re_snow = 10.0e-6_8
  f_ice_phy = 0.0_8
  f_rain_phy = 0.0_8

  warm_rain = .FALSE.
  is_cammgmp_used = .FALSE.
  f_qv = .TRUE.
  f_qc = .TRUE.
  f_qr = .TRUE.
  f_qi = .TRUE.
  f_qs = .TRUE.
  f_qg = .TRUE.

  ! Call RRTMG LW initialization
  CALL RRTMG_LWINIT(p_top=REAL(p_top, 8), allowed_to_read=.TRUE., &
                    ids=1, ide=num_cols, jds=1, jde=1, kds=1, kde=nz, &
                    ims=1, ime=num_cols, jms=1, jme=1, kms=1, kme=nz, &
                    its=1, ite=num_cols, jts=1, jte=1, kts=1, kte=nz)

  WRITE(*,*) 'tsk(1,1) =', tsk(1,1), ' p8w(1,1,1) =', p8w(1,1,1), ' p3d(1,1,1) =', p3d(1,1,1)
  WRITE(*,*) 't3d(1,1,1) =', t3d(1,1,1), ' qv3d(1,1,1) =', qv3d(1,1,1)
  WRITE(*,*) 'Calling RRTMG_LWRAD...'
  ! Execute RRTMG LWRAD
  CALL RRTMG_LWRAD(                                                     &
       RTHRATENLW=rthratenlw,                                           &
       RTHRATENLWC=rthratenlwc,                                         &
       LWUPT=lwupt, LWUPTC=lwuptc, LWUPTCLN=lwuptcln,                   &
       LWDNT=lwdnt, LWDNTC=lwdntc, LWDNTCLN=lwdntcln,                   &
       LWUPB=lwupb, LWUPBC=lwupbc, LWUPBCLN=lwupbcln,                   &
       LWDNB=lwdnb, LWDNBC=lwdnbc, LWDNBCLN=lwdnbcln,                   &
       GLW=glw, OLR=olr, LWCF=lwcf, EMISS=emiss,                        &
       P8W=p8w, P3D=p3d, PI3D=pi3d, DZ8W=dz8w, TSK=tsk,                 &
       T3D=t3d, T8W=t8w, RHO3D=rho3d, R=r_d, G=g,                       &
       ICLOUD=icloud, WARM_RAIN=warm_rain, CLDFRA3D=cldfra3d,           &
       cldovrlp=cldovrlp, idcor=idcor, XLAT=xlat,                       &
       IS_CAMMGMP_USED=is_cammgmp_used,                                 &
       F_ICE_PHY=f_ice_phy, F_RAIN_PHY=f_rain_phy,                       &
       XLAND=xland, XICE=xice, SNOW=snow,                               &
       QV3D=qv3d, QC3D=qc3d, QR3D=qr3d,                                 &
       QI3D=qi3d, QS3D=qs3d, QG3D=qg3d,                                 &
       O3INPUT=0, O33D=o33d,                                            &
       F_QV=f_qv, F_QC=f_qc, F_QR=f_qr,                                 &
       F_QI=f_qi, F_QS=f_qs, F_QG=f_qg,                                 &
       RE_CLOUD=re_cloud, RE_ICE=re_ice, RE_SNOW=re_snow,               &
       has_reqc=0, has_reqi=0, has_reqs=0,                              &
       calc_clean_atm_diag=0,                                           &
       YR=2000, JULIAN=1.0_8, GHG_INPUT=0,                              &
       MP_PHYSICS=mp_physics,                                           &
       IDS=1, IDE=num_cols, JDS=1, JDE=1, KDS=1, KDE=nz,               &
       IMS=1, IME=num_cols, JMS=1, JME=1, KMS=1, KME=nz,               &
       ITS=1, ITE=num_cols, JTS=1, JTE=1, KTS=1, KTE=nz                )

  WRITE(*,*) 'Finished RRTMG_LWRAD.'
  WRITE(*,*) 'Sample Col 1: GLW =', glw(1, 1), ' OLR =', olr(1, 1), ' LWCF =', lwcf(1, 1)

  ! Write output binary
  OPEN(NEWUNIT=uout, FILE=TRIM(outfile), ACCESS='STREAM', FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN', STATUS='REPLACE', IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'Error opening output file, iostat =', ios
    STOP 1
  END IF

  WRITE(uout) num_cols, nz

  ! 3D outputs: (nz, num_cols) in stream order
  DO k = 1, nz
    DO i = 1, num_cols
      WRITE(uout) rthratenlw(i, k, 1)
    END DO
  END DO
  DO k = 1, nz
    DO i = 1, num_cols
      WRITE(uout) rthratenlwc(i, k, 1)
    END DO
  END DO

  ! 2D outputs: (num_cols)
  DO i = 1, num_cols
    WRITE(uout) glw(i, 1)
  END DO
  DO i = 1, num_cols
    WRITE(uout) olr(i, 1)
  END DO
  DO i = 1, num_cols
    WRITE(uout) lwcf(i, 1)
  END DO
  DO i = 1, num_cols
    WRITE(uout) lwupt(i, 1)
  END DO
  DO i = 1, num_cols
    WRITE(uout) lwuptc(i, 1)
  END DO
  DO i = 1, num_cols
    WRITE(uout) lwdnt(i, 1)
  END DO
  DO i = 1, num_cols
    WRITE(uout) lwdntc(i, 1)
  END DO
  DO i = 1, num_cols
    WRITE(uout) lwupb(i, 1)
  END DO
  DO i = 1, num_cols
    WRITE(uout) lwupbc(i, 1)
  END DO
  DO i = 1, num_cols
    WRITE(uout) lwdnb(i, 1)
  END DO
  DO i = 1, num_cols
    WRITE(uout) lwdnbc(i, 1)
  END DO

  CLOSE(uout)
  WRITE(*,*) 'Fortran RRTMG LW run completed successfully.'
END PROGRAM driver_rrtmg_lw
