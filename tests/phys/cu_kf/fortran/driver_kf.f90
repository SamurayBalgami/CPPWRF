PROGRAM driver_kf
  USE module_cu_kfeta, ONLY: kf_eta_cps, kf_lutab
  IMPLICIT NONE

  INTEGER :: u_in, u_out
  INTEGER :: ncols, nlevels
  REAL*8 :: dt_d, dx_d, cudt_d
  INTEGER :: stepcu, ktau, trigger, warm_rain_int
  INTEGER :: flag_qi_int, flag_qs_int, flag_qr_int
  LOGICAL :: warm_rain_l, f_qi_l, f_qs_l, f_qr_l, f_qv_l, f_qc_l
  LOGICAL :: adapt_step_flag

  REAL*8, ALLOCATABLE :: u(:,:,:), v(:,:,:), th(:,:,:), t(:,:,:)
  REAL*8, ALLOCATABLE :: qv(:,:,:), dz8w(:,:,:), pcps(:,:,:), rho(:,:,:), pi(:,:,:)
  REAL*8, ALLOCATABLE :: w(:,:,:)
  REAL*8, ALLOCATABLE :: w0avg(:,:,:)
  REAL*8, ALLOCATABLE :: raincv(:,:), pratec(:,:), nca(:,:)
  REAL*8, ALLOCATABLE :: cutop(:,:), cubot(:,:)
  LOGICAL, ALLOCATABLE :: cu_act_flag(:,:)

  REAL*8, ALLOCATABLE :: rthcuten(:,:,:), rqvcuten(:,:,:), rqccuten(:,:,:)
  REAL*8, ALLOCATABLE :: rqrcuten(:,:,:), rqicuten(:,:,:), rqscuten(:,:,:)
  REAL*8, ALLOCATABLE :: cldfra_dp_kf(:,:,:), cldfra_sh_kf(:,:,:)
  REAL*8, ALLOCATABLE :: qc_kf(:,:,:), qi_kf(:,:,:)

  ! Dummy arrays for optional edrates
  REAL*8, ALLOCATABLE :: udr_kf(:,:,:), ddr_kf(:,:,:), uer_kf(:,:,:), der_kf(:,:,:)
  REAL*8, ALLOCATABLE :: timec_kf(:,:)
  REAL*8, ALLOCATABLE :: rqvften(:,:,:)
  REAL*8, ALLOCATABLE :: shall(:,:)

  REAL*8 :: xlv0, xlv1, xls0, xls1, cp, r_d, g, ep1, ep2
  REAL*8 :: svp1, svp2, svp3, svpt0
  CHARACTER(LEN=256) :: infile, outfile

  infile = "input_kf.bin"
  outfile = "output_fortran.bin"

  OPEN(NEWUNIT=u_in, FILE=TRIM(infile), FORM="unformatted", ACCESS="stream", STATUS="old")
  READ(u_in) ncols, nlevels, dt_d, dx_d, cudt_d, stepcu, ktau, trigger, &
             warm_rain_int, flag_qi_int, flag_qs_int, flag_qr_int

  warm_rain_l = (warm_rain_int /= 0)
  f_qi_l = (flag_qi_int /= 0)
  f_qs_l = (flag_qs_int /= 0)
  f_qr_l = (flag_qr_int /= 0)
  f_qv_l = .TRUE.
  f_qc_l = .TRUE.
  adapt_step_flag = .FALSE.

  ! Allocate 3D arrays with (ncols, nlevels + 1, 1)
  ALLOCATE(u(ncols, nlevels + 1, 1), v(ncols, nlevels + 1, 1))
  ALLOCATE(th(ncols, nlevels + 1, 1), t(ncols, nlevels + 1, 1))
  ALLOCATE(qv(ncols, nlevels + 1, 1), dz8w(ncols, nlevels + 1, 1))
  ALLOCATE(pcps(ncols, nlevels + 1, 1), rho(ncols, nlevels + 1, 1), pi(ncols, nlevels + 1, 1))
  ALLOCATE(w(ncols, nlevels + 1, 1))
  ALLOCATE(w0avg(ncols, nlevels + 1, 1))

  ALLOCATE(raincv(ncols, 1), pratec(ncols, 1), nca(ncols, 1))
  ALLOCATE(cutop(ncols, 1), cubot(ncols, 1))
  ALLOCATE(cu_act_flag(ncols, 1))

  ALLOCATE(rthcuten(ncols, nlevels + 1, 1), rqvcuten(ncols, nlevels + 1, 1), rqccuten(ncols, nlevels + 1, 1))
  ALLOCATE(rqrcuten(ncols, nlevels + 1, 1), rqicuten(ncols, nlevels + 1, 1), rqscuten(ncols, nlevels + 1, 1))
  ALLOCATE(cldfra_dp_kf(ncols, nlevels + 1, 1), cldfra_sh_kf(ncols, nlevels + 1, 1))
  ALLOCATE(qc_kf(ncols, nlevels + 1, 1), qi_kf(ncols, nlevels + 1, 1))

  ALLOCATE(udr_kf(ncols, nlevels + 1, 1), ddr_kf(ncols, nlevels + 1, 1))
  ALLOCATE(uer_kf(ncols, nlevels + 1, 1), der_kf(ncols, nlevels + 1, 1))
  ALLOCATE(timec_kf(ncols, 1), rqvften(ncols, nlevels + 1, 1), shall(ncols, 1))

  u = 0.0
  v = 0.0
  th = 0.0
  t = 0.0
  qv = 0.0
  dz8w = 0.0
  pcps = 0.0
  rho = 0.0
  pi = 0.0
  w = 0.0
  w0avg = 0.0

  ! Read 3D inputs
  READ(u_in) u(:, 1:nlevels, 1)
  READ(u_in) v(:, 1:nlevels, 1)
  READ(u_in) th(:, 1:nlevels, 1)
  READ(u_in) t(:, 1:nlevels, 1)
  READ(u_in) qv(:, 1:nlevels, 1)
  READ(u_in) dz8w(:, 1:nlevels, 1)
  READ(u_in) pcps(:, 1:nlevels, 1)
  READ(u_in) rho(:, 1:nlevels, 1)
  READ(u_in) pi(:, 1:nlevels, 1)
  READ(u_in) w0avg(:, 1:nlevels, 1)

  ! Read staggered w: nlevels + 1
  READ(u_in) w(:, 1:nlevels+1, 1)

  ! Read 2D inputs
  READ(u_in) raincv(:, 1)
  READ(u_in) pratec(:, 1)
  READ(u_in) nca(:, 1)
  READ(u_in) cutop(:, 1)
  READ(u_in) cubot(:, 1)
  CLOSE(u_in)

  ! Physical constants
  xlv0 = 3.15E6
  xlv1 = 2370.0
  xls0 = 2.905E6
  xls1 = 259.532
  cp = 1004.5
  r_d = 287.0
  g = 9.81
  ep1 = 461.6 / 287.0 - 1.0
  ep2 = 287.0 / 461.6
  svp1 = 0.6112
  svp2 = 17.67
  svp3 = 29.65
  svpt0 = 273.15

  ! Initialize lookup tables
  CALL kf_lutab(svp1, svp2, svp3, svpt0)

  ! Initialize tendency arrays
  rthcuten = 0.0
  rqvcuten = 0.0
  rqccuten = 0.0
  rqrcuten = 0.0
  rqicuten = 0.0
  rqscuten = 0.0
  cldfra_dp_kf = 0.0
  cldfra_sh_kf = 0.0
  qc_kf = 0.0
  qi_kf = 0.0
  rqvften = 0.0
  udr_kf = 0.0
  ddr_kf = 0.0
  uer_kf = 0.0
  der_kf = 0.0
  timec_kf = 0.0

  CALL KF_eta_CPS(                                              &
              ids=1, ide=ncols, jds=1, jde=1, kds=1, kde=nlevels &
             ,ims=1, ime=ncols, jms=1, jme=1, kms=1, kme=nlevels+1 &
             ,its=1, ite=ncols, jts=1, jte=1, kts=1, kte=nlevels &
             ,trigger=trigger                                   &
             ,DT=dt_d, KTAU=ktau, DX=dx_d, CUDT=cudt_d          &
             ,ADAPT_STEP_FLAG=adapt_step_flag                   &
             ,rho=rho, RAINCV=raincv, PRATEC=pratec, NCA=nca    &
             ,U=u, V=v, TH=th, T=t, W=w, dz8w=dz8w              &
             ,Pcps=pcps, pi=pi                                  &
             ,W0AVG=w0avg, XLV0=xlv0, XLV1=xlv1, XLS0=xls0      &
             ,XLS1=xls1, CP=cp, R=r_d, G=g, EP1=ep1, EP2=ep2    &
             ,SVP1=svp1, SVP2=svp2, SVP3=svp3, SVPT0=svpt0      &
             ,STEPCU=stepcu, CU_ACT_FLAG=cu_act_flag            &
             ,warm_rain=warm_rain_l, CUTOP=cutop, CUBOT=cubot   &
             ,QV=qv                                             &
             ,shall=shall                                       &
             ,F_QV=f_qv_l, F_QC=f_qc_l, F_QR=f_qr_l             &
             ,F_QI=f_qi_l, F_QS=f_qs_l                          &
             ,RTHCUTEN=rthcuten, RQVCUTEN=rqvcuten              &
             ,RQCCUTEN=rqccuten, RQRCUTEN=rqrcuten              &
             ,RQICUTEN=rqicuten, RQSCUTEN=rqscuten              &
             ,RQVFTEN=rqvften                                   &
             ,cldfra_dp_KF=cldfra_dp_kf                         &
             ,cldfra_sh_KF=cldfra_sh_kf                         &
             ,qc_KF=qc_kf, qi_KF=qi_kf                          &
             ,UDR_KF=udr_kf, DDR_KF=ddr_kf                      &
             ,UER_KF=uer_kf, DER_KF=der_kf                      &
             ,TIMEC_KF=timec_kf, KF_EDRATES=0                   )

  ! Write output file
  OPEN(NEWUNIT=u_out, FILE=TRIM(outfile), FORM="unformatted", ACCESS="stream", STATUS="replace")
  ! 2D outputs
  WRITE(u_out) raincv(:, 1)
  WRITE(u_out) pratec(:, 1)
  WRITE(u_out) nca(:, 1)
  WRITE(u_out) cutop(:, 1)
  WRITE(u_out) cubot(:, 1)

  ! 3D outputs: size (ncols, nlevels)
  WRITE(u_out) w0avg(:, 1:nlevels, 1)
  WRITE(u_out) rthcuten(:, 1:nlevels, 1)
  WRITE(u_out) rqvcuten(:, 1:nlevels, 1)
  WRITE(u_out) rqccuten(:, 1:nlevels, 1)
  WRITE(u_out) rqrcuten(:, 1:nlevels, 1)
  WRITE(u_out) rqicuten(:, 1:nlevels, 1)
  WRITE(u_out) rqscuten(:, 1:nlevels, 1)
  WRITE(u_out) cldfra_dp_kf(:, 1:nlevels, 1)
  WRITE(u_out) cldfra_sh_kf(:, 1:nlevels, 1)
  WRITE(u_out) qc_kf(:, 1:nlevels, 1)
  WRITE(u_out) qi_kf(:, 1:nlevels, 1)
  CLOSE(u_out)

  PRINT *, "Fortran KF driver finished successfully for ", ncols, " columns."
END PROGRAM driver_kf
