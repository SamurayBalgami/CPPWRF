PROGRAM driver_thompson
  USE module_mp_thompson, ONLY: thompson_init, mp_gt_driver
  IMPLICIT NONE

  INTEGER :: u_in, u_out
  INTEGER :: ncols, nlevels
  REAL*8  :: dt_d

  REAL*8, ALLOCATABLE :: th(:,:,:), qv(:,:,:), qc(:,:,:), qr(:,:,:)
  REAL*8, ALLOCATABLE :: qi(:,:,:), qs(:,:,:), qg(:,:,:), ni(:,:,:), nr(:,:,:)
  REAL*8, ALLOCATABLE :: pii(:,:,:), p(:,:,:), w(:,:,:), dz(:,:,:)
  REAL*8, ALLOCATABLE :: hgt(:,:,:)
  REAL*8, ALLOCATABLE :: refl_10cm(:,:,:), re_cloud(:,:,:), re_ice(:,:,:), re_snow(:,:,:)

  REAL*8, ALLOCATABLE :: rainnc(:,:), rainncv(:,:)
  REAL*8, ALLOCATABLE :: snownc(:,:), snowncv(:,:)
  REAL*8, ALLOCATABLE :: graupelnc(:,:), graupelncv(:,:)
  REAL*8, ALLOCATABLE :: sr(:,:)

  CHARACTER(LEN=256) :: infile, outfile

  infile = "input_thompson.bin"
  outfile = "output_fortran.bin"

  IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
     CALL GET_COMMAND_ARGUMENT(1, infile)
  END IF
  IF (COMMAND_ARGUMENT_COUNT() >= 2) THEN
     CALL GET_COMMAND_ARGUMENT(2, outfile)
  END IF

  OPEN(NEWUNIT=u_in, FILE=TRIM(infile), FORM="unformatted", ACCESS="stream", STATUS="old")
  READ(u_in) ncols, nlevels, dt_d

  ! Allocate 3D arrays: (ncols, nlevels, 1)
  ALLOCATE(th(ncols, nlevels, 1), qv(ncols, nlevels, 1), qc(ncols, nlevels, 1))
  ALLOCATE(qr(ncols, nlevels, 1), qi(ncols, nlevels, 1), qs(ncols, nlevels, 1))
  ALLOCATE(qg(ncols, nlevels, 1), ni(ncols, nlevels, 1), nr(ncols, nlevels, 1))
  ALLOCATE(pii(ncols, nlevels, 1), p(ncols, nlevels, 1), w(ncols, nlevels, 1), dz(ncols, nlevels, 1))
  ALLOCATE(hgt(ncols, nlevels, 1))
  ALLOCATE(refl_10cm(ncols, nlevels, 1))
  ALLOCATE(re_cloud(ncols, nlevels, 1), re_ice(ncols, nlevels, 1), re_snow(ncols, nlevels, 1))

  ! Allocate 2D arrays: (ncols, 1)
  ALLOCATE(rainnc(ncols, 1), rainncv(ncols, 1))
  ALLOCATE(snownc(ncols, 1), snowncv(ncols, 1))
  ALLOCATE(graupelnc(ncols, 1), graupelncv(ncols, 1))
  ALLOCATE(sr(ncols, 1))

  refl_10cm = 0.0
  re_cloud = 0.0
  re_ice = 0.0
  re_snow = 0.0
  hgt = 100.0

  ! Read 3D inputs
  READ(u_in) th(:, 1:nlevels, 1)
  READ(u_in) qv(:, 1:nlevels, 1)
  READ(u_in) qc(:, 1:nlevels, 1)
  READ(u_in) qr(:, 1:nlevels, 1)
  READ(u_in) qi(:, 1:nlevels, 1)
  READ(u_in) qs(:, 1:nlevels, 1)
  READ(u_in) qg(:, 1:nlevels, 1)
  READ(u_in) ni(:, 1:nlevels, 1)
  READ(u_in) nr(:, 1:nlevels, 1)
  READ(u_in) pii(:, 1:nlevels, 1)
  READ(u_in) p(:, 1:nlevels, 1)
  READ(u_in) w(:, 1:nlevels, 1)
  READ(u_in) dz(:, 1:nlevels, 1)

  ! Read 2D inputs
  READ(u_in) rainnc(:, 1)
  READ(u_in) rainncv(:, 1)
  READ(u_in) snownc(:, 1)
  READ(u_in) snowncv(:, 1)
  READ(u_in) graupelnc(:, 1)
  READ(u_in) graupelncv(:, 1)
  READ(u_in) sr(:, 1)
  CLOSE(u_in)

  ! Initialize tables
  CALL thompson_init(hgt=hgt, &
                     ids=1, ide=ncols+1, jds=1, jde=2, kds=1, kde=nlevels, &
                     ims=1, ime=ncols, jms=1, jme=1, kms=1, kme=nlevels, &
                     its=1, ite=ncols, jts=1, jte=1, kts=1, kte=nlevels)

  ! Run Thompson microphysics driver
  CALL mp_gt_driver(qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, &
                    ni=ni, nr=nr, th=th, &
                    pii=pii, p=p, w=w, dz=dz, &
                    dt_in=REAL(dt_d, KIND=8), itimestep=1, &
                    RAINNC=rainnc, RAINNCV=rainncv, &
                    SNOWNC=snownc, SNOWNCV=snowncv, &
                    GRAUPELNC=graupelnc, GRAUPELNCV=graupelncv, &
                    SR=sr, &
                    refl_10cm=refl_10cm, &
                    re_cloud=re_cloud, re_ice=re_ice, re_snow=re_snow, &
                    has_reqc=0, has_reqi=0, has_reqs=0, &
                    ids=1, ide=ncols+1, jds=1, jde=2, kds=1, kde=nlevels, &
                    ims=1, ime=ncols, jms=1, jme=1, kms=1, kme=nlevels, &
                    its=1, ite=ncols, jts=1, jte=1, kts=1, kte=nlevels)

  ! Write output binary
  OPEN(NEWUNIT=u_out, FILE=TRIM(outfile), FORM="unformatted", ACCESS="stream", STATUS="replace")

  ! 3D outputs: qv, qc, qr, qi, qs, qg, ni, nr, th
  WRITE(u_out) qv(:, 1:nlevels, 1)
  WRITE(u_out) qc(:, 1:nlevels, 1)
  WRITE(u_out) qr(:, 1:nlevels, 1)
  WRITE(u_out) qi(:, 1:nlevels, 1)
  WRITE(u_out) qs(:, 1:nlevels, 1)
  WRITE(u_out) qg(:, 1:nlevels, 1)
  WRITE(u_out) ni(:, 1:nlevels, 1)
  WRITE(u_out) nr(:, 1:nlevels, 1)
  WRITE(u_out) th(:, 1:nlevels, 1)

  ! 2D outputs: rainnc, rainncv, snownc, snowncv, graupelnc, graupelncv, sr
  WRITE(u_out) rainnc(:, 1)
  WRITE(u_out) rainncv(:, 1)
  WRITE(u_out) snownc(:, 1)
  WRITE(u_out) snowncv(:, 1)
  WRITE(u_out) graupelnc(:, 1)
  WRITE(u_out) graupelncv(:, 1)
  WRITE(u_out) sr(:, 1)
  CLOSE(u_out)

  PRINT *, "Fortran Thompson driver completed successfully."
END PROGRAM driver_thompson
