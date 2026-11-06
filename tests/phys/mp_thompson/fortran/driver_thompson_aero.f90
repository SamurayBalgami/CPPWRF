! SPDX-License-Identifier: GPL-3.0-or-later
! mp_physics = 28 (aerosol aware, aer_init_opt 0, wif_input_opt 1): the inputs of driver_thompson.f90 followed by
! nc, nwfa, nifa (3-D) and nwfa2d, nifa2d (2-D); thompson_init with the aerosol arguments (CCN_ACTIVATE.BIN in the
! working directory; nwfa > eps: no profile), mp_gt_driver as microphysics_driver calls it for THOMPSONAERO; the
! outputs of driver_thompson.f90 followed by nc, nwfa, nifa, nbca. Built with WRF's BYTESWAPIO (as wrf.exe): the input
! and output files are opened CONVERT='native'
PROGRAM driver_thompson_aero
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
  REAL*8, ALLOCATABLE :: nc(:,:,:), nwfa(:,:,:), nifa(:,:,:), nbca(:,:,:), orho(:,:,:)
  REAL*8, ALLOCATABLE :: nwfa2d(:,:), nifa2d(:,:), nbca2d(:,:)

  CHARACTER(LEN=256) :: infile, outfile

  infile = "input_thompson.bin"
  outfile = "output_fortran.bin"

  IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
     CALL GET_COMMAND_ARGUMENT(1, infile)
  END IF
  IF (COMMAND_ARGUMENT_COUNT() >= 2) THEN
     CALL GET_COMMAND_ARGUMENT(2, outfile)
  END IF

  OPEN(NEWUNIT=u_in, FILE=TRIM(infile), FORM="unformatted", ACCESS="stream", STATUS="old", CONVERT="native")
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
  ALLOCATE(nc(ncols, nlevels, 1), nwfa(ncols, nlevels, 1), nifa(ncols, nlevels, 1), nbca(ncols, nlevels, 1))
  ALLOCATE(orho(ncols, nlevels, 1), nwfa2d(ncols, 1), nifa2d(ncols, 1), nbca2d(ncols, 1))
  nbca = 0.0
  orho = 1.0
  nbca2d = 0.0

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
  READ(u_in) nc(:, 1:nlevels, 1)
  READ(u_in) nwfa(:, 1:nlevels, 1)
  READ(u_in) nifa(:, 1:nlevels, 1)
  READ(u_in) nwfa2d(:, 1)
  READ(u_in) nifa2d(:, 1)
  CLOSE(u_in)

  ! Initialize tables. write_thompson_tables is .true. by default (Registry); no namelist is read here, so it is set
  ! as WRF's: thompson_init writes qr_acr_qg_V4.dat, qr_acr_qsV2.dat and freezeH2O.dat, which the C++ driver reads
  ! as the model reads a run's table files
  CALL nl_set_write_thompson_tables(1, .TRUE.)
  ! ite / jte as start_domain passes them (the staggered ends): MAXVAL(nwfa(its:ite-1,:,jts:jte-1)) covers the columns
  CALL thompson_init(hgt=hgt, orho=orho, nwfa2d=nwfa2d, nbca2d=nbca2d, &
                     nwfa=nwfa, nifa=nifa, nbca=nbca, wif_input_opt=1, frc_urb2d=nbca2d, &
                     dx=15000.0d0, dy=15000.0d0, is_start=.TRUE., &
                     ids=1, ide=ncols+1, jds=1, jde=2, kds=1, kde=nlevels, &
                     ims=1, ime=ncols, jms=1, jme=1, kms=1, kme=nlevels, &
                     its=1, ite=ncols+1, jts=1, jte=2, kts=1, kte=nlevels)

  ! Run Thompson microphysics driver
  CALL mp_gt_driver(qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, &
                    ni=ni, nr=nr, nc=nc, &
                    nwfa=nwfa, nifa=nifa, nbca=nbca, nwfa2d=nwfa2d, nifa2d=nifa2d, nbca2d=nbca2d, &
                    aer_init_opt=0, wif_input_opt=1, th=th, &
                    pii=pii, p=p, w=w, dz=dz, &
                    dt_in=REAL(dt_d, KIND=8), itimestep=1, &
                    RAINNC=rainnc, RAINNCV=rainncv, &
                    SNOWNC=snownc, SNOWNCV=snowncv, &
                    GRAUPELNC=graupelnc, GRAUPELNCV=graupelncv, &
                    SR=sr, &
                    refl_10cm=refl_10cm, diagflag=.TRUE., ke_diag=nlevels, do_radar_ref=1, &
                    re_cloud=re_cloud, re_ice=re_ice, re_snow=re_snow, &
                    has_reqc=1, has_reqi=1, has_reqs=1, &
                    ids=1, ide=ncols+1, jds=1, jde=2, kds=1, kde=nlevels, &
                    ims=1, ime=ncols, jms=1, jme=1, kms=1, kme=nlevels, &
                    its=1, ite=ncols, jts=1, jte=1, kts=1, kte=nlevels)

  ! Write output binary
  OPEN(NEWUNIT=u_out, FILE=TRIM(outfile), FORM="unformatted", ACCESS="stream", STATUS="replace", CONVERT="native")

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

  ! 3D outputs: effective radii (has_reqc = has_reqi = has_reqs = 1)
  WRITE(u_out) re_cloud(:, 1:nlevels, 1)
  WRITE(u_out) re_ice(:, 1:nlevels, 1)
  WRITE(u_out) re_snow(:, 1:nlevels, 1)
  WRITE(u_out) refl_10cm(:, 1:nlevels, 1)  ! diagflag = .true., do_radar_ref = 1
  WRITE(u_out) nc(:, 1:nlevels, 1)
  WRITE(u_out) nwfa(:, 1:nlevels, 1)
  WRITE(u_out) nifa(:, 1:nlevels, 1)
  WRITE(u_out) nbca(:, 1:nlevels, 1)
  CLOSE(u_out)

  PRINT *, "Fortran aerosol-aware Thompson driver completed successfully."
END PROGRAM driver_thompson_aero
