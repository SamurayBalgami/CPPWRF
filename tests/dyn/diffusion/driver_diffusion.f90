! driver_diffusion.f90
! Differential testing harness for WRF v4.6.1 Diffusion and Damping kernels.
! Calls original WRF subroutines linked from libwrflib.a.

PROGRAM driver_diffusion
  USE module_configure, ONLY: grid_config_rec_type
  USE module_diffusion_em
  USE module_big_step_utilities_em
  IMPLICIT NONE

  INTEGER :: nx, ny, nz
  REAL(8) :: dx, dy, dt
  INTEGER :: ids, ide, jds, jde, kds, kde
  INTEGER :: ims, ime, jms, jme, kms, kme
  INTEGER :: its, ite, jts, jte, kts, kte
  INTEGER :: in_unit, out_unit, i, j, k, itf, jtf, ktf

  ! 1D arrays
  REAL(8), ALLOCATABLE :: dn(:), rdn(:), dnw(:), rdnw(:), fnm(:), fnp(:)
  REAL(8) :: cf1, cf2, cf3
  REAL(8), ALLOCATABLE :: c1_h(:), c2_h(:), c3_h(:), c4_h(:), c1f_h(:), c2f_h(:)
  REAL(8), ALLOCATABLE :: u_base(:), v_base(:), t_base(:), qv_base(:)

  ! 2D arrays
  REAL(8), ALLOCATABLE :: mut(:,:), muu(:,:), muv(:,:)
  REAL(8), ALLOCATABLE :: msftx(:,:), msfty(:,:), msfux(:,:), msfuy(:,:), msfvx(:,:), msfvy(:,:)
  REAL(8), ALLOCATABLE :: msfvx_inv(:,:)
  REAL(8), ALLOCATABLE :: ust(:,:), hfx(:,:), qfx(:,:)

  ! 3D input fields
  REAL(8), ALLOCATABLE :: phb(:,:,:), ph(:,:,:), u(:,:,:), v(:,:,:), w(:,:,:), ww(:,:,:)
  REAL(8), ALLOCATABLE :: thp(:,:,:), theta(:,:,:), t(:,:,:), pb(:,:,:), p(:,:,:), alt(:,:,:), rho(:,:,:)
  REAL(8), ALLOCATABLE :: moist(:,:,:,:), tke(:,:,:)
  REAL(8), ALLOCATABLE :: p8w(:,:,:), t8w(:,:,:)
  REAL(8), ALLOCATABLE :: dummy4d(:,:,:,:)
  REAL(8), ALLOCATABLE :: dummy_hpbl(:,:), dummy_dlk(:,:,:), dummy_meso(:,:,:)

  ! 3D outputs
  REAL(8), ALLOCATABLE :: z(:,:,:), rdz(:,:,:), rdzw(:,:,:), zx(:,:,:), zy(:,:,:)
  REAL(8), ALLOCATABLE :: defor11(:,:,:), defor22(:,:,:), defor33(:,:,:)
  REAL(8), ALLOCATABLE :: defor12(:,:,:), defor13(:,:,:), defor23(:,:,:)
  REAL(8), ALLOCATABLE :: div(:,:,:), bn2(:,:,:)
  REAL(8), ALLOCATABLE :: xkmh(:,:,:), xkmv(:,:,:), xkhh(:,:,:), xkhv(:,:,:)
  REAL(8), ALLOCATABLE :: xkmh_k1(:,:,:), xkmv_k1(:,:,:), xkhh_k1(:,:,:), xkhv_k1(:,:,:)
  REAL(8), ALLOCATABLE :: ru_tend(:,:,:), rv_tend(:,:,:), rw_tend(:,:,:), rt_tend(:,:,:)
  REAL(8), ALLOCATABLE :: moist_tend(:,:,:,:)
  REAL(8), ALLOCATABLE :: nba_rij(:,:,:,:), nba_mij(:,:,:,:)

  REAL(8) :: max_vert_cfl, max_horiz_cfl
  CHARACTER(LEN=256) :: in_path, out_path
  TYPE(grid_config_rec_type) :: config_flags

  IF (COMMAND_ARGUMENT_COUNT() .LT. 2) THEN
    in_path = "diffusion_input.bin"
    out_path = "diffusion_output_fortran.bin"
  ELSE
    CALL GET_COMMAND_ARGUMENT(1, in_path)
    CALL GET_COMMAND_ARGUMENT(2, out_path)
  END IF

  OPEN(NEWUNIT=in_unit, FILE=TRIM(in_path), ACCESS='STREAM', STATUS='OLD', FORM='UNFORMATTED')
  READ(in_unit) nx, ny, nz, dx, dy, dt

  ids = 1; ide = nx + 1; jds = 1; jde = ny + 1; kds = 1; kde = nz + 1
  ims = 0; ime = nx + 2; jms = 0; jme = ny + 2; kms = 1; kme = nz + 1
  its = 1; ite = nx; jts = 1; jte = ny; kts = 1; kte = kde

  ! Allocate 1D
  ALLOCATE(dn(nz), rdn(nz), dnw(nz), rdnw(nz), fnm(nz), fnp(nz))
  ALLOCATE(c1_h(nz), c2_h(nz), c3_h(nz), c4_h(nz), c1f_h(nz+1), c2f_h(nz+1))
  ALLOCATE(u_base(nz), v_base(nz), t_base(nz), qv_base(nz))

  READ(in_unit) dn
  READ(in_unit) rdn
  READ(in_unit) dnw
  READ(in_unit) rdnw
  READ(in_unit) fnm
  READ(in_unit) fnp
  READ(in_unit) cf1, cf2, cf3
  READ(in_unit) c1_h
  READ(in_unit) c2_h
  READ(in_unit) c3_h
  READ(in_unit) c4_h
  READ(in_unit) c1f_h
  READ(in_unit) c2f_h
  READ(in_unit) u_base
  READ(in_unit) v_base
  READ(in_unit) t_base
  READ(in_unit) qv_base

  ! Allocate 2D
  ALLOCATE(mut(ims:ime, jms:jme), muu(ims:ime, jms:jme), muv(ims:ime, jms:jme))
  ALLOCATE(msftx(ims:ime, jms:jme), msfty(ims:ime, jms:jme))
  ALLOCATE(msfux(ims:ime, jms:jme), msfuy(ims:ime, jms:jme))
  ALLOCATE(msfvx(ims:ime, jms:jme), msfvy(ims:ime, jms:jme), msfvx_inv(ims:ime, jms:jme))
  ALLOCATE(ust(ims:ime, jms:jme), hfx(ims:ime, jms:jme), qfx(ims:ime, jms:jme))

  mut = 0.0_8; muu = 0.0_8; muv = 0.0_8
  msftx = 1.0_8; msfty = 1.0_8; msfux = 1.0_8; msfuy = 1.0_8; msfvx = 1.0_8; msfvy = 1.0_8
  ust = 0.0_8; hfx = 0.0_8; qfx = 0.0_8

  READ(in_unit) mut(1:nx, 1:ny)
  READ(in_unit) muu(1:nx+1, 1:ny)
  READ(in_unit) muv(1:nx, 1:ny+1)
  READ(in_unit) msftx(1:nx, 1:ny)
  READ(in_unit) msfty(1:nx, 1:ny)
  READ(in_unit) msfux(1:nx+1, 1:ny)
  READ(in_unit) msfuy(1:nx+1, 1:ny)
  READ(in_unit) msfvx(1:nx, 1:ny+1)
  READ(in_unit) msfvy(1:nx, 1:ny+1)
  READ(in_unit) ust(1:nx, 1:ny)
  READ(in_unit) hfx(1:nx, 1:ny)
  READ(in_unit) qfx(1:nx, 1:ny)

  DO j = 1, ny+1
    DO i = 1, nx
      msfvx_inv(i, j) = 1.0_8 / msfvx(i, j)
    END DO
  END DO

  ! Allocate 3D fields
  ALLOCATE(phb(ims:ime, kms:kme, jms:jme), ph(ims:ime, kms:kme, jms:jme))
  ALLOCATE(u(ims:ime, kms:kme, jms:jme), v(ims:ime, kms:kme, jms:jme))
  ALLOCATE(w(ims:ime, kms:kme, jms:jme), ww(ims:ime, kms:kme, jms:jme))
  ALLOCATE(thp(ims:ime, kms:kme, jms:jme), theta(ims:ime, kms:kme, jms:jme), t(ims:ime, kms:kme, jms:jme))
  ALLOCATE(pb(ims:ime, kms:kme, jms:jme), p(ims:ime, kms:kme, jms:jme))
  ALLOCATE(alt(ims:ime, kms:kme, jms:jme), rho(ims:ime, kms:kme, jms:jme))
  ALLOCATE(moist(ims:ime, kms:kme, jms:jme, 2), tke(ims:ime, kms:kme, jms:jme))
  ALLOCATE(p8w(ims:ime, kms:kme, jms:jme), t8w(ims:ime, kms:kme, jms:jme))

  phb = 0.0_8; ph = 0.0_8; u = 0.0_8; v = 0.0_8; w = 0.0_8; ww = 0.0_8
  thp = 0.0_8; theta = 0.0_8; t = 0.0_8; pb = 0.0_8; p = 0.0_8; alt = 0.0_8; rho = 0.0_8
  moist = 0.0_8; tke = 0.0_8; p8w = 0.0_8; t8w = 0.0_8

  READ(in_unit) phb(1:nx, 1:nz+1, 1:ny)
  READ(in_unit) ph(1:nx, 1:nz+1, 1:ny)
  READ(in_unit) u(1:nx+1, 1:nz, 1:ny)
  READ(in_unit) v(1:nx, 1:nz, 1:ny+1)
  READ(in_unit) w(1:nx, 1:nz+1, 1:ny)
  READ(in_unit) ww(1:nx, 1:nz+1, 1:ny)
  READ(in_unit) thp(1:nx, 1:nz, 1:ny)
  READ(in_unit) theta(1:nx, 1:nz, 1:ny)
  READ(in_unit) t(1:nx, 1:nz, 1:ny)
  READ(in_unit) pb(1:nx, 1:nz, 1:ny)
  READ(in_unit) p(1:nx, 1:nz, 1:ny)
  READ(in_unit) alt(1:nx, 1:nz, 1:ny)
  READ(in_unit) rho(1:nx, 1:nz, 1:ny)
  READ(in_unit) moist(1:nx, 1:nz, 1:ny, 1)
  moist(:, :, :, 2) = moist(:, :, :, 1)
  READ(in_unit) tke(1:nx, 1:nz+1, 1:ny)
  CLOSE(in_unit)

  ! Compute p8w and t8w for calculate_N2
  DO j = 1, ny
    DO i = 1, nx
      p8w(i, 1, j) = pb(i, 1, j) + p(i, 1, j)
      t8w(i, 1, j) = theta(i, 1, j)
      DO k = 2, nz
        p8w(i, k, j) = 0.5_8 * (pb(i, k-1, j) + p(i, k-1, j) + pb(i, k, j) + p(i, k, j))
        t8w(i, k, j) = 0.5_8 * (theta(i, k-1, j) + theta(i, k, j))
      END DO
      p8w(i, nz+1, j) = pb(i, nz, j) + p(i, nz, j)
      t8w(i, nz+1, j) = theta(i, nz, j)
    END DO
  END DO

  ! Allocate output fields
  ALLOCATE(z(ims:ime, kms:kme, jms:jme), rdz(ims:ime, kms:kme, jms:jme), rdzw(ims:ime, kms:kme, jms:jme))
  ALLOCATE(zx(ims:ime, kms:kme, jms:jme), zy(ims:ime, kms:kme, jms:jme))
  ALLOCATE(defor11(ims:ime, kms:kme, jms:jme), defor22(ims:ime, kms:kme, jms:jme), defor33(ims:ime, kms:kme, jms:jme))
  ALLOCATE(defor12(ims:ime, kms:kme, jms:jme), defor13(ims:ime, kms:kme, jms:jme), defor23(ims:ime, kms:kme, jms:jme))
  ALLOCATE(div(ims:ime, kms:kme, jms:jme), bn2(ims:ime, kms:kme, jms:jme))
  ALLOCATE(xkmh(ims:ime, kms:kme, jms:jme), xkmv(ims:ime, kms:kme, jms:jme))
  ALLOCATE(xkhh(ims:ime, kms:kme, jms:jme), xkhv(ims:ime, kms:kme, jms:jme))
  ALLOCATE(xkmh_k1(ims:ime, kms:kme, jms:jme), xkmv_k1(ims:ime, kms:kme, jms:jme))
  ALLOCATE(xkhh_k1(ims:ime, kms:kme, jms:jme), xkhv_k1(ims:ime, kms:kme, jms:jme))
  ALLOCATE(ru_tend(ims:ime, kms:kme, jms:jme), rv_tend(ims:ime, kms:kme, jms:jme))
  ALLOCATE(rw_tend(ims:ime, kms:kme, jms:jme), rt_tend(ims:ime, kms:kme, jms:jme))
  ALLOCATE(moist_tend(ims:ime, kms:kme, jms:jme, 2))
  ALLOCATE(nba_rij(ims:ime, kms:kme, jms:jme, 1), nba_mij(ims:ime, kms:kme, jms:jme, 1))
  ALLOCATE(dummy4d(ims:ime, kms:kme, jms:jme, 1))
  ALLOCATE(dummy_hpbl(ims:ime, jms:jme), dummy_dlk(ims:ime, kms:kme, jms:jme), dummy_meso(ims:ime, kms:kme, jms:jme))

  dummy4d = 0.0_8; dummy_hpbl = 0.0_8; dummy_dlk = 0.0_8; dummy_meso = 0.0_8

  OPEN(NEWUNIT=out_unit, FILE=TRIM(out_path), ACCESS='STREAM', STATUS='REPLACE', FORM='UNFORMATTED')

  ! ============================================================================
  ! TEST 1: compute_diff_metrics
  ! Case A: Periodic boundary conditions
  ! ============================================================================
  config_flags%periodic_x = .TRUE.
  config_flags%periodic_y = .TRUE.
  config_flags%specified = .FALSE.
  config_flags%nested = .FALSE.
  config_flags%open_xs = .FALSE.
  config_flags%open_xe = .FALSE.
  config_flags%open_ys = .FALSE.
  config_flags%open_ye = .FALSE.

  z = 0.0_8; rdz = 0.0_8; rdzw = 0.0_8; zx = 0.0_8; zy = 0.0_8
  CALL fill_periodic_halos()
  CALL compute_diff_metrics(config_flags, ph, phb, z, rdz, rdzw, zx, zy, 1.0_8/dx, 1.0_8/dy, &
                            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kde)
  zx(nx+1, 1:nz, 1:ny) = zx(1, 1:nz, 1:ny)
  zy(1:nx, 1:nz, ny+1) = zy(1:nx, 1:nz, 1)

  CALL write_chunk(out_unit, "metrics_z_periodic", z(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "metrics_rdz_periodic", rdz(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "metrics_rdzw_periodic", rdzw(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "metrics_zx_periodic", zx(1:nx+1, 1:nz, 1:ny), (nx+1)*nz*ny)
  CALL write_chunk(out_unit, "metrics_zy_periodic", zy(1:nx, 1:nz, 1:ny+1), nx*nz*(ny+1))

  ! Case B: Specified boundary conditions
  config_flags%periodic_x = .FALSE.
  config_flags%periodic_y = .FALSE.
  config_flags%specified = .TRUE.

  z = 0.0_8; rdz = 0.0_8; rdzw = 0.0_8; zx = 0.0_8; zy = 0.0_8
  CALL compute_diff_metrics(config_flags, ph, phb, z, rdz, rdzw, zx, zy, 1.0_8/dx, 1.0_8/dy, &
                            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kde)
  zx(1, 1:nz, 1:ny) = 0.0_8
  zx(nx+1, 1:nz, 1:ny) = 0.0_8
  zy(1:nx, 1:nz, 1) = 0.0_8
  zy(1:nx, 1:nz, ny+1) = 0.0_8

  CALL write_chunk(out_unit, "metrics_z_specified", z(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "metrics_rdz_specified", rdz(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "metrics_rdzw_specified", rdzw(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "metrics_zx_specified", zx(1:nx+1, 1:nz, 1:ny), (nx+1)*nz*ny)
  CALL write_chunk(out_unit, "metrics_zy_specified", zy(1:nx, 1:nz, 1:ny+1), nx*nz*(ny+1))

  ! ============================================================================
  ! TEST 2: cal_deform_and_div
  ! Case A: Periodic
  ! ============================================================================
  config_flags%periodic_x = .TRUE.
  config_flags%periodic_y = .TRUE.
  config_flags%specified = .FALSE.
  config_flags%sfs_opt = 0

  defor11 = 0.0_8; defor22 = 0.0_8; defor33 = 0.0_8
  defor12 = 0.0_8; defor13 = 0.0_8; defor23 = 0.0_8; div = 0.0_8; nba_rij = 0.0_8

  CALL fill_periodic_halos()
  CALL cal_deform_and_div(config_flags, u, v, w, div, defor11, defor22, defor33, &
                          defor12, defor13, defor23, nba_rij, 1, u_base, v_base, &
                          msfux, msfuy, msfvx, msfvy, msftx, msfty, &
                          1.0_8/dx, 1.0_8/dy, dn, dnw, rdz, rdzw, fnm, fnp, cf1, cf2, cf3, zx, zy, &
                          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  defor12(nx+1, 1:nz, 1:ny+1) = defor12(1, 1:nz, 1:ny+1)
  defor12(1:nx+1, 1:nz, ny+1) = defor12(1:nx+1, 1:nz, 1)
  defor13(nx+1, 1:nz+1, 1:ny) = defor13(1, 1:nz+1, 1:ny)
  defor23(1:nx, 1:nz+1, ny+1) = defor23(1:nx, 1:nz+1, 1)

  CALL write_chunk(out_unit, "defor11_periodic", defor11(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "defor22_periodic", defor22(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "defor33_periodic", defor33(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "defor12_periodic", defor12(1:nx+1, 1:nz, 1:ny+1), (nx+1)*nz*(ny+1))
  CALL write_chunk(out_unit, "defor13_periodic", defor13(1:nx+1, 1:nz+1, 1:ny), (nx+1)*(nz+1)*ny)
  CALL write_chunk(out_unit, "defor23_periodic", defor23(1:nx, 1:nz+1, 1:ny+1), nx*(nz+1)*(ny+1))
  CALL write_chunk(out_unit, "div_periodic", div(1:nx, 1:nz, 1:ny), nx*nz*ny)

  ! Case B: Specified
  config_flags%periodic_x = .FALSE.
  config_flags%periodic_y = .FALSE.
  config_flags%specified = .TRUE.

  defor11 = 0.0_8; defor22 = 0.0_8; defor33 = 0.0_8
  defor12 = 0.0_8; defor13 = 0.0_8; defor23 = 0.0_8; div = 0.0_8; nba_rij = 0.0_8

  CALL cal_deform_and_div(config_flags, u, v, w, div, defor11, defor22, defor33, &
                          defor12, defor13, defor23, nba_rij, 1, u_base, v_base, &
                          msfux, msfuy, msfvx, msfvy, msftx, msfty, &
                          1.0_8/dx, 1.0_8/dy, dn, dnw, rdz, rdzw, fnm, fnp, cf1, cf2, cf3, zx, zy, &
                          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL write_chunk(out_unit, "defor11_specified", defor11(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "defor22_specified", defor22(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "defor33_specified", defor33(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "defor12_specified", defor12(1:nx+1, 1:nz, 1:ny+1), (nx+1)*nz*(ny+1))
  CALL write_chunk(out_unit, "defor13_specified", defor13(1:nx+1, 1:nz+1, 1:ny), (nx+1)*(nz+1)*ny)
  CALL write_chunk(out_unit, "defor23_specified", defor23(1:nx, 1:nz+1, 1:ny+1), nx*(nz+1)*(ny+1))
  CALL write_chunk(out_unit, "div_specified", div(1:nx, 1:nz, 1:ny), nx*nz*ny)

  ! ============================================================================
  ! TEST 3: calculate_km_kh
  ! Case A: km_opt = 1 (isotropic)
  ! ============================================================================
  config_flags%km_opt = 1
  config_flags%diff_opt = 1
  config_flags%c_s = 0.25_8
  xkmh = 0.0_8; xkmv = 0.0_8; xkhh = 0.0_8; xkhv = 0.0_8

  CALL calculate_km_kh(config_flags, dt, 0.0_8, 0.0_8, 0, xkmh, xkmv, xkhh, xkhv, bn2, &
                       12.0_8, 6.0_8, div, defor11, defor22, defor33, defor12, defor13, defor23, &
                       tke, p8w, t8w, theta, t, p, moist, dn, dnw, dx, dy, rdz, rdzw, 1, 1, &
                       cf1, cf2, cf3, .FALSE., 1000.0_8, msftx, msfty, zx, zy, &
                       dummy_hpbl, dummy_dlk, dummy_meso, &
                       ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL write_chunk(out_unit, "xkmh_kmopt1", xkmh(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "xkmv_kmopt1", xkmv(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "xkhh_kmopt1", xkhh(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "xkhv_kmopt1", xkhv(1:nx, 1:nz, 1:ny), nx*nz*ny)
  xkmh_k1 = xkmh; xkmv_k1 = xkmv; xkhh_k1 = xkhh; xkhv_k1 = xkhv

  ! Case B: km_opt = 4 (2D Smagorinsky with diff_opt = 2 slope reduction)
  config_flags%km_opt = 4
  config_flags%diff_opt = 2
  config_flags%c_s = 0.25_8
  xkmh = 0.0_8; xkmv = 0.0_8; xkhh = 0.0_8; xkhv = 0.0_8

  CALL calculate_km_kh(config_flags, dt, 0.0_8, 0.0_8, 0, xkmh, xkmv, xkhh, xkhv, bn2, &
                       12.0_8, 6.0_8, div, defor11, defor22, defor33, defor12, defor13, defor23, &
                       tke, p8w, t8w, theta, t, p, moist, dn, dnw, dx, dy, rdz, rdzw, 1, 1, &
                       cf1, cf2, cf3, .FALSE., 1000.0_8, msftx, msfty, zx, zy, &
                       dummy_hpbl, dummy_dlk, dummy_meso, &
                       ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL write_chunk(out_unit, "xkmh_kmopt4", xkmh(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "xkmv_kmopt4", xkmv(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "xkhh_kmopt4", xkhh(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "xkhv_kmopt4", xkhv(1:nx, 1:nz, 1:ny), nx*nz*ny)

  ! ============================================================================
  ! TEST 4: diff_opt = 1 (module_big_step_utilities_em)
  ! ============================================================================
  config_flags%diff_opt = 1
  config_flags%periodic_x = .FALSE.
  config_flags%periodic_y = .FALSE.
  config_flags%specified = .TRUE.

  ru_tend = 0.0_8; rv_tend = 0.0_8; rw_tend = 0.0_8; rt_tend = 0.0_8; moist_tend = 0.0_8

  CALL horizontal_diffusion('u', u, ru_tend, mut, c1_h, c2_h, config_flags, &
                            msfux, msfuy, msfvx, msfvx_inv, msfvy, msftx, msfty, &
                            10.0_8, xkmh, 1.0_8/dx, 1.0_8/dy, &
                            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL horizontal_diffusion('v', v, rv_tend, mut, c1_h, c2_h, config_flags, &
                            msfux, msfuy, msfvx, msfvx_inv, msfvy, msftx, msfty, &
                            10.0_8, xkmh, 1.0_8/dx, 1.0_8/dy, &
                            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL horizontal_diffusion('w', w, rw_tend, mut, c1f_h, c2f_h, config_flags, &
                            msfux, msfuy, msfvx, msfvx_inv, msfvy, msftx, msfty, &
                            10.0_8, xkmh, 1.0_8/dx, 1.0_8/dy, &
                            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL horizontal_diffusion('m', theta, rt_tend, mut, c1_h, c2_h, config_flags, &
                            msfux, msfuy, msfvx, msfvx_inv, msfvy, msftx, msfty, &
                            10.0_8, xkmh, 1.0_8/dx, 1.0_8/dy, &
                            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL horizontal_diffusion_3dmp('m', moist(:,:,:,1), moist_tend(:,:,:,1), mut, c1_h, c2_h, config_flags, &
                                 dummy4d(:,:,:,1), msfux, msfuy, msfvx, msfvx_inv, msfvy, msftx, msfty, &
                                 10.0_8, xkmh, 1.0_8/dx, 1.0_8/dy, &
                                 ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

  CALL write_chunk(out_unit, "u_hdiff1", ru_tend(1:nx+1, 1:nz, 1:ny), (nx+1)*nz*ny)
  CALL write_chunk(out_unit, "v_hdiff1", rv_tend(1:nx, 1:nz, 1:ny+1), nx*nz*(ny+1))
  CALL write_chunk(out_unit, "w_hdiff1", rw_tend(1:nx, 1:nz+1, 1:ny), nx*(nz+1)*ny)
  CALL write_chunk(out_unit, "t_hdiff1", rt_tend(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "m_hdiff1", moist_tend(1:nx, 1:nz, 1:ny, 1), nx*nz*ny)

  ! Vertical diffusion (diff_opt = 1)
  ru_tend = 0.0_8; rv_tend = 0.0_8; rw_tend = 0.0_8; rt_tend = 0.0_8; moist_tend = 0.0_8

  CALL vertical_diffusion_u(u, ru_tend, config_flags, u_base, c1_h, c2_h, alt, muu, rdn, rdnw, 5.0_8, &
                            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL vertical_diffusion_v(v, rv_tend, config_flags, v_base, c1_h, c2_h, alt, muv, rdn, rdnw, 5.0_8, &
                            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL vertical_diffusion('w', w, rw_tend, config_flags, c1f_h, c2f_h, alt, mut, rdn, rdnw, 5.0_8, &
                          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL vertical_diffusion('m', theta, rt_tend, config_flags, c1_h, c2_h, alt, mut, rdn, rdnw, 5.0_8, &
                          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL vertical_diffusion_3dmp(moist(:,:,:,1), moist_tend(:,:,:,1), config_flags, dummy4d(:,:,:,1), c1_h, c2_h, &
                               alt, mut, rdn, rdnw, 5.0_8, &
                               ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

  CALL write_chunk(out_unit, "u_vdiff1", ru_tend(1:nx+1, 1:nz, 1:ny), (nx+1)*nz*ny)
  CALL write_chunk(out_unit, "v_vdiff1", rv_tend(1:nx, 1:nz, 1:ny+1), nx*nz*(ny+1))
  CALL write_chunk(out_unit, "w_vdiff1", rw_tend(1:nx, 1:nz+1, 1:ny), nx*(nz+1)*ny)
  CALL write_chunk(out_unit, "t_vdiff1", rt_tend(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "m_vdiff1", moist_tend(1:nx, 1:nz, 1:ny, 1), nx*nz*ny)

  ! ============================================================================
  ! TEST 5: diff_opt = 2 (module_diffusion_em)
  ! ============================================================================
  config_flags%diff_opt = 2
  config_flags%km_opt = 4
  config_flags%isfflx = 0
  config_flags%tke_drag_coefficient = 0.0_8

  ru_tend = 0.0_8; rv_tend = 0.0_8; rw_tend = 0.0_8; rt_tend = 0.0_8; moist_tend = 0.0_8
  nba_mij = 0.0_8

  CALL horizontal_diffusion_2(rt_tend, ru_tend, rv_tend, rw_tend, tke(:,:,1:nz), &
                              moist_tend, 2, dummy4d, 0, dummy4d, 0, dummy4d, 0, &
                              thp, theta, tke, config_flags, defor11, defor22, defor12, defor13, defor23, &
                              nba_mij, 1, div, moist, dummy4d, dummy4d, dummy4d, &
                              msfux, msfuy, msfvx, msfvy, msftx, msfty, &
                              xkmh, xkmv, xkhh, config_flags%km_opt, 1.0_8/dx, 1.0_8/dy, rdz, rdzw, &
                              fnm, fnp, cf1, cf2, cf3, zx, zy, dn, dnw, rho, &
                              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

  CALL write_chunk(out_unit, "u_hdiff2", ru_tend(1:nx+1, 1:nz, 1:ny), (nx+1)*nz*ny)
  CALL write_chunk(out_unit, "v_hdiff2", rv_tend(1:nx, 1:nz, 1:ny+1), nx*nz*(ny+1))
  CALL write_chunk(out_unit, "w_hdiff2", rw_tend(1:nx, 1:nz+1, 1:ny), nx*(nz+1)*ny)
  CALL write_chunk(out_unit, "t_hdiff2", rt_tend(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "m_hdiff2", moist_tend(1:nx, 1:nz, 1:ny, 2), nx*nz*ny)

  ! Vertical diffusion 2
  ru_tend = 0.0_8; rv_tend = 0.0_8; rw_tend = 0.0_8; rt_tend = 0.0_8; moist_tend = 0.0_8
  nba_mij = 0.0_8

  CALL vertical_diffusion_2(ru_tend, rv_tend, rw_tend, rt_tend, tke(:,:,1:nz), &
                            moist_tend, 2, dummy4d, 0, dummy4d, 0, dummy4d, 0, &
                            u, v, thp, u_base, v_base, t_base, qv_base, tke, theta, &
                            config_flags, defor13, defor23, defor33, nba_mij, 1, div, &
                            moist, dummy4d, dummy4d, dummy4d, &
                            xkmv, xkhv, xkmh, config_flags%km_opt, fnm, fnp, dn, dnw, rdz, rdzw, &
                            hfx, qfx, ust, rho, &
                            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

  CALL write_chunk(out_unit, "u_vdiff2", ru_tend(1:nx+1, 1:nz, 1:ny), (nx+1)*nz*ny)
  CALL write_chunk(out_unit, "v_vdiff2", rv_tend(1:nx, 1:nz, 1:ny+1), nx*nz*(ny+1))
  CALL write_chunk(out_unit, "w_vdiff2", rw_tend(1:nx, 1:nz+1, 1:ny), nx*(nz+1)*ny)
  CALL write_chunk(out_unit, "t_vdiff2", rt_tend(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "m_vdiff2", moist_tend(1:nx, 1:nz, 1:ny, 2), nx*nz*ny)

  ! Case C: Vertical diffusion 2 with km_opt = 1 (nonzero xkhv to exercise theta and scalar vertical terms)
  config_flags%km_opt = 1
  ru_tend = 0.0_8; rv_tend = 0.0_8; rw_tend = 0.0_8; rt_tend = 0.0_8; moist_tend = 0.0_8
  nba_mij = 0.0_8

  CALL vertical_diffusion_2(ru_tend, rv_tend, rw_tend, rt_tend, tke(:,:,1:nz), &
                            moist_tend, 2, dummy4d, 0, dummy4d, 0, dummy4d, 0, &
                            u, v, thp, u_base, v_base, t_base, qv_base, tke, theta, &
                            config_flags, defor13, defor23, defor33, nba_mij, 1, div, &
                            moist, dummy4d, dummy4d, dummy4d, &
                            xkmv_k1, xkhv_k1, xkmh_k1, config_flags%km_opt, fnm, fnp, dn, dnw, rdz, rdzw, &
                            hfx, qfx, ust, rho, &
                            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

  CALL write_chunk(out_unit, "u_vdiff2_kmopt1", ru_tend(1:nx+1, 1:nz, 1:ny), (nx+1)*nz*ny)
  CALL write_chunk(out_unit, "v_vdiff2_kmopt1", rv_tend(1:nx, 1:nz, 1:ny+1), nx*nz*(ny+1))
  CALL write_chunk(out_unit, "w_vdiff2_kmopt1", rw_tend(1:nx, 1:nz+1, 1:ny), nx*(nz+1)*ny)
  CALL write_chunk(out_unit, "t_vdiff2_kmopt1", rt_tend(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "m_vdiff2_kmopt1", moist_tend(1:nx, 1:nz, 1:ny, 2), nx*nz*ny)

  ! ============================================================================
  ! TEST 6: sixth_order_diffusion
  ! Case A: diff_6th_opt = 1 (linear)
  ! ============================================================================
  ru_tend = 0.0_8; rv_tend = 0.0_8; rw_tend = 0.0_8; rt_tend = 0.0_8; moist_tend = 0.0_8
  config_flags%diff_6th_slopeopt = 0

  CALL sixth_order_diffusion('u', u, ru_tend, mut, dt, config_flags, c1_h, c2_h, &
                             1, 0.12_8, phb, ph, 1.0_8/dx, 1.0_8/dy, &
                             msftx, msfty, msfux, msfuy, msfvx, msfvy, &
                             ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL sixth_order_diffusion('v', v, rv_tend, mut, dt, config_flags, c1_h, c2_h, &
                             1, 0.12_8, phb, ph, 1.0_8/dx, 1.0_8/dy, &
                             msftx, msfty, msfux, msfuy, msfvx, msfvy, &
                             ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL sixth_order_diffusion('w', w, rw_tend, mut, dt, config_flags, c1f_h, c2f_h, &
                             1, 0.12_8, phb, ph, 1.0_8/dx, 1.0_8/dy, &
                             msftx, msfty, msfux, msfuy, msfvx, msfvy, &
                             ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL sixth_order_diffusion('m', theta, rt_tend, mut, dt, config_flags, c1_h, c2_h, &
                             1, 0.12_8, phb, ph, 1.0_8/dx, 1.0_8/dy, &
                             msftx, msfty, msfux, msfuy, msfvx, msfvy, &
                             ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL sixth_order_diffusion('m', moist(:,:,:,1), moist_tend(:,:,:,1), mut, dt, config_flags, c1_h, c2_h, &
                             1, 0.12_8, phb, ph, 1.0_8/dx, 1.0_8/dy, &
                             msftx, msfty, msfux, msfuy, msfvx, msfvy, &
                             ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

  CALL write_chunk(out_unit, "u_6th_opt1", ru_tend(1:nx+1, 1:nz, 1:ny), (nx+1)*nz*ny)
  CALL write_chunk(out_unit, "v_6th_opt1", rv_tend(1:nx, 1:nz, 1:ny+1), nx*nz*(ny+1))
  CALL write_chunk(out_unit, "w_6th_opt1", rw_tend(1:nx, 1:nz+1, 1:ny), nx*(nz+1)*ny)
  CALL write_chunk(out_unit, "t_6th_opt1", rt_tend(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "m_6th_opt1", moist_tend(1:nx, 1:nz, 1:ny, 1), nx*nz*ny)

  ! Case B: diff_6th_opt = 2 (monotonic)
  ru_tend = 0.0_8; rv_tend = 0.0_8; rw_tend = 0.0_8; rt_tend = 0.0_8; moist_tend = 0.0_8

  CALL sixth_order_diffusion('u', u, ru_tend, mut, dt, config_flags, c1_h, c2_h, &
                             2, 0.12_8, phb, ph, 1.0_8/dx, 1.0_8/dy, &
                             msftx, msfty, msfux, msfuy, msfvx, msfvy, &
                             ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL sixth_order_diffusion('v', v, rv_tend, mut, dt, config_flags, c1_h, c2_h, &
                             2, 0.12_8, phb, ph, 1.0_8/dx, 1.0_8/dy, &
                             msftx, msfty, msfux, msfuy, msfvx, msfvy, &
                             ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL sixth_order_diffusion('w', w, rw_tend, mut, dt, config_flags, c1f_h, c2f_h, &
                             2, 0.12_8, phb, ph, 1.0_8/dx, 1.0_8/dy, &
                             msftx, msfty, msfux, msfuy, msfvx, msfvy, &
                             ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL sixth_order_diffusion('m', theta, rt_tend, mut, dt, config_flags, c1_h, c2_h, &
                             2, 0.12_8, phb, ph, 1.0_8/dx, 1.0_8/dy, &
                             msftx, msfty, msfux, msfuy, msfvx, msfvy, &
                             ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
  CALL sixth_order_diffusion('m', moist(:,:,:,1), moist_tend(:,:,:,1), mut, dt, config_flags, c1_h, c2_h, &
                             2, 0.12_8, phb, ph, 1.0_8/dx, 1.0_8/dy, &
                             msftx, msfty, msfux, msfuy, msfvx, msfvy, &
                             ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

  CALL write_chunk(out_unit, "u_6th_opt2", ru_tend(1:nx+1, 1:nz, 1:ny), (nx+1)*nz*ny)
  CALL write_chunk(out_unit, "v_6th_opt2", rv_tend(1:nx, 1:nz, 1:ny+1), nx*nz*(ny+1))
  CALL write_chunk(out_unit, "w_6th_opt2", rw_tend(1:nx, 1:nz+1, 1:ny), nx*(nz+1)*ny)
  CALL write_chunk(out_unit, "t_6th_opt2", rt_tend(1:nx, 1:nz, 1:ny), nx*nz*ny)
  CALL write_chunk(out_unit, "m_6th_opt2", moist_tend(1:nx, 1:nz, 1:ny, 1), nx*nz*ny)

  ! ============================================================================
  ! TEST 7: w_damp
  ! ============================================================================
  config_flags%w_damping = 1
  config_flags%w_crit_cfl = 1.0_8
  config_flags%zadvect_implicit = 0
  config_flags%polar = .FALSE.
  rw_tend = 0.0_8

  CALL w_damp(rw_tend, max_vert_cfl, max_horiz_cfl, &
              u, v, ww, w, mut, c1f_h, c2f_h, rdnw, 1.0_8/dx, 1.0_8/dy, &
              msfux, msfuy, msfvx, msfvy, dt, config_flags, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

  CALL write_chunk(out_unit, "rw_tend_damp", rw_tend(1:nx, 1:nz+1, 1:ny), nx*(nz+1)*ny)
  CALL write_chunk(out_unit, "max_vert_cfl", (/ max_vert_cfl /), 1)
  CALL write_chunk(out_unit, "max_horiz_cfl", (/ max_horiz_cfl /), 1)

  CLOSE(out_unit)
  PRINT *, "Fortran driver finished successfully."

CONTAINS

  SUBROUTINE write_chunk(unit, name, arr, n)
    INTEGER, INTENT(IN) :: unit, n
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(8), INTENT(IN) :: arr(n)
    CHARACTER(LEN=32) :: tag

    tag = TRIM(name)
    WRITE(unit) tag
    WRITE(unit) n
    WRITE(unit) arr
  END SUBROUTINE write_chunk

  SUBROUTINE fill_periodic_halos()
    ! Map factors
    msftx(0, 1:ny) = msftx(nx, 1:ny); msftx(nx+1, 1:ny) = msftx(1, 1:ny)
    msftx(:, 0) = msftx(:, ny); msftx(:, ny+1) = msftx(:, 1)

    msfty(0, 1:ny) = msfty(nx, 1:ny); msfty(nx+1, 1:ny) = msfty(1, 1:ny)
    msfty(:, 0) = msfty(:, ny); msfty(:, ny+1) = msfty(:, 1)

    msfux(0, 1:ny) = msfux(nx, 1:ny); msfux(nx+1, 1:ny) = msfux(1, 1:ny); msfux(nx+2, 1:ny) = msfux(2, 1:ny)
    msfux(:, 0) = msfux(:, ny); msfux(:, ny+1) = msfux(:, 1)

    msfuy(0, 1:ny) = msfuy(nx, 1:ny); msfuy(nx+1, 1:ny) = msfuy(1, 1:ny); msfuy(nx+2, 1:ny) = msfuy(2, 1:ny)
    msfuy(:, 0) = msfuy(:, ny); msfuy(:, ny+1) = msfuy(:, 1)

    msfvx(0, 1:ny+1) = msfvx(nx, 1:ny+1); msfvx(nx+1, 1:ny+1) = msfvx(1, 1:ny+1)
    msfvx(:, 0) = msfvx(:, ny); msfvx(:, ny+1) = msfvx(:, 1); msfvx(:, ny+2) = msfvx(:, 2)

    msfvy(0, 1:ny+1) = msfvy(nx, 1:ny+1); msfvy(nx+1, 1:ny+1) = msfvy(1, 1:ny+1)
    msfvy(:, 0) = msfvy(:, ny); msfvy(:, ny+1) = msfvy(:, 1); msfvy(:, ny+2) = msfvy(:, 2)

    ! 3D fields
    ph(0, 1:nz+1, 1:ny) = ph(nx, 1:nz+1, 1:ny); ph(nx+1, 1:nz+1, 1:ny) = ph(1, 1:nz+1, 1:ny)
    ph(:, :, 0) = ph(:, :, ny); ph(:, :, ny+1) = ph(:, :, 1)

    phb(0, 1:nz+1, 1:ny) = phb(nx, 1:nz+1, 1:ny); phb(nx+1, 1:nz+1, 1:ny) = phb(1, 1:nz+1, 1:ny)
    phb(:, :, 0) = phb(:, :, ny); phb(:, :, ny+1) = phb(:, :, 1)

    u(0, 1:nz, 1:ny) = u(nx, 1:nz, 1:ny); u(nx+2, 1:nz, 1:ny) = u(2, 1:nz, 1:ny)
    u(:, :, 0) = u(:, :, ny); u(:, :, ny+1) = u(:, :, 1)

    v(0, 1:nz, 1:ny+1) = v(nx, 1:nz, 1:ny+1); v(nx+1, 1:nz, 1:ny+1) = v(1, 1:nz, 1:ny+1)
    v(:, :, 0) = v(:, :, ny); v(:, :, ny+2) = v(:, :, 2)

    w(0, 1:nz+1, 1:ny) = w(nx, 1:nz+1, 1:ny); w(nx+1, 1:nz+1, 1:ny) = w(1, 1:nz+1, 1:ny)
    w(:, :, 0) = w(:, :, ny); w(:, :, ny+1) = w(:, :, 1)

    theta(0, 1:nz, 1:ny) = theta(nx, 1:nz, 1:ny); theta(nx+1, 1:nz, 1:ny) = theta(1, 1:nz, 1:ny)
    theta(:, :, 0) = theta(:, :, ny); theta(:, :, ny+1) = theta(:, :, 1)

    rdz(0, 1:nz, 1:ny) = rdz(nx, 1:nz, 1:ny); rdz(nx+1, 1:nz, 1:ny) = rdz(1, 1:nz, 1:ny)
    rdz(:, :, 0) = rdz(:, :, ny); rdz(:, :, ny+1) = rdz(:, :, 1)

    rdzw(0, 1:nz, 1:ny) = rdzw(nx, 1:nz, 1:ny); rdzw(nx+1, 1:nz, 1:ny) = rdzw(1, 1:nz, 1:ny)
    rdzw(:, :, 0) = rdzw(:, :, ny); rdzw(:, :, ny+1) = rdzw(:, :, 1)

    zx(0, 1:nz+1, 1:ny) = zx(nx, 1:nz+1, 1:ny); zx(nx+2, 1:nz+1, 1:ny) = zx(2, 1:nz+1, 1:ny)
    zx(:, :, 0) = zx(:, :, ny); zx(:, :, ny+1) = zx(:, :, 1)

    zy(0, 1:nz+1, 1:ny+1) = zy(nx, 1:nz+1, 1:ny+1); zy(nx+1, 1:nz+1, 1:ny+1) = zy(1, 1:nz+1, 1:ny+1)
    zy(:, :, 0) = zy(:, :, ny); zy(:, :, ny+2) = zy(:, :, 2)
  END SUBROUTINE fill_periodic_halos

END PROGRAM driver_diffusion
