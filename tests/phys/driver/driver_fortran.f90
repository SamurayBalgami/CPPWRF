program driver_fortran
  use module_wrf_phys_driver_orig
  implicit none

  character(len=256) :: infile, outfile
  integer :: iunit, ounit
  integer :: nx, ny, nz, n_moist, use_theta_m
  real(8) :: p_top

  type(grid_config_rec_type) :: config_flags

  ! 1D vertical arrays
  real(8), allocatable :: znw(:), znu(:), dnw(:), dn(:), fnm(:), fnp(:)
  real(8), allocatable :: c1h(:), c2h(:), c1f(:), c2f(:)

  ! 2D arrays
  real(8), allocatable :: mut(:,:), muu(:,:), muv(:,:)
  real(8), allocatable :: msftx(:,:), msfty(:,:)
  real(8), allocatable :: msfux(:,:), msfuy(:,:)
  real(8), allocatable :: msfvx(:,:), msfvy(:,:), msfvx_inv(:,:)
  real(8), allocatable :: ht(:,:)

  ! 3D atmospheric arrays (ims:ime, kms:kme, jms:jme)
  real(8), allocatable :: u(:,:,:), v(:,:,:), t(:,:,:), p(:,:,:), pb(:,:,:)
  real(8), allocatable :: alt(:,:,:), ph(:,:,:), phb(:,:,:), moist(:,:,:,:)

  ! 3D uncoupled physics tendencies
  real(8), allocatable :: rthraten(:,:,:), rthblten(:,:,:), rublten(:,:,:), rvblten(:,:,:)
  real(8), allocatable :: rqvblten(:,:,:), rqcblten(:,:,:), rqiblten(:,:,:)
  real(8), allocatable :: rthcuten(:,:,:), rucuten(:,:,:), rvcuten(:,:,:)
  real(8), allocatable :: rqvcuten(:,:,:), rqccuten(:,:,:), rqrcuten(:,:,:), rqicuten(:,:,:), rqscuten(:,:,:)
  real(8), allocatable :: rthshten(:,:,:), rushten(:,:,:), rvshten(:,:,:)
  real(8), allocatable :: rqvshten(:,:,:), rqcshten(:,:,:), rqrshten(:,:,:), rqishten(:,:,:), rqsshten(:,:,:), rqgshten(:,:,:)

  ! Output arrays from phy_prep
  real(8), allocatable :: th_phy(:,:,:), th_phy_m_t0(:,:,:), p_phy(:,:,:), pi_phy(:,:,:)
  real(8), allocatable :: t_phy(:,:,:), rho(:,:,:), u_phy(:,:,:), v_phy(:,:,:)
  real(8), allocatable :: z(:,:,:), dz8w(:,:,:), z_at_w(:,:,:), p8w(:,:,:), t8w(:,:,:)
  real(8), allocatable :: p_hyd(:,:,:), p_hyd_w(:,:,:)

  ! Tendencies outputs
  real(8), allocatable :: ru_tendf(:,:,:), rv_tendf(:,:,:), t_tendf(:,:,:)
  real(8), allocatable :: moist_tendf(:,:,:,:)
  real(8), allocatable :: ru_tend(:,:,:), rv_tend(:,:,:), t_tend(:,:,:)

  ! Dummy arrays for verbatim original WRF signatures
  real(8), allocatable :: dummy_3d(:,:,:), dummy_2d(:,:)
  real(8), allocatable :: dummy_scalar(:,:,:,:), dummy_tracer(:,:,:,:)
  real(8), allocatable :: rw_tend(:,:,:), rw_tendf(:,:,:), ph_tend(:,:,:), ph_tendf(:,:,:)
  real(8), allocatable :: u_save(:,:,:), v_save(:,:,:), w_save(:,:,:), ph_save(:,:,:), t_save(:,:,:)
  real(8), allocatable :: mu_tend(:,:), mu_tendf(:,:)
  real(8), allocatable :: h_diabatic(:,:,:)

  integer :: i, j, k, m

  if (command_argument_count() >= 2) then
    call get_command_argument(1, infile)
    call get_command_argument(2, outfile)
  else
    infile = "driver_input.bin"
    outfile = "driver_output_fortran.bin"
  end if

  open(newunit=iunit, file=trim(infile), status='old', access='stream', form='unformatted')

  read(iunit) nx, ny, nz, n_moist
  read(iunit) p_top
  read(iunit) use_theta_m

  config_flags%use_theta_m = use_theta_m
  config_flags%ra_lw_physics = 1
  config_flags%ra_sw_physics = 1
  config_flags%bl_pbl_physics = 1
  config_flags%cu_physics = 11
  config_flags%shcu_physics = 2
  config_flags%periodic_x = .TRUE.
  config_flags%specified = .FALSE.
  config_flags%nested = .FALSE.

  ! Allocate 1D arrays
  allocate(znw(nz+1), znu(nz), dnw(nz), dn(nz+1), fnm(nz+1), fnp(nz+1))
  allocate(c1h(nz), c2h(nz), c1f(nz+1), c2f(nz+1))

  read(iunit) znw
  read(iunit) znu
  read(iunit) dnw
  read(iunit) dn
  read(iunit) fnm
  read(iunit) fnp
  read(iunit) c1h
  read(iunit) c2h
  read(iunit) c1f
  read(iunit) c2f

  ! Allocate 2D arrays
  allocate(mut(0:nx+1, 0:ny+1), muu(0:nx+1, 0:ny+1), muv(0:nx+1, 0:ny+1))
  allocate(msftx(0:nx+1, 0:ny+1), msfty(0:nx+1, 0:ny+1))
  allocate(msfux(0:nx+1, 0:ny+1), msfuy(0:nx+1, 0:ny+1))
  allocate(msfvx(0:nx+1, 0:ny+1), msfvy(0:nx+1, 0:ny+1), msfvx_inv(0:nx+1, 0:ny+1))
  allocate(ht(0:nx+1, 0:ny+1))

  mut = 0.0d0; muu = 0.0d0; muv = 0.0d0
  msftx = 1.0d0; msfty = 1.0d0; msfux = 1.0d0; msfuy = 1.0d0; msfvx = 1.0d0; msfvy = 1.0d0; msfvx_inv = 1.0d0
  ht = 0.0d0

  do j = 1, ny
    read(iunit) (mut(i, j), i = 1, nx)
  end do
  do j = 1, ny
    read(iunit) (muu(i, j), i = 1, nx + 1)
  end do
  do j = 1, ny + 1
    read(iunit) (muv(i, j), i = 1, nx)
  end do
  do j = 1, ny
    read(iunit) (msftx(i, j), i = 1, nx)
  end do
  do j = 1, ny
    read(iunit) (msfty(i, j), i = 1, nx)
  end do
  do j = 1, ny
    read(iunit) (msfux(i, j), i = 1, nx + 1)
  end do
  do j = 1, ny
    read(iunit) (msfuy(i, j), i = 1, nx + 1)
  end do
  do j = 1, ny + 1
    read(iunit) (msfvx(i, j), i = 1, nx)
  end do
  do j = 1, ny + 1
    read(iunit) (msfvy(i, j), i = 1, nx)
  end do
  do j = 1, ny + 1
    read(iunit) (msfvx_inv(i, j), i = 1, nx)
  end do
  do j = 1, ny
    read(iunit) (ht(i, j), i = 1, nx)
  end do

  ! Allocate 3D arrays: (0:nx+1, 1:nz+1, 0:ny+1)
  allocate(u(0:nx+1, 1:nz+1, 0:ny+1), v(0:nx+1, 1:nz+1, 0:ny+1), t(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(p(0:nx+1, 1:nz+1, 0:ny+1), pb(0:nx+1, 1:nz+1, 0:ny+1), alt(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(ph(0:nx+1, 1:nz+1, 0:ny+1), phb(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(moist(0:nx+1, 1:nz+1, 0:ny+1, n_moist))

  u = 0.0d0; v = 0.0d0; t = 0.0d0; p = 0.0d0; pb = 0.0d0; alt = 0.0d0; ph = 0.0d0; phb = 0.0d0; moist = 0.0d0

  ! Read u (nz, ny, nx+1)
  do k = 1, nz
    do j = 1, ny
      read(iunit) (u(i, k, j), i = 1, nx + 1)
    end do
  end do

  ! Read v (nz, ny+1, nx)
  do k = 1, nz
    do j = 1, ny + 1
      read(iunit) (v(i, k, j), i = 1, nx)
    end do
  end do

  ! Read t (nz, ny, nx)
  do k = 1, nz
    do j = 1, ny
      read(iunit) (t(i, k, j), i = 1, nx)
    end do
  end do

  ! Read p (nz, ny, nx)
  do k = 1, nz
    do j = 1, ny
      read(iunit) (p(i, k, j), i = 1, nx)
    end do
  end do

  ! Read pb (nz, ny, nx)
  do k = 1, nz
    do j = 1, ny
      read(iunit) (pb(i, k, j), i = 1, nx)
    end do
  end do

  ! Read alt (nz, ny, nx)
  do k = 1, nz
    do j = 1, ny
      read(iunit) (alt(i, k, j), i = 1, nx)
    end do
  end do

  ! Read ph (nz+1, ny, nx)
  do k = 1, nz + 1
    do j = 1, ny
      read(iunit) (ph(i, k, j), i = 1, nx)
    end do
  end do

  ! Read phb (nz+1, ny, nx)
  do k = 1, nz + 1
    do j = 1, ny
      read(iunit) (phb(i, k, j), i = 1, nx)
    end do
  end do

  ! Read moist (n_moist species, each nz, ny, nx)
  do m = 1, n_moist
    do k = 1, nz
      do j = 1, ny
        read(iunit) (moist(i, k, j, m), i = 1, nx)
      end do
    end do
  end do

  ! Allocate and read 3D tendencies
  allocate(rthraten(0:nx+1, 1:nz+1, 0:ny+1), rthblten(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(rublten(0:nx+1, 1:nz+1, 0:ny+1), rvblten(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(rqvblten(0:nx+1, 1:nz+1, 0:ny+1), rqcblten(0:nx+1, 1:nz+1, 0:ny+1), rqiblten(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(rthcuten(0:nx+1, 1:nz+1, 0:ny+1), rucuten(0:nx+1, 1:nz+1, 0:ny+1), rvcuten(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(rqvcuten(0:nx+1, 1:nz+1, 0:ny+1), rqccuten(0:nx+1, 1:nz+1, 0:ny+1), rqrcuten(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(rqicuten(0:nx+1, 1:nz+1, 0:ny+1), rqscuten(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(rthshten(0:nx+1, 1:nz+1, 0:ny+1), rushten(0:nx+1, 1:nz+1, 0:ny+1), rvshten(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(rqvshten(0:nx+1, 1:nz+1, 0:ny+1), rqcshten(0:nx+1, 1:nz+1, 0:ny+1), rqrshten(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(rqishten(0:nx+1, 1:nz+1, 0:ny+1), rqsshten(0:nx+1, 1:nz+1, 0:ny+1), rqgshten(0:nx+1, 1:nz+1, 0:ny+1))

  rthraten=0.0d0; rthblten=0.0d0; rublten=0.0d0; rvblten=0.0d0; rqvblten=0.0d0; rqcblten=0.0d0; rqiblten=0.0d0
  rthcuten=0.0d0; rucuten=0.0d0; rvcuten=0.0d0; rqvcuten=0.0d0; rqccuten=0.0d0; rqrcuten=0.0d0; rqicuten=0.0d0; rqscuten=0.0d0
  rthshten=0.0d0; rushten=0.0d0; rvshten=0.0d0; rqvshten=0.0d0; rqcshten=0.0d0; rqrshten=0.0d0; rqishten=0.0d0; rqsshten=0.0d0; rqgshten=0.0d0

  do k = 1, nz; do j = 1, ny; read(iunit) (rthraten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rthblten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rublten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rvblten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rqvblten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rqcblten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rqiblten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rthcuten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rucuten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rvcuten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rqvcuten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rqccuten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rqrcuten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rqicuten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rqscuten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rthshten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rushten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rvshten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rqvshten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rqcshten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rqrshten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rqishten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rqsshten(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; read(iunit) (rqgshten(i, k, j), i = 1, nx); end do; end do

  close(iunit)

  ! Periodic boundary condition on staggered faces
  do k = 1, nz
    do j = 1, ny
      u(nx + 1, k, j) = u(1, k, j)
    end do
  end do
  do k = 1, nz
    do i = 1, nx
      v(i, k, ny + 1) = v(i, k, 1)
    end do
  end do

  ! Allocate diagnostic outputs
  ! Allocate diagnostic outputs
  allocate(th_phy(0:nx+1, 1:nz+1, 0:ny+1), th_phy_m_t0(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(p_phy(0:nx+1, 1:nz+1, 0:ny+1), pi_phy(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(t_phy(0:nx+1, 1:nz+1, 0:ny+1), rho(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(u_phy(0:nx+1, 1:nz+1, 0:ny+1), v_phy(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(z(0:nx+1, 1:nz+1, 0:ny+1), dz8w(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(z_at_w(0:nx+1, 1:nz+1, 0:ny+1), p8w(0:nx+1, 1:nz+1, 0:ny+1), t8w(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(p_hyd(0:nx+1, 1:nz+1, 0:ny+1), p_hyd_w(0:nx+1, 1:nz+1, 0:ny+1))

  th_phy = 0.0d0; th_phy_m_t0 = 0.0d0; p_phy = 0.0d0; pi_phy = 0.0d0; t_phy = 0.0d0; rho = 0.0d0
  u_phy = 0.0d0; v_phy = 0.0d0; z = 0.0d0; dz8w = 0.0d0; z_at_w = 0.0d0; p8w = 0.0d0; t8w = 0.0d0
  p_hyd = 0.0d0; p_hyd_w = 0.0d0

  ! Allocate dummy arrays for original WRF argument lists
  allocate(dummy_3d(0:nx+1, 1:nz+1, 0:ny+1), dummy_2d(0:nx+1, 0:ny+1))
  allocate(dummy_scalar(0:nx+1, 1:nz+1, 0:ny+1, 1), dummy_tracer(0:nx+1, 1:nz+1, 0:ny+1, 1))
  allocate(rw_tend(0:nx+1, 1:nz+1, 0:ny+1), rw_tendf(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(ph_tend(0:nx+1, 1:nz+1, 0:ny+1), ph_tendf(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(u_save(0:nx+1, 1:nz+1, 0:ny+1), v_save(0:nx+1, 1:nz+1, 0:ny+1), w_save(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(ph_save(0:nx+1, 1:nz+1, 0:ny+1), t_save(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(mu_tend(0:nx+1, 0:ny+1), mu_tendf(0:nx+1, 0:ny+1))
  allocate(h_diabatic(0:nx+1, 1:nz+1, 0:ny+1))

  dummy_3d = 0.0d0; dummy_2d = 0.0d0; dummy_scalar = 0.0d0; dummy_tracer = 0.0d0
  rw_tend = 0.0d0; rw_tendf = 0.0d0; ph_tend = 0.0d0; ph_tendf = 0.0d0
  u_save = 0.0d0; v_save = 0.0d0; w_save = 0.0d0; ph_save = 0.0d0; t_save = 0.0d0
  mu_tend = 0.0d0; mu_tendf = 0.0d0; h_diabatic = 0.0d0

  ! 1. Call phy_prep with kte = nz + 1 (kde) as in WRF solve_em / module_first_rk_step_part1.F
  call phy_prep(config_flags,                                    &
                mut, muu, muv,                                   &
                c1h, c2h, c1f, c2f,                              &
                u, v, p, pb, alt, ph,                            &
                phb, t, moist, n_moist,                          &
                rho, th_phy, th_phy_m_t0,                        &
                p_phy, pi_phy,                                   &
                u_phy, v_phy, p8w, t_phy, t8w,                   &
                z, z_at_w, dz8w,                                 &
                p_hyd, p_hyd_w, dnw,                             &
                fnm, fnp, znw, p_top,                            &
                1, nx+1, 1, ny+1, 1, nz+1,                       &
                0, nx+1, 0, ny+1, 1, nz+1,                       &
                1, nx, 1, ny, 1, nz+1                            )

  ! 2. Call verbatim calculate_phy_tend (mu-coupling) with kte = nz + 1
  call calculate_phy_tend(config_flags, c1h, c2h,                &
                          mut, muu, muv, pi_phy,                 &
                          rthraten,                              &
                          rublten, rvblten, rthblten,            &
                          rqvblten, rqcblten, rqiblten,          &
                          rucuten, rvcuten, rthcuten,            &
                          rqvcuten, rqccuten, rqrcuten,          &
                          rqicuten, rqscuten,                    &
                          rushten, rvshten, rthshten,            &
                          rqvshten, rqcshten, rqrshten,          &
                          rqishten, rqsshten, rqgshten,          &
                          dummy_3d, dummy_3d, dummy_3d, dummy_3d,&
                          dummy_2d,                              &
                          dummy_scalar, dummy_scalar, 0,         &
                          dummy_tracer, dummy_tracer, 0,         &
                          1, nx+1, 1, ny+1, 1, nz+1,             &
                          0, nx+1, 0, ny+1, 1, nz+1,             &
                          1, nx, 1, ny, 1, nz+1                  )

  ! Periodic halo fill for tendencies before A-to-C stagger interpolation
  do k = 1, nz
    do j = 1, ny
      rublten(0, k, j) = rublten(nx, k, j)
      rublten(nx+1, k, j) = rublten(1, k, j)
      rucuten(0, k, j) = rucuten(nx, k, j)
      rucuten(nx+1, k, j) = rucuten(1, k, j)
      rushten(0, k, j) = rushten(nx, k, j)
      rushten(nx+1, k, j) = rushten(1, k, j)
    end do
  end do

  do k = 1, nz
    do i = 1, nx
      rvblten(i, k, 0) = rvblten(i, k, ny)
      rvblten(i, k, ny+1) = rvblten(i, k, 1)
      rvcuten(i, k, 0) = rvcuten(i, k, ny)
      rvcuten(i, k, ny+1) = rvcuten(i, k, 1)
      rvshten(i, k, 0) = rvshten(i, k, ny)
      rvshten(i, k, ny+1) = rvshten(i, k, 1)
    end do
  end do

  ! 3. Call verbatim update_phy_ten with kte = nz + 1
  allocate(ru_tendf(0:nx+1, 1:nz+1, 0:ny+1), rv_tendf(0:nx+1, 1:nz+1, 0:ny+1), t_tendf(0:nx+1, 1:nz+1, 0:ny+1))
  allocate(moist_tendf(0:nx+1, 1:nz+1, 0:ny+1, n_moist))
  ru_tendf = 0.0d0; rv_tendf = 0.0d0; t_tendf = 0.0d0; moist_tendf = 0.0d0

  call update_phy_ten(dummy_3d, t_tendf, ru_tendf, rv_tendf, moist_tendf, &
                      dummy_scalar, dummy_2d,                             &
                      rthraten, rthblten, rthcuten, rthshten,             &
                      rublten, rucuten, rushten,                          &
                      rvblten, rvcuten, rvshten,                          &
                      rqvblten, rqcblten, rqiblten, dummy_3d,             &
                      rqvcuten, rqccuten, rqrcuten, rqicuten, rqscuten,   &
                      dummy_3d, dummy_3d,                                 &
                      rqvshten, rqcshten, rqrshten, rqishten, rqsshten, rqgshten, &
                      dummy_3d, dummy_3d,                                 &
                      dummy_3d, dummy_3d, dummy_3d, dummy_3d,             &
                      dummy_3d, dummy_2d,                                 &
                      dummy_3d, dummy_3d,                                 &
                      n_moist, 0, config_flags, 1, .FALSE.,               &
                      1, nx+1, 1, ny+1, 1, nz+1,                          &
                      0, nx+1, 0, ny+1, 1, nz+1,                          &
                      1, nx, 1, ny, 1, nz+1                               )

  ! Periodic copy for staggered tendencies at domain edge
  do k = 1, nz
    do j = 1, ny
      ru_tendf(nx + 1, k, j) = ru_tendf(1, k, j)
    end do
  end do
  do k = 1, nz
    do i = 1, nx
      rv_tendf(i, k, ny + 1) = rv_tendf(i, k, 1)
    end do
  end do

  ! 4. Call verbatim rk_addtend_dry (map factor divisions) with kte = nz + 1
  allocate(ru_tend(0:nx+1, 1:nz+1, 0:ny+1), rv_tend(0:nx+1, 1:nz+1, 0:ny+1), t_tend(0:nx+1, 1:nz+1, 0:ny+1))
  ru_tend = 0.0d0; rv_tend = 0.0d0; t_tend = 0.0d0

  call rk_addtend_dry ( ru_tend, rv_tend, rw_tend, ph_tend, t_tend,      &
                        ru_tendf, rv_tendf, rw_tendf, ph_tendf, t_tendf, &
                        u_save, v_save, w_save, ph_save, t_save,         &
                        mu_tend, mu_tendf, 2, c1h, c2h,                  &
                        h_diabatic, mut, msftx, msfty, msfux, msfuy,     &
                        msfvx, msfvx_inv, msfvy,                         &
                        1, nx+1, 1, ny+1, 1, nz+1,                       &
                        0, nx+1, 0, ny+1, 1, nz+1,                       &
                        1, nx+1, 1, ny+1, 1, nz+1,                       &
                        1, nx, 1, ny, 1, nz+1                            )

  ! Periodic copy for large-step tendencies
  do k = 1, nz
    do j = 1, ny
      ru_tend(nx + 1, k, j) = ru_tend(1, k, j)
    end do
  end do
  do k = 1, nz
    do i = 1, nx
      rv_tend(i, k, ny + 1) = rv_tend(i, k, 1)
    end do
  end do

  ! Write output binary
  open(newunit=ounit, file=trim(outfile), status='replace', access='stream', form='unformatted')

  ! phy_prep outputs
  do k = 1, nz; do j = 1, ny; write(ounit) (th_phy(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; write(ounit) (th_phy_m_t0(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; write(ounit) (p_phy(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; write(ounit) (pi_phy(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; write(ounit) (t_phy(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; write(ounit) (rho(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; write(ounit) (u_phy(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; write(ounit) (v_phy(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; write(ounit) (z(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz + 1; do j = 1, ny; write(ounit) (dz8w(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz + 1; do j = 1, ny; write(ounit) (z_at_w(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz + 1; do j = 1, ny; write(ounit) (p8w(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz + 1; do j = 1, ny; write(ounit) (t8w(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; write(ounit) (p_hyd(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz + 1; do j = 1, ny; write(ounit) (p_hyd_w(i, k, j), i = 1, nx); end do; end do

  ! Tendency coupling outputs
  do k = 1, nz; do j = 1, ny; write(ounit) (ru_tendf(i, k, j), i = 1, nx + 1); end do; end do
  do k = 1, nz; do j = 1, ny + 1; write(ounit) (rv_tendf(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; write(ounit) (t_tendf(i, k, j), i = 1, nx); end do; end do
  do m = 1, n_moist
    do k = 1, nz; do j = 1, ny; write(ounit) (moist_tendf(i, k, j, m), i = 1, nx); end do; end do
  end do

  do k = 1, nz; do j = 1, ny; write(ounit) (ru_tend(i, k, j), i = 1, nx + 1); end do; end do
  do k = 1, nz; do j = 1, ny + 1; write(ounit) (rv_tend(i, k, j), i = 1, nx); end do; end do
  do k = 1, nz; do j = 1, ny; write(ounit) (t_tend(i, k, j), i = 1, nx); end do; end do

  close(ounit)

end program driver_fortran
