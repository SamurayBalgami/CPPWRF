! Reference driver for the P1 physics-suite differential test (tests/phys/suite).
!
! Links against an UNMODIFIED WRF v4.6.1 build configured with "-r8" (serial gfortran,
! -fdefault-real-8) and drives WRF's own code: the framework sets up the domain from
! namelist.input (wrf_init without the wrfinput read), the synthetic tile of suite_input.bin
! is copied into grid% in place of med_initialdata_input, WRF's start_domain runs (phy_init,
! table reads), and then, three times, the statements solve_em executes for the column
! physics of RK stage 1 and after the RK loop:
!     grid%itimestep = grid%itimestep + 1
!     first_rk_step_part1 (init_zero_tendency, phy_prep, radiation/surface/pbl/cumulus drivers)
!     first_rk_step_part2 (calculate_phy_tend, update_phy_ten, conv_t_tendf_to_moist, ...)
!     phy_prep_part2      (decoupling of the R*TEN fields with muts)
! No WRF statement is changed and no WRF module is stubbed: only the dynamics between the
! calls is replaced by a deterministic change of the atmospheric state, dumped as input.
!
! Output: suite_ref.bin = for each call n: records "IN<n>:<name>" (state before the call) and
! "OUT<n>:<name>" (after phy_prep_part2), same record format as suite_input.bin.
PROGRAM suite_ref
   USE module_domain
   USE module_configure
   USE module_state_description
   USE module_model_constants
   USE module_driver_constants
   USE module_timing
   USE module_wrf_error
   USE module_tiles, ONLY : set_tiles
   USE module_first_rk_step_part1, ONLY : first_rk_step_part1
   USE module_first_rk_step_part2, ONLY : first_rk_step_part2
   USE module_big_step_utilities_em, ONLY : phy_prep_part2
   USE module_wrf_top
   USE module_utility
   IMPLICIT NONE

   INTEGER :: ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, ips, ipe, jps, jpe, kps, kpe
   INTEGER :: imsx, imex, jmsx, jmex, kmsx, kmex, ipsx, ipex, jpsx, jpex, kpsx, kpex
   INTEGER :: imsy, imey, jmsy, jmey, kmsy, kmey, ipsy, ipey, jpsy, jpey, kpsy, kpey
   INTEGER :: k_start, k_end, i, j, k, n, ncall, ios, nd, dims(3), unit_out
   REAL, ALLOCATABLE :: buf(:)
   CHARACTER(LEN=32) :: name
   CHARACTER(LEN=256) :: mminlu_in
   TYPE(domain), POINTER :: grid
   REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: moist_tend, moist_old, chem_tend, tracer_tend, scalar_tend
   REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ru_tendf, rv_tendf, rw_tendf, t_tendf, ph_tendf, tke_tend
   REAL, ALLOCATABLE, DIMENSION(:,:,:) :: th_phy, pi_phy, p_phy, dz8w, p8w, t8w
   REAL, ALLOCATABLE, DIMENSION(:,:) :: mu_tendf, psim, psih, gz1oz0, chklowq, hol
   LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: cu_act_flag
   REAL :: curr_secs
   LOGICAL :: adapt_step_flag, f_flux, restart_flag

   ! ---------------- wrf_init without med_initialdata_input ----------------
   CALL init_modules(1)
   CALL WRFU_Initialize( defaultCalKind=WRFU_CAL_GREGORIAN )
   CALL init_modules(2)
   CALL initial_config
   CALL setup_physics_suite
   CALL set_derived_rconfigs
   CALL check_nml_consistency
   CALL set_physics_rconfigs
   NULLIFY( null_domain )
   CALL alloc_and_configure_domain ( domain_id = 1, active_this_task = domain_active_this_task(1), &
                                     grid = head_grid, parent = null_domain, kid = -1 )
   CALL model_to_grid_config_rec ( head_grid%id , model_config_rec , config_flags )
   CALL set_scalar_indices_from_config ( head_grid%id , ios, nd )
   CALL init_wrfio
   CALL Setup_Timekeeping (head_grid)
   grid => head_grid

   CALL get_ijk_from_grid ( grid, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                            ips, ipe, jps, jpe, kps, kpe, imsx, imex, jmsx, jmex, kmsx, kmex, &
                            ipsx, ipex, jpsx, jpex, kpsx, kpex, imsy, imey, jmsy, jmey, kmsy, kmey, &
                            ipsy, ipey, jpsy, jpey, kpsy, kpey )
   k_start = kps
   k_end   = kpe

   ! ---------------- synthetic wrfinput (what input_wrf would fill) ----------------
   ! global attributes of a MODIS (20 category) real-data input
   mminlu_in = 'MODIFIED_IGBP_MODIS_NOAH'
   CALL nl_set_mminlu ( 1, mminlu_in )
   CALL nl_set_iswater( 1, 17 )
   CALL nl_set_isice  ( 1, 15 )
   CALL nl_set_isurban( 1, 13 )
   CALL nl_set_isoilwater( 1, 14 )
   CALL nl_set_islake ( 1, -1 )
   CALL nl_set_map_proj ( 1, 1 )
   CALL nl_set_cen_lat ( 1, 39.3 )
   CALL nl_set_truelat1 ( 1, 35. )
   CALL nl_set_truelat2 ( 1, 42. )
   CALL nl_set_stand_lon ( 1, 34. )
   CALL model_to_grid_config_rec ( head_grid%id , model_config_rec , config_flags )
   grid%this_is_an_ideal_run = .FALSE.
   grid%p_top = 5000.
   grid%msftx = 1. ; grid%msfty = 1. ; grid%msfux = 1. ; grid%msfuy = 1.
   grid%msfvx = 1. ; grid%msfvy = 1. ; grid%msfvx_inv = 1.
   grid%msft = 1. ; grid%msfu = 1. ; grid%msfv = 1.
   grid%sina = 0. ; grid%cosa = 1. ; grid%f = 1.E-4 ; grid%e = 1.E-4
   grid%lakemask = 0. ; grid%imask_nostag = 1
   grid%c1h = 1. ; grid%c2h = 0. ; grid%c1f = 1. ; grid%c2f = 0. ; grid%c4h = 0. ; grid%c4f = 0.

   OPEN(10, FILE='suite_input.bin', ACCESS='stream', FORM='unformatted', STATUS='old', &
        CONVERT='little_endian')
   DO
      READ(10, IOSTAT=ios) name
      IF ( ios /= 0 ) EXIT
      READ(10) nd, dims
      ALLOCATE(buf(dims(1)*dims(2)*dims(3)))
      READ(10) buf
      CALL put_input(TRIM(name), dims, buf)
      DEALLOCATE(buf)
   END DO
   CLOSE(10)
   DO k = kds, kde
      grid%c3f(k) = grid%znw(k)
   END DO
   DO k = kds, kde-1
      grid%c3h(k) = grid%znu(k)
   END DO
   grid%rdx = 1./grid%dx
   grid%rdy = 1./grid%dy

   CALL start_domain ( grid , .TRUE. )

   ! ---------------- arrays that are locals of solve_em ----------------
   ALLOCATE( moist_tend(ims:ime,kms:kme,jms:jme,num_moist), moist_old(ims:ime,kms:kme,jms:jme,num_moist) )
   ALLOCATE( chem_tend(ims:ime,kms:kme,jms:jme,num_chem), tracer_tend(ims:ime,kms:kme,jms:jme,num_tracer) )
   ALLOCATE( scalar_tend(ims:ime,kms:kme,jms:jme,num_scalar) )
   ALLOCATE( ru_tendf(ims:ime,kms:kme,jms:jme), rv_tendf(ims:ime,kms:kme,jms:jme), rw_tendf(ims:ime,kms:kme,jms:jme) )
   ALLOCATE( t_tendf(ims:ime,kms:kme,jms:jme), ph_tendf(ims:ime,kms:kme,jms:jme), tke_tend(ims:ime,kms:kme,jms:jme) )
   ALLOCATE( th_phy(ims:ime,kms:kme,jms:jme), pi_phy(ims:ime,kms:kme,jms:jme), p_phy(ims:ime,kms:kme,jms:jme) )
   ALLOCATE( dz8w(ims:ime,kms:kme,jms:jme), p8w(ims:ime,kms:kme,jms:jme), t8w(ims:ime,kms:kme,jms:jme) )
   ALLOCATE( mu_tendf(ims:ime,jms:jme), psim(ims:ime,jms:jme), psih(ims:ime,jms:jme), gz1oz0(ims:ime,jms:jme) )
   ALLOCATE( chklowq(ims:ime,jms:jme), hol(ims:ime,jms:jme), cu_act_flag(ims:ime,jms:jme) )
   moist_tend = 0. ; moist_old = 0. ; chem_tend = 0. ; tracer_tend = 0. ; scalar_tend = 0.
   ru_tendf = 0. ; rv_tendf = 0. ; rw_tendf = 0. ; t_tendf = 0. ; ph_tendf = 0. ; tke_tend = 0.
   th_phy = 0. ; pi_phy = 0. ; p_phy = 0. ; dz8w = 0. ; p8w = 0. ; t8w = 0.
   mu_tendf = 0. ; psim = 0. ; psih = 0. ; gz1oz0 = 0. ; chklowq = 0. ; hol = 0. ; cu_act_flag = .FALSE.

   CALL set_tiles ( ZONE_SOLVE_EM, grid , ids , ide , jds , jde , ips , ipe , jps , jpe )
   f_flux = config_flags%do_avgflx_cugd .EQ. 1
   adapt_step_flag = .FALSE.
   restart_flag = .FALSE.

   OPEN(newunit=unit_out, FILE='suite_ref.bin', ACCESS='stream', FORM='unformatted', STATUS='replace', &
        CONVERT='little_endian')

   DO ncall = 1, 3
      ! ---- state before the call: the dynamics of the previous step is replaced by a
      !      deterministic change (call 1 = the start_domain state) ----
      IF ( ncall > 1 ) THEN
         DO j = jds, jde-1
         DO k = kds, kde-1
         DO i = ids, ide-1
            grid%t_2(i,k,j) = grid%t_2(i,k,j) + 0.4*REAL(ncall-1)*COS(0.3*i+0.2*j+0.15*k)
            grid%moist(i,k,j,P_QV) = grid%moist(i,k,j,P_QV)*(1.-0.03*REAL(ncall-1)*SIN(0.5*i-0.2*j))
         ENDDO
         ENDDO
         ENDDO
         DO j = jds, jde-1
         DO k = kds, kde-1
         DO i = ids, ide
            grid%u_2(i,k,j) = grid%u_2(i,k,j)*(1.+0.05*REAL(ncall-1)) + 0.3*SIN(0.4*j+0.1*k)
         ENDDO
         ENDDO
         ENDDO
         DO j = jds, jde
         DO k = kds, kde-1
         DO i = ids, ide-1
            grid%v_2(i,k,j) = grid%v_2(i,k,j)*(1.-0.04*REAL(ncall-1)) - 0.2*COS(0.3*i)
         ENDDO
         ENDDO
         ENDDO
         DO j = jds, jde-1
         DO i = ids, ide-1
            grid%mu_2(i,j) = grid%mu_2(i,j) + 15.*REAL(ncall-1)*SIN(0.2*i+0.3*j)
         ENDDO
         ENDDO
      ENDIF
      DO j = jds, jde-1
      DO i = ids, ide-1
         grid%mut(i,j)  = grid%mub(i,j) + grid%mu_2(i,j)
         grid%muts(i,j) = grid%mut(i,j) + 8.*COS(0.25*i-0.15*j)      ! the small steps' new mu
         ! precipitation of the step (microphysics/cumulus are off): exercises rain, snow and
         ! freezing rain on land, snow, ice and glacier points
         grid%rainncv(i,j) = MAX(0., 0.6*SIN(0.35*i+0.45*j+0.7*ncall))
         grid%raincv(i,j)  = 0.
      ENDDO
      ENDDO
      DO j = jds, jde-1
      DO k = kds, kde-1
      DO i = ids, ide-1
         grid%t_1(i,k,j) = grid%t_2(i,k,j)
      ENDDO
      ENDDO
      ENDDO
      ! solve_em: initialize_moist_old (RK stage 1: moist_old = moist for P_QV)
      DO j = jds, jde-1
      DO k = kds, kde-1
      DO i = ids, ide-1
         moist_old(i,k,j,P_QV) = grid%moist(i,k,j,P_QV)
      ENDDO
      ENDDO
      ENDDO
      CALL dump_state('IN', ncall)

      ! ---- solve_em ----
      grid%itimestep = grid%itimestep + 1
      curr_secs = REAL(grid%itimestep-1)*grid%dt
      CALL first_rk_step_part1 ( grid, config_flags, grid%moist, moist_tend, grid%chem, chem_tend, &
                                 grid%tracer, tracer_tend, grid%scalar, scalar_tend, grid%fdda3d, grid%fdda2d, &
                                 grid%aerod, ru_tendf, rv_tendf, rw_tendf, t_tendf, ph_tendf, mu_tendf, &
                                 tke_tend, config_flags%use_adaptive_time_step, curr_secs, &
                                 psim, psih, gz1oz0, chklowq, cu_act_flag, hol, th_phy, &
                                 pi_phy, p_phy, grid%t_phy, dz8w, p8w, t8w, &
                                 ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                                 ips, ipe, jps, jpe, kps, kpe, imsx, imex, jmsx, jmex, kmsx, kmex, &
                                 ipsx, ipex, jpsx, jpex, kpsx, kpex, imsy, imey, jmsy, jmey, kmsy, kmey, &
                                 ipsy, ipey, jpsy, jpey, kpsy, kpey, k_start, k_end, &
                                 f_flux=f_flux, aerocu=grid%aerocu, restart_flag=restart_flag )
      CALL dump_part1(ncall)
      CALL first_rk_step_part2 ( grid, config_flags, grid%moist, moist_old, moist_tend, grid%chem, chem_tend, &
                                 grid%tracer, tracer_tend, grid%scalar, scalar_tend, grid%fdda3d, grid%fdda2d, &
                                 ru_tendf, rv_tendf, rw_tendf, t_tendf, ph_tendf, mu_tendf, tke_tend, &
                                 adapt_step_flag, curr_secs, psim, psih, gz1oz0, chklowq, cu_act_flag, hol, &
                                 th_phy, pi_phy, p_phy, grid%t_phy, dz8w, p8w, t8w, &
                                 grid%nba_mij, num_nba_mij, grid%nba_rij, num_nba_rij, &
                                 ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                                 ips, ipe, jps, jpe, kps, kpe, imsx, imex, jmsx, jmex, kmsx, kmex, &
                                 ipsx, ipex, jpsx, jpex, kpsx, kpex, imsy, imey, jmsy, jmey, kmsy, kmey, &
                                 ipsy, ipey, jpsy, jpey, kpsy, kpey, k_start, k_end )
      CALL dump_part2(ncall)
      CALL phy_prep_part2 ( config_flags, grid%muts, grid%muus, grid%muvs, &
                            grid%c1h, grid%c2h, grid%c1f, grid%c2f, grid%rthraten, &
                            grid%rthblten, grid%rublten, grid%rvblten, &
                            grid%rqvblten, grid%rqcblten, grid%rqiblten, &
                            grid%rucuten,  grid%rvcuten,  grid%rthcuten, &
                            grid%rqvcuten, grid%rqccuten, grid%rqrcuten, &
                            grid%rqicuten, grid%rqscuten, &
                            grid%rushten,  grid%rvshten,  grid%rthshten, &
                            grid%rqvshten, grid%rqcshten, grid%rqrshten, &
                            grid%rqishten, grid%rqsshten, grid%rqgshten, &
                            grid%rthften,  grid%rqvften, &
                            grid%RUNDGDTEN, grid%RVNDGDTEN, grid%RTHNDGDTEN, &
                            grid%RPHNDGDTEN,grid%RQVNDGDTEN, grid%RMUNDGDTEN, &
                            grid%t_2, th_phy, grid%moist(:,:,:,P_QV), &
                            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                            grid%i_start(1), grid%i_end(1), grid%j_start(1), grid%j_end(1), k_start, k_end )
      CALL dump_after(ncall)
   END DO
   CLOSE(unit_out)
   PRINT *, 'suite_ref: 3 physics calls done, itimestep =', grid%itimestep, ' num_tiles =', grid%num_tiles

CONTAINS

   SUBROUTINE put_input(nm, d, b)
      CHARACTER(LEN=*), INTENT(IN) :: nm
      INTEGER, INTENT(IN) :: d(3)
      REAL, INTENT(IN) :: b(:)
      INTEGER :: ii, jj, kk, m
      SELECT CASE (nm)
      CASE ('ZNW') ; grid%znw(1:d(1)) = b
      CASE ('ZNU') ; grid%znu(1:d(1)) = b
      CASE ('DNW') ; grid%dnw(1:d(1)) = b
      CASE ('RDNW'); grid%rdnw(1:d(1)) = b
      CASE ('DN')  ; grid%dn(1:d(1)) = b
      CASE ('RDN') ; grid%rdn(1:d(1)) = b
      CASE ('FNM') ; grid%fnm(1:d(1)) = b
      CASE ('FNP') ; grid%fnp(1:d(1)) = b
      CASE ('ZS')  ; grid%zs(1:d(1)) = b
      CASE ('DZS') ; grid%dzs(1:d(1)) = b
      CASE DEFAULT
         m = 0
         DO jj = 1, d(3)
         DO kk = 1, d(2)
         DO ii = 1, d(1)
            m = m + 1
            CALL put_one(nm, ii, kk, jj, d, b(m))
         ENDDO
         ENDDO
         ENDDO
      END SELECT
   END SUBROUTINE put_input

   ! 2D fields arrive as (i,j,1): kk carries j, jj = 1
   SUBROUTINE put_one(nm, ii, kk, jj, d, v)
      CHARACTER(LEN=*), INTENT(IN) :: nm
      INTEGER, INTENT(IN) :: ii, kk, jj, d(3)
      REAL, INTENT(IN) :: v
      SELECT CASE (nm)
      CASE ('HT')      ; grid%ht(ii,kk) = v
      CASE ('MUB')     ; grid%mub(ii,kk) = v
      CASE ('MU_2')    ; grid%mu_2(ii,kk) = v
      CASE ('XLAND')   ; grid%xland(ii,kk) = v
      CASE ('LANDMASK'); grid%landmask(ii,kk) = v
      CASE ('XICE')    ; grid%xice(ii,kk) = v
      CASE ('IVGTYP')  ; grid%ivgtyp(ii,kk) = NINT(v)
      CASE ('ISLTYP')  ; grid%isltyp(ii,kk) = NINT(v)
      CASE ('LU_INDEX'); grid%lu_index(ii,kk) = v
      CASE ('VEGFRA')  ; grid%vegfra(ii,kk) = v
      CASE ('SHDMAX')  ; grid%shdmax(ii,kk) = v
      CASE ('SHDMIN')  ; grid%shdmin(ii,kk) = v
      CASE ('ALBBCK')  ; grid%albbck(ii,kk) = v
      CASE ('SNOALB')  ; grid%snoalb(ii,kk) = v
      CASE ('TMN')     ; grid%tmn(ii,kk) = v
      CASE ('TSK')     ; grid%tsk(ii,kk) = v
      CASE ('SST')     ; grid%sst(ii,kk) = v
      CASE ('SNOW')    ; grid%snow(ii,kk) = v
      CASE ('SNOWH')   ; grid%snowh(ii,kk) = v
      CASE ('CANWAT')  ; grid%canwat(ii,kk) = v
      CASE ('XLAT')    ; grid%xlat(ii,kk) = v
      CASE ('XLONG')   ; grid%xlong(ii,kk) = v
      CASE ('GLW')     ; grid%glw(ii,kk) = v
      CASE ('SWDOWN')  ; grid%swdown(ii,kk) = v
      CASE ('GSW')     ; grid%gsw(ii,kk) = v
      CASE ('U')       ; grid%u_2(ii,kk,jj) = v
      CASE ('V')       ; grid%v_2(ii,kk,jj) = v
      CASE ('W')       ; grid%w_2(ii,kk,jj) = v
      CASE ('T')       ; grid%t_2(ii,kk,jj) = v
      CASE ('QVAPOR')  ; grid%moist(ii,kk,jj,P_QV) = v
      CASE ('PH')      ; grid%ph_2(ii,kk,jj) = v
      CASE ('TSLB')    ; grid%tslb(ii,kk,jj) = v
      CASE ('SMOIS')   ; grid%smois(ii,kk,jj) = v
      CASE ('SH2O')    ; grid%sh2o(ii,kk,jj) = v
      CASE DEFAULT
         CALL wrf_error_fatal('suite_ref: unknown input field '//nm)
      END SELECT
   END SUBROUTINE put_one

   SUBROUTINE wrec(tag, n, nm, nd_, d, a)
      CHARACTER(LEN=*), INTENT(IN) :: tag, nm
      INTEGER, INTENT(IN) :: n, nd_, d(3)
      REAL, INTENT(IN) :: a(:)
      CHARACTER(LEN=32) :: full
      WRITE(full, '(A,I1,A,A)') tag, n, ':', nm
      WRITE(unit_out) full, nd_, d, a
   END SUBROUTINE wrec

   SUBROUTINE w2(tag, n, nm, a, ie, je)          ! mass 2D (ids:ie, jds:je)
      CHARACTER(LEN=*), INTENT(IN) :: tag, nm
      INTEGER, INTENT(IN) :: n, ie, je
      REAL, INTENT(IN) :: a(ims:ime,jms:jme)
      CALL wrec(tag, n, nm, 2, (/ie, je, 1/), RESHAPE(a(ids:ie,jds:je), (/ie*je/)))
   END SUBROUTINE w2

   SUBROUTINE w2i(tag, n, nm, a)
      CHARACTER(LEN=*), INTENT(IN) :: tag, nm
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: a(ims:ime,jms:jme)
      CALL w2(tag, n, nm, REAL(a), ide-1, jde-1)
   END SUBROUTINE w2i

   SUBROUTINE w3(tag, n, nm, a, ie, ke, je, k1)  ! 3D (ids:ie, k1:ke, jds:je) in (i,k,j) order
      CHARACTER(LEN=*), INTENT(IN) :: tag, nm
      INTEGER, INTENT(IN) :: n, ie, ke, je, k1
      REAL, INTENT(IN) :: a(:,:,:)
      INTEGER :: nk
      nk = ke - k1 + 1
      CALL wrec(tag, n, nm, 3, (/ie, nk, je/), RESHAPE(a(1:ie,1:nk,1:je), (/ie*nk*je/)))
   END SUBROUTINE w3

   SUBROUTINE m3(tag, n, nm, a)                   ! mass 3D (kds:kde-1)
      CHARACTER(LEN=*), INTENT(IN) :: tag, nm
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN) :: a(ims:ime,kms:kme,jms:jme)
      CALL w3(tag, n, nm, a(ids:ide-1,kds:kde-1,jds:jde-1), ide-1, kde-1, jde-1, kds)
   END SUBROUTINE m3

   SUBROUTINE wl3(tag, n, nm, a)                  ! w-level 3D (kds:kde)
      CHARACTER(LEN=*), INTENT(IN) :: tag, nm
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN) :: a(ims:ime,kms:kme,jms:jme)
      CALL w3(tag, n, nm, a(ids:ide-1,kds:kde,jds:jde-1), ide-1, kde, jde-1, kds)
   END SUBROUTINE wl3

   SUBROUTINE s3(tag, n, nm, a)                   ! soil 3D
      CHARACTER(LEN=*), INTENT(IN) :: tag, nm
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN) :: a(ims:ime,1:config_flags%num_soil_layers,jms:jme)
      CALL w3(tag, n, nm, a(ids:ide-1,:,jds:jde-1), ide-1, config_flags%num_soil_layers, jde-1, 1)
   END SUBROUTINE s3

   SUBROUTINE w1(tag, n, nm, a, nn)
      CHARACTER(LEN=*), INTENT(IN) :: tag, nm
      INTEGER, INTENT(IN) :: n, nn
      REAL, INTENT(IN) :: a(:)
      CALL wrec(tag, n, nm, 1, (/nn, 1, 1/), a(1:nn))
   END SUBROUTINE w1

   ! the land-surface / PBL state the suite carries between calls (inputs and outputs)
   SUBROUTINE dump_surface(tag, n)
      CHARACTER(LEN=*), INTENT(IN) :: tag
      INTEGER, INTENT(IN) :: n
      INTEGER :: ie, je
      ie = ide-1 ; je = jde-1
      CALL w2(tag,n,'TSK',grid%tsk,ie,je)       ; CALL w2(tag,n,'HFX',grid%hfx,ie,je)
      CALL w2(tag,n,'QFX',grid%qfx,ie,je)       ; CALL w2(tag,n,'LH',grid%lh,ie,je)
      CALL w2(tag,n,'UST',grid%ust,ie,je)       ; CALL w2(tag,n,'ZNT',grid%znt,ie,je)
      CALL w2(tag,n,'Z0',grid%z0,ie,je)         ; CALL w2(tag,n,'MOL',grid%mol,ie,je)
      CALL w2(tag,n,'RMOL',grid%rmol,ie,je)     ; CALL w2(tag,n,'PBLH',grid%pblh,ie,je)
      CALL w2i(tag,n,'KPBL',grid%kpbl)          ; CALL w2(tag,n,'QSFC',grid%qsfc,ie,je)
      CALL w2(tag,n,'FLHC',grid%flhc,ie,je)
      CALL w2(tag,n,'FLQC',grid%flqc,ie,je)     ; CALL w2(tag,n,'T2',grid%t2,ie,je)
      CALL w2(tag,n,'TH2',grid%th2,ie,je)       ; CALL w2(tag,n,'Q2',grid%q2,ie,je)
      CALL w2(tag,n,'U10',grid%u10,ie,je)       ; CALL w2(tag,n,'V10',grid%v10,ie,je)
      CALL w2(tag,n,'BR',grid%br,ie,je)         ; CALL w2(tag,n,'WSPD',grid%wspd,ie,je)
      CALL w2(tag,n,'FM',grid%fm,ie,je)         ; CALL w2(tag,n,'FH',grid%fh,ie,je)
      CALL w2(tag,n,'ZOL',grid%zol,ie,je)       ; CALL w2(tag,n,'REGIME',grid%regime,ie,je)
      CALL w2(tag,n,'PSFC',grid%psfc,ie,je)     ; CALL w2(tag,n,'GRDFLX',grid%grdflx,ie,je)
      CALL w2(tag,n,'ALBEDO',grid%albedo,ie,je) ; CALL w2(tag,n,'ALBBCK',grid%albbck,ie,je)
      CALL w2(tag,n,'EMISS',grid%emiss,ie,je)   ; CALL w2(tag,n,'EMBCK',grid%embck,ie,je)
      CALL w2(tag,n,'SNOW',grid%snow,ie,je)     ; CALL w2(tag,n,'SNOWH',grid%snowh,ie,je)
      CALL w2(tag,n,'SNOWC',grid%snowc,ie,je)   ; CALL w2(tag,n,'CANWAT',grid%canwat,ie,je)
      CALL w2(tag,n,'SMSTAV',grid%smstav,ie,je) ; CALL w2(tag,n,'SMSTOT',grid%smstot,ie,je)
      CALL w2(tag,n,'SFCRUNOFF',grid%sfcrunoff,ie,je) ; CALL w2(tag,n,'UDRUNOFF',grid%udrunoff,ie,je)
      CALL w2(tag,n,'ACSNOM',grid%acsnom,ie,je) ; CALL w2(tag,n,'ACSNOW',grid%acsnow,ie,je)
      CALL w2(tag,n,'SNOPCX',grid%snopcx,ie,je) ; CALL w2(tag,n,'POTEVP',grid%potevp,ie,je)
      CALL w2(tag,n,'NOAHRES',grid%noahres,ie,je) ; CALL w2(tag,n,'SFCEVP',grid%sfcevp,ie,je)
      CALL w2(tag,n,'SFCEXC',grid%sfcexc,ie,je) ; CALL w2(tag,n,'ACHFX',grid%achfx,ie,je)
      CALL w2(tag,n,'ACLHF',grid%aclhf,ie,je)   ; CALL w2(tag,n,'ACGRDFLX',grid%acgrdflx,ie,je)
      CALL w2(tag,n,'LAI',grid%lai,ie,je)       ; CALL w2(tag,n,'SNOTIME',grid%snotime,ie,je)
      CALL w2(tag,n,'QZ0',grid%qz0,ie,je)       ; CALL w2(tag,n,'RAINBL',grid%rainbl,ie,je)
      CALL w2(tag,n,'XLAND',grid%xland,ie,je)   ; CALL w2(tag,n,'XICE',grid%xice,ie,je)
      CALL w2(tag,n,'MAVAIL',grid%mavail,ie,je) ; CALL w2(tag,n,'TMN',grid%tmn,ie,je)
      CALL w2(tag,n,'VEGFRA',grid%vegfra,ie,je) ; CALL w2(tag,n,'SHDMAX',grid%shdmax,ie,je)
      CALL w2(tag,n,'SHDMIN',grid%shdmin,ie,je) ; CALL w2(tag,n,'SNOALB',grid%snoalb,ie,je)
      CALL w2(tag,n,'SST',grid%sst,ie,je)       ; CALL w2(tag,n,'LAKEMASK',grid%lakemask,ie,je)
      CALL w2(tag,n,'DX2D',grid%dx2d,ie,je)     ; CALL w2(tag,n,'UOCE',grid%uoce,ie,je)
      CALL w2(tag,n,'VOCE',grid%voce,ie,je)     ; CALL w2(tag,n,'CTOPO',grid%ctopo,ie,je)
      CALL w2(tag,n,'CTOPO2',grid%ctopo2,ie,je) ; CALL w2(tag,n,'WATER_DEPTH',grid%water_depth,ie,je)
      CALL w2(tag,n,'WSTAR_YSU',grid%wstar_ysu,ie,je) ; CALL w2(tag,n,'DELTA_YSU',grid%delta_ysu,ie,je)
      CALL w2(tag,n,'USTM',grid%ustm,ie,je)     ; CALL w2(tag,n,'CK',grid%ck,ie,je)
      CALL w2(tag,n,'CKA',grid%cka,ie,je)       ; CALL w2(tag,n,'CD',grid%cd,ie,je)
      CALL w2(tag,n,'CDA',grid%cda,ie,je)       ; CALL w2(tag,n,'SR',grid%sr,ie,je)
      CALL w2(tag,n,'HOL',hol,ie,je)
      CALL w2i(tag,n,'IVGTYP',grid%ivgtyp)      ; CALL w2i(tag,n,'ISLTYP',grid%isltyp)
      CALL s3(tag,n,'TSLB',grid%tslb)           ; CALL s3(tag,n,'SMOIS',grid%smois)
      CALL s3(tag,n,'SH2O',grid%sh2o)           ; CALL s3(tag,n,'SMCREL',grid%smcrel)
      CALL w1(tag,n,'DZS',grid%dzs,config_flags%num_soil_layers)
      CALL w1(tag,n,'ZS',grid%zs,config_flags%num_soil_layers)
      CALL m3(tag,n,'EXCH_H',grid%exch_h)       ; CALL m3(tag,n,'EXCH_M',grid%exch_m)
      CALL m3(tag,n,'RUBLTEN',grid%rublten)     ; CALL m3(tag,n,'RVBLTEN',grid%rvblten)
      CALL m3(tag,n,'RTHBLTEN',grid%rthblten)   ; CALL m3(tag,n,'RQVBLTEN',grid%rqvblten)
      CALL m3(tag,n,'RTHRATEN',grid%rthraten)
   END SUBROUTINE dump_surface

   SUBROUTINE dump_state(tag, n)
      CHARACTER(LEN=*), INTENT(IN) :: tag
      INTEGER, INTENT(IN) :: n
      INTEGER :: ie, je
      ie = ide-1 ; je = jde-1
      CALL dump_surface(tag, n)
      CALL w3(tag,n,'U',grid%u_2(ids:ide,kds:kde-1,jds:jde-1),ide,kde-1,je,1)
      CALL w3(tag,n,'V',grid%v_2(ids:ide-1,kds:kde-1,jds:jde),ie,kde-1,jde,1)
      CALL m3(tag,n,'T',grid%t_2)               ; CALL m3(tag,n,'T_1',grid%t_1)
      CALL m3(tag,n,'P',grid%p)                 ; CALL m3(tag,n,'PB',grid%pb)
      CALL m3(tag,n,'ALT',grid%alt)             ; CALL m3(tag,n,'QVAPOR',grid%moist(:,:,:,P_QV))
      CALL m3(tag,n,'QV_OLD',moist_old(:,:,:,P_QV))
      CALL wl3(tag,n,'PH',grid%ph_2)            ; CALL wl3(tag,n,'PHB',grid%phb)
      CALL w2(tag,n,'MUT',grid%mut,ie,je)       ; CALL w2(tag,n,'MUTS',grid%muts,ie,je)
      CALL w2(tag,n,'RAINNCV',grid%rainncv,ie,je) ; CALL w2(tag,n,'RAINCV',grid%raincv,ie,je)
      CALL w2(tag,n,'GLW',grid%glw,ie,je)       ; CALL w2(tag,n,'SWDOWN',grid%swdown,ie,je)
      CALL w2(tag,n,'GSW',grid%gsw,ie,je)       ; CALL w2(tag,n,'HT',grid%ht,ie,je)
      CALL w1(tag,n,'C1H',grid%c1h,kde-1)       ; CALL w1(tag,n,'C2H',grid%c2h,kde-1)
      CALL w1(tag,n,'FNM',grid%fnm,kde-1)       ; CALL w1(tag,n,'FNP',grid%fnp,kde-1)
      CALL w1(tag,n,'DNW',grid%dnw,kde-1)       ; CALL w1(tag,n,'ZNW',grid%znw,kde)
      CALL w1(tag,n,'P_TOP',(/grid%p_top/),1)   ; CALL w1(tag,n,'DT',(/grid%dt/),1)
      CALL w1(tag,n,'DX',(/grid%dx/),1)         ; CALL w1(tag,n,'ITIMESTEP',(/REAL(grid%itimestep)/),1)
   END SUBROUTINE dump_state

   SUBROUTINE dump_part1(n)                      ! phy_prep outputs (locals of solve_em)
      INTEGER, INTENT(IN) :: n
      CALL m3('P1_',n,'TH_PHY',th_phy) ; CALL m3('P1_',n,'P_PHY',p_phy) ; CALL m3('P1_',n,'PI_PHY',pi_phy)
      CALL m3('P1_',n,'T_PHY',grid%t_phy) ; CALL m3('P1_',n,'RHO',grid%rho)
      CALL m3('P1_',n,'U_PHY',grid%u_phy) ; CALL m3('P1_',n,'V_PHY',grid%v_phy)
      CALL m3('P1_',n,'Z',grid%z) ; CALL wl3('P1_',n,'Z_AT_W',grid%z_at_w)
      CALL wl3('P1_',n,'DZ8W',dz8w) ; CALL wl3('P1_',n,'P8W',p8w) ; CALL wl3('P1_',n,'T8W',t8w)
      CALL m3('P1_',n,'P_HYD',grid%p_hyd) ; CALL wl3('P1_',n,'P_HYD_W',grid%p_hyd_w)
      CALL w2('P1_',n,'PSIM',psim,ide-1,jde-1) ; CALL w2('P1_',n,'PSIH',psih,ide-1,jde-1)
      CALL w2('P1_',n,'GZ1OZ0',gz1oz0,ide-1,jde-1) ; CALL w2('P1_',n,'CHKLOWQ',chklowq,ide-1,jde-1)
      CALL dump_surface('P1_', n)
   END SUBROUTINE dump_part1

   SUBROUTINE dump_part2(n)                      ! coupled R*TEN and the *_tendf fields
      INTEGER, INTENT(IN) :: n
      CALL w3('P2_',n,'RU_TENDF',ru_tendf(ids:ide,kds:kde-1,jds:jde-1),ide,kde-1,jde-1,1)
      CALL w3('P2_',n,'RV_TENDF',rv_tendf(ids:ide-1,kds:kde-1,jds:jde),ide-1,kde-1,jde,1)
      CALL m3('P2_',n,'T_TENDF',t_tendf) ; CALL m3('P2_',n,'QV_TENDF',moist_tend(:,:,:,P_QV))
      CALL m3('P2_',n,'RUBLTEN',grid%rublten) ; CALL m3('P2_',n,'RVBLTEN',grid%rvblten)
      CALL m3('P2_',n,'RTHBLTEN',grid%rthblten) ; CALL m3('P2_',n,'RQVBLTEN',grid%rqvblten)
   END SUBROUTINE dump_part2

   SUBROUTINE dump_after(n)
      INTEGER, INTENT(IN) :: n
      CALL dump_surface('OUT', n)
   END SUBROUTINE dump_after

END PROGRAM suite_ref
