program driver_noah_ice_fortran
  ! Calls WRF's SFLX_GLACIAL (module_sf_noahlsm_glacial_only.F) and SFLX_SEAICE
  ! (module_sf_noah_seaice.F) column by column; seaice options as in WRF's defaults
  ! (seaice_albedo_opt = 0, seaice_albedo_default = 0.65, seaice_snowdepth_opt = 0,
  ! seaice_snowdepth_max = 1.E10, seaice_snowdepth_min = 0.001).
  use module_sf_noahlsm, only: lvcoef_data
  use module_sf_noahlsm_glacial_only, only: sflx_glacial
  use module_sf_noah_seaice, only: sflx_seaice
  implicit none

  integer, parameter :: nin = 24
  character(len=256) :: infile, outfile
  integer :: ncols, nsoil, i, k, iunit, ounit
  real(8) :: dt, lvcoef
  real(8), allocatable :: sldpth(:), vin(:,:), stcin(:,:), stc(:)
  real(8), allocatable :: og(:,:), os(:,:)
  real(8) :: ffrozp, lwdn, solnet, sfcprs, prcp, sfctmp, q2, th2, q2sat, dqsdt2, alb, snoalb, tbot, &
             z0brd, emissi, embrd, t1, snowh, sneqv, ch, sncovr, snotime1, ribb, sithick
  real(8) :: z0, albedo, eta, sheat, eta_kinematic, fdown, esnow, dew, etp, ssoil, flx1, flx2, flx3, &
             snomlt, runoff1, q1, albedosi, snowonsi
  real(8), parameter :: zlvl = 10.0
  integer, parameter :: isice = 15
  integer, parameter :: ng = 27, ns = 24

  call get_command_argument(1, infile)
  call get_command_argument(2, outfile)

  open(newunit=iunit, file=trim(infile), status='old', access='stream', form='unformatted')
  read(iunit) ncols, nsoil
  read(iunit) dt, lvcoef
  allocate(sldpth(nsoil), vin(ncols, nin), stcin(ncols, nsoil), stc(nsoil))
  allocate(og(ncols, ng + nsoil), os(ncols, ns + nsoil))
  read(iunit) sldpth
  read(iunit) ((vin(i, k), i = 1, ncols), k = 1, nin)
  read(iunit) ((stcin(i, k), i = 1, ncols), k = 1, nsoil)
  close(iunit)
  lvcoef_data = lvcoef

  do i = 1, ncols
    call load(i)
    stc(:) = stcin(i, :)
    call sflx_glacial(i, 1, isice, ffrozp, dt, zlvl, nsoil, sldpth, &
                      lwdn, solnet, sfcprs, prcp, sfctmp, q2, &
                      th2, q2sat, dqsdt2, &
                      alb, snoalb, tbot, z0brd, z0, emissi, embrd, &
                      t1, stc, snowh, sneqv, albedo, ch, &
                      eta, sheat, eta_kinematic, fdown, &
                      esnow, dew, &
                      etp, ssoil, &
                      flx1, flx2, flx3, &
                      snomlt, sncovr, &
                      runoff1, &
                      q1, &
                      snotime1, &
                      ribb)
    og(i, 1:ng) = (/ t1, snowh, sneqv, sncovr, snotime1, ribb, emissi, albedo, z0, eta, sheat, &
                     eta_kinematic, fdown, esnow, dew, etp, ssoil, flx1, flx2, flx3, snomlt, runoff1, q1, &
                     ch, alb, z0brd, embrd /)
    og(i, ng + 1:ng + nsoil) = stc(:)

    call load(i)
    stc(:) = stcin(i, :)
    albedosi = 0.0
    snowonsi = 0.0
    call sflx_seaice(i, 1, 0, 0.65, &
                     0, 1.E10, &
                     0.001, &
                     ffrozp, dt, zlvl, nsoil, &
                     sithick, &
                     lwdn, solnet, sfcprs, prcp, sfctmp, q2, &
                     th2, q2sat, dqsdt2, &
                     snoalb, tbot, z0brd, z0, emissi, &
                     t1, stc, snowh, sneqv, albedo, ch, &
                     albedosi, snowonsi, &
                     eta, sheat, eta_kinematic, fdown, &
                     esnow, dew, etp, ssoil, flx1, flx2, flx3, &
                     snomlt, sncovr, &
                     runoff1, q1, ribb)
    os(i, 1:ns) = (/ t1, snowh, sneqv, sncovr, ribb, emissi, albedo, z0, z0brd, eta, sheat, &
                     eta_kinematic, fdown, esnow, dew, etp, ssoil, flx1, flx2, flx3, snomlt, runoff1, q1, ch /)
    os(i, ns + 1:ns + nsoil) = stc(:)
  end do

  open(newunit=ounit, file=trim(outfile), status='replace', access='stream', form='unformatted')
  write(ounit) ncols, nsoil
  write(ounit) ((og(i, k), i = 1, ncols), k = 1, ng + nsoil)
  write(ounit) ((os(i, k), i = 1, ncols), k = 1, ns + nsoil)
  close(ounit)

contains

  subroutine load(ic)
    integer, intent(in) :: ic
    ffrozp = vin(ic, 1);  lwdn = vin(ic, 2);    solnet = vin(ic, 3);   sfcprs = vin(ic, 4)
    prcp = vin(ic, 5);    sfctmp = vin(ic, 6);  q2 = vin(ic, 7);       th2 = vin(ic, 8)
    q2sat = vin(ic, 9);   dqsdt2 = vin(ic, 10); alb = vin(ic, 11);     snoalb = vin(ic, 12)
    tbot = vin(ic, 13);   z0brd = vin(ic, 14);  emissi = vin(ic, 15);  embrd = vin(ic, 16)
    t1 = vin(ic, 17);     snowh = vin(ic, 18);  sneqv = vin(ic, 19);   ch = vin(ic, 20)
    sncovr = vin(ic, 21); snotime1 = vin(ic, 22); ribb = vin(ic, 23);  sithick = vin(ic, 24)
  end subroutine load

end program driver_noah_ice_fortran
