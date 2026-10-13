program read_tables_ref
  implicit none

  character(len=256) :: tbl_dir, outfile
  integer :: u_lu, u_veg, u_soil, u_gen, ounit
  character(len=256) :: lutype
  integer :: lucats, luseas, ls, lc, li, isn, n
  integer, parameter :: max_cats = 65, max_seas = 2
  real(8) :: albd(max_cats, max_seas), slmo(max_cats, max_seas)
  real(8) :: sfem(max_cats, max_seas), sfz0(max_cats, max_seas)
  real(8) :: therin(max_cats, max_seas), scfx(max_cats, max_seas), sfhc(max_cats, max_seas)

  ! Noah tables
  integer :: slpcats, slcats, veg_lucats, bare, natural
  real(8) :: slope_data(30), sbeta, fxexp, csoil, salp, refdk, refkdt, frzk, zbot, czil, smlow, smhigh, lvcoef
  real(8) :: bb(30), drysmc(30), f11(30), maxsmc(30), refsmc(30), satpsi(30), satdk(30), satdw(30), wltsmc(30), qtz(30)
  real(8) :: topt, cmcmax, cfactr, rsmax
  integer :: nrotbl(50)
  real(8) :: snuptbl(50), rstbl(50), rgltbl(50), hstbl(50), shdtbl(50), maxalb(50)
  real(8) :: emissmintbl(50), emissmaxtbl(50), laimintbl(50), laimaxtbl(50)
  real(8) :: z0mintbl(50), z0maxtbl(50), albedomintbl(50), albedomaxtbl(50)
  real(8) :: ztopvtbl(50), zbotvtbl(50)

  character(len=256) :: line, token
  integer :: ierr, i, k

  if (command_argument_count() >= 2) then
    call get_command_argument(1, tbl_dir)
    call get_command_argument(2, outfile)
  else
    tbl_dir = "."
    outfile = "tables_ref.bin"
  end if

  ! 1. Read LANDUSE.TBL for MODIFIED_IGBP_MODIS_NOAH
  open(newunit=u_lu, file=trim(tbl_dir)//'/LANDUSE.TBL', status='old', form='formatted')
  find_lu: do
    read(u_lu, '(A)', iostat=ierr) line
    if (ierr /= 0) stop "LANDUSE dataset not found"
    line = adjustl(line)
    if (trim(line) == 'MODIFIED_IGBP_MODIS_NOAH') exit find_lu
  end do find_lu

  read(u_lu, *) lucats, luseas
  albd = 0.0d0; slmo = 0.0d0; sfem = 0.0d0; sfz0 = 0.0d0
  therin = 0.0d0; scfx = 0.0d0; sfhc = 0.0d0

  do ls = 1, luseas
    read(u_lu, *) ! season line (SUMMER, WINTER)
    do lc = 1, lucats
      read(u_lu, *) li, albd(lc, ls), slmo(lc, ls), sfem(lc, ls), sfz0(lc, ls), therin(lc, ls), scfx(lc, ls), sfhc(lc, ls)
    end do
  end do
  close(u_lu)

  ! 2. Read GENPARM.TBL
  open(newunit=u_gen, file=trim(tbl_dir)//'/GENPARM.TBL', status='old', form='formatted')
  find_slope: do
    read(u_gen, '(A)', iostat=ierr) line
    if (ierr /= 0) stop "SLOPE_DATA not found"
    if (trim(adjustl(line)) == 'SLOPE_DATA') exit find_slope
  end do find_slope
  read(u_gen, *) slpcats
  slope_data = 0.0d0
  do i = 1, slpcats
    read(u_gen, *) slope_data(i)
  end do
  find_gen: do
    read(u_gen, '(A)', iostat=ierr) line
    if (ierr /= 0) exit find_gen
    line = adjustl(line)
    if (line(1:5) == 'SBETA') read(u_gen, *) sbeta
    if (line(1:5) == 'FXEXP') read(u_gen, *) fxexp
    if (line(1:5) == 'CSOIL') read(u_gen, *) csoil
    if (line(1:4) == 'SALP') read(u_gen, *) salp
    if (line(1:6) == 'REFKDT') then
      read(u_gen, *) refkdt
    else if (line(1:5) == 'REFDK') then
      read(u_gen, *) refdk
    end if
    if (line(1:4) == 'FRZK') read(u_gen, *) frzk
    if (line(1:4) == 'ZBOT') read(u_gen, *) zbot
    if (line(1:4) == 'CZIL') read(u_gen, *) czil
    if (line(1:5) == 'SMLOW') read(u_gen, *) smlow
    if (line(1:6) == 'SMHIGH') read(u_gen, *) smhigh
    if (line(1:6) == 'LVCOEF') read(u_gen, *) lvcoef
  end do find_gen
  close(u_gen)

  ! 3. Read SOILPARM.TBL (STAS)
  open(newunit=u_soil, file=trim(tbl_dir)//'/SOILPARM.TBL', status='old', form='formatted')
  find_soil: do
    read(u_soil, '(A)', iostat=ierr) line
    if (ierr /= 0) stop "Soil dataset STAS not found"
    if (trim(adjustl(line)) == 'Soil Parameters') then
      read(u_soil, '(A)') line
      if (trim(adjustl(line)) == 'STAS') exit find_soil
    end if
  end do find_soil
  read(u_soil, *) slcats
  bb = 0.0d0; drysmc = 0.0d0; f11 = 0.0d0; maxsmc = 0.0d0; refsmc = 0.0d0
  satpsi = 0.0d0; satdk = 0.0d0; satdw = 0.0d0; wltsmc = 0.0d0; qtz = 0.0d0
  do i = 1, slcats
    read(u_soil, *) k, bb(i), drysmc(i), f11(i), maxsmc(i), refsmc(i), satpsi(i), satdk(i), satdw(i), wltsmc(i), qtz(i)
  end do
  close(u_soil)

  ! 4. Read VEGPARM.TBL (MODIFIED_IGBP_MODIS_NOAH)
  open(newunit=u_veg, file=trim(tbl_dir)//'/VEGPARM.TBL', status='old', form='formatted')
  find_veg: do
    read(u_veg, '(A)', iostat=ierr) line
    if (ierr /= 0) stop "VEGPARM dataset MODIFIED_IGBP_MODIS_NOAH not found"
    if (trim(adjustl(line)) == 'MODIFIED_IGBP_MODIS_NOAH') exit find_veg
  end do find_veg
  read(u_veg, *) veg_lucats
  do i = 1, veg_lucats
    read(u_veg, *) k, shdtbl(i), nrotbl(i), rstbl(i), rgltbl(i), hstbl(i), &
                   snuptbl(i), maxalb(i), laimintbl(i), laimaxtbl(i), &
                   emissmintbl(i), emissmaxtbl(i), albedomintbl(i), albedomaxtbl(i), &
                   z0mintbl(i), z0maxtbl(i), ztopvtbl(i), zbotvtbl(i)
  end do
  find_veg_scalars: do
    read(u_veg, '(A)', iostat=ierr) line
    if (ierr /= 0) exit find_veg_scalars
    line = adjustl(line)
    if (trim(line) == 'TOPT_DATA') read(u_veg, *) topt
    if (trim(line) == 'CMCMAX_DATA') read(u_veg, *) cmcmax
    if (trim(line) == 'CFACTR_DATA') read(u_veg, *) cfactr
    if (trim(line) == 'RSMAX_DATA') read(u_veg, *) rsmax
    if (trim(line) == 'BARE') read(u_veg, *) bare
    if (trim(line) == 'NATURAL') then
      read(u_veg, *) natural
      exit find_veg_scalars
    end if
  end do find_veg_scalars
  close(u_veg)

  ! Write output binary
  open(newunit=ounit, file=trim(outfile), status='replace', access='stream', form='unformatted')
  ! Landuse
  write(ounit) lucats, luseas
  write(ounit) albd(1:lucats, 1:luseas)
  write(ounit) slmo(1:lucats, 1:luseas)
  write(ounit) sfem(1:lucats, 1:luseas)
  write(ounit) sfz0(1:lucats, 1:luseas)
  write(ounit) therin(1:lucats, 1:luseas)
  write(ounit) scfx(1:lucats, 1:luseas)
  write(ounit) sfhc(1:lucats, 1:luseas)

  ! Noah GENPARM
  write(ounit) slpcats
  write(ounit) slope_data
  write(ounit) sbeta, fxexp, csoil, salp, refdk, refkdt, frzk, zbot, czil, smlow, smhigh, lvcoef

  ! Noah SOILPARM
  write(ounit) slcats
  write(ounit) bb, drysmc, f11, maxsmc, refsmc, satpsi, satdk, satdw, wltsmc, qtz

  ! Noah VEGPARM
  write(ounit) veg_lucats, bare, natural
  write(ounit) topt, cmcmax, cfactr, rsmax
  write(ounit) nrotbl, snuptbl, rstbl, rgltbl, hstbl, shdtbl, maxalb
  write(ounit) emissmintbl, emissmaxtbl, laimintbl, laimaxtbl
  write(ounit) z0mintbl, z0maxtbl, albedomintbl, albedomaxtbl
  write(ounit) ztopvtbl, zbotvtbl

  close(ounit)
  print *, "Successfully wrote ", trim(outfile)

end program read_tables_ref
