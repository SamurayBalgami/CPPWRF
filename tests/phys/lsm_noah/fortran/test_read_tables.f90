program test_read_tables
  use module_sf_noahlsm
  implicit none

  character(len=256) :: vegfile, soilfile, genfile
  character(len=256) :: a_string
  integer :: i, lc, iindex, num_slope, ierr
  character(len=4) :: sltype_in, lutype_in

  vegfile = '../VEGPARM.TBL'
  soilfile = '../SOILPARM.TBL'
  genfile = '../GENPARM.TBL'

  ! Read GENPARM
  open(19, file=trim(genfile), status='old', action='read')
  read(19,*)
  read(19,*)
  read(19,*) num_slope
  slpcats = num_slope
  do lc = 1, slpcats
     read(19,*) slope_data(lc)
  end do
  read(19,*)
  read(19,*) sbeta_data
  read(19,*)
  read(19,*) fxexp_data
  read(19,*)
  read(19,*) csoil_data
  read(19,*)
  read(19,*) salp_data
  read(19,*)
  read(19,*) refdk_data
  read(19,*)
  read(19,*) refkdt_data
  read(19,*)
  read(19,*) frzk_data
  read(19,*)
  read(19,*) zbot_data
  read(19,*)
  read(19,*) czil_data
  read(19,*)
  read(19,*) smlow_data
  read(19,*)
  read(19,*) smhigh_data
  read(19,*)
  read(19,*) lvcoef_data
  close(19)

  print *, 'GENPARM read ok: slpcats=', slpcats, 'sbeta=', sbeta_data, 'csoil=', csoil_data

  ! Read SOILPARM for STAS
  open(19, file=trim(soilfile), status='old', action='read')
  read(19,*)
  read(19,'(A4)') sltype_in
  read(19,*) slcats, iindex
  print *, 'SOILPARM: sltype=', sltype_in, 'slcats=', slcats
  do lc = 1, slcats
     read(19,*) iindex, bb(lc), drysmc(lc), f11(lc), maxsmc(lc), &
                refsmc(lc), satpsi(lc), satdk(lc), satdw(lc),    &
                wltsmc(lc), qtz(lc)
  end do
  close(19)
  print *, 'SOILPARM read ok: cat 1 bb=', bb(1), 'maxsmc=', maxsmc(1), 'satdk=', satdk(1)

  ! Read VEGPARM for USGS
  open(19, file=trim(vegfile), status='old', action='read')
  read(19,'(A)') a_string
  read(19,*) lutype_in
  read(19,*) lucats, iindex
  print *, 'VEGPARM: lutype=', lutype_in, 'lucats=', lucats
  do lc = 1, lucats
     read(19,*) iindex, shdtbl(lc), nrotbl(lc), rstbl(lc), rgltbl(lc), hstbl(lc), &
                snuptbl(lc), maxalb(lc), laimintbl(lc), laimaxtbl(lc), emissmintbl(lc), &
                emissmaxtbl(lc), albedomintbl(lc), albedomaxtbl(lc), z0mintbl(lc), z0maxtbl(lc), &
                ztopvtbl(lc), zbotvtbl(lc)
  end do
  read(19,*)
  read(19,*) topt_data
  read(19,*)
  read(19,*) cmcmax_data
  read(19,*)
  read(19,*) cfactr_data
  read(19,*)
  read(19,*) rsmax_data
  read(19,*)
  read(19,*) bare
  read(19,*)
  read(19,*) natural
  close(19)

  print *, 'VEGPARM read ok: topt=', topt_data, 'bare=', bare, 'shd(2)=', shdtbl(2)

  ! Dump binary reference of tables
  open(unit=30, file='tables_ref.bin', access='stream', form='unformatted', status='replace')
  ! General
  write(30) slpcats
  write(30) slope_data
  write(30) sbeta_data, fxexp_data, csoil_data, salp_data, refdk_data, refkdt_data, &
            frzk_data, zbot_data, czil_data, smlow_data, smhigh_data, lvcoef_data

  ! Soil
  write(30) slcats
  write(30) bb, drysmc, f11, maxsmc, refsmc, satpsi, satdk, satdw, wltsmc, qtz

  ! Veg
  write(30) lucats, bare, natural
  write(30) topt_data, cmcmax_data, cfactr_data, rsmax_data
  write(30) nrotbl
  write(30) snuptbl, rstbl, rgltbl, hstbl, shdtbl, maxalb, emissmintbl, emissmaxtbl, &
            laimintbl, laimaxtbl, z0mintbl, z0maxtbl, albedomintbl, albedomaxtbl, &
            ztopvtbl, zbotvtbl
  close(30)
  print *, 'Dumped tables_ref.bin successfully'
end program test_read_tables
