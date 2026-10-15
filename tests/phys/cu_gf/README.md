# Grell-Freitas (cu_physics = 3) replay against WRF's GFDRV

1. Dump GFDRV's inputs and CPPWRF's outputs at step n (any CPU build):

       CPPWRF_GF_DUMP=gf.bin CPPWRF_GF_DUMP_STEP=n build/cppwrf configs/real_p8_r8.ini --set run.run_time=$((60*n))

2. Run WRF's own GFDRV (r8d8 objects) on them and compare RTHCUTEN, RQVCUTEN, RUCUTEN, RAINCV:

       W=~/wrf_real/WRF-4.6.1-r8d8/phys
       gfortran -O2 -fdefault-real-8 -fdefault-double-8 -I$W -o harness harness.f90 $W/module_cu_gf_wrfdrv.o \
         $W/module_cu_gf_deep.o $W/module_cu_gf_sh.o $W/module_gfs_physcons.o $W/module_gfs_machine.o
       ulimit -s unlimited; ./harness gf.bin [ishallow=1] [ichoice=0] [dt=60] [dx=15000]

   P8 steps 1, 5, 10 (deep and shallow+deep): every column identical (max relative difference 0).

Notes found on the way (all reproduced in src/phys/cu_gf):
- gfortran vectorises the zu/zd profile loops of get_zu_zd_pdf_fim and calls glibc's 2-lane vector pow
  `_ZGVbN2vv_pow` for pairs of iterations (scalar pow for an odd last one); it can differ from the scalar
  pow in the last bit. The C++ calls the same libmvec function on the host (gf_common.hpp, zu_profile).
- get_inversion_layers (shallow) keeps its scratch sec_deriv across the columns of a row, so the inversion
  chosen for a column can depend on the column before it: CPPWRF runs a per-row scan for it.
- That routine also reads sec_deriv(2), which nothing sets, when kbcon = 2; CPPWRF takes it as
  phys.gf_sec2 (default 1e9).
- GFDRV's constants are module_gfs_physcons (real*8): g = 9.80665, cp = 1004.6; the scheme's own
  are g = 9.81, cp = 1004, r_v = 461.
- phy_cu_ten does not add RUCUTEN/RVCUTEN for GF (momentum transport is computed but unused).
