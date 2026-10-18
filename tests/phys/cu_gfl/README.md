# Grell-Freitas-Li (`cu_physics = 3`, `GFLSCHEME`) replay against WRF v4.8.0 `cu_grell_freitas_li`

1. Dump `cu_grell_freitas_li`'s inputs and CPPWRF's outputs at step `n` (any CPU build):

       CPPWRF_GFL_DUMP=gfl.bin CPPWRF_GFL_DUMP_STEP=n build/cppwrf configs/real_p8_r8.ini --set run.run_time=$((60*n))

2. Run WRF v4.8.0's own `cu_grell_freitas_li` (`r8d8` objects) on them and compare `RTHCUTEN`, `RQVCUTEN`, `RQCCUTEN`, `RQICUTEN`, `RUCUTEN`, `RVCUTEN`, `RAINCV`, `PRATEC`, `HTOP`, and `HBOT`:

       W=$WRF_DIR/phys
       gfortran -O2 -fdefault-real-8 -fdefault-double-8 -I$W -o harness harness.f90 \
         $W/module_cu_gfl_wrf.o $W/module_cu_gfl_deep.o $W/module_cu_gfl_sh.o \
         $W/module_cu_gfl_common.o $W/ccpp_kind_types.o
       ulimit -s unlimited; ./harness gfl.bin [ishallow=1] [ichoice=0] [dt=60] [dx=15000] [itimestep=1]

   P8 steps 1, 5, 10 (deep and shallow+deep): every column identical (max absolute difference 0).

Notes on WRF v4.8.0 GFL (`src/phys/cu_gfl`):
- `gfortran` vectorises the `zu` profile loop of `get_zu_zd_pdf_fim` (`module_cu_gfl_deep.o`) and calls glibc's 2-lane vector pow `_ZGVbN2vv_pow` for pairs of iterations (scalar `pow` for an odd remainder iteration), evaluating `(((x**(a - 1.) * b) / beta) * (1. - x)**(b - 1.))`; it can differ from scalar `pow` or left-to-right multiplication by 1 ULP. The C++ calls the same `libmvec` function on the host (`gfl_common.hpp`, `zu_profile`) and preserves Fortran's evaluation order (`rates_up_pdf` `smooth_Heaviside` and `irainevap`'s `edt * evfact * (sig_init * sig_init)`).
- Unlike WRF v4.6.1 GF, `get_inversion_layers` in GFL v4.8.0 initializes `sec_deriv(:) = 0.` at entry, starts at `ilev = max(4, kstart(i)+1)`, and resets `sec_deriv(:) = 1.e9` inside the column loop `do i = its, itf`, so there is no row-carried state or uninitialized read of `sec_deriv(2)`.
- `cu_grell_freitas_li` receives WRF's model constants (`xlv = 2.5e6`, `cp = 1004.5`, `g = 9.81`, `r_v = 461.6`) as `real(kind=kind_phys)` for `dhdt`, `omeg`, and `mconv`, while `module_cu_gfl_deep` and `module_cu_gfl_sh` use their own module parameters (`g = 9.81`, `cp = 1004.`, `xlv = 2.5e6`, `r_v = 461.`, `tcrit = 258.`, `fluxtune = 1.5`).
- In `cu_grell_freitas_li`, shallow convection (`cu_gfl_sh_run`, `ichoice = 3`) runs after `cu_gfl_deep_run` + `neg_check` for columns where `deep.ierr == 0`.
- `phy_cu_ten` (`module_physics_addtendc.F`, `CASE (GDSCHEME, G3SCHEME, GFLSCHEME)`) still only adds `RTHCUTEN`, `RQVCUTEN`, `RQCCUTEN`, and `RQICUTEN` (momentum tendencies `RUCUTEN`/`RVCUTEN` are computed by GFL and used internally for ECMWF-style kinetic energy dissipation heating, but are not added to the model momentum equations).
