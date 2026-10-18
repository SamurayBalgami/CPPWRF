# RRTMG LW / SW checks

The RRTMG ports are verified against WRF in two ways:

1. **Whole model runs** (the main check): P3 (RRTMG every step) and P4 (cloudy RRTMG with WSM6
   radii) give wrfout and wrfrst files identical to WRF v4.8.0 real*8. This covers every radiation
   field (RTHRATLW/SW, GLW, GSW, SWDOWN, the LW/SW flux diagnostics and their accumulations, CLDFRA,
   O3RAD). Run `tests/verify/run_wrf_files_test.sh <cpu fp64 build> p3 p4`.
2. **Column checks against WRF dumps** (this directory and `../ra_rrtmg_sw`). A debug build of WRF
   writes the RRTMG_LWRAD / RRTMG_SWRAD inputs and outputs (whole memory arrays, big-endian) at a
   radiation call. The checks run the C++ code on the same inputs and compare every output:
   - `real_wrapper_check <in.bin> <out.bin> <tables.bin> <CAMtr_volume_mixing_ratio>`: the CPU column
     path (`rrtmg_lw_column` / `rrtmg_sw_column`), bit-identical on all columns;
   - `real_core_check <core.bin> <tables.bin>`: the rrtmg_lw / rrtmg_sw cores alone;
   - `batch_check <in.bin> <out.bin> <tables.bin> <CAMtr...> [chunk]`: the batch (GPU) path. It is
     bit-identical on the CPU and agrees to roundoff on the GPU (docs/PHYSICS_SUITE.md section 10).

`dump_tables.f90` (LW) and `../ra_rrtmg_sw/dump_sw_tables.f90` (SW) write the tables in
`data/rrtmg` from WRF's `rrtmg_lwinit` / `rrtmg_swinit` with the `RRTMG_*_DATA_DBL` files. Build them
against WRF's `phys/module_ra_rrtmg_lw.F` / `module_ra_rrtmg_sw.F` with WRF's real*8 flags. The
RRTMG data files and the table initialisation are the same in WRF v4.6.1 and v4.8.0.
