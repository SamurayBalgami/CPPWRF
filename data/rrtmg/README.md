# RRTMG tables (WRF v4.8.0)

`rrtmg_lw_tables.bin` / `rrtmg_sw_tables.bin` are the RRTMG longwave / shortwave tables exactly as
WRF v4.8.0 (real*8 build) holds them after `rrtmg_lwinit` / `rrtmg_swinit` read
`RRTMG_LW_DATA_DBL` / `RRTMG_SW_DATA_DBL` and reduce them to 140 / 112 g-points (cmbgb*), plus the
cloud property data and the Pade tables: raw memory images of `RrtmgLwTables` /
`RrtmgSwTables` (src/phys/ra_rrtmg_*/rrtmg_*_tables.hpp; sizes are checked by static_assert),
dumped by an instrumented WRF build (tests/phys/ra_rrtmg_sw/dump_sw_tables.f90 for the SW one).
With them the C++ RRTMG is bit-identical to WRF (docs/PHYSICS_SUITE.md sections 8-10).

The radiation driver looks for them in `phys.rrtmg_tables` (default: the WRF run directory of the
case) and then here (`CPPWRF_DATA_DIR`, set at build time). The ozone climatology
(`ozone*.formatted`) and `CAMtr_volume_mixing_ratio` are read from `phys.rrtmg_tables` as WRF reads
them from its run directory.
