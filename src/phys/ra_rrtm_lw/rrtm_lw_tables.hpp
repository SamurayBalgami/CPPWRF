// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// WRF 4.8.0 phys/module_ra_rrtm.F (RRTM longwave, ra_lw_physics = 1): the module's SAVE data the column code reads,
// as rrtminit / rrtm_lookuptable leave them: the Pade tables, the absorption coefficients of RRTM_DATA reduced from
// 256 to 140 g-points (CMBGB1..16), and the DATA tables (rrtm_lw_data.hpp) the column routines use. Every array has
// the Fortran element order (first index fastest) and is indexed with Fortran's bounds in rrtm_lw.hpp; a plain
// struct, copied to the device as it is.
#include <string>

namespace cppwrf::phys {

struct RrtmLwTables {
  // rrtm_lookuptable: TAU, TF, TRANS (0:5000), CORR1, CORR2 (0:200), BPADE; rrtminit: FLUXFAC; DATA: HEATFAC
  double tau[5001], tf[5001], trans[5001], corr1[201], corr2[201];
  double bpade, fluxfac, heatfac;
  int nspa[16], nspb[16], ngb[140], ixindx[35];
  double preflog[59], tref[59], totplnk[181 * 16], delwave[16];
  // band 1: ABSA1(65,NG1), ABSB1(235,NG1), SELFREFC1(10,NG1), FORREFC1(NG1), FRACREFAC1(NG1), FRACREFBC1(NG1)
  double absa1[65 * 8], absb1[235 * 8], selfrefc1[10 * 8], forrefc1[8], fracrefac1[8], fracrefbc1[8];
  double absa2[65 * 14], absb2[235 * 14], selfrefc2[10 * 14], forrefc2[14], fracrefac2[14 * 13], fracrefbc2[14];
  double absa3[650 * 16], absb3[1175 * 16], selfrefc3[10 * 16], forrefc3[16], absn2oac3[16], absn2obc3[16],
      fracrefac3[16 * 10], fracrefbc3[16 * 5];
  double absa4[585 * 14], absb4[1410 * 14], selfrefc4[10 * 14], fracrefac4[14 * 9], fracrefbc4[14 * 6];
  double absa5[585 * 16], absb5[1175 * 16], selfrefc5[10 * 16], ccl4c5[16], fracrefac5[16 * 9], fracrefbc5[16 * 5];
  double absa6[65 * 8], selfrefc6[10 * 8], absco2c6[8], cfc11adjc6[8], cfc12c6[8], fracrefac6[8];
  double absa7[585 * 12], absb7[235 * 12], selfrefc7[10 * 12], absco2c7[12], fracrefac7[12 * 9], fracrefbc7[12];
  double absa8[35 * 8], absb8[265 * 8], selfrefc8[10 * 8], absco2ac8[8], absco2bc8[8], absn2oac8[8],
      absn2obc8[8], cfc12c8[8], cfc22adjc8[8], fracrefac8[8], fracrefbc8[8];
  double absa9[715 * 12], absb9[235 * 12], selfrefc9[10 * 12], absn2oc9[3 * 12], fracrefac9[12 * 9],
      fracrefbc9[12];
  double absa10[65 * 6], absb10[235 * 6], fracrefac10[6], fracrefbc10[6];
  double absa11[65 * 8], absb11[235 * 8], selfrefc11[10 * 8], fracrefac11[8], fracrefbc11[8];
  double absa12[585 * 8], selfrefc12[10 * 8], fracrefac12[8 * 9];
  double absa13[585 * 4], selfrefc13[10 * 4], fracrefac13[4 * 9];
  double absa14[65 * 2], absb14[235 * 2], selfrefc14[10 * 2], fracrefac14[2], fracrefbc14[2];
  double absa15[585 * 2], selfrefc15[10 * 2], fracrefac15[2 * 9];
  double absa16[585 * 2], selfrefc16[10 * 2], fracrefac16[2 * 9];
  // the subroutines' DATA: TAUGB2 REFPARAM, TAUGB3 / TAUGB8 / TAUGB9 reference profiles, O3DATA's climatology,
  // MM5ATM's buffer-layer temperature profile
  double refparam2[13], etaref3[10], h2oref3[59], n2oref3[59], co2ref3[59], h2oref8[59], n2oref8[59], o3ref8[59],
      n2oref9[13], h2oref9[13], ch4ref9[13], etaref9[11];
  double o3sum[31], ppsum[31], o3win[31], ppwin[31], pprof[60], tprof[60];
};

// rrtm_lookuptable with RRTM_DATA at path (WRF's run directory file, RRTM_DATA_DBL in the real*8 build: big-endian,
// 4-byte record markers); false with the message WRF stops with when the file cannot be opened or a record is short
bool load_rrtm_lw_tables(const std::string& path, RrtmLwTables& t, std::string& error);

}  // namespace cppwrf::phys
