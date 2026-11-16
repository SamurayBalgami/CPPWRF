// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "phys/ra_rrtmg_common/rrtmg_real.hpp"

#include "core/backend.hpp"
#include <cstdint>
#include <string>

namespace cppwrf::phys {

#pragma pack(push, 1)

struct Band01Tables {
  rr fracrefa[10];
  rr fracrefb[10];
  rr ka[10 * 65];
  rr kb[10 * 235];
  rr ka_mn2[10 * 19];
  rr kb_mn2[10 * 19];
  rr selfref[10 * 10];
  rr forref[10 * 4];
};

struct Band02Tables {
  rr fracrefa[12];
  rr fracrefb[12];
  rr ka[12 * 65];
  rr kb[12 * 235];
  rr selfref[12 * 10];
  rr forref[12 * 4];
};

struct Band03Tables {
  rr fracrefa[9 * 16];
  rr fracrefb[5 * 16];
  rr ka[16 * 585];
  rr kb[16 * 1175];
  rr ka_mn2o[16 * 9 * 19];
  rr kb_mn2o[16 * 5 * 19];
  rr selfref[16 * 10];
  rr forref[16 * 4];
};

struct Band04Tables {
  rr fracrefa[9 * 14];
  rr fracrefb[5 * 14];
  rr ka[14 * 585];
  rr kb[14 * 1175];
  rr selfref[14 * 10];
  rr forref[14 * 4];
};

struct Band05Tables {
  rr fracrefa[9 * 16];
  rr fracrefb[5 * 16];
  rr ka[16 * 585];
  rr kb[16 * 1175];
  rr ka_mo3[16 * 9 * 19];
  rr ccl4[16];
  rr selfref[16 * 10];
  rr forref[16 * 4];
};

struct Band06Tables {
  rr fracrefa[8];
  rr ka[8 * 65];
  rr ka_mco2[8 * 19];
  rr cfc11adj[8];
  rr cfc12[8];
  rr selfref[8 * 10];
  rr forref[8 * 4];
};

struct Band07Tables {
  rr fracrefa[9 * 12];
  rr fracrefb[12];
  rr ka[12 * 585];
  rr kb[12 * 235];
  rr ka_mco2[12 * 9 * 19];
  rr kb_mco2[12 * 19];
  rr selfref[12 * 10];
  rr forref[12 * 4];
};

struct Band08Tables {
  rr fracrefa[8];
  rr fracrefb[8];
  rr ka[8 * 65];
  rr kb[8 * 235];
  rr ka_mco2[8 * 19];
  rr kb_mco2[8 * 19];
  rr ka_mn2o[8 * 19];
  rr kb_mn2o[8 * 19];
  rr ka_mo3[8 * 19];
  rr cfc12[8];
  rr cfc22adj[8];
  rr selfref[8 * 10];
  rr forref[8 * 4];
};

struct Band09Tables {
  rr fracrefa[9 * 12];
  rr fracrefb[12];
  rr ka[12 * 585];
  rr kb[12 * 235];
  rr ka_mn2o[12 * 9 * 19];
  rr kb_mn2o[12 * 19];
  rr selfref[12 * 10];
  rr forref[12 * 4];
};

struct Band10Tables {
  rr fracrefa[6];
  rr fracrefb[6];
  rr ka[6 * 65];
  rr kb[6 * 235];
  rr selfref[6 * 10];
  rr forref[6 * 4];
};

struct Band11Tables {
  rr fracrefa[8];
  rr fracrefb[8];
  rr ka[8 * 65];
  rr kb[8 * 235];
  rr ka_mo2[8 * 19];
  rr kb_mo2[8 * 19];
  rr selfref[8 * 10];
  rr forref[8 * 4];
};

struct Band12Tables {
  rr fracrefa[9 * 8];
  rr ka[8 * 585];
  rr selfref[8 * 10];
  rr forref[8 * 4];
};

struct Band13Tables {
  rr fracrefa[9 * 4];
  rr fracrefb[4];
  rr ka[4 * 585];
  rr ka_mco2[4 * 9 * 19];
  rr ka_mco[4 * 9 * 19];
  rr kb_mo3[4 * 19];
  rr selfref[4 * 10];
  rr forref[4 * 4];
};

struct Band14Tables {
  rr fracrefa[2];
  rr fracrefb[2];
  rr ka[2 * 65];
  rr kb[2 * 235];
  rr selfref[2 * 10];
  rr forref[2 * 4];
};

struct Band15Tables {
  rr fracrefa[9 * 2];
  rr ka[2 * 585];
  rr ka_mn2[2 * 9 * 19];
  rr selfref[2 * 10];
  rr forref[2 * 4];
};

struct Band16Tables {
  rr fracrefa[9 * 2];
  rr fracrefb[2];
  rr ka[2 * 585];
  rr kb[2 * 235];
  rr selfref[2 * 10];
  rr forref[2 * 4];
};

struct RrtmgLwTables {
  // Pade approximation tables (rrlw_tbl)
  rr tau_tbl[10001];
  rr exp_tbl[10001];
  rr tfn_tbl[10001];

  // Spectral / Planck tables (rrlw_wvn)
  rr totplnk[16 * 181];
  rr totplk16[181];
  rr delwave[16];
  rr wavenum1[16];
  rr wavenum2[16];

  int32_t ngc[16];
  int32_t ngs[16];
  int32_t ngn[140];
  int32_t ngb[140];
  int32_t ngm[256];
  rr wt[16];
  rr rwgt[256];

  // Cloud property tables (rrlw_cld)
  rr abscld1;
  rr absice0[2];
  rr absice1[2 * 5];
  rr absice2[43 * 16];
  rr absice3[46 * 16];
  rr absliq0;
  rr absliq1[58 * 16];

  // Reference profiles (rrlw_ref)
  rr preflog[59];
  rr tref[59];
  rr chi_mls[7 * 59];

  // 16 longwave absorption bands
  Band01Tables b01;
  Band02Tables b02;
  Band03Tables b03;
  Band04Tables b04;
  Band05Tables b05;
  Band06Tables b06;
  Band07Tables b07;
  Band08Tables b08;
  Band09Tables b09;
  Band10Tables b10;
  Band11Tables b11;
  Band12Tables b12;
  Band13Tables b13;
  Band14Tables b14;
  Band15Tables b15;
  Band16Tables b16;
};

#pragma pack(pop)


// Loads the tables dump (tests/phys/ra_rrtmg_lw/dump_tables.f90, data/rrtmg/rrtmg_lw_tables.bin)
bool load_rrtmg_lw_tables(const std::string& filepath, RrtmgLwTables& tables);

}  // namespace cppwrf::phys
