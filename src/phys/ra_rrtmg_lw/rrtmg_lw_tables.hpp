#pragma once

#include "core/backend.hpp"
#include <cstdint>
#include <string>

namespace cppwrf::phys {

#pragma pack(push, 1)

struct Band01Tables {
  double fracrefa[10];
  double fracrefb[10];
  double ka[10 * 65];
  double kb[10 * 235];
  double ka_mn2[10 * 19];
  double kb_mn2[10 * 19];
  double selfref[10 * 10];
  double forref[10 * 4];
};

struct Band02Tables {
  double fracrefa[12];
  double fracrefb[12];
  double ka[12 * 65];
  double kb[12 * 235];
  double selfref[12 * 10];
  double forref[12 * 4];
};

struct Band03Tables {
  double fracrefa[9 * 16];
  double fracrefb[5 * 16];
  double ka[16 * 585];
  double kb[16 * 1175];
  double ka_mn2o[16 * 9 * 19];
  double kb_mn2o[16 * 5 * 19];
  double selfref[16 * 10];
  double forref[16 * 4];
};

struct Band04Tables {
  double fracrefa[9 * 14];
  double fracrefb[5 * 14];
  double ka[14 * 585];
  double kb[14 * 1175];
  double selfref[14 * 10];
  double forref[14 * 4];
};

struct Band05Tables {
  double fracrefa[9 * 16];
  double fracrefb[5 * 16];
  double ka[16 * 585];
  double kb[16 * 1175];
  double ka_mo3[16 * 9 * 19];
  double ccl4[16];
  double selfref[16 * 10];
  double forref[16 * 4];
};

struct Band06Tables {
  double fracrefa[8];
  double ka[8 * 65];
  double ka_mco2[8 * 19];
  double cfc11adj[8];
  double cfc12[8];
  double selfref[8 * 10];
  double forref[8 * 4];
};

struct Band07Tables {
  double fracrefa[9 * 12];
  double fracrefb[12];
  double ka[12 * 585];
  double kb[12 * 235];
  double ka_mco2[12 * 9 * 19];
  double kb_mco2[12 * 19];
  double selfref[12 * 10];
  double forref[12 * 4];
};

struct Band08Tables {
  double fracrefa[8];
  double fracrefb[8];
  double ka[8 * 65];
  double kb[8 * 235];
  double ka_mco2[8 * 19];
  double kb_mco2[8 * 19];
  double ka_mn2o[8 * 19];
  double kb_mn2o[8 * 19];
  double ka_mo3[8 * 19];
  double cfc12[8];
  double cfc22adj[8];
  double selfref[8 * 10];
  double forref[8 * 4];
};

struct Band09Tables {
  double fracrefa[9 * 12];
  double fracrefb[12];
  double ka[12 * 585];
  double kb[12 * 235];
  double ka_mn2o[12 * 9 * 19];
  double kb_mn2o[12 * 19];
  double selfref[12 * 10];
  double forref[12 * 4];
};

struct Band10Tables {
  double fracrefa[6];
  double fracrefb[6];
  double ka[6 * 65];
  double kb[6 * 235];
  double selfref[6 * 10];
  double forref[6 * 4];
};

struct Band11Tables {
  double fracrefa[8];
  double fracrefb[8];
  double ka[8 * 65];
  double kb[8 * 235];
  double ka_mo2[8 * 19];
  double kb_mo2[8 * 19];
  double selfref[8 * 10];
  double forref[8 * 4];
};

struct Band12Tables {
  double fracrefa[9 * 8];
  double ka[8 * 585];
  double selfref[8 * 10];
  double forref[8 * 4];
};

struct Band13Tables {
  double fracrefa[9 * 4];
  double fracrefb[4];
  double ka[4 * 585];
  double ka_mco2[4 * 9 * 19];
  double ka_mco[4 * 9 * 19];
  double kb_mo3[4 * 19];
  double selfref[4 * 10];
  double forref[4 * 4];
};

struct Band14Tables {
  double fracrefa[2];
  double fracrefb[2];
  double ka[2 * 65];
  double kb[2 * 235];
  double selfref[2 * 10];
  double forref[2 * 4];
};

struct Band15Tables {
  double fracrefa[9 * 2];
  double ka[2 * 585];
  double ka_mn2[2 * 9 * 19];
  double selfref[2 * 10];
  double forref[2 * 4];
};

struct Band16Tables {
  double fracrefa[9 * 2];
  double fracrefb[2];
  double ka[2 * 585];
  double kb[2 * 235];
  double selfref[2 * 10];
  double forref[2 * 4];
};

struct RrtmgLwTables {
  // Pade approximation tables (rrlw_tbl)
  double tau_tbl[10001];
  double exp_tbl[10001];
  double tfn_tbl[10001];

  // Spectral / Planck tables (rrlw_wvn)
  double totplnk[16 * 181];
  double totplk16[181];
  double delwave[16];
  double wavenum1[16];
  double wavenum2[16];

  int32_t ngc[16];
  int32_t ngs[16];
  int32_t ngn[140];
  int32_t ngb[140];
  int32_t ngm[256];
  double wt[16];
  double rwgt[256];

  // Cloud property tables (rrlw_cld)
  double abscld1;
  double absice0[2];
  double absice1[2 * 5];
  double absice2[43 * 16];
  double absice3[46 * 16];
  double absliq0;
  double absliq1[58 * 16];

  // Reference profiles (rrlw_ref)
  double preflog[59];
  double tref[59];
  double chi_mls[7 * 59];

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

static_assert(sizeof(RrtmgLwTables) == 1435576, "RrtmgLwTables size must be exactly 1435576 bytes");

// Loads the tables dump (tests/phys/ra_rrtmg_lw/dump_tables.f90, data/rrtmg/rrtmg_lw_tables.bin)
bool load_rrtmg_lw_tables(const std::string& filepath, RrtmgLwTables& tables);

}  // namespace cppwrf::phys
