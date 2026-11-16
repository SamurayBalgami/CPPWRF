// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "phys/ra_rrtmg_common/rrtmg_real.hpp"

#include "core/backend.hpp"
#include <cstdint>
#include <string>

namespace cppwrf::phys {

#pragma pack(push, 1)

struct SwBand16Tables {
  rr ka[9 * 5 * 13 * 6];       // (9,5,13,6)
  rr kb[5 * 47 * 6];           // (5,13:59,6)
  rr selfref[10 * 6];          // (10,6)
  rr forref[3 * 6];            // (3,6)
  rr sfluxref[6];              // (6)
  rr rayl;
  rr strrat1;
  int32_t layreffr;
};

struct Band17Tables {
  rr ka[9 * 5 * 13 * 12];      // (9,5,13,12)
  rr kb[5 * 5 * 47 * 12];      // (5,5,13:59,12)
  rr selfref[10 * 12];         // (10,12)
  rr forref[4 * 12];           // (4,12)
  rr sfluxref[12 * 5];         // (12,5)
  rr rayl;
  rr strrat;
  int32_t layreffr;
};

struct Band18Tables {
  rr ka[9 * 5 * 13 * 8];       // (9,5,13,8)
  rr kb[5 * 47 * 8];           // (5,13:59,8)
  rr selfref[10 * 8];          // (10,8)
  rr forref[3 * 8];            // (3,8)
  rr sfluxref[8 * 9];          // (8,9)
  rr rayl;
  rr strrat;
  int32_t layreffr;
};

struct Band19Tables {
  rr ka[9 * 5 * 13 * 8];       // (9,5,13,8)
  rr kb[5 * 47 * 8];           // (5,13:59,8)
  rr selfref[10 * 8];          // (10,8)
  rr forref[3 * 8];            // (3,8)
  rr sfluxref[8 * 9];          // (8,9)
  rr rayl;
  rr strrat;
  int32_t layreffr;
};

struct Band20Tables {
  rr ka[5 * 13 * 10];          // (5,13,10)
  rr kb[5 * 47 * 10];          // (5,13:59,10)
  rr selfref[10 * 10];         // (10,10)
  rr forref[4 * 10];           // (4,10)
  rr sfluxref[10];             // (10)
  rr absch4[10];               // (10)
  rr rayl;
  int32_t layreffr;
};

struct Band21Tables {
  rr ka[9 * 5 * 13 * 10];      // (9,5,13,10)
  rr kb[5 * 5 * 47 * 10];      // (5,5,13:59,10)
  rr selfref[10 * 10];         // (10,10)
  rr forref[4 * 10];           // (4,10)
  rr sfluxref[10 * 9];         // (10,9)
  rr rayl;
  rr strrat;
  int32_t layreffr;
};

struct Band22Tables {
  rr ka[9 * 5 * 13 * 2];       // (9,5,13,2)
  rr kb[5 * 47 * 2];           // (5,13:59,2)
  rr selfref[10 * 2];          // (10,2)
  rr forref[3 * 2];            // (3,2)
  rr sfluxref[2 * 9];          // (2,9)
  rr rayl;
  rr strrat;
  int32_t layreffr;
};

struct Band23Tables {
  rr ka[5 * 13 * 10];          // (5,13,10)
  rr selfref[10 * 10];         // (10,10)
  rr forref[3 * 10];           // (3,10)
  rr sfluxref[10];             // (10)
  rr rayl[10];                 // (10)
  rr givfac;
  int32_t layreffr;
};

struct Band24Tables {
  rr ka[9 * 5 * 13 * 8];       // (9,5,13,8)
  rr kb[5 * 47 * 8];           // (5,13:59,8)
  rr selfref[10 * 8];          // (10,8)
  rr forref[3 * 8];            // (3,8)
  rr sfluxref[8 * 9];          // (8,9)
  rr abso3a[8];                // (8)
  rr abso3b[8];                // (8)
  rr rayla[8 * 9];             // (8,9)
  rr raylb[8];                 // (8)
  rr strrat;
  int32_t layreffr;
};

struct Band25Tables {
  rr ka[5 * 13 * 6];           // (5,13,6)
  rr sfluxref[6];              // (6)
  rr abso3a[6];                // (6)
  rr abso3b[6];                // (6)
  rr rayl[6];                  // (6)
  int32_t layreffr;
};

struct Band26Tables {
  rr sfluxref[6];              // (6)
  rr rayl[6];                  // (6)
};

struct Band27Tables {
  rr ka[5 * 13 * 8];           // (5,13,8)
  rr kb[5 * 47 * 8];           // (5,13:59,8)
  rr sfluxref[8];              // (8)
  rr rayl[8];                  // (8)
  rr scalekur;
  int32_t layreffr;
};

struct Band28Tables {
  rr ka[9 * 5 * 13 * 6];       // (9,5,13,6)
  rr kb[5 * 5 * 47 * 6];       // (5,5,13:59,6)
  rr sfluxref[6 * 5];          // (6,5)
  rr rayl;
  rr strrat;
  int32_t layreffr;
};

struct Band29Tables {
  rr ka[5 * 13 * 12];          // (5,13,12)
  rr kb[5 * 47 * 12];          // (5,13:59,12)
  rr selfref[10 * 12];         // (10,12)
  rr forref[4 * 12];           // (4,12)
  rr sfluxref[12];             // (12)
  rr absh2o[12];               // (12)
  rr absco2[12];               // (12)
  rr rayl;
  int32_t layreffr;
};

struct RrtmgSwTables {
  // 1. Pade approximation lookup table (rrsw_tbl)
  rr exp_tbl[10001];

  // 2. Spectral info (rrsw_wvn)
  rr wavenum1[14];
  rr wavenum2[14];
  rr delwave[14];
  int32_t ngc[14];
  int32_t ngs[14];
  int32_t ngn[112];
  int32_t ngb[112];
  int32_t ngm[224];
  rr wt[16];
  rr rwgt[224];

  // 3. Cloud property tables (rrsw_cld)
  rr extliq1[58 * 14];
  rr ssaliq1[58 * 14];
  rr asyliq1[58 * 14];
  rr extice2[43 * 14];
  rr ssaice2[43 * 14];
  rr asyice2[43 * 14];
  rr extice3[46 * 14];
  rr ssaice3[46 * 14];
  rr asyice3[46 * 14];
  rr fdlice3[46 * 14];
  rr abari[5];
  rr bbari[5];
  rr cbari[5];
  rr dbari[5];
  rr ebari[5];
  rr fbari[5];

  // 4. Reference profiles (rrsw_ref)
  rr preflog[59];
  rr tref[59];

  // 5. Aerosol properties (rrsw_aer)
  rr rsrtaua[14 * 6];
  rr rsrpiza[14 * 6];
  rr rsrasya[14 * 6];

  // 6. 14 shortwave absorption bands
  SwBand16Tables b16;
  Band17Tables b17;
  Band18Tables b18;
  Band19Tables b19;
  Band20Tables b20;
  Band21Tables b21;
  Band22Tables b22;
  Band23Tables b23;
  Band24Tables b24;
  Band25Tables b25;
  Band26Tables b26;
  Band27Tables b27;
  Band28Tables b28;
  Band29Tables b29;
};

#pragma pack(pop)


// Loads the tables dump (tests/phys/ra_rrtmg_sw/dump_sw_tables.f90, data/rrtmg/rrtmg_sw_tables.bin)
bool load_rrtmg_sw_tables(const std::string& filepath, RrtmgSwTables& tables);

}  // namespace cppwrf::phys
