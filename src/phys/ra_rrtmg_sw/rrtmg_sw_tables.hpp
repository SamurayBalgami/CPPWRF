#pragma once

#include "core/backend.hpp"
#include <cstdint>
#include <string>

namespace cppwrf::phys {

#pragma pack(push, 1)

struct Band16Tables {
  double ka[9 * 5 * 13 * 6];       // (9,5,13,6)
  double kb[5 * 47 * 6];           // (5,13:59,6)
  double selfref[10 * 6];          // (10,6)
  double forref[3 * 6];            // (3,6)
  double sfluxref[6];              // (6)
  double rayl;
  double strrat1;
  int32_t layreffr;
};

struct Band17Tables {
  double ka[9 * 5 * 13 * 12];      // (9,5,13,12)
  double kb[5 * 5 * 47 * 12];      // (5,5,13:59,12)
  double selfref[10 * 12];         // (10,12)
  double forref[4 * 12];           // (4,12)
  double sfluxref[12 * 5];         // (12,5)
  double rayl;
  double strrat;
  int32_t layreffr;
};

struct Band18Tables {
  double ka[9 * 5 * 13 * 8];       // (9,5,13,8)
  double kb[5 * 47 * 8];           // (5,13:59,8)
  double selfref[10 * 8];          // (10,8)
  double forref[3 * 8];            // (3,8)
  double sfluxref[8 * 9];          // (8,9)
  double rayl;
  double strrat;
  int32_t layreffr;
};

struct Band19Tables {
  double ka[9 * 5 * 13 * 8];       // (9,5,13,8)
  double kb[5 * 47 * 8];           // (5,13:59,8)
  double selfref[10 * 8];          // (10,8)
  double forref[3 * 8];            // (3,8)
  double sfluxref[8 * 9];          // (8,9)
  double rayl;
  double strrat;
  int32_t layreffr;
};

struct Band20Tables {
  double ka[5 * 13 * 10];          // (5,13,10)
  double kb[5 * 47 * 10];          // (5,13:59,10)
  double selfref[10 * 10];         // (10,10)
  double forref[4 * 10];           // (4,10)
  double sfluxref[10];             // (10)
  double absch4[10];               // (10)
  double rayl;
  int32_t layreffr;
};

struct Band21Tables {
  double ka[9 * 5 * 13 * 10];      // (9,5,13,10)
  double kb[5 * 5 * 47 * 10];      // (5,5,13:59,10)
  double selfref[10 * 10];         // (10,10)
  double forref[4 * 10];           // (4,10)
  double sfluxref[10 * 9];         // (10,9)
  double rayl;
  double strrat;
  int32_t layreffr;
};

struct Band22Tables {
  double ka[9 * 5 * 13 * 2];       // (9,5,13,2)
  double kb[5 * 47 * 2];           // (5,13:59,2)
  double selfref[10 * 2];          // (10,2)
  double forref[3 * 2];            // (3,2)
  double sfluxref[2 * 9];          // (2,9)
  double rayl;
  double strrat;
  int32_t layreffr;
};

struct Band23Tables {
  double ka[5 * 13 * 10];          // (5,13,10)
  double selfref[10 * 10];         // (10,10)
  double forref[3 * 10];           // (3,10)
  double sfluxref[10];             // (10)
  double rayl[10];                 // (10)
  double givfac;
  int32_t layreffr;
};

struct Band24Tables {
  double ka[9 * 5 * 13 * 8];       // (9,5,13,8)
  double kb[5 * 47 * 8];           // (5,13:59,8)
  double selfref[10 * 8];          // (10,8)
  double forref[3 * 8];            // (3,8)
  double sfluxref[8 * 9];          // (8,9)
  double abso3a[8];                // (8)
  double abso3b[8];                // (8)
  double rayla[8 * 9];             // (8,9)
  double raylb[8];                 // (8)
  double strrat;
  int32_t layreffr;
};

struct Band25Tables {
  double ka[5 * 13 * 6];           // (5,13,6)
  double sfluxref[6];              // (6)
  double abso3a[6];                // (6)
  double abso3b[6];                // (6)
  double rayl[6];                  // (6)
  int32_t layreffr;
};

struct Band26Tables {
  double sfluxref[6];              // (6)
  double rayl[6];                  // (6)
};

struct Band27Tables {
  double ka[5 * 13 * 8];           // (5,13,8)
  double kb[5 * 47 * 8];           // (5,13:59,8)
  double sfluxref[8];              // (8)
  double rayl[8];                  // (8)
  double scalekur;
  int32_t layreffr;
};

struct Band28Tables {
  double ka[9 * 5 * 13 * 6];       // (9,5,13,6)
  double kb[5 * 5 * 47 * 6];       // (5,5,13:59,6)
  double sfluxref[6 * 5];          // (6,5)
  double rayl;
  double strrat;
  int32_t layreffr;
};

struct Band29Tables {
  double ka[5 * 13 * 12];          // (5,13,12)
  double kb[5 * 47 * 12];          // (5,13:59,12)
  double selfref[10 * 12];         // (10,12)
  double forref[4 * 12];           // (4,12)
  double sfluxref[12];             // (12)
  double absh2o[12];               // (12)
  double absco2[12];               // (12)
  double rayl;
  int32_t layreffr;
};

struct RrtmgSwTables {
  // 1. Pade approximation lookup table (rrsw_tbl)
  double exp_tbl[10001];

  // 2. Spectral info (rrsw_wvn)
  double wavenum1[14];
  double wavenum2[14];
  double delwave[14];
  int32_t ngc[14];
  int32_t ngs[14];
  int32_t ngn[112];
  int32_t ngb[112];
  int32_t ngm[224];
  double wt[16];
  double rwgt[224];

  // 3. Cloud property tables (rrsw_cld)
  double extliq1[58 * 14];
  double ssaliq1[58 * 14];
  double asyliq1[58 * 14];
  double extice2[43 * 14];
  double ssaice2[43 * 14];
  double asyice2[43 * 14];
  double extice3[46 * 14];
  double ssaice3[46 * 14];
  double asyice3[46 * 14];
  double fdlice3[46 * 14];
  double abari[5];
  double bbari[5];
  double cbari[5];
  double dbari[5];
  double ebari[5];
  double fbari[5];

  // 4. Reference profiles (rrsw_ref)
  double preflog[59];
  double tref[59];

  // 5. Aerosol properties (rrsw_aer)
  double rsrtaua[14 * 6];
  double rsrpiza[14 * 6];
  double rsrasya[14 * 6];

  // 6. 14 shortwave absorption bands
  Band16Tables b16;
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

// Loads precomputed tables from binary dump file
bool load_rrtmg_sw_tables(const std::string& filepath, RrtmgSwTables& tables);

// Global default tables instance (host)
RrtmgSwTables& get_default_rrtmg_sw_tables();

}  // namespace cppwrf::phys
