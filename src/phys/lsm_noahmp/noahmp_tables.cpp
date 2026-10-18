#include "phys/lsm_noahmp/noahmp_tables.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cppwrf::phys::noahmp {

void validate_noahmp_options(const NoahmpOptions& opt) {
  auto check = [](const char* name, int got, int expected) {
    if (got != expected) {
      throw std::runtime_error(std::string("Noah-MP: unsupported ") + name + "=" +
                               std::to_string(got) + " (only " + std::to_string(expected) +
                               " is ported)");
    }
  };
  check("dveg", opt.dveg, 4);
  check("opt_crs", opt.opt_crs, 1);
  check("opt_btr", opt.opt_btr, 1);
  check("opt_runsrf", opt.opt_runsrf, 3);
  check("opt_runsub", opt.opt_runsub, 3);
  check("opt_sfc", opt.opt_sfc, 1);
  check("opt_frz", opt.opt_frz, 1);
  check("opt_inf", opt.opt_inf, 1);
  check("opt_rad", opt.opt_rad, 3);
  check("opt_alb", opt.opt_alb, 2);
  check("opt_snf", opt.opt_snf, 1);
  check("opt_tbot", opt.opt_tbot, 2);
  check("opt_stc", opt.opt_stc, 1);
  check("opt_gla", opt.opt_gla, 1);
  check("opt_rsf", opt.opt_rsf, 1);
  check("opt_soil", opt.opt_soil, 1);
  check("opt_pedo", opt.opt_pedo, 1);
  check("opt_crop", opt.opt_crop, 0);
  check("opt_irr", opt.opt_irr, 0);
  check("opt_irrm", opt.opt_irrm, 0);
  check("opt_infdv", opt.opt_infdv, 0);
  check("opt_tdrn", opt.opt_tdrn, 0);
  check("opt_tksno", opt.opt_tksno, 1);
  check("opt_compact", opt.opt_compact, 2);
  check("opt_scf", opt.opt_scf, 1);
  check("opt_wetland", opt.opt_wetland, 0);
  check("sf_urban_physics", opt.sf_urban_physics, 0);
  if (opt.soiltstep != 0.0) {
    throw std::runtime_error("Noah-MP: unsupported soiltstep=" + std::to_string(opt.soiltstep) +
                             " (only 0.0 is ported)");
  }
}

namespace {

std::string trim(const std::string& s) {
  std::size_t a = 0;
  while (a < s.size() && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
  std::size_t b = s.size();
  while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
  return s.substr(a, b - a);
}

std::string to_upper(std::string s) {
  for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

std::string to_lower(std::string s) {
  for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

std::string strip_quotes(const std::string& s) {
  std::string t = trim(s);
  if (t.size() >= 2 && ((t.front() == '"' && t.back() == '"') || (t.front() == '\'' && t.back() == '\''))) {
    return trim(t.substr(1, t.size() - 2));
  }
  return t;
}

real parse_real(const std::string& s) {
  if constexpr (sizeof(real) == 8) {
    return static_cast<real>(std::strtod(s.c_str(), nullptr));
  } else {
    return static_cast<real>(std::strtof(s.c_str(), nullptr));
  }
}

int parse_int(const std::string& s) {
  return static_cast<int>(std::strtol(s.c_str(), nullptr, 10));
}

using SectionMap = std::map<std::string, std::vector<std::string>>;
using TableData = std::map<std::string, SectionMap>;

bool parse_table_file(const std::string& path, TableData& out) {
  std::ifstream in(path);
  if (!in.is_open()) return false;

  std::string cur_sec;
  std::string cur_key;
  std::string line;
  while (std::getline(in, line)) {
    // Strip comment outside quotes
    bool in_q = false;
    char q_char = 0;
    std::size_t comment_pos = std::string::npos;
    for (std::size_t i = 0; i < line.size(); ++i) {
      char c = line[i];
      if (in_q) {
        if (c == q_char) in_q = false;
      } else {
        if (c == '"' || c == '\'') {
          in_q = true;
          q_char = c;
        } else if (c == '!') {
          comment_pos = i;
          break;
        }
      }
    }
    if (comment_pos != std::string::npos) line = line.substr(0, comment_pos);
    line = trim(line);
    if (line.empty()) continue;

    if (line[0] == '&') {
      cur_sec = to_lower(trim(line.substr(1)));
      cur_key.clear();
      continue;
    }
    if (line == "/") {
      cur_sec.clear();
      cur_key.clear();
      continue;
    }
    if (cur_sec.empty()) continue;

    std::size_t eq = line.find('=');
    std::string rhs;
    if (eq != std::string::npos) {
      cur_key = to_upper(trim(line.substr(0, eq)));
      rhs = line.substr(eq + 1);
    } else {
      rhs = line;
    }
    if (cur_key.empty()) continue;

    // Split rhs by comma or whitespace (respecting quotes)
    auto& vec = out[cur_sec][cur_key];
    std::string tok;
    in_q = false;
    q_char = 0;
    for (std::size_t i = 0; i < rhs.size(); ++i) {
      char c = rhs[i];
      if (in_q) {
        tok.push_back(c);
        if (c == q_char) in_q = false;
      } else {
        if (c == '"' || c == '\'') {
          in_q = true;
          q_char = c;
          tok.push_back(c);
        } else if (c == ',' || std::isspace(static_cast<unsigned char>(c))) {
          std::string t = trim(tok);
          if (!t.empty()) vec.push_back(t);
          tok.clear();
        } else {
          tok.push_back(c);
        }
      }
    }
    std::string t = trim(tok);
    if (!t.empty()) vec.push_back(t);
  }
  return true;
}

void init_undefined(NoahmpTables& t) {
  constexpr real u = real(-9999.0);
  auto fill1 = [&](real* a, int n) {
    for (int i = 0; i <= n; ++i) a[i] = u;
  };
  fill1(t.ch2op, kMaxVegCats);
  fill1(t.dleaf, kMaxVegCats);
  fill1(t.z0mvt, kMaxVegCats);
  fill1(t.hvt, kMaxVegCats);
  fill1(t.hvb, kMaxVegCats);
  fill1(t.den, kMaxVegCats);
  fill1(t.rc, kMaxVegCats);
  fill1(t.mfsno, kMaxVegCats);
  fill1(t.scffac, kMaxVegCats);
  fill1(t.cbiom, kMaxVegCats);
  fill1(t.xl, kMaxVegCats);
  fill1(t.cwpvt, kMaxVegCats);
  fill1(t.c3psn, kMaxVegCats);
  fill1(t.kc25, kMaxVegCats);
  fill1(t.akc, kMaxVegCats);
  fill1(t.ko25, kMaxVegCats);
  fill1(t.ako, kMaxVegCats);
  fill1(t.avcmx, kMaxVegCats);
  fill1(t.aqe, kMaxVegCats);
  fill1(t.ltovrc, kMaxVegCats);
  fill1(t.dilefc, kMaxVegCats);
  fill1(t.dilefw, kMaxVegCats);
  fill1(t.rmf25, kMaxVegCats);
  fill1(t.sla, kMaxVegCats);
  fill1(t.fragr, kMaxVegCats);
  fill1(t.tmin, kMaxVegCats);
  fill1(t.vcmx25, kMaxVegCats);
  fill1(t.tdlef, kMaxVegCats);
  fill1(t.bp, kMaxVegCats);
  fill1(t.mp, kMaxVegCats);
  fill1(t.qe25, kMaxVegCats);
  fill1(t.rms25, kMaxVegCats);
  fill1(t.rmr25, kMaxVegCats);
  fill1(t.arm, kMaxVegCats);
  fill1(t.folnmx, kMaxVegCats);
  fill1(t.wdpool, kMaxVegCats);
  fill1(t.wrrat, kMaxVegCats);
  fill1(t.mrp, kMaxVegCats);
  fill1(t.nroot, kMaxVegCats);
  fill1(t.rgl, kMaxVegCats);
  fill1(t.rs, kMaxVegCats);
  fill1(t.hs, kMaxVegCats);
  fill1(t.topt, kMaxVegCats);
  fill1(t.rsmax, kMaxVegCats);
  fill1(t.rtovrc, kMaxVegCats);
  fill1(t.rswoodc, kMaxVegCats);
  fill1(t.bf, kMaxVegCats);
  fill1(t.wstrc, kMaxVegCats);
  fill1(t.laimin, kMaxVegCats);
  fill1(t.xsamin, kMaxVegCats);

  for (int i = 0; i <= kMaxVegCats; ++i) {
    for (int m = 0; m <= kNumMonths; ++m) {
      t.saim[i][m] = u;
      t.laim[i][m] = u;
    }
    for (int b = 0; b <= kNumRadBands; ++b) {
      t.rhol[i][b] = u;
      t.rhos[i][b] = u;
      t.taul[i][b] = u;
      t.taus[i][b] = u;
    }
  }

  fill1(t.bexp, kMaxSoilCats);
  fill1(t.smcdry, kMaxSoilCats);
  fill1(t.smcmax, kMaxSoilCats);
  fill1(t.smcref, kMaxSoilCats);
  fill1(t.psisat, kMaxSoilCats);
  fill1(t.dksat, kMaxSoilCats);
  fill1(t.dwsat, kMaxSoilCats);
  fill1(t.smcwlt, kMaxSoilCats);
  fill1(t.quartz, kMaxSoilCats);
  fill1(t.bvic, kMaxSoilCats);
  fill1(t.axaj, kMaxSoilCats);
  fill1(t.bxaj, kMaxSoilCats);
  fill1(t.xxaj, kMaxSoilCats);
  fill1(t.bdvic, kMaxSoilCats);
  fill1(t.gdvic, kMaxSoilCats);
  fill1(t.bbvic, kMaxSoilCats);

  fill1(t.slope, kMaxSlopeCats);

  for (int c = 0; c <= kMaxSoilColors; ++c) {
    for (int b = 0; b <= kNumRadBands; ++b) {
      t.albsat[c][b] = u;
      t.albdry[c][b] = u;
    }
  }
  fill1(t.albice, kNumRadBands);
  fill1(t.alblak, kNumRadBands);
  fill1(t.omegas, kNumRadBands);
  fill1(t.eg, 2);

  for (int c = 0; c <= kNumCropCats; ++c) {
    t.pltday[c] = -9999;
    t.hsday[c] = -9999;
  }
  fill1(t.plantpop, kNumCropCats);
  fill1(t.irri, kNumCropCats);
  fill1(t.gddtbase, kNumCropCats);
  fill1(t.gddtcut, kNumCropCats);
  fill1(t.gdds1, kNumCropCats);
  fill1(t.gdds2, kNumCropCats);
  fill1(t.gdds3, kNumCropCats);
  fill1(t.gdds4, kNumCropCats);
  fill1(t.gdds5, kNumCropCats);
  fill1(t.c3psni, kNumCropCats);
  fill1(t.kc25i, kNumCropCats);
  fill1(t.akci, kNumCropCats);
  fill1(t.ko25i, kNumCropCats);
  fill1(t.akoi, kNumCropCats);
  fill1(t.avcmxi, kNumCropCats);
  fill1(t.vcmx25i, kNumCropCats);
  fill1(t.bpi, kNumCropCats);
  fill1(t.mpi, kNumCropCats);
  fill1(t.folnmxi, kNumCropCats);
  fill1(t.qe25i, kNumCropCats);
  fill1(t.aref, kNumCropCats);
  fill1(t.psnrf, kNumCropCats);
  fill1(t.i2par, kNumCropCats);
  fill1(t.tassim0, kNumCropCats);
  fill1(t.tassim1, kNumCropCats);
  fill1(t.tassim2, kNumCropCats);
  fill1(t.k_crop, kNumCropCats);
  fill1(t.epsi, kNumCropCats);
  fill1(t.q10mr, kNumCropCats);
  fill1(t.lefreez, kNumCropCats);
  fill1(t.fra_gr, kNumCropCats);
  fill1(t.lfmr25, kNumCropCats);
  fill1(t.stmr25, kNumCropCats);
  fill1(t.rtmr25, kNumCropCats);
  fill1(t.grainmr25, kNumCropCats);
  fill1(t.bio2lai, kNumCropCats);
  for (int c = 0; c <= kNumCropCats; ++c) {
    for (int s = 0; s <= kNumCropStages; ++s) {
      t.dile_fc[c][s] = u;
      t.dile_fw[c][s] = u;
      t.lf_ovrc[c][s] = u;
      t.st_ovrc[c][s] = u;
      t.rt_ovrc[c][s] = u;
      t.lfpt[c][s] = u;
      t.stpt[c][s] = u;
      t.rtpt[c][s] = u;
      t.grainpt[c][s] = u;
      t.lfct[c][s] = u;
      t.stct[c][s] = u;
      t.rtct[c][s] = u;
    }
  }

  for (int s = 0; s <= kMaxSoilCats; ++s) t.td_depth[s] = -9999;
  fill1(t.tdsmc_fac, kMaxSoilCats);
  fill1(t.td_dc, kMaxSoilCats);
  fill1(t.td_dcoef, kMaxSoilCats);
  fill1(t.td_d, kMaxSoilCats);
  fill1(t.td_adepth, kMaxSoilCats);
  fill1(t.td_radi, kMaxSoilCats);
  fill1(t.td_spac, kMaxSoilCats);
  fill1(t.td_ddrain, kMaxSoilCats);
  fill1(t.klat_fac, kMaxSoilCats);
}

int get_int(const SectionMap& sec, const char* key, int def = -9999) {
  auto it = sec.find(key);
  if (it == sec.end() || it->second.empty()) return def;
  return parse_int(it->second[0]);
}

real get_real(const SectionMap& sec, const char* key, real def = real(-9999.0)) {
  auto it = sec.find(key);
  if (it == sec.end() || it->second.empty()) return def;
  return parse_real(it->second[0]);
}

void get_real_arr(const SectionMap& sec, const char* key, real* dst, int n) {
  auto it = sec.find(key);
  if (it == sec.end()) return;
  int cnt = std::min(n, static_cast<int>(it->second.size()));
  for (int i = 0; i < cnt; ++i) dst[i + 1] = parse_real(it->second[i]);
}

void get_int_arr_as_real(const SectionMap& sec, const char* key, real* dst, int n) {
  auto it = sec.find(key);
  if (it == sec.end()) return;
  int cnt = std::min(n, static_cast<int>(it->second.size()));
  for (int i = 0; i < cnt; ++i) dst[i + 1] = static_cast<real>(parse_int(it->second[i]));
}

void get_int_arr(const SectionMap& sec, const char* key, int* dst, int n) {
  auto it = sec.find(key);
  if (it == sec.end()) return;
  int cnt = std::min(n, static_cast<int>(it->second.size()));
  for (int i = 0; i < cnt; ++i) dst[i + 1] = parse_int(it->second[i]);
}

void get_real_col(const SectionMap& sec, const char* key, real* col_base, int stride, int col, int n) {
  auto it = sec.find(key);
  if (it == sec.end()) return;
  int cnt = std::min(n, static_cast<int>(it->second.size()));
  for (int i = 0; i < cnt; ++i) col_base[(i + 1) * stride + col] = parse_real(it->second[i]);
}

}  // namespace

bool read_noahmp_tables(const std::string& table_dir,
                        const std::string& veg_dataset,
                        NoahmpTables& t) {
  init_undefined(t);

  std::string path = table_dir;
  if (path.size() < 4 || path.substr(path.size() - 4) != ".TBL") {
    path += "/NoahmpTable.TBL";
  }

  TableData data;
  if (!parse_table_file(path, data)) return false;

  // 1. Vegetation parameters
  std::string vcat_sec, vparm_sec;
  std::string vds = trim(veg_dataset);
  if (vds == "USGS") {
    vcat_sec = "noahmp_usgs_veg_categories";
    vparm_sec = "noahmp_usgs_parameters";
  } else if (vds == "MODIFIED_IGBP_MODIS_NOAH") {
    vcat_sec = "noahmp_modis_veg_categories";
    vparm_sec = "noahmp_modis_parameters";
  } else {
    return false;
  }

  const auto& vcat = data[vcat_sec];
  const auto& vp = data[vparm_sec];
  const int nveg = get_int(vcat, "NVEG", 0);
  if (nveg <= 0 || nveg > kMaxVegCats) return false;
  t.nveg = nveg;

  t.isurban = get_int(vp, "ISURBAN");
  t.iswater = get_int(vp, "ISWATER");
  t.isbarren = get_int(vp, "ISBARREN");
  t.isice = get_int(vp, "ISICE");
  t.iscrop = get_int(vp, "ISCROP");
  t.eblforest = get_int(vp, "EBLFOREST");
  t.natural = get_int(vp, "NATURAL");
  t.urbtype_beg = get_int(vp, "URBTYPE_BEG");
  t.lcz_1 = get_int(vp, "LCZ_1");
  t.lcz_2 = get_int(vp, "LCZ_2");
  t.lcz_3 = get_int(vp, "LCZ_3");
  t.lcz_4 = get_int(vp, "LCZ_4");
  t.lcz_5 = get_int(vp, "LCZ_5");
  t.lcz_6 = get_int(vp, "LCZ_6");
  t.lcz_7 = get_int(vp, "LCZ_7");
  t.lcz_8 = get_int(vp, "LCZ_8");
  t.lcz_9 = get_int(vp, "LCZ_9");
  t.lcz_10 = get_int(vp, "LCZ_10");
  t.lcz_11 = get_int(vp, "LCZ_11");

  get_real_arr(vp, "CH2OP", t.ch2op, nveg);
  get_real_arr(vp, "DLEAF", t.dleaf, nveg);
  get_real_arr(vp, "Z0MVT", t.z0mvt, nveg);
  get_real_arr(vp, "HVT", t.hvt, nveg);
  get_real_arr(vp, "HVB", t.hvb, nveg);
  get_real_arr(vp, "DEN", t.den, nveg);
  get_real_arr(vp, "RC", t.rc, nveg);
  get_real_arr(vp, "MFSNO", t.mfsno, nveg);
  get_real_arr(vp, "SCFFAC", t.scffac, nveg);
  get_real_arr(vp, "CBIOM", t.cbiom, nveg);
  get_real_arr(vp, "XL", t.xl, nveg);
  get_real_arr(vp, "CWPVT", t.cwpvt, nveg);
  get_real_arr(vp, "C3PSN", t.c3psn, nveg);
  get_real_arr(vp, "KC25", t.kc25, nveg);
  get_real_arr(vp, "AKC", t.akc, nveg);
  get_real_arr(vp, "KO25", t.ko25, nveg);
  get_real_arr(vp, "AKO", t.ako, nveg);
  get_real_arr(vp, "AVCMX", t.avcmx, nveg);
  get_real_arr(vp, "AQE", t.aqe, nveg);
  get_real_arr(vp, "LTOVRC", t.ltovrc, nveg);
  get_real_arr(vp, "DILEFC", t.dilefc, nveg);
  get_real_arr(vp, "DILEFW", t.dilefw, nveg);
  get_real_arr(vp, "RMF25", t.rmf25, nveg);
  get_real_arr(vp, "SLA", t.sla, nveg);
  get_real_arr(vp, "FRAGR", t.fragr, nveg);
  get_real_arr(vp, "TMIN", t.tmin, nveg);
  get_real_arr(vp, "VCMX25", t.vcmx25, nveg);
  get_real_arr(vp, "TDLEF", t.tdlef, nveg);
  get_real_arr(vp, "BP", t.bp, nveg);
  get_real_arr(vp, "MP", t.mp, nveg);
  get_real_arr(vp, "QE25", t.qe25, nveg);
  get_real_arr(vp, "RMS25", t.rms25, nveg);
  get_real_arr(vp, "RMR25", t.rmr25, nveg);
  get_real_arr(vp, "ARM", t.arm, nveg);
  get_real_arr(vp, "FOLNMX", t.folnmx, nveg);
  get_real_arr(vp, "WDPOOL", t.wdpool, nveg);
  get_real_arr(vp, "WRRAT", t.wrrat, nveg);
  get_real_arr(vp, "MRP", t.mrp, nveg);
  get_int_arr_as_real(vp, "NROOT", t.nroot, nveg);
  get_real_arr(vp, "RGL", t.rgl, nveg);
  get_real_arr(vp, "RS", t.rs, nveg);
  get_real_arr(vp, "HS", t.hs, nveg);
  get_real_arr(vp, "TOPT", t.topt, nveg);
  get_real_arr(vp, "RSMAX", t.rsmax, nveg);
  get_real_arr(vp, "RTOVRC", t.rtovrc, nveg);
  get_real_arr(vp, "RSWOODC", t.rswoodc, nveg);
  get_real_arr(vp, "BF", t.bf, nveg);
  get_real_arr(vp, "WSTRC", t.wstrc, nveg);
  get_real_arr(vp, "LAIMIN", t.laimin, nveg);
  get_real_arr(vp, "XSAMIN", t.xsamin, nveg);

  const char* months[12] = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN",
                            "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};
  for (int m = 0; m < 12; ++m) {
    std::string ks = std::string("SAI_") + months[m];
    std::string kl = std::string("LAI_") + months[m];
    get_real_col(vp, ks.c_str(), &t.saim[0][0], kNumMonths + 1, m + 1, nveg);
    get_real_col(vp, kl.c_str(), &t.laim[0][0], kNumMonths + 1, m + 1, nveg);
  }
  get_real_col(vp, "RHOL_VIS", &t.rhol[0][0], kNumRadBands + 1, 1, nveg);
  get_real_col(vp, "RHOL_NIR", &t.rhol[0][0], kNumRadBands + 1, 2, nveg);
  get_real_col(vp, "RHOS_VIS", &t.rhos[0][0], kNumRadBands + 1, 1, nveg);
  get_real_col(vp, "RHOS_NIR", &t.rhos[0][0], kNumRadBands + 1, 2, nveg);
  get_real_col(vp, "TAUL_VIS", &t.taul[0][0], kNumRadBands + 1, 1, nveg);
  get_real_col(vp, "TAUL_NIR", &t.taul[0][0], kNumRadBands + 1, 2, nveg);
  get_real_col(vp, "TAUS_VIS", &t.taus[0][0], kNumRadBands + 1, 1, nveg);
  get_real_col(vp, "TAUS_NIR", &t.taus[0][0], kNumRadBands + 1, 2, nveg);

  // 2. Soil parameters
  const auto& scat = data["noahmp_stas_soil_categories"];
  std::string sltype = "STAS";
  if (auto it = scat.find("SLTYPE"); it != scat.end() && !it->second.empty()) {
    sltype = strip_quotes(it->second[0]);
  }
  const int slcats = get_int(scat, "SLCATS", 0);
  if (slcats <= 0 || slcats > kMaxSoilCats) return false;
  t.slcats = slcats;

  std::string sparm_sec = (sltype == "STAS_RUC") ? "noahmp_soil_stas_ruc_parameters"
                                                 : "noahmp_soil_stas_parameters";
  const auto& sp = data[sparm_sec];
  get_real_arr(sp, "BB", t.bexp, slcats);
  get_real_arr(sp, "DRYSMC", t.smcdry, slcats);
  get_real_arr(sp, "MAXSMC", t.smcmax, slcats);
  get_real_arr(sp, "REFSMC", t.smcref, slcats);
  get_real_arr(sp, "SATPSI", t.psisat, slcats);
  get_real_arr(sp, "SATDK", t.dksat, slcats);
  get_real_arr(sp, "SATDW", t.dwsat, slcats);
  get_real_arr(sp, "WLTSMC", t.smcwlt, slcats);
  get_real_arr(sp, "QTZ", t.quartz, slcats);
  get_real_arr(sp, "BVIC", t.bvic, slcats);
  get_real_arr(sp, "AXAJ", t.axaj, slcats);
  get_real_arr(sp, "BXAJ", t.bxaj, slcats);
  get_real_arr(sp, "XXAJ", t.xxaj, slcats);
  get_real_arr(sp, "BDVIC", t.bdvic, slcats);
  get_real_arr(sp, "GDVIC", t.gdvic, slcats);
  get_real_arr(sp, "BBVIC", t.bbvic, slcats);

  // 3. General parameters
  const auto& gp = data["noahmp_general_parameters"];
  get_real_arr(gp, "SLOPE_DATA", t.slope, kMaxSlopeCats);
  t.csoil = get_real(gp, "CSOIL_DATA");
  t.refdk = get_real(gp, "REFDK_DATA");
  t.refkdt = get_real(gp, "REFKDT_DATA");
  t.frzk = get_real(gp, "FRZK_DATA");
  t.zbot = get_real(gp, "ZBOT_DATA");
  t.czil = get_real(gp, "CZIL_DATA");

  // 4. Radiation parameters
  const auto& rp = data["noahmp_rad_parameters"];
  get_real_col(rp, "ALBSAT_VIS", &t.albsat[0][0], kNumRadBands + 1, 1, kMaxSoilColors);
  get_real_col(rp, "ALBSAT_NIR", &t.albsat[0][0], kNumRadBands + 1, 2, kMaxSoilColors);
  get_real_col(rp, "ALBDRY_VIS", &t.albdry[0][0], kNumRadBands + 1, 1, kMaxSoilColors);
  get_real_col(rp, "ALBDRY_NIR", &t.albdry[0][0], kNumRadBands + 1, 2, kMaxSoilColors);
  get_real_arr(rp, "ALBICE", t.albice, kNumRadBands);
  get_real_arr(rp, "ALBLAK", t.alblak, kNumRadBands);
  get_real_arr(rp, "OMEGAS", t.omegas, kNumRadBands);
  t.betads = get_real(rp, "BETADS");
  t.betais = get_real(rp, "BETAIS");
  get_real_arr(rp, "EG", t.eg, 2);
  t.eice = get_real(rp, "EICE");

  // 5. Global parameters
  const auto& glp = data["noahmp_global_parameters"];
  t.co2 = get_real(glp, "CO2");
  t.o2 = get_real(glp, "O2");
  t.timean = get_real(glp, "TIMEAN");
  t.fsatmx = get_real(glp, "FSATMX");
  t.z0sno = get_real(glp, "Z0SNO");
  t.ssi = get_real(glp, "SSI");
  t.snow_ret_fac = get_real(glp, "SNOW_RET_FAC");
  t.snow_emis = get_real(glp, "SNOW_EMIS");
  t.swemx = get_real(glp, "SWEMX");
  t.tau0 = get_real(glp, "TAU0");
  t.grain_growth = get_real(glp, "GRAIN_GROWTH");
  t.extra_growth = get_real(glp, "EXTRA_GROWTH");
  t.dirt_soot = get_real(glp, "DIRT_SOOT");
  t.bats_cosz = get_real(glp, "BATS_COSZ");
  t.bats_vis_new = get_real(glp, "BATS_VIS_NEW");
  t.bats_nir_new = get_real(glp, "BATS_NIR_NEW");
  t.bats_vis_age = get_real(glp, "BATS_VIS_AGE");
  t.bats_nir_age = get_real(glp, "BATS_NIR_AGE");
  t.bats_vis_dir = get_real(glp, "BATS_VIS_DIR");
  t.bats_nir_dir = get_real(glp, "BATS_NIR_DIR");
  t.rsurf_snow = get_real(glp, "RSURF_SNOW");
  t.rsurf_exp = get_real(glp, "RSURF_EXP");
  t.c2_snowcompact = get_real(glp, "C2_SNOWCOMPACT");
  t.c3_snowcompact = get_real(glp, "C3_SNOWCOMPACT");
  t.c4_snowcompact = get_real(glp, "C4_SNOWCOMPACT");
  t.c5_snowcompact = get_real(glp, "C5_SNOWCOMPACT");
  t.dm_snowcompact = get_real(glp, "DM_SNOWCOMPACT");
  t.eta0_snowcompact = get_real(glp, "ETA0_SNOWCOMPACT");
  t.snowcompactm_ar24 = get_real(glp, "SNOWCOMPACTM_AR24");
  t.snowcompactb_ar24 = get_real(glp, "SNOWCOMPACTB_AR24");
  t.snowcompact_p1_ar24 = get_real(glp, "SNOWCOMPACT_P1_AR24");
  t.snowcompact_p2_ar24 = get_real(glp, "SNOWCOMPACT_P2_AR24");
  t.snowcompact_p3_ar24 = get_real(glp, "SNOWCOMPACT_P3_AR24");
  t.snowcompact_up_ar24 = get_real(glp, "SNOWCOMPACT_UP_AR24");
  t.scfm1_ar25 = get_real(glp, "SCFM1_AR25");
  t.scfm2_ar25 = get_real(glp, "SCFM2_AR25");
  t.scfac1_ar25 = get_real(glp, "SCFAC1_AR25");
  t.scfac2_ar25 = get_real(glp, "SCFAC2_AR25");
  t.snliqmaxfrac = get_real(glp, "SNLIQMAXFRAC");
  t.swemaxgla = get_real(glp, "SWEMAXGLA");
  t.wslmax = get_real(glp, "WSLMAX");
  t.rous = get_real(glp, "ROUS");
  t.cmic = get_real(glp, "CMIC");
  t.snowden_max = get_real(glp, "SNOWDEN_MAX");
  t.class_alb_ref = get_real(glp, "CLASS_ALB_REF");
  t.class_sno_age = get_real(glp, "CLASS_SNO_AGE");
  t.class_alb_new = get_real(glp, "CLASS_ALB_NEW");
  t.psiwlt = get_real(glp, "PSIWLT");
  t.z0soil = get_real(glp, "Z0SOIL");
  t.z0lake = get_real(glp, "Z0LAKE");
  t.wcap = get_real(glp, "WCAP");

  // 6. Irrigation parameters
  const auto& ip = data["noahmp_irrigation_parameters"];
  t.irr_frac = get_real(ip, "IRR_FRAC");
  t.irr_har = get_int(ip, "IRR_HAR");
  t.irr_lai = get_real(ip, "IRR_LAI");
  t.irr_mad = get_real(ip, "IRR_MAD");
  t.filoss = get_real(ip, "FILOSS");
  t.sprir_rate = get_real(ip, "SPRIR_RATE");
  t.micir_rate = get_real(ip, "MICIR_RATE");
  t.firtfac = get_real(ip, "FIRTFAC");
  t.ir_rain = get_real(ip, "IR_RAIN");

  // 7. Crop parameters
  const auto& cp = data["noahmp_crop_parameters"];
  t.default_crop = get_int(cp, "DEFAULT_CROP");
  get_int_arr(cp, "PLTDAY", t.pltday, kNumCropCats);
  get_int_arr(cp, "HSDAY", t.hsday, kNumCropCats);
  get_real_arr(cp, "PLANTPOP", t.plantpop, kNumCropCats);
  get_real_arr(cp, "GDDTBASE", t.gddtbase, kNumCropCats);
  get_real_arr(cp, "GDDTCUT", t.gddtcut, kNumCropCats);
  get_real_arr(cp, "GDDS1", t.gdds1, kNumCropCats);
  get_real_arr(cp, "GDDS2", t.gdds2, kNumCropCats);
  get_real_arr(cp, "GDDS3", t.gdds3, kNumCropCats);
  get_real_arr(cp, "GDDS4", t.gdds4, kNumCropCats);
  get_real_arr(cp, "GDDS5", t.gdds5, kNumCropCats);
  get_real_arr(cp, "C3PSNI", t.c3psni, kNumCropCats);
  get_real_arr(cp, "KC25I", t.kc25i, kNumCropCats);
  get_real_arr(cp, "AKCI", t.akci, kNumCropCats);
  get_real_arr(cp, "KO25I", t.ko25i, kNumCropCats);
  get_real_arr(cp, "AKOI", t.akoi, kNumCropCats);
  get_real_arr(cp, "AVCMXI", t.avcmxi, kNumCropCats);
  get_real_arr(cp, "VCMX25I", t.vcmx25i, kNumCropCats);
  get_real_arr(cp, "BPI", t.bpi, kNumCropCats);
  get_real_arr(cp, "MPI", t.mpi, kNumCropCats);
  get_real_arr(cp, "FOLNMXI", t.folnmxi, kNumCropCats);
  get_real_arr(cp, "QE25I", t.qe25i, kNumCropCats);
  get_real_arr(cp, "AREF", t.aref, kNumCropCats);
  get_real_arr(cp, "PSNRF", t.psnrf, kNumCropCats);
  get_real_arr(cp, "I2PAR", t.i2par, kNumCropCats);
  get_real_arr(cp, "TASSIM0", t.tassim0, kNumCropCats);
  get_real_arr(cp, "TASSIM1", t.tassim1, kNumCropCats);
  get_real_arr(cp, "TASSIM2", t.tassim2, kNumCropCats);
  get_real_arr(cp, "K", t.k_crop, kNumCropCats);
  get_real_arr(cp, "EPSI", t.epsi, kNumCropCats);
  get_real_arr(cp, "Q10MR", t.q10mr, kNumCropCats);
  get_real_arr(cp, "LEFREEZ", t.lefreez, kNumCropCats);
  get_real_arr(cp, "FRA_GR", t.fra_gr, kNumCropCats);
  get_real_arr(cp, "LFMR25", t.lfmr25, kNumCropCats);
  get_real_arr(cp, "STMR25", t.stmr25, kNumCropCats);
  get_real_arr(cp, "RTMR25", t.rtmr25, kNumCropCats);
  get_real_arr(cp, "GRAINMR25", t.grainmr25, kNumCropCats);
  get_real_arr(cp, "BIO2LAI", t.bio2lai, kNumCropCats);
  for (int s = 1; s <= kNumCropStages; ++s) {
    std::string sfx = "_S" + std::to_string(s);
    get_real_col(cp, ("DILE_FC" + sfx).c_str(), &t.dile_fc[0][0], kNumCropStages + 1, s, kNumCropCats);
    get_real_col(cp, ("DILE_FW" + sfx).c_str(), &t.dile_fw[0][0], kNumCropStages + 1, s, kNumCropCats);
    get_real_col(cp, ("LF_OVRC" + sfx).c_str(), &t.lf_ovrc[0][0], kNumCropStages + 1, s, kNumCropCats);
    get_real_col(cp, ("ST_OVRC" + sfx).c_str(), &t.st_ovrc[0][0], kNumCropStages + 1, s, kNumCropCats);
    get_real_col(cp, ("RT_OVRC" + sfx).c_str(), &t.rt_ovrc[0][0], kNumCropStages + 1, s, kNumCropCats);
    get_real_col(cp, ("LFPT" + sfx).c_str(), &t.lfpt[0][0], kNumCropStages + 1, s, kNumCropCats);
    get_real_col(cp, ("STPT" + sfx).c_str(), &t.stpt[0][0], kNumCropStages + 1, s, kNumCropCats);
    get_real_col(cp, ("RTPT" + sfx).c_str(), &t.rtpt[0][0], kNumCropStages + 1, s, kNumCropCats);
    get_real_col(cp, ("GRAINPT" + sfx).c_str(), &t.grainpt[0][0], kNumCropStages + 1, s, kNumCropCats);
    get_real_col(cp, ("LFCT" + sfx).c_str(), &t.lfct[0][0], kNumCropStages + 1, s, kNumCropCats);
    get_real_col(cp, ("STCT" + sfx).c_str(), &t.stct[0][0], kNumCropStages + 1, s, kNumCropCats);
    get_real_col(cp, ("RTCT" + sfx).c_str(), &t.rtct[0][0], kNumCropStages + 1, s, kNumCropCats);
  }

  // 8. Tile drainage parameters
  const auto& tdp = data["noahmp_tiledrain_parameters"];
  const int nsoiltype = get_int(tdp, "NSOILTYPE", 0);
  t.drain_layer_opt = get_int(tdp, "DRAIN_LAYER_OPT");
  if (nsoiltype > 0 && nsoiltype <= kMaxSoilCats) {
    get_real_arr(tdp, "TDSMC_FAC", t.tdsmc_fac, nsoiltype);
    get_int_arr(tdp, "TD_DEPTH", t.td_depth, nsoiltype);
    get_real_arr(tdp, "TD_DC", t.td_dc, nsoiltype);
    get_real_arr(tdp, "TD_DCOEF", t.td_dcoef, nsoiltype);
    get_real_arr(tdp, "TD_D", t.td_d, nsoiltype);
    get_real_arr(tdp, "TD_ADEPTH", t.td_adepth, nsoiltype);
    get_real_arr(tdp, "TD_RADI", t.td_radi, nsoiltype);
    get_real_arr(tdp, "TD_SPAC", t.td_spac, nsoiltype);
    get_real_arr(tdp, "TD_DDRAIN", t.td_ddrain, nsoiltype);
    get_real_arr(tdp, "KLAT_FAC", t.klat_fac, nsoiltype);
  }

  // 9. Optional pedotransfer parameters
  const auto& op = data["noahmp_optional_parameters"];
  t.sr2006_theta_1500t_a = get_real(op, "SR2006_THETA_1500T_A");
  t.sr2006_theta_1500t_b = get_real(op, "SR2006_THETA_1500T_B");
  t.sr2006_theta_1500t_c = get_real(op, "SR2006_THETA_1500T_C");
  t.sr2006_theta_1500t_d = get_real(op, "SR2006_THETA_1500T_D");
  t.sr2006_theta_1500t_e = get_real(op, "SR2006_THETA_1500T_E");
  t.sr2006_theta_1500t_f = get_real(op, "SR2006_THETA_1500T_F");
  t.sr2006_theta_1500t_g = get_real(op, "SR2006_THETA_1500T_G");
  t.sr2006_theta_1500_a = get_real(op, "SR2006_THETA_1500_A");
  t.sr2006_theta_1500_b = get_real(op, "SR2006_THETA_1500_B");
  t.sr2006_theta_33t_a = get_real(op, "SR2006_THETA_33T_A");
  t.sr2006_theta_33t_b = get_real(op, "SR2006_THETA_33T_B");
  t.sr2006_theta_33t_c = get_real(op, "SR2006_THETA_33T_C");
  t.sr2006_theta_33t_d = get_real(op, "SR2006_THETA_33T_D");
  t.sr2006_theta_33t_e = get_real(op, "SR2006_THETA_33T_E");
  t.sr2006_theta_33t_f = get_real(op, "SR2006_THETA_33T_F");
  t.sr2006_theta_33t_g = get_real(op, "SR2006_THETA_33T_G");
  t.sr2006_theta_33_a = get_real(op, "SR2006_THETA_33_A");
  t.sr2006_theta_33_b = get_real(op, "SR2006_THETA_33_B");
  t.sr2006_theta_33_c = get_real(op, "SR2006_THETA_33_C");
  t.sr2006_theta_s33t_a = get_real(op, "SR2006_THETA_S33T_A");
  t.sr2006_theta_s33t_b = get_real(op, "SR2006_THETA_S33T_B");
  t.sr2006_theta_s33t_c = get_real(op, "SR2006_THETA_S33T_C");
  t.sr2006_theta_s33t_d = get_real(op, "SR2006_THETA_S33T_D");
  t.sr2006_theta_s33t_e = get_real(op, "SR2006_THETA_S33T_E");
  t.sr2006_theta_s33t_f = get_real(op, "SR2006_THETA_S33T_F");
  t.sr2006_theta_s33t_g = get_real(op, "SR2006_THETA_S33T_G");
  t.sr2006_theta_s33_a = get_real(op, "SR2006_THETA_S33_A");
  t.sr2006_theta_s33_b = get_real(op, "SR2006_THETA_S33_B");
  t.sr2006_psi_et_a = get_real(op, "SR2006_PSI_ET_A");
  t.sr2006_psi_et_b = get_real(op, "SR2006_PSI_ET_B");
  t.sr2006_psi_et_c = get_real(op, "SR2006_PSI_ET_C");
  t.sr2006_psi_et_d = get_real(op, "SR2006_PSI_ET_D");
  t.sr2006_psi_et_e = get_real(op, "SR2006_PSI_ET_E");
  t.sr2006_psi_et_f = get_real(op, "SR2006_PSI_ET_F");
  t.sr2006_psi_et_g = get_real(op, "SR2006_PSI_ET_G");
  t.sr2006_psi_e_a = get_real(op, "SR2006_PSI_E_A");
  t.sr2006_psi_e_b = get_real(op, "SR2006_PSI_E_B");
  t.sr2006_psi_e_c = get_real(op, "SR2006_PSI_E_C");
  t.sr2006_smcmax_a = get_real(op, "SR2006_SMCMAX_A");
  t.sr2006_smcmax_b = get_real(op, "SR2006_SMCMAX_B");

  // 10. SNICAR parameters
  const auto& snp = data["noahmp_snicar_parameters"];
  t.depbchydropho = get_real(snp, "DEPBCHYDROPHO");
  t.depbchydrophi = get_real(snp, "DEPBCHYDROPHI");
  t.depochydropho = get_real(snp, "DEPOCHYDROPHO");
  t.depochydrophi = get_real(snp, "DEPOCHYDROPHI");
  t.depdust1 = get_real(snp, "DEPDUST1");
  t.depdust2 = get_real(snp, "DEPDUST2");
  t.depdust3 = get_real(snp, "DEPDUST3");
  t.depdust4 = get_real(snp, "DEPDUST4");
  t.depdust5 = get_real(snp, "DEPDUST5");
  t.snowradiusmin = get_real(snp, "SNOWRADIUSMIN");
  t.freshsnowradiusmax = get_real(snp, "FRESHSNOWRADIUSMAX");
  t.snowradiusrefrz = get_real(snp, "SNOWRADIUSREFRZ");
  t.scaveffmeltscale = get_real(snp, "SCAVEFFMELTSCALE");
  t.scaveffmeltbcphi = get_real(snp, "SCAVEFFMELTBCPHI");
  t.scaveffmeltbcpho = get_real(snp, "SCAVEFFMELTBCPHO");
  t.scaveffmeltocphi = get_real(snp, "SCAVEFFMELTOCPHI");
  t.scaveffmeltocpho = get_real(snp, "SCAVEFFMELTOCPHO");
  t.scaveffmeltdust1 = get_real(snp, "SCAVEFFMELTDUST1");
  t.scaveffmeltdust2 = get_real(snp, "SCAVEFFMELTDUST2");
  t.scaveffmeltdust3 = get_real(snp, "SCAVEFFMELTDUST3");
  t.scaveffmeltdust4 = get_real(snp, "SCAVEFFMELTDUST4");
  t.scaveffmeltdust5 = get_real(snp, "SCAVEFFMELTDUST5");
  t.snowradiusmax = get_real(snp, "SNOWRADIUSMAX");
  t.snowwetagec1brun89 = get_real(snp, "SNOWWETAGEC1BRUN89");
  t.snowwetagec2brun89 = get_real(snp, "SNOWWETAGEC2BRUN89");
  t.snowagescalefac = get_real(snp, "SNOWAGESCALEFAC");

  return true;
}

}  // namespace cppwrf::phys::noahmp
