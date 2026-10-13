#include "phys/lsm_noah/noah_tables.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace cppwrf::phys::noah {

namespace {

std::string trim(const std::string& s) {
  auto start = s.find_first_not_of(" \t\r\n");
  if (start == std::string::npos) return "";
  auto end = s.find_last_not_of(" \t\r\n");
  return s.substr(start, end - start + 1);
}

// Tokenize a line separated by spaces or commas
std::vector<std::string> tokenize(const std::string& line) {
  std::vector<std::string> tokens;
  std::string cur;
  for (char c : line) {
    if (c == ',' || std::isspace(static_cast<unsigned char>(c))) {
      if (!cur.empty()) {
        tokens.push_back(cur);
        cur.clear();
      }
    } else {
      cur += c;
    }
  }
  if (!cur.empty()) tokens.push_back(cur);
  return tokens;
}

bool read_genparm(const std::string& path, NoahTables& tables) {
  std::ifstream f(path);
  if (!f.is_open()) {
    std::cerr << "Cannot open GENPARM table: " << path << std::endl;
    return false;
  }
  std::string line;
  // Look for SLOPE_DATA
  while (std::getline(f, line)) {
    if (trim(line) == "SLOPE_DATA") break;
  }
  if (!std::getline(f, line)) return false;
  int n_slope = std::stoi(trim(line));
  tables.slpcats = n_slope;
  for (int i = 0; i < n_slope && i < kMaxSlopeCats; ++i) {
    if (!std::getline(f, line)) return false;
    tables.slope_data[i] = static_cast<real>(std::stod(trim(line)));
  }

  while (std::getline(f, line)) {
    std::string key = trim(line);
    if (key.empty()) continue;
    if (!std::getline(f, line)) break;
    double val = std::stod(trim(line));
    if (key == "SBETA_DATA") tables.sbeta_data = static_cast<real>(val);
    else if (key == "FXEXP_DATA") tables.fxexp_data = static_cast<real>(val);
    else if (key == "CSOIL_DATA") tables.csoil_data = static_cast<real>(val);
    else if (key == "SALP_DATA") tables.salp_data = static_cast<real>(val);
    else if (key == "REFDK_DATA") tables.refdk_data = static_cast<real>(val);
    else if (key == "REFKDT_DATA") tables.refkdt_data = static_cast<real>(val);
    else if (key == "FRZK_DATA") tables.frzk_data = static_cast<real>(val);
    else if (key == "ZBOT_DATA") tables.zbot_data = static_cast<real>(val);
    else if (key == "CZIL_DATA") tables.czil_data = static_cast<real>(val);
    else if (key == "SMLOW_DATA") tables.smlow_data = static_cast<real>(val);
    else if (key == "SMHIGH_DATA") tables.smhigh_data = static_cast<real>(val);
    else if (key == "LVCOEF_DATA") tables.lvcoef_data = static_cast<real>(val);
  }
  return true;
}

bool read_soilparm(const std::string& path, const std::string& dataset, NoahTables& tables) {
  std::ifstream f(path);
  if (!f.is_open()) {
    std::cerr << "Cannot open SOILPARM table: " << path << std::endl;
    return false;
  }
  std::string line;
  bool found = false;
  while (std::getline(f, line)) {
    std::string t = trim(line);
    if (t == dataset) {
      found = true;
      break;
    }
  }
  if (!found) {
    std::cerr << "Soil dataset " << dataset << " not found in " << path << std::endl;
    return false;
  }
  if (!std::getline(f, line)) return false;
  auto tokens = tokenize(line);
  if (tokens.empty()) return false;
  int slcats = std::stoi(tokens[0]);
  tables.slcats = slcats;

  for (int lc = 0; lc < slcats; ++lc) {
    if (!std::getline(f, line)) return false;
    auto toks = tokenize(line);
    if (toks.size() < 11) return false;
    // Format: IINDEX, BB, DRYSMC, F11, MAXSMC, REFSMC, SATPSI, SATDK, SATDW, WLTSMC, QTZ
    tables.bb[lc] = static_cast<real>(std::stod(toks[1]));
    tables.drysmc[lc] = static_cast<real>(std::stod(toks[2]));
    tables.f11[lc] = static_cast<real>(std::stod(toks[3]));
    tables.maxsmc[lc] = static_cast<real>(std::stod(toks[4]));
    tables.refsmc[lc] = static_cast<real>(std::stod(toks[5]));
    tables.satpsi[lc] = static_cast<real>(std::stod(toks[6]));
    tables.satdk[lc] = static_cast<real>(std::stod(toks[7]));
    tables.satdw[lc] = static_cast<real>(std::stod(toks[8]));
    tables.wltsmc[lc] = static_cast<real>(std::stod(toks[9]));
    tables.qtz[lc] = static_cast<real>(std::stod(toks[10]));
  }
  return true;
}

bool read_vegparm(const std::string& path, const std::string& dataset, NoahTables& tables) {
  std::ifstream f(path);
  if (!f.is_open()) {
    std::cerr << "Cannot open VEGPARM table: " << path << std::endl;
    return false;
  }
  std::string line;
  bool found = false;
  while (std::getline(f, line)) {
    std::string t = trim(line);
    if (t == dataset) {
      found = true;
      break;
    }
  }
  if (!found) {
    std::cerr << "Vegetation dataset " << dataset << " not found in " << path << std::endl;
    return false;
  }
  if (!std::getline(f, line)) return false;
  auto tokens = tokenize(line);
  if (tokens.empty()) return false;
  int lucats = std::stoi(tokens[0]);
  tables.lucats = lucats;

  for (int lc = 0; lc < lucats; ++lc) {
    if (!std::getline(f, line)) return false;
    auto toks = tokenize(line);
    if (toks.size() < 17) return false;
    // IINDEX, SHDTBL, NROTBL, RSTBL, RGLTBL, HSTBL, SNUPTBL, MAXALB,
    // LAIMINTBL, LAIMAXTBL, EMISSMINTBL, EMISSMAXTBL, ALBEDOMINTBL, ALBEDOMAXTBL,
    // Z0MINTBL, Z0MAXTBL, ZTOPVTBL, ZBOTVTBL
    tables.shdtbl[lc] = static_cast<real>(std::stod(toks[1]));
    tables.nrotbl[lc] = std::stoi(toks[2]);
    tables.rstbl[lc] = static_cast<real>(std::stod(toks[3]));
    tables.rgltbl[lc] = static_cast<real>(std::stod(toks[4]));
    tables.hstbl[lc] = static_cast<real>(std::stod(toks[5]));
    tables.snuptbl[lc] = static_cast<real>(std::stod(toks[6]));
    tables.maxalb[lc] = static_cast<real>(std::stod(toks[7]));
    tables.laimintbl[lc] = static_cast<real>(std::stod(toks[8]));
    tables.laimaxtbl[lc] = static_cast<real>(std::stod(toks[9]));
    tables.emissmintbl[lc] = static_cast<real>(std::stod(toks[10]));
    tables.emissmaxtbl[lc] = static_cast<real>(std::stod(toks[11]));
    tables.albedomintbl[lc] = static_cast<real>(std::stod(toks[12]));
    tables.albedomaxtbl[lc] = static_cast<real>(std::stod(toks[13]));
    tables.z0mintbl[lc] = static_cast<real>(std::stod(toks[14]));
    tables.z0maxtbl[lc] = static_cast<real>(std::stod(toks[15]));
    if (toks.size() >= 18) {
      tables.ztopvtbl[lc] = static_cast<real>(std::stod(toks[16]));
      tables.zbotvtbl[lc] = static_cast<real>(std::stod(toks[17]));
    }
  }

  while (std::getline(f, line)) {
    std::string key = trim(line);
    if (key.empty()) continue;
    if (key == "TOPT_DATA") {
      if (std::getline(f, line)) tables.topt_data = static_cast<real>(std::stod(trim(line)));
    } else if (key == "CMCMAX_DATA") {
      if (std::getline(f, line)) tables.cmcmax_data = static_cast<real>(std::stod(trim(line)));
    } else if (key == "CFACTR_DATA") {
      if (std::getline(f, line)) tables.cfactr_data = static_cast<real>(std::stod(trim(line)));
    } else if (key == "RSMAX_DATA") {
      if (std::getline(f, line)) tables.rsmax_data = static_cast<real>(std::stod(trim(line)));
    } else if (key == "BARE") {
      if (std::getline(f, line)) tables.bare = std::stoi(trim(line));
    } else if (key == "NATURAL") {
      if (std::getline(f, line)) tables.natural = std::stoi(trim(line));
    } else if (key.rfind("LCZ_", 0) == 0) {
      int n = std::stoi(key.substr(4));
      if (std::getline(f, line) && n >= 1 && n <= 11) tables.lcz[n - 1] = std::stoi(trim(line));
      if (n == 11) break;
    } else if (key.rfind("Vegetation Parameters", 0) == 0) {
      break;
    }
  }
  return true;
}

}  // namespace

bool read_noah_tables(const std::string& table_dir,
                      const std::string& veg_dataset,
                      const std::string& soil_dataset,
                      NoahTables& tables) {
  std::string gen_path = table_dir + "/GENPARM.TBL";
  std::string soil_path = table_dir + "/SOILPARM.TBL";
  std::string veg_path = table_dir + "/VEGPARM.TBL";

  if (!read_genparm(gen_path, tables)) return false;
  if (!read_soilparm(soil_path, soil_dataset, tables)) return false;
  if (!read_vegparm(veg_path, veg_dataset, tables)) return false;
  return true;
}

}  // namespace cppwrf::phys::noah
