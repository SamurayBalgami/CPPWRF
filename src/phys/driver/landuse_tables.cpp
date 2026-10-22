#include "phys/driver/landuse_tables.hpp"
#include "core/strings.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace cppwrf::phys::driver {

namespace {

// Tokenize a line separated by spaces or commas, respecting single quotes
std::vector<std::string> tokenize(const std::string& line) {
  std::vector<std::string> tokens;
  std::string cur;
  bool in_quotes = false;

  for (std::size_t i = 0; i < line.size(); ++i) {
    char c = line[i];
    if (c == '\'') {
      in_quotes = !in_quotes;
      continue;
    }
    if (!in_quotes && (c == ',' || std::isspace(static_cast<unsigned char>(c)))) {
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

}  // namespace

bool read_landuse_table(const std::string& filepath,
                        const std::string& dataset_name,
                        LanduseTables& tables) {
  std::ifstream f(filepath);
  if (!f.is_open()) {
    std::cerr << "Cannot open LANDUSE table: " << filepath << "\n";
    return false;
  }

  // landuse_init (module_physics_init.F): records read in order; a dataset is a name record (list-directed:
  // first token, compared case-sensitively with MMINLU), LUCATS,LUSEAS, then per season a title record and
  // LUCATS category records (LC .NE. LI is fatal); other datasets' records are skipped one by one
  const std::string target = trim(dataset_name);
  std::string line;
  auto next_nonblank = [&](std::string& out) {  // list-directed reads skip blank records
    while (std::getline(f, out))
      if (!tokenize(out).empty()) return true;
    return false;
  };
  while (next_nonblank(line)) {
    const std::string lutype = tokenize(line)[0];
    if (!next_nonblank(line)) break;
    auto hdr_tokens = tokenize(line);
    if (hdr_tokens.size() < 2) return false;
    const int n_cats = std::stoi(hdr_tokens[0]);
    const int n_seasons = std::stoi(hdr_tokens[1]);
    const bool found = lutype == target;
    if (found) {
      tables.dataset_name = dataset_name;
      tables.n_cats = n_cats;
      tables.n_seasons = n_seasons;
    }
    for (int sn = 0; sn < n_seasons; ++sn) {
      if (!std::getline(f, line)) return false;  // season title (READ (landuse_unit,*): one record)
      for (int c = 1; c <= n_cats; ++c) {
        if (!found) {
          if (!std::getline(f, line)) return false;
          continue;
        }
        if (!next_nonblank(line)) return false;
        auto tokens = tokenize(line);
        if (tokens.size() < 8) return false;
        const int cat_idx = std::stoi(tokens[0]);
        if (cat_idx != c) {
          std::cerr << "module_start: MISSING LANDUSE UNIT (" << filepath << ")\n";
          return false;
        }
        if (cat_idx >= kMaxLanduseCats) continue;
        LanduseCategory& cat = tables.categories[sn][cat_idx];
        cat.albd = static_cast<real>(std::stod(tokens[1]));
        cat.slmo = static_cast<real>(std::stod(tokens[2]));
        cat.sfem = static_cast<real>(std::stod(tokens[3]));
        cat.sfz0 = static_cast<real>(std::stod(tokens[4]));
        cat.therin = static_cast<real>(std::stod(tokens[5]));
        cat.scfx = static_cast<real>(std::stod(tokens[6]));
        cat.sfhc = static_cast<real>(std::stod(tokens[7]));
        if (tokens.size() >= 9) cat.name = tokens[8];
      }
    }
    if (found) return true;
  }
  std::cerr << "LANDUSE dataset '" << dataset_name << "' not found in " << filepath << "\n";
  return false;
  return true;
}

}  // namespace cppwrf::phys::driver
