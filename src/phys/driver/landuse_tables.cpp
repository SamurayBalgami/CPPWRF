#include "phys/driver/landuse_tables.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace cppwrf::phys::driver {

namespace {

std::string trim(const std::string& s) {
  auto start = s.find_first_not_of(" \t\r\n");
  if (start == std::string::npos) return "";
  auto end = s.find_last_not_of(" \t\r\n");
  return s.substr(start, end - start + 1);
}

std::string to_upper(std::string s) {
  for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

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

  const std::string target = to_upper(trim(dataset_name));
  std::string line;
  bool found = false;

  while (std::getline(f, line)) {
    std::string s = to_upper(trim(line));
    if (s.empty()) continue;
    if (s == target) {
      found = true;
      break;
    }
  }

  if (!found) {
    std::cerr << "LANDUSE dataset '" << dataset_name << "' not found in " << filepath << "\n";
    return false;
  }

  tables.dataset_name = dataset_name;

  // Header line: LUCATS, LUSEAS, 'ALBD SLMO ...'
  if (!std::getline(f, line)) return false;
  auto hdr_tokens = tokenize(line);
  if (hdr_tokens.size() < 2) return false;

  int n_cats = std::stoi(hdr_tokens[0]);
  int n_seasons = std::stoi(hdr_tokens[1]);

  tables.n_cats = n_cats;
  tables.n_seasons = n_seasons;

  for (int s = 0; s < n_seasons; ++s) {
    // Season name line (e.g. SUMMER, WINTER)
    if (!std::getline(f, line)) return false;

    for (int c = 1; c <= n_cats; ++c) {
      if (!std::getline(f, line)) return false;
      auto tokens = tokenize(line);
      if (tokens.size() < 8) return false;

      int cat_idx = std::stoi(tokens[0]);
      if (cat_idx >= kMaxLanduseCats) continue;

      LanduseCategory& cat = tables.categories[s][cat_idx];
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

  return true;
}

}  // namespace cppwrf::phys::driver
