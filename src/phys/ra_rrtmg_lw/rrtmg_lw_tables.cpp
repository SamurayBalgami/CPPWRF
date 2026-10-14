#include "phys/ra_rrtmg_lw/rrtmg_lw_tables.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace cppwrf::phys {

bool load_rrtmg_lw_tables(const std::string& filepath, RrtmgLwTables& tables) {
  std::vector<std::string> search_paths;
  if (!filepath.empty()) {
    search_paths.push_back(filepath);
  }
  search_paths.push_back("rrtmg_lw_tables.bin");
  search_paths.push_back("tests/phys/ra_rrtmg_lw/rrtmg_lw_tables.bin");
  search_paths.push_back("tests/phys/ra_rrtmg_lw/build/rrtmg_lw_tables.bin");
  search_paths.push_back("../tests/phys/ra_rrtmg_lw/rrtmg_lw_tables.bin");
  search_paths.push_back("../../tests/phys/ra_rrtmg_lw/rrtmg_lw_tables.bin");

  for (const auto& path : search_paths) {
    std::ifstream file(path, std::ios::binary);
    if (file.is_open()) {
      file.read(reinterpret_cast<char*>(&tables), sizeof(RrtmgLwTables));
      if (file.gcount() == static_cast<std::streamsize>(sizeof(RrtmgLwTables))) {
        return true;
      }
    }
  }

  std::cerr << "Error: Could not locate or read valid rrtmg_lw_tables.bin in search paths.\n";
  return false;
}

RrtmgLwTables& get_default_rrtmg_lw_tables() {
  static RrtmgLwTables default_tables;
  static bool initialized = false;
  if (!initialized) {
    if (!load_rrtmg_lw_tables("", default_tables)) {
      throw std::runtime_error("Failed to load default RRTMG LW tables");
    }
    initialized = true;
  }
  return default_tables;
}

}  // namespace cppwrf::phys
