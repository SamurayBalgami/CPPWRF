#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace cppwrf::phys {

bool load_rrtmg_sw_tables(const std::string& filepath, RrtmgSwTables& tables) {
  std::vector<std::string> search_paths;
  if (!filepath.empty()) {
    search_paths.push_back(filepath);
  }
  search_paths.push_back("rrtmg_sw_tables.bin");
  search_paths.push_back("tests/phys/ra_rrtmg_sw/rrtmg_sw_tables.bin");
  search_paths.push_back("tests/phys/ra_rrtmg_sw/build/rrtmg_sw_tables.bin");
  search_paths.push_back("../tests/phys/ra_rrtmg_sw/rrtmg_sw_tables.bin");
  search_paths.push_back("../../tests/phys/ra_rrtmg_sw/rrtmg_sw_tables.bin");

  for (const auto& path : search_paths) {
    std::ifstream file(path, std::ios::binary);
    if (file.is_open()) {
      file.read(reinterpret_cast<char*>(&tables), sizeof(RrtmgSwTables));
      if (file.gcount() == static_cast<std::streamsize>(sizeof(RrtmgSwTables))) {
        return true;
      }
    }
  }

  std::cerr << "Error: Could not locate or read valid rrtmg_sw_tables.bin in search paths.\n";
  return false;
}

RrtmgSwTables& get_default_rrtmg_sw_tables() {
  static RrtmgSwTables default_tables;
  static bool initialized = false;
  if (!initialized) {
    if (!load_rrtmg_sw_tables("", default_tables)) {
      throw std::runtime_error("Failed to load default RRTMG SW tables");
    }
    initialized = true;
  }
  return default_tables;
}

}  // namespace cppwrf::phys
