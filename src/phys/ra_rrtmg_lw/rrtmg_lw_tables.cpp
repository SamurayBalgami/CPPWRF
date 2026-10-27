#include "phys/ra_rrtmg_lw/rrtmg_lw_tables.hpp"
#include <fstream>
#include <iostream>

namespace cppwrf::phys {

bool load_rrtmg_lw_tables(const std::string& filepath, RrtmgLwTables& tables) {
  std::ifstream file(filepath, std::ios::binary);
  if (file.is_open()) {
    file.read(reinterpret_cast<char*>(&tables), sizeof(RrtmgLwTables));
    if (file.gcount() == static_cast<std::streamsize>(sizeof(RrtmgLwTables))) return true;
  }
  std::cerr << "Error: could not read rrtmg_lw_tables.bin from " << filepath << "\n";
  return false;
}

}  // namespace cppwrf::phys
