#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"
#include <fstream>
#include <iostream>

namespace cppwrf::phys {

bool load_rrtmg_sw_tables(const std::string& filepath, RrtmgSwTables& tables) {
  std::ifstream file(filepath, std::ios::binary);
  if (file.is_open()) {
    file.read(reinterpret_cast<char*>(&tables), sizeof(RrtmgSwTables));
    if (file.gcount() == static_cast<std::streamsize>(sizeof(RrtmgSwTables))) return true;
  }
  std::cerr << "Error: could not read rrtmg_sw_tables.bin from " << filepath << "\n";
  return false;
}

}  // namespace cppwrf::phys
