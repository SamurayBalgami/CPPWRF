#pragma once

#include <string>
#include <vector>
#include "core/precision.hpp"

namespace cppwrf::phys::driver {

static constexpr int kMaxLanduseCats = 65;
static constexpr int kMaxSeasons = 2; // Summer = 0, Winter = 1

struct LanduseCategory {
  real albd = 0;    // Albedo (%)
  real slmo = 0;    // Soil moisture availability (0-1)
  real sfem = 0;    // Surface emissivity (0-1)
  real sfz0 = 0;    // Roughness length (cm)
  real therin = 0;  // Thermal inertia
  real scfx = 0;    // Snow cover factor
  real sfhc = 0;    // Surface heat capacity
  std::string name;
};

struct LanduseTables {
  std::string dataset_name;
  int n_cats = 0;
  int n_seasons = 0;
  // Indexed as [season][category_1based]
  LanduseCategory categories[kMaxSeasons][kMaxLanduseCats];
};

// Reads LANDUSE.TBL for the specified dataset (e.g. "USGS", "MODIFIED_IGBP_MODIS_NOAH")
bool read_landuse_table(const std::string& filepath,
                        const std::string& dataset_name,
                        LanduseTables& tables);

}  // namespace cppwrf::phys::driver
