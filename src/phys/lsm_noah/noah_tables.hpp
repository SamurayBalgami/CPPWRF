#pragma once

#include <string>
#include "core/precision.hpp"

namespace cppwrf::phys::noah {

static constexpr int kMaxVegCats = 50;
static constexpr int kMaxSoilCats = 30;
static constexpr int kMaxSlopeCats = 30;

struct NoahTables {
  // General parameters
  int slpcats = 0;
  real slope_data[kMaxSlopeCats] = {0};
  real sbeta_data = 0;
  real fxexp_data = 0;
  real csoil_data = 0;
  real salp_data = 0;
  real refdk_data = 0;
  real refkdt_data = 0;
  real frzk_data = 0;
  real zbot_data = 0;
  real czil_data = 0;
  real smlow_data = 0;
  real smhigh_data = 0;
  real lvcoef_data = 0;

  // Soil parameters
  int slcats = 0;
  real bb[kMaxSoilCats] = {0};
  real drysmc[kMaxSoilCats] = {0};
  real f11[kMaxSoilCats] = {0};
  real maxsmc[kMaxSoilCats] = {0};
  real refsmc[kMaxSoilCats] = {0};
  real satpsi[kMaxSoilCats] = {0};
  real satdk[kMaxSoilCats] = {0};
  real satdw[kMaxSoilCats] = {0};
  real wltsmc[kMaxSoilCats] = {0};
  real qtz[kMaxSoilCats] = {0};

  // Vegetation parameters
  int lucats = 0;
  int bare = 0;
  int natural = 0;
  real topt_data = 0;
  real cmcmax_data = 0;
  real cfactr_data = 0;
  real rsmax_data = 0;
  int nrotbl[kMaxVegCats] = {0};
  real snuptbl[kMaxVegCats] = {0};
  real rstbl[kMaxVegCats] = {0};
  real rgltbl[kMaxVegCats] = {0};
  real hstbl[kMaxVegCats] = {0};
  real shdtbl[kMaxVegCats] = {0};
  real maxalb[kMaxVegCats] = {0};
  real emissmintbl[kMaxVegCats] = {0};
  real emissmaxtbl[kMaxVegCats] = {0};
  real laimintbl[kMaxVegCats] = {0};
  real laimaxtbl[kMaxVegCats] = {0};
  real z0mintbl[kMaxVegCats] = {0};
  real z0maxtbl[kMaxVegCats] = {0};
  real albedomintbl[kMaxVegCats] = {0};
  real albedomaxtbl[kMaxVegCats] = {0};
  real ztopvtbl[kMaxVegCats] = {0};
  real zbotvtbl[kMaxVegCats] = {0};
};

bool read_noah_tables(const std::string& table_dir,
                      const std::string& veg_dataset,
                      const std::string& soil_dataset,
                      NoahTables& tables);

}  // namespace cppwrf::phys::noah
