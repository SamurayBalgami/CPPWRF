#pragma once

#include <string>
#include <vector>

#include "core/field.hpp"
#include "core/grid.hpp"
#include "dyn/state.hpp"
#include "phys/cu_kf/kf_tables.hpp"
#include "phys/driver/landuse_tables.hpp"
#include "phys/lsm_noah/noah.hpp"
#include "phys/lsm_noah/noah_tables.hpp"
#include "phys/mp_wsm6/wsm6.hpp"
#include "phys/sf_sfclayrev/sfclayrev.hpp"

namespace cppwrf::phys::driver {

// Container holding all static physics parameters and tables
struct PhysicsStaticData {
  LanduseTables landuse;
  noah::NoahTables noah_tables;
  sfclayrev::SFClayRevTables<real> sfclay_tables;
  cu_kf::KFTables<real> kf_tables;
  wsm6::Wsm6Params<real> wsm6_params;
  bool initialized = false;
};

// Initializes all physics parameter tables from files in table_dir
bool init_physics_tables(const std::string& table_dir,
                         const std::string& lu_dataset,
                         PhysicsStaticData& data);

// Initializes 2D surface characteristics from LU_INDEX and tables (WRF LANDUSE_INIT)
void init_surface_properties(const Field2& lu_index,
                             const Field2& xice,
                             int julday,
                             real cen_lat,
                             int is_water,
                             int is_ice,
                             const LanduseTables& tables,
                             Field2& albedo,
                             Field2& albbck,
                             Field2& thc,
                             Field2& z0,
                             Field2& znt,
                             Field2& embck,
                             Field2& emiss,
                             Field2& mavail,
                             Field2& xland);

// Initializes Noah LSM soil properties: SNOALB from IVGTYP, and SH2O (liquid soil moisture)
// from TSLB and SMOIS using FRH2O for frozen soil (WRF lsminit)
void init_noah_soil(const Field2& ivgtyp,
                    const Field2& isltyp,
                    const Field3& tslb,
                    const Field3& smois,
                    const noah::NoahTables& tables,
                    int num_soil_layers,
                    Field2& snoalb,
                    Field3& sh2o);

}  // namespace cppwrf::phys::driver
