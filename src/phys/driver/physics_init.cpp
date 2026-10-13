#include "phys/driver/physics_init.hpp"

#include <cmath>
#include <iostream>

#include "core/backend.hpp"
#include "core/halo.hpp"

namespace cppwrf::phys::driver {

namespace {

using V2 = View2<real>;
using C2 = View2<const real>;
using V3 = View3<real>;
using C3 = View3<const real>;

struct KInitSurface {
  C2 lu_index;
  C2 xice;
  V2 albedo;
  V2 albbck;
  V2 thc;
  V2 z0;
  V2 znt;
  V2 embck;
  V2 emiss;
  V2 mavail;
  V2 xland;
  int isn;
  int is_water;
  int is_ice;
  LanduseTables tables;

  CPPWRF_HD void operator()(int i, int j) const {
    int is = static_cast<int>(std::round(lu_index(i, j)));
    if (is <= 0 || is >= kMaxLanduseCats) is = is_water;

    const LanduseCategory& cat = tables.categories[isn][is];
    albbck(i, j) = cat.albd * real(0.01);
    albedo(i, j) = albbck(i, j);
    thc(i, j) = cat.therin * real(0.01);
    z0(i, j) = cat.sfz0 * real(0.01);
    znt(i, j) = z0(i, j);
    embck(i, j) = cat.sfem;
    emiss(i, j) = embck(i, j);
    mavail(i, j) = cat.slmo;
    xland(i, j) = (is == is_water) ? real(2.0) : real(1.0);

    if (xice(i, j) >= real(0.02)) {
      xland(i, j) = real(1.0);
      const LanduseCategory& ice = tables.categories[isn][is_ice];
      albbck(i, j) = ice.albd * real(0.01);
      albedo(i, j) = albbck(i, j);
      embck(i, j) = ice.sfem;
      emiss(i, j) = embck(i, j);
      thc(i, j) = ice.therin * real(0.01);
      z0(i, j) = ice.sfz0 * real(0.01);
      znt(i, j) = z0(i, j);
      mavail(i, j) = ice.slmo;
    }
  }
};

struct KInitNoahSoil {
  C2 ivgtyp;
  C2 isltyp;
  C3 tslb;
  C3 smois;
  V2 snoalb;
  V3 sh2o;
  noah::NoahTables tables;
  int num_soil_layers;

  CPPWRF_HD void operator()(int i, int j) const {
    int iv = static_cast<int>(std::round(ivgtyp(i, j)));
    if (iv >= 1 && iv <= 50) {
      snoalb(i, j) = tables.maxalb[iv - 1] * real(0.01);
    } else {
      snoalb(i, j) = real(0.75);
    }

    int is = static_cast<int>(std::round(isltyp(i, j)));
    real bx = (is >= 1 && is <= tables.slcats) ? tables.bb[is - 1] : real(0);
    real smcmax = (is >= 1 && is <= tables.slcats) ? tables.maxsmc[is - 1] : real(0);
    real psisat = (is >= 1 && is <= tables.slcats) ? tables.satpsi[is - 1] : real(0);

    for (int k = 0; k < num_soil_layers; ++k) {
      real t = tslb(i, j, k);
      real sm = smois(i, j, k);
      if (t < real(273.149) && bx > real(0) && smcmax > real(0) && psisat > real(0)) {
        real free_w = real(0);
        noah::NoahKernel::frh2o(free_w, t, sm, sm, smcmax, bx, psisat);
        sh2o(i, j, k) = free_w;
      } else {
        sh2o(i, j, k) = sm;
      }
    }
  }
};

}  // namespace

bool init_physics_tables(const std::string& table_dir,
                         const std::string& lu_dataset,
                         PhysicsStaticData& data) {
  // 1. LANDUSE.TBL
  std::string lu_file = table_dir + "/LANDUSE.TBL";
  if (!read_landuse_table(lu_file, lu_dataset, data.landuse)) {
    std::cerr << "Error reading LANDUSE.TBL from " << lu_file << "\n";
    return false;
  }

  // 2. Noah tables (VEGPARM.TBL, SOILPARM.TBL, GENPARM.TBL)
  if (!noah::read_noah_tables(table_dir, lu_dataset, "STAS", data.noah_tables)) {
    std::cerr << "Error reading Noah tables from " << table_dir << "\n";
    return false;
  }

  // 3. sfclayrev tables
  sfclayrev::init_sfclayrev_tables(data.sfclay_tables);

  // 4. KF tables
  cu_kf::init_kf_tables(data.kf_tables);

  // 5. WSM6 parameters
  data.wsm6_params.init(real(1.28), real(1000.0), real(100.0), real(4190.0), real(1846.4), 0);

  data.initialized = true;
  return true;
}

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
                             Field2& xland) {
  int isn = (julday < 105 || julday > 288) ? 1 : 0;
  if (cen_lat < real(0)) isn = 1 - isn;

  const auto& L = lu_index.layout();
  column_for(Range2D{0, L.nx, 0, L.ny},
             KInitSurface{lu_index.cview(), xice.cview(), albedo.view(), albbck.view(),
                          thc.view(), z0.view(), znt.view(), embck.view(), emiss.view(),
                          mavail.view(), xland.view(), isn, is_water, is_ice, tables});
}

void init_noah_soil(const Field2& ivgtyp,
                    const Field2& isltyp,
                    const Field3& tslb,
                    const Field3& smois,
                    const noah::NoahTables& tables,
                    int num_soil_layers,
                    Field2& snoalb,
                    Field3& sh2o) {
  const auto& L = ivgtyp.layout();
  column_for(Range2D{0, L.nx, 0, L.ny},
             KInitNoahSoil{ivgtyp.cview(), isltyp.cview(), tslb.cview(), smois.cview(),
                           snoalb.view(), sh2o.view(), tables, num_soil_layers});
}

}  // namespace cppwrf::phys::driver
