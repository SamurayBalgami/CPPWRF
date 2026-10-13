#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "phys/driver/landuse_tables.hpp"
#include "phys/driver/physics_init.hpp"

using namespace cppwrf;
using namespace cppwrf::phys::driver;

int main(int argc, char** argv) {
  std::string tbl_dir = ".";
  std::string ref_bin = "tables_ref.bin";
  if (argc >= 2) tbl_dir = argv[1];
  if (argc >= 3) ref_bin = argv[2];

  std::cout << "Testing physics table readers against Fortran reference...\n";
  std::cout << "Table dir: " << tbl_dir << "\n";
  std::cout << "Ref bin:   " << ref_bin << "\n";

  // 1. Read tables via C++
  PhysicsStaticData data;
  if (!init_physics_tables(tbl_dir, "MODIFIED_IGBP_MODIS_NOAH", data)) {
    std::cerr << "FAILED: init_physics_tables failed.\n";
    return 1;
  }

  // 2. Read Fortran binary reference
  std::ifstream rf(ref_bin, std::ios::binary);
  if (!rf.is_open()) {
    std::cerr << "FAILED: could not open " << ref_bin << "\n";
    return 1;
  }

  int ref_lucats = 0, ref_luseas = 0;
  rf.read(reinterpret_cast<char*>(&ref_lucats), sizeof(int));
  rf.read(reinterpret_cast<char*>(&ref_luseas), sizeof(int));

  if (ref_lucats != data.landuse.n_cats || ref_luseas != data.landuse.n_seasons) {
    std::cerr << "Mismatch in LANDUSE cats/seasons: C++(" << data.landuse.n_cats << ", "
              << data.landuse.n_seasons << ") vs Fortran(" << ref_lucats << ", " << ref_luseas << ")\n";
    return 1;
  }

  // Fortran arrays: albd(lucats, luseas), slmo, sfem, sfz0, therin, scfx(lucats), sfhc
  std::vector<double> ref_albd(ref_lucats * ref_luseas);
  std::vector<double> ref_slmo(ref_lucats * ref_luseas);
  std::vector<double> ref_sfem(ref_lucats * ref_luseas);
  std::vector<double> ref_sfz0(ref_lucats * ref_luseas);
  std::vector<double> ref_therin(ref_lucats * ref_luseas);
  std::vector<double> ref_scfx(ref_lucats * ref_luseas);
  std::vector<double> ref_sfhc(ref_lucats * ref_luseas);

  rf.read(reinterpret_cast<char*>(ref_albd.data()), ref_albd.size() * sizeof(double));
  rf.read(reinterpret_cast<char*>(ref_slmo.data()), ref_slmo.size() * sizeof(double));
  rf.read(reinterpret_cast<char*>(ref_sfem.data()), ref_sfem.size() * sizeof(double));
  rf.read(reinterpret_cast<char*>(ref_sfz0.data()), ref_sfz0.size() * sizeof(double));
  rf.read(reinterpret_cast<char*>(ref_therin.data()), ref_therin.size() * sizeof(double));
  rf.read(reinterpret_cast<char*>(ref_scfx.data()), ref_scfx.size() * sizeof(double));
  rf.read(reinterpret_cast<char*>(ref_sfhc.data()), ref_sfhc.size() * sizeof(double));

  double max_diff_lu = 0.0;
  for (int s = 0; s < ref_luseas; ++s) {
    for (int c = 1; c <= ref_lucats; ++c) {
      int idx = (c - 1) + s * ref_lucats;
      const auto& cat = data.landuse.categories[s][c];
      double d_alb = std::abs(double(cat.albd) - ref_albd[idx]);
      double d_slmo = std::abs(double(cat.slmo) - ref_slmo[idx]);
      double d_sfem = std::abs(double(cat.sfem) - ref_sfem[idx]);
      double d_sfz0 = std::abs(double(cat.sfz0) - ref_sfz0[idx]);
      double d_therin = std::abs(double(cat.therin) - ref_therin[idx]);
      double d_scfx = std::abs(double(cat.scfx) - ref_scfx[idx]);
      double d_sfhc = std::abs(double(cat.sfhc) - ref_sfhc[idx]);
      if (d_alb > 1e-4) std::cout << "ALBD diff at s=" << s << " c=" << c << ": cpp=" << cat.albd << " ref=" << ref_albd[idx] << "\n";
      if (d_slmo > 1e-4) std::cout << "SLMO diff at s=" << s << " c=" << c << ": cpp=" << cat.slmo << " ref=" << ref_slmo[idx] << "\n";
      if (d_sfem > 1e-4) std::cout << "SFEM diff at s=" << s << " c=" << c << ": cpp=" << cat.sfem << " ref=" << ref_sfem[idx] << "\n";
      if (d_sfz0 > 1e-4) std::cout << "SFZ0 diff at s=" << s << " c=" << c << ": cpp=" << cat.sfz0 << " ref=" << ref_sfz0[idx] << "\n";
      if (d_therin > 1e-4) std::cout << "THERIN diff at s=" << s << " c=" << c << ": cpp=" << cat.therin << " ref=" << ref_therin[idx] << "\n";
      if (d_scfx > 1e-4) std::cout << "SCFX diff at s=" << s << " c=" << c << ": cpp=" << cat.scfx << " ref=" << ref_scfx[idx] << "\n";
      if (d_sfhc > 1e-4) std::cout << "SFHC diff at s=" << s << " c=" << c << ": cpp=" << cat.sfhc << " ref=" << ref_sfhc[idx] << "\n";
      max_diff_lu = std::max(max_diff_lu, d_alb);
      max_diff_lu = std::max(max_diff_lu, d_slmo);
      max_diff_lu = std::max(max_diff_lu, d_sfem);
      max_diff_lu = std::max(max_diff_lu, d_sfz0);
      max_diff_lu = std::max(max_diff_lu, d_therin);
      max_diff_lu = std::max(max_diff_lu, d_scfx);
      max_diff_lu = std::max(max_diff_lu, d_sfhc);
    }
  }

  std::cout << "LANDUSE.TBL max diff vs Fortran: " << max_diff_lu << "\n";
  if (max_diff_lu > 1e-12) {
    std::cerr << "FAILED: LANDUSE.TBL diff too large!\n";
    return 1;
  }

  // Noah GENPARM
  int ref_slpcats = 0;
  rf.read(reinterpret_cast<char*>(&ref_slpcats), sizeof(int));
  std::vector<double> ref_slope(30);
  rf.read(reinterpret_cast<char*>(ref_slope.data()), 30 * sizeof(double));
  double sbeta, fxexp, csoil, salp, refdk, refkdt, frzk, zbot, czil, smlow, smhigh, lvcoef;
  rf.read(reinterpret_cast<char*>(&sbeta), sizeof(double));
  rf.read(reinterpret_cast<char*>(&fxexp), sizeof(double));
  rf.read(reinterpret_cast<char*>(&csoil), sizeof(double));
  rf.read(reinterpret_cast<char*>(&salp), sizeof(double));
  rf.read(reinterpret_cast<char*>(&refdk), sizeof(double));
  rf.read(reinterpret_cast<char*>(&refkdt), sizeof(double));
  rf.read(reinterpret_cast<char*>(&frzk), sizeof(double));
  rf.read(reinterpret_cast<char*>(&zbot), sizeof(double));
  rf.read(reinterpret_cast<char*>(&czil), sizeof(double));
  rf.read(reinterpret_cast<char*>(&smlow), sizeof(double));
  rf.read(reinterpret_cast<char*>(&smhigh), sizeof(double));
  rf.read(reinterpret_cast<char*>(&lvcoef), sizeof(double));

  double max_diff_gen = 0.0;
  max_diff_gen = std::max(max_diff_gen, std::abs(double(data.noah_tables.sbeta_data) - sbeta));
  max_diff_gen = std::max(max_diff_gen, std::abs(double(data.noah_tables.fxexp_data) - fxexp));
  max_diff_gen = std::max(max_diff_gen, std::abs(double(data.noah_tables.csoil_data) - csoil));
  max_diff_gen = std::max(max_diff_gen, std::abs(double(data.noah_tables.salp_data) - salp));
  max_diff_gen = std::max(max_diff_gen, std::abs(double(data.noah_tables.refdk_data) - refdk));
  max_diff_gen = std::max(max_diff_gen, std::abs(double(data.noah_tables.refkdt_data) - refkdt));
  max_diff_gen = std::max(max_diff_gen, std::abs(double(data.noah_tables.frzk_data) - frzk));
  max_diff_gen = std::max(max_diff_gen, std::abs(double(data.noah_tables.zbot_data) - zbot));
  max_diff_gen = std::max(max_diff_gen, std::abs(double(data.noah_tables.czil_data) - czil));
  max_diff_gen = std::max(max_diff_gen, std::abs(double(data.noah_tables.smlow_data) - smlow));
  max_diff_gen = std::max(max_diff_gen, std::abs(double(data.noah_tables.smhigh_data) - smhigh));
  max_diff_gen = std::max(max_diff_gen, std::abs(double(data.noah_tables.lvcoef_data) - lvcoef));

  std::cout << "GENPARM.TBL max diff vs Fortran: " << max_diff_gen << "\n";
  if (max_diff_gen > 1e-12) {
    std::cerr << "FAILED: GENPARM.TBL diff too large!\n";
    return 1;
  }

  // Noah SOILPARM
  int ref_slcats = 0;
  rf.read(reinterpret_cast<char*>(&ref_slcats), sizeof(int));
  std::vector<double> bb(30), drysmc(30), f11(30), maxsmc(30), refsmc(30), satpsi(30), satdk(30), satdw(30), wltsmc(30), qtz(30);
  rf.read(reinterpret_cast<char*>(bb.data()), 30 * sizeof(double));
  rf.read(reinterpret_cast<char*>(drysmc.data()), 30 * sizeof(double));
  rf.read(reinterpret_cast<char*>(f11.data()), 30 * sizeof(double));
  rf.read(reinterpret_cast<char*>(maxsmc.data()), 30 * sizeof(double));
  rf.read(reinterpret_cast<char*>(refsmc.data()), 30 * sizeof(double));
  rf.read(reinterpret_cast<char*>(satpsi.data()), 30 * sizeof(double));
  rf.read(reinterpret_cast<char*>(satdk.data()), 30 * sizeof(double));
  rf.read(reinterpret_cast<char*>(satdw.data()), 30 * sizeof(double));
  rf.read(reinterpret_cast<char*>(wltsmc.data()), 30 * sizeof(double));
  rf.read(reinterpret_cast<char*>(qtz.data()), 30 * sizeof(double));

  double max_diff_soil = 0.0;
  for (int i = 0; i < ref_slcats; ++i) {
    max_diff_soil = std::max(max_diff_soil, std::abs(double(data.noah_tables.bb[i]) - bb[i]));
    max_diff_soil = std::max(max_diff_soil, std::abs(double(data.noah_tables.drysmc[i]) - drysmc[i]));
    max_diff_soil = std::max(max_diff_soil, std::abs(double(data.noah_tables.f11[i]) - f11[i]));
    max_diff_soil = std::max(max_diff_soil, std::abs(double(data.noah_tables.maxsmc[i]) - maxsmc[i]));
    max_diff_soil = std::max(max_diff_soil, std::abs(double(data.noah_tables.refsmc[i]) - refsmc[i]));
    max_diff_soil = std::max(max_diff_soil, std::abs(double(data.noah_tables.satpsi[i]) - satpsi[i]));
    max_diff_soil = std::max(max_diff_soil, std::abs(double(data.noah_tables.satdk[i]) - satdk[i]));
    max_diff_soil = std::max(max_diff_soil, std::abs(double(data.noah_tables.satdw[i]) - satdw[i]));
    max_diff_soil = std::max(max_diff_soil, std::abs(double(data.noah_tables.wltsmc[i]) - wltsmc[i]));
    max_diff_soil = std::max(max_diff_soil, std::abs(double(data.noah_tables.qtz[i]) - qtz[i]));
  }

  std::cout << "SOILPARM.TBL max diff vs Fortran: " << max_diff_soil << "\n";
  if (max_diff_soil > 1e-12) {
    std::cerr << "FAILED: SOILPARM.TBL diff too large!\n";
    return 1;
  }

  // Noah VEGPARM
  int veg_lucats = 0, bare = 0, natural = 0;
  rf.read(reinterpret_cast<char*>(&veg_lucats), sizeof(int));
  rf.read(reinterpret_cast<char*>(&bare), sizeof(int));
  rf.read(reinterpret_cast<char*>(&natural), sizeof(int));

  double topt, cmcmax, cfactr, rsmax;
  rf.read(reinterpret_cast<char*>(&topt), sizeof(double));
  rf.read(reinterpret_cast<char*>(&cmcmax), sizeof(double));
  rf.read(reinterpret_cast<char*>(&cfactr), sizeof(double));
  rf.read(reinterpret_cast<char*>(&rsmax), sizeof(double));

  std::vector<int> nrotbl(50);
  std::vector<double> snuptbl(50), rstbl(50), rgltbl(50), hstbl(50), shdtbl(50), maxalb(50);
  std::vector<double> emissmintbl(50), emissmaxtbl(50), laimintbl(50), laimaxtbl(50);
  std::vector<double> z0mintbl(50), z0maxtbl(50), albedomintbl(50), albedomaxtbl(50);
  std::vector<double> ztopvtbl(50), zbotvtbl(50);

  rf.read(reinterpret_cast<char*>(nrotbl.data()), 50 * sizeof(int));
  rf.read(reinterpret_cast<char*>(snuptbl.data()), 50 * sizeof(double));
  rf.read(reinterpret_cast<char*>(rstbl.data()), 50 * sizeof(double));
  rf.read(reinterpret_cast<char*>(rgltbl.data()), 50 * sizeof(double));
  rf.read(reinterpret_cast<char*>(hstbl.data()), 50 * sizeof(double));
  rf.read(reinterpret_cast<char*>(shdtbl.data()), 50 * sizeof(double));
  rf.read(reinterpret_cast<char*>(maxalb.data()), 50 * sizeof(double));
  rf.read(reinterpret_cast<char*>(emissmintbl.data()), 50 * sizeof(double));
  rf.read(reinterpret_cast<char*>(emissmaxtbl.data()), 50 * sizeof(double));
  rf.read(reinterpret_cast<char*>(laimintbl.data()), 50 * sizeof(double));
  rf.read(reinterpret_cast<char*>(laimaxtbl.data()), 50 * sizeof(double));
  rf.read(reinterpret_cast<char*>(z0mintbl.data()), 50 * sizeof(double));
  rf.read(reinterpret_cast<char*>(z0maxtbl.data()), 50 * sizeof(double));
  rf.read(reinterpret_cast<char*>(albedomintbl.data()), 50 * sizeof(double));
  rf.read(reinterpret_cast<char*>(albedomaxtbl.data()), 50 * sizeof(double));
  rf.read(reinterpret_cast<char*>(ztopvtbl.data()), 50 * sizeof(double));
  rf.read(reinterpret_cast<char*>(zbotvtbl.data()), 50 * sizeof(double));

  double max_diff_veg = 0.0;
  max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.topt_data) - topt));
  max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.cmcmax_data) - cmcmax));
  max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.cfactr_data) - cfactr));
  max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.rsmax_data) - rsmax));

  for (int i = 0; i < veg_lucats; ++i) {
    max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.nrotbl[i]) - nrotbl[i]));
    max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.snuptbl[i]) - snuptbl[i]));
    max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.rstbl[i]) - rstbl[i]));
    max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.rgltbl[i]) - rgltbl[i]));
    max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.hstbl[i]) - hstbl[i]));
    max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.shdtbl[i]) - shdtbl[i]));
    max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.maxalb[i]) - maxalb[i]));
    max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.emissmintbl[i]) - emissmintbl[i]));
    max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.emissmaxtbl[i]) - emissmaxtbl[i]));
    max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.laimintbl[i]) - laimintbl[i]));
    max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.laimaxtbl[i]) - laimaxtbl[i]));
    max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.z0mintbl[i]) - z0mintbl[i]));
    max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.z0maxtbl[i]) - z0maxtbl[i]));
    max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.albedomintbl[i]) - albedomintbl[i]));
    max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.albedomaxtbl[i]) - albedomaxtbl[i]));
    max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.ztopvtbl[i]) - ztopvtbl[i]));
    max_diff_veg = std::max(max_diff_veg, std::abs(double(data.noah_tables.zbotvtbl[i]) - zbotvtbl[i]));
  }

  std::cout << "VEGPARM.TBL max diff vs Fortran: " << max_diff_veg << "\n";
  if (max_diff_veg > 1e-12) {
    std::cerr << "FAILED: VEGPARM.TBL diff too large!\n";
    return 1;
  }

  std::cout << "ALL PHYSICS TABLES MATCH FORTRAN REFERENCE EXACTLY (max diff <= 1e-12)!\n";
  return 0;
}
