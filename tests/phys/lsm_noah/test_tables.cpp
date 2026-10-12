#include <cmath>
#include <cstdio>
#include <iostream>
#include "phys/lsm_noah/noah_tables.hpp"

using namespace cppwrf::phys::noah;

int main() {
  NoahTables tables;
  if (!read_noah_tables(".", "USGS", "STAS", tables)) {
    std::cerr << "Failed to read Noah tables from .\n";
    return 1;
  }

  FILE* f = std::fopen("fortran/tables_ref.bin", "rb");
  if (!f) {
    std::cerr << "Cannot open fortran/tables_ref.bin\n";
    return 1;
  }

  int slpcats = 0;
  double slope_data[30] = {0};
  double sbeta, fxexp, csoil, salp, refdk, refkdt, frzk, zbot, czil, smlow, smhigh, lvcoef;

  std::fread(&slpcats, sizeof(int), 1, f);
  std::fread(slope_data, sizeof(double), 30, f);
  std::fread(&sbeta, sizeof(double), 1, f);
  std::fread(&fxexp, sizeof(double), 1, f);
  std::fread(&csoil, sizeof(double), 1, f);
  std::fread(&salp, sizeof(double), 1, f);
  std::fread(&refdk, sizeof(double), 1, f);
  std::fread(&refkdt, sizeof(double), 1, f);
  std::fread(&frzk, sizeof(double), 1, f);
  std::fread(&zbot, sizeof(double), 1, f);
  std::fread(&czil, sizeof(double), 1, f);
  std::fread(&smlow, sizeof(double), 1, f);
  std::fread(&smhigh, sizeof(double), 1, f);
  std::fread(&lvcoef, sizeof(double), 1, f);

  int slcats = 0;
  double bb[30], drysmc[30], f11[30], maxsmc[30], refsmc[30], satpsi[30], satdk[30], satdw[30], wltsmc[30], qtz[30];
  std::fread(&slcats, sizeof(int), 1, f);
  std::fread(bb, sizeof(double), 30, f);
  std::fread(drysmc, sizeof(double), 30, f);
  std::fread(f11, sizeof(double), 30, f);
  std::fread(maxsmc, sizeof(double), 30, f);
  std::fread(refsmc, sizeof(double), 30, f);
  std::fread(satpsi, sizeof(double), 30, f);
  std::fread(satdk, sizeof(double), 30, f);
  std::fread(satdw, sizeof(double), 30, f);
  std::fread(wltsmc, sizeof(double), 30, f);
  std::fread(qtz, sizeof(double), 30, f);

  int lucats = 0, bare = 0, natural = 0;
  double topt, cmcmax, cfactr, rsmax;
  int nrotbl[50];
  double snuptbl[50], rstbl[50], rgltbl[50], hstbl[50], shdtbl[50], maxalb[50];
  double emissmintbl[50], emissmaxtbl[50], laimintbl[50], laimaxtbl[50];
  double z0mintbl[50], z0maxtbl[50], albedomintbl[50], albedomaxtbl[50];
  double ztopvtbl[50], zbotvtbl[50];

  std::fread(&lucats, sizeof(int), 1, f);
  std::fread(&bare, sizeof(int), 1, f);
  std::fread(&natural, sizeof(int), 1, f);
  std::fread(&topt, sizeof(double), 1, f);
  std::fread(&cmcmax, sizeof(double), 1, f);
  std::fread(&cfactr, sizeof(double), 1, f);
  std::fread(&rsmax, sizeof(double), 1, f);
  std::fread(nrotbl, sizeof(int), 50, f);
  std::fread(snuptbl, sizeof(double), 50, f);
  std::fread(rstbl, sizeof(double), 50, f);
  std::fread(rgltbl, sizeof(double), 50, f);
  std::fread(hstbl, sizeof(double), 50, f);
  std::fread(shdtbl, sizeof(double), 50, f);
  std::fread(maxalb, sizeof(double), 50, f);
  std::fread(emissmintbl, sizeof(double), 50, f);
  std::fread(emissmaxtbl, sizeof(double), 50, f);
  std::fread(laimintbl, sizeof(double), 50, f);
  std::fread(laimaxtbl, sizeof(double), 50, f);
  std::fread(z0mintbl, sizeof(double), 50, f);
  std::fread(z0maxtbl, sizeof(double), 50, f);
  std::fread(albedomintbl, sizeof(double), 50, f);
  std::fread(albedomaxtbl, sizeof(double), 50, f);
  std::fread(ztopvtbl, sizeof(double), 50, f);
  std::fread(zbotvtbl, sizeof(double), 50, f);
  std::fclose(f);

  double max_diff = 0.0;
  auto check = [&](double a, double b, const char* name) {
    double diff = std::abs(a - b);
    if (diff > max_diff) max_diff = diff;
    if (diff > 1e-14) {
      std::cout << "Mismatch in " << name << ": C++=" << a << ", Fortran=" << b << "\n";
    }
  };

  check(tables.slpcats, slpcats, "slpcats");
  for (int i = 0; i < tables.slpcats; ++i) check(tables.slope_data[i], slope_data[i], "slope_data");
  check(tables.sbeta_data, sbeta, "sbeta");
  check(tables.csoil_data, csoil, "csoil");
  check(tables.lvcoef_data, lvcoef, "lvcoef");

  check(tables.slcats, slcats, "slcats");
  for (int i = 0; i < tables.slcats; ++i) {
    check(tables.bb[i], bb[i], "bb");
    check(tables.satdk[i], satdk[i], "satdk");
    check(tables.maxsmc[i], maxsmc[i], "maxsmc");
  }

  check(tables.lucats, lucats, "lucats");
  check(tables.bare, bare, "bare");
  check(tables.natural, natural, "natural");
  check(tables.topt_data, topt, "topt");
  for (int i = 0; i < tables.lucats; ++i) {
    check(tables.nrotbl[i], nrotbl[i], "nrotbl");
    check(tables.shdtbl[i], shdtbl[i], "shdtbl");
    check(tables.rstbl[i], rstbl[i], "rstbl");
  }

  std::cout << "Table comparison finished. Max diff = " << max_diff << "\n";
  return (max_diff > 1e-14) ? 1 : 0;
}
