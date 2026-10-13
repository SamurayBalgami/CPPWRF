#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "core/field.hpp"
#include "core/grid.hpp"
#include "phys/driver/landuse_tables.hpp"
#include "phys/driver/physics_init.hpp"

using namespace cppwrf;
using namespace cppwrf::phys::driver;

int main(int argc, char** argv) {
  std::string tbl_dir = "/home/ostriquetrum/wrf_real/wrf_run";
  std::string bin_file = "wrfinput_extracted.bin";
  if (argc >= 2) tbl_dir = argv[1];
  if (argc >= 3) bin_file = argv[2];

  std::cout << "======================================================================\n";
  std::cout << " Real Data Physics Initialization Smoke Test (R0 wrfinput)\n";
  std::cout << " Table dir: " << tbl_dir << "\n";
  std::cout << " Input bin: " << bin_file << "\n";
  std::cout << "======================================================================\n";

  // 1. Read physics tables
  PhysicsStaticData static_data;
  if (!init_physics_tables(tbl_dir, "MODIFIED_IGBP_MODIS_NOAH", static_data)) {
    std::cerr << "FAILED: init_physics_tables failed.\n";
    return 1;
  }
  std::cout << "Physics tables loaded successfully.\n";

  // 2. Open extracted wrfinput binary
  std::ifstream f(bin_file, std::ios::binary);
  if (!f.is_open()) {
    std::cerr << "FAILED: could not open " << bin_file << "\n";
    return 1;
  }

  int hdr[4] = {0};
  f.read(reinterpret_cast<char*>(hdr), 4 * sizeof(int));
  int nx = hdr[0];
  int ny = hdr[1];
  int num_soil_layers = hdr[2];
  int julday = hdr[3];

  double cen_lat_d = 0;
  f.read(reinterpret_cast<char*>(&cen_lat_d), sizeof(double));
  real cen_lat = static_cast<real>(cen_lat_d);

  std::cout << "wrfinput domain: nx=" << nx << ", ny=" << ny
            << ", num_soil_layers=" << num_soil_layers
            << ", julday=" << julday << ", cen_lat=" << cen_lat << "\n";

  if (nx <= 0 || ny <= 0 || num_soil_layers <= 0) {
    std::cerr << "FAILED: invalid domain dimensions.\n";
    return 1;
  }

  std::size_t n2d = std::size_t(nx) * ny;
  std::size_t n3d_soil = std::size_t(nx) * ny * num_soil_layers;

  std::vector<double> buf_lu(n2d), buf_ivg(n2d), buf_isl(n2d), buf_xland(n2d);
  std::vector<double> buf_seaice(n2d), buf_tsk(n2d);
  std::vector<double> buf_tslb(n3d_soil), buf_smois(n3d_soil);
  std::vector<double> buf_albbck_ref(n2d), buf_snoalb_ref(n2d);

  f.read(reinterpret_cast<char*>(buf_lu.data()), n2d * sizeof(double));
  f.read(reinterpret_cast<char*>(buf_ivg.data()), n2d * sizeof(double));
  f.read(reinterpret_cast<char*>(buf_isl.data()), n2d * sizeof(double));
  f.read(reinterpret_cast<char*>(buf_xland.data()), n2d * sizeof(double));
  f.read(reinterpret_cast<char*>(buf_seaice.data()), n2d * sizeof(double));
  f.read(reinterpret_cast<char*>(buf_tsk.data()), n2d * sizeof(double));
  f.read(reinterpret_cast<char*>(buf_tslb.data()), n3d_soil * sizeof(double));
  f.read(reinterpret_cast<char*>(buf_smois.data()), n3d_soil * sizeof(double));
  f.read(reinterpret_cast<char*>(buf_albbck_ref.data()), n2d * sizeof(double));
  f.read(reinterpret_cast<char*>(buf_snoalb_ref.data()), n2d * sizeof(double));

  // Create CPPWRF Grid and Fields
  Grid g;
  g.nx = nx;
  g.ny = ny;
  g.nz = num_soil_layers;

  Field2 lu_index(g, Stagger::M, "lu_index");
  Field2 ivgtyp(g, Stagger::M, "ivgtyp");
  Field2 isltyp(g, Stagger::M, "isltyp");
  Field2 seaice(g, Stagger::M, "seaice");
  Field2 tsk(g, Stagger::M, "tsk");
  Field3 tslb(g, Stagger::M, "tslb");
  Field3 smois(g, Stagger::M, "smois");

  Field2 albedo(g, Stagger::M, "albedo");
  Field2 albbck(g, Stagger::M, "albbck");
  Field2 thc(g, Stagger::M, "thc");
  Field2 z0(g, Stagger::M, "z0");
  Field2 znt(g, Stagger::M, "znt");
  Field2 embck(g, Stagger::M, "embck");
  Field2 emiss(g, Stagger::M, "emiss");
  Field2 mavail(g, Stagger::M, "mavail");
  Field2 xland(g, Stagger::M, "xland");
  Field2 snoalb(g, Stagger::M, "snoalb");
  Field3 sh2o(g, Stagger::M, "sh2o");

  auto v_lu = lu_index.view();
  auto v_ivg = ivgtyp.view();
  auto v_isl = isltyp.view();
  auto v_seaice = seaice.view();
  auto v_tsk = tsk.view();
  auto v_tslb = tslb.view();
  auto v_smois = smois.view();

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      std::size_t idx2 = std::size_t(j) * nx + i;
      v_lu(i, j) = static_cast<real>(buf_lu[idx2]);
      v_ivg(i, j) = static_cast<real>(buf_ivg[idx2]);
      v_isl(i, j) = static_cast<real>(buf_isl[idx2]);
      v_seaice(i, j) = static_cast<real>(buf_seaice[idx2]);
      v_tsk(i, j) = static_cast<real>(buf_tsk[idx2]);

      for (int k = 0; k < num_soil_layers; ++k) {
        std::size_t idx3 = std::size_t(k) * (nx * ny) + std::size_t(j) * nx + i;
        v_tslb(i, j, k) = static_cast<real>(buf_tslb[idx3]);
        v_smois(i, j, k) = static_cast<real>(buf_smois[idx3]);
      }
    }
  }

  // Run surface initialisation
  std::cout << "Executing init_surface_properties...\n";
  init_surface_properties(lu_index, seaice, julday, cen_lat, 17, 15, static_data.landuse,
                          albedo, albbck, thc, z0, znt, embck, emiss, mavail, xland);

  // Run Noah soil initialisation
  std::cout << "Executing init_noah_soil...\n";
  init_noah_soil(ivgtyp, isltyp, tslb, smois, static_data.noah_tables, num_soil_layers,
                 snoalb, sh2o);

  // Statistics and Validation
  auto v_alb = albedo.view();
  auto v_thc = thc.view();
  auto v_z0 = z0.view();
  auto v_emiss = emiss.view();
  auto v_mavail = mavail.view();
  auto v_xland = xland.view();
  auto v_snoalb = snoalb.view();
  auto v_sh2o = sh2o.view();

  real min_alb = 1e9, max_alb = -1e9;
  real min_z0 = 1e9, max_z0 = -1e9;
  real min_emiss = 1e9, max_emiss = -1e9;
  real min_thc = 1e9, max_thc = -1e9;
  real min_mav = 1e9, max_mav = -1e9;
  real min_sno = 1e9, max_sno = -1e9;
  int n_water = 0, n_land = 0;

  bool has_nan = false;
  int sh2o_violations = 0;

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      real a = v_alb(i, j);
      real z = v_z0(i, j);
      real e = v_emiss(i, j);
      real t = v_thc(i, j);
      real m = v_mavail(i, j);
      real xl = v_xland(i, j);
      real sn = v_snoalb(i, j);

      if (std::isnan(a) || std::isnan(z) || std::isnan(e) || std::isnan(t) ||
          std::isnan(m) || std::isnan(xl) || std::isnan(sn)) {
        has_nan = true;
      }

      min_alb = std::min(min_alb, a); max_alb = std::max(max_alb, a);
      min_z0 = std::min(min_z0, z); max_z0 = std::max(max_z0, z);
      min_emiss = std::min(min_emiss, e); max_emiss = std::max(max_emiss, e);
      min_thc = std::min(min_thc, t); max_thc = std::max(max_thc, t);
      min_mav = std::min(min_mav, m); max_mav = std::max(max_mav, m);
      min_sno = std::min(min_sno, sn); max_sno = std::max(max_sno, sn);

      if (xl == real(2.0)) n_water++;
      else if (xl == real(1.0)) n_land++;

      for (int k = 0; k < num_soil_layers; ++k) {
        real sh = v_sh2o(i, j, k);
        real sm = v_smois(i, j, k);
        if (std::isnan(sh)) has_nan = true;
        if (sh > sm + real(1e-6)) sh2o_violations++;
      }
    }
  }

  std::cout << "\n======================================================================\n";
  std::cout << " Physics Initialization Diagnostic Report:\n";
  std::cout << "----------------------------------------------------------------------\n";
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "  ALBEDO:    [" << min_alb << ", " << max_alb << "]\n";
  std::cout << "  EMISS:     [" << min_emiss << ", " << max_emiss << "]\n";
  std::cout << "  Z0 (m):    [" << min_z0 << ", " << max_z0 << "]\n";
  std::cout << "  THC:       [" << min_thc << ", " << max_thc << "]\n";
  std::cout << "  MAVAIL:    [" << min_mav << ", " << max_mav << "]\n";
  std::cout << "  SNOALB:    [" << min_sno << ", " << max_sno << "]\n";
  std::cout << "  XLAND:     Land points = " << n_land << ", Water points = " << n_water
            << " (Total = " << (n_land + n_water) << " / " << (nx * ny) << ")\n";
  std::cout << "  SH2O <= SMOIS violations: " << sh2o_violations << "\n";
  std::cout << "======================================================================\n";

  if (has_nan) {
    std::cerr << "FAILED: NaN detected in physics initialization fields!\n";
    return 1;
  }
  if (sh2o_violations > 0) {
    std::cerr << "FAILED: Liquid water content SH2O exceeds total soil moisture SMOIS!\n";
    return 1;
  }
  if (min_alb < real(0.01) || max_alb > real(0.95)) {
    std::cerr << "FAILED: Albedo outside physical bounds [0.01, 0.95]!\n";
    return 1;
  }
  if (min_emiss < real(0.80) || max_emiss > real(1.00)) {
    std::cerr << "FAILED: Emissivity outside physical bounds [0.80, 1.00]!\n";
    return 1;
  }
  if (n_land + n_water != nx * ny) {
    std::cerr << "FAILED: XLAND points do not partition into 1.0 (land) and 2.0 (water)!\n";
    return 1;
  }

  std::cout << "REAL DATA SMOKE TEST PASSED: All fields physically valid and consistent.\n";
  return 0;
}
