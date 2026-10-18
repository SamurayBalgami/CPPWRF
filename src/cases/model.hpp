#pragma once
// Model assembly from a config (host side): grid, case, vertical levels, base and initial
// state, scalars, numerical parameters and the RK3 integrator. Shared by main and tests.
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "cases/cases.hpp"
#include "dyn/lateral_bc.hpp"
#include "dyn/diffusion_driver.hpp"
#include "dyn/rk3.hpp"
#include "io/config.hpp"
#include "phys/driver/suite.hpp"

namespace cppwrf::cases {

struct Model {
  Grid g;
  CaseDef cd;
  std::vector<double> eta;
  std::vector<double> xlat;
  std::vector<double> xlong;
  std::unique_ptr<dyn::BaseState> base;
  std::unique_ptr<dyn::State> st;
  std::unique_ptr<dyn::Tendencies> tend;
  dyn::LargeStepParams lp;
  dyn::AcousticParams ap;
  dyn::ScalarParams sp;
  double run_time = 0;
  long nsteps = 0;
  int n_sound = 0;
  int diff_opt = 0;
  int km_opt = 1;
  int hybrid_opt = 0;
  int map_proj = 0;
  double cen_lat = 0.0, cen_lon = 0.0;
  double truelat1 = 0.0, truelat2 = 0.0;
  double moad_cen_lat = 0.0, stand_lon = 0.0;
  std::string start_date;
  std::unique_ptr<dyn::RK3> rk;  // null when nsteps == 0
  std::unique_ptr<dyn::LateralBC> bc;  // real data: specified lateral boundaries
  std::unique_ptr<dyn::DiffusionDriver> diff;
  std::unique_ptr<phys::driver::PhysicsSuite> suite;  // column physics (sfclay/LSM/PBL), null: off  // WRF diff_opt 1 / sixth-order diffusion (null: off)
  // Called before step n (1-based) of the time loop: boundary interval switching (real data).
  std::function<void(long)> before_step;
  static Model build(const io::Config& cfg);
};

// case = real: grid, base and initial state, boundaries from wrfinput/wrfbdy (src/cases/real_case.cpp)
void init_real_case(const io::Config& cfg, Model& m);
// Column physics (real data, P1 suite): true if any of phys.sf_sfclay_physics,
// sf_surface_physics, bl_pbl_physics is set; init_real_physics builds the suite, fills the
// initial surface state as WRF start_domain/phy_init and attaches it to m.rk (real_physics.cpp).
bool wants_column_physics(const io::Config& cfg);
void init_real_physics(const io::Config& cfg, Model& m);
// phy_init's landuse_init (LANDUSE.TBL for the season of julday) on the wrfinput LU_INDEX, SEAICE, SNOWC
struct LanduseInit {
  std::vector<double> albbck, albedo, z0, embck, emiss, mavail, xland, thc;
  int isn = 1;                  // season (1 summer, 2 winter)
  int lucats = 0, luseas = 0;   // LANDUSE.TBL categories and seasons
  std::vector<double> lu_state;  // WRF's LU_STATE: per category SCFX, then per season ALBD SLMO SFEM SFZ0 SFHC THERIN
};
LanduseInit landuse_init(const std::string& table_dir, const std::string& mminlu, int julday, double cen_lat,
                         int iswater, int isice, const std::vector<double>& lu, const std::vector<double>& xice,
                         const std::vector<double>& snowc);

// Comma-separated list (whitespace trimmed, empty entries dropped).
std::vector<std::string> split_list(const std::string& s);

}  // namespace cppwrf::cases
