#pragma once
// WRF phys/module_radiation_driver.F radiation_driver for ra_lw_physics = 4 / ra_sw_physics = 4
// (RRTMG) with o3input = 2 (CAM ozone), ghg_input = 1 (CAMtr_volume_mixing_ratio), aer_opt = 0,
// swint_opt = 0, icloud = 1, cldovrlp, idcor, slope_rad = topo_shading = 0, single domain.
// First version: the columns run on the host (OpenMP); in CUDA builds the inputs are downloaded
// and the outputs uploaded once per radiation call.
#include <memory>
#include <string>
#include <vector>

#include "core/field.hpp"
#include "core/grid.hpp"

namespace cppwrf::phys::driver {

struct RadiationOptions {
  int ra_lw_physics = 4, ra_sw_physics = 4;
  double radt = 0;          // minutes; stepra = max(nint(radt*60/dt), 1)
  int ra_call_offset = 0;
  int icloud = 1, cldovrlp = 2, idcor = 0, o3input = 2, ghg_input = 1, aer_opt = 0;
  int julyr = 0, julday = 1;  // start (wrfinput JULYR/JULDAY), gmt start hour
  double gmt = 0;
  double p_top = 0;
  std::string table_dir = ".";  // RRTMG_*_DATA tables dump, ozone*.formatted, CAMtr_volume_mixing_ratio
};

struct RadiationInputs {  // device fields (mass points) read by the radiation
  const Field3 *t_phy, *p_hyd, *p_hyd_w, *pi_phy, *rho, *dz8w, *t8w, *qv;
  const Field3* cldfra;  // nullptr: no clouds (cldfra = 0)
  const Field2 *tsk, *albedo, *emiss, *xland, *xice, *snow;
};

struct RadiationOutputs {
  Field3* rthraten;  // uncoupled theta tendency (K/s), LW + SW
  Field2 *glw, *gsw, *swdown, *coszen;
};

class RadiationDriver {
 public:
  RadiationDriver(const Grid& g, const RadiationOptions& o, const std::vector<double>& xlat,
                  const std::vector<double>& xlong);
  ~RadiationDriver();
  int stepra() const { return stepra_; }
  // radiation_driver for time step itimestep (1-based) at xtime minutes after the start and
  // julian (day of year - 1 + fraction); returns false when the radiation is not called
  bool run(long itimestep, double xtime, double julian, const RadiationInputs& in, const RadiationOutputs& out);

 private:
  struct Impl;
  std::unique_ptr<Impl> p_;
  Grid g_;
  RadiationOptions o_;
  int stepra_ = 1;
};

}  // namespace cppwrf::phys::driver
