#pragma once
// WRF phys/module_ra_clWRF_support.F read_CAMgases (ghg_input = 1): yearly trace-gas table
// CAMtr_volume_mixing_ratio (read once), interpolated in time (year, julian day) at every call.
#include <string>
#include <vector>

namespace cppwrf::phys {

struct CamGases {
  std::vector<int> yr;
  std::vector<double> jul, co2, n2o, ch4, cfc11, cfc12;
  // read the table (WRF reads it in phy_init); throws if missing
  void load(const std::string& path);
  // read_CAMgases(yr, julian, .false., "RRTMG", ...): volume mixing ratios
  void interpolate(int year, double julian, double& co2v, double& n2ov, double& ch4v, double& cfc11v,
                   double& cfc12v) const;
};

}  // namespace cppwrf::phys
