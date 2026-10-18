#pragma once
// Global attributes of WRF's share/output_wrf.F for the input, boundary, history and restart streams,
// in WRF's order.
#include <map>
#include <string>

#include "wrf_io/nc_out.hpp"
#include "wrf_io/wrf_schema.hpp"

namespace cppwrf::wrf_io {

enum class Stream { input, boundary, history, restart };

// Domain values output_wrf takes from the grid rather than the namelist.
struct DomainAtts {
  std::string title;                  // TITLE, e.g. " OUTPUT FROM WRF V4.8.0 MODEL"
  std::string start_date;             // START_DATE (this domain's start time)
  std::string simulation_start_date;  // SIMULATION_START_DATE of history/restart files
  int ide = 0, jde = 0, kde = 0;      // staggered grid dimensions
  int grid_id = 1, parent_id = 0, i_parent_start = 1, j_parent_start = 1, parent_grid_ratio = 1;
  double dt = 0;                      // grid%dt
  std::map<std::string, double> proj; // CEN_LAT, CEN_LON, TRUELAT1, TRUELAT2, MOAD_CEN_LAT, STAND_LON, POLE_LAT, POLE_LON
  int map_proj = 0;
  std::string mminlu;
  int num_land_cat = 21, iswater = 17, islake = -1, isice = 15, isurban = 13, isoilwater = 14;
  bool ideal = false;                 // SIMULATION_INITIALIZATION_TYPE
  int real_bytes = 4;                 // WRF's REAL kind (4, or 8 for -r8)
};

void write_global_atts(NcOut& nc, const WrfConfig& c, const DomainAtts& d, Stream s);

// geth_julgmt for a date yyyy-mm-dd_hh:mm:ss
void julgmt(const std::string& date, int& julyr, int& julday, double& gmt);

}  // namespace cppwrf::wrf_io
