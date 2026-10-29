#pragma once
// Global attributes of WRF's share/output_wrf.F for the input, boundary, history and restart streams,
// in WRF's order.
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "wrf_io/nc_out.hpp"
#include "wrf_io/wrf_schema.hpp"

namespace cppwrf::wrf_io {

enum class Stream { input, boundary, history, restart, auxinput };  // auxinput: e.g. wrflowinp (auxinput4)

// Domain values output_wrf takes from the grid rather than the namelist.
struct DomainAtts {
  std::string title;                  // TITLE, e.g. " OUTPUT FROM WRF V4.8.0 MODEL"
  std::string start_date;             // START_DATE (this domain's start time)
  std::string simulation_start_date;  // SIMULATION_START_DATE of history/restart files
  std::string gmt_date;               // GMT/JULYR/JULDAY: the domain's start (empty: start_date); a
                                      // resumed run keeps it while START_DATE is the restart time
  // config_flags%gmt/julyr/julday as WRF has them (wrf.exe: read from wrfinput; ideal: the namelist);
  // without them (real.exe) geth_julgmt of gmt_date
  bool have_julgmt = false;
  double gmt = 0;
  int julyr = 0, julday = 1;
  int ide = 0, jde = 0, kde = 0;      // staggered grid dimensions
  int grid_id = 1, parent_id = 0, i_parent_start = 1, j_parent_start = 1, parent_grid_ratio = 1;
  double dt = 0;                      // grid%dt
  // use_adaptive_time_step: ADAPT_DT_START / _MAX / _MIN of history files (grid%starting_time_step, ...)
  bool adaptive = false;
  double adapt_dt_start = 0, adapt_dt_max = 0, adapt_dt_min = 0;
  std::map<std::string, double> proj; // CEN_LAT, CEN_LON, TRUELAT1, TRUELAT2, MOAD_CEN_LAT, STAND_LON, POLE_LAT, POLE_LON
  int map_proj = 0;
  std::string mminlu;
  int num_land_cat = 21, iswater = 17, islake = -1, isice = 15, isurban = 13, isoilwater = 14;
  bool ideal = false;                 // SIMULATION_INITIALIZATION_TYPE
  int real_bytes = 4;                 // WRF's REAL kind (4, or 8 for -r8)
  // restart files: FLAG_RESTART, MAX_WRF_ALARMS and the alarm state (name, value) in WRF's order
  std::vector<std::pair<std::string, int>> restart_alarms;
};

void write_global_atts(NcOut& nc, const WrfConfig& c, const DomainAtts& d, Stream s);

// geth_julgmt (share/module_date_time.F) for a date yyyy-mm-dd_hh:mm:ss; as WRF's, February keeps 29 days
// once a leap-year date was seen
void julgmt(const std::string& date, int& julyr, int& julday, double& gmt);
// date + seconds (geth_newdate), yyyy-mm-dd_hh:mm:ss
std::string add_seconds(const std::string& date, long seconds);
// construct_filename2a (share/module_io_domain.F): a stream's file name from its namelist name (history_outname,
// rst_outname; quotes and trailing blanks removed): <domain> -> the 2-digit domain id, <date>, <year>, <month>,
// <day>, <hour> -> parts of date (nocolons is not ported)
std::string stream_filename(std::string basename, int domain_id, const std::string& date);

}  // namespace cppwrf::wrf_io
