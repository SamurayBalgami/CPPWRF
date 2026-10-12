#pragma once

#include <string>
#include <vector>

#include "io_netcdf/netcdf_common.hpp"

namespace cppwrf::io_netcdf {

struct RestartMeta {
  long step = 0;
  double time_sec = 0.0;
  int nx = 0, ny = 0, nz = 0;
  double dx = 0.0, dy = 0.0, dt = 0.0;
  std::string case_name;
  std::string date_str = "0001-01-01_00:00:00";
};

// Writes a restart file containing all prognostic fields (*_1 and *_2, plus scalars)
// to a NetCDF file with exact precision (NC_FLOAT for F32, NC_DOUBLE for F64).
void write_restart(const std::string& path, const RestartMeta& meta,
                   const std::vector<VarView>& vars);

// Reads a restart file. Fills meta and populates data into the provided vars buffers.
// The provided vars must have their data pointers pre-allocated to receive the data.
void read_restart(const std::string& path, RestartMeta& meta,
                  std::vector<VarView>& vars);

}  // namespace cppwrf::io_netcdf
