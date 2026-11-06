// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once

#include <string>
#include <vector>

#include "io_netcdf/netcdf_common.hpp"

namespace cppwrf::io_netcdf {

// the format of CPPWRF's own restart file (the global attribute CPPWRF_RESTART_VERSION): raise it when the file's
// contents change; a file of another version is refused (WRF's wrfrst files are the cross-version restart)
inline constexpr int kRestartVersion = 1;

struct RestartMeta {
  int version = kRestartVersion;  // read: 0 when the file has none (written before the version was recorded)
  std::string commit;             // CPPWRF_COMMIT: the commit of the sources that wrote it
  std::string options;            // CPPWRF_OPTIONS: the run's options, "section.key = value" lines
  long step = 0;
  double time_sec = 0.0;
  int nx = 0, ny = 0, nz = 0;
  double dx = 0.0, dy = 0.0, dt = 0.0;
  std::string case_name;
  std::string date_str = "0001-01-01_00:00:00";
};

// CPPWRF's own restart file (restart_NNNNNN.nc; WRF's wrfrst_d0N files are WrfHistory::write_restart):
// the given fields (the caller's list of state and carried physics fields) at their exact precision
// (NC_FLOAT for F32, NC_DOUBLE for F64).
void write_restart(const std::string& path, const RestartMeta& meta,
                   const std::vector<VarView>& vars);

// Reads such a file: fills meta and every given field into its pre-allocated buffer; a field missing from
// the file or of another size is an error.
void read_restart(const std::string& path, RestartMeta& meta,
                  std::vector<VarView>& vars);

// The global attributes alone (version, commit, options, step, time, grid), to check a file before reading it.
RestartMeta read_restart_meta(const std::string& path);

}  // namespace cppwrf::io_netcdf
