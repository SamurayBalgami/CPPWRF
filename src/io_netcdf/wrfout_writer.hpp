#pragma once

#include <memory>
#include <string>
#include <vector>

#include "io_netcdf/netcdf_common.hpp"

namespace cppwrf::io_netcdf {

class WrfoutWriter {
 public:
  explicit WrfoutWriter(std::string out_dir, WrfoutMeta meta);
  ~WrfoutWriter();

  WrfoutWriter(const WrfoutWriter&) = delete;
  WrfoutWriter& operator=(const WrfoutWriter&) = delete;
  WrfoutWriter(WrfoutWriter&&) noexcept;
  WrfoutWriter& operator=(WrfoutWriter&&) noexcept;

  // Writes one time frame to a WRF-format NetCDF file.
  // If target_file is empty, defaults to "${out_dir}/wrfout_d01_${date_str}".
  void write_frame(long step, double time_sec, const std::vector<VarView>& vars,
                   const std::string& target_file = "");

  const std::string& last_written_file() const { return last_file_; }

 private:
  class Impl;
  std::unique_ptr<Impl> impl_;
  std::string last_file_;
};

}  // namespace cppwrf::io_netcdf
