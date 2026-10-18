#pragma once
// wrfout_d0N history files as WRF v4.8.0 writes them: the history stream of the WRF I/O schema
// (data/wrf_schema) for this configuration, with WRF's dimensions, variables, attributes and global
// attributes. Values come from the model state (by WRF name); a variable the model does not carry keeps
// its wrfinput value (WRF reads it from the input stream and leaves it unchanged) or WRF's initial 0.
#ifdef CPPWRF_HAVE_NETCDF
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "cases/model.hpp"
#include "io/config.hpp"

namespace cppwrf::wrf_io {
class StreamFile;
class WrfSchema;
class WrfConfig;
class Namelist;
}  // namespace cppwrf::wrf_io

namespace cppwrf::cases {

class WrfHistory {
 public:
  // fields: model fields written under their (WRF lower-case) names
  WrfHistory(const io::Config& cfg, const Model& m, std::string out_dir, std::string domain_tag,
             std::string start_date, std::vector<const Field2*> fields2, std::vector<const Field3*> fields3);
  ~WrfHistory();
  // one history frame at model time time_sec after the start (step = WRF's itimestep)
  void write(long step, double time_sec);
  const std::string& last_file() const { return last_file_; }

 private:
  bool value(const std::string& name, const std::string& varname, const std::vector<std::size_t>& dims,
             std::vector<double>& out);
  bool from_wrfinput(const std::string& name, std::size_t n, std::vector<double>& out);

  const Model& m_;
  std::string out_dir_, domain_tag_, start_date_, wrfinput_, outname_;
  std::vector<const Field2*> f2_;
  std::vector<const Field3*> f3_;
  std::unique_ptr<wrf_io::Namelist> nl_;
  std::unique_ptr<wrf_io::WrfSchema> schema_;
  std::unique_ptr<wrf_io::WrfConfig> wc_;
  std::unique_ptr<wrf_io::StreamFile> file_;
  int frames_per_outfile_ = 10;
  int wrfinput_id_ = -1;
  std::vector<double> p_hyd0_;  // P_HYD of the initial state (start_em)
  long step_ = 0;
  double time_sec_ = 0;
  std::string last_file_;
};

}  // namespace cppwrf::cases
#endif
