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
#include "phys/driver/phy_prep.hpp"

namespace cppwrf::wrf_io {
struct DomainAtts;
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
  // wrfrst_d0N_<date>: the restart stream at model time time_sec (alarms from the intervals in seconds)
  void write_restart(long step, double time_sec, double history_interval_s, double restart_interval_s);
  const std::string& last_restart() const { return last_restart_; }
  // before step n of a frame that will be written: runs without column physics evaluate phy_prep as
  // first_rk_step_part1 does every step (P_HYD of the frame)
  void before_step(long n);
  // phy_prep of that step where first_rk_step_part1 runs it (RK stage 1, after rk_step_prep)
  void first_rk_prep();
  // after the solve of step n, before nest feedback changes the state: WRF's T (th_phy_m_t0) of the frame
  void after_solve(long n);
  // med_nest_feedback for the history fields the model does not carry (Registry u= functions)
  void feedback_from(WrfHistory& child, int i_parent_start, int j_parent_start, int ratio);
  const std::string& last_file() const { return last_file_; }
  // a run resumed from a restart file at time_sec: START_DATE is the restart time
  void set_restart_time(double time_sec) { restart_sec_ = time_sec; }

 private:
  bool value(const std::string& name, const std::string& varname, const std::vector<std::size_t>& dims,
             std::vector<double>& out);
  wrf_io::DomainAtts domain_atts();
  void start_em_slopes();
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
  double auxinput4_interval_s_ = 0;  // sst_update: the AUXINPUT4 alarm
  int wrfinput_id_ = -1;
  std::vector<double> p_hyd0_;  // P_HYD of the initial state (start_em)
  std::unique_ptr<phys::driver::PhyPrepFields> prep_;  // runs without column physics
  std::map<long, std::vector<double>> p_hyd_snap_;  // P_HYD of phy_prep at the start of the step
  std::map<long, std::vector<double>> t_snap_;  // T at the end of the solve, by step
  bool same_m_ = false, same_u_ = false, same_v_ = false;
  bool model_only_ = false;
  std::map<std::string, std::vector<double>> fed_;  // nest-feedback values of such fields
  std::map<std::string, std::vector<double>> phy_init_;  // phy_init fields of runs without column physics
  long step_ = 0;
  long prep_step_ = -1;
  double restart_sec_ = 0;
  double time_sec_ = 0;
  std::string last_file_, last_restart_;
};

}  // namespace cppwrf::cases
#endif
