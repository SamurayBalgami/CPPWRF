#pragma once

#include <chrono>
#include <memory>
#include <string>
#include <vector>

#include "cases/model.hpp"
#include "dyn/diagnostics.hpp"
#include "io/config.hpp"
#include "io/output.hpp"

#ifdef CPPWRF_HAVE_NETCDF
#include "io_netcdf/restart.hpp"
#include "cases/wrf_history.hpp"
#endif

namespace cppwrf::cases {

class DomainRunner {
 public:
  DomainRunner(const io::Config& cfg, cases::Model& m, const std::string& out_dir,
               std::string domain_tag = "d01");
  ~DomainRunner();

  DomainRunner(const DomainRunner&) = delete;
  DomainRunner& operator=(const DomainRunner&) = delete;

  // Initializes domain state, loads restarts if needed, writes run.json, diag header, first diag/snap
  bool start();
  // output field lists and the wrfout writer (idempotent; start() calls it)
  void prepare_output();

  // Called before step n (1-based): boundary interval updates, start timing
  void before_step(long n);

  // Called after step n (1-based): step timing, diag, snapshot, wrfout, restart
  bool after_step(long n);
  // Called right after the solve of step n, before nests run and feed back (WRF's end of solve_em)
  void after_solve(long n);
  // Nested domains: WRF writes history in med_before_solve_io, i.e. a frame for step n at the start
  // of step n+1 (after the parent's next med_nest_force, whose couple/uncouple round trip changes
  // the state by an ulp); the last frame after the last step (med_last_solve_io).
  void set_defer_output(bool d) { defer_output_ = d; }
  // two-way nests: the WRF-field view (WrfHistory) carries the physics feedback, also when no wrfout is written
  void set_need_history(bool h) { need_history_ = h; }

  // Prints completion statistics
  void finish();

  long start_step = 0;
  double step_time = 0.0;

  const std::string& out_dir() const { return out_dir_; }
#ifdef CPPWRF_HAVE_NETCDF
  WrfHistory* history() { return wrfout_.get(); }
#endif
  const std::string& domain_tag() const { return domain_tag_; }

 private:
  bool defer_output_ = false;
  bool need_history_ = false;
  bool output_prepared_ = false;
  long pending_ = -1;  // deferred diag/snap step
  bool do_diag(long n);
  void do_snap(long n);

#ifdef CPPWRF_HAVE_NETCDF
  void do_restart(long n, bool ring);
  void save_restart_file(const std::string& path, long step, double time_sec);
  void load_restart_file(const std::string& path);
#endif

  const io::Config& cfg_;
  cases::Model& m_;
  std::string out_dir_;
  std::string domain_tag_;

  long diag_iv_ = 10;
  long out_iv_ = 1;
  double rst_iv_s_ = 0.0;  // restart interval (s), WRF's RESTART_ALARM
  double rst_prev_ = 0.0;  // the alarm's PrevRingTime (s)
  bool rst_due(long n) const;
  long rst_end_ = 0;  // restart.end: last step with a restart file (0: no end)
  bool stop_on_nonfinite_ = true;
  std::string start_date_;
  std::string cname_;
  std::string out_fmt_;
  bool write_native_ = true;
  bool write_wrfout_ = false;

  std::chrono::steady_clock::time_point t0_;
  std::chrono::steady_clock::time_point step_start_;

  std::unique_ptr<io::DiagWriter> diag_;
  std::unique_ptr<io::SnapshotWriter> snap_;
  std::unique_ptr<dyn::Diagnostics> dg_;
  std::vector<std::vector<real>> store_;
  std::vector<const Field2*> extra_2d_;
  std::vector<const Field3*> extra_3d_;

#ifdef CPPWRF_HAVE_NETCDF
  std::unique_ptr<WrfHistory> wrfout_;
#endif
};

}  // namespace cppwrf::cases
