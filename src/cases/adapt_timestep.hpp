#pragma once
// WRF's adaptive time step (use_adaptive_time_step): dyn_em/adapt_timestep_em.F adapt_timestep and calc_dt, and
// the adaptive block of dyn_em/start_em.F, for a single domain (the nest branches are not ported).
#include "core/precision.hpp"
#include "core/wrf_time.hpp"

namespace cppwrf::cases {

struct AdaptOptions {  // &domains (Registry defaults) and the &time_control / &bdy_control values adapt_timestep reads
  real target_cfl = real(1.2), target_hcfl = real(0.84);
  int max_step_increase_pct = 5;
  int starting_time_step = -1, starting_time_step_den = 0;
  int max_time_step = -1, max_time_step_den = 0;
  int min_time_step = -1, min_time_step_den = 0;
  bool step_to_output_time = true;
  int history_interval = 60, history_interval_d = 0, history_interval_h = 0, history_interval_m = 0,
      history_interval_s = 0;
  int interval_seconds = 21600;
};

class AdaptTimeStep {
 public:
  explicit AdaptTimeStep(const AdaptOptions& o) : o_(o) {}

  // start_em.F:932-1000 before phy_init: the unset starting / max / min steps from dx, dy, grid%dt = the starting
  // step, the CFL memory zeroed, max_msftx / max_msfty (MAXVAL of the map factors); then (not a restart) the first
  // adapt_timestep. itimestep: grid%itimestep (0, or the restart file's).
  void start_em(long itimestep, double dx, double dy, real max_msftx, real max_msfty, bool restart);
  // solve_em.F:352: at the start of step itimestep + 1; now = current time - this run's start time,
  // dtbc = grid%dtbc (before this step's dt is added)
  void adapt(long itimestep, const wrf_time::Interval& now, real dtbc) { adapt_timestep(itimestep, now, dtbc, false); }
  // solve_em.F:919-934: the maxima of w_damp's vertical and horizontal CFL numbers of one RK stage
  void add_cfl(real vert, real horiz) {
    if (horiz > max_horiz_cfl) max_horiz_cfl = horiz;
    if (vert > max_vert_cfl) max_vert_cfl = vert;
  }
  // solve_em.F:484-486, after adapt_timestep
  void reset_cfl() { max_horiz_cfl = 0; max_vert_cfl = 0; }

  const AdaptOptions& options() const { return o_; }
  real dt() const { return dt_; }                                  // grid%dt
  const wrf_time::Interval& time_step() const { return timestep_; }  // the clock's TimeStep
  // wrfout / wrfrst attributes ADAPT_DT_START / _MAX / _MIN (output_wrf.F:176-178)
  int starting_time_step() const { return o_.starting_time_step; }
  int max_time_step() const { return o_.max_time_step; }
  int min_time_step() const { return o_.min_time_step; }
  real max_msftx() const { return max_msftx_; }
  real max_msfty() const { return max_msfty_; }

  // restart fields (Registry.EM_COMMON:366-374)
  int last_dt_sec = 0, last_dt_sec_num = 0, last_dt_sec_den = 0;
  int last_step_updated = -1;
  bool stepping_to_time = false;
  // domain-type members (frame/module_domain_type.F:225-230), not in the restart file
  real max_vert_cfl = 0, max_horiz_cfl = 0, last_max_vert_cfl = 0, last_max_horiz_cfl = 0;

 private:
  void adapt_timestep(long itimestep, const wrf_time::Interval& now, real dtbc, bool first);
  AdaptOptions o_;
  real dt_ = 0;
  real max_msftx_ = 0, max_msfty_ = 0;
  wrf_time::Interval timestep_;
};

}  // namespace cppwrf::cases
