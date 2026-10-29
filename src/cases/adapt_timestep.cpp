#include "cases/adapt_timestep.hpp"

#include <algorithm>
#include <cmath>

namespace cppwrf::cases {

using wrf_time::Interval;

namespace {

constexpr int precision = 100;

// Fortran NINT / INT of a REAL
long nint(double x) { return std::lround(x); }
long fint(double x) { return long(x); }

// adapt_timestep_em.F calc_dt
Interval calc_dt(real max_cfl, real max_increase_factor, const Interval& last_dtInterval, real target_cfl) {
  int num, den;
  if (max_cfl < real(0.001)) {
    // If the max_cfl is small, then we increase dtInterval the maximum amount allowable.
    num = int(fint(max_increase_factor * real(precision) + real(0.5)));
    den = precision;
    return last_dtInterval * num / den;
  }
  real factor;
  if (max_cfl > target_cfl) {
    // reducing: below the target cfl by half the difference between max and target
    factor = (target_cfl - real(0.5) * (max_cfl - target_cfl)) / max_cfl;
    factor = std::max(real(0.1), factor);
  } else {
    factor = target_cfl / max_cfl;
  }
  num = int(fint(factor * real(precision) + real(0.5)));
  den = precision;
  return last_dtInterval * num / den;
}

}  // namespace

void AdaptTimeStep::start_em(long itimestep, double dx, double dy, real max_msftx, real max_msfty, bool restart) {
  if (last_step_updated != itimestep) {
    const real dmin = std::min(real(dx), real(dy));
    if (o_.starting_time_step == -1) o_.starting_time_step = int(nint(real(4) * dmin / real(1000)));
    if (o_.max_time_step == -1) o_.max_time_step = int(nint(real(8) * dmin / real(1000)));
    if (o_.min_time_step == -1) o_.min_time_step = int(nint(real(3) * dmin / real(1000)));
    // Set a starting timestep.
    dt_ = real(o_.starting_time_step);
    timestep_ = wrf_time::make(o_.starting_time_step);  // setup_timekeeping: the namelist time_step (replaced below)
    last_max_vert_cfl = 0;
    last_max_horiz_cfl = 0;
    max_msftx_ = max_msftx;
    max_msfty_ = max_msfty;
  }
  // This first call just initializes variables; a restart takes them from the restart file.
  if (!restart) adapt_timestep(itimestep, Interval{}, real(0), true);
}

void AdaptTimeStep::adapt_timestep(long itimestep, const Interval& now, real dtbc, bool first) {
  bool use_last2 = false;
  // Assign last_dtInterval type values from restart file
  Interval grid_last_dtInterval = wrf_time::make(last_dt_sec, last_dt_sec_num, last_dt_sec_den);
  // If this step has already been adapted, no need to do it again.
  if (last_step_updated == itimestep) return;
  last_step_updated = int(itimestep);

  Interval last_dtInterval = grid_last_dtInterval;
  // time since the beginning of this run (current time - start time)
  real curr_secs = wrf_time::real_time<real>(now);
  const real max_increase_factor = real(1.) + real(o_.max_step_increase_pct) / real(100.);

  Interval dtInterval;
  if (first) {  // domain_get_advanceCount == 0 .AND. .NOT. restart
    dtInterval = o_.starting_time_step_den == 0 ? wrf_time::make(0, o_.starting_time_step, 1)
                                                : wrf_time::make(0, o_.starting_time_step, o_.starting_time_step_den);
    curr_secs = 0;
    last_dtInterval = wrf_time::make(0, 0, 1);
  } else {
    real mv, mh;
    if (stepping_to_time) {
      mv = last_max_vert_cfl;
      mh = last_max_horiz_cfl;
    } else {
      mv = max_vert_cfl;
      mh = max_horiz_cfl;
    }
    const Interval dtInterval_vert = calc_dt(mv, max_increase_factor, last_dtInterval, o_.target_cfl);
    const Interval dtInterval_horiz = calc_dt(mh, max_increase_factor, last_dtInterval, o_.target_hcfl);
    dtInterval = dtInterval_vert < dtInterval_horiz ? dtInterval_vert : dtInterval_horiz;
  }

  // Limit the increase of dtInterval to the specified input limit
  int num = int(nint(max_increase_factor * real(precision)));
  int den = precision;
  Interval tmpTimeInterval = last_dtInterval * num / den;
  if (now != Interval{} && dtInterval > tmpTimeInterval) dtInterval = tmpTimeInterval;

  // Round dtInterval to the nearest 1/100 (dt is a double precision local)
  double dt = double(wrf_time::real_time<real>(dtInterval));
  num = int(nint(dt * precision));
  den = precision;
  dtInterval = wrf_time::make(0, num, den);

  // Limit the maximum and the minimum dtInterval based on user input
  tmpTimeInterval = o_.max_time_step_den == 0 ? wrf_time::make(0, o_.max_time_step, 1)
                                              : wrf_time::make(0, o_.max_time_step, o_.max_time_step_den);
  if (dtInterval > tmpTimeInterval) dtInterval = tmpTimeInterval;
  tmpTimeInterval = o_.min_time_step_den == 0 ? wrf_time::make(0, o_.min_time_step, 1)
                                              : wrf_time::make(0, o_.min_time_step, o_.min_time_step_den);
  if (dtInterval < tmpTimeInterval) dtInterval = tmpTimeInterval;

  // Assure that we fall on a BC time (only the step of the domain with lateral boundary times)
  stepping_to_time = false;
  bool stepping_to_bc;
  const real time_to_bc = real(o_.interval_seconds) - dtbc;
  num = int(fint(time_to_bc * real(precision) + real(0.5)));
  den = precision;
  tmpTimeInterval = wrf_time::make(0, num, den);
  if (tmpTimeInterval < dtInterval * 2 && tmpTimeInterval > dtInterval) {
    dtInterval = tmpTimeInterval / 2;
    use_last2 = true;
    stepping_to_bc = true;
    stepping_to_time = true;
  } else if (tmpTimeInterval <= dtInterval) {
    const long bc_time = nint((curr_secs + time_to_bc) / real(o_.interval_seconds)) * o_.interval_seconds;
    tmpTimeInterval = wrf_time::make(bc_time);
    dtInterval = tmpTimeInterval - now;
    use_last2 = true;
    stepping_to_bc = true;
    stepping_to_time = true;
  } else {
    stepping_to_bc = false;
  }

  // step_to_output_time: fall on an output time, looking out two time steps (no very short step)
  if (o_.step_to_output_time && !stepping_to_bc) {
    if (o_.history_interval_m == 0) o_.history_interval_m = o_.history_interval;
    const int history_interval_sec = o_.history_interval_s + o_.history_interval_m * 60 +
                                     o_.history_interval_h * 3600 + o_.history_interval_d * 86400;
    const real time_to_output = real(history_interval_sec) - std::fmod(curr_secs, real(history_interval_sec));
    num = int(fint(time_to_output * real(precision) + real(0.5)));
    den = precision;
    tmpTimeInterval = wrf_time::make(0, num, den);
    if (tmpTimeInterval < dtInterval * 2 && tmpTimeInterval > dtInterval) {
      dtInterval = tmpTimeInterval / 2;
      use_last2 = true;
      stepping_to_time = true;
    } else if (tmpTimeInterval <= dtInterval) {
      // round to the nearest history time to avoid rounding errors
      const long output_time =
          nint((curr_secs + time_to_output) / real(history_interval_sec)) * history_interval_sec;
      tmpTimeInterval = wrf_time::make(output_time);
      dtInterval = tmpTimeInterval - now;
      use_last2 = true;
      stepping_to_time = true;
    }
  }

  // (adapt_step_using_child: adaptation_domain with max_dom = 2 only)
  if (use_last2) {
    grid_last_dtInterval = last_dtInterval;
  } else {
    grid_last_dtInterval = dtInterval;
    last_max_vert_cfl = max_vert_cfl;
    last_max_horiz_cfl = max_horiz_cfl;
  }
  dt_ = wrf_time::real_time<real>(dtInterval);
  last_max_vert_cfl = max_vert_cfl;
  // Update the clock based on the new time step
  timestep_ = dtInterval;
  // Update last timestep info for restart file
  last_dt_sec = int(grid_last_dtInterval.S);
  last_dt_sec_num = int(grid_last_dtInterval.Sn);
  last_dt_sec_den = int(grid_last_dtInterval.Sd);
}

}  // namespace cppwrf::cases
