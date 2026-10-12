#include "dyn/rk3.hpp"

#include <algorithm>
#include <stdexcept>

#ifdef CPPWRF_HAVE_SMALL_STEP
#include "dyn/small_step.hpp"
#endif

namespace cppwrf::dyn {

#ifdef CPPWRF_HAVE_SMALL_STEP
class SmallStepHolder {
 public:
  SmallStepHolder(const Grid& g, const AcousticParams& a)
      : ss(g, SmallStepParams{a.epssm, a.smdiv, a.emdiv}) {}
  SmallStep ss;
};
bool RK3::small_step_available() { return true; }
#else
class SmallStepHolder {};
bool RK3::small_step_available() { return false; }
#endif

RK3::RK3(const Grid& g, const LargeStepParams& lp, const AcousticParams& ap, int tss, const ScalarParams& sp)
    : g_(g), large_(g, lp), scalars_(g, sp), n_sound_(tss > 0 ? tss : auto_sound_steps(g.dt, g.dx, g.dy, g.two_d())) {
  if (n_sound_ < 2 || n_sound_ % 2) throw std::runtime_error("time_step_sound must be even and >= 2");
#ifdef CPPWRF_HAVE_SMALL_STEP
  small_ = std::make_unique<SmallStepHolder>(g, ap);
#else
  (void)ap;
#endif
}

RK3::~RK3() = default;

int RK3::auto_sound_steps(double dt, double dx, double dy, bool two_d) {
  double spacing = two_d ? dx : std::min(dx, dy);
  return std::max(2 * (int(300.0 * dt / spacing - 0.01) + 1), 4);
}

void RK3::step(State& s, Tendencies& t, const BaseState& b) {
#ifndef CPPWRF_HAVE_SMALL_STEP
  (void)s; (void)t; (void)b;
  throw std::runtime_error("acoustic small step (src/dyn/small_step.*) not available in this build");
#else
  // time-t state (WRF small_step_prep, rk_step == 1)
  s.u_1.copy_from(s.u_2); s.v_1.copy_from(s.v_2); s.w_1.copy_from(s.w_2);
  s.t_1.copy_from(s.t_2); s.ph_1.copy_from(s.ph_2); s.mu_1.copy_from(s.mu_2);
  for (int n = 0; n < s.n_scalars(); ++n) s.q_1[n].copy_from(s.q_2[n]);  // rk_update_scalar, rk_step 1
  const double dt = g_.dt, dts = dt / n_sound_;
  for (int rk = 1; rk <= 3; ++rk) {
    double dt_rk = rk == 1 ? dt / 3.0 : (rk == 2 ? 0.5 * dt : dt);
    double dts_rk = rk == 1 ? dt_rk : dts;
    int ns = rk == 1 ? 1 : (rk == 2 ? n_sound_ / 2 : n_sound_);
    large_.prep(s, b);
    large_.tendencies(s, t, b, rk);
    small_->ss.run_stage(s, t, b, real(dt_rk), real(dts_rk), ns, rk);
    scalars_.run_stage(s, t, b, rk, real(dt_rk));
    small_->ss.calc_p_rho_phi(s, b);
  }
  if (phys_ && phys_->params().mp_physics != 0) {  // moist_physics_prep/microphysics/finish, calc_p_rho_phi
    phys_->step(s, b, real(dt));
    small_->ss.calc_p_rho_phi(s, b);
  }
#endif
}

}  // namespace cppwrf::dyn
