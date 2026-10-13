#include "dyn/rk3.hpp"

#include <algorithm>
#include <stdexcept>

#ifdef CPPWRF_HAVE_SMALL_STEP
#include "core/halo.hpp"
#include "dyn/lateral_bc.hpp"
#include "dyn/small_step.hpp"
#endif

namespace cppwrf::dyn {

#ifdef CPPWRF_HAVE_SMALL_STEP
class SmallStepHolder {
 public:
  SmallStepHolder(const Grid& g, const AcousticParams& a)
      : ss(g, SmallStepParams{.epssm = a.epssm, .smdiv = a.smdiv, .emdiv = a.emdiv, .top_lid = a.top_lid,
                              .damp_opt = a.damp_opt, .hypsometric_opt = a.hypsometric_opt, .zdamp = a.zdamp,
                              .dampcoef = a.dampcoef}) {}
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

void RK3::set_lateral_bc(LateralBC* bc) {
  bc_ = bc;
  large_.set_lateral_bc(bc);
  scalars_.set_lateral_bc(bc);
#ifdef CPPWRF_HAVE_SMALL_STEP
  small_->ss.set_lateral_bc(bc);
#endif
}

int RK3::auto_sound_steps(double dt, double dx, double dy, bool two_d) {
  double spacing = two_d ? dx : std::min(dx, dy);
  return std::max(2 * (int(300.0 * dt / spacing - 0.01) + 1), 4);
}

void RK3::step(State& s, Tendencies& t, const BaseState& b) {
#ifndef CPPWRF_HAVE_SMALL_STEP
  (void)s; (void)t; (void)b;
  throw std::runtime_error("acoustic small step (src/dyn/small_step.*) not available in this build");
#else
  if (bc_) bc_->advance_time(g_.dt);  // solve_em: dtbc = dtbc + dt
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
    if (rk == 1 && cp_) {  // first_rk_step_part1: init_zero_tendency, then the column physics
      large_.zero_tendf(t);
      cp_->part1(s, b, t);
    }
    large_.tendencies(s, t, b, rk);  // part2 (diffusion, update_phy_ten) + rk_tendency
    small_->ss.run_stage(s, t, b, real(dt_rk), real(dts_rk), ns, rk);
    scalars_.run_stage(s, t, b, rk, real(dt_rk));
    small_->ss.calc_p_rho_phi(s, b);
  }
  if (cp_) cp_->after_rk(s, b);  // phy_prep_part2
  if (phys_ && phys_->params().mp_physics != 0) {  // moist_physics_prep/microphysics/finish, calc_p_rho_phi
    phys_->step(s, b, real(dt));
    small_->ss.calc_p_rho_phi(s, b);
  }
  if (bc_) {  // spec_bdy_final: exact boundary values at t + dt
    bc_->final(s, b);
    fill_halos({&s.u_2, &s.v_2, &s.t_2, &s.ph_2}, g_);
    fill_halo(s.mu_2, g_);
  }
  // solve_em: "reset surface w for consistency" (kinematic condition from the new u, v)
  if (!b.ht.uniform()) set_w_surface(s.w_2, s.u_2, s.v_2, b, g_, false);
#endif
}

namespace {
struct KSetWSurface {  // WRF module_bc_em set_w_surface
  View3<real> w; View3<const real> u, v; BView2<const real> ht, msftx, msfty; View1<const real> znw;
  real cf1, cf2, cf3, rdx, rdy; int nz; bool fill;
  CPPWRF_HD void operator()(int i, int j) const {
    const real ws =
        msfty(i, j) * real(0.5) * rdy *
            ((ht(i, j + 1) - ht(i, j)) * (cf1 * v(i, j + 1, 0) + cf2 * v(i, j + 1, 1) + cf3 * v(i, j + 1, 2)) +
             (ht(i, j) - ht(i, j - 1)) * (cf1 * v(i, j, 0) + cf2 * v(i, j, 1) + cf3 * v(i, j, 2))) +
        msftx(i, j) * real(0.5) * rdx *
            ((ht(i + 1, j) - ht(i, j)) * (cf1 * u(i + 1, j, 0) + cf2 * u(i + 1, j, 1) + cf3 * u(i + 1, j, 2)) +
             (ht(i, j) - ht(i - 1, j)) * (cf1 * u(i, j, 0) + cf2 * u(i, j, 1) + cf3 * u(i, j, 2)));
    w(i, j, 0) = ws;
    if (fill)
      for (int k = 1; k <= nz; ++k) w(i, j, k) = ws * znw(k) * znw(k);
  }
};
}  // namespace

// Periodic lateral boundaries only (the halo supplies i-1, i+1); specified boundaries clamp
// the stencil at the domain edge in WRF (TODO with the lateral boundary integration).
void set_w_surface(Field3& w, Field3& u, Field3& v, const BaseState& b, const Grid& g, bool fill) {
  fill_halos({&u, &v}, g);
  column_for(g.columns(), KSetWSurface{w.view(), u.cview(), v.cview(), b.ht.view(), b.msftx.view(), b.msfty.view(),
                                       b.znw.view(), real(b.cf1), real(b.cf2), real(b.cf3), real(b.rdx),
                                       real(b.rdy), g.nz, fill});
  fill_halo(w, g);
}

}  // namespace cppwrf::dyn
