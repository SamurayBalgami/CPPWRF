#include "dyn/rk3.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <initializer_list>
#include <stdexcept>
#include <vector>

#ifdef CPPWRF_HAVE_SMALL_STEP
#include "core/halo.hpp"
#include "core/vm_dump.hpp"
#include "dyn/diffusion_driver.hpp"
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
    : g_(g), large_(g, lp), scalars_(g, sp), n_sound_(tss > 0 ? tss : auto_sound_steps(g.dt, g.dx, g.dy)) {
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

// In WRF's REAL, as adaptive_sound_steps: 300*dt/spacing - 0.01 can land on either side of an
// integer (dx 1000 m, dt 6.7 s: 1.99999... in double -> 4 steps, 2 in float -> 6 steps)
int RK3::auto_sound_steps(double dt, double dx, double dy) {
  const real spacing = std::min(real(dx), real(dy));  // solve_em: also in 2-D (dy of the namelist)
  return std::max(2 * (int(real(300.) * real(dt) / spacing - real(0.01)) + 1), 4);
}

int RK3::adaptive_sound_steps(real dt, real dx, real dy, real max_msft) {
  const real spacing = std::min(dx, dy);
  return std::max(2 * (int(real(300.) * dt / (spacing / max_msft) - real(0.01)) + 1), 4);
}

void RK3::set_dt(double dt, int n_sound) {
  if (n_sound < 2 || n_sound % 2) throw std::runtime_error("time_step_sound must be even and >= 2");
  g_.dt = dt;
  n_sound_ = n_sound;
  large_.set_dt(dt);
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
  fill_halo(s.mu_1, g_);  // WRF: MU_1 = MU_2 on the domain, then set_physical_bc2d(mu_1) fills its halo
  for (int n = 0; n < s.n_scalars(); ++n) s.q_1[n].copy_from(s.q_2[n]);  // rk_update_scalar, rk_step 1
  const double dt = g_.dt, dts = dt / n_sound_;
  for (int rk = 1; rk <= 3; ++rk) {
    double dt_rk = rk == 1 ? dt / 3.0 : (rk == 2 ? 0.5 * dt : dt);
    double dts_rk = rk == 1 ? dt_rk : dts;
    int ns = rk == 1 ? 1 : (rk == 2 ? n_sound_ / 2 : n_sound_);
    large_.prep(s, b);
    if (rk == 1 && on_first_rk_prep) on_first_rk_prep();
    if (rk == 1 && cp_) {  // first_rk_step_part1: init_zero_tendency, then the column physics
      large_.zero_tendf(t);
      { const char* e = std::getenv("CPPWRF_VM_DUMP_STEP"); vmdump::on() = vm_step_ == (e ? std::atol(e) - 1 : 0); }
      cp_->part1(s, b, t);
      vmdump::on() = false;
    }
    if (rk == 1) { const char* e = std::getenv("CPPWRF_VM_DUMP_STEP"); vmdump::on() = vm_step_ == (e ? std::atol(e) - 1 : 0); }
    if (rk == 1 && dd_) dd_->prepare(s, b);  // first_rk_step_part2: metrics, deformation, K
    large_.tendencies(s, t, b, rk);  // part2 (diffusion, update_phy_ten) + rk_tendency
    if (on_cfl) {  // w_damp's CFL maxima of this stage (solve_em.F:919-934)
      real vert, horiz;
      large_.cfl_max(s, b, g_.dt, vert, horiz);
      on_cfl(vert, horiz);
    }
    vmdump::on() = false;
    { const char* e = std::getenv("CPPWRF_VM_DUMP_STEP"); vmdump::on() = vm_step_ == (e ? std::atol(e) - 1 : 0); }
    if (rk == 1 && bc_ && vmdump::on()) {  // debug: boundary slabs (values, tendencies) of T and U, all sides
      if (std::FILE* f = vmdump::open("vm_bdy", 0)) {
        auto* bc = const_cast<LateralBC*>(bc_);
        for (LateralBC::Var v : {LateralBC::T, LateralBC::U})
          for (int sd = 0; sd < 4; ++sd)
            for (bool tend : {false, true}) {
              const auto side = LateralBC::Side(sd);
              const int hdr[3] = {bc->slab_npos(v, side), bc->slab_nk(v), bc->width()};
              std::fwrite(hdr, sizeof(int), 3, f);
              const std::size_t n = std::size_t(hdr[0]) * hdr[1] * hdr[2];
              std::vector<real> h(n);
              copy_to_host(h.data(), bc->slab(v, side, tend), n * sizeof(real));
              for (real x : h) { const double d = double(x); std::fwrite(&d, sizeof(double), 1, f); }
            }
        const double dt = bc_->dtbc;
        std::fwrite(&dt, sizeof(double), 1, f);
        std::fclose(f);
      }
    }
    vmdump::dump("vm_tend_rk", rk, {&t.ru_tend, &t.rv_tend, &t.rw_tend, &t.ph_tend, &t.t_tend, &s.ru, &s.rv, &s.php, &s.alt, &s.ww}, {&s.muu, &s.mut});
    small_->ss.run_stage(s, t, b, real(dt_rk), real(dts_rk), ns, rk);
    scalars_.run_stage(s, t, b, rk, real(dt_rk));
    vmdump::dump("vm_state_rk", rk, {&s.u_2, &s.v_2, &s.w_2, &s.t_2, &s.ph_2, s.P_QV >= 0 ? &s.q_2[s.P_QV] : nullptr, &s.ww, &s.ru_m, &s.rv_m, &s.ww_m}, {&s.mu_2, &s.muts});
    small_->ss.calc_p_rho_phi(s, b);
    vmdump::dump("vm_p_rk", rk, {&s.p, &s.al});
    vmdump::on() = false;
  }
  ++vm_step_;
  if (cp_) cp_->after_rk(s, b);  // phy_prep_part2
  if (phys_ && phys_->params().mp_physics != 0) {  // moist_physics_prep/microphysics/finish, calc_p_rho_phi
    { const char* e = std::getenv("CPPWRF_VM_DUMP_STEP"); vmdump::on() = vm_step_ - 1 == (e ? std::atol(e) - 1 : 0); }
    phys_->step(s, b, real(dt));
    vmdump::on() = false;
    // solve_em scalar_tile_loop_1ba: calc_p_rho_phi on the tile without the spec zone
    if (bc_) small_->ss.calc_p_rho_phi_interior(s, b, bc_->spec_zone());
    else small_->ss.calc_p_rho_phi(s, b);
  }
  if (bc_) {  // spec_bdy_final: exact boundary values at t + dt
    bc_->final(s, b);
    // WRF (serial) does not refresh any halo after spec_bdy_final: until the next refresh the halos
    // keep the values from before it (read e.g. by calc_cq for q and by calc_ww_cp for mu_2)
    small_->ss.calc_p_rho_phi_interior(s, b, bc_->spec_zone());  // solve_em after spec_bdy_final
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

// The ht halo supplies i-1, i+1. WRF clamps the stencil at non-periodic edges (im1 = ids,
// ip1 = ide-1), i.e. ht(i) - ht(im1) = 0 there; the flat (specified) and mirror (symmetric)
// mass-point halos of ht give exactly the same value.
void set_w_surface(Field3& w, Field3& u, Field3& v, const BaseState& b, const Grid& g, bool fill) {
  fill_halos({&u, &v}, g);
  column_for(g.columns(), KSetWSurface{w.view(), u.cview(), v.cview(), b.ht.view(), b.msftx.view(), b.msfty.view(),
                                       b.znw.view(), real(b.cf1), real(b.cf2), real(b.cf3), real(b.rdx),
                                       real(b.rdy), g.nz, fill});
  fill_halo(w, g);
}

}  // namespace cppwrf::dyn
