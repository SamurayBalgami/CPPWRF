#include "dyn/scalars.hpp"

#include <stdexcept>
#include <vector>

#include "core/backend.hpp"
#include "core/halo.hpp"
#include "dyn/diffusion_driver.hpp"
#include "dyn/kernels_common.hpp"
#include "dyn/lateral_bc.hpp"

// WRF dyn_em/module_advect_em.F (advect_scalar, advect_scalar_pd) and module_em.F
// (rk_scalar_tend, rk_update_scalar, rk_update_scalar_pd). C++ k = Fortran k - 1.
// With map factors (msftx, msfty as in advect_scalar*) and periodic, symmetric or specified
// lateral boundaries (degraded flux orders and the spec zone as WRF for specified).

namespace cppwrf::dyn {
namespace {

using V3 = View3<real>;
using C3 = View3<const real>;
using C2 = View2<const real>;
using B2 = BView2<const real>;
constexpr real H = real(0.5);

// High-order (advect_scalar) and first-order upwind (advect_scalar_pd) face fluxes.
struct Fluxes {
  C3 q, qo, ru, rv, ww; C2 mut; VCoef c;
  real rdx, rdy, dt, hup, vup; int nz; bool do_y;
  B2 msftx, msfty; kern::Spec sp;
  CPPWRF_HD real fx(int i, int j, int k) const {  // west face of cell i
    const real vel = ru(i, j, k);
    const int o = kern::ord_face(i, sp.nx, sp.x);
    if (o == 2) return H * vel * (q(i, j, k) + q(i - 1, j, k));
    if (o == 3) return vel * kern::flux3(q(i - 2, j, k), q(i - 1, j, k), q(i, j, k), q(i + 1, j, k), vel, hup);
    return vel * kern::flux5(q(i - 3, j, k), q(i - 2, j, k), q(i - 1, j, k), q(i, j, k), q(i + 1, j, k), q(i + 2, j, k), vel, hup);
  }
  CPPWRF_HD real fy(int i, int j, int k) const {
    const real vel = rv(i, j, k);
    const int o = kern::ord_face(j, sp.ny, sp.y);
    if (o == 2) return H * vel * (q(i, j, k) + q(i, j - 1, k));
    if (o == 3) return vel * kern::flux3(q(i, j - 2, k), q(i, j - 1, k), q(i, j, k), q(i, j + 1, k), vel, hup);
    return vel * kern::flux5(q(i, j - 3, k), q(i, j - 2, k), q(i, j - 1, k), q(i, j, k), q(i, j + 1, k), q(i, j + 2, k), vel, hup);
  }
  CPPWRF_HD real fz(int i, int j, int k) const {  // w level k
    if (k == 0 || k == nz) return 0;
    const real vel = ww(i, j, k);
    if (k == 1 || k == nz - 1) return vel * (c.fnm(k) * q(i, j, k) + c.fnp(k) * q(i, j, k - 1));
    return vel * kern::flux3(q(i, j, k - 2), q(i, j, k - 1), q(i, j, k), q(i, j, k + 1), -vel, vup);
  }
  CPPWRF_HD real mh(int i, int j, int k) const { return c.c1h(k) * mut(i, j) + c.c2h(k); }
  CPPWRF_HD real fxl(int i, int j, int k) const {
    const real dx = real(2) / (msfty(i, j) + msfty(i - 1, j)) / rdx, mu = H * (mh(i, j, k) + mh(i - 1, j, k)),
               vel = ru(i, j, k);
    const real cr = vel * dt / dx / mu;
    return mu * (dx / dt) * kern::flux_upwind(qo(i - 1, j, k), qo(i, j, k), cr);
  }
  CPPWRF_HD real fyl(int i, int j, int k) const {
    const real dy = real(2) / (msftx(i, j) + msftx(i, j - 1)) / rdy, mu = H * (mh(i, j, k) + mh(i, j - 1, k)),
               vel = rv(i, j, k);
    const real cr = vel * dt / dy / mu;
    return mu * (dy / dt) * kern::flux_upwind(qo(i, j - 1, k), qo(i, j, k), cr);
  }
  CPPWRF_HD real fzl(int i, int j, int k) const {
    if (k == 0 || k == nz) return 0;
    const real dz = real(2) / (c.rdnw(k) + c.rdnw(k - 1)), mu = mh(i, j, k), vel = ww(i, j, k);
    const real cr = vel * dt / dz / mu;
    return mu * (dz / dt) * kern::flux_upwind(qo(i, j, k - 1), qo(i, j, k), cr);
  }
  // spec zone of a specified boundary: WRF advect_scalar(_pd) does not compute it (tendency 0)
  CPPWRF_HD bool spec_cell(int i, int j) const {
    return (sp.x && (i == 0 || i == sp.nx - 1)) || (sp.y && (j == 0 || j == sp.ny - 1));
  }
  CPPWRF_HD real adv(int i, int j, int k) const {  // advect_scalar: horizontal terms times msftx
    if (spec_cell(i, j)) return 0;
    real r = -(msftx(i, j) * rdx) * (fx(i + 1, j, k) - fx(i, j, k));
    if (do_y) r -= (msftx(i, j) * rdy) * (fy(i, j + 1, k) - fy(i, j, k));
    return r - c.rdnw(k) * (fz(i, j, k + 1) - fz(i, j, k));
  }
};

CPPWRF_HD inline real m_full(const VCoef& c, int k, real mub, real mu) { return c.c1h(k) * (mub + mu) + c.c2h(k); }

// rk_update_scalar (non-PD stages / options): q_new = (m_old q_1 + dt (adv + tendf)) / m_new
struct KUpdate {
  Fluxes f; C3 q1, tendf; C2 mu1, mu2; B2 mub; V3 qn; real dt_rk; bool add_tendf;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    real tend = f.adv(i, j, k) * f.msfty(i, j);  // rk_update_scalar: advect_tend * msfty + sc_tend
    if (add_tendf) tend += tendf(i, j, k);
    qn(i, j, k) = (m_full(f.c, k, mub(i, j), mu1(i, j)) * q1(i, j, k) + dt_rk * tend) /
                  m_full(f.c, k, mub(i, j), mu2(i, j));
  }
};

// rk_update_scalar_pd (last stage, before advection): q_1 += dt tendf / m_old
struct KPreUpdatePD {
  V3 q1; C3 tendf; C2 mu1; B2 mub; VCoef c; real dt;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const real m = m_full(c, k, mub(i, j), mu1(i, j));
    q1(i, j, k) = (m * q1(i, j, k) + dt * tendf(i, j, k)) / m;
  }
};

// advect_scalar_pd, step 1: every face flux once. dh = high-order minus first-order upwind
// flux, lo = upwind flux, on x (U points), y (V points) and eta (W levels) faces. The limiter and
// the update read them (the same expressions as evaluating Fluxes per cell: bitwise identical).
struct KFacePDX {
  Fluxes f; V3 dh, lo;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const real l = f.fxl(i, j, k);
    dh(i, j, k) = f.fx(i, j, k) - l;
    lo(i, j, k) = l;
  }
};
struct KFacePDY {
  Fluxes f; V3 dh, lo;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const real l = f.fyl(i, j, k);
    dh(i, j, k) = f.fy(i, j, k) - l;
    lo(i, j, k) = l;
  }
};
struct KFacePDZ {
  Fluxes f; V3 dh, lo;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const real l = f.fzl(i, j, k);
    dh(i, j, k) = f.fz(i, j, k) - l;
    lo(i, j, k) = l;
  }
};

struct FacesPD {
  C3 dhx, lox, dhy, loy, dhz, loz;
};

// advect_scalar_pd limiter: scale factor of the outgoing high-minus-low fluxes of each cell
struct KScalePD {
  Fluxes f; FacesPD a; C2 mu_old; B2 mub; V3 scale;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    if (f.spec_cell(i, j)) { scale(i, j, k) = 1; return; }  // WRF limits interior cells only
    const VCoef& c = f.c;
    const real mxy = f.msftx(i, j) * f.msfty(i, j);
    real hlow = f.rdx * (a.lox(i + 1, j, k) - a.lox(i, j, k));
    if (f.do_y) hlow += f.rdy * (a.loy(i, j + 1, k) - a.loy(i, j, k));
    const real ph_low = ((c.c1h(k) * mub(i, j) + c.c2h(k)) + c.c1h(k) * mu_old(i, j)) * f.qo(i, j, k) -
                        f.dt * (mxy * hlow + f.msfty(i, j) * c.rdnw(k) * (a.loz(i, j, k + 1) - a.loz(i, j, k)));
    const real fxp = a.dhx(i + 1, j, k), fxm = a.dhx(i, j, k);
    real hout = f.rdx * ((fxp > 0 ? fxp : real(0)) - (fxm < 0 ? fxm : real(0)));
    if (f.do_y) {
      const real fyp = a.dhy(i, j + 1, k), fym = a.dhy(i, j, k);
      hout += f.rdy * ((fyp > 0 ? fyp : real(0)) - (fym < 0 ? fym : real(0)));
    }
    const real fzp = a.dhz(i, j, k + 1), fzm = a.dhz(i, j, k);
    const real out = mxy * hout + f.msfty(i, j) * c.rdnw(k) * ((fzp < 0 ? fzp : real(0)) - (fzm > 0 ? fzm : real(0)));
    const real flux_out = f.dt * out;
    real s = 1;
    if (flux_out > ph_low) {
      s = ph_low / (flux_out + real(1e-20));
      if (s < 0) s = 0;
    }
    scale(i, j, k) = s;
  }
};

// PD-limited flux divergence + rk_update_scalar (tendf already folded into q_1)
struct KUpdatePD {
  Fluxes f; FacesPD a; C3 scale, q1; C2 mu1, mu2; B2 mub; V3 qn; real dt_rk;
  // limited high-minus-low flux on the west face of (i,j,k): scaled by the upwind (outflow) cell
  CPPWRF_HD real dx(int i, int j, int k) const {
    const real d = a.dhx(i, j, k);
    return d * (d > 0 ? scale(i - 1, j, k) : (d < 0 ? scale(i, j, k) : real(1)));
  }
  CPPWRF_HD real dy(int i, int j, int k) const {
    const real d = a.dhy(i, j, k);
    return d * (d > 0 ? scale(i, j - 1, k) : (d < 0 ? scale(i, j, k) : real(1)));
  }
  CPPWRF_HD real dz(int i, int j, int k) const {  // eta decreases upward: fz > 0 leaves cell k
    if (k == 0 || k == f.nz) return 0;
    const real d = a.dhz(i, j, k);
    return d * (d > 0 ? scale(i, j, k) : (d < 0 ? scale(i, j, k - 1) : real(1)));
  }
  CPPWRF_HD void operator()(int i, int j, int k) const {  // WRF advect_scalar_pd final tendency
    if (f.spec_cell(i, j)) {
      qn(i, j, k) = m_full(f.c, k, mub(i, j), mu1(i, j)) * q1(i, j, k) / m_full(f.c, k, mub(i, j), mu2(i, j));
      return;
    }
    real tend = -f.c.rdnw(k) * (dz(i, j, k + 1) - dz(i, j, k) + a.loz(i, j, k + 1) - a.loz(i, j, k));
    tend -= f.msftx(i, j) * (f.rdx * (dx(i + 1, j, k) - dx(i, j, k) + a.lox(i + 1, j, k) - a.lox(i, j, k)));
    if (f.do_y) tend -= f.msftx(i, j) * (f.rdy * (dy(i, j + 1, k) - dy(i, j, k) + a.loy(i, j + 1, k) - a.loy(i, j, k)));
    qn(i, j, k) = (m_full(f.c, k, mub(i, j), mu1(i, j)) * q1(i, j, k) + dt_rk * (tend * f.msfty(i, j))) /
                  m_full(f.c, k, mub(i, j), mu2(i, j));
  }
};

}  // namespace

ScalarTransport::ScalarTransport(const Grid& g, const ScalarParams& p)
    : g_(g), p_(p), q_new_(g, Stagger::M, "q_new"), scale_(g, Stagger::M, "pd_scale"),
      dhx_(g, Stagger::U, "pd_dhx"), lox_(g, Stagger::U, "pd_lox"), dhy_(g, Stagger::V, "pd_dhy"),
      loy_(g, Stagger::V, "pd_loy"), dhz_(g, Stagger::W, "pd_dhz"), loz_(g, Stagger::W, "pd_loz") {}

void ScalarTransport::run_stage(State& s, Tendencies& t, const BaseState& b, int rk, real dt_rk) {
  const Grid& g = g_;
  const int ns = s.n_scalars();
  if (ns == 0) return;
  if (int(t.q_tendf.size()) != ns) throw std::runtime_error("Tendencies: q_tendf not registered for all scalars");
  const VCoef c = VCoef::from(b);
  const Range3D mr = g.mass_cells();
  const bool do_y = !g.two_d();
  const real hup = p_.h_adv_order == 5 ? real(1) : real(0), vup = p_.v_adv_order == 3 ? real(1) : real(0);
  const bool diff = !dd_ && (p_.khdif > 0 || p_.kvdif > 0);  // legacy constant-K kernels
  if (b.mapped && diff) throw std::runtime_error("scalar diffusion with map factors is not implemented yet");
  for (int n = 0; n < ns; ++n) {
    Field3 &q1 = s.q_1[n], &q2 = s.q_2[n], &tf = t.q_tendf[n];
    const bool moist = n < s.n_moist;
    const bool pd = rk == 3 && (moist ? p_.moist_adv_opt : p_.scalar_adv_opt) == 1;
    if (rk == 1) {  // rk_scalar_tend, rk_step == 1: diffusion into the held tendency (diff_opt = 1)
      // q_tendf was zeroed by zero_tendf and holds the physics tendencies (update_phy_ten)
      if (diff) {
        const bool is_qv = n == s.P_QV;  // vertical_diffusion_mp diffuses qv - qv_base
        parallel_for(mr, kern::KDiffM{q2.cview(), s.alt.cview(), b.qv_base.view(), s.muts.cview(), c, tf.view(),
                                      p_.khdif * p_.pr_inv, p_.kvdif * p_.pr_inv, real(b.rdx), real(b.rdy), g.nz,
                                      do_y, false, is_qv});
      }
      if (dd_) dd_->scalar(s, b, n, tf, dt_rk);  // diff_opt 1 / sixth order (rk_scalar_tend, rk_step 1)
      if (bc_ && n == s.P_QV) bc_->relax_scalar(tf, q2, s, b);  // relax_bdy_scalar (rk 1)
    }
    if (bc_ && n == s.P_QV) bc_->spec_scalar(tf);  // spec_bdy_scalar (every stage)
    Fluxes fl{q2.cview(), q1.cview(), s.ru_m.cview(), s.rv_m.cview(), s.ww_m.cview(), s.muts.cview(), c,
              real(b.rdx), real(b.rdy), dt_rk, hup, vup, g.nz, do_y, b.msftx.view(), b.msfty.view(),
              kern::Spec{g.bc_x == BC::Specified, g.bc_y == BC::Specified, g.nx, g.ny}};
    if (!pd) {
      parallel_for(mr, KUpdate{fl, q1.cview(), tf.cview(), s.mu_1.cview(), s.mu_2.cview(), b.mub.view(),
                               q_new_.view(), dt_rk, true});
    } else {
      parallel_for(mr, KPreUpdatePD{q1.view(), tf.cview(), s.mu_1.cview(), b.mub.view(), c, dt_rk});
      fill_halo(q1, g);
      parallel_for(Range3D{0, g.nx + 1, 0, g.ny, 0, g.nz}, KFacePDX{fl, dhx_.view(), lox_.view()});
      if (do_y) parallel_for(Range3D{0, g.nx, 0, g.ny + 1, 0, g.nz}, KFacePDY{fl, dhy_.view(), loy_.view()});
      parallel_for(Range3D{0, g.nx, 0, g.ny, 0, g.nz + 1}, KFacePDZ{fl, dhz_.view(), loz_.view()});
      const FacesPD fa{dhx_.cview(), lox_.cview(), dhy_.cview(), loy_.cview(), dhz_.cview(), loz_.cview()};
      parallel_for(mr, KScalePD{fl, fa, s.mu_1.cview(), b.mub.view(), scale_.view()});
      fill_halo(scale_, g);
      parallel_for(mr, KUpdatePD{fl, fa, scale_.cview(), q1.cview(), s.mu_1.cview(), s.mu_2.cview(), b.mub.view(),
                                 q_new_.view(), dt_rk});
    }
    swap(q2, q_new_);
    if (bc_ && moist && n != s.P_QV) bc_->flow_dep(q2, s);  // have_bcs_moist = .false.
  }
  std::vector<Field3*> qs;  // halos of the new scalars, batched (one launch per 8 fields)
  for (int n = 0; n < ns; ++n) qs.push_back(&s.q_2[n]);
  fill_halos(qs, g);
}

}  // namespace cppwrf::dyn
