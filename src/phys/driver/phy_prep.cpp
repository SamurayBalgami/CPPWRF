#include "phys/driver/phy_prep.hpp"

#include <cmath>
#include <stdexcept>

#include "phys/driver/wrf_constants.hpp"

namespace cppwrf::phys::driver {

namespace {

constexpr int kMaxMoist = 12;

// phy_prep (module_big_step_utilities_em.F), one column per thread; each WRF loop nest is a
// loop over k here (the loops are independent across columns).
struct KPhyPrep {
  View3<const real> u, v, p, alt, ph, t;
  BView3<const real> pb, phb;
  View3<const real> q[kMaxMoist];
  int n_moist, p_qv;
  View2<const real> mut;
  View1<const real> c1h, c2h, fnm, fnp, dnw;
  real p_top;
  int nz, use_theta_m;
  View3<real> th_phy, th_phy_m_t0, p_phy, pi_phy, t_phy, rho, u_phy, v_phy, z, z_at_w, dz8w, p8w, t8w, p_hyd,
      p_hyd_w;

  CPPWRF_HD void operator()(int i, int j) const {
    using namespace wrfc;
    const int kte = nz;  // WRF kte = kde (w levels 1..kde -> 0..nz), k_end = kde-1 -> nz-1
    if (use_theta_m == 1 && p_qv >= 0) {
      for (int k = 0; k < nz; ++k) th_phy(i, j, k) = (t(i, j, k) + t0) / (real(1.) + r_v / r_d * q[p_qv](i, j, k));
    } else {
      for (int k = 0; k < nz; ++k) th_phy(i, j, k) = t(i, j, k) + t0;
    }
    for (int k = 0; k < nz; ++k) {
      th_phy_m_t0(i, j, k) = th_phy(i, j, k) - t0;
      p_phy(i, j, k) = p(i, j, k) + pb(i, j, k);
      pi_phy(i, j, k) = rpow(p_phy(i, j, k) / p1000mb, rcp);
      t_phy(i, j, k) = th_phy(i, j, k) * pi_phy(i, j, k);
      rho(i, j, k) = real(1.) / alt(i, j, k) * (real(1.) + q[p_qv](i, j, k));
      u_phy(i, j, k) = real(0.5) * (u(i, j, k) + u(i + 1, j, k));
      v_phy(i, j, k) = real(0.5) * (v(i, j, k) + v(i, j + 1, k));
    }
    for (int k = 0; k <= kte; ++k) z_at_w(i, j, k) = (phb(i, j, k) + ph(i, j, k)) / g;
    for (int k = 0; k <= kte - 1; ++k) dz8w(i, j, k) = z_at_w(i, j, k + 1) - z_at_w(i, j, k);
    dz8w(i, j, kte) = real(0.);
    for (int k = 0; k < nz; ++k) z(i, j, k) = real(0.5) * (z_at_w(i, j, k) + z_at_w(i, j, k + 1));
    for (int k = 1; k < nz; ++k) {
      p8w(i, j, k) = fnm(k) * p_phy(i, j, k) + fnp(k) * p_phy(i, j, k - 1);
      t8w(i, j, k) = fnm(k) * t_phy(i, j, k) + fnp(k) * t_phy(i, j, k - 1);
    }
    {
      real z0 = z_at_w(i, j, 0);
      real z1 = z(i, j, 0);
      real z2 = z(i, j, 1);
      real w1 = (z0 - z2) / (z1 - z2);
      real w2 = real(1.) - w1;
      p8w(i, j, 0) = w1 * p_phy(i, j, 0) + w2 * p_phy(i, j, 1);
      t8w(i, j, 0) = w1 * t_phy(i, j, 0) + w2 * t_phy(i, j, 1);
      z0 = z_at_w(i, j, kte);
      z1 = z(i, j, nz - 1);
      z2 = z(i, j, nz - 2);
      w1 = (z0 - z2) / (z1 - z2);
      w2 = real(1.) - w1;
      p8w(i, j, nz) = std::exp(w1 * rlog(p_phy(i, j, nz - 1)) + w2 * rlog(p_phy(i, j, nz - 2)));
      t8w(i, j, nz) = w1 * t_phy(i, j, nz - 1) + w2 * t_phy(i, j, nz - 2);
    }
    p_hyd_w(i, j, kte) = p_top;
    for (int k = kte - 1; k >= 0; --k) {
      real qtot = real(0.);
      for (int n = 0; n < n_moist; ++n) qtot = qtot + q[n](i, j, k);
      p_hyd_w(i, j, k) = p_hyd_w(i, j, k + 1) - (real(1.) + qtot) * (c1h(k) * mut(i, j) + c2h(k)) * dnw(k);
    }
    for (int k = 0; k < nz; ++k) p_hyd(i, j, k) = real(0.5) * (p_hyd_w(i, j, k) + p_hyd_w(i, j, k + 1));
  }
};

}  // namespace

PhyPrepFields::PhyPrepFields(const Grid& g)
    : th_phy(g, Stagger::M, "th_phy"), th_phy_m_t0(g, Stagger::M, "th_phy_m_t0"), p_phy(g, Stagger::M, "p_phy"),
      pi_phy(g, Stagger::M, "pi_phy"), t_phy(g, Stagger::M, "t_phy"), rho(g, Stagger::M, "rho"),
      u_phy(g, Stagger::M, "u_phy"), v_phy(g, Stagger::M, "v_phy"), z(g, Stagger::M, "z"),
      p_hyd(g, Stagger::M, "p_hyd"), z_at_w(g, Stagger::W, "z_at_w"), dz8w(g, Stagger::W, "dz8w"),
      p8w(g, Stagger::W, "p8w"), t8w(g, Stagger::W, "t8w"), p_hyd_w(g, Stagger::W, "p_hyd_w") {}

void run_phy_prep(const dyn::State& s, const dyn::BaseState& b, PhyPrepFields& p, const PhyPrepParams& params) {
  if (s.P_QV < 0) throw std::runtime_error("phy_prep: the physics suite needs qv");
  if (s.n_moist > kMaxMoist) throw std::runtime_error("phy_prep: too many moist species");
  KPhyPrep k{};
  k.u = s.u_2.cview(); k.v = s.v_2.cview(); k.p = s.p.cview(); k.alt = s.alt.cview();
  k.ph = s.ph_2.cview(); k.t = s.t_2.cview();
  k.pb = b.pb.view(); k.phb = b.phb.view();
  for (int n = 0; n < s.n_moist; ++n) k.q[n] = s.q_2[n].cview();
  k.n_moist = s.n_moist; k.p_qv = s.P_QV;
  k.mut = s.mut.cview();
  k.c1h = b.c1h.view(); k.c2h = b.c2h.view(); k.fnm = b.fnm.view(); k.fnp = b.fnp.view(); k.dnw = b.dnw.view();
  k.p_top = real(b.p_top);
  k.nz = s.t_2.layout().nz;
  k.use_theta_m = params.use_theta_m;
  k.th_phy = p.th_phy.view(); k.th_phy_m_t0 = p.th_phy_m_t0.view(); k.p_phy = p.p_phy.view();
  k.pi_phy = p.pi_phy.view(); k.t_phy = p.t_phy.view(); k.rho = p.rho.view(); k.u_phy = p.u_phy.view();
  k.v_phy = p.v_phy.view(); k.z = p.z.view(); k.z_at_w = p.z_at_w.view(); k.dz8w = p.dz8w.view();
  k.p8w = p.p8w.view(); k.t8w = p.t8w.view(); k.p_hyd = p.p_hyd.view(); k.p_hyd_w = p.p_hyd_w.view();
  column_for(Range2D{0, s.t_2.layout().nx, 0, s.t_2.layout().ny}, k);
}

}  // namespace cppwrf::phys::driver
