#include "phys/driver/phy_prep.hpp"

#include <cmath>
#include "core/backend.hpp"
#include "core/halo.hpp"
#include "phys/phys_constants.hpp"

namespace cppwrf::phys::driver {

namespace {

using V3 = View3<real>;
using C3 = View3<const real>;
using V2 = View2<real>;
using C2 = View2<const real>;
using B3 = BView3<const real>;
using B2 = BView2<const real>;
using A1 = View1<const double>;

constexpr real G = real(cst::g);
constexpr real T0 = real(cst::t0);
constexpr real P1000 = real(cst::p1000mb);
constexpr real RCP = real(cst::rcp);
constexpr real RVRD = real(cst::rvovrd);

// 1. Mass-point diagnostic quantities: th_phy, th_phy_m_t0, p_phy, pi_phy, t_phy, rho, u_phy, v_phy
struct KMassPointDiags {
  C3 t, p, alt, u, v, qv;
  B3 pb;
  V3 th_phy, th_phy_m_t0, p_phy, pi_phy, t_phy, rho, u_phy, v_phy;
  bool use_theta_m, has_qv;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    const real q = has_qv ? qv(i, j, k) : real(0);
    real th;
    if (use_theta_m && has_qv) {
      th = (t(i, j, k) + T0) / (real(1) + RVRD * q);
    } else {
      th = t(i, j, k) + T0;
    }
    th_phy(i, j, k) = th;
    th_phy_m_t0(i, j, k) = th - T0;

    const real p_full = p(i, j, k) + pb(i, j, k);
    p_phy(i, j, k) = p_full;

    const real pi = std::pow(p_full / P1000, RCP);
    pi_phy(i, j, k) = pi;
    t_phy(i, j, k) = th * pi;

    rho(i, j, k) = (real(1) + q) / alt(i, j, k);

    u_phy(i, j, k) = real(0.5) * (u(i, j, k) + u(i + 1, j, k));
    v_phy(i, j, k) = real(0.5) * (v(i, j, k) + v(i, j + 1, k));
  }
};

// 2. W-point diagnostic quantities: z_at_w
struct KWPointDiags {
  C3 ph;
  B3 phb;
  V3 z_at_w;

  CPPWRF_HD void operator()(int i, int j, int k) const {
    z_at_w(i, j, k) = (phb(i, j, k) + ph(i, j, k)) / G;
  }
};

// 3. Column-vertical operations: dz8w, z, p8w, t8w, and hydrostatic pressure p_hyd_w, p_hyd
struct KColumnVertical {
  C3 z_at_w, p_phy, t_phy, qtot;
  C2 mut;
  V3 dz8w, z, p8w, t8w, p_hyd, p_hyd_w;
  A1 fnm, fnp, c1h, c2h, dnw;
  real p_top;
  int nz;
  bool has_qtot;

  CPPWRF_HD void operator()(int i, int j) const {
    // A. dz8w and mass-level z
    for (int k = 0; k < nz; ++k) {
      const real zw0 = z_at_w(i, j, k);
      const real zw1 = z_at_w(i, j, k + 1);
      dz8w(i, j, k) = zw1 - zw0;
      z(i, j, k) = real(0.5) * (zw0 + zw1);
    }
    dz8w(i, j, nz) = real(0);

    // B. Interior p8w and t8w (k = 1 .. nz - 1)
    for (int k = 1; k < nz; ++k) {
      const real fm = real(fnm(k));
      const real fp = real(fnp(k));
      p8w(i, j, k) = fm * p_phy(i, j, k) + fp * p_phy(i, j, k - 1);
      t8w(i, j, k) = fm * t_phy(i, j, k) + fp * t_phy(i, j, k - 1);
    }

    // C. Bottom surface extrapolation (k = 0)
    {
      const real z0 = z_at_w(i, j, 0);
      const real z1 = z(i, j, 0);
      const real z2 = z(i, j, 1);
      const real w1 = (z0 - z2) / (z1 - z2);
      const real w2 = real(1) - w1;
      p8w(i, j, 0) = w1 * p_phy(i, j, 0) + w2 * p_phy(i, j, 1);
      t8w(i, j, 0) = w1 * t_phy(i, j, 0) + w2 * t_phy(i, j, 1);
    }

    // D. Top extrapolation (k = nz)
    {
      const real z0 = z_at_w(i, j, nz);
      const real z1 = z(i, j, nz - 1);
      const real z2 = z(i, j, nz - 2);
      const real w1 = (z0 - z2) / (z1 - z2);
      const real w2 = real(1) - w1;
      p8w(i, j, nz) = std::exp(w1 * std::log(p_phy(i, j, nz - 1)) + w2 * std::log(p_phy(i, j, nz - 2)));
      t8w(i, j, nz) = w1 * t_phy(i, j, nz - 1) + w2 * t_phy(i, j, nz - 2);
    }

    // E. Hydrostatic pressure calculation (full levels, top-down)
    p_hyd_w(i, j, nz) = p_top;
    const real mu_ij = mut(i, j);
    for (int k = nz - 1; k >= 0; --k) {
      const real qt = has_qtot ? qtot(i, j, k) : real(0);
      const real col_mass = real(c1h(k)) * mu_ij + real(c2h(k));
      const real dn = real(dnw(k));
      // WRF: p_hyd_w(k) = p_hyd_w(k+1) - (1. + qtot)*(c1(k)*MUT(i,j) + c2(k))*dnw(k)
      p_hyd_w(i, j, k) = p_hyd_w(i, j, k + 1) - (real(1) + qt) * col_mass * dn;
    }

    // F. Hydrostatic pressure at half levels (mass points)
    for (int k = 0; k < nz; ++k) {
      p_hyd(i, j, k) = real(0.5) * (p_hyd_w(i, j, k) + p_hyd_w(i, j, k + 1));
    }
  }
};

}  // namespace

PhyPrepFields::PhyPrepFields(const Grid& g)
    : th_phy(g, Stagger::M, "th_phy"),
      th_phy_m_t0(g, Stagger::M, "th_phy_m_t0"),
      p_phy(g, Stagger::M, "p_phy"),
      pi_phy(g, Stagger::M, "pi_phy"),
      t_phy(g, Stagger::M, "t_phy"),
      rho(g, Stagger::M, "rho"),
      u_phy(g, Stagger::M, "u_phy"),
      v_phy(g, Stagger::M, "v_phy"),
      z(g, Stagger::M, "z"),
      p_hyd(g, Stagger::M, "p_hyd"),
      z_at_w(g, Stagger::W, "z_at_w"),
      dz8w(g, Stagger::W, "dz8w"),
      p8w(g, Stagger::W, "p8w"),
      t8w(g, Stagger::W, "t8w"),
      p_hyd_w(g, Stagger::W, "p_hyd_w") {}

void run_phy_prep(const dyn::State& s, const dyn::BaseState& b,
                  PhyPrepFields& p, const PhyPrepParams& params) {
  const Grid g = p.th_phy.layout().ny == 1 ? Grid{s.u_2.layout().nx - 1, 1, s.t_2.layout().nz}
                                           : Grid{s.u_2.layout().nx - 1, s.v_2.layout().ny - 1, s.t_2.layout().nz};

  const bool has_qv = (s.P_QV >= 0 && s.P_QV < int(s.q_2.size()));
  C3 qv_view = has_qv ? s.q_2[s.P_QV].cview() : s.t_2.cview();

  // 1. Diagnostics at mass points
  parallel_for(
      g.mass_cells(),
      KMassPointDiags{s.t_2.cview(), s.p.cview(), s.alt.cview(), s.u_2.cview(), s.v_2.cview(), qv_view,
                      b.pb.view(), p.th_phy.view(), p.th_phy_m_t0.view(), p.p_phy.view(), p.pi_phy.view(),
                      p.t_phy.view(), p.rho.view(), p.u_phy.view(), p.v_phy.view(),
                      params.use_theta_m == 1, has_qv});

  // 2. Diagnostics at w points
  parallel_for(g.cells(Stagger::W), KWPointDiags{s.ph_2.cview(), b.phb.view(), p.z_at_w.view()});

  // 3. Vertical column operations
  const bool has_qtot = (s.n_moist > 0);
  C3 qtot_view = has_qtot ? s.qtot.cview() : s.t_2.cview();

  column_for(
      g.columns(),
      KColumnVertical{p.z_at_w.cview(), p.p_phy.cview(), p.t_phy.cview(), qtot_view,
                      s.mut.cview(), p.dz8w.view(), p.z.view(), p.p8w.view(), p.t8w.view(),
                      p.p_hyd.view(), p.p_hyd_w.view(),
                      b.fnm.view(), b.fnp.view(), b.c1h.view(), b.c2h.view(), b.dnw.view(),
                      real(b.p_top), g.nz, has_qtot});

  // Halos
  fill_halo(p.th_phy, g);
  fill_halo(p.th_phy_m_t0, g);
  fill_halo(p.p_phy, g);
  fill_halo(p.pi_phy, g);
  fill_halo(p.t_phy, g);
  fill_halo(p.rho, g);
  fill_halo(p.u_phy, g);
  fill_halo(p.v_phy, g);
  fill_halo(p.z, g);
  fill_halo(p.p_hyd, g);
  fill_halo(p.z_at_w, g);
  fill_halo(p.dz8w, g);
  fill_halo(p.p8w, g);
  fill_halo(p.t8w, g);
  fill_halo(p.p_hyd_w, g);
}

}  // namespace cppwrf::phys::driver
