#include "phys/driver/tendency_coupling.hpp"

#include "core/backend.hpp"
#include "core/halo.hpp"

namespace cppwrf::phys::driver {

namespace {

using V3 = View3<real>;
using C3 = View3<const real>;
using V2 = View2<real>;
using C2 = View2<const real>;
using B2 = BView2<const real>;
using A1 = View1<const double>;

struct KAddA2A {
  V3 lvar;
  C3 rvar;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    lvar(i, j, k) += rvar(i, j, k);
  }
};

struct KAddA2CU {
  V3 lvar;
  C3 rvar;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    lvar(i, j, k) += real(0.5) * (rvar(i, j, k) + rvar(i - 1, j, k));
  }
};

struct KAddA2CV {
  V3 lvar;
  C3 rvar;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    lvar(i, j, k) += real(0.5) * (rvar(i, j, k) + rvar(i, j - 1, k));
  }
};

struct KMuCouple {
  V3 var;
  C2 mut;
  A1 c1h, c2h;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const real col_mass = real(c1h(k)) * mut(i, j) + real(c2h(k));
    var(i, j, k) *= col_mass;
  }
};

struct KMapFactorU {
  V3 ru_tend;
  C3 ru_tendf;
  B2 msfuy;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    ru_tend(i, j, k) += ru_tendf(i, j, k) / msfuy(i, j);
  }
};

struct KMapFactorV {
  V3 rv_tend;
  C3 rv_tendf;
  B2 msfvx_inv;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    rv_tend(i, j, k) += rv_tendf(i, j, k) * msfvx_inv(i, j);
  }
};

struct KMapFactorT {
  V3 t_tend;
  C3 t_tendf;
  B2 msfty;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    t_tend(i, j, k) += t_tendf(i, j, k) / msfty(i, j);
  }
};

}  // namespace

PhysicsTendencies::PhysicsTendencies(const Grid& g)
    : rthraten(g, Stagger::M, "rthraten"),
      rthblten(g, Stagger::M, "rthblten"),
      rublten(g, Stagger::M, "rublten"),
      rvblten(g, Stagger::M, "rvblten"),
      rqvblten(g, Stagger::M, "rqvblten"),
      rqcblten(g, Stagger::M, "rqcblten"),
      rqiblten(g, Stagger::M, "rqiblten"),
      rthcuten(g, Stagger::M, "rthcuten"),
      rucuten(g, Stagger::M, "rucuten"),
      rvcuten(g, Stagger::M, "rvcuten"),
      rqvcuten(g, Stagger::M, "rqvcuten"),
      rqccuten(g, Stagger::M, "rqccuten"),
      rqrcuten(g, Stagger::M, "rqrcuten"),
      rqicuten(g, Stagger::M, "rqicuten"),
      rqscuten(g, Stagger::M, "rqscuten"),
      rthshten(g, Stagger::M, "rthshten"),
      rushten(g, Stagger::M, "rushten"),
      rvshten(g, Stagger::M, "rvshten"),
      rqvshten(g, Stagger::M, "rqvshten"),
      rqcshten(g, Stagger::M, "rqcshten"),
      rqrshten(g, Stagger::M, "rqrshten"),
      rqishten(g, Stagger::M, "rqishten"),
      rqsshten(g, Stagger::M, "rqsshten"),
      rqgshten(g, Stagger::M, "rqgshten") {}

void add_a2a(View3<real> lvar, View3<const real> rvar, Range3D range) {
  parallel_for(range, KAddA2A{lvar, rvar});
}

void add_a2a(Field3& lvar, const Field3& rvar, const Grid& g) {
  add_a2a(lvar.view(), rvar.cview(), g.mass_cells());
}

void add_a2c_u(View3<real> lvar, View3<const real> rvar, Range3D range_u) {
  parallel_for(range_u, KAddA2CU{lvar, rvar});
}

void add_a2c_u(Field3& lvar, const Field3& rvar, const Grid& g) {
  const Range3D ur{g.u_i0(), g.u_i1(), 0, g.ny, 0, g.nz};
  add_a2c_u(lvar.view(), rvar.cview(), ur);
}

void add_a2c_v(View3<real> lvar, View3<const real> rvar, Range3D range_v) {
  parallel_for(range_v, KAddA2CV{lvar, rvar});
}

void add_a2c_v(Field3& lvar, const Field3& rvar, const Grid& g) {
  const Range3D vr{0, g.nx, g.v_j0(), g.v_j1(), 0, g.nz};
  add_a2c_v(lvar.view(), rvar.cview(), vr);
}

void calculate_phy_tend(PhysicsTendencies& pt, const dyn::State& s, const dyn::BaseState& b, const Grid& g) {
  const Range3D mr = g.mass_cells();
  Field3* tends[] = {
      &pt.rthraten, &pt.rthblten, &pt.rublten,  &pt.rvblten,  &pt.rqvblten, &pt.rqcblten,
      &pt.rqiblten, &pt.rthcuten, &pt.rucuten,  &pt.rvcuten,  &pt.rqvcuten, &pt.rqccuten,
      &pt.rqrcuten, &pt.rqicuten, &pt.rqscuten, &pt.rthshten, &pt.rushten,  &pt.rvshten,
      &pt.rqvshten, &pt.rqcshten, &pt.rqrshten, &pt.rqishten, &pt.rqsshten, &pt.rqgshten};

  for (Field3* t : tends) {
    parallel_for(mr, KMuCouple{t->view(), s.mut.cview(), b.c1h.view(), b.c2h.view()});
    fill_halo(*t, g);
  }
}

void update_phy_ten(dyn::Tendencies& t, const PhysicsTendencies& pt, const dyn::State& s, const Grid& g) {
  const Range3D mr = g.mass_cells();
  const Range3D ur{g.u_i0(), g.u_i1(), 0, g.ny, 0, g.nz};
  const Range3D vr{0, g.nx, g.v_j0(), g.v_j1(), 0, g.nz};

  // Temperature / Theta tendency: RTHRATEN, RTHBLTEN, RTHCUTEN, RTHSHTEN
  add_a2a(t.t_tendf.view(), pt.rthraten.cview(), mr);
  add_a2a(t.t_tendf.view(), pt.rthblten.cview(), mr);
  add_a2a(t.t_tendf.view(), pt.rthcuten.cview(), mr);
  add_a2a(t.t_tendf.view(), pt.rthshten.cview(), mr);

  // Momentum U tendency: RUBLTEN, RUCUTEN, RUSHTEN
  add_a2c_u(t.ru_tendf.view(), pt.rublten.cview(), ur);
  add_a2c_u(t.ru_tendf.view(), pt.rucuten.cview(), ur);
  add_a2c_u(t.ru_tendf.view(), pt.rushten.cview(), ur);

  // Momentum V tendency: RVBLTEN, RVCUTEN, RVSHTEN
  add_a2c_v(t.rv_tendf.view(), pt.rvblten.cview(), vr);
  add_a2c_v(t.rv_tendf.view(), pt.rvcuten.cview(), vr);
  add_a2c_v(t.rv_tendf.view(), pt.rvshten.cview(), vr);

  // Moisture species tendencies into t.q_tendf
  auto add_scalar_tend = [&](int idx, const Field3& tend) {
    if (idx >= 0 && idx < int(t.q_tendf.size())) {
      add_a2a(t.q_tendf[idx].view(), tend.cview(), mr);
    }
  };

  add_scalar_tend(s.P_QV, pt.rqvblten);
  add_scalar_tend(s.P_QV, pt.rqvcuten);
  add_scalar_tend(s.P_QV, pt.rqvshten);

  add_scalar_tend(s.P_QC, pt.rqcblten);
  add_scalar_tend(s.P_QC, pt.rqccuten);
  add_scalar_tend(s.P_QC, pt.rqcshten);

  add_scalar_tend(s.P_QR, pt.rqrcuten);
  add_scalar_tend(s.P_QR, pt.rqrshten);

  add_scalar_tend(s.P_QI, pt.rqiblten);
  add_scalar_tend(s.P_QI, pt.rqicuten);
  add_scalar_tend(s.P_QI, pt.rqishten);

  add_scalar_tend(s.P_QS, pt.rqscuten);
  add_scalar_tend(s.P_QS, pt.rqsshten);

  add_scalar_tend(s.P_QG, pt.rqgshten);

  fill_halo(t.t_tendf, g);
  fill_halo(t.ru_tendf, g);
  fill_halo(t.rv_tendf, g);
  for (Field3& qf : t.q_tendf) fill_halo(qf, g);
}

void apply_tendf_map_factors(dyn::Tendencies& t, const dyn::BaseState& b, const Grid& g) {
  const Range3D ur{g.u_i0(), g.u_i1(), 0, g.ny, 0, g.nz};
  const Range3D vr{0, g.nx, g.v_j0(), g.v_j1(), 0, g.nz};
  const Range3D mr = g.mass_cells();

  parallel_for(ur, KMapFactorU{t.ru_tend.view(), t.ru_tendf.cview(), b.msfuy.view()});
  parallel_for(vr, KMapFactorV{t.rv_tend.view(), t.rv_tendf.cview(), b.msfvx_inv.view()});
  parallel_for(mr, KMapFactorT{t.t_tend.view(), t.t_tendf.cview(), b.msfty.view()});

  fill_halo(t.ru_tend, g);
  fill_halo(t.rv_tend, g);
  fill_halo(t.t_tend, g);
}

void couple_physics_tendencies(dyn::Tendencies& t, PhysicsTendencies& pt,
                               const dyn::State& s, const dyn::BaseState& b, const Grid& g,
                               bool apply_map_factors) {
  calculate_phy_tend(pt, s, b, g);
  update_phy_ten(t, pt, s, g);
  if (apply_map_factors) {
    apply_tendf_map_factors(t, b, g);
  }
}

}  // namespace cppwrf::phys::driver
