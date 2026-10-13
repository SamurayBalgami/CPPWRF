#include "phys/driver/tendency_coupling.hpp"

#include "core/halo.hpp"
#include "phys/driver/wrf_constants.hpp"

namespace cppwrf::phys::driver {

namespace {

using V3 = View3<real>;
using C3 = View3<const real>;
using C2 = View2<const real>;
using A1 = View1<const real>;

// ten = (c1(k)*mut + c2(k)) * ten
struct KCouple {
  V3 ten;
  C2 mut;
  A1 c1, c2;
  CPPWRF_HD void operator()(int i, int j, int k) const { ten(i, j, k) = (c1(k) * mut(i, j) + c2(k)) * ten(i, j, k); }
};
// ten = ten / (c1(k)*MUT + c2(k))
struct KDecouple {
  V3 ten;
  C2 mut;
  A1 c1, c2;
  CPPWRF_HD void operator()(int i, int j, int k) const { ten(i, j, k) = ten(i, j, k) / (c1(k) * mut(i, j) + c2(k)); }
};
struct KAddA2A {
  V3 lvar;
  C3 rvar;
  CPPWRF_HD void operator()(int i, int j, int k) const { lvar(i, j, k) = lvar(i, j, k) + rvar(i, j, k); }
};
struct KAddA2CU {
  V3 lvar;
  C3 rvar;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    lvar(i, j, k) = lvar(i, j, k) + real(0.5) * (rvar(i, j, k) + rvar(i - 1, j, k));
  }
};
struct KAddA2CV {
  V3 lvar;
  C3 rvar;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    lvar(i, j, k) = lvar(i, j, k) + real(0.5) * (rvar(i, j, k) + rvar(i, j - 1, k));
  }
};
struct KConvTMoist {
  V3 t_tendf;
  C3 t_1, moist_old, moist_tend;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    using namespace wrfc;
    t_tendf(i, j, k) = (real(1.) + (r_v / r_d) * moist_old(i, j, k)) * t_tendf(i, j, k) +
                       (r_v / r_d) * (t_1(i, j, k) + t0) / (real(1.) + r_v / r_d * moist_old(i, j, k)) *
                           moist_tend(i, j, k);
  }
};

// add_a2a / add_a2c_u / add_a2c_v loop bounds (config_flags%specified: outer row excluded;
// periodic_x: full rows). C++ k = WRF k - 1; ktf = kde-1 -> nz-1.
Range3D a2a_range(const Grid& g) {
  if (g.specified()) return {1, g.nx - 1, 1, g.ny - 1, 0, g.nz};
  return {0, g.nx, 0, g.ny, 0, g.nz};
}
Range3D a2c_u_range(const Grid& g) {
  if (g.specified()) return {1, g.nx, 1, g.ny - 1, 0, g.nz};
  return {0, g.nx + 1, 0, g.ny, 0, g.nz};
}
Range3D a2c_v_range(const Grid& g) {  // WRF loops k = kts..kte here; level kde holds rvblten = 0
  if (g.specified()) return {1, g.nx - 1, 1, g.ny, 0, g.nz};
  return {0, g.nx, 0, g.ny + (g.two_d() ? 0 : 1), 0, g.nz};
}

}  // namespace

PhysicsTendencies::PhysicsTendencies(const Grid& g)
    : rthraten(g, Stagger::M, "rthraten"), rublten(g, Stagger::M, "rublten"), rvblten(g, Stagger::M, "rvblten"),
      rthblten(g, Stagger::M, "rthblten"), rqvblten(g, Stagger::M, "rqvblten"), rqcblten(g, Stagger::M, "rqcblten"),
      rqiblten(g, Stagger::M, "rqiblten") {}

void calculate_phy_tend(PhysicsTendencies& pt, const dyn::State& s, const dyn::BaseState& b, const Grid& g,
                        const CouplingFlags& f) {
  const Range3D mr = g.mass_cells();
  auto couple = [&](Field3& x) { parallel_for(mr, KCouple{x.view(), s.mut.cview(), b.c1h.view(), b.c2h.view()}); };
  if (f.radiation) couple(pt.rthraten);
  if (f.pbl) {
    couple(pt.rublten);
    couple(pt.rvblten);
    couple(pt.rthblten);
    if (s.P_QV >= 0) couple(pt.rqvblten);
    if (s.P_QC >= 0) couple(pt.rqcblten);
    if (s.P_QI >= 0) couple(pt.rqiblten);
  }
}

void update_phy_ten(dyn::Tendencies& t, PhysicsTendencies& pt, const dyn::State& s, const Grid& g,
                    const CouplingFlags& f) {
  const Range3D ar = a2a_range(g);
  if (f.radiation) parallel_for(ar, KAddA2A{t.t_tendf.view(), pt.rthraten.cview()});  // phy_ra_ten
  if (f.pbl) {  // phy_bl_ten, YSUSCHEME
    if (!g.specified()) fill_halos({&pt.rublten, &pt.rvblten}, g);  // periodic: i-1 / j-1 neighbours
    parallel_for(ar, KAddA2A{t.t_tendf.view(), pt.rthblten.cview()});
    parallel_for(a2c_u_range(g), KAddA2CU{t.ru_tendf.view(), pt.rublten.cview()});
    parallel_for(a2c_v_range(g), KAddA2CV{t.rv_tendf.view(), pt.rvblten.cview()});
    if (s.P_QV >= 0) parallel_for(ar, KAddA2A{t.q_tendf[s.P_QV].view(), pt.rqvblten.cview()});
    if (s.P_QC >= 0) parallel_for(ar, KAddA2A{t.q_tendf[s.P_QC].view(), pt.rqcblten.cview()});
    if (s.P_QI >= 0) parallel_for(ar, KAddA2A{t.q_tendf[s.P_QI].view(), pt.rqiblten.cview()});
  }
}

void conv_t_tendf_to_moist(dyn::Tendencies& t, const dyn::State& s, const Grid& g) {
  parallel_for(g.mass_cells(),
               KConvTMoist{t.t_tendf.view(), s.t_1.cview(), s.q_1[s.P_QV].cview(), t.q_tendf[s.P_QV].cview()});
}

void phy_prep_part2(PhysicsTendencies& pt, const dyn::State& s, const dyn::BaseState& b, const Grid& g,
                    const CouplingFlags& f) {
  const Range3D mr = g.mass_cells();
  auto decouple = [&](Field3& x) {
    parallel_for(mr, KDecouple{x.view(), s.muts.cview(), b.c1h.view(), b.c2h.view()});
  };
  if (f.radiation) decouple(pt.rthraten);
  if (f.pbl) {
    decouple(pt.rublten);
    decouple(pt.rvblten);
    decouple(pt.rthblten);
    if (s.P_QV >= 0) decouple(pt.rqvblten);
    if (s.P_QC >= 0) decouple(pt.rqcblten);
    if (s.P_QI >= 0) decouple(pt.rqiblten);
  }
}

}  // namespace cppwrf::phys::driver
