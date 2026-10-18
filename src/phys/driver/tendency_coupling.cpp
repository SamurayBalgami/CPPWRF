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

// phy_prep_part2, use_theta_m = 1: RTHFTEN of theta_m -> RTHFTEN of dry theta
struct KDryThetaForcing {
  V3 rthften;
  C3 rqvften, t_new, qv;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    using namespace wrfc;
    const real th_phy = (t_new(i, j, k) + t0) / (real(1.) + (r_v / r_d) * qv(i, j, k));
    rthften(i, j, k) =
        th_phy / (t_new(i, j, k) + t0) * (rthften(i, j, k) - (r_v / r_d) * th_phy * rqvften(i, j, k));
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
      rqiblten(g, Stagger::M, "rqiblten"), rucuten(g, Stagger::M, "rucuten"), rvcuten(g, Stagger::M, "rvcuten"),
      rthcuten(g, Stagger::M, "rthcuten"), rqvcuten(g, Stagger::M, "rqvcuten"), rqccuten(g, Stagger::M, "rqccuten"),
      rqrcuten(g, Stagger::M, "rqrcuten"), rqicuten(g, Stagger::M, "rqicuten"), rqscuten(g, Stagger::M, "rqscuten"),
      rthften(g, Stagger::M, "rthften"), rqvften(g, Stagger::M, "rqvften") {}

void calculate_phy_tend(PhysicsTendencies& pt, const dyn::State& s, const dyn::BaseState& b, const Grid& g,
                        const CouplingFlags& f) {
  const Range3D mr = g.mass_cells();
  auto couple = [&](Field3& x) { parallel_for(mr, KCouple{x.view(), s.mut.cview(), b.c1h.view(), b.c2h.view()}); };
  if (f.radiation) couple(pt.rthraten);
  if (f.cu_physics > 0) {
    couple(pt.rucuten);
    couple(pt.rvcuten);
    couple(pt.rthcuten);
    couple(pt.rqvcuten);
    if (s.P_QC >= 0) couple(pt.rqccuten);
    if (s.P_QR >= 0) couple(pt.rqrcuten);
    if (s.P_QI >= 0) couple(pt.rqicuten);
    if (s.P_QS >= 0) couple(pt.rqscuten);
  }
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
  if (f.cu_physics == 1) {  // phy_cu_ten, KFETASCHEME
    parallel_for(ar, KAddA2A{t.t_tendf.view(), pt.rthcuten.cview()});
    if (s.P_QV >= 0) parallel_for(ar, KAddA2A{t.q_tendf[s.P_QV].view(), pt.rqvcuten.cview()});
    if (s.P_QC >= 0) parallel_for(ar, KAddA2A{t.q_tendf[s.P_QC].view(), pt.rqccuten.cview()});
    if (s.P_QR >= 0) parallel_for(ar, KAddA2A{t.q_tendf[s.P_QR].view(), pt.rqrcuten.cview()});
    if (s.P_QI >= 0) parallel_for(ar, KAddA2A{t.q_tendf[s.P_QI].view(), pt.rqicuten.cview()});
    if (s.P_QS >= 0) parallel_for(ar, KAddA2A{t.q_tendf[s.P_QS].view(), pt.rqscuten.cview()});
  }
  if (f.cu_physics == 3) {  // phy_cu_ten, GFLSCHEME: no momentum (RUCUTEN/RVCUTEN are not added)
    parallel_for(ar, KAddA2A{t.t_tendf.view(), pt.rthcuten.cview()});
    if (s.P_QV >= 0) parallel_for(ar, KAddA2A{t.q_tendf[s.P_QV].view(), pt.rqvcuten.cview()});
    if (s.P_QC >= 0) parallel_for(ar, KAddA2A{t.q_tendf[s.P_QC].view(), pt.rqccuten.cview()});
    if (s.P_QI >= 0) parallel_for(ar, KAddA2A{t.q_tendf[s.P_QI].view(), pt.rqicuten.cview()});
  }
  if (f.cu_physics == 16) {  // phy_cu_ten, NTIEDTKESCHEME
    if (!g.specified()) fill_halos({&pt.rucuten, &pt.rvcuten}, g);
    parallel_for(a2c_u_range(g), KAddA2CU{t.ru_tendf.view(), pt.rucuten.cview()});
    parallel_for(a2c_v_range(g), KAddA2CV{t.rv_tendf.view(), pt.rvcuten.cview()});
    parallel_for(ar, KAddA2A{t.t_tendf.view(), pt.rthcuten.cview()});
    if (s.P_QV >= 0) parallel_for(ar, KAddA2A{t.q_tendf[s.P_QV].view(), pt.rqvcuten.cview()});
    if (s.P_QC >= 0) parallel_for(ar, KAddA2A{t.q_tendf[s.P_QC].view(), pt.rqccuten.cview()});
    if (s.P_QI >= 0) parallel_for(ar, KAddA2A{t.q_tendf[s.P_QI].view(), pt.rqicuten.cview()});
  }
}

void conv_t_tendf_to_moist(dyn::Tendencies& t, const Field3& t_1, const dyn::State& s, const Grid& g) {
  parallel_for(g.mass_cells(),
               KConvTMoist{t.t_tendf.view(), t_1.cview(), s.q_1[s.P_QV].cview(), t.q_tendf[s.P_QV].cview()});
}

void phy_prep_part2(PhysicsTendencies& pt, const dyn::State& s, const dyn::BaseState& b, const Grid& g,
                    const CouplingFlags& f) {
  const Range3D mr = g.mass_cells();
  auto decouple = [&](Field3& x) {
    parallel_for(mr, KDecouple{x.view(), s.muts.cview(), b.c1h.view(), b.c2h.view()});
  };
  if (f.radiation) decouple(pt.rthraten);
  if (f.cu_physics > 0) {
    decouple(pt.rucuten);
    decouple(pt.rvcuten);
    decouple(pt.rthcuten);
    if (s.P_QV >= 0) decouple(pt.rqvcuten);
    if (s.P_QC >= 0) decouple(pt.rqccuten);
    if (s.P_QR >= 0) decouple(pt.rqrcuten);
    if (s.P_QI >= 0) decouple(pt.rqicuten);
    if (s.P_QS >= 0) decouple(pt.rqscuten);
  }
  if (f.pbl) {
    decouple(pt.rublten);
    decouple(pt.rvblten);
    decouple(pt.rthblten);
    if (s.P_QV >= 0) decouple(pt.rqvblten);
    if (s.P_QC >= 0) decouple(pt.rqcblten);
    if (s.P_QI >= 0) decouple(pt.rqiblten);
  }
  if (f.cu_physics == 1 && s.P_QV >= 0) decouple(pt.rqvften);  // phy_prep_part2: KFETASCHEME (RQVFTEN only)
  if (f.cu_physics == 16 || f.cu_physics == 3) {  // advective forcing of the cumulus scheme
    if (s.P_QV >= 0) decouple(pt.rqvften);
    decouple(pt.rthften);
    if (f.theta_m)
      parallel_for(mr, KDryThetaForcing{pt.rthften.view(), pt.rqvften.cview(), s.t_2.cview(), s.q_2[s.P_QV].cview()});
  }
}

}  // namespace cppwrf::phys::driver
