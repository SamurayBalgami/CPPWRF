// SPDX-License-Identifier: GPL-3.0-or-later
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
// add_a2a / add_a2c_u / add_a2c_v of one tendency (the kernels above), for a batched launch
struct KAdd {
  enum Kind : int { A2A, A2CU, A2CV };
  V3 lvar;
  C3 rvar;
  Kind kind;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    if (kind == A2A) KAddA2A{lvar, rvar}(i, j, k);
    else if (kind == A2CU) KAddA2CU{lvar, rvar}(i, j, k);
    else KAddA2CV{lvar, rvar}(i, j, k);
  }
};
// The additions of update_phy_ten in launches of independent additions: an addition goes into the launch after the
// last one that adds to the same tendency, so every tendency gets its terms in WRF's order (the sources are only
// read); additions to different tendencies are independent.
class AddLaunches {
 public:
  void add(const Range3D& r, KAdd::Kind kind, Field3& lvar, const Field3& rvar) {
    if constexpr (!kFuseMulAdd) {  // CUDA FP32: each addition as its own kernel, in WRF's order (backend.hpp)
      if (kind == KAdd::A2A) parallel_for(r, KAddA2A{lvar.view(), rvar.cview()});
      else if (kind == KAdd::A2CU) parallel_for(r, KAddA2CU{lvar.view(), rvar.cview()});
      else parallel_for(r, KAddA2CV{lvar.view(), rvar.cview()});
      return;
    }
    int w = 0;
    for (const Op& o : ops_)
      if (o.k.lvar.p == lvar.data() && o.wave >= w) w = o.wave + 1;
    ops_.push_back(Op{r, KAdd{lvar.view(), rvar.cview(), kind}, w});
    if (w + 1 > nwave_) nwave_ = w + 1;
  }
  void run() {
    for (int w = 0; w < nwave_; ++w) {
      LaunchBatch<KAdd, 16> b;
      for (const Op& o : ops_)
        if (o.wave == w) b.add(o.r, o.k);
      b.run();
    }
  }

 private:
  struct Op {
    Range3D r;
    KAdd k;
    int wave;
  };
  std::vector<Op> ops_;
  int nwave_ = 0;
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
  LaunchBatch<KCouple, 16> batch;  // every tendency once: independent, one launch
  auto couple = [&](Field3& x) { batch.add(mr, KCouple{x.view(), s.mut.cview(), b.c1h.view(), b.c2h.view()}); };
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
  batch.run();
}

void update_phy_ten(dyn::Tendencies& t, PhysicsTendencies& pt, const dyn::State& s, const Grid& g,
                    const CouplingFlags& f) {
  const Range3D ar = a2a_range(g);
  // periodic: i-1 / j-1 neighbours of the momentum tendencies (before the additions, which only read them)
  if (f.pbl && !g.specified()) fill_halos({&pt.rublten, &pt.rvblten}, g);
  if ((f.cu_physics == 16 || f.cu_physics == 6) && !g.specified()) fill_halos({&pt.rucuten, &pt.rvcuten}, g);
  AddLaunches add;
  if (f.radiation) add.add(ar, KAdd::A2A, t.t_tendf, pt.rthraten);  // phy_ra_ten
  if (f.pbl) {  // phy_bl_ten, YSUSCHEME
    add.add(ar, KAdd::A2A, t.t_tendf, pt.rthblten);
    add.add(a2c_u_range(g), KAdd::A2CU, t.ru_tendf, pt.rublten);
    add.add(a2c_v_range(g), KAdd::A2CV, t.rv_tendf, pt.rvblten);
    if (s.P_QV >= 0) add.add(ar, KAdd::A2A, t.q_tendf[s.P_QV], pt.rqvblten);
    if (s.P_QC >= 0) add.add(ar, KAdd::A2A, t.q_tendf[s.P_QC], pt.rqcblten);
    if (s.P_QI >= 0) add.add(ar, KAdd::A2A, t.q_tendf[s.P_QI], pt.rqiblten);
  }
  if (f.cu_physics == 1) {  // phy_cu_ten, KFETASCHEME
    add.add(ar, KAdd::A2A, t.t_tendf, pt.rthcuten);
    if (s.P_QV >= 0) add.add(ar, KAdd::A2A, t.q_tendf[s.P_QV], pt.rqvcuten);
    if (s.P_QC >= 0) add.add(ar, KAdd::A2A, t.q_tendf[s.P_QC], pt.rqccuten);
    if (s.P_QR >= 0) add.add(ar, KAdd::A2A, t.q_tendf[s.P_QR], pt.rqrcuten);
    if (s.P_QI >= 0) add.add(ar, KAdd::A2A, t.q_tendf[s.P_QI], pt.rqicuten);
    if (s.P_QS >= 0) add.add(ar, KAdd::A2A, t.q_tendf[s.P_QS], pt.rqscuten);
  }
  if (f.cu_physics == 3) {  // phy_cu_ten, GFLSCHEME: no momentum (RUCUTEN/RVCUTEN are not added)
    add.add(ar, KAdd::A2A, t.t_tendf, pt.rthcuten);
    if (s.P_QV >= 0) add.add(ar, KAdd::A2A, t.q_tendf[s.P_QV], pt.rqvcuten);
    if (s.P_QC >= 0) add.add(ar, KAdd::A2A, t.q_tendf[s.P_QC], pt.rqccuten);
    if (s.P_QI >= 0) add.add(ar, KAdd::A2A, t.q_tendf[s.P_QI], pt.rqicuten);
  }
  if (f.cu_physics == 16 || f.cu_physics == 6) {  // phy_cu_ten, TIEDTKESCHEME and NTIEDTKESCHEME
    add.add(a2c_u_range(g), KAdd::A2CU, t.ru_tendf, pt.rucuten);
    add.add(a2c_v_range(g), KAdd::A2CV, t.rv_tendf, pt.rvcuten);
    add.add(ar, KAdd::A2A, t.t_tendf, pt.rthcuten);
    if (s.P_QV >= 0) add.add(ar, KAdd::A2A, t.q_tendf[s.P_QV], pt.rqvcuten);
    if (s.P_QC >= 0) add.add(ar, KAdd::A2A, t.q_tendf[s.P_QC], pt.rqccuten);
    if (s.P_QI >= 0) add.add(ar, KAdd::A2A, t.q_tendf[s.P_QI], pt.rqicuten);
  }
  add.run();
}

void conv_t_tendf_to_moist(dyn::Tendencies& t, const Field3& t_1, const dyn::State& s, const Grid& g) {
  parallel_for(g.mass_cells(),
               KConvTMoist{t.t_tendf.view(), t_1.cview(), s.q_1[s.P_QV].cview(), t.q_tendf[s.P_QV].cview()});
}

void phy_prep_part2(PhysicsTendencies& pt, const dyn::State& s, const dyn::BaseState& b, const Grid& g,
                    const CouplingFlags& f) {
  const Range3D mr = g.mass_cells();
  LaunchBatch<KDecouple, 16> batch;  // every tendency once: independent, one launch per 16
  auto decouple = [&](Field3& x) { batch.add(mr, KDecouple{x.view(), s.muts.cview(), b.c1h.view(), b.c2h.view()}); };
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
  if ((f.cu_physics == 1 || f.cu_physics == 6) && s.P_QV >= 0) decouple(pt.rqvften);  // KFETASCHEME, TIEDTKESCHEME: RQVFTEN only
  if (f.cu_physics == 16 || f.cu_physics == 3) {  // advective forcing of the cumulus scheme
    if (s.P_QV >= 0) decouple(pt.rqvften);
    decouple(pt.rthften);
  }
  batch.run();
  if (f.cu_physics == 16 || f.cu_physics == 3) {
    if (f.theta_m)
      parallel_for(mr, KDryThetaForcing{pt.rthften.view(), pt.rqvften.cview(), s.t_2.cview(), s.q_2[s.P_QV].cview()});
  }
}

}  // namespace cppwrf::phys::driver
