#include "dyn/physics_coupling.hpp"

#include <stdexcept>

#include "core/backend.hpp"
#include "core/constants.hpp"
#include "core/halo.hpp"
#include "phys/mp_kessler/kessler.hpp"
#include "phys/phys_constants.hpp"

// WRF dyn_em/module_big_step_utilities_em.F: moist_physics_prep_em, moist_physics_finish_em
// (use_theta_m = 1 branch, no DFI, no WRFPLUS). C++ k = Fortran k - 1.

namespace cppwrf::dyn {
namespace {

using V3 = View3<real>;
using C3 = View3<const real>;
using V2 = View2<real>;
using B3 = BView3<const real>;
constexpr real RVRD = real(cst::rvovrd);
constexpr real T0 = real(cst::t0);

// moist_physics_prep_em: z_at_w, dz8w, th_phy (dry theta), h_diabatic = th_phy, qv_diabatic,
// rho, pii, z (mass-level height)
struct KPrep {
  C3 t, qv, ph, al, p; B3 phb, alb, pb;
  V3 th_phy, hd, qvd, rho, pii, z, dz8w;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const real zw0 = (ph(i, j, k) + phb(i, j, k)) / real(cst::g);
    const real zw1 = (ph(i, j, k + 1) + phb(i, j, k + 1)) / real(cst::g);
    dz8w(i, j, k) = zw1 - zw0;
    z(i, j, k) = real(0.5) * (zw0 + zw1);
    const real q = qv(i, j, k);
    const real th = (t(i, j, k) + T0) / (real(1) + RVRD * q);
    th_phy(i, j, k) = th;
    hd(i, j, k) = th;
    qvd(i, j, k) = q;
    rho(i, j, k) = real(1) / (al(i, j, k) + alb(i, j, k));
    pii(i, j, k) = real(pow(double(p(i, j, k) + pb(i, j, k)) / cst::p1000mb, cst::rcp));
  }
};

// moist_physics_finish_em, use_theta_m = 1: new theta_m - t0 and the heating rate h_diabatic
struct KFinish {
  V3 t, hd; C3 th_phy, qv, qvd; real dt, lim;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const real th = th_phy(i, j, k), th_old = hd(i, j, k), q = qv(i, j, k), q_old = qvd(i, j, k);
    real mpten = th - th_old;
    mpten = mpten < lim * dt ? mpten : lim * dt;   // min(mp_tend_lim*dt, mpten)
    mpten = mpten > -lim * dt ? mpten : -lim * dt; // max(-mp_tend_lim*dt, mpten)
    const real qvten = q - q_old;
    t(i, j, k) = th_old * (real(1) + RVRD * q_old) + mpten * (real(1) + RVRD * q) + RVRD * qvten * th - T0;
    hd(i, j, k) = (mpten * (real(1) + RVRD * q) + RVRD * qvten * th) / dt;
  }
};

}  // namespace

PhysicsCoupling::PhysicsCoupling(const Grid& g, const PhysicsParams& p)
    : rainnc(g, Stagger::M, "rainnc"), rainncv(g, Stagger::M, "rainncv"), g_(g), p_(p),
      th_phy(g, Stagger::M, "th_phy"), rho(g, Stagger::M, "rho"), pii(g, Stagger::M, "pii"), z(g, Stagger::M, "z"),
      dz8w(g, Stagger::M, "dz8w"), qv_diabatic(g, Stagger::M, "qv_diabatic") {
  if (p_.mp_physics != 0 && p_.mp_physics != 1) throw std::runtime_error("phys.mp_physics must be 0 or 1");
  if (p_.mp_physics == 1 && g.nz > phys::kessler::kMaxLevels) throw std::runtime_error("Kessler: nz > kMaxLevels");
}

void PhysicsCoupling::step(State& s, const BaseState& b, real dt) {
  if (p_.mp_physics == 0) return;
  if (s.P_QV < 0 || s.P_QC < 0 || s.P_QR < 0) throw std::runtime_error("mp_physics = 1 needs qv, qc, qr");
  const Grid& g = g_;
  Field3& qv = s.q_2[s.P_QV];
  parallel_for(g.mass_cells(), KPrep{s.t_2.cview(), qv.cview(), s.ph_2.cview(), s.al.cview(), s.p.cview(),
                                     b.phb.view(), b.alb.view(), b.pb.view(), th_phy.view(), s.h_diabatic.view(),
                                     qv_diabatic.view(), rho.view(), pii.view(), z.view(), dz8w.view()});
  namespace pc = phys::cst;
  column_for(g.columns(),
             phys::kessler::KesslerKernel{th_phy.view(), qv.view(), s.q_2[s.P_QC].view(), s.q_2[s.P_QR].view(),
                                          rho.cview(), pii.cview(), dz8w.cview(), z.cview(), rainnc.view(),
                                          rainncv.view(), dt, real(pc::xlv), real(pc::cp), real(pc::ep_2),
                                          real(pc::svp1), real(pc::svp2), real(pc::svp3), real(pc::svpt0),
                                          real(pc::rhowater), g.nz});
  parallel_for(g.mass_cells(), KFinish{s.t_2.view(), s.h_diabatic.view(), th_phy.cview(), qv.cview(),
                                       qv_diabatic.cview(), dt, p_.mp_tend_lim});
  fill_halo(s.t_2, g);
  fill_halo(s.h_diabatic, g);
  for (int n = 0; n < s.n_moist; ++n) fill_halo(s.q_2[n], g);
}

}  // namespace cppwrf::dyn
