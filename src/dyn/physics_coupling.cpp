// SPDX-License-Identifier: GPL-3.0-or-later
#include "dyn/physics_coupling.hpp"

#include <map>
#include <mutex>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "core/backend.hpp"
#include "core/constants.hpp"
#include "core/halo.hpp"
#include "core/limits.hpp"
#include "core/vm_dump.hpp"
#include "phys/mp_kessler/kessler.hpp"
#include "phys/mp_thompson/thompson.hpp"
#include "phys/mp_wsm6/wsm6.hpp"
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
  V3 th_phy, hd, qvd, rho, pii, z, dz8w, p_phy;
  int i0, ni;  // WRF's i loop (its..ite = i0..i0+ni-1): gfortran pairs it into libmvec's 2-lane pow (cppwrf::vm)
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
    p_phy(i, j, k) = p(i, j, k) + pb(i, j, k);
    pii(i, j, k) = vm::pow((p(i, j, k) + pb(i, j, k)) / real(cst::p1000mb), real(cst::rcp), i - i0, ni);
  }
};

// moist_physics_finish_em, use_theta_m = 1: new theta_m - t0 and the heating rate h_diabatic
struct KFinish {
  V3 t, hd; C3 th_phy, qv; V3 qvd; C3 qc; V3 qcd; real dt, lim; bool has_qc;
  CPPWRF_HD void operator()(int i, int j, int k) const {
    const real th = th_phy(i, j, k), th_old = hd(i, j, k), q = qv(i, j, k), q_old = qvd(i, j, k);
    real mpten = th - th_old;
    mpten = mpten < lim * dt ? mpten : lim * dt;   // min(mp_tend_lim*dt, mpten)
    mpten = mpten > -lim * dt ? mpten : -lim * dt; // max(-mp_tend_lim*dt, mpten)
    const real qvten = q - q_old;
    t(i, j, k) = th_old * (real(1) + RVRD * q_old) + mpten * (real(1) + RVRD * q) + RVRD * qvten * th - T0;
    hd(i, j, k) = (mpten * (real(1) + RVRD * q) + RVRD * qvten * th) / dt;
    qvd(i, j, k) = qvten / dt;  // qv_diabatic, qc_diabatic: the microphysics tendencies of qv and qc
    if (has_qc) qcd(i, j, k) = (qc(i, j, k) - qcd(i, j, k)) / dt;
  }
};

struct KSaveQc {  // moist_physics_prep_em: qc_diabatic = qc
  C3 qc; V3 qcd;
  CPPWRF_HD void operator()(int i, int j, int k) const { qcd(i, j, k) = qc(i, j, k); }
};

// module_mp_thompson.F thompson_init computes (or reads) the tables only while tcg_racg is not allocated: the
// first domain's call makes them and every nest uses that one set (~380 MB in FP64)
std::shared_ptr<const phys::mp_thompson::ThompsonDeviceTables<real>> thompson_tables(const std::string& dir,
                                                                                     bool force_read, bool aero) {
  using Tables = phys::mp_thompson::ThompsonDeviceTables<real>;
  static std::mutex m;
  static std::map<std::tuple<std::string, bool, bool>, std::weak_ptr<const Tables>> made;
  std::lock_guard<std::mutex> lock(m);
  std::weak_ptr<const Tables>& w = made[{dir, force_read, aero}];
  if (auto t = w.lock()) return t;
  phys::mp_thompson::ThompsonHostTables<real> h;
  phys::mp_thompson::init_thompson_tables(h, dir, force_read, aero);
  auto t = std::make_shared<Tables>();
  t->upload(h);
  w = t;
  return t;
}

}  // namespace

PhysicsCoupling::PhysicsCoupling(const Grid& g, const PhysicsParams& p)
    : rainnc(g, Stagger::M, "rainnc"), rainncv(g, Stagger::M, "rainncv"), g_(g), p_(p),
      th_phy(g, Stagger::M, "th_phy"), rho(g, Stagger::M, "rho"), pii(g, Stagger::M, "pii"), z(g, Stagger::M, "z"),
      dz8w(g, Stagger::M, "dz8w"), qv_diabatic(g, Stagger::M, "qv_diabatic"), qc_diabatic(g, Stagger::M, "qc_diabatic"),
      p_phy(g, Stagger::M, "p_phy") {
  if (p_.mp_physics != 0 && p_.mp_physics != 1 && p_.mp_physics != 6 && p_.mp_physics != 8 && p_.mp_physics != 28)
    throw std::runtime_error("phys.mp_physics must be 0, 1, 6, 8 or 28");
  if (p_.mp_physics == 1) check_limit("Kessler", "nz", g.nz, phys::kessler::kMaxLevels);
  if (p_.do_radar_ref != 0 && p_.do_radar_ref != 1) throw std::runtime_error("phys.do_radar_ref must be 0 or 1");
  if (p_.microp_re) {  // phy_init: background values RE_QC_BG, RE_QI_BG, RE_QS_BG
    re_cloud = Field3(g, Stagger::M, "re_cloud");
    re_ice = Field3(g, Stagger::M, "re_ice");
    re_snow = Field3(g, Stagger::M, "re_snow");
    re_cloud.fill_host([](int, int, int) { return real(2.49E-6); });
    re_ice.fill_host([](int, int, int) { return real(4.99E-6); });
    re_snow.fill_host([](int, int, int) { return real(9.99E-6); });
  }
  if (p_.mp_physics == 6) {
    check_limit("WSM6", "nz", g.nz, phys::wsm6::kMaxLevels);
    snownc = Field2(g, Stagger::M, "snownc");
    snowncv = Field2(g, Stagger::M, "snowncv");
    graupelnc = Field2(g, Stagger::M, "graupelnc");
    graupelncv = Field2(g, Stagger::M, "graupelncv");
    sr = Field2(g, Stagger::M, "sr");
    namespace pc = phys::cst;
    wp_.g = real(cst::g); wp_.cpd = real(pc::cp); wp_.cpv = real(4.0 * pc::r_v); wp_.rd = real(pc::r_d);
    wp_.rv = real(pc::r_v); wp_.t0c = real(pc::svpt0); wp_.ep1 = real(pc::ep_1); wp_.ep2 = real(pc::ep_2);
    wp_.qmin = real(1e-15); wp_.xls = real(pc::xls); wp_.xlv0 = real(pc::xlv); wp_.xlf0 = real(pc::xlf);
    wp_.cliq = real(pc::cliq); wp_.cice = real(pc::cice); wp_.psat = real(610.78);
    // wsm6init (WRF microphysics_driver arguments); once, not every step
    wp_.init(real(pc::rhoair0), real(pc::rhowater), real(pc::rhosnow), real(pc::cliq), real(4.0 * pc::r_v),
             p_.hail_opt);
#if defined(CPPWRF_DOUBLE)
    {  // the work arrays of wsm6_stages for step()'s mp_cols
      const int sx = g.bc_x == BC::Specified ? 1 : 0, sy = g.bc_y == BC::Specified ? 1 : 0;
      wsm6_work_ = Buffer(phys::wsm6::Wsm6Work<real>::size(Range2D{sx, g.nx - sx, sy, g.ny - sy}, g.nz));
    }
#endif
    if (p_.do_radar_ref == 1)  // mp_wsm6_init: m(D) = pi*den/6 D**3 of rain, snow and graupel (hail), mu = 0
      radar_ = std::make_shared<const phys::radar::RadarDeviceTables<real>>(phys::radar::radar_init<real>(
          wp_.pi * wp_.denr / real(6.), real(3.), real(0.), wp_.pi * wp_.dens / real(6.), real(3.), real(0.),
          wp_.pi * wp_.deng / real(6.), real(3.), real(0.)));
  }
  if (p_.do_radar_ref == 1 && (p_.mp_physics == 6 || thompson())) refl_10cm = Field3(g, Stagger::M, "refl_10cm");
  if (thompson()) {
    check_limit("Thompson", "nz", g.nz, phys::mp_thompson::kMaxLevels);
    snownc = Field2(g, Stagger::M, "snownc");
    snowncv = Field2(g, Stagger::M, "snowncv");
    graupelnc = Field2(g, Stagger::M, "graupelnc");
    graupelncv = Field2(g, Stagger::M, "graupelncv");
    sr = Field2(g, Stagger::M, "sr");
    if (p_.mp_physics == 28) {  // the input's QNWFA2D / QNIFA2D (real.exe with aer_init_opt 0: zero)
      qnwfa2d = Field2(g, Stagger::M, "qnwfa2d");
      qnifa2d = Field2(g, Stagger::M, "qnifa2d");
    }
    tp_ = thompson_tables(p_.table_dir, p_.force_read_thompson, p_.mp_physics == 28);
    {  // the column arrays' work buffer on the GPU, for step()'s mp_cols
      const int sx = g.bc_x == BC::Specified ? 1 : 0, sy = g.bc_y == BC::Specified ? 1 : 0;
      phys::mp_thompson::thompson_reserve<real>(Range2D{sx, g.nx - sx, sy, g.ny - sy}, g.nz);
    }
    if (p_.do_radar_ref == 1) {  // thompson_init (calc_refl10cm reads only the bins, K_w, PI5, lamda4, m_w_0, m_i_0)
      namespace mt = phys::mp_thompson;
      radar_ = std::make_shared<const phys::radar::RadarDeviceTables<real>>(phys::radar::radar_init<real>(
          real(mt::am_r), real(mt::bm_r), real(mt::mu_r), real(mt::am_s), real(mt::bm_s), real(mt::mu_s),
          real(mt::am_g5), real(mt::bm_g), real(mt::mu_g)));
    }
  }
}

void PhysicsCoupling::step(State& s, const BaseState& b, real dt) {
  if (p_.mp_physics == 0) return;
  if (s.P_QV < 0 || s.P_QC < 0 || s.P_QR < 0) throw std::runtime_error("microphysics needs qv, qc, qr");
  if (p_.mp_physics == 6 && (s.P_QI < 0 || s.P_QS < 0 || s.P_QG < 0))
    throw std::runtime_error("mp_physics = 6 (WSM6) needs qv, qc, qr, qi, qs, qg");
  if (thompson() && (s.P_QI < 0 || s.P_QS < 0 || s.P_QG < 0 || s.P_QNI < 0 || s.P_QNR < 0))
    throw std::runtime_error("mp_physics = 8 / 28 (Thompson) needs qv, qc, qr, qi, qs, qg, qni, qnr");
  if (p_.mp_physics == 28 && (s.P_QNC < 0 || s.P_QNWFA < 0 || s.P_QNIFA < 0 || s.P_QNBCA < 0))
    throw std::runtime_error("mp_physics = 28 (Thompson aerosol aware) also needs qnc, qnwfa, qnifa, qnbca");
  const Grid& g = g_;
  Field3& qv = s.q_2[s.P_QV];
  // solve_em calls moist_physics_prep_em on the tile without the spec zone (sz) in x (unless periodic_x)
  // and y; only the lanes of its i loop matter here (the mp columns all lie inside it)
  const int pz = g.bc_x == BC::Specified || g.bc_y == BC::Specified ? 1 : 0;
  const int pi0 = g.bc_x == BC::Periodic ? 0 : pz;
  parallel_for(g.mass_cells(), KPrep{s.t_2.cview(), qv.cview(), s.ph_2.cview(), s.al.cview(), s.p.cview(),
                                     b.phb.view(), b.alb.view(), b.pb.view(), th_phy.view(), s.h_diabatic.view(),
                                     qv_diabatic.view(), rho.view(), pii.view(), z.view(), dz8w.view(), p_phy.view(), pi0, g.nx - 2 * pi0});
  if (s.P_QC >= 0) parallel_for(g.mass_cells(), KSaveQc{s.q_2[s.P_QC].cview(), qc_diabatic.view()});
  namespace pc = phys::cst;
  // microphysics_driver: with specified boundaries the spec zone is skipped (sz = spec_zone)
  const int sx = g.bc_x == BC::Specified ? 1 : 0, sy = g.bc_y == BC::Specified ? 1 : 0;
  const Range2D mp_cols{sx, g.nx - sx, sy, g.ny - sy};
  if (p_.mp_physics == 6) {
    // debug (CPPWRF_VM_DUMP): WSM6 inputs and outputs, for a replay through WRF's mp_wsm6
    vmdump::dump("vm_mp_in", 1, {&th_phy, &qv, &s.q_2[s.P_QC], &s.q_2[s.P_QR], &s.q_2[s.P_QI], &s.q_2[s.P_QS],
                                 &s.q_2[s.P_QG], &rho, &pii, &p_phy, &dz8w},
                 {&rainnc, &rainncv, &snownc, &snowncv, &graupelnc, &graupelncv, &sr});
    // the column arrays sized to the next tier above nz (CUDA reserves a kernel's stack frame for every resident
    // thread of the GPU: FP64 49 KB per thread with 128 levels, 19 KB with 48)
    auto wsm6 = [&](auto nl) {
#if defined(CPPWRF_DOUBLE)  // FP64 builds: the level loops in parallel, the sedimentation per column (wsm6.hpp)
      phys::wsm6::wsm6_stages<real, decltype(nl)::value>(
          wp_, dt, g.nz, mp_cols, wsm6_work_.get(), th_phy.view(), qv.view(), s.q_2[s.P_QC].view(),
          s.q_2[s.P_QR].view(), s.q_2[s.P_QI].view(), s.q_2[s.P_QS].view(), s.q_2[s.P_QG].view(), rho.cview(),
          pii.cview(), p_phy.cview(), dz8w.cview(), rainnc.view(), rainncv.view(), snownc.view(), snowncv.view(),
          graupelnc.view(), graupelncv.view(), sr.view());
#else
      column_for(mp_cols, phys::wsm6::Wsm6Kernel<real, decltype(nl)::value>{
                              wp_, dt, g.nz, th_phy.view(), qv.view(), s.q_2[s.P_QC].view(), s.q_2[s.P_QR].view(),
                              s.q_2[s.P_QI].view(), s.q_2[s.P_QS].view(), s.q_2[s.P_QG].view(), rho.cview(),
                              pii.cview(), p_phy.cview(), dz8w.cview(), rainnc.view(), rainncv.view(), snownc.view(),
                              snowncv.view(), graupelnc.view(), graupelncv.view(), sr.view()});
#endif
    };
    if (g.nz <= 48) wsm6(std::integral_constant<int, 48>{});
    else if (g.nz <= 64) wsm6(std::integral_constant<int, 64>{});
    else wsm6(std::integral_constant<int, phys::wsm6::kMaxLevels>{});
    if (p_.microp_re)
      parallel_for(Range3D{mp_cols.i0, mp_cols.i1, mp_cols.j0, mp_cols.j1, 0, g.nz},
                   phys::wsm6::Wsm6EffRadKernel<real>{
                       wp_.pidnc, wp_.pidn0s, wp_.alpha, wp_.n0s, wp_.n0smax, wp_.t0c, wp_.qmin, real(2.49E-6),
                       real(4.99E-6), real(9.99E-6), real(50.E-6), real(125.E-6), real(999.E-6), th_phy.cview(),
                       pii.cview(), s.q_2[s.P_QC].cview(), s.q_2[s.P_QI].cview(), s.q_2[s.P_QS].cview(), rho.cview(),
                       re_cloud.view(), re_ice.view(), re_snow.view()});
    if (radar_ && diagflag) {  // module_mp_wsm6.F: refl10cm_wsm6 after the scheme (diagflag, do_radar_ref = 1)
      auto refl = [&](auto nl) {
        column_for(mp_cols, phys::wsm6::Wsm6ReflKernel<real, decltype(nl)::value>{
                                radar_->get(), wp_.n0r, wp_.n0s, wp_.n0smax, wp_.alpha, wp_.n0g, wp_.pi, g.nz,
                                th_phy.cview(), pii.cview(), p_phy.cview(), qv.cview(), s.q_2[s.P_QR].cview(),
                                s.q_2[s.P_QS].cview(), s.q_2[s.P_QG].cview(), refl_10cm.view()});
      };
      if (g.nz <= 48) refl(std::integral_constant<int, 48>{});
      else refl(std::integral_constant<int, phys::wsm6::kMaxLevels>{});
    }
    vmdump::dump("vm_mp_out", 1, {&th_phy, &qv, &s.q_2[s.P_QC], &s.q_2[s.P_QR], &s.q_2[s.P_QI], &s.q_2[s.P_QS],
                                  &s.q_2[s.P_QG]},
                 {&rainnc, &rainncv, &snownc, &snowncv, &graupelnc, &graupelncv, &sr});
  } else if (thompson()) {
    // debug (CPPWRF_VM_DUMP): Thompson inputs and outputs, for a replay through WRF's mp_gt_driver
    vmdump::dump("vm_mp_in", 1, {&th_phy, &qv, &s.q_2[s.P_QC], &s.q_2[s.P_QR], &s.q_2[s.P_QI], &s.q_2[s.P_QS],
                                 &s.q_2[s.P_QG], &s.q_2[s.P_QNI], &s.q_2[s.P_QNR], &pii, &p_phy, &s.w_2, &dz8w},
                 {&rainnc, &rainncv, &snownc, &snowncv, &graupelnc, &graupelncv, &sr});
    phys::mp_thompson::ThompsonAeroFields<real> aero{};
    if (p_.mp_physics == 28)
      aero = {s.q_2[s.P_QNC].view(), s.q_2[s.P_QNWFA].view(), s.q_2[s.P_QNIFA].view(), s.q_2[s.P_QNBCA].view(),
              qnwfa2d.cview(), qnifa2d.cview()};
    phys::mp_thompson::run_thompson<real>(
        th_phy.view(), qv.view(), s.q_2[s.P_QC].view(), s.q_2[s.P_QR].view(),
        s.q_2[s.P_QI].view(), s.q_2[s.P_QS].view(), s.q_2[s.P_QG].view(),
        s.q_2[s.P_QNI].view(), s.q_2[s.P_QNR].view(), pii.cview(), p_phy.cview(),
        s.w_2.cview(), dz8w.cview(), rainnc.view(), rainncv.view(), snownc.view(),
        snowncv.view(), graupelnc.view(), graupelncv.view(), sr.view(),
        dt, mp_cols, g.nz, tp_->view(), p_.microp_re ? re_cloud.view() : View3<real>{},
        p_.microp_re ? re_ice.view() : View3<real>{}, p_.microp_re ? re_snow.view() : View3<real>{},
        p_.mp_physics == 28 ? &aero : nullptr, radar_ && diagflag ? refl_10cm.view() : View3<real>{},
        radar_ ? radar_->get() : nullptr);
    vmdump::dump("vm_mp_out", 1, {&th_phy, &qv, &s.q_2[s.P_QC], &s.q_2[s.P_QR], &s.q_2[s.P_QI], &s.q_2[s.P_QS],
                                  &s.q_2[s.P_QG], &s.q_2[s.P_QNI], &s.q_2[s.P_QNR]},
                 {&rainnc, &rainncv, &snownc, &snowncv, &graupelnc, &graupelncv, &sr});
  } else {
    column_for(mp_cols,
               phys::kessler::KesslerKernel{th_phy.view(), qv.view(), s.q_2[s.P_QC].view(), s.q_2[s.P_QR].view(),
                                            rho.cview(), pii.cview(), dz8w.cview(), z.cview(), rainnc.view(),
                                            rainncv.view(), dt, real(pc::xlv), real(pc::cp), real(pc::ep_2),
                                            real(pc::svp1), real(pc::svp2), real(pc::svp3), real(pc::svpt0),
                                            real(pc::rhowater), g.nz, mp_cols.i0, mp_cols.i1 - mp_cols.i0});
  }
  parallel_for(g.mass_cells(), KFinish{s.t_2.view(), s.h_diabatic.view(), th_phy.cview(), qv.cview(),
                                       qv_diabatic.view(), s.P_QC >= 0 ? s.q_2[s.P_QC].cview() : qv.cview(),
                                       qc_diabatic.view(), dt, p_.mp_tend_lim, s.P_QC >= 0});
  std::vector<Field3*> fs{&s.t_2, &s.h_diabatic};
  for (int n = 0; n < s.n_moist; ++n) fs.push_back(&s.q_2[n]);
  for (int n : {s.P_QNI, s.P_QNR, s.P_QNC, s.P_QNWFA, s.P_QNIFA, s.P_QNBCA})
    if (n >= 0) fs.push_back(&s.q_2[n]);
  fill_halos(fs, g);
}

// module_mp_thompson.F thompson_init (is_aerosol_aware, wif_input_opt /= 2), on the host as a setup step. hgt is
// start_em's z_at_q = (ph_2 + phb)/g at the w levels k = 1..kde of the whole patch; the profile loops k = 2..kte
// (kte = kde) are vectorised by gfortran (libmvec exp in pairs from k = 2); C++ k = Fortran k - 1, the level kde
// (WRF's unused top memory level) is only a lane of that loop.
void PhysicsCoupling::init_aerosols(State& s, const BaseState& b) {
  if (p_.mp_physics != 28) return;
  const Grid& g = g_;
  const int nx = g.nx, ny = g.ny, nz = g.nz;
  std::vector<real> ph, phb, nw, w2d;
  s.ph_2.download(ph);
  b.phb.download(phb);
  const Layout &LW = s.ph_2.layout(), &LB = b.phb.layout(), &LM = s.q_2[s.P_QNWFA].layout();
  const Layout& L2 = qnwfa2d.layout();
  auto hgt = [&](int i, int k, int j) {
    return (double(ph[LW.index(i, j, k)]) + double(phb[LB.index(i, j, k)])) / cst::g;
  };
  // log(naCCN1/naCCN0), log(naIN1/naIN0): constants that gfortran folds correctly rounded
  const double lccn = double(logl((long double)(phys::mp_thompson::naCCN1 / phys::mp_thompson::naCCN0)));
  const double lin = double(logl((long double)(phys::mp_thompson::naIN1 / phys::mp_thompson::naIN0)));
  auto profile = [&](int n, double na0, double na1, double lna, bool emission) {
    Field3& f = s.q_2[n];
    f.download(nw);
    double max_test = double(nw[LM.index(0, 0, 0)]);
    for (int k = 0; k < nz; ++k)  // MAXVAL(nwfa(its:ite-1,:,jts:jte-1)): every mass point (WRF's kde level is 0)
      for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) max_test = std::max(max_test, double(nw[LM.index(i, j, k)]));
    if (!(max_test < phys::mp_thompson::eps_val)) return;
    if (emission) qnwfa2d.download(w2d);
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        const double h1 = hgt(i, 0, j);
        double h_01;
        if (h1 <= 1000.0) h_01 = 0.8;
        else if (h1 >= 2500.0) h_01 = 0.01;
        else h_01 = 0.8 * libm::cos(h1 * 0.001 - 1.0);
        const double ni3 = -1.0 * lna / h_01;
        nw[LM.index(i, j, 0)] = real(na1 + na0 * libm::exp(-((hgt(i, 1, j) - h1) / 1000.) * ni3));
        if (emission) {
          const double z1 = hgt(i, 1, j) - h1;
          w2d[L2.index(i, j, 0)] = real(double(nw[LM.index(i, j, 0)]) * 0.000196 * (50. / z1));
        }
        for (int k = 1; k <= nz; ++k) {  // Fortran k = 2..kde
          const double v = na1 + na0 * vm::exp(-((hgt(i, k, j) - h1) / 1000.) * ni3, k - 1, nz);
          if (k < nz) nw[LM.index(i, j, k)] = real(v);
        }
      }
    f.fill_host([&](int i, int j, int k) { return nw[LM.index(i, j, k)]; });
    fill_halo(f, g);
    s.q_1[n].copy_from(f);
    if (emission) qnwfa2d.fill_host([&](int i, int j, int) { return w2d[L2.index(i, j, 0)]; });
  };
  profile(s.P_QNWFA, phys::mp_thompson::naCCN0, phys::mp_thompson::naCCN1, lccn, true);
  profile(s.P_QNIFA, phys::mp_thompson::naIN0, phys::mp_thompson::naIN1, lin, false);
}

std::vector<const Field3*> PhysicsCoupling::fields() const {
  std::vector<const Field3*> f{&rho, &qv_diabatic, &qc_diabatic};
  if (refl_10cm.layout().total > 0) f.push_back(&refl_10cm);
  return f;
}

}  // namespace cppwrf::dyn
