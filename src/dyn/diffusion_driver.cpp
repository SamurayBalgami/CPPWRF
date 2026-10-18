#include "dyn/diffusion_driver.hpp"
#include "core/vm_dump.hpp"

#include <stdexcept>

#include "core/halo.hpp"
#include "dyn/diffusion/diffusion.hpp"
#include "phys/driver/phy_prep.hpp"
#include "phys/driver/suite.hpp"

namespace cppwrf::dyn {

namespace dif = diffusion;
using S = Stagger;

namespace {
Grid plus_one_level(const Grid& g) {
  Grid p = g;
  p.nz = g.nz + 1;
  return p;
}
}  // namespace

DiffusionDriver::DiffusionDriver(const Grid& g, const DiffusionConfig& c)
    : g_(g), gp_(plus_one_level(g)), c_(c),
      z_(g, S::M, "dif_z"), rdz_(g, S::M, "dif_rdz"), rdzw_(g, S::M, "dif_rdzw"),
      zx_(gp_, S::U, "dif_zx"), zy_(gp_, S::V, "dif_zy"),
      d11_(g, S::M, "defor11"), d22_(g, S::M, "defor22"), d33_(g, S::M, "defor33"),
      d12_(g, S::U, "defor12"),  // (nx+1) x (ny+1): row ny lives in the y halo (never halo-filled)
      d13_(gp_, S::U, "defor13"), d23_(gp_, S::V, "defor23"), div_(g, S::M, "div"),
      xkmh_(g, S::M, "xkmh"), xkmv_(g, S::M, "xkmv"), xkhh_(g, S::M, "xkhh"), xkhv_(g, S::M, "xkhv") {
  if (c.diff_opt < 0 || c.diff_opt > 2) throw std::runtime_error("dyn.diff_opt must be 0, 1 or 2");
  if (c.diff_opt != 0 && (c.km_opt < 1 || c.km_opt > 4)) throw std::runtime_error("dyn.km_opt must be 1..4");
  if (c.diff_6th_opt < 0 || c.diff_6th_opt > 2) throw std::runtime_error("dyn.diff_6th_opt must be 0, 1 or 2");
  if (g.nz + 1 >= 128) throw std::runtime_error("diffusion kernels support nz <= 126");
  p_.diff_opt = c.diff_opt;
  p_.km_opt = c.km_opt;
  p_.damp_opt = c.damp_opt;
  p_.zdamp = c.zdamp;
  p_.dampcoef = c.dampcoef;
  p_.khdif = c.khdif;
  p_.kvdif = c.kvdif;
  p_.c_s = c.c_s;
  p_.c_k = c.c_k;
  p_.isotropic = c.mix_isotropic;
  p_.mix_full_fields = c.mix_full_fields;
  p_.mix_upper_bound = c.mix_upper_bound;
  p_.prandtl = c.prandtl;
  p_.tke_drag_coefficient = c.tke_drag_coefficient;
  p_.tke_heat_flux = c.tke_heat_flux;
  p_.tke_upper_bound = c.tke_upper_bound;
  p_.moist_mix2_off = c.moist_mix2_off;
  p_.tke_mix2_off = c.tke_mix2_off;
  p_.tke_mix6_off = c.tke_mix6_off;
  p_.diff_6th_opt = c.diff_6th_opt;
  p_.diff_6th_factor = c.diff_6th_factor;
  p_.bl_pbl_physics = c.pbl ? 1 : 0;
  p_.isfflx = c.isfflx;
  p_.use_theta_m = c.use_theta_m;
  p_.fft_filter_lat = c.fft_filter_lat;
  p_.periodic_x = g.bc_x == BC::Periodic;
  p_.periodic_y = g.bc_y == BC::Periodic;
  p_.specified = g.specified();

  if (c.km_opt == 2 || c.km_opt == 3 || (c.n2_diagnostic && c.diff_opt != 0)) {
    bn2_ = Field3(g, S::M, "bn2");
  }
  if (c.diff_opt == 2) {
    titau11_ = Field3(g, S::M, "titau11");
    titau22_ = Field3(g, S::M, "titau22");
    titau33_ = Field3(g, S::M, "titau33");
    titau12_ = Field3(g, S::U, "titau12");
    titau13_ = Field3(gp_, S::U, "titau13");
    titau23_ = Field3(gp_, S::V, "titau23");
  }
  if (c.diff_opt == 2 || c.km_opt == 2 || c.km_opt == 3 || (c.n2_diagnostic && c.diff_opt != 0)) {
    prep_ = std::make_unique<phys::driver::PhyPrepFields>(g);
    ustm_ = Field2(g, S::M, "dif_ustm");
    hfx_  = Field2(g, S::M, "dif_hfx");
    qfx_  = Field2(g, S::M, "dif_qfx");
  }
}

DiffusionDriver::~DiffusionDriver() = default;

// first_rk_step_part2 (rk_step 1): metrics, deformation, eddy viscosities
void DiffusionDriver::prepare(const State& s, const BaseState& b) {
  if (c_.diff_opt != 1 && c_.diff_opt != 2) return;
  const Grid& g = g_;
  const real rdx = real(b.rdx), rdy = real(b.rdy);
  if (!suite_ && prep_) {
    phys::driver::run_phy_prep(s, b, *prep_, phys::driver::PhyPrepParams{s.P_QV >= 0 ? c_.use_theta_m : 0});
  }
  dif::compute_diff_metrics(p_, s.ph_2.cview(), b.phb.view(), z_.view(), rdz_.view(), rdzw_.view(), zx_.view(),
                            zy_.view(), rdx, rdy, g.nx, g.ny, g.nz, true);
  fill_halos({&z_, &rdz_, &rdzw_}, g);
  fill_halo(zx_, gp_, Parity::Even, Parity::Even);
  fill_halo(zy_, gp_, Parity::Even, Parity::Even);
  if (suite_) {
    fill_halos({&suite_->ustm, &suite_->ust}, g);
  }
  dif::cal_deform_and_div(p_, s.u_2.cview(), s.v_2.cview(), s.w_2.cview(), b.u_base.view(), b.v_base.view(),
                          zx_.cview(), zy_.cview(), rdz_.cview(), rdzw_.cview(), b.msftx.view(), b.msfty.view(),
                          b.msfux.view(), b.msfuy.view(), b.msfvx.view(), b.msfvy.view(), b.fnm.view(), b.fnp.view(),
                          b.dn.view(), b.dnw.view(), real(b.cf1), real(b.cf2), real(b.cf3), rdx, rdy, d11_.view(),
                          d22_.view(), d33_.view(), d12_.view(), d13_.view(), d23_.view(), div_.view(), g.nx, g.ny,
                          g.nz);
  if (c_.km_opt == 1 || c_.km_opt == 4) {
    if (c_.n2_diagnostic) {  // calculate_km_kh: calculate_N2 for every km_opt
      const phys::driver::PhyPrepFields& pr = suite_ ? suite_->prep : *prep_;
      dif::calculate_N2(p_, pr.th_phy.cview(), pr.t_phy.cview(), pr.p_phy.cview(), pr.p8w.cview(), pr.t8w.cview(),
                        rdz_.cview(), rdzw_.cview(), s.P_QV >= 0 ? s.q_2[s.P_QV].cview() : View3<const real>{},
                        s.P_QV >= 0, s.P_QC >= 0 ? s.q_2[s.P_QC].cview() : View3<const real>{}, s.P_QC >= 0,
                        s.P_QI >= 0 ? s.q_2[s.P_QI].cview() : View3<const real>{}, s.P_QI >= 0, real(b.cf1),
                        real(b.cf2), real(b.cf3), bn2_.view(), g.nx, g.ny, g.nz);
    }
    dif::calculate_km_kh(p_, real(g.dx), real(g.dy), d11_.cview(), d22_.cview(), d12_.cview(), rdzw_.cview(),
                         zx_.cview(), zy_.cview(), b.msftx.view(), b.msfty.view(), xkmh_.view(), xkmv_.view(),
                         xkhh_.view(), xkhv_.view(), g.nx, g.ny, g.nz);
    if (c_.damp_opt == 1) {
      dif::cal_dampkm(p_, real(g.dx), real(g.dy), real(g.dt), rdz_.cview(), rdzw_.cview(), b.msftx.view(),
                      b.msfty.view(), xkmh_.view(), xkhh_.view(), xkmv_.view(), xkhv_.view(), g.nx, g.ny, g.nz);
    }
  } else {
    const phys::driver::PhyPrepFields& pr = suite_ ? suite_->prep : *prep_;
    dif::calculate_km_kh(p_, real(g.dx), real(g.dy), real(g.dt), real(b.cf1), real(b.cf2), real(b.cf3),
                         d11_.cview(), d22_.cview(), d33_.cview(), d12_.cview(), d13_.cview(), d23_.cview(),
                         s.P_TKE >= 0 ? s.q_2[s.P_TKE].cview() : View3<const real>{},
                         pr.p8w.cview(), pr.t8w.cview(), pr.th_phy.cview(), pr.t_phy.cview(), pr.p_phy.cview(),
                         s.P_QV >= 0 ? s.q_2[s.P_QV].cview() : View3<const real>{}, s.P_QV >= 0,
                         s.P_QC >= 0 ? s.q_2[s.P_QC].cview() : View3<const real>{}, s.P_QC >= 0,
                         s.P_QI >= 0 ? s.q_2[s.P_QI].cview() : View3<const real>{}, s.P_QI >= 0,
                         rdz_.cview(), rdzw_.cview(), zx_.cview(), zy_.cview(), b.msftx.view(), b.msfty.view(),
                         bn2_.view(), xkmh_.view(), xkmv_.view(), xkhh_.view(), xkhv_.view(), g.nx, g.ny, g.nz);
  }
  fill_halos({&xkmh_, &xkmv_, &xkhh_, &xkhv_}, g);
  if (c_.diff_opt == 2) {
    Field3& rho_f = suite_ ? suite_->prep.rho : prep_->rho;
    fill_halos({&div_, &d11_, &d22_, &d33_, &rho_f}, g);
    fill_halo(d13_, gp_, Parity::Even, Parity::Even);
    fill_halo(d23_, gp_, Parity::Even, Parity::Even);
  }
  vmdump::dump("vm_xk", 0, {&xkhh_, &xkmh_});
}

// first_rk_step_part2 (diff_opt 2) and rk_tendency, rk_step 1 (module_em.F 800-944)
void DiffusionDriver::dynamics(const State& s, const BaseState& b, Tendencies& t) const {
  const Grid& g = g_;
  const int nx = g.nx, ny = g.ny, nz = g.nz;
  const bool do_y = !g.two_d();
  const real rdx = real(b.rdx), rdy = real(b.rdy);
  if (c_.diff_opt == 1) {
    dif::horizontal_diffusion('u', p_, s.u_2.cview(), t.ru_tendf.view(), s.mut.cview(), b.msftx.view(),
                              b.msfty.view(), b.msfux.view(), b.msfuy.view(), b.msfvx.view(), b.msfvy.view(),
                              b.msfvx_inv.view(), b.c1h.view(), b.c2h.view(), xkmh_.cview(), rdx, rdy, nx, ny, nz);
    if (do_y)
      dif::horizontal_diffusion('v', p_, s.v_2.cview(), t.rv_tendf.view(), s.mut.cview(), b.msftx.view(),
                                b.msfty.view(), b.msfux.view(), b.msfuy.view(), b.msfvx.view(), b.msfvy.view(),
                                b.msfvx_inv.view(), b.c1h.view(), b.c2h.view(), xkmh_.cview(), rdx, rdy, nx, ny, nz);
    dif::horizontal_diffusion('w', p_, s.w_2.cview(), t.rw_tendf.view(), s.mut.cview(), b.msftx.view(),
                              b.msfty.view(), b.msfux.view(), b.msfuy.view(), b.msfvx.view(), b.msfvy.view(),
                              b.msfvx_inv.view(), b.c1f.view(), b.c2f.view(), xkmh_.cview(), rdx, rdy, nx, ny, nz);
    dif::horizontal_diffusion_3dmp(p_, s.t_2.cview(), b.t_init.view(), t.t_tendf.view(), s.mut.cview(),
                                   b.msftx.view(), b.msfty.view(), b.msfux.view(), b.msfuy.view(), b.msfvy.view(),
                                   b.msfvx_inv.view(), b.c1h.view(), b.c2h.view(), xkhh_.cview(), rdx, rdy, nx, ny,
                                   nz);
    if (!c_.pbl) {
      dif::vertical_diffusion_u(p_, s.u_2.cview(), t.ru_tendf.view(), b.u_base.view(), s.muu.cview(), b.c1h.view(),
                                b.c2h.view(), s.alt.cview(), b.rdn.view(), b.rdnw.view(), c_.kvdif, nx, ny, nz);
      if (do_y)
        dif::vertical_diffusion_v(p_, s.v_2.cview(), t.rv_tendf.view(), b.v_base.view(), s.muv.cview(),
                                  b.c1h.view(), b.c2h.view(), s.alt.cview(), b.rdn.view(), b.rdnw.view(), c_.kvdif,
                                  nx, ny, nz);
      dif::vertical_diffusion('w', s.w_2.cview(), t.rw_tendf.view(), s.mut.cview(), b.c1f.view(), b.c2f.view(),
                              s.alt.cview(), b.rdn.view(), b.rdnw.view(), c_.kvdif, nx, ny, nz);
      // vertical_diffusion_3dmp(t, t_init, kvdq = kvdif / prandtl)
      dif::vertical_diffusion_mp(s.t_2.cview(), t.t_tendf.view(), b.t_init.view(), s.mut.cview(), b.c1h.view(),
                                 b.c2h.view(), s.alt.cview(), b.rdn.view(), b.rdnw.view(), c_.kvdif / p_.prandtl, nx,
                                 ny, nz);
    }
  } else if (c_.diff_opt == 2) {
    const phys::driver::PhyPrepFields& pr = suite_ ? suite_->prep : *prep_;
    const Field2& ustm = suite_ ? suite_->ustm : ustm_;
    const Field2& hfx  = suite_ ? suite_->hfx  : hfx_;
    const Field2& qfx  = suite_ ? suite_->qfx  : qfx_;
    const bool has_qv  = s.P_QV >= 0;
    const bool has_tke = s.P_TKE >= 0;
    if (c_.km_opt == 2 && has_tke) {
      dif::tke_rhs(p_, real(g.dx), real(g.dy), real(g.dt), d11_.cview(), d22_.cview(), d33_.cview(),
                   d12_.cview(), d13_.cview(), d23_.cview(), s.u_2.cview(), s.v_2.cview(), s.q_2[s.P_TKE].cview(),
                   bn2_.cview(), pr.th_phy.cview(), pr.rho.cview(), rdzw_.cview(), xkmh_.cview(), xkmv_.cview(),
                   xkhv_.cview(), has_qv ? s.q_2[s.P_QV].cview() : View3<const real>{}, has_qv,
                   s.mut.cview(), b.c1h.view(), b.c2h.view(), b.msftx.view(), b.msfty.view(),
                   ustm.cview(), true, hfx.cview(), true, t.q_tendf[s.P_TKE].view(), nx, ny, nz);
    }
    const bool bdy_xs = p_.open_xs || p_.specified || p_.nested;
    const bool bdy_xe = p_.open_xe || p_.specified || p_.nested;
    const bool bdy_ys = p_.open_ys || p_.specified || p_.nested;
    const bool bdy_ye = p_.open_ye || p_.specified || p_.nested;
    const int is_s = (!p_.periodic_x && bdy_xs) ? 1 : 0;
    const int is_e = (!p_.periodic_x && bdy_xe) ? nx - 1 : nx;
    const int js_s = bdy_ys ? 1 : 0;
    const int js_e = bdy_ye ? ny - 1 : ny;
    if (!c_.pbl) {
      dif::cal_titau_vert(pr.rho.cview(), xkmh_.cview(), xkmv_.cview(), d13_.cview(), d23_.cview(), d33_.cview(),
                          b.fnm.view(), b.fnp.view(), titau13_.view(), titau23_.view(), titau33_.view(),
                          p_.specified, nx, ny, nz, p_);
      dif::vertical_diffusion_2(p_, titau13_.cview(), titau23_.cview(), titau33_.cview(), s.t_2.cview(),
                                b.t_base.view(), xkhv_.cview(), pr.rho.cview(), rdz_.cview(),
                                b.fnm.view(), b.fnp.view(), b.dn.view(), b.dnw.view(),
                                t.ru_tendf.view(), t.rv_tendf.view(), t.rw_tendf.view(), t.t_tendf.view(),
                                has_qv ? s.q_2[s.P_QV].cview() : View3<const real>{},
                                has_qv ? t.q_tendf[s.P_QV].view() : View3<real>{}, has_qv,
                                nx, ny, nz,
                                has_tke ? s.q_2[s.P_TKE].cview() : View3<const real>{},
                                has_tke ? t.q_tendf[s.P_TKE].view() : View3<real>{},
                                xkmv_.cview(), has_tke,
                                s.u_2.cview(), s.v_2.cview(), pr.th_phy.cview(),
                                b.qv_base_1d.view(), true, ustm.cview(), true, hfx.cview(), true, qfx.cview(), true);
      for (int n = 0; n < s.n_scalars(); ++n) {
        if (n == s.P_QV || n == s.P_TKE) continue;
        column_for(Range2D{is_s, is_e, js_s, js_e},
                   dif::KVertDiffS2{s.q_2[n].cview(), xkhv_.cview(), pr.rho.cview(), rdz_.cview(),
                                    b.fnm.view(), b.fnp.view(), b.dnw.view(), t.q_tendf[n].view(), nz, false});
      }
    }
    dif::horizontal_diffusion_2(p_, d11_.cview(), d22_.cview(), d33_.cview(), d12_.cview(), d13_.cview(),
                                d23_.cview(), pr.rho.cview(), xkmh_.cview(), xkmv_.cview(), xkhh_.cview(),
                                zx_.cview(), zy_.cview(), rdz_.cview(), rdzw_.cview(),
                                b.msftx.view(), b.msfty.view(), b.msfux.view(), b.msfuy.view(),
                                b.msfvx.view(), b.msfvy.view(), b.fnm.view(), b.fnp.view(),
                                b.dn.view(), b.dnw.view(), real(b.cf1), real(b.cf2), real(b.cf3), rdx, rdy,
                                titau11_.view(), titau22_.view(), titau33_.view(), titau12_.view(),
                                titau13_.view(), titau23_.view(),
                                t.ru_tendf.view(), t.rv_tendf.view(), t.rw_tendf.view(),
                                s.t_2.cview(), t.t_tendf.view(),
                                has_qv ? s.q_2[s.P_QV].cview() : View3<const real>{},
                                has_qv ? t.q_tendf[s.P_QV].view() : View3<real>{}, has_qv,
                                nx, ny, nz,
                                has_tke ? s.q_2[s.P_TKE].cview() : View3<const real>{},
                                has_tke ? t.q_tendf[s.P_TKE].view() : View3<real>{}, has_tke);
    for (int n = 0; n < s.n_scalars(); ++n) {
      if (n == s.P_QV || n == s.P_TKE) continue;
      const bool mix_off = n < s.n_moist ? c_.moist_mix2_off : c_.scalar_mix2_off;
      if (!mix_off) {
        column_for(Range2D{is_s, is_e, js_s, js_e},
                   dif::KHorizDiffS2{s.q_2[n].cview(), xkhh_.cview(), pr.rho.cview(), zx_.cview(), zy_.cview(),
                                     rdz_.cview(), rdzw_.cview(), b.msftx.view(), b.msfty.view(),
                                     b.msfux.view(), b.msfvy.view(), b.fnm.view(), b.fnp.view(),
                                     b.dn.view(), b.dnw.view(), real(b.cf1), real(b.cf2), real(b.cf3),
                                     t.q_tendf[n].view(), rdx, rdy, nz, false});
      }
    }
  }
  if (c_.diff_6th_opt != 0) {
    const real dt = real(g.dt);
    dif::sixth_order_diffusion('u', p_, s.u_2.cview(), t.ru_tendf.view(), s.mut.cview(), dt, b.c1h.view(),
                               b.c2h.view(), b.msftx.view(), b.msfty.view(), b.msfux.view(), b.msfuy.view(),
                               b.msfvx.view(), b.msfvy.view(), nx, ny, nz);
    if (do_y)
      dif::sixth_order_diffusion('v', p_, s.v_2.cview(), t.rv_tendf.view(), s.mut.cview(), dt, b.c1h.view(),
                                 b.c2h.view(), b.msftx.view(), b.msfty.view(), b.msfux.view(), b.msfuy.view(),
                                 b.msfvx.view(), b.msfvy.view(), nx, ny, nz);
    dif::sixth_order_diffusion('w', p_, s.w_2.cview(), t.rw_tendf.view(), s.mut.cview(), dt, b.c1f.view(),
                               b.c2f.view(), b.msftx.view(), b.msfty.view(), b.msfux.view(), b.msfuy.view(),
                               b.msfvx.view(), b.msfvy.view(), nx, ny, nz);
    dif::sixth_order_diffusion('m', p_, s.t_2.cview(), t.t_tendf.view(), s.mut.cview(), dt, b.c1h.view(),
                               b.c2h.view(), b.msftx.view(), b.msfty.view(), b.msfux.view(), b.msfuy.view(),
                               b.msfvx.view(), b.msfvy.view(), nx, ny, nz);
  }
  if (c_.damp_opt == 2) {
    dif::rk_rayleigh_damp(p_, s.u_2.cview(), s.v_2.cview(), s.w_2.cview(), s.t_2.cview(), s.ph_2.cview(),
                          b.phb.view(), s.mut.cview(), s.muu.cview(), s.muv.cview(),
                          b.c1h.view(), b.c2h.view(), b.c1f.view(), b.c2f.view(),
                          b.u_base.view(), b.v_base.view(), b.t_base.view(), b.z_base.view(),
                          t.ru_tendf.view(), t.rv_tendf.view(), t.rw_tendf.view(), t.t_tendf.view(),
                          nx, ny, nz);
  }
}

// rk_scalar_tend, rk_step 1 (module_em.F 1378-1440): mut there is muts
void DiffusionDriver::scalar(const State& s, const BaseState& b, int n, Field3& tend, real dt_rk) const {
  const Grid& g = g_;
  const int nx = g.nx, ny = g.ny, nz = g.nz;
  const Field3& q = s.q_2[n];
  const bool is_tke = (n == s.P_TKE);
  const bool moist = (n < s.n_moist);
  const bool mix2_off = is_tke ? c_.tke_mix2_off : (moist ? c_.moist_mix2_off : c_.scalar_mix2_off);
  const bool mix6_off = is_tke ? c_.tke_mix6_off : (moist ? c_.moist_mix6_off : c_.scalar_mix6_off);
  if (c_.diff_opt == 1 && !mix2_off) {
    dif::horizontal_diffusion('m', p_, q.cview(), tend.view(), s.muts.cview(), b.msftx.view(), b.msfty.view(),
                              b.msfux.view(), b.msfuy.view(), b.msfvx.view(), b.msfvy.view(), b.msfvx_inv.view(),
                              b.c1h.view(), b.c2h.view(), xkhh_.cview(), real(b.rdx), real(b.rdy), nx, ny, nz);
    if (!c_.pbl) {
      const real kvdq = c_.kvdif / p_.prandtl;  // WRF kvdq = kvdif / prandtl
      if (moist && n == s.P_QV)
        dif::vertical_diffusion_mp(q.cview(), tend.view(), b.qv_base.view(), s.muts.cview(), b.c1h.view(),
                                   b.c2h.view(), s.alt.cview(), b.rdn.view(), b.rdnw.view(), kvdq, nx, ny, nz);
      else
        dif::vertical_diffusion('m', q.cview(), tend.view(), s.muts.cview(), b.c1h.view(), b.c2h.view(),
                                s.alt.cview(), b.rdn.view(), b.rdnw.view(), kvdq, nx, ny, nz);
    }
  }
  if (c_.diff_6th_opt != 0 && !mix6_off)
    dif::sixth_order_diffusion('m', p_, q.cview(), tend.view(), s.muts.cview(), dt_rk, b.c1h.view(), b.c2h.view(),
                               b.msftx.view(), b.msfty.view(), b.msfux.view(), b.msfuy.view(), b.msfvx.view(),
                               b.msfvy.view(), nx, ny, nz);
}

void DiffusionDriver::bound_tke(Field3& tke) const {
  const real ub = p_.tke_upper_bound;
  View3<real> v = tke.view();
  parallel_for(g_.mass_cells(), [=] CPPWRF_HD(int i, int j, int k) {
    real x = v(i, j, k);
    if (x < real(0)) x = real(0);
    if (x > ub) x = ub;
    v(i, j, k) = x;
  });
}

std::vector<const Field3*> DiffusionDriver::fields() const {
  std::vector<const Field3*> v = {&d11_, &d22_, &d33_, &d12_, &d13_, &d23_, &div_, &xkmh_, &xkmv_, &xkhh_, &xkhv_};
  if (bn2_.layout().total > 0) v.push_back(&bn2_);
  if (prep_)
    for (const Field3* f : {&prep_->th_phy, &prep_->p_phy, &prep_->pi_phy, &prep_->t_phy, &prep_->rho, &prep_->u_phy,
                            &prep_->v_phy, &prep_->z, &prep_->p_hyd, &prep_->z_at_w, &prep_->dz8w, &prep_->p8w,
                            &prep_->t8w, &prep_->p_hyd_w})
      v.push_back(f);
  return v;
}

}  // namespace cppwrf::dyn
