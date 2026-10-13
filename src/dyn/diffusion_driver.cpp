#include "dyn/diffusion_driver.hpp"

#include <stdexcept>

#include "core/halo.hpp"
#include "dyn/diffusion/diffusion.hpp"

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
  if (c.diff_opt != 0 && c.diff_opt != 1) throw std::runtime_error("dyn.diff_opt: only 0 and 1 are wired");
  if (c.diff_opt == 1 && c.km_opt != 1 && c.km_opt != 4) throw std::runtime_error("dyn.km_opt must be 1 or 4");
  if (c.diff_6th_opt < 0 || c.diff_6th_opt > 2) throw std::runtime_error("dyn.diff_6th_opt must be 0, 1 or 2");
  if (g.nz + 1 >= 128) throw std::runtime_error("diffusion kernels support nz <= 126");
  p_.diff_opt = c.diff_opt;
  p_.km_opt = c.km_opt;
  p_.khdif = c.khdif;
  p_.kvdif = c.kvdif;
  p_.c_s = c.c_s;
  p_.diff_6th_opt = c.diff_6th_opt;
  p_.diff_6th_factor = c.diff_6th_factor;
  p_.periodic_x = g.bc_x == BC::Periodic;
  p_.periodic_y = g.bc_y == BC::Periodic;
  p_.specified = g.specified();
}

// first_rk_step_part2 (rk_step 1): metrics, deformation, eddy viscosities
void DiffusionDriver::prepare(const State& s, const BaseState& b) {
  if (c_.diff_opt != 1) return;
  const Grid& g = g_;
  const real rdx = real(b.rdx), rdy = real(b.rdy);
  dif::compute_diff_metrics(p_, s.ph_2.cview(), b.phb.view(), z_.view(), rdz_.view(), rdzw_.view(), zx_.view(),
                            zy_.view(), rdx, rdy, g.nx, g.ny, g.nz);
  fill_halos({&rdz_, &rdzw_}, g);
  fill_halo(zx_, gp_, Parity::Even, Parity::Even);
  fill_halo(zy_, gp_, Parity::Even, Parity::Even);
  dif::cal_deform_and_div(p_, s.u_2.cview(), s.v_2.cview(), s.w_2.cview(), b.u_base.view(), b.v_base.view(),
                          zx_.cview(), zy_.cview(), rdz_.cview(), rdzw_.cview(), b.msftx.view(), b.msfty.view(),
                          b.msfux.view(), b.msfuy.view(), b.msfvx.view(), b.msfvy.view(), b.fnm.view(), b.fnp.view(),
                          b.dn.view(), b.dnw.view(), real(b.cf1), real(b.cf2), real(b.cf3), rdx, rdy, d11_.view(),
                          d22_.view(), d33_.view(), d12_.view(), d13_.view(), d23_.view(), div_.view(), g.nx, g.ny,
                          g.nz);
  dif::calculate_km_kh(p_, real(g.dx), real(g.dy), d11_.cview(), d22_.cview(), d12_.cview(), rdzw_.cview(),
                       zx_.cview(), zy_.cview(), b.msftx.view(), b.msfty.view(), xkmh_.view(), xkmv_.view(),
                       xkhh_.view(), xkhv_.view(), g.nx, g.ny, g.nz);
  fill_halos({&xkmh_, &xkmv_, &xkhh_, &xkhv_}, g);
}

// rk_tendency, rk_step 1 (module_em.F 800-935)
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
      // vertical_diffusion_3dmp(t, t_init, kvdq = 3 kvdif)
      dif::vertical_diffusion_mp(s.t_2.cview(), t.t_tendf.view(), b.t_init.view(), s.mut.cview(), b.c1h.view(),
                                 b.c2h.view(), s.alt.cview(), b.rdn.view(), b.rdnw.view(), real(3) * c_.kvdif, nx,
                                 ny, nz);
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
}

// rk_scalar_tend, rk_step 1 (module_em.F 1378-1440): mut there is muts
void DiffusionDriver::scalar(const State& s, const BaseState& b, int n, Field3& tend, real dt_rk) const {
  const Grid& g = g_;
  const int nx = g.nx, ny = g.ny, nz = g.nz;
  const Field3& q = s.q_2[n];
  if (c_.diff_opt == 1) {
    dif::horizontal_diffusion('m', p_, q.cview(), tend.view(), s.muts.cview(), b.msftx.view(), b.msfty.view(),
                              b.msfux.view(), b.msfuy.view(), b.msfvx.view(), b.msfvy.view(), b.msfvx_inv.view(),
                              b.c1h.view(), b.c2h.view(), xkhh_.cview(), real(b.rdx), real(b.rdy), nx, ny, nz);
    if (!c_.pbl) {
      const real kvdq = real(3) * c_.kvdif;
      if (n < s.n_moist && n == s.P_QV)
        dif::vertical_diffusion_mp(q.cview(), tend.view(), b.qv_base.view(), s.muts.cview(), b.c1h.view(),
                                   b.c2h.view(), s.alt.cview(), b.rdn.view(), b.rdnw.view(), kvdq, nx, ny, nz);
      else
        dif::vertical_diffusion('m', q.cview(), tend.view(), s.muts.cview(), b.c1h.view(), b.c2h.view(),
                                s.alt.cview(), b.rdn.view(), b.rdnw.view(), kvdq, nx, ny, nz);
    }
  }
  if (c_.diff_6th_opt != 0)
    dif::sixth_order_diffusion('m', p_, q.cview(), tend.view(), s.muts.cview(), dt_rk, b.c1h.view(), b.c2h.view(),
                               b.msftx.view(), b.msfty.view(), b.msfux.view(), b.msfuy.view(), b.msfvx.view(),
                               b.msfvy.view(), nx, ny, nz);
}

}  // namespace cppwrf::dyn
