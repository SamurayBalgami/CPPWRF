// PhysicsSuite::cumulus_driver, convtrans_prep and advance_ppt (WRF phys/module_cumulus_driver.F,
// dyn_em/module_convtrans_prep.F, phys/module_physics_addtendc.F)
#include "phys/driver/suite.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

#include "phys/cu_gfl/gfl_driver.hpp"
#include "phys/cu_ntiedtke/ntiedtke.hpp"
#include "phys/cu_tiedtke/tiedtke.hpp"
#include "phys/driver/wrf_constants.hpp"

#include "phys/driver/suite_internal.hpp"

namespace cppwrf::phys::driver {

using namespace suite_detail;

namespace {

// ------------------------------------------------------------ Grell-Freitas-Li
// module_cu_gfl_wrf.F90 cu_grell_freitas_li (GFLSCHEME, cu_physics = 3).
struct KGrellFreitas {
  int nz, ishallow, ichoice;
  real dt, dx, xlv, cp, g, r_v;
  C3 u, v, w, t, pi, rho, dz8w, pcps, p8w, qv;
  C2 xland, hfx, qfx;
  CI2 kpbl;
  BView2<const real> ht;
  C3 rthraten, rthblten, rqvblten, rthften, rqvften;
  V2 raincv, pratec, htop, hbot;
  V3 rthcuten, rqvcuten, rqccuten, rqicuten, rucuten, rvcuten, gd_cloud, gd_cloud2;
  bool flag_qc, flag_qi;
  CPPWRF_HD void column(int i, int j, cu_gfl::GflDrvColumn<real>& c) const {
    c.nz = nz; c.dt = dt; c.dx = dx; c.ht = ht(i, j); c.hfx = hfx(i, j); c.qfx = qfx(i, j); c.xland = xland(i, j);
    c.xlv = xlv; c.cp = cp; c.g = g; c.r_v = r_v;
    c.kpbl = kpbl(i, j);
    for (int k = 0; k < cu_gfl::KM; ++k) {
      c.u[k] = 0; c.v[k] = 0; c.t[k] = 0; c.q[k] = 0; c.p[k] = 0; c.pi[k] = 0; c.rho[k] = 0; c.dz8w[k] = 0;
      c.w[k] = 0; c.p8w[k] = 0; c.rthften[k] = 0; c.rthraten[k] = 0; c.rthblten[k] = 0; c.rqvften[k] = 0;
      c.rqvblten[k] = 0;
    }
    for (int k = 0; k < nz; ++k) {
      c.u[k] = u(i, j, k); c.v[k] = v(i, j, k); c.t[k] = t(i, j, k); c.q[k] = qv(i, j, k); c.p[k] = pcps(i, j, k);
      c.pi[k] = pi(i, j, k); c.rho[k] = rho(i, j, k); c.dz8w[k] = dz8w(i, j, k);
      c.rthften[k] = rthften(i, j, k); c.rthraten[k] = rthraten(i, j, k); c.rthblten[k] = rthblten(i, j, k);
      c.rqvften[k] = rqvften(i, j, k); c.rqvblten[k] = rqvblten(i, j, k);
    }
    for (int k = 0; k <= nz; ++k) {
      c.w[k] = w(i, j, k);
      c.p8w[k] = p8w(i, j, k);
    }
  }
  CPPWRF_HD void operator()(int i, int j) const {
    cu_gfl::GflDrvColumn<real> c;
    column(i, j, c);
    cu_gfl::GflDrvOut<real> o;
    cu_gfl::gfl_column(c, ishallow, ichoice, o);
    raincv(i, j) = o.raincv;
    pratec(i, j) = o.pratec;
    htop(i, j) = o.htop;
    hbot(i, j) = o.hbot;
    for (int k = 0; k < nz; ++k) {
      rthcuten(i, j, k) = o.rthcuten[k];
      rqvcuten(i, j, k) = o.rqvcuten[k];
      rqccuten(i, j, k) = o.rqccuten[k];
      rqicuten(i, j, k) = o.rqicuten[k];
      rucuten(i, j, k) = o.rucuten[k];
      rvcuten(i, j, k) = o.rvcuten[k];
      gd_cloud(i, j, k) = o.gdc[k];
      gd_cloud2(i, j, k) = o.gdc2[k];
      (void)flag_qc;
      (void)flag_qi;
    }
  }
};

// ------------------------------------------------------------ convtrans_prep
// dyn_em/module_convtrans_prep.F convtrans_prep (cu_diag == 1).
struct KConvtransPrep {
  int nz, ktau, stepave_count, cu_phys;
  bool avg_end_flag, first_period_flag, flag_qc, flag_qi;
  C2 raincv;
  V2 raincv_a, raincv_b;
  C3 gd_cloud, gd_cloud2, t_phy, p_phy, qv, qc, qi;
  V3 gd_cloud_a, gd_cloud2_a, qc_cu, qi_cu, gd_cldfr;
  CPPWRF_HD void operator()(int i, int j) const {
    if (ktau <= 1) {
      raincv_a(i, j) = real(0.);
      raincv_b(i, j) = real(0.);
      for (int k = 0; k < nz; ++k) {
        gd_cloud2_a(i, j, k) = real(0.);
        gd_cldfr(i, j, k) = real(0.);
      }
    }
    if (first_period_flag) {
      raincv_b(i, j) = raincv(i, j);
    }
    raincv_a(i, j) = raincv_a(i, j) + raincv(i, j);
    if (avg_end_flag) {
      raincv_b(i, j) = raincv_a(i, j) / real(stepave_count);
      raincv_a(i, j) = real(0.);
    }
    if (cu_phys == 3 || cu_phys == 5 || cu_phys == 93) {
      if (first_period_flag) {
        for (int k = 0; k < nz; ++k) {
          qc_cu(i, j, k) = gd_cloud(i, j, k);
          qi_cu(i, j, k) = gd_cloud2(i, j, k);
        }
      }
      for (int k = 0; k < nz; ++k) {
        gd_cloud_a(i, j, k) = gd_cloud_a(i, j, k) + gd_cloud(i, j, k);
        gd_cloud2_a(i, j, k) = gd_cloud2_a(i, j, k) + gd_cloud2(i, j, k);
      }
      if (avg_end_flag) {
        for (int k = 0; k < nz; ++k) {
          qc_cu(i, j, k) = real(.1) * gd_cloud_a(i, j, k) / real(stepave_count);
          gd_cloud_a(i, j, k) = real(0.);
          qi_cu(i, j, k) = real(.1) * gd_cloud2_a(i, j, k) / real(stepave_count);
          gd_cloud2_a(i, j, k) = real(0.);
        }
      }
    }
    if (first_period_flag || avg_end_flag) {
      for (int k = 0; k < nz; ++k) {
        gd_cldfr(i, j, k) = real(0.);
        if (flag_qc && !flag_qi) {
          real satvp = real(3.80) * std::exp(real(17.27) * (t_phy(i, j, k) - real(273.)) /
                                             (t_phy(i, j, k) - real(36.))) /
                       (real(.01) * p_phy(i, j, k));
          real rhgrid = std::max(real(.1), std::min(real(.95), qv(i, j, k) / satvp));
          real h2oliq = real(1000.) * (qc_cu(i, j, k) + qc(i, j, k));
          satvp = real(1000.) * satvp;
          real c = (real(1.) - std::exp((-real(100.) * h2oliq) /
                                        rpow((real(1.) - rhgrid) * satvp, real(0.49)))) *
                   rpow(rhgrid, real(0.25));
          c = std::max(real(0.), std::min(real(1.), c));
          if (qc(i, j, k) == real(0.)) c = c * real(.2);
          gd_cldfr(i, j, k) = c;
        }
        if (flag_qc && flag_qi) {
          real satvp = real(3.80) * std::exp(real(17.27) * (t_phy(i, j, k) - real(273.)) /
                                             (t_phy(i, j, k) - real(36.))) /
                       (real(.01) * p_phy(i, j, k));
          real rhgrid = std::max(real(.1), std::min(real(.95), qv(i, j, k) / satvp));
          real h2oliq = real(1000.) * (qc_cu(i, j, k) + qc(i, j, k) + qi_cu(i, j, k) + qi(i, j, k));
          satvp = real(1000.) * satvp;
          real c = (real(1.) - std::exp((-real(100.) * h2oliq) /
                                        rpow((real(1.) - rhgrid) * satvp, real(0.49)))) *
                   rpow(rhgrid, real(0.25));
          c = std::max(real(0.), std::min(real(1.), c));
          if (qc(i, j, k) == real(0.) && qi(i, j, k) == real(0.)) c = c * real(.2);
          gd_cldfr(i, j, k) = c;
        }
      }
    }
  }
};

// ---------------------------------------------------------------- New Tiedtke
// module_cumulus_driver.F (G3/NTIEDTKE forcing: RTHFTEN = (RTHFTEN + RTHRATEN + RTHBLTEN)*pi,
// RQVFTEN += RQVBLTEN) and module_cu_ntiedtke.F cu_ntiedtke_driver (pre_run, cu_ntiedtke_run,
// post_run) for one column. C++ k = WRF k - 1; the scheme's arrays run top (1) to bottom (nz).
struct KNTiedtke {
  cu_ntiedtke::NtConst<real> c;
  int nz; long itimestep; int stepcu; real dt;
  C3 u, v, w, t, pi, rho, dz8w, pcps, p8w, qv, qc, qi;
  C2 xland, dx2d, hfx, qfx;
  C3 rthraten, rthblten, rqvblten;
  V3 rthften, rqvften;
  V2 raincv, pratec;
  V3 rthcuten, rqvcuten, rqccuten, rqicuten, rucuten, rvcuten;
  bool flag_qc, flag_qi;
  void* dump;  // debug (CPU, one thread): FILE* receiving inputs and outputs of cu_ntiedtke_run
  CPPWRF_HD void operator()(int i, int j) const {
    using cu_ntiedtke::KM;
    using namespace wrfc;
    for (int k = 0; k < nz; ++k) {
      rthften(i, j, k) = (rthften(i, j, k) + rthraten(i, j, k) + rthblten(i, j, k)) * pi(i, j, k);
      rqvften(i, j, k) = rqvften(i, j, k) + rqvblten(i, j, k);
    }
    // cu_ntiedtke_pre_run
    const real delt = dt * real(stepcu);
    const int slimsk = int(std::fabs(xland(i, j) - real(2.)));
    real zi[KM], zl[KM], dot[KM];
    zi[0] = real(0.);
    for (int k = 0; k < nz; ++k) zi[k + 1] = zi[k] + dz8w(i, j, k);
    for (int k = 0; k < nz; ++k) {
      zl[k] = real(0.5) * (zi[k] + zi[k + 1]);
      dot[k] = -real(0.5) * g * rho(i, j, k) * (w(i, j, k) + w(i, j, k + 1));
    }
    real ghti[KM], prsi[KM], ghtl[KM], omg[KM], prsl[KM];
    real tf[KM], qvf[KM], qcf[KM], qif[KM], uf[KM], vf[KM], qvftenz[KM], thftenz[KM];
    for (int k = 0; k <= nz; ++k) {
      ghti[nz + 1 - k] = zi[k];
      prsi[nz + 1 - k] = p8w(i, j, k);
    }
    for (int k = 0; k < nz; ++k) {
      const int zz = nz - k;
      ghtl[zz] = zl[k]; omg[zz] = dot[k]; prsl[zz] = pcps(i, j, k);
      tf[zz] = t(i, j, k); qvf[zz] = qv(i, j, k); qcf[zz] = qc(i, j, k); qif[zz] = qi(i, j, k);
      uf[zz] = u(i, j, k); vf[zz] = v(i, j, k);
      qvftenz[zz] = itimestep == 1 ? real(0.) : rqvften(i, j, k);
      thftenz[zz] = itimestep == 1 ? real(0.) : rthften(i, j, k);
    }
    real rn = real(0.);
#if !defined(__CUDA_ARCH__)
    auto wr = [&](const real* a, int n) { std::fwrite(a + 1, sizeof(real), n, static_cast<std::FILE*>(dump)); };
    if (dump) {
      const double sc[6] = {double(qfx(i, j)), double(hfx(i, j)), double(slimsk), double(delt), double(dx2d(i, j)),
                            double(nz)};
      std::fwrite(sc, sizeof(double), 6, static_cast<std::FILE*>(dump));
      for (const real* a : {uf, vf, tf, qvf, qcf, qif, qvftenz, thftenz, ghtl, omg, prsl}) wr(a, nz);
      wr(ghti, nz + 1); wr(prsi, nz + 1);
    }
#endif
    cu_ntiedtke::cu_ntiedtke_run<real>(c, nz, uf, vf, tf, qvf, qcf, qif, qvftenz, thftenz, ghtl, ghti, omg, prsl,
                                       prsi, qfx(i, j), hfx(i, j), rn, slimsk, delt, dx2d(i, j));
#if !defined(__CUDA_ARCH__)
    if (dump) {
      for (const real* a : {uf, vf, tf, qvf, qcf, qif}) wr(a, nz);
      std::fwrite(&rn, sizeof(real), 1, static_cast<std::FILE*>(dump));
    }
#endif
    // cu_ntiedtke_post_run
    const real rdelt = real(1.) / delt;
    raincv(i, j) = rn / real(stepcu);
    pratec(i, j) = rn / (real(stepcu) * dt);
    for (int k = 0; k < nz; ++k) {
      const int zz = nz - k;
      rthcuten(i, j, k) = (tf[zz] - t(i, j, k)) / pi(i, j, k) * rdelt;
      rqvcuten(i, j, k) = (qvf[zz] - qv(i, j, k)) * rdelt;
      if (flag_qc) rqccuten(i, j, k) = (qcf[zz] - qc(i, j, k)) * rdelt;
      if (flag_qi) rqicuten(i, j, k) = (qif[zz] - qi(i, j, k)) * rdelt;
      rucuten(i, j, k) = (uf[zz] - u(i, j, k)) * rdelt;
      rvcuten(i, j, k) = (vf[zz] - v(i, j, k)) * rdelt;
    }
  }
};

// ---------------------------------------------------------------- Tiedtke
// module_cu_tiedtke.F cu_tiedtke (TIEDTKESCHEME, cu_physics = 6): the column inputs of tiecnv. Pass 1
// (first = true) stores the column's icum after the first cloud ascent in icum1; pass 2 runs tiecnv with
// the row's icum (row_icum) and writes the tendencies.
struct KTiedtke {
  cu_tiedtke::TdConst<real> c;
  int nz, nx; long itimestep; int stepcu; real dt;
  bool first;
  C3 u, v, w, t, pi, rho, dz8w, pcps, p8w, qv, qc, qi, rqvften, rqvblten;
  View1<const real> znu;
  C2 xland, qfx;
  I2 icum1;
  CI2 row_icum;
  V2 raincv, pratec;
  V3 rthcuten, rqvcuten, rqccuten, rqicuten, rucuten, rvcuten;
  bool flag_qc, flag_qi;
  CPPWRF_HD void operator()(int i, int j) const {
    using cu_tiedtke::KM;
    using namespace wrfc;
    const real delt = dt * real(stepcu);
    const real rdelt = real(1.) / delt;
    real zi[KM], zl[KM];
    zi[0] = real(0.0);  // zi(kts) = 0, zi(k) = zi(k-1) + dz8w(k-1), k = kts+1 .. kte
    for (int k = 1; k < nz; ++k) zi[k] = zi[k - 1] + dz8w(i, j, k - 1);
    for (int k = 1; k < nz; ++k) zl[k - 1] = (zi[k] + zi[k - 1]) * real(0.5);
    zl[nz - 1] = real(2.) * zi[nz - 1] - zl[nz - 2];
    const int slimsk = int(std::fabs(xland(i, j) - real(2.)));
    cu_tiedtke::TdInput<real> in;
    for (int k = 0; k < nz; ++k) {
      const int zz = nz - k;
      in.pu[zz] = u(i, j, k);
      in.pv[zz] = v(i, j, k);
      in.pt[zz] = t(i, j, k);
      in.pqv[zz] = qv(i, j, k);
      in.pqvf[zz] = itimestep == 1 ? real(0.) : rqvften(i, j, k);
      in.pqvbl[zz] = itimestep == 1 ? real(0.) : rqvblten(i, j, k);
      in.pqc[zz] = qc(i, j, k);
      in.pqi[zz] = qi(i, j, k);
      in.pomg[zz] = -real(0.5) * g * rho(i, j, k) * (w(i, j, k) + w(i, j, k + 1));
      in.poz[zz] = zl[k];
      in.pap[zz] = pcps(i, j, k);
      in.sig1[zz] = znu(k);
    }
    for (int k = 0; k <= nz; ++k) in.paph[nz + 1 - k] = p8w(i, j, k);
    const real evap = qfx(i, j);
    if (first) {
      icum1(i, j) = cu_tiedtke::tiedtke_first_ascent<real>(c, nz, in, evap, slimsk, delt, i, nx);
      return;
    }
    real u1[KM], v1[KM], t1[KM], q1[KM], q2[KM], q3[KM];
    for (int k = 1; k <= nz; ++k) {
      u1[k] = in.pu[k]; v1[k] = in.pv[k]; t1[k] = in.pt[k]; q1[k] = in.pqv[k]; q2[k] = in.pqc[k]; q3[k] = in.pqi[k];
    }
    const real rn = cu_tiedtke::tiecnv<real>(c, nz, in, u1, v1, t1, q1, q2, q3, evap, slimsk, delt, i, nx, row_icum(0, j));
    raincv(i, j) = rn / real(stepcu);
    pratec(i, j) = rn / (real(stepcu) * dt);
    for (int k = 0; k < nz; ++k) {
      const int zz = nz - k;
      rthcuten(i, j, k) = (t1[zz] - t(i, j, k)) / pi(i, j, k) * rdelt;
      rqvcuten(i, j, k) = (q1[zz] - qv(i, j, k)) * rdelt;
      rucuten(i, j, k) = (u1[zz] - u(i, j, k)) * rdelt;
      rvcuten(i, j, k) = (v1[zz] - v(i, j, k)) * rdelt;
      if (flag_qc) rqccuten(i, j, k) = (q2[zz] - qc(i, j, k)) * rdelt;
      if (flag_qi) rqicuten(i, j, k) = (q3[zz] - qi(i, j, k)) * rdelt;
    }
  }
};

// row_icum(0, j) = OR of icum1(i, j) over the row (cumastr_new's icum for the row)
struct KRowAny {
  CI2 icum1;
  I2 row_icum;
  int nx;
  CPPWRF_HD void operator()(int, int j) const {
    int any = 0;
    for (int i = 0; i < nx; ++i) any = any || icum1(i, j) != 0;
    row_icum(0, j) = any;
  }
};

// ---------------------------------------------------------------- advance_ppt
// phys/module_physics_addtendc.F advance_ppt (CASE KFETASCHEME)
struct KAdvancePpt {
  V3 rthcuten, rqvcuten, rqccuten, rqrcuten, rqicuten, rqscuten;
  V2 rainc, rainsh, pratec, pratesh, nca, htop, hbot, cuppt;
  C2 cutop, cubot;
  real dt;
  int nz, kde;
  bool flag_qi, flag_qs, kfeta;  // kfeta: CASE (KFETASCHEME); otherwise CASE DEFAULT (accumulation only)
  CPPWRF_HD static int nint_f(real x) { return int(std::round(x)); }  // Fortran NINT (lround: half away from 0)
  CPPWRF_HD void operator()(int i, int j) const {
    rainc(i, j) = rainc(i, j) + pratec(i, j) * dt;
    rainsh(i, j) = rainsh(i, j) + pratesh(i, j) * dt;
    cuppt(i, j) = cuppt(i, j) + (pratec(i, j) + pratesh(i, j)) * dt / real(1000.);
    if (!kfeta) return;
    const int ncutop = nint_f(cutop(i, j));
    const int ncubot = nint_f(cubot(i, j));
    if (ncutop > 1 && ncutop < kde) {
      if (cutop(i, j) > htop(i, j)) htop(i, j) = cutop(i, j);
    }
    if (ncubot > 0 && ncubot < kde) {
      if (cubot(i, j) < hbot(i, j)) hbot(i, j) = cubot(i, j);
    }
    if (nca(i, j) > real(0.)) {
      if (nint_f(nca(i, j) / dt) <= 1) {
        for (int k = 0; k < nz; ++k) {
          rthcuten(i, j, k) = real(0.);
          rqvcuten(i, j, k) = real(0.);
          rqccuten(i, j, k) = real(0.);
          rqrcuten(i, j, k) = real(0.);
          if (flag_qi) rqicuten(i, j, k) = real(0.);
          if (flag_qs) rqscuten(i, j, k) = real(0.);
        }
      }
      nca(i, j) = nca(i, j) - dt;
    }
  }
};

}  // namespace

// module_cumulus_driver.F cumulus_driver for KFETASCHEME (cu_physics = 1)
void PhysicsSuite::cumulus_driver(const dyn::State& s, const dyn::BaseState& b) {
  if (o_.cu_physics == 0) return;
  const bool run_param = itimestep == 1 || o_.cudt == 0 || stepcu == 1 || itimestep % stepcu == 0;
  if (!run_param) return;
  using namespace wrfc;
  if (o_.cu_physics == 3) {
    KGrellFreitas k{};
    k.nz = g_.nz; k.ishallow = o_.ishallow; k.ichoice = o_.gf_ichoice;
    k.dt = real(g_.dt); k.dx = real(g_.dx);
    k.xlv = xlv; k.cp = cp; k.g = g; k.r_v = r_v;
    k.u = prep.u_phy.cview(); k.v = prep.v_phy.cview(); k.w = s.w_2.cview(); k.t = prep.t_phy.cview();
    k.pi = prep.pi_phy.cview(); k.rho = prep.rho.cview(); k.dz8w = prep.dz8w.cview(); k.pcps = prep.p_hyd.cview();
    k.p8w = prep.p_hyd_w.cview(); k.qv = s.q_2[s.P_QV].cview();
    k.xland = xland.cview(); k.hfx = hfx.cview(); k.qfx = qfx.cview(); k.kpbl = kpbl.cview(); k.ht = b.ht.view();
    k.rthraten = pt.rthraten.cview(); k.rthblten = pt.rthblten.cview(); k.rqvblten = pt.rqvblten.cview();
    k.rthften = pt.rthften.cview(); k.rqvften = pt.rqvften.cview();
    // first_rk_step_part1 passes HTOP=grid%cutop, HBOT=grid%cubot to cumulus_driver; advance_ppt
    // hits CASE DEFAULT for GFLSCHEME, so grid%htop and grid%hbot remain 0.
    k.raincv = raincv.view(); k.pratec = pratec.view(); k.htop = cutop.view(); k.hbot = cubot.view();
    k.rthcuten = pt.rthcuten.view(); k.rqvcuten = pt.rqvcuten.view(); k.rqccuten = pt.rqccuten.view();
    k.rqicuten = pt.rqicuten.view(); k.rucuten = pt.rucuten.view(); k.rvcuten = pt.rvcuten.view();
    k.gd_cloud = gd_cloud.view(); k.gd_cloud2 = gd_cloud2.view();
    k.flag_qc = s.P_QC >= 0; k.flag_qi = s.P_QI >= 0;
    column_for(g_.columns(), k);
    // debug: CPPWRF_GFL_DUMP=file CPPWRF_GFL_DUMP_STEP=n (or CPPWRF_GF_DUMP / CPPWRF_GF_DUMP_STEP)
    // writes cu_grell_freitas_li's inputs and outputs of step n
    // (tests/phys/cu_gfl/harness.f90 runs WRF's cu_grell_freitas_li on them)
    const char* d = std::getenv("CPPWRF_GFL_DUMP");
    if (!d) d = std::getenv("CPPWRF_GF_DUMP");
    const char* ds = std::getenv("CPPWRF_GFL_DUMP_STEP");
    if (!ds) ds = std::getenv("CPPWRF_GF_DUMP_STEP");
    if (d && ds && itimestep == std::atol(ds)) {
      std::FILE* f = std::fopen(d, "wb");
      const int hdr[3] = {g_.nx, g_.ny, g_.nz};
      std::fwrite(hdr, sizeof(int), 3, f);
      auto w3 = [&](const Field3& x, int nk) {
        std::vector<real> h;
        x.download(h);
        const Layout& L = x.layout();
        for (int j = 0; j < g_.ny; ++j)
          for (int kk = 0; kk < nk; ++kk)
            for (int i = 0; i < g_.nx; ++i) {
              const double v = double(h[L.index(i, j, kk)]);
              std::fwrite(&v, sizeof(double), 1, f);
            }
      };
      auto w2 = [&](const Field2& x) {
        std::vector<real> h;
        x.download(h);
        const Layout& L = x.layout();
        for (int j = 0; j < g_.ny; ++j)
          for (int i = 0; i < g_.nx; ++i) {
            const double v = double(h[L.index(i, j, 0)]);
            std::fwrite(&v, sizeof(double), 1, f);
          }
      };
      const int nz = g_.nz;
      w3(prep.u_phy, nz); w3(prep.v_phy, nz); w3(prep.t_phy, nz); w3(s.q_2[s.P_QV], nz); w3(prep.p_hyd, nz);
      w3(prep.pi_phy, nz); w3(prep.rho, nz); w3(prep.dz8w, nz); w3(s.w_2, nz + 1); w3(prep.p_hyd_w, nz + 1);
      w3(pt.rthften, nz); w3(pt.rthraten, nz); w3(pt.rthblten, nz); w3(pt.rqvften, nz); w3(pt.rqvblten, nz);
      {
        std::vector<real> h;
        b.ht.download(h);
        const Layout& L = b.ht.layout();
        for (int j = 0; j < g_.ny; ++j)
          for (int i = 0; i < g_.nx; ++i) {
            const double v = b.ht.uniform() ? double(h[L.index(0, 0, 0)]) : double(h[L.index(i, j, 0)]);
            std::fwrite(&v, sizeof(double), 1, f);
          }
      }
      w2(hfx); w2(qfx); w2(xland);
      {
        std::vector<int> h;
        kpbl.download(h);
        const Layout& L = kpbl.layout();
        for (int j = 0; j < g_.ny; ++j)
          for (int i = 0; i < g_.nx; ++i) {
            const double v = double(h[L.index(i, j, 0)]);
            std::fwrite(&v, sizeof(double), 1, f);
          }
      }
      w3(pt.rthcuten, nz); w3(pt.rqvcuten, nz); w3(pt.rqccuten, nz); w3(pt.rqicuten, nz);
      w3(pt.rucuten, nz); w3(pt.rvcuten, nz);
      w2(raincv); w2(pratec); w2(cutop); w2(cubot);
      std::fclose(f);
    }
    return;
  }
  if (o_.cu_physics == 6) {
    KTiedtke k{};
    k.c = cu_tiedtke::make_tdconst<real>(r_d, r_v, cp, xlv, xls, xlf, g);
    k.nz = g_.nz; k.nx = g_.nx; k.itimestep = itimestep; k.stepcu = stepcu; k.dt = real(g_.dt);
    k.u = prep.u_phy.cview(); k.v = prep.v_phy.cview(); k.w = s.w_2.cview(); k.t = prep.t_phy.cview();
    k.pi = prep.pi_phy.cview(); k.rho = prep.rho.cview(); k.dz8w = prep.dz8w.cview(); k.pcps = prep.p_hyd.cview();
    k.p8w = prep.p_hyd_w.cview();
    k.qv = s.q_2[s.P_QV].cview();
    k.qc = s.P_QC >= 0 ? s.q_2[s.P_QC].cview() : zero3_.cview();
    k.qi = s.P_QI >= 0 ? s.q_2[s.P_QI].cview() : zero3_.cview();
    k.rqvften = pt.rqvften.cview(); k.rqvblten = pt.rqvblten.cview();
    k.znu = b.znu.view();
    k.xland = xland.cview(); k.qfx = qfx.cview();
    k.icum1 = td_icum.view(); k.row_icum = td_row_icum.cview();
    k.raincv = raincv.view(); k.pratec = pratec.view();
    k.rthcuten = pt.rthcuten.view(); k.rqvcuten = pt.rqvcuten.view(); k.rqccuten = pt.rqccuten.view();
    k.rqicuten = pt.rqicuten.view(); k.rucuten = pt.rucuten.view(); k.rvcuten = pt.rvcuten.view();
    k.flag_qc = s.P_QC >= 0; k.flag_qi = s.P_QI >= 0;
    k.first = true;
    column_for(g_.columns(), k);
    column_for(Range2D{0, 1, 0, g_.ny}, KRowAny{td_icum.cview(), td_row_icum.view(), g_.nx});
    k.first = false;
    column_for(g_.columns(), k);
    return;
  }
  if (o_.cu_physics == 16) {
    KNTiedtke k{};
    k.c = cu_ntiedtke::make_ntconst<real>(cp, r_d, r_v, xlv, xls, xlf, g);
    k.nz = g_.nz; k.itimestep = itimestep; k.stepcu = stepcu; k.dt = real(g_.dt);
    k.u = prep.u_phy.cview(); k.v = prep.v_phy.cview(); k.w = s.w_2.cview(); k.t = prep.t_phy.cview();
    k.pi = prep.pi_phy.cview(); k.rho = prep.rho.cview(); k.dz8w = prep.dz8w.cview(); k.pcps = prep.p_hyd.cview();
    k.p8w = prep.p_hyd_w.cview();
    k.qv = s.q_2[s.P_QV].cview();
    k.qc = s.P_QC >= 0 ? s.q_2[s.P_QC].cview() : zero3_.cview();
    k.qi = s.P_QI >= 0 ? s.q_2[s.P_QI].cview() : zero3_.cview();
    k.xland = xland.cview(); k.dx2d = dx2d.cview(); k.hfx = hfx.cview(); k.qfx = qfx.cview();
    k.rthraten = pt.rthraten.cview(); k.rthblten = pt.rthblten.cview(); k.rqvblten = pt.rqvblten.cview();
    k.rthften = pt.rthften.view(); k.rqvften = pt.rqvften.view();
    k.raincv = raincv.view(); k.pratec = pratec.view();
    k.rthcuten = pt.rthcuten.view(); k.rqvcuten = pt.rqvcuten.view(); k.rqccuten = pt.rqccuten.view();
    k.rqicuten = pt.rqicuten.view(); k.rucuten = pt.rucuten.view(); k.rvcuten = pt.rvcuten.view();
    k.flag_qc = s.P_QC >= 0; k.flag_qi = s.P_QI >= 0;
    k.dump = nullptr;
#if !defined(CPPWRF_USE_CUDA)
    // debug: CPPWRF_NTIEDTKE_DUMP=file CPPWRF_NTIEDTKE_DUMP_STEP=n (OMP_NUM_THREADS=1) writes the scheme's
    // column inputs and outputs of step n (tests/phys/cu_ntiedtke/harness.f90 replays them through WRF)
    const char* d = std::getenv("CPPWRF_NTIEDTKE_DUMP");
    const char* ds = std::getenv("CPPWRF_NTIEDTKE_DUMP_STEP");
    if (d && ds && itimestep == std::atol(ds)) k.dump = std::fopen(d, "wb");
#endif
    column_for(g_.columns(), k);
#if !defined(CPPWRF_USE_CUDA)
    if (k.dump) std::fclose(static_cast<std::FILE*>(k.dump));
#endif
    return;
  }
  cu_kf::KFKernel<real> kf{};
  kf.shall = shall.view();
  kf.nx = g_.nx; kf.ny = g_.ny; kf.nz = g_.nz;
  kf.dt = real(g_.dt); kf.dx = real(g_.dx); kf.cudt = real(o_.cudt);
  kf.stepcu = stepcu; kf.ktau = int(itimestep); kf.trigger = o_.kfeta_trigger;
  kf.warm_rain = o_.warm_rain;
  kf.flag_qi = s.P_QI >= 0; kf.flag_qs = s.P_QS >= 0; kf.flag_qr = s.P_QR >= 0;
  kf.xlv0 = xlv0; kf.xlv1 = xlv1; kf.xls0 = xls0; kf.xls1 = xls1;
  kf.cp = cp; kf.r = r_d; kf.g = g; kf.ep1 = ep_1; kf.ep2 = ep_2;
  kf.svp1 = svp1; kf.svp2 = svp2; kf.svp3 = svp3; kf.svpt0 = svpt0;
  kf.tables = d_kf_tables_;
  kf.u = prep.u_phy.cview(); kf.v = prep.v_phy.cview(); kf.w = s.w_2.cview();
  kf.th = prep.th_phy.cview(); kf.t = prep.t_phy.cview(); kf.qv = s.q_2[s.P_QV].cview();
  kf.dz8w = prep.dz8w.cview(); kf.pcps = prep.p_hyd.cview();  // P=grid%p_hyd in first_rk_step_part1
  kf.rho = prep.rho.cview(); kf.pi = prep.pi_phy.cview();
  kf.w0avg = w0avg.view();
  kf.raincv = raincv.view(); kf.pratec = pratec.view(); kf.nca = nca.view();
  kf.cutop = cutop.view(); kf.cubot = cubot.view();
  kf.rthcuten = pt.rthcuten.view(); kf.rqvcuten = pt.rqvcuten.view();
  kf.rqccuten = pt.rqccuten.view(); kf.rqrcuten = pt.rqrcuten.view();
  kf.rqicuten = pt.rqicuten.view(); kf.rqscuten = pt.rqscuten.view();
  kf.cldfra_dp_kf = cldfra_dp.view(); kf.cldfra_sh_kf = cldfra_sh.view();
  kf.qc_kf = qc_cu.view(); kf.qi_kf = qi_cu.view();
  cu_kf::run_kf(kf, g_.columns());
}

namespace {
// module_convtrans_prep.F's INTEGER, SAVE :: stepave_count: one counter for all domains, not in the restart
// file (0 in a new process), reset by any domain's first step (ktau <= 1)
int stepave_count = 0;
}  // namespace

// dyn_em/module_convtrans_prep.F convtrans_prep (called in first_rk_step_part1 when cu_diag == 1)
void PhysicsSuite::convtrans_prep(const dyn::State& s) {
  const int ktau = int(itimestep);
  const int stepave = int(real(o_.convtrans_avglen_m) * real(60.) / real(g_.dt));
  if (stepave < 1)  // WRF: MOD(ktau, stepave) divides by zero
    throw std::runtime_error("convtrans_prep: convtrans_avglen_m is shorter than the time step");
  const bool avg_end_flag = (ktau % stepave == 0);
  const bool first_period_flag = (ktau <= stepave);
  if (ktau <= 1) stepave_count = 0;
  stepave_count = stepave_count + 1;
  KConvtransPrep k{};
  k.nz = g_.nz; k.ktau = ktau; k.stepave_count = stepave_count; k.cu_phys = o_.cu_physics;
  k.avg_end_flag = avg_end_flag; k.first_period_flag = first_period_flag;
  k.flag_qc = s.P_QC >= 0; k.flag_qi = s.P_QI >= 0;
  k.raincv = raincv.cview();
  k.raincv_a = raincv_a.view(); k.raincv_b = raincv_b.view();
  k.gd_cloud = gd_cloud.cview(); k.gd_cloud2 = gd_cloud2.cview();
  k.t_phy = prep.t_phy.cview(); k.p_phy = prep.p_phy.cview();
  k.qv = s.q_2[s.P_QV].cview();
  k.qc = s.P_QC >= 0 ? s.q_2[s.P_QC].cview() : zero3_.cview();
  k.qi = s.P_QI >= 0 ? s.q_2[s.P_QI].cview() : zero3_.cview();
  k.gd_cloud_a = gd_cloud_a.view(); k.gd_cloud2_a = gd_cloud2_a.view();
  k.qc_cu = qc_cu.view(); k.qi_cu = qi_cu.view(); k.gd_cldfr = gd_cldfr.view();
  column_for(g_.columns(), k);
  if (avg_end_flag) stepave_count = 0;
}

// module_physics_addtendc.F advance_ppt (called in solve_em before phy_prep_part2)
void PhysicsSuite::advance_ppt(const dyn::State& s) {
  if (o_.cu_physics == 0) return;
  column_for(g_.columns(), KAdvancePpt{pt.rthcuten.view(), pt.rqvcuten.view(), pt.rqccuten.view(), pt.rqrcuten.view(),
                                       pt.rqicuten.view(), pt.rqscuten.view(), rainc.view(), rainsh.view(),
                                       pratec.view(), pratesh.view(), nca.view(), htop.view(), hbot.view(),
                                       cuppt.view(), cutop.cview(), cubot.cview(), real(g_.dt), g_.nz, g_.nz + 1,
                                       s.P_QI >= 0, s.P_QS >= 0, o_.cu_physics == 1});
}

}  // namespace cppwrf::phys::driver
