#pragma once
// Grell-Freitas-Li WRF driver: WRF 4.8.0 phys/GFL/WRF/module_cu_gfl_wrf.F90
// cu_grell_freitas_li for one column (imid = 0, sub3d = 0, dicycle = 0,
// nranflag = 0, do_capsuppress = 0, ccnclean = 250., ccn = 1500., tcrit = 258.).
//
// In module_cumulus_driver.F, cu_grell_freitas_li receives WRF's model constants
// (xlv = 2.5e6, cp = 1004.5, g = 9.81, r_v = 461.6) as real(kind=kind_phys).

#include "phys/cu_gfl/gfl_deep.hpp"
#include "phys/cu_gfl/gfl_shallow.hpp"

namespace cppwrf::phys::cu_gfl {

#define R(x) static_cast<real>(x)

// Column inputs of cu_grell_freitas_li in C++ level order (k = 0..nz-1, w and p8w at 0..nz interfaces)
template <typename real>
struct GflDrvColumn {
  int nz;
  real dt, dx, ht, hfx, qfx, xland;
  real xlv, cp, g, r_v;
  int kpbl;
  real u[KM], v[KM], t[KM], q[KM], p[KM], pi[KM], rho[KM], dz8w[KM], w[KM], p8w[KM];
  real rthften[KM], rthraten[KM], rthblten[KM], rqvften[KM], rqvblten[KM];
};

// Prepared 1-based column arrays in cu_grell_freitas_li
template <typename real>
struct GflPrep {
  real zo[KM], po[KM], rhoi[KM], us[KM], vs[KM], t2d[KM], q2d[KM], qcheck[KM];
  real qc[KM], qi[KM], tn[KM], qo[KM], tshall[KM], qshall[KM], dhdt[KM], omeg[KM];
  real ter11, psur, mconv;
};

template <typename real>
CPPWRF_HOST_DEVICE inline void gfl_prep(const GflDrvColumn<real>& c,
                                        GflPrep<real>& p) {
  const int ktf = c.nz;
  for (int k = 0; k < KM; ++k) {
    p.zo[k] = R(0.0);
    p.po[k] = R(0.0);
    p.rhoi[k] = R(0.0);
    p.us[k] = R(0.0);
    p.vs[k] = R(0.0);
    p.t2d[k] = R(0.0);
    p.q2d[k] = R(0.0);
    p.qcheck[k] = R(0.0);
    p.qc[k] = R(0.0);
    p.qi[k] = R(0.0);
    p.tn[k] = R(0.0);
    p.qo[k] = R(0.0);
    p.tshall[k] = R(0.0);
    p.qshall[k] = R(0.0);
    p.dhdt[k] = R(0.0);
    p.omeg[k] = R(0.0);
  }
  p.ter11 = max(R(0.0), c.ht);
  p.zo[1] = p.ter11 + R(0.5) * c.dz8w[0];
  for (int k = 2; k <= ktf; ++k) {
    p.zo[k] = p.zo[k - 1] + R(0.5) * (c.dz8w[k - 2] + c.dz8w[k - 1]);
  }
  p.psur = c.p8w[0] * R(0.01);
  for (int k = 1; k <= ktf; ++k) {
    const int m = k - 1;
    p.us[k] = c.u[m];
    p.vs[k] = c.v[m];
    p.rhoi[k] = c.rho[m];
    p.t2d[k] = c.t[m];
    p.q2d[k] = c.q[m];
    p.qcheck[k] = c.q[m];
    p.qc[k] = R(0.0);
    p.qi[k] = R(0.0);
    if (p.q2d[k] < R(1.e-08)) p.q2d[k] = R(1.e-08);
    p.tn[k] = p.t2d[k] +
              (c.rthften[m] + c.rthraten[m] + c.rthblten[m]) * c.pi[m] * c.dt;
    p.qo[k] = p.q2d[k] + (c.rqvften[m] + c.rqvblten[m]) * c.dt;
    if (p.tn[k] < R(200.0)) p.tn[k] = p.t2d[k];
    if (p.qo[k] < R(1.e-08)) p.qo[k] = R(1.e-08);
    p.po[k] = c.p[m] * R(0.01);
    p.tshall[k] = p.t2d[k] + c.rthblten[m] * c.pi[m] * c.dt;
    p.qshall[k] = p.q2d[k] + c.rqvblten[m] * c.dt;
    p.dhdt[k] = c.cp * c.rthblten[m] * c.pi[m] + c.xlv * c.rqvblten[m];
  }
  p.omeg[1] = R(0.0);
  for (int k = 2; k <= ktf; ++k) {
    p.omeg[k] = -c.g * R(0.5) * (c.rho[k - 1] + c.rho[k - 2]) * c.w[k - 1];
  }
  p.mconv = R(0.0);
  for (int k = 2; k <= ktf; ++k) {
    const real dq = p.q2d[k] - p.q2d[k - 1];
    p.mconv = p.mconv + p.omeg[k] * dq / c.g;
  }
}

// Outputs of cu_grell_freitas_li for one column (C++ levels 0..nz-1)
template <typename real>
struct GflDrvOut {
  real rthcuten[KM], rqvcuten[KM], rqccuten[KM], rqicuten[KM], rucuten[KM], rvcuten[KM];
  real raincv, pratec, htop, hbot;
};

template <typename real>
CPPWRF_HOST_DEVICE inline void gfl_column(const GflDrvColumn<real>& c,
                                          int ishallow, int ichoice_deep,
                                          GflDrvOut<real>& o) {
  using C = GflConstants<real>;
  const int kte = c.nz, ktf = c.nz;
  for (int k = 0; k < KM; ++k) {
    o.rthcuten[k] = R(0.0);
    o.rqvcuten[k] = R(0.0);
    o.rqccuten[k] = R(0.0);
    o.rqicuten[k] = R(0.0);
    o.rucuten[k] = R(0.0);
    o.rvcuten[k] = R(0.0);
  }
  o.hbot = static_cast<real>(kte);
  o.htop = R(1.0);
  o.raincv = R(0.0);
  o.pratec = R(0.0);

  GflPrep<real> p;
  gfl_prep(c, p);

  real outt[KM] = {}, outq[KM] = {}, outqc[KM] = {}, outu[KM] = {}, outv[KM] = {};
  real cupclw[KM] = {}, cnvwt[KM] = {};
  real subten_h[KM] = {}, subten_q[KM] = {}, subten_t[KM] = {};
  real subten_u[KM] = {}, subten_v[KM] = {};
  real outts[KM] = {}, outqs[KM] = {}, outqcs[KM] = {}, outus[KM] = {}, outvs[KM] = {};
  real cupclws[KM] = {}, cnvwts[KM] = {}, zus[KM] = {};

  real pret = R(0.0), prets = R(0.0);
  real xmb = R(0.0), xmbs = R(0.0);
  real frhs = R(0.0);

  GflColumnIn<real> d{};
  d.nz = c.nz;
  d.t = p.t2d;
  d.q = p.q2d;
  d.tn = p.tn;
  d.qo = p.qo;
  d.po = p.po;
  d.us = p.us;
  d.vs = p.vs;
  d.rho = p.rhoi;
  d.omeg = p.omeg;
  d.zo = p.zo;
  d.dhdt = p.dhdt;
  d.z1 = p.ter11;
  d.psur = p.psur;
  d.hfx = c.hfx;
  d.qfx = c.qfx;
  d.xland = c.xland;
  d.dx = c.dx;
  d.dtime = c.dt;
  d.kpbl = c.kpbl;
  d.ichoice = ichoice_deep;
  d.xmbs_in = R(0.0);

  const GflDeepOut r =
      cu_gfl_deep_run(d, outu, outv, outt, outq, outqc, subten_h, subten_q,
                      subten_t, subten_u, subten_v, cupclw, cnvwt, pret, xmb);

  neg_check(false, c.dt, p.qcheck, outq, subten_q, outt, subten_t, outu,
            subten_u, outv, subten_v, outqc, pret, r.ktop);

  GflShallowOut sh{0, 0, 0, 0};
  if (ishallow == 1) {
    if (r.ierr != 0) {
      sh.ierr = 999;
    } else {
      GflShallowIn<real> s{};
      s.nz = c.nz;
      s.us = p.us;
      s.vs = p.vs;
      s.zo = p.zo;
      s.t = p.t2d;
      s.q = p.q2d;
      s.tn = p.tshall;
      s.qo = p.qshall;
      s.po = p.po;
      s.dhdt = p.dhdt;
      s.rho = p.rhoi;
      s.z1 = p.ter11;
      s.psur = p.psur;
      s.hfx = c.hfx;
      s.qfx = c.qfx;
      s.xland = c.xland;
      s.dx = c.dx;
      s.tcrit = C::tcrit;
      s.dtime = c.dt;
      s.kpbl = c.kpbl;
      s.ichoice = 3;
      sh = cu_gfl_sh_run(s, frhs, zus, xmbs, outts, outqs, outqcs, outus, outvs,
                         cnvwts, prets, cupclws);
    }
  }

  real factor = R(1.0);
  if (pret > R(0.0) || prets > R(0.0)) {
    o.pratec = factor * (pret + prets);
    o.raincv = factor * (pret + prets) * c.dt;
    const int ktopmax = max(r.ktop, sh.ktop);
    const int kbconmax = max(r.kbcon, sh.kbcon);
    if (static_cast<real>(ktopmax) > o.htop)
      o.htop = static_cast<real>(ktopmax) + R(0.001);
    if (static_cast<real>(kbconmax) < o.hbot)
      o.hbot = static_cast<real>(kbconmax) + R(0.001);
  } else if (r.ierr > 0) {
    factor = R(0.0);
  }

  for (int k = 1; k <= kte; ++k) {
    o.rthcuten[k - 1] = factor * (outts[k] + outt[k]) / c.pi[k - 1];
    o.rqvcuten[k - 1] = factor * (outqs[k] + outq[k]);
    o.rucuten[k - 1] = factor * (outu[k] + outus[k]);
    o.rvcuten[k - 1] = factor * (outv[k] + outvs[k]);
    if (p.t2d[k] < C::tcrit) {
      o.rqccuten[k - 1] = R(0.0);
      o.rqicuten[k - 1] = factor * (outqcs[k] + outqc[k]);
    } else {
      o.rqicuten[k - 1] = R(0.0);
      o.rqccuten[k - 1] = factor * (outqcs[k] + outqc[k]);
    }
  }
}

#undef R

}  // namespace cppwrf::phys::cu_gfl
