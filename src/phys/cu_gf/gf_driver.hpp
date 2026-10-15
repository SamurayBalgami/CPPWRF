#pragma once
// module_cu_gf_wrfdrv.F GFDRV for one column (ideep = 1, imid_gf = 0, ichoice from clos_choice,
// ishallow_g3 from ishallow). The driver's constants come from module_gfs_physcons (kind_phys = real*8
// even in a single-precision WRF): g = 9.80665, cp = 1004.6, xlv = 2.5e6; expressions with them are
// evaluated in double, as gfortran does, and stored in T.
#include "phys/cu_gf/gf_deep.hpp"
#include "phys/cu_gf/gf_shallow.hpp"

namespace cppwrf::phys::cu_gf {

inline constexpr double gfs_g = 9.80665, gfs_cp = 1004.6, gfs_xlv = 2.5e6;

// Column inputs of GFDRV in WRF order (k = 0..nz-1 C++ levels, w and p8w at nz+1 interfaces)
template <typename T>
struct GfDrvColumn {
  int nz;
  T dt, dx, ht, hfx, qfx, xland;
  int kpbl;
  // accessors (C++ level index k = 0..nz-1, or 0..nz for w/p8w)
  T u[KM], v[KM], t[KM], q[KM], p[KM], pi[KM], rho[KM], dz8w[KM], w[KM], p8w[KM];
  T rthften[KM], rthraten[KM], rthblten[KM], rqvften[KM], rqvblten[KM];
};

// GFDRV's arrays for one column (1-based levels)
template <typename T>
struct GfPrep {
  T zo[KM], po[KM], rhoi[KM], us[KM], vs[KM], t2d[KM], q2d[KM], tn[KM], qo[KM], tshall[KM], qshall[KM];
  T dhdt[KM], omeg[KM];
  T ter11, psur;
};

template <typename T>
CPPWRF_HD inline void gf_prep(const GfDrvColumn<T>& c, GfPrep<T>& p) {
  const int ktf = c.nz;
  for (int k = 0; k < KM; ++k) {
    p.zo[k] = 0; p.po[k] = 0; p.rhoi[k] = 0; p.us[k] = 0; p.vs[k] = 0; p.t2d[k] = 0; p.q2d[k] = 0;
    p.tn[k] = 0; p.qo[k] = 0; p.tshall[k] = 0; p.qshall[k] = 0; p.dhdt[k] = 0; p.omeg[k] = 0;
  }
  p.psur = c.p8w[0] * T(.01);
  p.ter11 = fmx(T(0.), c.ht);
  p.zo[1] = p.ter11 + T(.5) * c.dz8w[0];
  for (int k = 2; k <= ktf; ++k) p.zo[k] = p.zo[k - 1] + T(.5) * (c.dz8w[k - 2] + c.dz8w[k - 1]);
  for (int k = 1; k <= ktf; ++k) {
    const int m = k - 1;
    p.po[k] = c.p[m] * T(.01);
    p.rhoi[k] = c.rho[m];
    p.us[k] = c.u[m];
    p.vs[k] = c.v[m];
    p.t2d[k] = c.t[m];
    p.q2d[k] = c.q[m];
    if (p.q2d[k] < T(1.E-08)) p.q2d[k] = T(1.E-08);
    p.tn[k] = p.t2d[k] + (c.rthften[m] + c.rthraten[m] + c.rthblten[m]) * c.pi[m] * c.dt;
    p.qo[k] = p.q2d[k] + (c.rqvften[m] + c.rqvblten[m]) * c.dt;
    p.tshall[k] = p.t2d[k] + c.rthblten[m] * c.pi[m] * c.dt;
    p.dhdt[k] = T(gfs_cp * double(c.rthblten[m]) * double(c.pi[m]) + gfs_xlv * double(c.rqvblten[m]));
    p.qshall[k] = p.q2d[k] + c.rqvblten[m] * c.dt;
    if (p.tn[k] < T(200.)) p.tn[k] = p.t2d[k];
    if (p.qo[k] < T(1.E-08)) p.qo[k] = T(1.E-08);
  }
  for (int k = 1; k <= ktf; ++k) p.omeg[k] = T(-(gfs_g * double(c.rho[k - 1]) * double(c.w[k - 1])));
}

template <typename T>
CPPWRF_HD inline GfShallowIn<T> gf_shallow_in(const GfDrvColumn<T>& c, const GfPrep<T>& p) {
  GfShallowIn<T> s{};
  s.nz = c.nz;
  s.t = p.t2d; s.q = p.q2d; s.tn = p.tshall; s.qo = p.qshall; s.po = p.po; s.rho = p.rhoi; s.zo = p.zo;
  s.dhdt = p.dhdt;
  s.z1 = p.ter11; s.psur = p.psur; s.hfx = c.hfx; s.qfx = c.qfx; s.xland = c.xland; s.dtime = c.dt;
  s.kpbl = c.kpbl;
  return s;
}

// GFDRV outputs for one column (C++ levels 0..nz-1)
template <typename T>
struct GfDrvOut {
  T rthcuten[KM], rqvcuten[KM], rqccuten[KM], rqicuten[KM], rucuten[KM], rvcuten[KM];
  T raincv, pratec, htop, hbot;
};

// The rest of GFDRV for one column, after the row scan of the shallow inversion layers.
// active = the column lies in ibegc..iendc x jbegc..jendc (else WRF leaves zero tendencies).
template <typename T>
CPPWRF_HD inline void gf_column(const GfDrvColumn<T>& c, int ishallow, int ichoice, GfInvInfo<T>& inv,
                                T sec2_garbage, bool active, GfDrvOut<T>& o, int dbg = 0) {
  const int kte = c.nz, ktf = c.nz;
  for (int k = 0; k < KM; ++k) {
    o.rthcuten[k] = 0; o.rqvcuten[k] = 0; o.rqccuten[k] = 0; o.rqicuten[k] = 0; o.rucuten[k] = 0; o.rvcuten[k] = 0;
  }
  o.hbot = T(kte);
  o.htop = T(1);
  o.raincv = T(0.);
  o.pratec = T(0.);
  GfPrep<T> p;
  gf_prep(c, p);
  T outt[KM] = {}, outq[KM] = {}, outqc[KM] = {}, outu[KM] = {}, outv[KM] = {}, cupclw[KM] = {}, cnvwt[KM] = {};
  T outts[KM] = {}, outqs[KM] = {}, outqcs[KM] = {}, outus[KM] = {}, outvs[KM] = {}, cupclws[KM] = {};
  T cnvwts[KM] = {}, zus[KM] = {};
  T pret = T(0.), prets = T(0.), pretm = T(0.);
  T cuten = T(0.), cutenm = T(0.), cutens = T(1.);
  if (ishallow == 0) cutens = T(0.);
  T xmbs = T(0.);
  if (ishallow == 1) {
    const GfShallowIn<T> s = gf_shallow_in(c, p);
    cup_gf_sh(s, 1, inv, sec2_garbage, zus, xmbs, outts, outqs, outqcs, cnvwts, prets, cupclws);
    if (xmbs <= T(0.)) cutens = T(0.);
    neg_check(true, outqs, outts, outus, outvs, outqcs, prets, ktf);
  }
  GfColumnIn<T> d{};
  d.nz = c.nz;
  d.t = p.t2d; d.q = p.q2d; d.tn = p.tn; d.qo = p.qo; d.po = p.po; d.us = p.us; d.vs = p.vs; d.rho = p.rhoi;
  d.omeg = p.omeg; d.zo = p.zo; d.dhdt = p.dhdt;
  d.z1 = p.ter11; d.psur = p.psur; d.hfx = c.hfx; d.qfx = c.qfx; d.xland = c.xland; d.dx = c.dx; d.dtime = c.dt;
  d.kpbl = c.kpbl; d.ichoice = ichoice; d.xmbs_in = xmbs; d.dbg = dbg;
  T xmb = T(0);
  const GfDeepOut r = cup_gf_deep(d, outu, outv, outt, outq, outqc, cupclw, cnvwt, pret, xmb);
  neg_check(false, outq, outt, outu, outv, outqc, pret, ktf);
  if (!active) return;
  int kbcon = r.kbcon, ktop = r.ktop;
  if (pret > T(0.)) {
    cuten = T(1.);
  } else {
    cuten = T(0.);
    kbcon = 0;
    ktop = 0;
  }
  const T outtm = T(0.), outqm = T(0.), outum = T(0.), outqcm = T(0.);  // imid_gf = 0
  for (int k = 1; k <= ktf; ++k) {
    o.rthcuten[k - 1] = (cutens * outts[k] + cutenm * outtm + cuten * outt[k]) / c.pi[k - 1];
    o.rqvcuten[k - 1] = cuten * outq[k] + cutens * outqs[k] + cutenm * outqm;
    o.rucuten[k - 1] = outum * cutenm + outu[k] * cuten;
    o.rvcuten[k - 1] = outum * cutenm + outv[k] * cuten;
  }
  if (pret > T(0.) || pretm > T(0.) || prets > T(0.)) {
    o.pratec = cuten * pret + cutenm * pretm + cutens * prets;
    o.raincv = o.pratec * c.dt;
    if (T(ktop) > o.htop) o.htop = T(ktop) + T(.001);
    if (T(kbcon) < o.hbot) o.hbot = T(kbcon) + T(.001);
  }
  for (int k = 1; k <= ktf; ++k) {
    const T qc = outqcm + outqcs[k] + outqc[k] * cuten;
    if (p.t2d[k] < T(258.)) {
      o.rqicuten[k - 1] = qc;
      o.rqccuten[k - 1] = T(0.);
    } else {
      o.rqicuten[k - 1] = T(0.);
      o.rqccuten[k - 1] = qc;
    }
  }
}

}  // namespace cppwrf::phys::cu_gf
