// PhysicsSuite::pbl_driver (WRF phys/module_pbl_driver.F)
#include "phys/driver/suite.hpp"

#include "phys/bl_gwdo/gwdo.hpp"
#include "phys/bl_myj/myjpbl.hpp"
#include "phys/bl_mynn/mynnedmf.hpp"
#include "phys/bl_ysu/ysu.hpp"
#include "phys/driver/wrf_constants.hpp"
#include "phys/fogdes/fogdes.hpp"

#include "phys/driver/suite_internal.hpp"

namespace cppwrf::phys::driver {

using namespace suite_detail;

namespace {

// ---------------------------------------------------------------- pbl_driver
struct KPblPrep {
  C3 u_phy, v_phy, p8w;
  V3 u_phytmp, v_phytmp, rthblten, rublten, rvblten, rqcblten, rqvblten, rqiblten;
  V2 psfc;
  int nz;
  bool flag_qi;
  CPPWRF_HD void operator()(int i, int j) const {
    constexpr real u_frame = real(0.), v_frame = real(0.);
    for (int k = 0; k < nz; ++k) {
      v_phytmp(i, j, k) = v_phy(i, j, k) + v_frame;
      u_phytmp(i, j, k) = u_phy(i, j, k) + u_frame;
    }
    psfc(i, j) = p8w(i, j, 0);
    for (int k = 0; k < nz; ++k) {  // kts..min(kte+1,kde): the field has nz levels, level kde is never read
      rthblten(i, j, k) = real(0.);
      rublten(i, j, k) = real(0.);
      rvblten(i, j, k) = real(0.);
      rqcblten(i, j, k) = real(0.);
      rqvblten(i, j, k) = real(0.);
    }
    if (flag_qi)
      for (int k = 0; k < nz; ++k) rqiblten(i, j, k) = real(0.);
  }
};

// module_pbl_driver.F gwdo (gwd_opt = 1) after the PBL scheme: module_bl_gwdo.F -> bl_gwdo_run
struct KGwdo {
  int nz;
  real dt, dx_factor, cp, g, rd, ep1, pi;
  bool if_nonhyd;
  C3 u, v, t, qv, p, p8w, pi3d, z;
  V3 rublten, rvblten, dtaux3d, dtauy3d;
  C2 var, oc1, oa1, oa2, oa3, oa4, ol1, ol2, ol3, ol4, omax, dx;
  BView2<const real> sina, cosa;
  V2 dusfcg, dvsfcg;
  CPPWRF_HD void operator()(int i, int j) const {
    bl_gwdo::GwdoColumn<real> c;
    for (int k = 0; k <= nz; ++k) c.prsi[k + 1] = p8w(i, j, k);
    for (int k = 0; k < nz; ++k) {
      c.rublten[k + 1] = rublten(i, j, k);
      c.rvblten[k + 1] = rvblten(i, j, k);
      c.uproj[k + 1] = u(i, j, k);
      c.vproj[k + 1] = v(i, j, k);
      c.t1[k + 1] = t(i, j, k);
      c.q1[k + 1] = qv(i, j, k);
      c.prsl[k + 1] = p(i, j, k);
      c.prslk[k + 1] = pi3d(i, j, k);
      c.zl[k + 1] = z(i, j, k);
    }
    c.dxmeter = dx(i, j); c.var = var(i, j); c.oc1 = oc1(i, j); c.sina = sina(i, j); c.cosa = cosa(i, j);
    c.oa2d1 = oa1(i, j); c.oa2d2 = oa2(i, j); c.oa2d3 = oa3(i, j); c.oa2d4 = oa4(i, j);
    c.ol2d1 = ol1(i, j); c.ol2d2 = ol2(i, j); c.ol2d3 = ol3(i, j); c.ol2d4 = ol4(i, j);
    c.omax = omax(i, j);
    bl_gwdo::bl_gwdo_run<real>(c, nz, dx_factor, if_nonhyd, g, cp, rd, ep1, pi, dt);
    for (int k = 0; k < nz; ++k) {
      rublten(i, j, k) = c.rublten[k + 1];
      rvblten(i, j, k) = c.rvblten[k + 1];
      dtaux3d(i, j, k) = c.dtaux3d[k + 1];
      dtauy3d(i, j, k) = c.dtauy3d[k + 1];
    }
    dusfcg(i, j) = c.dusfcg;
    dvsfcg(i, j) = c.dvsfcg;
  }
};

// module_bl_mynnedmf_driver.F90 mynnedmf_driver for one column (j = row): the surface-flux limits,
// mynnedmf_pre_run, mynnedmf and mynnedmf_post_run (flag_qc = flag_qi = flag_qs = .true.); the model levels
// k = 0 .. nz-1 are the scheme's k = 1 .. kte
struct KMynnPbl {
  int nz, initflag;
  mynnedmf::Options op;
  int tke_budget = 0, bl_mynn_output = 0;
  V3 edmf_a, edmf_w, edmf_thl, edmf_qt, edmf_ent, edmf_qc;  // bl_mynn_output > 0
  V3 qshear, qbuoy, qdiss, qwt, dqke;                        // tke_budget = 1
  real delt, dx;
  C3 dz8w, u, v, w, th, t, p, pi, rho, qv, qc, qi, qs;
  C2 xland, tsk, psfc, ust, hfx, qfx, wspd, uoce, voce;
  V3 qke, el, sh3d, sm3d, exch_h, exch_m, tsq, qsq, cov, qc_bl, qi_bl, cldfra_bl;
  V3 rublten, rvblten, rthblten, rqvblten, rqcblten, rqiblten;
  V2 pblh, maxwidth, maxmf, ztop_plume, excess_h, excess_q;
  View2<int> kpbl;
  V3 qke_adv_1, qke_adv_2;  // the scalar qke_adv (State q_1 / q_2), when present
  bool has_qke_adv = false;
  CPPWRF_HD void operator()(int i, int j) const {
    using namespace mynnedmf;
    Column<kMaxLevels> c;
    real qv1[kMaxLevels + 2], qc1[kMaxLevels + 2], qi1[kMaxLevels + 2], qs1[kMaxLevels + 2];
    c.kte = nz;
    for (int k = 1; k <= nz; ++k) {
      const int kk = k - 1;
      c.dz[k] = dz8w(i, j, kk); c.u[k] = u(i, j, kk); c.v[k] = v(i, j, kk); c.w[k] = w(i, j, kk);
      c.th[k] = th(i, j, kk); c.tk[k] = t(i, j, kk); c.pres[k] = p(i, j, kk); c.ex[k] = pi(i, j, kk);
      c.rho[k] = rho(i, j, kk);
      qv1[k] = qv(i, j, kk); qc1[k] = qc(i, j, kk); qi1[k] = qi(i, j, kk); qs1[k] = qs(i, j, kk);
    }
    real hfx1 = hfx(i, j);
    if (hfx1 > real(1200.)) hfx1 = real(1200.);
    if (hfx1 < real(-600.)) hfx1 = real(-600.);
    real qfx1 = qfx(i, j);
    if (qfx1 > real(9e-4)) qfx1 = real(9e-4);
    if (qfx1 < real(-3e-4)) qfx1 = real(-3e-4);
    c.pblh = pblh(i, j);
    c.kpbl = kpbl(i, j);
    if (initflag == 0) {  // .or. restart
      for (int k = 1; k <= nz; ++k) {
        const int kk = k - 1;
        c.qc_bl[k] = qc_bl(i, j, kk); c.qi_bl[k] = qi_bl(i, j, kk); c.cldfra_bl[k] = cldfra_bl(i, j, kk);
        c.qke[k] = qke(i, j, kk); c.qsq[k] = qsq(i, j, kk); c.tsq[k] = tsq(i, j, kk); c.cov[k] = cov(i, j, kk);
        c.sh[k] = sh3d(i, j, kk); c.sm[k] = sm3d(i, j, kk); c.kh[k] = exch_h(i, j, kk); c.km[k] = exch_m(i, j, kk);
        c.el[k] = el(i, j, kk);
        // bl_mynn_tkeadvect: the advected scalar qke_adv (grid%scalar), else qke
        c.qke_adv[k] = op.tkeadvect ? qke_adv_2(i, j, kk) : qke(i, j, kk);
      }
    }
    for (int k = 1; k <= nz; ++k) {  // mynnedmf_pre_run
      c.sqv[k] = qv1[k] / (one + qv1[k]);
      c.sqc[k] = qc1[k] / (one + qv1[k]);
      c.sqi[k] = qi1[k] / (one + qv1[k]);
      c.sqs[k] = qs1[k] / (one + qv1[k]);
    }
    mynnedmf::mynnedmf<kMaxLevels>(c, op, tke_budget, initflag, delt, dx, xland(i, j), tsk(i, j), psfc(i, j),
                                   ust(i, j), hfx1, qfx1, wspd(i, j), uoce(i, j), voce(i, j));
    for (int k = 1; k <= nz; ++k) {  // mynnedmf_post_run
      const real sq = qv1[k] / (one + qv1[k]);
      const real sqv = sq + c.dqv[k] * delt;
      const real rq = sqv / (one - sqv);
      c.dqv[k] = (rq - qv1[k]) / delt;
      real sqx = qc1[k] / (one + qv1[k]);
      real sqn = sqx + c.dqc[k] * delt;
      c.dqc[k] = (sqn * (one + sqv) - qc1[k]) / delt;
      sqx = qi1[k] / (one + qv1[k]);
      sqn = sqx + c.dqi[k] * delt;
      c.dqi[k] = (sqn * (one + sqv) - qi1[k]) / delt;
    }
    for (int k = 1; k <= nz; ++k) {
      const int kk = k - 1;
      qke(i, j, kk) = c.qke[k]; el(i, j, kk) = c.el[k]; sh3d(i, j, kk) = c.sh[k]; sm3d(i, j, kk) = c.sm[k];
      exch_h(i, j, kk) = c.kh[k]; exch_m(i, j, kk) = c.km[k]; tsq(i, j, kk) = c.tsq[k]; qsq(i, j, kk) = c.qsq[k];
      cov(i, j, kk) = c.cov[k];
      rublten(i, j, kk) = c.du[k]; rvblten(i, j, kk) = c.dv[k]; rthblten(i, j, kk) = c.dth[k];
      rqvblten(i, j, kk) = c.dqv[k]; rqcblten(i, j, kk) = c.dqc[k]; rqiblten(i, j, kk) = c.dqi[k];
      qc_bl(i, j, kk) = c.qc_bl[k] / (one - c.sqv[k]);
      qi_bl(i, j, kk) = c.qi_bl[k] / (one - c.sqv[k]);
      cldfra_bl(i, j, kk) = c.cldfra_bl[k];
      if (tke_budget == 1) {
        dqke(i, j, kk) = c.dqke[k]; qwt(i, j, kk) = c.qwt[k]; qshear(i, j, kk) = c.qshear[k];
        qbuoy(i, j, kk) = c.qbuoy[k]; qdiss(i, j, kk) = c.qdiss[k];
      }
      if (bl_mynn_output > 0) {
        edmf_a(i, j, kk) = c.edmf_a[k]; edmf_w(i, j, kk) = c.edmf_w[k]; edmf_qt(i, j, kk) = c.edmf_qt[k];
        edmf_thl(i, j, kk) = c.edmf_thl[k]; edmf_ent(i, j, kk) = c.edmf_ent[k]; edmf_qc(i, j, kk) = c.edmf_qc[k];
      }
      // qke_adv(i,k,j) = qke_adv1(k) into grid%scalar; rk_update_scalar (rk_step 1) then copies it to the old level
      if (has_qke_adv) qke_adv_2(i, j, kk) = qke_adv_1(i, j, kk) = c.qke_adv[k];
    }
    kpbl(i, j) = c.kpbl;
    pblh(i, j) = c.pblh;
    if (op.edmf > 0) {
      maxwidth(i, j) = c.maxwidth; maxmf(i, j) = c.maxmf; ztop_plume(i, j) = c.ztop_plume;
      excess_h(i, j) = c.excess_h; excess_q(i, j) = c.excess_q;
    }
  }
};


// module_pbl_driver.F bl_fogdes (grav_settling >= 1) after the PBL scheme and the gravity wave drag:
// the settling of qc_curr with the deposition velocity VDFG of the surface driver, added to RQCBLTEN
struct KBlFogdes {
  C3 qc, dz8w;
  C2 vdfg;
  V3 rqcblten;
  int grav_settling, nz;
  CPPWRF_HD void operator()(int i, int j) const {
    auto q = [&](int k) { return qc(i, j, k); };
    auto dz = [&](int k) { return dz8w(i, j, k); };
    auto dqc = [&](int k) -> real& { return rqcblten(i, j, k); };
    fogdes::bl_fogdes(grav_settling, vdfg(i, j), nz, q, dz, dqc);
  }
};
}  // namespace

// module_pbl_driver.F pbl_driver for YSUSCHEME / MYJPBLSCHEME
void PhysicsSuite::pbl_driver(dyn::State& s, const dyn::BaseState& b) {
  if (o_.bl_pbl_physics == 0) return;  // module_pbl_driver: no PBL scheme (surface fluxes reach diff_opt 2)
  const Range2D cols = g_.columns();
  const bool run_param = itimestep == 1 || o_.bldt == 0 || stepbl == 1 || itimestep % stepbl == 0;
  if (!run_param) return;
  const real dt = real(g_.dt);
  const real dtbl = o_.bldt == 0 ? dt : dt * real(stepbl);
  const bool flag_qc = s.P_QC >= 0, flag_qi = s.P_QI >= 0;
  column_for(cols, KPblPrep{prep.u_phy.cview(), prep.v_phy.cview(), prep.p_hyd_w.cview(), u_phytmp.view(),
                            v_phytmp.view(), pt.rthblten.view(), pt.rublten.view(), pt.rvblten.view(),
                            pt.rqcblten.view(), pt.rqvblten.view(), pt.rqiblten.view(), psfc.view(), g_.nz,
                            flag_qi || o_.bl_pbl_physics == 2});
  if (o_.bl_pbl_physics == 1) {
    using namespace wrfc;
    // YSU: U3D=u_phytmp, T3D=t_phy, P3D=p_phy (grid%p_hyd), P3DI=p8w (grid%p_hyd_w), PI3D=pi_phy,
    // PSIM=fm, PSIH=fhh, DT=dtbl
    ysu::run_ysu(u_phytmp.cview(), v_phytmp.cview(), prep.t_phy.cview(), s.q_2[s.P_QV].cview(),
                 flag_qc ? s.q_2[s.P_QC].cview() : zero3_.cview(), flag_qi ? s.q_2[s.P_QI].cview() : zero3_.cview(),
                 prep.p_hyd.cview(), prep.p_hyd_w.cview(), prep.pi_phy.cview(), prep.dz8w.cview(),
                 pt.rthraten.cview(), psfc.cview(), znt.cview(), ust.cview(), fm.cview(), fh.cview(), xland.cview(),
                 hfx.cview(), qfx.cview(), wspd.cview(), br.cview(), u10.view(), v10.view(), uoce.cview(),
                 voce.cview(), ctopo.cview(), ctopo2.cview(), pt.rublten.view(), pt.rvblten.view(),
                 pt.rthblten.view(), pt.rqvblten.view(), pt.rqcblten.view(), pt.rqiblten.view(), exch_h.view(),
                 exch_m.view(), pblh.view(), kpbl.view(), wstar.view(), delta.view(), hol.view(), hol.view(),
                 hol.view(), hol.view(), dtbl, cols, g_.nz, flag_qc, flag_qi, o_.ysu_topdown_pblmix == 1,
                 /*have_surface_diag=*/false, cp, g, rcp, r_d, rovg, ep_1, ep_2, karman,
                 xlv, r_v);
  } else if (o_.bl_pbl_physics == 2) {
    // MYJPBL (module_pbl_driver.F:1447-1466: DT=dt, STEPBL=stepbl, PMID=p_phy, PINT=p8w, QCI=qi_curr always PRESENT)
    const C2 ht_v = b.ht.uniform() ? slice0(zero3_) : C2{b.ht.data(), int(b.ht.layout().sj)};
    myjpbl::run_myjpbl(
        ht_v, prep.dz8w.cview(), prep.p_hyd.cview(), prep.p_hyd_w.cview(), prep.th_phy.cview(), prep.t_phy.cview(),
        prep.pi_phy.cview(), s.q_2[s.P_QV].cview(), flag_qc ? s.q_2[s.P_QC].cview() : zero3_.cview(),
        flag_qi ? s.q_2[s.P_QI].cview() : zero3_.cview(), prep.u_phy.cview(), prep.v_phy.cview(), prep.rho.cview(),
        tsk.cview(), chklowq.cview(), thz0.cview(), qz0.cview(), uz0.cview(), vz0.cview(), lowlyr.cview(),
        xland.cview(), xice.cview(), snow.cview(), ust.cview(), znt.cview(), akhs.cview(), akms.cview(), lh.cview(),
        qsfc.view(), thz0.view(), qz0.view(), ct.view(), tke_pbl.view(), exch_h.view(), el_pbl.view(), pblh.view(),
        kpbl.view(), mixht.view(), pt.rublten.view(), pt.rvblten.view(), pt.rthblten.view(), pt.rqvblten.view(),
        pt.rqcblten.view(), pt.rqiblten.view(), cols, dt, stepbl, g_.nz, /*has_qi=*/true);
  } else if (o_.bl_pbl_physics == 5) {
    // MYNNPBLSCHEME (module_pbl_driver.F: initflag = 1 at itimestep 1, delt = dtbl, dz = dz8w, dxc = dx,
    // u = u_phy, v = v_phy, w = w_2, p = p_phy, exner = pi_phy, T3D = t_phy, ts = tsk, ps = psfc)
    if (!flag_qc || !flag_qi || s.P_QS < 0)
      throw std::runtime_error("bl_pbl_physics=5: the microphysics must carry qc, qi and qs");
    KMynnPbl k{};
    k.nz = g_.nz; k.initflag = itimestep == 1 ? 1 : 0; k.delt = dtbl; k.dx = real(g_.dx);
    k.op = o_.mynn; k.tke_budget = o_.tke_budget; k.bl_mynn_output = o_.bl_mynn_output;
    if (o_.tke_budget == 1) {
      k.qshear = qshear.view(); k.qbuoy = qbuoy.view(); k.qdiss = qdiss.view(); k.qwt = qwt.view();
      k.dqke = dqke.view();
    }
    if (o_.bl_mynn_output > 0) {
      k.edmf_a = edmf_a.view(); k.edmf_w = edmf_w.view(); k.edmf_thl = edmf_thl.view(); k.edmf_qt = edmf_qt.view();
      k.edmf_ent = edmf_ent.view(); k.edmf_qc = edmf_qc.view();
    }
    k.dz8w = prep.dz8w.cview(); k.u = prep.u_phy.cview(); k.v = prep.v_phy.cview(); k.w = s.w_2.cview();
    k.th = prep.th_phy.cview(); k.t = prep.t_phy.cview(); k.p = prep.p_hyd.cview(); k.pi = prep.pi_phy.cview();
    k.rho = prep.rho.cview(); k.qv = s.q_2[s.P_QV].cview(); k.qc = s.q_2[s.P_QC].cview();
    k.qi = s.q_2[s.P_QI].cview(); k.qs = s.q_2[s.P_QS].cview();
    k.xland = xland.cview(); k.tsk = tsk.cview(); k.psfc = psfc.cview(); k.ust = ust.cview(); k.hfx = hfx.cview();
    k.qfx = qfx.cview(); k.wspd = wspd.cview(); k.uoce = uoce.cview(); k.voce = voce.cview();
    k.qke = qke.view(); k.el = el_pbl.view(); k.sh3d = sh3d.view(); k.sm3d = sm3d.view(); k.exch_h = exch_h.view();
    k.exch_m = exch_m.view(); k.tsq = tsq.view(); k.qsq = qsq.view(); k.cov = cov.view(); k.qc_bl = qc_bl.view();
    k.qi_bl = qi_bl.view(); k.cldfra_bl = cldfra_bl.view();
    k.rublten = pt.rublten.view(); k.rvblten = pt.rvblten.view(); k.rthblten = pt.rthblten.view();
    k.rqvblten = pt.rqvblten.view(); k.rqcblten = pt.rqcblten.view(); k.rqiblten = pt.rqiblten.view();
    k.pblh = pblh.view(); k.maxwidth = maxwidth.view(); k.maxmf = maxmf.view(); k.ztop_plume = ztop_plume.view();
    k.excess_h = excess_h.view(); k.excess_q = excess_q.view(); k.kpbl = kpbl.view();
    k.has_qke_adv = s.P_QKE_ADV >= 0;
    if (o_.mynn.tkeadvect && !k.has_qke_adv)
      throw std::runtime_error("bl_mynn_tkeadvect: the scalar qke_adv is missing");
    if (k.has_qke_adv) { k.qke_adv_1 = s.q_1[s.P_QKE_ADV].view(); k.qke_adv_2 = s.q_2[s.P_QKE_ADV].view(); }
    column_for(cols, k);
  }
  if (o_.gwd_opt == 1) {
    using namespace wrfc;
    KGwdo k{};
    k.nz = g_.nz; k.dt = dtbl; k.dx_factor = real(o_.gwd_dx_factor); k.if_nonhyd = o_.gwd_if_nonhyd;
    k.cp = cp; k.g = g; k.rd = r_d; k.ep1 = ep_1; k.pi = real(3.141592653);
    k.u = prep.u_phy.cview(); k.v = prep.v_phy.cview(); k.t = prep.t_phy.cview(); k.qv = s.q_2[s.P_QV].cview();
    k.p = prep.p_hyd.cview(); k.p8w = prep.p_hyd_w.cview(); k.pi3d = prep.pi_phy.cview(); k.z = prep.z.cview();
    k.rublten = pt.rublten.view(); k.rvblten = pt.rvblten.view();
    k.dtaux3d = dtaux3d.view(); k.dtauy3d = dtauy3d.view();
    k.var = var2d.cview(); k.oc1 = oc12d.cview(); k.oa1 = oa1.cview(); k.oa2 = oa2.cview(); k.oa3 = oa3.cview();
    k.oa4 = oa4.cview(); k.ol1 = ol1.cview(); k.ol2 = ol2.cview(); k.ol3 = ol3.cview(); k.ol4 = ol4.cview();
    k.omax = elvmax.cview(); k.dx = dx2d.cview(); k.sina = b.sina.view(); k.cosa = b.cosa.view();
    k.dusfcg = dusfcg.view(); k.dvsfcg = dvsfcg.view();
    column_for(cols, k);
  }
  if (o_.grav_settling >= 1) {
    if (!flag_qc) throw std::runtime_error("grav_settling: the microphysics must carry qc (qc_curr)");
    column_for(cols, KBlFogdes{s.q_2[s.P_QC].cview(), prep.dz8w.cview(), vdfg.cview(), pt.rqcblten.view(),
                               o_.grav_settling, g_.nz});
  }
}

}  // namespace cppwrf::phys::driver
