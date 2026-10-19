// PhysicsSuite::pbl_driver (WRF phys/module_pbl_driver.F)
#include "phys/driver/suite.hpp"

#include "phys/bl_gwdo/gwdo.hpp"
#include "phys/bl_myj/myjpbl.hpp"
#include "phys/bl_ysu/ysu.hpp"
#include "phys/driver/wrf_constants.hpp"

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

}  // namespace

// module_pbl_driver.F pbl_driver for YSUSCHEME / MYJPBLSCHEME
void PhysicsSuite::pbl_driver(const dyn::State& s, const dyn::BaseState& b) {
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
                 /*have_ctopo=*/true, /*have_surface_diag=*/false, cp, g, rcp, r_d, rovg, ep_1, ep_2, karman,
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
}

}  // namespace cppwrf::phys::driver
