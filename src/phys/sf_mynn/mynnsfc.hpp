#pragma once
// MYNN surface layer driver: WRF v4.8.0 phys/MYNN-SFC/WRF/module_sf_mynnsfc_driver.F90 (mynnsfc_driver) as
// called by module_surface_driver.F (sf_sfclay_physics = 5, fractional_seaice = 0): per point the land
// (no snow), water and ice (snow >= 5 cm on land, partial sea ice) modules in WRF's order, each starting
// from the point state the previous one left; qv, qsfc, qgh, qfx, lh, q2 as mixing ratios outside the modules.
#include "core/field.hpp"
#include "phys/sf_mynn/mynnsfc_ice.hpp"

namespace cppwrf::phys::mynnsfc {

struct MynnSfcKernel {
  // level 0 (kts) and 1 (kts+1) of the 3-D inputs
  View3<const real> u3d, v3d, t3d, qv3d, p3d, dz8w, rho3d;
  View2<const real> psfcpa, mavail, xland, tsk, snowh, dx, xice, lakemask;
  View2<real> pblh;  // INOUT (set at the first step)
  // INOUT / OUT
  View2<real> chs, chs2, cqs, cqs2, cpm, znt, ust, ustm, zol, mol, rmol, psim, psih, qgh, hfx, qfx, lh, flhc, flqc,
      qsfc, u10, v10, th2, t2, q2, gz1oz0, wspd, br, ch, ck, cka, cd, cda;
  const PsiTables* tb = nullptr;  // psi_init tables (device memory in CUDA builds)
  int itimestep = 1;
  bool restart = false;
  int isfflx = 1;
  real xice_threshold = real(0.5);

  CPPWRF_HD void operator()(int i, int j) const {
    constexpr real snow_thresh = real(0.05);
    if (itimestep == 1) {  // .not. restart .or. .not. cycling: always at the first step
      if (ust(i, j) <= real(1e-4) || ust(i, j) > real(4.0)) {
        const real u_1 = u3d(i, j, 0), v_1 = v3d(i, j, 0);
        ust(i, j) = c_max(real(0.04) * std::sqrt(u_1 * u_1 + v_1 * v_1), real(0.001));
      }
      pblh(i, j) = c_max(ten, real(700.) * ust(i, j));
      qfx(i, j) = zero;
      hfx(i, j) = zero;
      mol(i, j) = zero;
    }
    const real qvmr = qv3d(i, j, 0);
    SfcIn in;
    in.itimestep = itimestep;
    in.flag_restart = restart;
    in.isfflx = isfflx;
    in.dz8w_1 = dz8w(i, j, 0);
    in.qv_1 = c_max(real(1e-8), qvmr / (one + qvmr));
    in.p_1 = p3d(i, j, 0);
    in.t_1 = t3d(i, j, 0);
    in.rho_1 = rho3d(i, j, 0);
    in.dz8w_2 = dz8w(i, j, 1);
    in.u_2 = u3d(i, j, 1);
    in.v_2 = v3d(i, j, 1);
    in.u_1 = u3d(i, j, 0);
    in.v_1 = v3d(i, j, 0);
    in.mavail = mavail(i, j);
    in.pblh = pblh(i, j);
    in.psfcpa = psfcpa(i, j);
    in.snowh = snowh(i, j);
    in.dx = dx(i, j);
    in.tskin = tsk(i, j);
    in.tsurf = in.tskin;
    SfcPoint o;
    o.ust = ust(i, j);
    o.stress = o.ust * o.ust;
    o.mol = mol(i, j);
    o.qsfc = qsfc(i, j) / (one + qvmr);
    o.qgh = qgh(i, j) / (one + qvmr);
    o.qstar = zero;  // read by no default option before it is set
    o.hfx = hfx(i, j);
    o.qfx = qfx(i, j) / (one + qvmr);
    o.lh = lh(i, j);
    o.rmol = rmol(i, j);
    o.znt = znt(i, j);
    o.zol = zol(i, j);
    o.cpm = cpm(i, j);
    o.chs2 = chs2(i, j);
    o.cqs2 = cqs2(i, j);
    o.cqs = cqs(i, j);
    o.chs = chs(i, j);
    o.ch = ch(i, j);
    o.flhc = flhc(i, j);
    o.flqc = flqc(i, j);
    o.gz1oz0 = gz1oz0(i, j);
    o.wspd = wspd(i, j);
    o.rb = br(i, j);
    o.psim = psim(i, j);
    o.psih = psih(i, j);
    o.u10 = zero; o.v10 = zero; o.th2 = zero; o.t2 = zero; o.q2 = zero; o.wstar = zero;
    o.ustm = ustm(i, j);
    o.cm = zero; o.psix = zero; o.psit = zero; o.psix10 = zero; o.psit2 = zero;
    o.ck = zero; o.cka = zero; o.cd = zero; o.cda = zero;
    const real xland_1 = xland(i, j), xice_1 = xice(i, j);
    if (((xland_1 - real(1.5)) < zero) && (in.snowh < snow_thresh)) mynnsfc_land(in, o, *tb);
    if ((xland_1 - real(1.5)) >= zero) mynnsfc_water(in, lakemask(i, j), o, *tb);
    if ((((xland_1 - real(1.5)) < zero) && (in.snowh >= snow_thresh)) || (xice_1 >= xice_threshold && xice_1 < one))
      mynnsfc_ice(in, o, *tb);
    hfx(i, j) = o.hfx;
    qfx(i, j) = o.qfx / (one - in.qv_1);
    lh(i, j) = o.lh / (one - in.qv_1);
    mol(i, j) = o.mol;
    rmol(i, j) = o.rmol;
    qsfc(i, j) = o.qsfc / (one - in.qv_1);
    qgh(i, j) = o.qgh / (one - in.qv_1);
    znt(i, j) = o.znt;
    zol(i, j) = o.zol;
    ust(i, j) = o.ust;
    cpm(i, j) = o.cpm;
    chs2(i, j) = o.chs2;
    cqs2(i, j) = o.cqs2;
    cqs(i, j) = o.cqs;
    chs(i, j) = o.chs;
    ch(i, j) = o.ch;
    flhc(i, j) = o.flhc;
    flqc(i, j) = o.flqc;
    gz1oz0(i, j) = o.gz1oz0;
    wspd(i, j) = o.wspd;
    br(i, j) = o.rb;
    psim(i, j) = o.psim;
    psih(i, j) = o.psih;
    u10(i, j) = o.u10;
    v10(i, j) = o.v10;
    th2(i, j) = o.th2;
    t2(i, j) = o.t2;
    q2(i, j) = o.q2 / (one - in.qv_1);
    ck(i, j) = o.ck;
    cka(i, j) = o.cka;
    cd(i, j) = o.cd;
    cda(i, j) = o.cda;
    ustm(i, j) = o.ustm;
  }
};

}  // namespace cppwrf::phys::mynnsfc
