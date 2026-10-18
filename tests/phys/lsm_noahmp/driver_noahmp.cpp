#include "phys/lsm_noahmp/noahmp.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using cppwrf::real;
using cppwrf::View2;
using cppwrf::View3;

namespace {

void read_exact(std::ifstream& in, void* buf, std::size_t bytes) {
  in.read(reinterpret_cast<char*>(buf), static_cast<std::streamsize>(bytes));
  if (!in) {
    std::cerr << "Error reading input binary file\n";
    std::exit(1);
  }
}

std::vector<real> to_real_vec(const std::vector<double>& src) {
  std::vector<real> dst(src.size());
  for (std::size_t i = 0; i < src.size(); ++i) {
    dst[i] = static_cast<real>(src[i]);
  }
  return dst;
}

void write_vec2d(std::ofstream& out, const std::vector<real>& src, int n_cols) {
  std::vector<double> tmp(n_cols);
  for (int i = 0; i < n_cols; ++i) {
    tmp[i] = static_cast<double>(src[i]);
  }
  out.write(reinterpret_cast<const char*>(tmp.data()), static_cast<std::streamsize>(n_cols * sizeof(double)));
}

void write_vec3d(std::ofstream& out, const std::vector<real>& src, int n_cols, int nz) {
  std::vector<double> tmp(static_cast<std::size_t>(n_cols) * nz);
  for (int i = 0; i < n_cols * nz; ++i) {
    tmp[i] = static_cast<double>(src[i]);
  }
  out.write(reinterpret_cast<const char*>(tmp.data()), static_cast<std::streamsize>(tmp.size() * sizeof(double)));
}

}  // namespace

int main(int argc, char** argv) {
  std::string infile = (argc >= 2) ? argv[1] : "noahmp_in.bin";
  std::string outfile = (argc >= 3) ? argv[2] : "noahmp_cpp.bin";
  std::string tbl_path = (argc >= 4) ? argv[3] : "NoahmpTable.TBL";

  std::ifstream in(infile, std::ios::binary);
  if (!in) {
    std::cerr << "Cannot open input file: " << infile << "\n";
    return 1;
  }

  int n_cols = 0, nsoil = 0, n_steps = 0, yr = 0;
  double dt_d = 0.0, dx_d = 0.0, julian_start_d = 0.0, xice_threshold_d = 0.0;
  read_exact(in, &n_cols, sizeof(int));
  read_exact(in, &nsoil, sizeof(int));
  read_exact(in, &n_steps, sizeof(int));
  read_exact(in, &yr, sizeof(int));
  read_exact(in, &dt_d, sizeof(double));
  read_exact(in, &dx_d, sizeof(double));
  read_exact(in, &julian_start_d, sizeof(double));
  read_exact(in, &xice_threshold_d, sizeof(double));

  std::vector<double> dzs_d(nsoil);
  read_exact(in, dzs_d.data(), nsoil * sizeof(double));

  std::vector<int> ivgtyp(n_cols), isltyp(n_cols);
  read_exact(in, ivgtyp.data(), n_cols * sizeof(int));
  read_exact(in, isltyp.data(), n_cols * sizeof(int));

  auto read_d2 = [&]() {
    std::vector<double> v(n_cols);
    read_exact(in, v.data(), n_cols * sizeof(double));
    return to_real_vec(v);
  };
  auto read_d3 = [&](int nz) {
    std::vector<double> v(static_cast<std::size_t>(n_cols) * nz);
    read_exact(in, v.data(), v.size() * sizeof(double));
    return to_real_vec(v);
  };

  auto xlat = read_d2();
  auto xlong = read_d2();
  auto coszen = read_d2();
  auto vegfra = read_d2();
  auto shdmax = read_d2();
  auto tmn = read_d2();
  auto xland = read_d2();
  auto xice = read_d2();
  auto swdown = read_d2();
  auto glw = read_d2();
  auto rainbl = read_d2();
  auto sr = read_d2();
  auto mp_rainc = read_d2();
  auto mp_rainnc = read_d2();
  auto mp_shcv = read_d2();
  auto mp_snow = read_d2();
  auto mp_graup = read_d2();
  auto mp_hail = read_d2();
  auto dz8w_1 = read_d2();
  auto t_phy_1 = read_d2();
  auto qv_curr_1 = read_d2();
  auto u_phy_1 = read_d2();
  auto v_phy_1 = read_d2();
  auto p8w_1 = read_d2();
  auto p8w_2 = read_d2();
  auto tsk = read_d2();
  auto snow = read_d2();
  auto snowh = read_d2();
  auto canwat = read_d2();
  auto lai = read_d2();

  auto tslb = read_d3(nsoil);
  auto smois = read_d3(nsoil);
  auto sh2o = read_d3(nsoil);
  in.close();

  const int nsnow = 3;
  std::vector<real> p8w(static_cast<std::size_t>(n_cols) * 2, real(0.0));
  for (int i = 0; i < n_cols; ++i) {
    p8w[i] = p8w_1[i];
    p8w[n_cols + i] = p8w_2[i];
  }

  // Allocate state and output vectors matching driver_noahmp.f90
  std::vector<int> isnowxy(n_cols, 0), cropcat(n_cols, 0), pgsxy(n_cols, 0);
  std::vector<real> hfx(n_cols, real(0.0)), qfx(n_cols, real(0.0)), lh(n_cols, real(0.0)), grdflx(n_cols, real(0.0));
  std::vector<real> smstav(n_cols, real(0.0)), smstot(n_cols, real(0.0));
  std::vector<real> sfcrunoff(n_cols, real(0.0)), udrunoff(n_cols, real(0.0));
  std::vector<real> albedo(n_cols, real(0.18)), snowc(n_cols, real(0.0));
  std::vector<real> acsnom(n_cols, real(0.0)), acsnow(n_cols, real(0.0));
  std::vector<real> emiss(n_cols, real(0.98)), qsfc(n_cols, real(0.005));
  std::vector<real> z0(n_cols, real(0.05)), znt(n_cols, real(0.05)), rs(n_cols, real(0.0));
  std::vector<real> tvxy(n_cols, real(0.0)), tgxy(n_cols, real(0.0));
  std::vector<real> canicexy(n_cols, real(0.0)), canliqxy(n_cols, real(0.0));
  std::vector<real> eahxy(n_cols, real(0.0)), tahxy(n_cols, real(0.0));
  std::vector<real> cmxy(n_cols, real(0.0)), chxy(n_cols, real(0.0));
  std::vector<real> fwetxy(n_cols, real(0.0)), sneqvoxy(n_cols, real(0.0)), alboldxy(n_cols, real(0.0));
  std::vector<real> qsnowxy(n_cols, real(0.0)), qrainxy(n_cols, real(0.0));
  std::vector<real> wslakexy(n_cols, real(0.0)), zwtxy(n_cols, real(0.0));
  std::vector<real> waxy(n_cols, real(0.0)), wtxy(n_cols, real(0.0));
  std::vector<real> lfmassxy(n_cols, real(0.0)), rtmassxy(n_cols, real(0.0));
  std::vector<real> stmassxy(n_cols, real(0.0)), woodxy(n_cols, real(0.0));
  std::vector<real> stblcpxy(n_cols, real(0.0)), fastcpxy(n_cols, real(0.0));
  std::vector<real> xsaixy(n_cols, real(0.0)), taussxy(n_cols, real(0.0));
  std::vector<real> smcwtdxy(n_cols, real(0.0)), deeprechxy(n_cols, real(0.0)), rechxy(n_cols, real(0.0));
  std::vector<real> grainxy(n_cols, real(0.0)), gddxy(n_cols, real(0.0)), qtdrain(n_cols, real(0.0));
  std::vector<real> acc_ssoilxy(n_cols, real(0.0)), acc_qinsurxy(n_cols, real(0.0)), acc_qsevaxy(n_cols, real(0.0));
  std::vector<real> acc_dwaterxy(n_cols, real(0.0)), acc_prcpxy(n_cols, real(0.0));
  std::vector<real> acc_ecanxy(n_cols, real(0.0)), acc_etranxy(n_cols, real(0.0)), acc_edirxy(n_cols, real(0.0));
  std::vector<real> chstarxy(n_cols, real(0.0));

  // 2D OUT diagnostics
  std::vector<real> t2mvxy(n_cols, real(0.0)), t2mbxy(n_cols, real(0.0));
  std::vector<real> q2mvxy(n_cols, real(0.0)), q2mbxy(n_cols, real(0.0));
  std::vector<real> tradxy(n_cols, real(0.0)), neexy(n_cols, real(0.0)), gppxy(n_cols, real(0.0)), nppxy(n_cols, real(0.0));
  std::vector<real> fvegxy(n_cols, real(0.0)), runsfxy(n_cols, real(0.0)), runsbxy(n_cols, real(0.0));
  std::vector<real> ecanxy(n_cols, real(0.0)), edirxy(n_cols, real(0.0)), etranxy(n_cols, real(0.0));
  std::vector<real> fsaxy(n_cols, real(0.0)), firaxy(n_cols, real(0.0)), aparxy(n_cols, real(0.0)), psnxy(n_cols, real(0.0));
  std::vector<real> savxy(n_cols, real(0.0)), sagxy(n_cols, real(0.0)), rssunxy(n_cols, real(0.0)), rsshaxy(n_cols, real(0.0));
  std::vector<real> bgapxy(n_cols, real(0.0)), wgapxy(n_cols, real(0.0)), tgvxy(n_cols, real(0.0)), tgbxy(n_cols, real(0.0));
  std::vector<real> chvxy(n_cols, real(0.0)), chbxy(n_cols, real(0.0));
  std::vector<real> shgxy(n_cols, real(0.0)), shcxy(n_cols, real(0.0)), shbxy(n_cols, real(0.0));
  std::vector<real> evgxy(n_cols, real(0.0)), evbxy(n_cols, real(0.0));
  std::vector<real> ghvxy(n_cols, real(0.0)), ghbxy(n_cols, real(0.0));
  std::vector<real> irgxy(n_cols, real(0.0)), ircxy(n_cols, real(0.0)), irbxy(n_cols, real(0.0));
  std::vector<real> trxy(n_cols, real(0.0)), evcxy(n_cols, real(0.0));
  std::vector<real> chleafxy(n_cols, real(0.0)), chucxy(n_cols, real(0.0)), chv2xy(n_cols, real(0.0)), chb2xy(n_cols, real(0.0));
  std::vector<real> pahxy(n_cols, real(0.0)), pahgxy(n_cols, real(0.0)), pahvxy(n_cols, real(0.0)), pahbxy(n_cols, real(0.0));
  std::vector<real> qintsxy(n_cols, real(0.0)), qintrxy(n_cols, real(0.0)), qdripsxy(n_cols, real(0.0)), qdriprxy(n_cols, real(0.0));
  std::vector<real> qthrosxy(n_cols, real(0.0)), qthrorxy(n_cols, real(0.0)), qsnsubxy(n_cols, real(0.0)), qsnfroxy(n_cols, real(0.0));
  std::vector<real> qsubcxy(n_cols, real(0.0)), qfrocxy(n_cols, real(0.0)), qevacxy(n_cols, real(0.0)), qdewcxy(n_cols, real(0.0));
  std::vector<real> qfrzcxy(n_cols, real(0.0)), qmeltcxy(n_cols, real(0.0)), qsnbotxy(n_cols, real(0.0)), qmeltxy(n_cols, real(0.0));
  std::vector<real> pondingxy(n_cols, real(0.0)), fpicexy(n_cols, real(0.0)), rainlsm(n_cols, real(0.0)), snowlsm(n_cols, real(0.0));
  std::vector<real> forctlsm(n_cols, real(0.0)), forcqlsm(n_cols, real(0.0)), forcplsm(n_cols, real(0.0));
  std::vector<real> forczlsm(n_cols, real(0.0)), forcwlsm(n_cols, real(0.0));
  std::vector<real> eflxbxy(n_cols, real(0.0)), soilenergy(n_cols, real(0.0)), snowenergy(n_cols, real(0.0)), canhsxy(n_cols, real(0.0));

  // 3D INOUT state
  std::vector<real> smoiseq(static_cast<std::size_t>(n_cols) * nsoil, real(0.0));
  std::vector<real> acc_etranixy(static_cast<std::size_t>(n_cols) * nsoil, real(0.0));
  std::vector<real> tsnoxy(static_cast<std::size_t>(n_cols) * nsnow, real(0.0));
  std::vector<real> snicexy(static_cast<std::size_t>(n_cols) * nsnow, real(0.0));
  std::vector<real> snliqxy(static_cast<std::size_t>(n_cols) * nsnow, real(0.0));
  std::vector<real> zsnsoxy(static_cast<std::size_t>(n_cols) * (nsnow + nsoil), real(0.0));
  std::vector<real> albsoildirxy(static_cast<std::size_t>(n_cols) * 2, real(0.0));
  std::vector<real> albsoildifxy(static_cast<std::size_t>(n_cols) * 2, real(0.0));

  const std::string mminlu = "MODIFIED_IGBP_MODIS_NOAH";
  cppwrf::phys::noahmp::NoahmpOptions opts{};
  cppwrf::phys::noahmp::validate_noahmp_options(opts);
  cppwrf::phys::noahmp::NoahmpTables tables{};
  if (!cppwrf::phys::noahmp::read_noahmp_tables(tbl_path, mminlu, tables)) {
    std::cerr << "Failed to read NoahmpTable.TBL from " << tbl_path << "\n";
    return 1;
  }

  auto v2i = [&](std::vector<int>& v) { return View2<int>{v.data(), n_cols}; };
  auto v2ci = [&](const std::vector<int>& v) { return View2<const int>{v.data(), n_cols}; };
  auto v2r = [&](std::vector<real>& v) { return View2<real>{v.data(), n_cols}; };
  auto v2cr = [&](const std::vector<real>& v) { return View2<const real>{v.data(), n_cols}; };
  auto v3r = [&](std::vector<real>& v, int) { return View3<real>{v.data(), n_cols, n_cols}; };
  auto v3cr = [&](const std::vector<real>& v, int) { return View3<const real>{v.data(), n_cols, n_cols}; };
  const cppwrf::Range2D range{0, n_cols, 0, 1};

  // 1. Run Noah-MP initialization
  cppwrf::phys::noahmp::NoahmpInitKernel init_k{};
  init_k.ivgtyp = v2ci(ivgtyp);
  init_k.isltyp = v2ci(isltyp);
  init_k.tsk = v2cr(tsk);
  init_k.xice = v2cr(xice);
  init_k.tmn = v2r(tmn);
  init_k.isnowxy = v2i(isnowxy);
  init_k.cropcat = v2i(cropcat);
  init_k.smois = v3r(smois, nsoil);
  init_k.sh2o = v3r(sh2o, nsoil);
  init_k.tslb = v3r(tslb, nsoil);
  init_k.snow = v2r(snow);
  init_k.snowh = v2r(snowh);
  init_k.canwat = v2r(canwat);
  init_k.canicexy = v2r(canicexy);
  init_k.canliqxy = v2r(canliqxy);
  init_k.zsnsoxy = v3r(zsnsoxy, nsnow + nsoil);
  init_k.tsnoxy = v3r(tsnoxy, nsnow);
  init_k.snicexy = v3r(snicexy, nsnow);
  init_k.snliqxy = v3r(snliqxy, nsnow);
  init_k.tvxy = v2r(tvxy);
  init_k.tgxy = v2r(tgxy);
  init_k.eahxy = v2r(eahxy);
  init_k.tahxy = v2r(tahxy);
  init_k.cmxy = v2r(cmxy);
  init_k.chxy = v2r(chxy);
  init_k.fwetxy = v2r(fwetxy);
  init_k.sneqvoxy = v2r(sneqvoxy);
  init_k.alboldxy = v2r(alboldxy);
  init_k.qsnowxy = v2r(qsnowxy);
  init_k.qrainxy = v2r(qrainxy);
  init_k.wslakexy = v2r(wslakexy);
  init_k.zwtxy = v2r(zwtxy);
  init_k.waxy = v2r(waxy);
  init_k.wtxy = v2r(wtxy);
  init_k.lfmassxy = v2r(lfmassxy);
  init_k.rtmassxy = v2r(rtmassxy);
  init_k.stmassxy = v2r(stmassxy);
  init_k.woodxy = v2r(woodxy);
  init_k.stblcpxy = v2r(stblcpxy);
  init_k.fastcpxy = v2r(fastcpxy);
  init_k.lai = v2r(lai);
  init_k.xsaixy = v2r(xsaixy);
  init_k.t2mvxy = v2r(t2mvxy);
  init_k.t2mbxy = v2r(t2mbxy);
  init_k.grainxy = v2r(grainxy);
  init_k.gddxy = v2r(gddxy);
  init_k.qtdrain = v2r(qtdrain);
  init_k.chstarxy = v2r(chstarxy);
  init_k.albsoildirxy = v3r(albsoildirxy, 2);
  init_k.albsoildifxy = v3r(albsoildifxy, 2);
  init_k.nsoil = nsoil;
  init_k.nsnow = nsnow;
  init_k.fndsnowh = true;
  init_k.restart = false;
  init_k.iopt_runsub = opts.opt_runsub;
  init_k.iopt_crop = opts.opt_crop;
  init_k.sf_urban_physics = opts.sf_urban_physics;
  for (int k = 0; k < nsoil; ++k) {
    init_k.dzs[k] = static_cast<real>(dzs_d[k]);
  }
  init_k.tables = &tables;

  cppwrf::phys::noahmp::run_noahmp_init(init_k, range);

  // 2. Run Noah-MP per-step kernel for n_steps
  cppwrf::phys::noahmp::NoahmpKernel step_k{};
  step_k.ivgtyp = v2ci(ivgtyp);
  step_k.isltyp = v2ci(isltyp);
  step_k.cropcat = v2ci(cropcat);
  step_k.coszen = v2cr(coszen);
  step_k.xlat = v2cr(xlat);
  step_k.xlong = v2cr(xlong);
  step_k.dz8w_0 = v2cr(dz8w_1);
  step_k.vegfra = v2cr(vegfra);
  step_k.shdmax = v2cr(shdmax);
  step_k.tmn = v2cr(tmn);
  step_k.xland = v2cr(xland);
  step_k.xice = v2cr(xice);
  step_k.t_phy_0 = v2cr(t_phy_1);
  step_k.qv_curr_0 = v2cr(qv_curr_1);
  step_k.u_phy_0 = v2cr(u_phy_1);
  step_k.v_phy_0 = v2cr(v_phy_1);
  step_k.swdown = v2cr(swdown);
  step_k.glw = v2cr(glw);
  step_k.p8w = v3cr(p8w, 2);
  step_k.rainbl = v2cr(rainbl);
  step_k.sr = v2cr(sr);
  step_k.raincv = v2cr(mp_rainc);
  step_k.rainncv = v2cr(mp_rainnc);
  step_k.rainshv = v2cr(mp_shcv);
  step_k.snowncv = v2cr(mp_snow);
  step_k.graupelncv = v2cr(mp_graup);
  step_k.hailncv = v2cr(mp_hail);

  step_k.tsk = v2r(tsk);
  step_k.hfx = v2r(hfx);
  step_k.qfx = v2r(qfx);
  step_k.lh = v2r(lh);
  step_k.grdflx = v2r(grdflx);
  step_k.smstav = v2r(smstav);
  step_k.smstot = v2r(smstot);
  step_k.sfcrunoff = v2r(sfcrunoff);
  step_k.udrunoff = v2r(udrunoff);
  step_k.albedo = v2r(albedo);
  step_k.snowc = v2r(snowc);
  step_k.smois = v3r(smois, nsoil);
  step_k.sh2o = v3r(sh2o, nsoil);
  step_k.tslb = v3r(tslb, nsoil);
  step_k.snow = v2r(snow);
  step_k.snowh = v2r(snowh);
  step_k.canwat = v2r(canwat);
  step_k.acsnom = v2r(acsnom);
  step_k.acsnow = v2r(acsnow);
  step_k.emiss = v2r(emiss);
  step_k.qsfc = v2r(qsfc);
  step_k.z0 = v2r(z0);
  step_k.znt = v2r(znt);
  step_k.rs = v2r(rs);

  step_k.isnowxy = v2i(isnowxy);
  step_k.tvxy = v2r(tvxy);
  step_k.tgxy = v2r(tgxy);
  step_k.canicexy = v2r(canicexy);
  step_k.canliqxy = v2r(canliqxy);
  step_k.eahxy = v2r(eahxy);
  step_k.tahxy = v2r(tahxy);
  step_k.cmxy = v2r(cmxy);
  step_k.chxy = v2r(chxy);
  step_k.fwetxy = v2r(fwetxy);
  step_k.sneqvoxy = v2r(sneqvoxy);
  step_k.alboldxy = v2r(alboldxy);
  step_k.qsnowxy = v2r(qsnowxy);
  step_k.qrainxy = v2r(qrainxy);
  step_k.wslakexy = v2r(wslakexy);
  step_k.zwtxy = v2r(zwtxy);
  step_k.waxy = v2r(waxy);
  step_k.wtxy = v2r(wtxy);
  step_k.tsnoxy = v3r(tsnoxy, nsnow);
  step_k.zsnsoxy = v3r(zsnsoxy, nsnow + nsoil);
  step_k.snicexy = v3r(snicexy, nsnow);
  step_k.snliqxy = v3r(snliqxy, nsnow);
  step_k.lfmassxy = v2r(lfmassxy);
  step_k.rtmassxy = v2r(rtmassxy);
  step_k.stmassxy = v2r(stmassxy);
  step_k.woodxy = v2r(woodxy);
  step_k.stblcpxy = v2r(stblcpxy);
  step_k.fastcpxy = v2r(fastcpxy);
  step_k.lai = v2r(lai);
  step_k.xsaixy = v2r(xsaixy);
  step_k.taussxy = v2r(taussxy);
  step_k.smoiseq = v3r(smoiseq, nsoil);
  step_k.smcwtdxy = v2r(smcwtdxy);
  step_k.deeprechxy = v2r(deeprechxy);
  step_k.rechxy = v2r(rechxy);
  step_k.grainxy = v2r(grainxy);
  step_k.gddxy = v2r(gddxy);
  step_k.pgsxy = v2i(pgsxy);
  step_k.qtdrain = v2r(qtdrain);
  step_k.acc_ssoilxy = v2r(acc_ssoilxy);
  step_k.acc_qinsurxy = v2r(acc_qinsurxy);
  step_k.acc_qsevaxy = v2r(acc_qsevaxy);
  step_k.acc_etranixy = v3r(acc_etranixy, nsoil);
  step_k.acc_dwaterxy = v2r(acc_dwaterxy);
  step_k.acc_prcpxy = v2r(acc_prcpxy);
  step_k.acc_ecanxy = v2r(acc_ecanxy);
  step_k.acc_etranxy = v2r(acc_etranxy);
  step_k.acc_edirxy = v2r(acc_edirxy);
  step_k.albsoildirxy = v3r(albsoildirxy, 2);
  step_k.albsoildifxy = v3r(albsoildifxy, 2);

  step_k.t2mvxy = v2r(t2mvxy);
  step_k.t2mbxy = v2r(t2mbxy);
  step_k.q2mvxy = v2r(q2mvxy);
  step_k.q2mbxy = v2r(q2mbxy);
  step_k.fvegxy = v2r(fvegxy);
  step_k.tradxy = v2r(tradxy);
  step_k.neexy = v2r(neexy);
  step_k.gppxy = v2r(gppxy);
  step_k.nppxy = v2r(nppxy);
  step_k.runsfxy = v2r(runsfxy);
  step_k.runsbxy = v2r(runsbxy);
  step_k.ecanxy = v2r(ecanxy);
  step_k.edirxy = v2r(edirxy);
  step_k.etranxy = v2r(etranxy);
  step_k.fsaxy = v2r(fsaxy);
  step_k.firaxy = v2r(firaxy);
  step_k.aparxy = v2r(aparxy);
  step_k.psnxy = v2r(psnxy);
  step_k.savxy = v2r(savxy);
  step_k.sagxy = v2r(sagxy);
  step_k.rssunxy = v2r(rssunxy);
  step_k.rsshaxy = v2r(rsshaxy);
  step_k.bgapxy = v2r(bgapxy);
  step_k.wgapxy = v2r(wgapxy);
  step_k.tgvxy = v2r(tgvxy);
  step_k.tgbxy = v2r(tgbxy);
  step_k.chvxy = v2r(chvxy);
  step_k.chbxy = v2r(chbxy);
  step_k.shgxy = v2r(shgxy);
  step_k.shcxy = v2r(shcxy);
  step_k.shbxy = v2r(shbxy);
  step_k.evgxy = v2r(evgxy);
  step_k.evbxy = v2r(evbxy);
  step_k.ghvxy = v2r(ghvxy);
  step_k.ghbxy = v2r(ghbxy);
  step_k.irgxy = v2r(irgxy);
  step_k.ircxy = v2r(ircxy);
  step_k.irbxy = v2r(irbxy);
  step_k.trxy = v2r(trxy);
  step_k.evcxy = v2r(evcxy);
  step_k.chleafxy = v2r(chleafxy);
  step_k.chucxy = v2r(chucxy);
  step_k.chv2xy = v2r(chv2xy);
  step_k.chb2xy = v2r(chb2xy);
  step_k.qintsxy = v2r(qintsxy);
  step_k.qintrxy = v2r(qintrxy);
  step_k.qdripsxy = v2r(qdripsxy);
  step_k.qdriprxy = v2r(qdriprxy);
  step_k.qthrosxy = v2r(qthrosxy);
  step_k.qthrorxy = v2r(qthrorxy);
  step_k.qsnsubxy = v2r(qsnsubxy);
  step_k.qsnfroxy = v2r(qsnfroxy);
  step_k.qsubcxy = v2r(qsubcxy);
  step_k.qfrocxy = v2r(qfrocxy);
  step_k.qevacxy = v2r(qevacxy);
  step_k.qdewcxy = v2r(qdewcxy);
  step_k.qfrzcxy = v2r(qfrzcxy);
  step_k.qmeltcxy = v2r(qmeltcxy);
  step_k.qsnbotxy = v2r(qsnbotxy);
  step_k.qmeltxy = v2r(qmeltxy);
  step_k.pondingxy = v2r(pondingxy);
  step_k.pahxy = v2r(pahxy);
  step_k.pahgxy = v2r(pahgxy);
  step_k.pahvxy = v2r(pahvxy);
  step_k.pahbxy = v2r(pahbxy);
  step_k.fpicexy = v2r(fpicexy);
  step_k.rainlsm = v2r(rainlsm);
  step_k.snowlsm = v2r(snowlsm);
  step_k.forctlsm = v2r(forctlsm);
  step_k.forcqlsm = v2r(forcqlsm);
  step_k.forcplsm = v2r(forcplsm);
  step_k.forczlsm = v2r(forczlsm);
  step_k.forcwlsm = v2r(forcwlsm);
  step_k.eflxbxy = v2r(eflxbxy);
  step_k.soilenergy = v2r(soilenergy);
  step_k.snowenergy = v2r(snowenergy);
  step_k.canhsxy = v2r(canhsxy);

  step_k.yr = yr;
  step_k.nsoil = nsoil;
  step_k.nsnow = nsnow;
  step_k.dtbl = static_cast<real>(dt_d);
  step_k.dx = static_cast<real>(dx_d);
  step_k.xice_threshold = static_cast<real>(xice_threshold_d);
  for (int k = 0; k < nsoil; ++k) {
    step_k.dzs[k] = static_cast<real>(dzs_d[k]);
  }
  step_k.opts = opts;
  step_k.tables = &tables;

  for (int step = 1; step <= n_steps; ++step) {
    step_k.itimestep = step;
    step_k.julian_in = static_cast<real>(julian_start_d + static_cast<double>(step - 1) * dt_d / 86400.0);
    cppwrf::phys::noahmp::run_noahmp(step_k, range);
  }

  std::printf("C++ Col 1: tsk = %.15f hfx = %.15f lh = %.15f\n",
              static_cast<double>(tsk[0]), static_cast<double>(hfx[0]), static_cast<double>(lh[0]));

  std::ofstream out(outfile, std::ios::binary | std::ios::trunc);
  if (!out) {
    std::cerr << "Cannot open output file: " << outfile << "\n";
    return 1;
  }
  out.write(reinterpret_cast<const char*>(&n_cols), sizeof(int));
  out.write(reinterpret_cast<const char*>(&nsoil), sizeof(int));

  std::vector<real> isnowxy_r(n_cols);
  for (int i = 0; i < n_cols; ++i) {
    isnowxy_r[i] = static_cast<real>(isnowxy[i]);
  }

  // 62 2D INOUT state arrays
  const std::vector<real>* inout2d[] = {
      &tsk, &hfx, &qfx, &lh, &grdflx, &smstav, &smstot, &sfcrunoff, &udrunoff,
      &albedo, &snowc, &snow, &snowh, &canwat, &acsnom, &acsnow, &emiss, &qsfc,
      &z0, &znt, &rs, &isnowxy_r, &tvxy, &tgxy, &canicexy, &canliqxy, &eahxy,
      &tahxy, &cmxy, &chxy, &fwetxy, &sneqvoxy, &alboldxy, &qsnowxy, &qrainxy,
      &wslakexy, &zwtxy, &waxy, &wtxy, &lfmassxy, &rtmassxy, &stmassxy, &woodxy,
      &stblcpxy, &fastcpxy, &lai, &xsaixy, &taussxy, &smcwtdxy, &deeprechxy,
      &rechxy, &grainxy, &gddxy, &qtdrain, &acc_ssoilxy, &acc_qinsurxy,
      &acc_qsevaxy, &acc_dwaterxy, &acc_prcpxy, &acc_ecanxy, &acc_etranxy, &acc_edirxy};
  for (const auto* vec : inout2d) {
    write_vec2d(out, *vec, n_cols);
  }

  // 77 2D OUT diagnostic arrays
  const std::vector<real>* out2d[] = {
      &t2mvxy, &t2mbxy, &q2mvxy, &q2mbxy, &tradxy, &neexy, &gppxy, &nppxy,
      &fvegxy, &runsfxy, &runsbxy, &ecanxy, &edirxy, &etranxy, &fsaxy, &firaxy,
      &aparxy, &psnxy, &savxy, &sagxy, &rssunxy, &rsshaxy, &bgapxy, &wgapxy,
      &tgvxy, &tgbxy, &chvxy, &chbxy, &shgxy, &shcxy, &shbxy, &evgxy, &evbxy,
      &ghvxy, &ghbxy, &irgxy, &ircxy, &irbxy, &trxy, &evcxy, &chleafxy, &chucxy,
      &chv2xy, &chb2xy, &pahxy, &pahgxy, &pahvxy, &pahbxy, &qintsxy, &qintrxy,
      &qdripsxy, &qdriprxy, &qthrosxy, &qthrorxy, &qsnsubxy, &qsnfroxy,
      &qsubcxy, &qfrocxy, &qevacxy, &qdewcxy, &qfrzcxy, &qmeltcxy, &qsnbotxy,
      &qmeltxy, &pondingxy, &fpicexy, &rainlsm, &snowlsm, &forctlsm, &forcqlsm,
      &forcplsm, &forczlsm, &forcwlsm, &eflxbxy, &soilenergy, &snowenergy, &canhsxy};
  for (const auto* vec : out2d) {
    write_vec2d(out, *vec, n_cols);
  }

  // 11 3D arrays
  write_vec3d(out, smois, n_cols, nsoil);
  write_vec3d(out, sh2o, n_cols, nsoil);
  write_vec3d(out, tslb, n_cols, nsoil);
  write_vec3d(out, smoiseq, n_cols, nsoil);
  write_vec3d(out, acc_etranixy, n_cols, nsoil);
  write_vec3d(out, tsnoxy, n_cols, nsnow);
  write_vec3d(out, snicexy, n_cols, nsnow);
  write_vec3d(out, snliqxy, n_cols, nsnow);
  write_vec3d(out, zsnsoxy, n_cols, nsnow + nsoil);
  write_vec3d(out, albsoildirxy, n_cols, 2);
  write_vec3d(out, albsoildifxy, n_cols, 2);

  out.close();
  return 0;
}
