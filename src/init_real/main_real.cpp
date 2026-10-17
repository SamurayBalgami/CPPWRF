// cppwrf_real: WRF real.exe in C++. Reads namelist.input and met_em.d01.<date>.nc, writes wrfinput_d01
// (and wrfbdy_d01 for more than one time).
//   cppwrf_real [--namelist namelist.input] [--met-dir dir] [--out-dir dir] [--precision double|float]
#include <cstdio>
#include <cstring>
#include <memory>
#include <iostream>
#include <string>

#include "init_real/met_em.hpp"
#include "init_real/namelist.hpp"
#include "init_real/real_domain.hpp"
#include "init_real/real_output.hpp"
#include "init_real/wrf_schema.hpp"

using namespace cppwrf::init_real;

template <typename T>
static int run(const RealOptions& o, const Namelist& nl, const WrfSchema& sc, const std::string& met_dir,
               const std::string& out_dir) {
  RealDomain<T> dom(o);
  WrfConfig cfg(nl, sc.defaults(), o.dom);
  char dn[8];
  std::snprintf(dn, sizeof dn, "d%02d", o.dom);
  // setup_physics_suite (physics_suite = 'none'): unset physics options are 0
  for (const char* p : {"mp_physics", "cu_physics", "ra_lw_physics", "ra_sw_physics", "sf_sfclay_physics",
                        "sf_surface_physics", "bl_pbl_physics"})
    if (cfg.integer(p) < 0) cfg.set(p, "0");
  // set_physics_rconfigs / check_nml_consistency values used by the package conditions
  cfg.set("cu_used", cfg.integer("cu_physics") != 0 ? "1" : "0");
  cfg.set("shcu_used", cfg.integer("shcu_physics") != 0 ? "1" : "0");
  cfg.set("gwd_used", cfg.integer("gwd_opt") == 1 ? "1" : cfg.integer("gwd_opt") == 3 ? "3" : "0");
  cfg.set("aercu_used", cfg.integer("aercu_opt") > 0 ? "1" : "0");
  cfg.set("fft_used", "0");
  cfg.set("cam_used", "0");
  cfg.set("use_wps_input", "1");
  const Dims& d = dom.d;
  // compute_si_start_and_end
  int nsteps = 0;
  for (std::string t = o.start_date; t <= o.end_date; t = add_seconds(t, o.interval_seconds)) ++nsteps;
  // nests: only the initial time (no grid FDDA, sst_update or qna_update)
  const int time_loop_max = o.dom > 1 ? 1 : std::max(1, nsteps);
  const int w = o.spec_bdy_width;
  const T dt_bdy = T(o.interval_seconds);
  F3<T> ub1(1, d.ide, 1, d.kde, 1, d.jde), vb1 = ub1, tb1 = ub1, pb1 = ub1, qb1 = ub1, mb1(1, d.ide, 1, 1, 1, d.jde);
  F3<T> ub2 = ub1, vb2 = ub1, tb2 = ub1, pb2 = ub1, qb2 = ub1, mb2 = mb1;
  BdyVar<T> bu, bv, bt, bp, bq, bm, tu, tv, tt, tp, tq, tm;
  for (auto [b, k1] : {std::pair<BdyVar<T>*, int>{&bu, d.kde}, {&bv, d.kde}, {&bt, d.kde}, {&bp, d.kde}, {&bq, d.kde},
                       {&tu, d.kde}, {&tv, d.kde}, {&tt, d.kde}, {&tp, d.kde}, {&tq, d.kde}})
    b->alloc(d, 1, k1, w);
  bm.alloc(d, 1, 1, w);
  tm.alloc(d, 1, 1, w);
  const Dims d1{1, d.ide, 1, d.jde, 1, 1};
  auto do_couple = [&](F3<T>& u, F3<T>& v, F3<T>& t, F3<T>& p, F3<T>& q) {
    couple(dom.mu_2, dom.mub, u, dom.u_2, 'u', dom.msfuy, dom.c1h, dom.c2h, dom.c1f, dom.c2f, d);
    couple(dom.mu_2, dom.mub, v, dom.v_2, 'v', dom.msfvx, dom.c1h, dom.c2h, dom.c1f, dom.c2f, d);
    couple(dom.mu_2, dom.mub, t, dom.t_2, 't', dom.msfty, dom.c1h, dom.c2h, dom.c1h, dom.c2h, d);
    couple(dom.mu_2, dom.mub, p, dom.ph_2, 'h', dom.msfty, dom.c1h, dom.c2h, dom.c1f, dom.c2f, d);
    couple(dom.mu_2, dom.mub, q, dom.qv, 't', dom.msfty, dom.c1h, dom.c2h, dom.c1h, dom.c2h, d);
  };
  auto do_stuff = [&]() {
    stuff(bu, 'U', w, 1, d.kde, d, [&](int i, int k, int j) { return ub1(i, k, j); });
    stuff(bv, 'V', w, 1, d.kde, d, [&](int i, int k, int j) { return vb1(i, k, j); });
    stuff(bt, 'T', w, 1, d.kde, d, [&](int i, int k, int j) { return tb1(i, k, j); });
    stuff(bp, 'W', w, 1, d.kde, d, [&](int i, int k, int j) { return pb1(i, k, j); });
    stuff(bq, 'T', w, 1, d.kde, d, [&](int i, int k, int j) { return qb1(i, k, j); });
    stuff(bm, 'M', w, 1, 1, d, [&](int i, int k, int j) { return mb1(i, k, j); });
  };
  std::unique_ptr<WrfBdyWriter<T>> bdy;
  std::string date = o.start_date, prev = o.start_date;
  for (int loop = 1; loop <= time_loop_max; ++loop) {
    const MetEm m = MetEm::read(met_dir + "/met_em." + dn + "." + date + ".nc");
    dom.read_input(m);
    if (loop == 1) cfg.set("num_land_cat", std::to_string(m.att_int("NUM_LAND_CAT", 21)));
    dom.init_domain(loop, date);
    if (loop == 1) {
      write_wrfinput(dom, m, sc, cfg, out_dir + "/wrfinput_" + dn, date);
      std::cout << "cppwrf_real: wrote " << out_dir << "/wrfinput_" << dn << "\n";
      if (time_loop_max > 1) {
        do_couple(ub1, vb1, tb1, pb1, qb1);
        for (int j = 1; j <= d.jde - 1; ++j)
          for (int i = 1; i <= d.ide - 1; ++i) mb1(i, 1, j) = dom.mu_2(i, j);
        do_stuff();
      }
    } else {
      if (loop == 2) bdy = std::make_unique<WrfBdyWriter<T>>(out_dir + "/wrfbdy_d01", dom, m, sc, cfg, o.start_date, w);
      do_couple(ub2, vb2, tb2, pb2, qb2);
      for (int j = 1; j <= d.jde; ++j)
        for (int i = 1; i <= d.ide; ++i) mb2(i, 1, j) = dom.mu_2(i, j);
      auto tend = [&](BdyVar<T>& b, char st, int k1, const F3<T>& n2, const F3<T>& n1) {
        stuff(b, st, w, 1, k1, d, [&](int i, int k, int j) { return (n2(i, k, j) - n1(i, k, j)) / dt_bdy; });
      };
      tend(tu, 'U', d.kde, ub2, ub1); tend(tv, 'V', d.kde, vb2, vb1); tend(tt, 'T', d.kde, tb2, tb1);
      tend(tp, 'W', d.kde, pb2, pb1); tend(tq, 'T', d.kde, qb2, qb1); tend(tm, 'M', 1, mb2, mb1);
      (void)d1;
      bdy->write(std::size_t(loop - 2), prev, date,
                 {{"U", {&bu, &tu}}, {"V", {&bv, &tv}}, {"T", {&bt, &tt}}, {"PH", {&bp, &tp}}, {"QVAPOR", {&bq, &tq}},
                  {"MU", {&bm, &tm}}});
      if (loop < time_loop_max) {
        ub1 = ub2; vb1 = vb2; tb1 = tb2; pb1 = pb2; qb1 = qb2; mb1 = mb2;
        do_stuff();
      }
    }
    prev = date;
    date = add_seconds(date, o.interval_seconds);
  }
  if (bdy) std::cout << "cppwrf_real: wrote " << out_dir << "/wrfbdy_d01 (" << time_loop_max - 1 << " times)\n";
  return 0;
}

int main(int argc, char** argv) {
  std::string nml = "namelist.input", met_dir = ".", out_dir = ".", prec = "double";
  std::string schema = std::string(CPPWRF_DATA_DIR) + "/wrf_schema/wrf-4.8.0.tsv";
  for (int a = 1; a < argc; ++a) {
    std::string s = argv[a];
    auto next = [&]() -> std::string {
      if (a + 1 >= argc) throw std::runtime_error("missing value for " + s);
      return argv[++a];
    };
    if (s == "--namelist") nml = next();
    else if (s == "--met-dir") met_dir = next();
    else if (s == "--out-dir") out_dir = next();
    else if (s == "--precision") prec = next();
    else if (s == "--schema") schema = next();
    else {
      std::cerr << "usage: cppwrf_real [--namelist f] [--met-dir d] [--out-dir d] [--precision double|float]\n";
      return 2;
    }
  }
  try {
    const Namelist nl = Namelist::read(nml);
    const WrfSchema sc = WrfSchema::load(schema);
    const int max_dom = nl.get_int("domains", "max_dom", 1);
    for (int id = 1; id <= max_dom; ++id) {
      const RealOptions o = RealOptions::from_namelist(nl, id);
      const int rc = prec == "float" ? run<float>(o, nl, sc, met_dir, out_dir) : run<double>(o, nl, sc, met_dir, out_dir);
      if (rc != 0) return rc;
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "cppwrf_real: " << e.what() << "\n";
    return 1;
  }
}
