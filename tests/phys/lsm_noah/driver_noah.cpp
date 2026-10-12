#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "phys/lsm_noah/noah.hpp"
#include "phys/lsm_noah/noah_tables.hpp"

using namespace cppwrf;

template <typename T>
struct DevBuffer {
  T* p = nullptr;
  std::size_t n = 0;
  DevBuffer() = default;
  explicit DevBuffer(std::size_t count)
      : n(count), p(static_cast<T*>(dev_alloc(count * sizeof(T)))) {}
  ~DevBuffer() {
    if (p) dev_free(p);
  }
  DevBuffer(DevBuffer&& o) noexcept : p(o.p), n(o.n) {
    o.p = nullptr;
    o.n = 0;
  }
  DevBuffer& operator=(DevBuffer&& o) noexcept {
    if (this != &o) {
      if (p) dev_free(p);
      p = o.p;
      n = o.n;
      o.p = nullptr;
      o.n = 0;
    }
    return *this;
  }
  DevBuffer(const DevBuffer&) = delete;
  DevBuffer& operator=(const DevBuffer&) = delete;
  T* get() const { return p; }
};

int main(int argc, char** argv) {
  std::string infile = "tests/phys/lsm_noah/noah_input.bin";
  std::string outfile = "tests/phys/lsm_noah/noah_output_cpp.bin";
  std::string table_dir = "tests/phys/lsm_noah";
  bool bench_mode = false;
  int bench_iters = 10;

  for (int a = 1; a < argc; ++a) {
    std::string arg = argv[a];
    if (arg == "--bench") {
      bench_mode = true;
      if (a + 1 < argc && argv[a + 1][0] != '-') {
        bench_iters = std::atoi(argv[++a]);
      }
    } else if (a == 1 && arg[0] != '-') {
      infile = arg;
    } else if (a == 2 && arg[0] != '-') {
      outfile = arg;
    } else if (a == 3 && arg[0] != '-') {
      table_dir = arg;
    }
  }

  // Read parameter tables
  phys::noah::NoahTables tables;
  if (!phys::noah::read_noah_tables(table_dir, "USGS", "STAS", tables)) {
    std::fprintf(stderr, "Failed reading Noah tables from %s\n", table_dir.c_str());
    return 1;
  }

  // Open input file
  FILE* f = std::fopen(infile.c_str(), "rb");
  if (!f) {
    std::fprintf(stderr, "Cannot open input file: %s\n", infile.c_str());
    return 1;
  }

  int32_t n_cols = 0, nsoil = 0, isurban = 0;
  double dt_in = 0.0, zlvl_in = 0.0;
  if (std::fread(&n_cols, sizeof(int32_t), 1, f) != 1 ||
      std::fread(&nsoil, sizeof(int32_t), 1, f) != 1 ||
      std::fread(&dt_in, sizeof(double), 1, f) != 1 ||
      std::fread(&zlvl_in, sizeof(double), 1, f) != 1 ||
      std::fread(&isurban, sizeof(int32_t), 1, f) != 1) {
    std::fprintf(stderr, "Failed reading scalar inputs from %s\n", infile.c_str());
    std::fclose(f);
    return 1;
  }

  std::vector<double> h_sldpth(nsoil);
  if (std::fread(h_sldpth.data(), sizeof(double), nsoil, f) != std::size_t(nsoil)) {
    std::fprintf(stderr, "Failed reading sldpth from %s\n", infile.c_str());
    std::fclose(f);
    return 1;
  }

  std::size_t n2d = std::size_t(n_cols);
  std::size_t n3d = std::size_t(n_cols) * std::size_t(nsoil);

  // Read 2D integer fields
  std::vector<int32_t> h_vegtyp(n2d), h_soiltyp(n2d), h_slopetyp(n2d);
  if (std::fread(h_vegtyp.data(), sizeof(int32_t), n2d, f) != n2d ||
      std::fread(h_soiltyp.data(), sizeof(int32_t), n2d, f) != n2d ||
      std::fread(h_slopetyp.data(), sizeof(int32_t), n2d, f) != n2d) {
    std::fprintf(stderr, "Failed reading integer fields from %s\n", infile.c_str());
    std::fclose(f);
    return 1;
  }

  // Read 2D double fields (surface + forcing)
  std::vector<double> h_shdfac(n2d), h_shdmin(n2d), h_shdmax(n2d);
  std::vector<double> h_alb(n2d), h_snoalb(n2d), h_tbot(n2d), h_z0brd(n2d), h_embrd(n2d);
  std::vector<double> h_lwdn(n2d), h_soldn(n2d), h_solnet(n2d), h_sfcprs(n2d), h_prcp(n2d);
  std::vector<double> h_sfctmp(n2d), h_th2(n2d), h_q2(n2d), h_sfcspd(n2d), h_q2sat(n2d), h_dqsdt2(n2d), h_ffrozp(n2d);

  auto read_2d = [&](std::vector<double>& v, const char* name) {
    if (std::fread(v.data(), sizeof(double), n2d, f) != n2d) {
      std::fprintf(stderr, "Failed reading %s from %s\n", name, infile.c_str());
      std::fclose(f);
      std::exit(1);
    }
  };

  read_2d(h_shdfac, "shdfac");
  read_2d(h_shdmin, "shdmin");
  read_2d(h_shdmax, "shdmax");
  read_2d(h_alb, "alb");
  read_2d(h_snoalb, "snoalb");
  read_2d(h_tbot, "tbot");
  read_2d(h_z0brd, "z0brd");
  read_2d(h_embrd, "embrd");
  read_2d(h_lwdn, "lwdn");
  read_2d(h_soldn, "soldn");
  read_2d(h_solnet, "solnet");
  read_2d(h_sfcprs, "sfcprs");
  read_2d(h_prcp, "prcp");
  read_2d(h_sfctmp, "sfctmp");
  read_2d(h_th2, "th2");
  read_2d(h_q2, "q2");
  read_2d(h_sfcspd, "sfcspd");
  read_2d(h_q2sat, "q2sat");
  read_2d(h_dqsdt2, "dqsdt2");
  read_2d(h_ffrozp, "ffrozp");

  // Read state fields (inout)
  std::vector<double> h_cmc(n2d), h_t1(n2d), h_snowh(n2d), h_sneqv(n2d), h_albedo(n2d);
  std::vector<double> h_ch(n2d), h_cm(n2d), h_emissi(n2d), h_z0(n2d), h_xlai(n2d), h_snotime1(n2d), h_ribb(n2d);

  read_2d(h_cmc, "cmc");
  read_2d(h_t1, "t1");
  read_2d(h_snowh, "snowh");
  read_2d(h_sneqv, "sneqv");
  read_2d(h_albedo, "albedo");
  read_2d(h_ch, "ch");
  read_2d(h_cm, "cm");
  read_2d(h_emissi, "emissi");
  read_2d(h_z0, "z0");
  read_2d(h_xlai, "xlai");
  read_2d(h_snotime1, "snotime1");
  read_2d(h_ribb, "ribb");

  // Read 3D state fields (inout)
  std::vector<double> h_stc(n3d), h_smc(n3d), h_sh2o(n3d);
  auto read_3d = [&](std::vector<double>& v, const char* name) {
    if (std::fread(v.data(), sizeof(double), n3d, f) != n3d) {
      std::fprintf(stderr, "Failed reading %s from %s\n", name, infile.c_str());
      std::fclose(f);
      std::exit(1);
    }
  };

  read_3d(h_stc, "stc");
  read_3d(h_smc, "smc");
  read_3d(h_sh2o, "sh2o");
  std::fclose(f);

  // Allocate device memory
  DevBuffer<int> d_vegtyp(n2d), d_soiltyp(n2d), d_slopetyp(n2d);
  copy_to_dev(d_vegtyp.get(), h_vegtyp.data(), n2d * sizeof(int));
  copy_to_dev(d_soiltyp.get(), h_soiltyp.data(), n2d * sizeof(int));
  copy_to_dev(d_slopetyp.get(), h_slopetyp.data(), n2d * sizeof(int));

  auto make_real_2d = [&](const std::vector<double>& hd, DevBuffer<real>& d) {
    std::vector<real> hr(hd.begin(), hd.end());
    copy_to_dev(d.get(), hr.data(), n2d * sizeof(real));
  };
  auto make_real_3d = [&](const std::vector<double>& hd, DevBuffer<real>& d) {
    std::vector<real> hr(hd.begin(), hd.end());
    copy_to_dev(d.get(), hr.data(), n3d * sizeof(real));
  };

  DevBuffer<real> d_shdfac(n2d), d_shdmin(n2d), d_shdmax(n2d);
  DevBuffer<real> d_alb(n2d), d_snoalb(n2d), d_tbot(n2d), d_z0brd(n2d), d_embrd(n2d);
  DevBuffer<real> d_lwdn(n2d), d_soldn(n2d), d_solnet(n2d), d_sfcprs(n2d), d_prcp(n2d);
  DevBuffer<real> d_sfctmp(n2d), d_th2(n2d), d_q2(n2d), d_sfcspd(n2d), d_q2sat(n2d), d_dqsdt2(n2d), d_ffrozp(n2d);

  make_real_2d(h_shdfac, d_shdfac);
  make_real_2d(h_shdmin, d_shdmin);
  make_real_2d(h_shdmax, d_shdmax);
  make_real_2d(h_alb, d_alb);
  make_real_2d(h_snoalb, d_snoalb);
  make_real_2d(h_tbot, d_tbot);
  make_real_2d(h_z0brd, d_z0brd);
  make_real_2d(h_embrd, d_embrd);
  make_real_2d(h_lwdn, d_lwdn);
  make_real_2d(h_soldn, d_soldn);
  make_real_2d(h_solnet, d_solnet);
  make_real_2d(h_sfcprs, d_sfcprs);
  make_real_2d(h_prcp, d_prcp);
  make_real_2d(h_sfctmp, d_sfctmp);
  make_real_2d(h_th2, d_th2);
  make_real_2d(h_q2, d_q2);
  make_real_2d(h_sfcspd, d_sfcspd);
  make_real_2d(h_q2sat, d_q2sat);
  make_real_2d(h_dqsdt2, d_dqsdt2);
  make_real_2d(h_ffrozp, d_ffrozp);

  DevBuffer<real> d_cmc(n2d), d_t1(n2d), d_snowh(n2d), d_sneqv(n2d), d_albedo(n2d);
  DevBuffer<real> d_ch(n2d), d_cm(n2d), d_emissi(n2d), d_z0(n2d), d_xlai(n2d), d_snotime1(n2d), d_ribb(n2d);
  DevBuffer<real> d_stc(n3d), d_smc(n3d), d_sh2o(n3d);

  make_real_2d(h_cmc, d_cmc);
  make_real_2d(h_t1, d_t1);
  make_real_2d(h_snowh, d_snowh);
  make_real_2d(h_sneqv, d_sneqv);
  make_real_2d(h_albedo, d_albedo);
  make_real_2d(h_ch, d_ch);
  make_real_2d(h_cm, d_cm);
  make_real_2d(h_emissi, d_emissi);
  make_real_2d(h_z0, d_z0);
  make_real_2d(h_xlai, d_xlai);
  make_real_2d(h_snotime1, d_snotime1);
  make_real_2d(h_ribb, d_ribb);
  make_real_3d(h_stc, d_stc);
  make_real_3d(h_smc, d_smc);
  make_real_3d(h_sh2o, d_sh2o);

  // Diagnostic outputs
  DevBuffer<real> d_eta(n2d), d_sheat(n2d), d_eta_kinematic(n2d), d_fdown(n2d), d_ec(n2d), d_edir(n2d);
  DevBuffer<real> d_ett(n2d), d_esnow(n2d), d_drip(n2d), d_dew(n2d), d_beta(n2d), d_etp(n2d), d_ssoil(n2d);
  DevBuffer<real> d_flx1(n2d), d_flx2(n2d), d_flx3(n2d), d_snomlt(n2d), d_sncovr(n2d);
  DevBuffer<real> d_runoff1(n2d), d_runoff2(n2d), d_runoff3(n2d), d_rc(n2d), d_pc(n2d), d_rsmin(n2d);
  DevBuffer<real> d_soilw(n2d), d_soilm(n2d), d_q1(n2d);
  DevBuffer<real> d_et(n3d), d_smav(n3d);

  // Set up views
  int sj2d = n_cols;
  int sj3d = n_cols * nsoil;
  int sk3d = n_cols;

  auto v2c_int = [&](const DevBuffer<int>& b) { return View2<const int>{b.get(), sj2d}; };
  auto v2c = [&](const DevBuffer<real>& b) { return View2<const real>{b.get(), sj2d}; };
  auto v2 = [&](DevBuffer<real>& b) { return View2<real>{b.get(), sj2d}; };
  auto v3 = [&](DevBuffer<real>& b) { return View3<real>{b.get(), sj3d, sk3d}; };

  phys::noah::NoahKernel kernel;
  kernel.lwdn = v2c(d_lwdn);
  kernel.soldn = v2c(d_soldn);
  kernel.solnet = v2c(d_solnet);
  kernel.sfcprs = v2c(d_sfcprs);
  kernel.prcp = v2c(d_prcp);
  kernel.sfctmp = v2c(d_sfctmp);
  kernel.th2 = v2c(d_th2);
  kernel.q2 = v2c(d_q2);
  kernel.sfcspd = v2c(d_sfcspd);
  kernel.q2sat = v2c(d_q2sat);
  kernel.dqsdt2 = v2c(d_dqsdt2);
  kernel.ffrozp = v2c(d_ffrozp);

  kernel.vegtyp = v2c_int(d_vegtyp);
  kernel.soiltyp = v2c_int(d_soiltyp);
  kernel.slopetyp = v2c_int(d_slopetyp);
  kernel.shdfac = v2c(d_shdfac);
  kernel.shdmin = v2c(d_shdmin);
  kernel.shdmax = v2c(d_shdmax);
  kernel.alb = v2c(d_alb);
  kernel.snoalb = v2c(d_snoalb);
  kernel.tbot = v2c(d_tbot);
  kernel.z0brd = v2(d_z0brd);
  kernel.embrd = v2(d_embrd);

  kernel.cmc = v2(d_cmc);
  kernel.t1 = v2(d_t1);
  kernel.stc = v3(d_stc);
  kernel.smc = v3(d_smc);
  kernel.sh2o = v3(d_sh2o);
  kernel.snowh = v2(d_snowh);
  kernel.sneqv = v2(d_sneqv);
  kernel.albedo = v2(d_albedo);
  kernel.ch = v2(d_ch);
  kernel.cm = v2(d_cm);
  kernel.emissi = v2(d_emissi);
  kernel.z0 = v2(d_z0);
  kernel.xlai = v2(d_xlai);
  kernel.snotime1 = v2(d_snotime1);
  kernel.ribb = v2(d_ribb);

  kernel.eta = v2(d_eta);
  kernel.sheat = v2(d_sheat);
  kernel.eta_kinematic = v2(d_eta_kinematic);
  kernel.fdown = v2(d_fdown);
  kernel.ec = v2(d_ec);
  kernel.edir = v2(d_edir);
  kernel.et = v3(d_et);
  kernel.ett = v2(d_ett);
  kernel.esnow = v2(d_esnow);
  kernel.drip = v2(d_drip);
  kernel.dew = v2(d_dew);
  kernel.beta = v2(d_beta);
  kernel.etp = v2(d_etp);
  kernel.ssoil = v2(d_ssoil);
  kernel.flx1 = v2(d_flx1);
  kernel.flx2 = v2(d_flx2);
  kernel.flx3 = v2(d_flx3);
  kernel.snomlt = v2(d_snomlt);
  kernel.sncovr = v2(d_sncovr);
  kernel.runoff1 = v2(d_runoff1);
  kernel.runoff2 = v2(d_runoff2);
  kernel.runoff3 = v2(d_runoff3);
  kernel.rc = v2(d_rc);
  kernel.pc = v2(d_pc);
  kernel.rsmin = v2(d_rsmin);
  kernel.soilw = v2(d_soilw);
  kernel.soilm = v2(d_soilm);
  kernel.q1 = v2(d_q1);
  kernel.smav = v3(d_smav);

  kernel.dt = real(dt_in);
  kernel.zlvl = real(zlvl_in);
  kernel.nsoil = nsoil;
  for (int k = 0; k < nsoil; ++k) {
    kernel.sldpth[k] = real(h_sldpth[k]);
  }
  kernel.isurban = isurban;
  kernel.tables = tables;

  Range2D range{0, n_cols, 0, 1};

  // Run Noah scheme once
  phys::noah::run_noah(kernel, range);
  dev_sync();

  if (bench_mode) {
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < bench_iters; ++iter) {
      phys::noah::run_noah(kernel, range);
    }
    dev_sync();
    auto t1 = std::chrono::high_resolution_clock::now();
    double total_sec = std::chrono::duration<double>(t1 - t0).count();
    double ns_per_col = (total_sec * 1e9) / (double(bench_iters) * double(n_cols));
    std::cout << "Benchmark (" << backend_name() << "): "
              << n_cols << " columns x " << nsoil << " soil layers, "
              << bench_iters << " iterations: "
              << total_sec * 1e3 / bench_iters << " ms/iter, "
              << ns_per_col << " ns/column" << std::endl;
  }

  // Copy outputs back to host
  auto get_2d = [&](DevBuffer<real>& d, std::vector<double>& h) {
    std::vector<real> r(n2d);
    copy_to_host(r.data(), d.get(), n2d * sizeof(real));
    for (std::size_t i = 0; i < n2d; ++i) h[i] = double(r[i]);
  };
  auto get_3d = [&](DevBuffer<real>& d, std::vector<double>& h) {
    std::vector<real> r(n3d);
    copy_to_host(r.data(), d.get(), n3d * sizeof(real));
    for (std::size_t i = 0; i < n3d; ++i) h[i] = double(r[i]);
  };

  std::vector<double> out_eta(n2d), out_sheat(n2d), out_eta_kinematic(n2d), out_fdown(n2d);
  std::vector<double> out_ec(n2d), out_edir(n2d), out_ett(n2d), out_esnow(n2d), out_drip(n2d), out_dew(n2d);
  std::vector<double> out_beta(n2d), out_etp(n2d), out_ssoil(n2d), out_flx1(n2d), out_flx2(n2d), out_flx3(n2d);
  std::vector<double> out_snomlt(n2d), out_sncovr(n2d), out_runoff1(n2d), out_runoff2(n2d), out_runoff3(n2d);
  std::vector<double> out_rc(n2d), out_pc(n2d), out_rsmin(n2d), out_soilw(n2d), out_soilm(n2d), out_q1(n2d);

  std::vector<double> out_cmc(n2d), out_t1(n2d), out_snowh(n2d), out_sneqv(n2d), out_albedo(n2d);
  std::vector<double> out_emissi(n2d), out_z0(n2d), out_xlai(n2d), out_snotime1(n2d), out_ribb(n2d);

  std::vector<double> out_stc(n3d), out_smc(n3d), out_sh2o(n3d), out_et(n3d), out_smav(n3d);

  get_2d(d_eta, out_eta);
  get_2d(d_sheat, out_sheat);
  get_2d(d_eta_kinematic, out_eta_kinematic);
  get_2d(d_fdown, out_fdown);
  get_2d(d_ec, out_ec);
  get_2d(d_edir, out_edir);
  get_2d(d_ett, out_ett);
  get_2d(d_esnow, out_esnow);
  get_2d(d_drip, out_drip);
  get_2d(d_dew, out_dew);
  get_2d(d_beta, out_beta);
  get_2d(d_etp, out_etp);
  get_2d(d_ssoil, out_ssoil);
  get_2d(d_flx1, out_flx1);
  get_2d(d_flx2, out_flx2);
  get_2d(d_flx3, out_flx3);
  get_2d(d_snomlt, out_snomlt);
  get_2d(d_sncovr, out_sncovr);
  get_2d(d_runoff1, out_runoff1);
  get_2d(d_runoff2, out_runoff2);
  get_2d(d_runoff3, out_runoff3);
  get_2d(d_rc, out_rc);
  get_2d(d_pc, out_pc);
  get_2d(d_rsmin, out_rsmin);
  get_2d(d_soilw, out_soilw);
  get_2d(d_soilm, out_soilm);
  get_2d(d_q1, out_q1);

  get_2d(d_cmc, out_cmc);
  get_2d(d_t1, out_t1);
  get_2d(d_snowh, out_snowh);
  get_2d(d_sneqv, out_sneqv);
  get_2d(d_albedo, out_albedo);
  get_2d(d_emissi, out_emissi);
  get_2d(d_z0, out_z0);
  get_2d(d_xlai, out_xlai);
  get_2d(d_snotime1, out_snotime1);
  get_2d(d_ribb, out_ribb);

  get_3d(d_stc, out_stc);
  get_3d(d_smc, out_smc);
  get_3d(d_sh2o, out_sh2o);
  get_3d(d_et, out_et);
  get_3d(d_smav, out_smav);

  // Write outputs in exact Fortran driver binary order
  FILE* fo = std::fopen(outfile.c_str(), "wb");
  if (!fo) {
    std::fprintf(stderr, "Cannot open output file: %s\n", outfile.c_str());
    return 1;
  }

  std::fwrite(&n_cols, sizeof(int32_t), 1, fo);
  std::fwrite(&nsoil, sizeof(int32_t), 1, fo);

  auto write_2d = [&](const std::vector<double>& v) {
    std::fwrite(v.data(), sizeof(double), n2d, fo);
  };
  auto write_3d = [&](const std::vector<double>& v) {
    std::fwrite(v.data(), sizeof(double), n3d, fo);
  };

  // 2D outputs
  write_2d(out_eta);
  write_2d(out_sheat);
  write_2d(out_eta_kinematic);
  write_2d(out_fdown);
  write_2d(out_ec);
  write_2d(out_edir);
  write_2d(out_ett);
  write_2d(out_esnow);
  write_2d(out_drip);
  write_2d(out_dew);
  write_2d(out_beta);
  write_2d(out_etp);
  write_2d(out_ssoil);
  write_2d(out_flx1);
  write_2d(out_flx2);
  write_2d(out_flx3);
  write_2d(out_snomlt);
  write_2d(out_sncovr);
  write_2d(out_runoff1);
  write_2d(out_runoff2);
  write_2d(out_runoff3);
  write_2d(out_rc);
  write_2d(out_pc);
  write_2d(out_rsmin);
  write_2d(out_soilw);
  write_2d(out_soilm);
  write_2d(out_q1);

  // Inout updated states
  write_2d(out_cmc);
  write_2d(out_t1);
  write_2d(out_snowh);
  write_2d(out_sneqv);
  write_2d(out_albedo);
  write_2d(out_emissi);
  write_2d(out_z0);
  write_2d(out_xlai);
  write_2d(out_snotime1);
  write_2d(out_ribb);

  // 3D outputs (nsoil x n_cols)
  write_3d(out_stc);
  write_3d(out_smc);
  write_3d(out_sh2o);
  write_3d(out_et);
  write_3d(out_smav);

  std::fclose(fo);
  std::cout << "C++ Noah run complete: " << n_cols << " columns processed into " << outfile << std::endl;
  return 0;
}
