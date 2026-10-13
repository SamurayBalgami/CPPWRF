#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "core/field.hpp"
#include "phys/sf_sfclayrev/sfclayrev.hpp"

using namespace cppwrf;

int main(int argc, char** argv) {
  std::string infile = "sfclayrev_input.bin";
  std::string outfile = "sfclayrev_output_cpp.bin";
  bool bench_mode = false;
  int bench_iters = 10;

  bool infile_set = false;
  for (int a = 1; a < argc; ++a) {
    std::string arg = argv[a];
    if (arg == "--bench") {
      bench_mode = true;
      if (a + 1 < argc && argv[a + 1][0] != '-') {
        bench_iters = std::atoi(argv[++a]);
      }
    } else if (!infile_set) {
      infile = arg;
      infile_set = true;
    } else {
      outfile = arg;
    }
  }

  FILE* f = std::fopen(infile.c_str(), "rb");
  if (!f) {
    std::fprintf(stderr, "Cannot open input file: %s\n", infile.c_str());
    return 1;
  }

  int32_t n_cols = 0, nz = 0;
  if (std::fread(&n_cols, sizeof(int32_t), 1, f) != 1 ||
      std::fread(&nz, sizeof(int32_t), 1, f) != 1) {
    std::fprintf(stderr, "Failed reading dimensions from %s\n", infile.c_str());
    std::fclose(f);
    return 1;
  }

  int32_t isfflx_i = 1, shalwater_z0_i = 0, isftcflx_i = 0, iz0tlnd_i = 0, scm_force_flux_i = 0;
  if (std::fread(&isfflx_i, sizeof(int32_t), 1, f) != 1 ||
      std::fread(&shalwater_z0_i, sizeof(int32_t), 1, f) != 1 ||
      std::fread(&isftcflx_i, sizeof(int32_t), 1, f) != 1 ||
      std::fread(&iz0tlnd_i, sizeof(int32_t), 1, f) != 1 ||
      std::fread(&scm_force_flux_i, sizeof(int32_t), 1, f) != 1) {
    std::fprintf(stderr, "Failed reading option flags from %s\n", infile.c_str());
    std::fclose(f);
    return 1;
  }

  double cp, g, rovcp, r, xlv, svp1, svp2, svp3, svpt0, ep1, ep2, karman, p1000mb;
  if (std::fread(&cp, sizeof(double), 1, f) != 1 ||
      std::fread(&g, sizeof(double), 1, f) != 1 ||
      std::fread(&rovcp, sizeof(double), 1, f) != 1 ||
      std::fread(&r, sizeof(double), 1, f) != 1 ||
      std::fread(&xlv, sizeof(double), 1, f) != 1 ||
      std::fread(&svp1, sizeof(double), 1, f) != 1 ||
      std::fread(&svp2, sizeof(double), 1, f) != 1 ||
      std::fread(&svp3, sizeof(double), 1, f) != 1 ||
      std::fread(&svpt0, sizeof(double), 1, f) != 1 ||
      std::fread(&ep1, sizeof(double), 1, f) != 1 ||
      std::fread(&ep2, sizeof(double), 1, f) != 1 ||
      std::fread(&karman, sizeof(double), 1, f) != 1 ||
      std::fread(&p1000mb, sizeof(double), 1, f) != 1) {
    std::fprintf(stderr, "Failed reading constants from %s\n", infile.c_str());
    std::fclose(f);
    return 1;
  }

  std::size_t n2d = std::size_t(n_cols);

  // Read FP64 double arrays from file
  std::vector<double> h_u1d(n2d), h_v1d(n2d), h_t1d(n2d), h_qv1d(n2d), h_p1d(n2d), h_dz8w1d(n2d);
  std::vector<double> h_psfcpa(n2d), h_tsk(n2d), h_xland(n2d), h_dx(n2d), h_mavail(n2d), h_pblh(n2d);
  std::vector<double> h_lakemask(n2d), h_water_depth(n2d);

  std::vector<double> h_regime(n2d), h_hfx(n2d), h_qfx(n2d), h_qsfc(n2d), h_mol(n2d), h_rmol(n2d);
  std::vector<double> h_gz1oz0(n2d), h_wspd(n2d), h_br(n2d), h_psim(n2d), h_psih(n2d), h_fm(n2d), h_fh(n2d);
  std::vector<double> h_znt(n2d), h_zol(n2d), h_ust(n2d), h_cpm(n2d), h_chs2(n2d), h_cqs2(n2d), h_chs(n2d);
  std::vector<double> h_flhc(n2d), h_flqc(n2d), h_qgh(n2d), h_ustm(n2d);

  auto read_arr = [&](std::vector<double>& arr) {
    if (std::fread(arr.data(), sizeof(double), n2d, f) != n2d) {
      std::fprintf(stderr, "Failed reading array data\n");
      std::exit(1);
    }
  };

  read_arr(h_u1d); read_arr(h_v1d); read_arr(h_t1d); read_arr(h_qv1d); read_arr(h_p1d); read_arr(h_dz8w1d);
  read_arr(h_psfcpa); read_arr(h_tsk); read_arr(h_xland); read_arr(h_dx); read_arr(h_mavail); read_arr(h_pblh);
  read_arr(h_lakemask); read_arr(h_water_depth);

  read_arr(h_regime); read_arr(h_hfx); read_arr(h_qfx); read_arr(h_qsfc); read_arr(h_mol); read_arr(h_rmol);
  read_arr(h_gz1oz0); read_arr(h_wspd); read_arr(h_br); read_arr(h_psim); read_arr(h_psih); read_arr(h_fm); read_arr(h_fh);
  read_arr(h_znt); read_arr(h_zol); read_arr(h_ust); read_arr(h_cpm); read_arr(h_chs2); read_arr(h_cqs2); read_arr(h_chs);
  read_arr(h_flhc); read_arr(h_flqc); read_arr(h_qgh); read_arr(h_ustm);
  std::fclose(f);

  // Convert to working precision `real`
  auto to_real = [](const std::vector<double>& in) {
    return std::vector<real>(in.begin(), in.end());
  };

  auto r_u1d = to_real(h_u1d), r_v1d = to_real(h_v1d), r_t1d = to_real(h_t1d), r_qv1d = to_real(h_qv1d), r_p1d = to_real(h_p1d), r_dz8w1d = to_real(h_dz8w1d);
  auto r_psfcpa = to_real(h_psfcpa), r_tsk = to_real(h_tsk), r_xland = to_real(h_xland), r_dx = to_real(h_dx), r_mavail = to_real(h_mavail), r_pblh = to_real(h_pblh);
  auto r_lakemask = to_real(h_lakemask), r_water_depth = to_real(h_water_depth);

  auto r_regime = to_real(h_regime), r_hfx = to_real(h_hfx), r_qfx = to_real(h_qfx), r_qsfc = to_real(h_qsfc), r_mol = to_real(h_mol), r_rmol = to_real(h_rmol);
  auto r_gz1oz0 = to_real(h_gz1oz0), r_wspd = to_real(h_wspd), r_br = to_real(h_br), r_psim = to_real(h_psim), r_psih = to_real(h_psih), r_fm = to_real(h_fm), r_fh = to_real(h_fh);
  auto r_znt = to_real(h_znt), r_zol = to_real(h_zol), r_ust = to_real(h_ust), r_cpm = to_real(h_cpm), r_chs2 = to_real(h_chs2), r_cqs2 = to_real(h_cqs2), r_chs = to_real(h_chs);
  auto r_flhc = to_real(h_flhc), r_flqc = to_real(h_flqc), r_qgh = to_real(h_qgh), r_ustm = to_real(h_ustm);

  std::vector<real> r_lh(n2d, real(0)), r_u10(n2d, real(0)), r_v10(n2d, real(0)), r_th2(n2d, real(0)), r_t2(n2d, real(0)), r_q2(n2d, real(0));
  std::vector<real> r_ck(n2d, real(0)), r_cka(n2d, real(0)), r_cd(n2d, real(0)), r_cda(n2d, real(0));

  // Initialize lookup tables
  phys::sfclayrev::SFClayRevTables<real> h_tables;
  phys::sfclayrev::init_sfclayrev_tables(h_tables);

  Buffer d_tables(sizeof(phys::sfclayrev::SFClayRevTables<real>) / sizeof(real) + 1);
  copy_to_dev(d_tables.get(), &h_tables, sizeof(phys::sfclayrev::SFClayRevTables<real>));

  // Allocate device memory using Buffer
  Buffer d_u1d(n2d), d_v1d(n2d), d_t1d(n2d), d_qv1d(n2d), d_p1d(n2d), d_dz8w1d(n2d);
  Buffer d_psfcpa(n2d), d_tsk(n2d), d_xland(n2d), d_dx(n2d), d_mavail(n2d), d_pblh(n2d);
  Buffer d_lakemask(n2d), d_water_depth(n2d);

  Buffer d_regime(n2d), d_hfx(n2d), d_qfx(n2d), d_qsfc(n2d), d_mol(n2d), d_rmol(n2d);
  Buffer d_gz1oz0(n2d), d_wspd(n2d), d_br(n2d), d_psim(n2d), d_psih(n2d), d_fm(n2d), d_fh(n2d);
  Buffer d_znt(n2d), d_zol(n2d), d_ust(n2d), d_cpm(n2d), d_chs2(n2d), d_cqs2(n2d), d_chs(n2d);
  Buffer d_flhc(n2d), d_flqc(n2d), d_qgh(n2d), d_ustm(n2d);

  Buffer d_lh(n2d), d_u10(n2d), d_v10(n2d), d_th2(n2d), d_t2(n2d), d_q2(n2d);
  Buffer d_ck(n2d), d_cka(n2d), d_cd(n2d), d_cda(n2d);

  auto copy_up = [&](Buffer& b, const std::vector<real>& v) {
    copy_to_dev(b.get(), v.data(), n2d * sizeof(real));
  };

  copy_up(d_u1d, r_u1d); copy_up(d_v1d, r_v1d); copy_up(d_t1d, r_t1d); copy_up(d_qv1d, r_qv1d); copy_up(d_p1d, r_p1d); copy_up(d_dz8w1d, r_dz8w1d);
  copy_up(d_psfcpa, r_psfcpa); copy_up(d_tsk, r_tsk); copy_up(d_xland, r_xland); copy_up(d_dx, r_dx); copy_up(d_mavail, r_mavail); copy_up(d_pblh, r_pblh);
  copy_up(d_lakemask, r_lakemask); copy_up(d_water_depth, r_water_depth);

  copy_up(d_regime, r_regime); copy_up(d_hfx, r_hfx); copy_up(d_qfx, r_qfx); copy_up(d_qsfc, r_qsfc); copy_up(d_mol, r_mol); copy_up(d_rmol, r_rmol);
  copy_up(d_gz1oz0, r_gz1oz0); copy_up(d_wspd, r_wspd); copy_up(d_br, r_br); copy_up(d_psim, r_psim); copy_up(d_psih, r_psih); copy_up(d_fm, r_fm); copy_up(d_fh, r_fh);
  copy_up(d_znt, r_znt); copy_up(d_zol, r_zol); copy_up(d_ust, r_ust); copy_up(d_cpm, r_cpm); copy_up(d_chs2, r_chs2); copy_up(d_cqs2, r_cqs2); copy_up(d_chs, r_chs);
  copy_up(d_flhc, r_flhc); copy_up(d_flqc, r_flqc); copy_up(d_qgh, r_qgh); copy_up(d_ustm, r_ustm);

  copy_up(d_lh, r_lh); copy_up(d_u10, r_u10); copy_up(d_v10, r_v10); copy_up(d_th2, r_th2); copy_up(d_t2, r_t2); copy_up(d_q2, r_q2);
  copy_up(d_ck, r_ck); copy_up(d_cka, r_cka); copy_up(d_cd, r_cd); copy_up(d_cda, r_cda);

  View2<const real> v_u1d{d_u1d.get(), n_cols};
  View2<const real> v_v1d{d_v1d.get(), n_cols};
  View2<const real> v_t1d{d_t1d.get(), n_cols};
  View2<const real> v_qv1d{d_qv1d.get(), n_cols};
  View2<const real> v_p1d{d_p1d.get(), n_cols};
  View2<const real> v_dz8w1d{d_dz8w1d.get(), n_cols};
  View2<const real> v_psfcpa{d_psfcpa.get(), n_cols};
  View2<const real> v_tsk{d_tsk.get(), n_cols};
  View2<const real> v_xland{d_xland.get(), n_cols};
  View2<const real> v_dx{d_dx.get(), n_cols};
  View2<const real> v_mavail{d_mavail.get(), n_cols};
  View2<const real> v_pblh{d_pblh.get(), n_cols};
  View2<const real> v_lakemask{d_lakemask.get(), n_cols};
  View2<const real> v_water_depth{d_water_depth.get(), n_cols};

  View2<real> v_regime{d_regime.get(), n_cols};
  View2<real> v_hfx{d_hfx.get(), n_cols};
  View2<real> v_qfx{d_qfx.get(), n_cols};
  View2<real> v_qsfc{d_qsfc.get(), n_cols};
  View2<real> v_mol{d_mol.get(), n_cols};
  View2<real> v_rmol{d_rmol.get(), n_cols};
  View2<real> v_gz1oz0{d_gz1oz0.get(), n_cols};
  View2<real> v_wspd{d_wspd.get(), n_cols};
  View2<real> v_br{d_br.get(), n_cols};
  View2<real> v_psim{d_psim.get(), n_cols};
  View2<real> v_psih{d_psih.get(), n_cols};
  View2<real> v_fm{d_fm.get(), n_cols};
  View2<real> v_fh{d_fh.get(), n_cols};
  View2<real> v_znt{d_znt.get(), n_cols};
  View2<real> v_zol{d_zol.get(), n_cols};
  View2<real> v_ust{d_ust.get(), n_cols};
  View2<real> v_cpm{d_cpm.get(), n_cols};
  View2<real> v_chs2{d_chs2.get(), n_cols};
  View2<real> v_cqs2{d_cqs2.get(), n_cols};
  View2<real> v_chs{d_chs.get(), n_cols};
  View2<real> v_flhc{d_flhc.get(), n_cols};
  View2<real> v_flqc{d_flqc.get(), n_cols};
  View2<real> v_qgh{d_qgh.get(), n_cols};
  View2<real> v_ustm{d_ustm.get(), n_cols};

  View2<real> v_lh{d_lh.get(), n_cols};
  View2<real> v_u10{d_u10.get(), n_cols};
  View2<real> v_v10{d_v10.get(), n_cols};
  View2<real> v_th2{d_th2.get(), n_cols};
  View2<real> v_t2{d_t2.get(), n_cols};
  View2<real> v_q2{d_q2.get(), n_cols};
  View2<real> v_ck{d_ck.get(), n_cols};
  View2<real> v_cka{d_cka.get(), n_cols};
  View2<real> v_cd{d_cd.get(), n_cols};
  View2<real> v_cda{d_cda.get(), n_cols};

  const auto* tbl_ptr = reinterpret_cast<const phys::sfclayrev::SFClayRevTables<real>*>(d_tables.get());
  Range2D range{0, n_cols, 0, 1};

  bool isfflx = (isfflx_i == 1);
  bool shalwater_z0 = (shalwater_z0_i == 1);
  bool scm_force_flux = (scm_force_flux_i == 1);

  // Run kernel once
  phys::sfclayrev::run_sfclayrev(
      v_u1d, v_v1d, v_t1d, v_qv1d, v_p1d, v_dz8w1d, v_psfcpa, v_tsk, v_xland, v_dx,
      v_mavail, v_pblh, v_lakemask, v_water_depth,
      v_regime, v_hfx, v_qfx, v_qsfc, v_mol, v_rmol, v_gz1oz0, v_wspd, v_br,
      v_psim, v_psih, v_fm, v_fh, v_znt, v_zol, v_ust, v_cpm, v_chs2, v_cqs2,
      v_chs, v_flhc, v_flqc, v_qgh, v_ustm,
      v_lh, v_u10, v_v10, v_th2, v_t2, v_q2, v_ck, v_cka, v_cd, v_cda,
      tbl_ptr, range,
      isfflx, shalwater_z0, isftcflx_i, iz0tlnd_i, scm_force_flux,
      real(cp), real(g), real(rovcp), real(r), real(xlv),
      real(svp1), real(svp2), real(svp3), real(svpt0), real(ep1), real(ep2), real(karman), real(p1000mb));
  dev_sync();

  if (bench_mode) {
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < bench_iters; ++iter) {
      phys::sfclayrev::run_sfclayrev(
          v_u1d, v_v1d, v_t1d, v_qv1d, v_p1d, v_dz8w1d, v_psfcpa, v_tsk, v_xland, v_dx,
          v_mavail, v_pblh, v_lakemask, v_water_depth,
          v_regime, v_hfx, v_qfx, v_qsfc, v_mol, v_rmol, v_gz1oz0, v_wspd, v_br,
          v_psim, v_psih, v_fm, v_fh, v_znt, v_zol, v_ust, v_cpm, v_chs2, v_cqs2,
          v_chs, v_flhc, v_flqc, v_qgh, v_ustm,
          v_lh, v_u10, v_v10, v_th2, v_t2, v_q2, v_ck, v_cka, v_cd, v_cda,
          tbl_ptr, range,
          isfflx, shalwater_z0, isftcflx_i, iz0tlnd_i, scm_force_flux,
          real(cp), real(g), real(rovcp), real(r), real(xlv),
          real(svp1), real(svp2), real(svp3), real(svpt0), real(ep1), real(ep2), real(karman), real(p1000mb));
    }
    dev_sync();
    auto t1 = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    double avg_ms = total_ms / bench_iters;
    double ns_per_col = (avg_ms * 1e6) / n_cols;
    std::printf("BENCHMARK: %d cols x 1 level | %d iters | avg %.3f ms | %.1f ns/col\n",
                n_cols, bench_iters, avg_ms, ns_per_col);
  }

  // Copy outputs back to host
  auto copy_down = [&](std::vector<real>& v, Buffer& b) {
    copy_to_host(v.data(), b.get(), n2d * sizeof(real));
  };

  copy_down(r_regime, d_regime); copy_down(r_hfx, d_hfx); copy_down(r_qfx, d_qfx); copy_down(r_lh, d_lh);
  copy_down(r_qsfc, d_qsfc); copy_down(r_mol, d_mol); copy_down(r_rmol, d_rmol); copy_down(r_gz1oz0, d_gz1oz0);
  copy_down(r_wspd, d_wspd); copy_down(r_br, d_br); copy_down(r_psim, d_psim); copy_down(r_psih, d_psih);
  copy_down(r_fm, d_fm); copy_down(r_fh, d_fh); copy_down(r_znt, d_znt); copy_down(r_zol, d_zol);
  copy_down(r_ust, d_ust); copy_down(r_cpm, d_cpm); copy_down(r_chs2, d_chs2); copy_down(r_cqs2, d_cqs2);
  copy_down(r_chs, d_chs); copy_down(r_flhc, d_flhc); copy_down(r_flqc, d_flqc); copy_down(r_qgh, d_qgh);
  copy_down(r_ustm, d_ustm); copy_down(r_u10, d_u10); copy_down(r_v10, d_v10); copy_down(r_th2, d_th2);
  copy_down(r_t2, d_t2); copy_down(r_q2, d_q2); copy_down(r_ck, d_ck); copy_down(r_cka, d_cka);
  copy_down(r_cd, d_cd); copy_down(r_cda, d_cda);

  // Convert back to double for output file
  auto to_double = [](const std::vector<real>& in) {
    return std::vector<double>(in.begin(), in.end());
  };

  auto out_regime = to_double(r_regime), out_hfx = to_double(r_hfx), out_qfx = to_double(r_qfx), out_lh = to_double(r_lh);
  auto out_qsfc = to_double(r_qsfc), out_mol = to_double(r_mol), out_rmol = to_double(r_rmol), out_gz1oz0 = to_double(r_gz1oz0);
  auto out_wspd = to_double(r_wspd), out_br = to_double(r_br), out_psim = to_double(r_psim), out_psih = to_double(r_psih);
  auto out_fm = to_double(r_fm), out_fh = to_double(r_fh), out_znt = to_double(r_znt), out_zol = to_double(r_zol);
  auto out_ust = to_double(r_ust), out_cpm = to_double(r_cpm), out_chs2 = to_double(r_chs2), out_cqs2 = to_double(r_cqs2);
  auto out_chs = to_double(r_chs), out_flhc = to_double(r_flhc), out_flqc = to_double(r_flqc), out_qgh = to_double(r_qgh);
  auto out_ustm = to_double(r_ustm), out_u10 = to_double(r_u10), out_v10 = to_double(r_v10), out_th2 = to_double(r_th2);
  auto out_t2 = to_double(r_t2), out_q2 = to_double(r_q2), out_ck = to_double(r_ck), out_cka = to_double(r_cka);
  auto out_cd = to_double(r_cd), out_cda = to_double(r_cda);

  FILE* fo = std::fopen(outfile.c_str(), "wb");
  if (!fo) {
    std::fprintf(stderr, "Cannot open output file: %s\n", outfile.c_str());
    return 1;
  }

  auto write_arr = [&](const std::vector<double>& arr) {
    std::fwrite(arr.data(), sizeof(double), n2d, fo);
  };

  write_arr(out_regime); write_arr(out_hfx); write_arr(out_qfx); write_arr(out_lh);
  write_arr(out_qsfc); write_arr(out_mol); write_arr(out_rmol); write_arr(out_gz1oz0);
  write_arr(out_wspd); write_arr(out_br); write_arr(out_psim); write_arr(out_psih);
  write_arr(out_fm); write_arr(out_fh); write_arr(out_znt); write_arr(out_zol);
  write_arr(out_ust); write_arr(out_cpm); write_arr(out_chs2); write_arr(out_cqs2);
  write_arr(out_chs); write_arr(out_flhc); write_arr(out_flqc); write_arr(out_qgh);
  write_arr(out_ustm); write_arr(out_u10); write_arr(out_v10); write_arr(out_th2);
  write_arr(out_t2); write_arr(out_q2); write_arr(out_ck); write_arr(out_cka);
  write_arr(out_cd); write_arr(out_cda);

  std::fclose(fo);
  return 0;
}
