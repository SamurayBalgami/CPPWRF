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
#include "phys/sf_myjsfc/myjsfc.hpp"

using namespace cppwrf;
using namespace cppwrf::phys::sf_myjsfc;

int main(int argc, char** argv) {
  std::string infile = "sf_myjsfc_input.bin";
  std::string outfile = "sf_myjsfc_output_cpp.bin";
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

  int32_t n_cols = 0, nz = 0, itimestep = 1, isurban = 1, iz0tlnd = 0;
  if (std::fread(&n_cols, sizeof(int32_t), 1, f) != 1 ||
      std::fread(&nz, sizeof(int32_t), 1, f) != 1 ||
      std::fread(&itimestep, sizeof(int32_t), 1, f) != 1 ||
      std::fread(&isurban, sizeof(int32_t), 1, f) != 1 ||
      std::fread(&iz0tlnd, sizeof(int32_t), 1, f) != 1) {
    std::fprintf(stderr, "Failed reading dimensions and flags from %s\n", infile.c_str());
    std::fclose(f);
    return 1;
  }

  double p1000mb_d = 100000.0;
  if (std::fread(&p1000mb_d, sizeof(double), 1, f) != 1) {
    std::fprintf(stderr, "Failed reading p1000mb from %s\n", infile.c_str());
    std::fclose(f);
    return 1;
  }
  real p1000mb = static_cast<real>(p1000mb_d);

  std::size_t n2d = std::size_t(n_cols);
  std::size_t n3d = std::size_t(n_cols) * std::size_t(nz);
  std::size_t n3d_stag = std::size_t(n_cols) * std::size_t(nz + 1);

  std::vector<double> h_dz(n3d), h_pmid(n3d), h_pint(n3d_stag), h_th(n3d), h_t(n3d);
  std::vector<double> h_qv(n3d), h_qc(n3d), h_u(n3d), h_v(n3d), h_q2(n3d);

  std::vector<double> h_ht(n2d), h_tsk(n2d), h_xland(n2d), h_z0base(n2d), h_mavail(n2d);
  std::vector<int32_t> h_lowlyr(n2d), h_ivgtyp(n2d);

  std::vector<double> h_qsfc(n2d), h_thz0(n2d), h_qz0(n2d), h_uz0(n2d), h_vz0(n2d);
  std::vector<double> h_ustar(n2d), h_znt(n2d), h_pblh(n2d), h_rmol(n2d), h_akhs(n2d), h_akms(n2d);

  auto read_3d = [&](std::vector<double>& arr, std::size_t count) {
    if (std::fread(arr.data(), sizeof(double), count, f) != count) {
      std::fprintf(stderr, "Failed reading 3D array data\n");
      std::exit(1);
    }
  };
  auto read_2d = [&](std::vector<double>& arr) {
    if (std::fread(arr.data(), sizeof(double), n2d, f) != n2d) {
      std::fprintf(stderr, "Failed reading 2D array data\n");
      std::exit(1);
    }
  };
  auto read_2d_int = [&](std::vector<int32_t>& arr) {
    if (std::fread(arr.data(), sizeof(int32_t), n2d, f) != n2d) {
      std::fprintf(stderr, "Failed reading 2D int array data\n");
      std::exit(1);
    }
  };

  read_3d(h_dz, n3d);
  read_3d(h_pmid, n3d);
  read_3d(h_pint, n3d_stag);
  read_3d(h_th, n3d);
  read_3d(h_t, n3d);
  read_3d(h_qv, n3d);
  read_3d(h_qc, n3d);
  read_3d(h_u, n3d);
  read_3d(h_v, n3d);
  read_3d(h_q2, n3d);

  read_2d(h_ht);
  read_2d(h_tsk);
  read_2d(h_xland);
  read_2d(h_z0base);
  read_2d(h_mavail);
  read_2d_int(h_lowlyr);
  read_2d_int(h_ivgtyp);

  read_2d(h_qsfc);
  read_2d(h_thz0);
  read_2d(h_qz0);
  read_2d(h_uz0);
  read_2d(h_vz0);
  read_2d(h_ustar);
  read_2d(h_znt);
  read_2d(h_pblh);
  read_2d(h_rmol);
  read_2d(h_akhs);
  read_2d(h_akms);

  std::fclose(f);

  Buffer d_dz(n3d), d_pmid(n3d), d_pint(n3d_stag), d_th(n3d), d_t(n3d);
  Buffer d_qv(n3d), d_qc(n3d), d_u(n3d), d_v(n3d), d_q2(n3d);

  Buffer d_ht(n2d), d_tsk(n2d), d_xland(n2d), d_z0base(n2d), d_mavail(n2d);
  int* d_lowlyr = static_cast<int*>(dev_alloc(n2d * sizeof(int)));
  int* d_ivgtyp = static_cast<int*>(dev_alloc(n2d * sizeof(int)));

  Buffer d_qsfc(n2d), d_thz0(n2d), d_qz0(n2d), d_uz0(n2d), d_vz0(n2d);
  Buffer d_ustar(n2d), d_znt(n2d), d_pblh(n2d), d_rmol(n2d), d_akhs(n2d), d_akms(n2d);

  Buffer d_rib(n2d), d_chs(n2d), d_chs2(n2d), d_cqs2(n2d), d_hfx(n2d), d_qfx(n2d);
  Buffer d_flx_lh(n2d), d_flhc(n2d), d_flqc(n2d), d_qgh(n2d), d_cpm(n2d), d_ct(n2d);
  Buffer d_u10(n2d), d_v10(n2d), d_t02(n2d), d_th02(n2d), d_tshltr(n2d);
  Buffer d_th10(n2d), d_q02(n2d), d_qshltr(n2d), d_q10(n2d), d_pshltr(n2d);
  Buffer d_u10e(n2d), d_v10e(n2d);

  auto upload_real = [](Buffer& buf, const std::vector<double>& h, std::size_t sz) {
    std::vector<real> r(sz);
    for (std::size_t i = 0; i < sz; ++i) r[i] = static_cast<real>(h[i]);
    copy_to_dev(buf.get(), r.data(), sz * sizeof(real));
  };

  upload_real(d_dz, h_dz, n3d);
  upload_real(d_pmid, h_pmid, n3d);
  upload_real(d_pint, h_pint, n3d_stag);
  upload_real(d_th, h_th, n3d);
  upload_real(d_t, h_t, n3d);
  upload_real(d_qv, h_qv, n3d);
  upload_real(d_qc, h_qc, n3d);
  upload_real(d_u, h_u, n3d);
  upload_real(d_v, h_v, n3d);
  upload_real(d_q2, h_q2, n3d);

  upload_real(d_ht, h_ht, n2d);
  upload_real(d_tsk, h_tsk, n2d);
  upload_real(d_xland, h_xland, n2d);
  upload_real(d_z0base, h_z0base, n2d);
  upload_real(d_mavail, h_mavail, n2d);

  copy_to_dev(d_lowlyr, h_lowlyr.data(), n2d * sizeof(int));
  copy_to_dev(d_ivgtyp, h_ivgtyp.data(), n2d * sizeof(int));

  upload_real(d_qsfc, h_qsfc, n2d);
  upload_real(d_thz0, h_thz0, n2d);
  upload_real(d_qz0, h_qz0, n2d);
  upload_real(d_uz0, h_uz0, n2d);
  upload_real(d_vz0, h_vz0, n2d);
  upload_real(d_ustar, h_ustar, n2d);
  upload_real(d_znt, h_znt, n2d);
  upload_real(d_pblh, h_pblh, n2d);
  upload_real(d_rmol, h_rmol, n2d);
  upload_real(d_akhs, h_akhs, n2d);
  upload_real(d_akms, h_akms, n2d);

  MYJSFCTableManager tbl_mgr;
  MYJSFCTables tables = tbl_mgr.get_tables();

  View3<const real> v_dz{d_dz.get(), 0, n_cols};
  View3<const real> v_pmid{d_pmid.get(), 0, n_cols};
  View3<const real> v_pint{d_pint.get(), 0, n_cols};
  View3<const real> v_th{d_th.get(), 0, n_cols};
  View3<const real> v_t{d_t.get(), 0, n_cols};
  View3<const real> v_qv{d_qv.get(), 0, n_cols};
  View3<const real> v_qc{d_qc.get(), 0, n_cols};
  View3<const real> v_u{d_u.get(), 0, n_cols};
  View3<const real> v_v{d_v.get(), 0, n_cols};
  View3<const real> v_q2{d_q2.get(), 0, n_cols};

  View2<const real> v_ht{d_ht.get(), 0};
  View2<const real> v_tsk{d_tsk.get(), 0};
  View2<const real> v_xland{d_xland.get(), 0};
  View2<const real> v_z0base{d_z0base.get(), 0};
  View2<const real> v_mavail{d_mavail.get(), 0};
  View2<const int> v_lowlyr{d_lowlyr, 0};
  View2<const int> v_ivgtyp{d_ivgtyp, 0};

  View2<real> v_qsfc{d_qsfc.get(), 0};
  View2<real> v_thz0{d_thz0.get(), 0};
  View2<real> v_qz0{d_qz0.get(), 0};
  View2<real> v_uz0{d_uz0.get(), 0};
  View2<real> v_vz0{d_vz0.get(), 0};
  View2<real> v_ustar{d_ustar.get(), 0};
  View2<real> v_znt{d_znt.get(), 0};
  View2<real> v_pblh{d_pblh.get(), 0};
  View2<real> v_rmol{d_rmol.get(), 0};
  View2<real> v_akhs{d_akhs.get(), 0};
  View2<real> v_akms{d_akms.get(), 0};

  View2<real> v_rib{d_rib.get(), 0};
  View2<real> v_chs{d_chs.get(), 0};
  View2<real> v_chs2{d_chs2.get(), 0};
  View2<real> v_cqs2{d_cqs2.get(), 0};
  View2<real> v_hfx{d_hfx.get(), 0};
  View2<real> v_qfx{d_qfx.get(), 0};
  View2<real> v_flx_lh{d_flx_lh.get(), 0};
  View2<real> v_flhc{d_flhc.get(), 0};
  View2<real> v_flqc{d_flqc.get(), 0};
  View2<real> v_qgh{d_qgh.get(), 0};
  View2<real> v_cpm{d_cpm.get(), 0};
  View2<real> v_ct{d_ct.get(), 0};
  View2<real> v_u10{d_u10.get(), 0};
  View2<real> v_v10{d_v10.get(), 0};
  View2<real> v_t02{d_t02.get(), 0};
  View2<real> v_th02{d_th02.get(), 0};
  View2<real> v_tshltr{d_tshltr.get(), 0};
  View2<real> v_th10{d_th10.get(), 0};
  View2<real> v_q02{d_q02.get(), 0};
  View2<real> v_qshltr{d_qshltr.get(), 0};
  View2<real> v_q10{d_q10.get(), 0};
  View2<real> v_pshltr{d_pshltr.get(), 0};
  View2<real> v_u10e{d_u10e.get(), 0};
  View2<real> v_v10e{d_v10e.get(), 0};

  Range2D range{0, n_cols, 0, 1};

  run_myjsfc(
      v_ht, v_dz, v_pmid, v_pint, v_th, v_t, v_qv, v_qc, v_u, v_v, v_q2,
      v_tsk, v_xland, v_ivgtyp, v_z0base, v_mavail, v_lowlyr,
      v_qsfc, v_thz0, v_qz0, v_uz0, v_vz0, v_ustar, v_znt, v_pblh, v_rmol, v_akhs, v_akms,
      v_rib, v_chs, v_chs2, v_cqs2, v_hfx, v_qfx, v_flx_lh, v_flhc, v_flqc, v_qgh, v_cpm, v_ct,
      v_u10, v_v10, v_t02, v_th02, v_tshltr, v_th10, v_q02, v_qshltr, v_q10, v_pshltr, v_u10e, v_v10e,
      tables, range, itimestep, isurban, iz0tlnd, p1000mb, nz);

  dev_sync();

  auto download_real = [](const Buffer& buf, std::vector<double>& h, std::size_t sz) {
    std::vector<real> r(sz);
    copy_to_host(r.data(), buf.get(), sz * sizeof(real));
    for (std::size_t i = 0; i < sz; ++i) h[i] = static_cast<double>(r[i]);
  };

  download_real(d_qsfc, h_qsfc, n2d);
  download_real(d_thz0, h_thz0, n2d);
  download_real(d_qz0, h_qz0, n2d);
  download_real(d_uz0, h_uz0, n2d);
  download_real(d_vz0, h_vz0, n2d);
  download_real(d_ustar, h_ustar, n2d);
  download_real(d_znt, h_znt, n2d);
  download_real(d_pblh, h_pblh, n2d);
  download_real(d_rmol, h_rmol, n2d);
  download_real(d_akhs, h_akhs, n2d);
  download_real(d_akms, h_akms, n2d);

  std::vector<double> h_rib(n2d), h_chs(n2d), h_chs2(n2d), h_cqs2(n2d), h_hfx(n2d), h_qfx(n2d);
  std::vector<double> h_flx_lh(n2d), h_flhc(n2d), h_flqc(n2d), h_qgh(n2d), h_cpm(n2d), h_ct(n2d);
  std::vector<double> h_u10(n2d), h_v10(n2d), h_t02(n2d), h_th02(n2d), h_tshltr(n2d);
  std::vector<double> h_th10(n2d), h_q02(n2d), h_qshltr(n2d), h_q10(n2d), h_pshltr(n2d);
  std::vector<double> h_u10e(n2d), h_v10e(n2d);

  download_real(d_rib, h_rib, n2d);
  download_real(d_chs, h_chs, n2d);
  download_real(d_chs2, h_chs2, n2d);
  download_real(d_cqs2, h_cqs2, n2d);
  download_real(d_hfx, h_hfx, n2d);
  download_real(d_qfx, h_qfx, n2d);
  download_real(d_flx_lh, h_flx_lh, n2d);
  download_real(d_flhc, h_flhc, n2d);
  download_real(d_flqc, h_flqc, n2d);
  download_real(d_qgh, h_qgh, n2d);
  download_real(d_cpm, h_cpm, n2d);
  download_real(d_ct, h_ct, n2d);
  download_real(d_u10, h_u10, n2d);
  download_real(d_v10, h_v10, n2d);
  download_real(d_t02, h_t02, n2d);
  download_real(d_th02, h_th02, n2d);
  download_real(d_tshltr, h_tshltr, n2d);
  download_real(d_th10, h_th10, n2d);
  download_real(d_q02, h_q02, n2d);
  download_real(d_qshltr, h_qshltr, n2d);
  download_real(d_q10, h_q10, n2d);
  download_real(d_pshltr, h_pshltr, n2d);
  download_real(d_u10e, h_u10e, n2d);
  download_real(d_v10e, h_v10e, n2d);

  FILE* fo = std::fopen(outfile.c_str(), "wb");
  if (!fo) {
    std::fprintf(stderr, "Cannot open output file: %s\n", outfile.c_str());
    return 1;
  }

  auto write_2d = [&](const std::vector<double>& arr) {
    std::fwrite(arr.data(), sizeof(double), n2d, fo);
  };

  write_2d(h_qsfc);
  write_2d(h_thz0);
  write_2d(h_qz0);
  write_2d(h_uz0);
  write_2d(h_vz0);
  write_2d(h_ustar);
  write_2d(h_znt);
  write_2d(h_pblh);
  write_2d(h_rmol);
  write_2d(h_akhs);
  write_2d(h_akms);

  write_2d(h_rib);
  write_2d(h_chs);
  write_2d(h_chs2);
  write_2d(h_cqs2);
  write_2d(h_hfx);
  write_2d(h_qfx);
  write_2d(h_flx_lh);
  write_2d(h_flhc);
  write_2d(h_flqc);
  write_2d(h_qgh);
  write_2d(h_cpm);
  write_2d(h_ct);
  write_2d(h_u10);
  write_2d(h_v10);
  write_2d(h_t02);
  write_2d(h_th02);
  write_2d(h_tshltr);
  write_2d(h_th10);
  write_2d(h_q02);
  write_2d(h_qshltr);
  write_2d(h_q10);
  write_2d(h_pshltr);
  write_2d(h_u10e);
  write_2d(h_v10e);

  std::fclose(fo);

  dev_free(d_lowlyr);
  dev_free(d_ivgtyp);

  return 0;
}
