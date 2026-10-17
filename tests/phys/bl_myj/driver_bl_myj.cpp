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
#include "phys/bl_myj/myjpbl.hpp"

using namespace cppwrf;
using namespace cppwrf::phys::myjpbl;

int main(int argc, char** argv) {
  std::string infile = "bl_myj_input.bin";
  std::string outfile = "bl_myj_output_cpp.bin";
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

  int32_t n_cols = 0, nz = 0, stepbl = 1;
  double dt_d = 12.0;

  if (std::fread(&n_cols, sizeof(int32_t), 1, f) != 1 ||
      std::fread(&nz, sizeof(int32_t), 1, f) != 1 ||
      std::fread(&stepbl, sizeof(int32_t), 1, f) != 1 ||
      std::fread(&dt_d, sizeof(double), 1, f) != 1) {
    std::fprintf(stderr, "Failed reading header from %s\n", infile.c_str());
    std::fclose(f);
    return 1;
  }
  real dt = static_cast<real>(dt_d);

  std::size_t n2d = std::size_t(n_cols);
  std::size_t n3d = std::size_t(n_cols) * std::size_t(nz);
  std::size_t n3d_stag = std::size_t(n_cols) * std::size_t(nz + 1);

  std::vector<double> h_dz(n3d), h_pmid(n3d), h_pint(n3d_stag), h_th(n3d), h_t(n3d);
  std::vector<double> h_exner(n3d), h_qv(n3d), h_qc(n3d), h_qi(n3d);
  std::vector<double> h_u(n3d), h_v(n3d), h_rho(n3d), h_tke(n3d), h_exch_h(n3d);

  std::vector<double> h_ht(n2d), h_tsk(n2d), h_xland(n2d), h_sice(n2d), h_snow(n2d);
  std::vector<double> h_chklowq(n2d), h_ustar(n2d), h_znt(n2d), h_akhs(n2d), h_akms(n2d), h_elflx(n2d);
  std::vector<int32_t> h_lowlyr(n2d);

  std::vector<double> h_qsfc(n2d), h_thz0(n2d), h_qz0(n2d), h_uz0(n2d), h_vz0(n2d), h_ct(n2d);

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
  read_3d(h_exner, n3d);
  read_3d(h_qv, n3d);
  read_3d(h_qc, n3d);
  read_3d(h_qi, n3d);
  read_3d(h_u, n3d);
  read_3d(h_v, n3d);
  read_3d(h_rho, n3d);
  read_3d(h_tke, n3d);
  read_3d(h_exch_h, n3d);

  read_2d(h_ht);
  read_2d(h_tsk);
  read_2d(h_xland);
  read_2d(h_sice);
  read_2d(h_snow);
  read_2d(h_chklowq);
  read_2d_int(h_lowlyr);
  read_2d(h_ustar);
  read_2d(h_znt);
  read_2d(h_akhs);
  read_2d(h_akms);
  read_2d(h_elflx);

  read_2d(h_qsfc);
  read_2d(h_thz0);
  read_2d(h_qz0);
  read_2d(h_uz0);
  read_2d(h_vz0);
  read_2d(h_ct);

  std::fclose(f);

  Buffer d_dz(n3d), d_pmid(n3d), d_pint(n3d_stag), d_th(n3d), d_t(n3d);
  Buffer d_exner(n3d), d_qv(n3d), d_qc(n3d), d_qi(n3d), d_u(n3d), d_v(n3d), d_rho(n3d);
  Buffer d_tke_myj(n3d), d_exch_h(n3d);

  Buffer d_ht(n2d), d_tsk(n2d), d_xland(n2d), d_sice(n2d), d_snow(n2d);
  Buffer d_chklowq(n2d), d_ustar(n2d), d_znt(n2d), d_akhs(n2d), d_akms(n2d), d_elflx(n2d);
  int* d_lowlyr = static_cast<int*>(dev_alloc(n2d * sizeof(int)));

  Buffer d_qsfc(n2d), d_thz0(n2d), d_thz0_in(n2d), d_qz0(n2d), d_qz0_in(n2d);
  Buffer d_uz0(n2d), d_vz0(n2d), d_ct(n2d);

  Buffer d_pblh(n2d), d_mixht(n2d);
  int* d_kpbl = static_cast<int*>(dev_alloc(n2d * sizeof(int)));

  Buffer d_el_myj(n3d), d_rublten(n3d), d_rvblten(n3d), d_rthblten(n3d);
  Buffer d_rqvblten(n3d), d_rqcblten(n3d), d_rqiblten(n3d);

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
  upload_real(d_exner, h_exner, n3d);
  upload_real(d_qv, h_qv, n3d);
  upload_real(d_qc, h_qc, n3d);
  upload_real(d_qi, h_qi, n3d);
  upload_real(d_u, h_u, n3d);
  upload_real(d_v, h_v, n3d);
  upload_real(d_rho, h_rho, n3d);
  upload_real(d_tke_myj, h_tke, n3d);
  upload_real(d_exch_h, h_exch_h, n3d);

  upload_real(d_ht, h_ht, n2d);
  upload_real(d_tsk, h_tsk, n2d);
  upload_real(d_xland, h_xland, n2d);
  upload_real(d_sice, h_sice, n2d);
  upload_real(d_snow, h_snow, n2d);
  upload_real(d_chklowq, h_chklowq, n2d);
  copy_to_dev(d_lowlyr, h_lowlyr.data(), n2d * sizeof(int));
  upload_real(d_ustar, h_ustar, n2d);
  upload_real(d_znt, h_znt, n2d);
  upload_real(d_akhs, h_akhs, n2d);
  upload_real(d_akms, h_akms, n2d);
  upload_real(d_elflx, h_elflx, n2d);

  upload_real(d_qsfc, h_qsfc, n2d);
  upload_real(d_thz0, h_thz0, n2d);
  upload_real(d_thz0_in, h_thz0, n2d);
  upload_real(d_qz0, h_qz0, n2d);
  upload_real(d_qz0_in, h_qz0, n2d);
  upload_real(d_uz0, h_uz0, n2d);
  upload_real(d_vz0, h_vz0, n2d);
  upload_real(d_ct, h_ct, n2d);

  Range2D range{0, n_cols, 0, 1};

  View2<const real> v_ht{d_ht.get(), 0};
  View3<const real> v_dz{d_dz.get(), 0, n_cols};
  View3<const real> v_pmid{d_pmid.get(), 0, n_cols};
  View3<const real> v_pint{d_pint.get(), 0, n_cols};
  View3<const real> v_th{d_th.get(), 0, n_cols};
  View3<const real> v_t{d_t.get(), 0, n_cols};
  View3<const real> v_exner{d_exner.get(), 0, n_cols};
  View3<const real> v_qv{d_qv.get(), 0, n_cols};
  View3<const real> v_qc{d_qc.get(), 0, n_cols};
  View3<const real> v_qi{d_qi.get(), 0, n_cols};
  View3<const real> v_u{d_u.get(), 0, n_cols};
  View3<const real> v_v{d_v.get(), 0, n_cols};
  View3<const real> v_rho{d_rho.get(), 0, n_cols};

  View2<const real> v_tsk{d_tsk.get(), 0};
  View2<const real> v_chklowq{d_chklowq.get(), 0};
  View2<const real> v_thz0_in{d_thz0_in.get(), 0};
  View2<const real> v_qz0_in{d_qz0_in.get(), 0};
  View2<const real> v_uz0{d_uz0.get(), 0};
  View2<const real> v_vz0{d_vz0.get(), 0};
  View2<const int> v_lowlyr{d_lowlyr, 0};
  View2<const real> v_xland{d_xland.get(), 0};
  View2<const real> v_sice{d_sice.get(), 0};
  View2<const real> v_snow{d_snow.get(), 0};
  View2<const real> v_ustar{d_ustar.get(), 0};
  View2<const real> v_znt{d_znt.get(), 0};
  View2<const real> v_akhs{d_akhs.get(), 0};
  View2<const real> v_akms{d_akms.get(), 0};
  View2<const real> v_elflx{d_elflx.get(), 0};

  View2<real> v_qsfc{d_qsfc.get(), 0};
  View2<real> v_thz0{d_thz0.get(), 0};
  View2<real> v_qz0{d_qz0.get(), 0};
  View2<real> v_ct{d_ct.get(), 0};
  View3<real> v_tke_myj{d_tke_myj.get(), 0, n_cols};
  View3<real> v_exch_h{d_exch_h.get(), 0, n_cols};

  View3<real> v_el_myj{d_el_myj.get(), 0, n_cols};
  View2<real> v_pblh{d_pblh.get(), 0};
  View2<int> v_kpbl{d_kpbl, 0};
  View2<real> v_mixht{d_mixht.get(), 0};

  View3<real> v_rublten{d_rublten.get(), 0, n_cols};
  View3<real> v_rvblten{d_rvblten.get(), 0, n_cols};
  View3<real> v_rthblten{d_rthblten.get(), 0, n_cols};
  View3<real> v_rqvblten{d_rqvblten.get(), 0, n_cols};
  View3<real> v_rqcblten{d_rqcblten.get(), 0, n_cols};
  View3<real> v_rqiblten{d_rqiblten.get(), 0, n_cols};

  run_myjpbl(
      v_ht, v_dz, v_pmid, v_pint, v_th, v_t, v_exner, v_qv, v_qc, v_qi,
      v_u, v_v, v_rho, v_tsk, v_chklowq, v_thz0_in, v_qz0_in, v_uz0, v_vz0,
      v_lowlyr, v_xland, v_sice, v_snow, v_ustar, v_znt, v_akhs, v_akms, v_elflx,
      v_qsfc, v_thz0, v_qz0, v_ct, v_tke_myj, v_exch_h, v_el_myj,
      v_pblh, v_kpbl, v_mixht,
      v_rublten, v_rvblten, v_rthblten, v_rqvblten, v_rqcblten, v_rqiblten,
      range, dt, stepbl, nz, true);
  dev_sync();

  auto download_real = [](const Buffer& buf, std::vector<double>& h, std::size_t sz) {
    std::vector<real> r(sz);
    copy_to_host(r.data(), buf.get(), sz * sizeof(real));
    for (std::size_t i = 0; i < sz; ++i) h[i] = static_cast<double>(r[i]);
  };

  download_real(d_qsfc, h_qsfc, n2d);
  download_real(d_thz0, h_thz0, n2d);
  download_real(d_qz0, h_qz0, n2d);
  download_real(d_ct, h_ct, n2d);
  std::vector<double> h_pblh(n2d), h_mixht(n2d);
  download_real(d_pblh, h_pblh, n2d);
  download_real(d_mixht, h_mixht, n2d);

  std::vector<int> temp_kpbl(n2d);
  copy_to_host(temp_kpbl.data(), d_kpbl, n2d * sizeof(int));
  std::vector<double> h_kpbl(n2d);
  for (std::size_t i = 0; i < n2d; ++i) h_kpbl[i] = static_cast<double>(temp_kpbl[i]);

  std::vector<double> out_tke(n3d), out_exch_h(n3d), out_el_myj(n3d);
  std::vector<double> out_rublten(n3d), out_rvblten(n3d), out_rthblten(n3d);
  std::vector<double> out_rqvblten(n3d), out_rqcblten(n3d), out_rqiblten(n3d);

  download_real(d_tke_myj, out_tke, n3d);
  download_real(d_exch_h, out_exch_h, n3d);
  download_real(d_el_myj, out_el_myj, n3d);
  download_real(d_rublten, out_rublten, n3d);
  download_real(d_rvblten, out_rvblten, n3d);
  download_real(d_rthblten, out_rthblten, n3d);
  download_real(d_rqvblten, out_rqvblten, n3d);
  download_real(d_rqcblten, out_rqcblten, n3d);
  download_real(d_rqiblten, out_rqiblten, n3d);

  FILE* fo = std::fopen(outfile.c_str(), "wb");
  if (!fo) {
    std::fprintf(stderr, "Cannot open output file: %s\n", outfile.c_str());
    return 1;
  }

  auto write_2d = [&](const std::vector<double>& arr) {
    std::fwrite(arr.data(), sizeof(double), n2d, fo);
  };
  auto write_3d = [&](const std::vector<double>& arr) {
    std::fwrite(arr.data(), sizeof(double), n3d, fo);
  };

  write_2d(h_qsfc);
  write_2d(h_thz0);
  write_2d(h_qz0);
  write_2d(h_ct);
  write_2d(h_pblh);
  write_2d(h_kpbl);
  write_2d(h_mixht);

  write_3d(out_tke);
  write_3d(out_exch_h);
  write_3d(out_el_myj);
  write_3d(out_rublten);
  write_3d(out_rvblten);
  write_3d(out_rthblten);
  write_3d(out_rqvblten);
  write_3d(out_rqcblten);
  write_3d(out_rqiblten);

  std::fclose(fo);

  dev_free(d_lowlyr);
  dev_free(d_kpbl);

  return 0;
}
