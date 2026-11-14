// SPDX-License-Identifier: GPL-3.0-or-later
// BMJ (cu_physics = 2) column test, the C++ side: cu_bmj::init_bmj_tables and cu_bmj::bmjdrv_column on the columns of
// generate_inputs.py, the output as driver_bmj.f90 writes it (FP64: compared bit for bit; FP32: a smoke run).
//   driver_bmj_cpu IN OUT
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>

#include "phys/cu_bmj/bmj.hpp"
#include "phys/cu_bmj/bmj_tables.hpp"
#include "phys/driver/wrf_constants.hpp"

using namespace cppwrf;
using cppwrf::real;

namespace {
// the views of one row of columns (i = column, k = level; j unused)
struct Cols {
  int nz;
  std::vector<real>*t, *qv, *p, *pi, *rho, *dz8w, *p8w, *rthcuten, *rqvcuten, *ccldfra, *qcconv, *qiconv;
  std::vector<real>*xland, *cldefi, *raincv, *pratec, *cutop, *cubot, *convcld;
  std::vector<int>* lowlyr_;
  struct A3 {
    std::vector<real>* v;
    int n;
    real& operator()(int i, int, int k) const { return (*v)[std::size_t(i) * n + k]; }
  };
  struct A2 {
    std::vector<real>* v;
    real& operator()(int i, int) const { return (*v)[std::size_t(i)]; }
  };
  struct I2 {
    std::vector<int>* v;
    int operator()(int i, int) const { return (*v)[std::size_t(i)]; }
  };
  A3 t3(std::vector<real>* v) const { return A3{v, nz}; }
};
struct F {
  Cols::A3 t, qv, p, pi, rho, dz8w, p8w, rthcuten, rqvcuten, ccldfra, qcconv, qiconv;
  Cols::A2 xland, cldefi, raincv, pratec, cutop, cubot, convcld;
  Cols::I2 lowlyr;
};
}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s IN OUT\n", argv[0]);
    return 2;
  }
  std::FILE* in = std::fopen(argv[1], "rb");
  if (!in) return 1;
  int32_t h[4];
  double dt;
  if (std::fread(h, 4, 4, in) != 4 || std::fread(&dt, 8, 1, in) != 1) return 1;
  const int ncols = h[0], nz = h[1], stepcu = h[2];
  const bool fb = h[3] != 0;
  auto rd = [&](std::size_t n) {
    std::vector<double> d(n);
    if (std::fread(d.data(), 8, n, in) != n) std::exit(1);
    return std::vector<real>(d.begin(), d.end());
  };
  const std::size_t n3 = std::size_t(nz) * ncols, n3w = std::size_t(nz + 1) * ncols;
  std::vector<real> t = rd(n3), qv = rd(n3), p = rd(n3), pi = rd(n3), rho = rd(n3), dz8w = rd(n3), p8w = rd(n3w);
  std::vector<real> xland = rd(ncols), cldefi = rd(ncols);
  std::fclose(in);
  std::vector<real> rthcuten(n3), rqvcuten(n3), ccldfra(n3), qcconv(n3), qiconv(n3);
  std::vector<real> raincv(ncols), pratec(ncols), cutop(ncols), cubot(ncols), convcld(ncols);
  std::vector<int> lowlyr(ncols, 1);  // bmjinit
  using namespace phys::wrfc;
  auto tb = std::make_unique<phys::cu_bmj::BmjTables<real>>();
  phys::cu_bmj::init_bmj_tables<real>(*tb, cp, r_d);
  F f{{&t, nz}, {&qv, nz}, {&p, nz}, {&pi, nz}, {&rho, nz}, {&dz8w, nz}, {&p8w, nz + 1}, {&rthcuten, nz},
      {&rqvcuten, nz}, {&ccldfra, nz}, {&qcconv, nz}, {&qiconv, nz}, {&xland}, {&cldefi}, {&raincv}, {&pratec},
      {&cutop}, {&cubot}, {&convcld}, {&lowlyr}};
  for (int i = 0; i < ncols; ++i)
    phys::cu_bmj::bmjdrv_column<real>(f, *tb, i, 0, nz, real(dt), stepcu, cp, r_d, xlv, g, svpt0, ep_1, fb,
                                      ColumnSlots{});
  std::FILE* out = std::fopen(argv[2], "wb");
  if (!out) return 1;
  auto wr = [&](const real* a, std::size_t n) {
    for (std::size_t q = 0; q < n; ++q) {
      const double v = double(a[q]);
      std::fwrite(&v, 8, 1, out);
    }
  };
  wr(tb->qs0, phys::cu_bmj::JTB);
  wr(tb->sqs, phys::cu_bmj::JTB);
  wr(&tb->ptbl[0][0], std::size_t(phys::cu_bmj::JTB) * phys::cu_bmj::ITB);
  for (const auto* a : {&rthcuten, &rqvcuten, &ccldfra, &qcconv, &qiconv}) wr(a->data(), n3);
  for (const auto* a : {&raincv, &pratec, &cutop, &cubot, &cldefi, &convcld}) wr(a->data(), std::size_t(ncols));
  std::fclose(out);
  return 0;
}
