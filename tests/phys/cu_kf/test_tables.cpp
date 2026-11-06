// SPDX-License-Identifier: GPL-3.0-or-later
// Kain-Fritsch lookup tables: init_kf_tables (kf_tables.cpp) against WRF's kf_lutab, whose tables
// dump_tables_f.f90 writes from module_cu_kfeta.o (tables_fortran.bin: rdpr, rdthk, plutop, the0k, alu, ttab,
// qstab). Every value must have the same bits (tolerance 0; -0 differs from +0).
//   test_tables [tables_fortran.bin]
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "src/phys/cu_kf/kf_tables.hpp"

using namespace cppwrf::phys::cu_kf;

int main(int argc, char** argv) {
  const std::string path = argc > 1 ? argv[1] : "tables_fortran.bin";
  KFTables<double> t{};
  init_kf_tables(t);

  std::ifstream in(path, std::ios::binary);
  if (!in) {
    std::cerr << "Failed to open " << path << "\n";
    return 1;
  }
  const std::vector<char> f((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

  // the same values in the order of the Fortran file (ttab / qstab: Fortran (KFNT, KFNP), C++ [KFNP][KFNT])
  struct Part {
    const char* name;
    const double* p;
    std::size_t n;
  };
  const Part parts[] = {{"rdpr", &t.rdpr, 1},
                        {"rdthk", &t.rdthk, 1},
                        {"plutop", &t.plutop, 1},
                        {"the0k", t.the0k, KFNP},
                        {"alu", t.alu, KFNALU},
                        {"ttab", &t.ttab[0][0], std::size_t(KFNP) * KFNT},
                        {"qstab", &t.qstab[0][0], std::size_t(KFNP) * KFNT}};
  std::size_t total = 0;
  for (const Part& q : parts) total += q.n * sizeof(double);
  if (f.size() != total) {
    std::cerr << "FAIL: " << path << " holds " << f.size() << " bytes, the tables " << total << "\n";
    return 1;
  }
  bool pass = true;
  std::size_t off = 0;
  for (const Part& q : parts) {
    std::size_t bad = 0, first = 0;
    for (std::size_t i = 0; i < q.n; ++i) {
      if (std::memcmp(f.data() + off + i * sizeof(double), q.p + i, sizeof(double)) != 0 && bad++ == 0) first = i;
    }
    if (bad) {
      double w;
      std::memcpy(&w, f.data() + off + first * sizeof(double), sizeof(double));
      std::printf("FAIL %s: %zu/%zu differ (first at %zu: C++ %.17g WRF %.17g)\n", q.name, bad, q.n, first,
                  q.p[first], w);
      pass = false;
    } else {
      std::printf("%s: %zu values identical\n", q.name, q.n);
    }
    off += q.n * sizeof(double);
  }
  std::puts(pass ? "KF tables: PASS (bit-identical)" : "KF tables: FAIL");
  return pass ? 0 : 1;
}
