// cppwrf: ./cppwrf <config.ini> [--set section.key=value ...]
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

#include "cases/domain_runner.hpp"
#include "cases/model.hpp"
#include "core/backend.hpp"
#include "io/config.hpp"

using namespace cppwrf;

static int run(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <config.ini> [--set section.key=value ...]\n", argv[0]);
    return 2;
  }
  std::vector<std::string> ov;
  for (int a = 2; a < argc; ++a) {
    std::string s = argv[a];
    if (s == "--set" && a + 1 < argc) ov.push_back(argv[++a]);
    else if (s.rfind("--set=", 0) == 0) ov.push_back(s.substr(6));
    else throw std::runtime_error("unknown argument: " + s);
  }
  io::Config cfg = io::Config::load(argv[1], ov);
  const std::string cname = cfg.get_string("run.case");
  const std::string out_dir = cfg.get_string("output.dir", "out/" + cname);

  cases::Model m = cases::Model::build(cfg);
  cases::DomainRunner r(cfg, m, out_dir, "d01");
  if (!r.start()) return 3;

  for (long n = r.start_step + 1; n <= m.nsteps; ++n) {
    r.before_step(n);
    m.rk->step(*m.st, *m.tend, *m.base);
    dev_sync();
    if (!r.after_step(n)) return 3;
  }
  r.finish();
  return 0;
}

int main(int argc, char** argv) {
  try {
    return run(argc, argv);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "cppwrf: error: %s\n", e.what());
    return 2;
  }
}
