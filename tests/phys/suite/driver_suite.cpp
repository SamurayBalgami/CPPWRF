// C++ side of the P1 physics-suite differential test (tests/phys/suite).
//
// Reads suite_ref.bin (written by fortran/suite_ref.F90 from the ORIGINAL WRF v4.6.1 code):
//   IN1:*  the state after WRF's start_domain (phy_init) = the suite's initial surface state
//   INn:*  the atmospheric state before call n (the dynamics in between is a fixed change)
// and runs the C++ PhysicsSuite (part1, part2, after_rk) on the same inputs, three calls with
// the surface state carried over by the C++ suite itself. Writes suite_cpp.bin with the same
// record names (P1_n after part1, P2_n after part2, OUTn after after_rk) for compare.py.
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/field.hpp"
#include "core/grid.hpp"
#include "dyn/state.hpp"
#include "phys/driver/suite.hpp"

using namespace cppwrf;
using phys::driver::IntField2;
using phys::driver::PhysicsSuite;
using phys::driver::SuiteOptions;

struct Rec {
  int nd = 0, n1 = 1, n2 = 1, n3 = 1;
  std::vector<double> a;
  double at(int i, int k, int j) const { return a[i + std::size_t(n1) * (k + std::size_t(n2) * j)]; }
  double at2(int i, int j) const { return a[i + std::size_t(n1) * j]; }
};

static std::map<std::string, Rec> read_records(const std::string& path) {
  std::ifstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("cannot open " + path);
  std::map<std::string, Rec> m;
  char name[33];
  while (f.read(name, 32)) {
    name[32] = 0;
    std::string n(name);
    n = n.substr(0, n.find_last_not_of(' ') + 1);
    Rec r;
    int hdr[4];
    f.read(reinterpret_cast<char*>(hdr), 16);
    r.nd = hdr[0]; r.n1 = hdr[1]; r.n2 = hdr[2]; r.n3 = hdr[3];
    r.a.resize(std::size_t(r.n1) * r.n2 * r.n3);
    f.read(reinterpret_cast<char*>(r.a.data()), r.a.size() * 8);
    m[n] = std::move(r);
  }
  return m;
}

class Writer {
 public:
  explicit Writer(const std::string& p) : f_(p, std::ios::binary) {}
  void rec(const std::string& name, int nd, int n1, int n2, int n3, const std::vector<double>& a) {
    char buf[32];
    std::memset(buf, ' ', 32);
    std::memcpy(buf, name.data(), name.size() < 32 ? name.size() : 32);
    f_.write(buf, 32);
    int hdr[4] = {nd, n1, n2, n3};
    f_.write(reinterpret_cast<const char*>(hdr), 16);
    f_.write(reinterpret_cast<const char*>(a.data()), a.size() * 8);
  }

 private:
  std::ofstream f_;
};

// Field (k = 0..nk-1 at the given stagger extents) <-> record in WRF (i,k,j) order
static void put3(Field3& f, const Rec& r) {
  f.fill_host([&](int i, int j, int k) { return (i < r.n1 && k < r.n2 && j < r.n3) ? r.at(i, k, j) : 0.0; });
}
static void put2(Field2& f, const Rec& r) {
  f.fill_host([&](int i, int j, int) { return (i < r.n1 && j < r.n2) ? r.at2(i, j) : 0.0; });
}
static std::vector<double> get3(const Field3& f, int n1, int n2, int n3) {
  std::vector<real> h;
  f.download(h);
  std::vector<double> a(std::size_t(n1) * n2 * n3);
  for (int j = 0; j < n3; ++j)
    for (int k = 0; k < n2; ++k)
      for (int i = 0; i < n1; ++i) a[i + std::size_t(n1) * (k + std::size_t(n2) * j)] = h[f.layout().index(i, j, k)];
  return a;
}
static std::vector<double> get2(const Field2& f, int n1, int n2) {
  std::vector<real> h;
  f.download(h);
  std::vector<double> a(std::size_t(n1) * n2);
  for (int j = 0; j < n2; ++j)
    for (int i = 0; i < n1; ++i) a[i + std::size_t(n1) * j] = h[f.layout().index(i, j, 0)];
  return a;
}
static void put_int(IntField2& f, const Rec& r) {
  std::vector<int> h(f.layout().total, 0);
  for (int j = 0; j < r.n2; ++j)
    for (int i = 0; i < r.n1; ++i) h[f.layout().index(i, j, 0)] = int(r.at2(i, j) + (r.at2(i, j) >= 0 ? 0.5 : -0.5));
  f.upload(h);
}
static std::vector<double> get_int(const IntField2& f, int n1, int n2) {
  std::vector<int> h;
  f.download(h);
  std::vector<double> a(std::size_t(n1) * n2);
  for (int j = 0; j < n2; ++j)
    for (int i = 0; i < n1; ++i) a[i + std::size_t(n1) * j] = h[f.layout().index(i, j, 0)];
  return a;
}

int main(int argc, char** argv) {
  const std::string ref = argc > 1 ? argv[1] : "suite_ref.bin";
  const std::string out = argc > 2 ? argv[2] : "suite_cpp.bin";
  const std::string tables = argc > 3 ? argv[3] : ".";
  auto R = read_records(ref);
  auto rec = [&](const std::string& n) -> const Rec& {
    auto it = R.find(n);
    if (it == R.end()) throw std::runtime_error("missing record " + n);
    return it->second;
  };
  const Rec& t1 = rec("IN1:T");
  Grid g;
  g.nx = t1.n1; g.nz = t1.n2; g.ny = t1.n3;
  g.dx = g.dy = rec("IN1:DX").a[0];
  g.dt = rec("IN1:DT").a[0];
  g.p_top = rec("IN1:P_TOP").a[0];
  g.bc_x = g.bc_y = BC::Specified;
  const int nx = g.nx, ny = g.ny, nz = g.nz;
  const int nsoil = rec("IN1:TSLB").n2;

  dyn::State s(g);
  dyn::Tendencies t(g);
  s.add_scalar(g, "qv", true);
  t.add_scalar(g, "qv");
  dyn::BaseState b(g, /*flat=*/false, /*mapped=*/false);
  b.p_top = g.p_top;
  {
    const Rec &c1 = rec("IN1:C1H"), &c2 = rec("IN1:C2H"), &fm = rec("IN1:FNM"), &fp = rec("IN1:FNP"),
              &dn = rec("IN1:DNW"), &zw = rec("IN1:ZNW");
    for (int k = 0; k < nz; ++k) {
      b.c1h[k] = c1.a[k]; b.c2h[k] = c2.a[k]; b.fnm[k] = fm.a[k]; b.fnp[k] = fp.a[k]; b.dnw[k] = dn.a[k];
    }
    for (int k = 0; k <= nz; ++k) b.znw[k] = zw.a[k];
    b.upload_levels();
  }

  SuiteOptions o;
  o.table_dir = tables;
  o.num_soil_layers = nsoil;
  PhysicsSuite su(g, o);
  for (int k = 0; k < nsoil; ++k) su.dzs[k] = rec("IN1:DZS").a[k];

  // ---- initial surface state: WRF after start_domain ----
  std::map<std::string, Field2*> f2 = {
      {"TSK", &su.tsk}, {"HFX", &su.hfx}, {"QFX", &su.qfx}, {"LH", &su.lh}, {"UST", &su.ust}, {"ZNT", &su.znt},
      {"Z0", &su.z0}, {"MOL", &su.mol}, {"RMOL", &su.rmol}, {"PBLH", &su.pblh}, {"QSFC", &su.qsfc},
      {"FLHC", &su.flhc}, {"FLQC", &su.flqc}, {"T2", &su.t2}, {"TH2", &su.th2}, {"Q2", &su.q2}, {"U10", &su.u10},
      {"V10", &su.v10}, {"BR", &su.br}, {"WSPD", &su.wspd}, {"FM", &su.fm}, {"FH", &su.fh}, {"ZOL", &su.zol},
      {"REGIME", &su.regime}, {"PSFC", &su.psfc}, {"GRDFLX", &su.grdflx}, {"ALBEDO", &su.albedo},
      {"ALBBCK", &su.albbck}, {"EMISS", &su.emiss}, {"EMBCK", &su.embck}, {"SNOW", &su.snow}, {"SNOWH", &su.snowh},
      {"SNOWC", &su.snowc}, {"CANWAT", &su.canwat}, {"SMSTAV", &su.smstav}, {"SMSTOT", &su.smstot},
      {"SFCRUNOFF", &su.sfcrunoff}, {"UDRUNOFF", &su.udrunoff}, {"ACSNOM", &su.acsnom}, {"ACSNOW", &su.acsnow},
      {"SNOPCX", &su.snopcx}, {"POTEVP", &su.potevp}, {"NOAHRES", &su.noahres}, {"SFCEVP", &su.sfcevp},
      {"SFCEXC", &su.sfcexc}, {"ACHFX", &su.achfx}, {"ACLHF", &su.aclhf}, {"ACGRDFLX", &su.acgrdflx},
      {"LAI", &su.lai}, {"SNOTIME", &su.snotime}, {"QZ0", &su.qz0}, {"RAINBL", &su.rainbl}, {"XLAND", &su.xland},
      {"XICE", &su.xice}, {"MAVAIL", &su.mavail}, {"TMN", &su.tmn}, {"VEGFRA", &su.vegfra},
      {"SHDMAX", &su.shdmax}, {"SHDMIN", &su.shdmin}, {"SNOALB", &su.snoalb}, {"LAKEMASK", &su.lakemask},
      {"DX2D", &su.dx2d}, {"UOCE", &su.uoce}, {"VOCE", &su.voce}, {"CTOPO", &su.ctopo}, {"CTOPO2", &su.ctopo2},
      {"WATER_DEPTH", &su.water_depth}, {"WSTAR_YSU", &su.wstar}, {"DELTA_YSU", &su.delta}, {"USTM", &su.ustm}, {"CK", &su.ck}, {"CKA", &su.cka},
      {"CD", &su.cd}, {"CDA", &su.cda}, {"SR", &su.sr}};
  std::map<std::string, Field3*> fsoil = {
      {"TSLB", &su.tslb}, {"SMOIS", &su.smois}, {"SH2O", &su.sh2o}, {"SMCREL", &su.smcrel}};
  std::map<std::string, Field3*> fm3 = {{"EXCH_H", &su.exch_h},       {"EXCH_M", &su.exch_m},
                                        {"RUBLTEN", &su.pt.rublten},  {"RVBLTEN", &su.pt.rvblten},
                                        {"RTHBLTEN", &su.pt.rthblten}, {"RQVBLTEN", &su.pt.rqvblten},
                                        {"RTHRATEN", &su.pt.rthraten}};
  std::map<std::string, IntField2*> fi = {{"KPBL", &su.kpbl}, {"IVGTYP", &su.ivgtyp}, {"ISLTYP", &su.isltyp}};
  for (auto& [n, f] : f2) put2(*f, rec("IN1:" + n));
  for (auto& [n, f] : fsoil) put3(*f, rec("IN1:" + n));
  for (auto& [n, f] : fm3) put3(*f, rec("IN1:" + n));
  for (auto& [n, f] : fi) put_int(*f, rec("IN1:" + n));

  Writer w(out);
  auto dump_surface = [&](const std::string& tag) {
    for (auto& [n, f] : f2) w.rec(tag + n, 2, nx, ny, 1, get2(*f, nx, ny));
    for (auto& [n, f] : fsoil) w.rec(tag + n, 3, nx, nsoil, ny, get3(*f, nx, nsoil, ny));
    for (auto& [n, f] : fm3) w.rec(tag + n, 3, nx, nz, ny, get3(*f, nx, nz, ny));
    for (auto& [n, f] : fi) w.rec(tag + n, 2, nx, ny, 1, get_int(*f, nx, ny));
  };

  for (int n = 1; n <= 3; ++n) {
    const std::string in = "IN" + std::to_string(n) + ":";
    // atmospheric state and step inputs before the call
    put3(s.u_2, rec(in + "U")); put3(s.v_2, rec(in + "V")); put3(s.t_2, rec(in + "T")); put3(s.t_1, rec(in + "T_1"));
    put3(s.p, rec(in + "P")); put3(s.alt, rec(in + "ALT")); put3(s.ph_2, rec(in + "PH"));
    put3(s.q_2[s.P_QV], rec(in + "QVAPOR")); put3(s.q_1[s.P_QV], rec(in + "QV_OLD"));
    put2(s.mut, rec(in + "MUT")); put2(s.muts, rec(in + "MUTS"));
    {
      const Rec &pb = rec(in + "PB"), &phb = rec(in + "PHB");
      b.pb.fill_host([&](int i, int j, int k) { return pb.at(i, k, j); });
      b.phb.fill_host([&](int i, int j, int k) { return phb.at(i, k, j); });
    }
    put2(su.rainncv, rec(in + "RAINNCV")); put2(su.raincv, rec(in + "RAINCV"));
    put2(su.glw, rec(in + "GLW")); put2(su.swdown, rec(in + "SWDOWN")); put2(su.gsw, rec(in + "GSW"));
    // init_zero_tendency (the dycore does it before ColumnPhysics::part1)
    for (Field3* f : {&t.ru_tendf, &t.rv_tendf, &t.rw_tendf, &t.t_tendf, &t.ph_tendf, &t.q_tendf[0]})
      f->fill_host([](int, int, int) { return 0.0; });

    su.part1(s, b, t);
    const std::string p1 = "P1_" + std::to_string(n) + ":";
    w.rec(p1 + "TH_PHY", 3, nx, nz, ny, get3(su.prep.th_phy, nx, nz, ny));
    w.rec(p1 + "P_PHY", 3, nx, nz, ny, get3(su.prep.p_phy, nx, nz, ny));
    w.rec(p1 + "PI_PHY", 3, nx, nz, ny, get3(su.prep.pi_phy, nx, nz, ny));
    w.rec(p1 + "T_PHY", 3, nx, nz, ny, get3(su.prep.t_phy, nx, nz, ny));
    w.rec(p1 + "RHO", 3, nx, nz, ny, get3(su.prep.rho, nx, nz, ny));
    w.rec(p1 + "U_PHY", 3, nx, nz, ny, get3(su.prep.u_phy, nx, nz, ny));
    w.rec(p1 + "V_PHY", 3, nx, nz, ny, get3(su.prep.v_phy, nx, nz, ny));
    w.rec(p1 + "Z", 3, nx, nz, ny, get3(su.prep.z, nx, nz, ny));
    w.rec(p1 + "Z_AT_W", 3, nx, nz + 1, ny, get3(su.prep.z_at_w, nx, nz + 1, ny));
    w.rec(p1 + "DZ8W", 3, nx, nz + 1, ny, get3(su.prep.dz8w, nx, nz + 1, ny));
    w.rec(p1 + "P8W", 3, nx, nz + 1, ny, get3(su.prep.p8w, nx, nz + 1, ny));
    w.rec(p1 + "T8W", 3, nx, nz + 1, ny, get3(su.prep.t8w, nx, nz + 1, ny));
    w.rec(p1 + "P_HYD", 3, nx, nz, ny, get3(su.prep.p_hyd, nx, nz, ny));
    w.rec(p1 + "P_HYD_W", 3, nx, nz + 1, ny, get3(su.prep.p_hyd_w, nx, nz + 1, ny));
    w.rec(p1 + "PSIM", 2, nx, ny, 1, get2(su.psim, nx, ny));
    w.rec(p1 + "PSIH", 2, nx, ny, 1, get2(su.psih, nx, ny));
    w.rec(p1 + "GZ1OZ0", 2, nx, ny, 1, get2(su.gz1oz0, nx, ny));
    w.rec(p1 + "CHKLOWQ", 2, nx, ny, 1, get2(su.chklowq, nx, ny));
    dump_surface(p1);

    su.part2(s, b, t);
    const std::string p2 = "P2_" + std::to_string(n) + ":";
    w.rec(p2 + "RU_TENDF", 3, nx + 1, nz, ny, get3(t.ru_tendf, nx + 1, nz, ny));
    w.rec(p2 + "RV_TENDF", 3, nx, nz, ny + 1, get3(t.rv_tendf, nx, nz, ny + 1));
    w.rec(p2 + "T_TENDF", 3, nx, nz, ny, get3(t.t_tendf, nx, nz, ny));
    w.rec(p2 + "QV_TENDF", 3, nx, nz, ny, get3(t.q_tendf[0], nx, nz, ny));
    w.rec(p2 + "RUBLTEN", 3, nx, nz, ny, get3(su.pt.rublten, nx, nz, ny));
    w.rec(p2 + "RVBLTEN", 3, nx, nz, ny, get3(su.pt.rvblten, nx, nz, ny));
    w.rec(p2 + "RTHBLTEN", 3, nx, nz, ny, get3(su.pt.rthblten, nx, nz, ny));
    w.rec(p2 + "RQVBLTEN", 3, nx, nz, ny, get3(su.pt.rqvblten, nx, nz, ny));

    su.after_rk(s, b);
    dump_surface("OUT" + std::to_string(n) + ":");
  }
  std::printf("driver_suite: 3 calls done (%s, real = %s), wrote %s\n", backend_name(), precision_name, out.c_str());
  return 0;
}
