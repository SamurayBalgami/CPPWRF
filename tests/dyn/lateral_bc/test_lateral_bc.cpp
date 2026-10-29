// Differential test of src/dyn/lateral_bc.cpp (LateralBC, specified boundaries) against WRF v4.8.0's
// own routines (lbc_bridge.f90): every LateralBC call against the WRF calls solve_em makes at that
// point, from the same random fields and boundary data, compared bit for bit (-0 differs from +0).
// The grid is not square (x/y mix-ups show) and the map factors vary in the horizontal.
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "core/field.hpp"
#include "core/grid.hpp"
#include "dyn/lateral_bc.hpp"
#include "dyn/state.hpp"

using namespace cppwrf;
using namespace cppwrf::dyn;

extern "C" {
void lbc_set_nested_c(int nested);
void lbc_fcx_gcx_c(int nx, int ny, int nz, double dt, double* fcx, double* gcx);
void relax_bdy_dry_c(int nx, int ny, int nz, double* u_save, double* v_save, double* ph_save, double* t_save,
                     double* w_save, double* mu_tend, const double* c1h, const double* c2h, const double* c1f,
                     const double* c2f, const double* ru, const double* rv, const double* ph, const double* t,
                     const double* w, const double* mu, const double* mut, const double* ub, const double* vb,
                     const double* phb, const double* tb, const double* wb, const double* mub, const double* fcx,
                     const double* gcx, double dtbc);
void spec_bdy_dry_c(int nx, int ny, int nz, double* ru_tend, double* rv_tend, double* ph_tend, double* t_tend,
                    double* rw_tend, double* mu_tend, const double* ub, const double* vb, const double* phb,
                    const double* tb, const double* wb, const double* mub);
void spec_bdyupdate_c(int nx, int ny, int nz, double* field, const double* tend, double dts, char var, int two_d);
void update_ph_w_c(int nx, int ny, int nz, const double* ph_save, double* ph, const double* ph_tend,
                   const double* mu_tend, const double* muts, const double* c1f, const double* c2f, double dts,
                   double* w, const double* rw_tend);
void relax_spec_scalar_c(int nx, int ny, int nz, double* q_tend, const double* q, const double* mut,
                         const double* c1h, const double* c2h, const double* qb, const double* fcx, const double* gcx,
                         double dtbc, int relax, int spec);
void flow_dep_bdy_c(int nx, int ny, int nz, double* q, const double* ru_m, const double* rv_m);
void spec_bdy_final_c(int nx, int ny, int nz, double* field, const double* mu, const double* c1, const double* c2,
                      const double* msf, const double* fb, char var, double dtbc, int two_d);
}

namespace {

constexpr int NX = 13, NY = 11, NZ = 7, H = 3, WIDTH = 5;
constexpr int NI = NX + 1 + 2 * H, NJ = NY + 1 + 2 * H, NK = NZ + 1;  // WRF memory dims (ims:ime, kms:kme, jms:jme)
constexpr double DT = 60.0, DTBC = 1234.5, DTS = 10.0;

// WRF memory index of C++ point (i, j, k): Fortran (i+1, k+1, j+1), ims = jms = 1-H, kms = 1
std::size_t f3(int i, int j, int k) { return std::size_t(i + H) + NI * (std::size_t(k) + NK * std::size_t(j + H)); }
std::size_t f2(int i, int j) { return std::size_t(i + H) + NI * std::size_t(j + H); }

std::uint64_t rng = 0x9E3779B97F4A7C15ull;
double rnd(double lo, double hi) {
  rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
  return lo + (hi - lo) * double(rng >> 11) / double(1ull << 53);
}

void fill(Field3& f, double lo, double hi) { f.fill_host([&](int, int, int) { return rnd(lo, hi); }); }
void fill(Field2& f, double lo, double hi) { f.fill_host([&](int, int, int) { return rnd(lo, hi); }); }

// WRF memory array of a field (halo points CPPWRF does not store are 0, as CPPWRF's own halo)
template <int R> std::vector<double> wrf(const FieldT<R>& f) {
  std::vector<real> h;
  f.download(h);
  const Layout& L = f.layout();
  std::vector<double> a(R == 3 ? std::size_t(NI) * NK * NJ : std::size_t(NI) * NJ, 0.0);
  for (int j = -H; j <= NY + H; ++j)
    for (int k = 0; k < (R == 3 ? NK : 1); ++k)
      for (int i = -H; i <= NX + H; ++i) {
        if (k >= L.nz || j < -L.hy || j >= L.ny + L.hy || i >= L.nx + H) continue;
        a[R == 3 ? f3(i, j, k) : f2(i, j)] = h[L.index(i, j, k)];
      }
  return a;
}
std::vector<double> wrf(const BaseField2& f) {
  std::vector<real> h;
  f.download(h);
  const Layout& L = f.layout();
  std::vector<double> a(std::size_t(NI) * NJ, 0.0);
  for (int j = -H; j <= NY + H; ++j)
    for (int i = -H; i <= NX + H; ++i) {
      if (L.si == 0) { a[f2(i, j)] = h[0]; continue; }
      if (j < -L.hy || j >= L.ny + L.hy || i >= L.nx + H) continue;
      a[f2(i, j)] = h[L.index(i, j, 0)];
    }
  return a;
}
std::vector<double> wrf(const Array1& c) {  // (kms:kme)
  std::vector<double> a(NK, 0.0);
  for (int k = 0; k < c.size(); ++k) a[k] = c[k];
  return a;
}

// Boundary data of one variable: CPPWRF slabs (pos, k, b) per side, and WRF's packed arrays
// xs, xe, ys, ye, btxs, btxe, btys, btye with k extent kext (nz + 1, or 1 for mu)
struct Bdy {
  std::vector<double> v[4], t[4];
  std::vector<double> packed;
};
Bdy make_bdy(int nposx, int nposy, int nk, int kext, double lo, double hi) {
  Bdy d;
  const std::size_t sx = std::size_t(NJ) * kext * WIDTH, sy = std::size_t(NI) * kext * WIDTH;
  d.packed.assign(4 * sx + 4 * sy, 0.0);
  std::size_t off = 0;
  for (int half = 0; half < 2; ++half)
    for (int s = 0; s < 4; ++s) {
      const int npos = s < 2 ? nposx : nposy, nm = s < 2 ? NJ : NI;
      std::vector<double>& c = half == 0 ? d.v[s] : d.t[s];
      c.resize(std::size_t(npos) * nk * WIDTH);
      for (int b = 0; b < WIDTH; ++b)
        for (int k = 0; k < nk; ++k)
          for (int p = 0; p < npos; ++p) {
            const double x = half == 0 ? rnd(lo, hi) : rnd(lo, hi) * 1e-3;
            c[p + std::size_t(k) * npos + std::size_t(b) * npos * nk] = x;
            d.packed[off + (p + H) + std::size_t(nm) * (k + std::size_t(kext) * b)] = x;
          }
      off += s < 2 ? sx : sy;
    }
  return d;
}

int failures = 0;
bool same(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }
// CPPWRF's field over its logical points against WRF's array
template <int R> void check(const char* what, const FieldT<R>& f, const std::vector<double>& ref) {
  std::vector<real> h;
  f.download(h);
  const Layout& L = f.layout();
  long bad = 0, n = 0;
  for (int k = 0; k < L.nz; ++k)
    for (int j = 0; j < L.ny; ++j)
      for (int i = 0; i < L.nx; ++i, ++n) {
        const double c = h[L.index(i, j, k)], w = ref[R == 3 ? f3(i, j, k) : f2(i, j)];
        if (same(c, w)) continue;
        if (++bad <= 3) std::printf("    %s(i=%d, j=%d, k=%d): CPPWRF %a  WRF %a\n", what, i, j, k, c, w);
      }
  std::printf("  %-34s %6ld points  %s\n", what, n, bad ? (std::to_string(bad) + " differ  FAIL").c_str() : "identical");
  if (bad) ++failures;
}

// One pass over every LateralBC call: a real-data domain (specified) or a nest (nested: w relaxed,
// specified, updated and finalised from boundary data too)
void run(bool nested) {
  Grid g;
  g.nx = NX; g.ny = NY; g.nz = NZ; g.dt = DT;
  g.bc_x = g.bc_y = BC::Specified;
  State s(g);
  Tendencies t(g);
  BaseState b(g, true, true);  // map factors vary in the horizontal
  const int qv = s.add_scalar(g, "QVAPOR", true);
  t.add_scalar(g, "QVAPOR");

  auto msf = [&](BaseField2& m) { m.fill_host_halo([&](int, int, int) { return rnd(0.9, 1.1); }); };
  msf(b.msftx); msf(b.msfty); msf(b.msfux); msf(b.msfuy); msf(b.msfvx); msf(b.msfvy);
  b.c1h = Array1(NZ); b.c2h = Array1(NZ); b.c1f = Array1(NZ + 1); b.c2f = Array1(NZ + 1);
  for (int k = 0; k <= NZ; ++k) {
    if (k < NZ) { b.c1h[k] = rnd(0.3, 1.0); b.c2h[k] = rnd(0.0, 5000.0); }
    b.c1f[k] = rnd(0.3, 1.0); b.c2f[k] = rnd(0.0, 5000.0);
  }
  b.c1h.upload(); b.c2h.upload(); b.c1f.upload(); b.c2f.upload();

  lbc_set_nested_c(nested);
  LateralBC lbc(g, DT, WIDTH, 1, 4, nested);
  lbc.dtbc = DTBC;
  // boundary data as wrfbdy holds it: u, v coupled (mu * u / msf), t and ph coupled, mu
  Bdy bu = make_bdy(NY, NX + 1, NZ, NK, -2e6, 2e6), bv = make_bdy(NY + 1, NX, NZ, NK, -2e6, 2e6);
  Bdy bt = make_bdy(NY, NX, NZ, NK, 1e5, 2e6), bph = make_bdy(NY, NX, NZ + 1, NK, -5e7, 5e7);
  Bdy bmu = make_bdy(NY, NX, 1, 1, -2e3, 2e3), bq = make_bdy(NY, NX, NZ, NK, 0.0, 1e3);
  Bdy bwd = make_bdy(NY, NX, NZ + 1, NK, -5e6, 5e6);  // w: nests only (zero on a real-data domain)
  const std::vector<double> bw = nested ? bwd.packed : std::vector<double>(bwd.packed.size(), 0.0);
  const LateralBC::Var vars[] = {LateralBC::U, LateralBC::V, LateralBC::T, LateralBC::PH, LateralBC::MU, LateralBC::W};
  const Bdy* data[] = {&bu, &bv, &bt, &bph, &bmu, &bwd};
  for (int v = 0; v < (nested ? 6 : 5); ++v)
    for (int sd = 0; sd < 4; ++sd) lbc.set(vars[v], LateralBC::Side(sd), data[v]->v[sd], data[v]->t[sd]);
  for (int sd = 0; sd < 4; ++sd) lbc.set_scalar(qv, LateralBC::Side(sd), bq.v[sd], bq.t[sd]);
  std::vector<double> fcx(WIDTH), gcx(WIDTH);
  lbc_fcx_gcx_c(NX, NY, NZ, DT, fcx.data(), gcx.data());

  std::vector<double> c1h = wrf(b.c1h), c2h = wrf(b.c2h), c1f = wrf(b.c1f), c2f = wrf(b.c2f);
  std::vector<double> msfuy = wrf(b.msfuy), msfvx = wrf(b.msfvx), msfty = wrf(b.msfty);

  std::printf("LateralBC vs WRF v4.8.0, %s (nx %d, ny %d, nz %d, dtbc %g s)\n",
              nested ? "nest (nested)" : "real-data domain (specified)", NX, NY, NZ, DTBC);

  {  // relax_bdy_dry -> *_save, then rk_addtend_dry (rk 1): *_tendf + *_save (* msf for u, v, w); mu_tend direct
    fill(s.ru, -3e6, 3e6); fill(s.rv, -3e6, 3e6); fill(s.ph_2, -5e3, 5e3); fill(s.t_2, -20, 40); fill(s.w_2, -5, 5);
    fill(s.mu_2, -2e3, 2e3); fill(s.mut, 5e4, 9e4);
    fill(t.ru_tendf, -50, 50); fill(t.rv_tendf, -50, 50); fill(t.ph_tendf, -50, 50); fill(t.t_tendf, -50, 50);
    fill(t.rw_tendf, -50, 50);
    t.mu_tend.fill_host([](int, int, int) { return 0.0; });  // WRF zero_tend before relax_bdy_dry
    const std::vector<double> nul(std::size_t(NI) * NK * NJ, 0.0);
    std::vector<double> u_save = nul, v_save = nul, ph_save = nul, t_save = nul, w_save = nul, mu_tend = wrf(t.mu_tend);
    relax_bdy_dry_c(NX, NY, NZ, u_save.data(), v_save.data(), ph_save.data(), t_save.data(), w_save.data(),
                    mu_tend.data(), c1h.data(), c2h.data(), c1f.data(), c2f.data(), wrf(s.ru).data(), wrf(s.rv).data(),
                    wrf(s.ph_2).data(), wrf(s.t_2).data(), wrf(s.w_2).data(), wrf(s.mu_2).data(), wrf(s.mut).data(),
                    bu.packed.data(), bv.packed.data(), bph.packed.data(), bt.packed.data(), bw.data(),
                    bmu.packed.data(), fcx.data(), gcx.data(), DTBC);
    std::vector<double> ru = wrf(t.ru_tendf), rv = wrf(t.rv_tendf), ph = wrf(t.ph_tendf), tt = wrf(t.t_tendf),
                        rw = wrf(t.rw_tendf);
    for (int j = -H; j <= NY + H; ++j)
      for (int k = 0; k < NK; ++k)
        for (int i = -H; i <= NX + H; ++i) {
          const std::size_t n = f3(i, j, k), m = f2(i, j);
          ru[n] = ru[n] + u_save[n] * msfuy[m];
          rv[n] = rv[n] + v_save[n] * msfvx[m];
          rw[n] = rw[n] + w_save[n] * msfty[m];
          ph[n] = ph[n] + ph_save[n];
          tt[n] = tt[n] + t_save[n];
        }
    lbc.relax_dry(s, t, b);
    check("relax_dry ru_tendf", t.ru_tendf, ru);
    check("relax_dry rv_tendf", t.rv_tendf, rv);
    check("relax_dry rw_tendf", t.rw_tendf, rw);
    check("relax_dry ph_tendf", t.ph_tendf, ph);
    check("relax_dry t_tendf", t.t_tendf, tt);
    check("relax_dry mu_tend", t.mu_tend, mu_tend);
  }
  {  // spec_bdy_dry
    fill(t.ru_tend, -50, 50); fill(t.rv_tend, -50, 50); fill(t.ph_tend, -50, 50); fill(t.t_tend, -50, 50);
    fill(t.rw_tend, -50, 50); fill(t.mu_tend, -5, 5);
    std::vector<double> ru = wrf(t.ru_tend), rv = wrf(t.rv_tend), ph = wrf(t.ph_tend), tt = wrf(t.t_tend),
                        rw = wrf(t.rw_tend), mu = wrf(t.mu_tend);
    spec_bdy_dry_c(NX, NY, NZ, ru.data(), rv.data(), ph.data(), tt.data(), rw.data(), mu.data(), bu.packed.data(),
                   bv.packed.data(), bph.packed.data(), bt.packed.data(), bw.data(), bmu.packed.data());
    lbc.spec_dry(t);
    check("spec_dry ru_tend", t.ru_tend, ru);
    check("spec_dry rv_tend", t.rv_tend, rv);
    check("spec_dry rw_tend", t.rw_tend, rw);
    check("spec_dry ph_tend", t.ph_tend, ph);
    check("spec_dry t_tend", t.t_tend, tt);
    check("spec_dry mu_tend", t.mu_tend, mu);
  }
  {  // spec_bdyupdate u, v (after advance_uv)
    fill(s.u_2, -30, 30); fill(s.v_2, -30, 30); fill(t.ru_tend, -50, 50); fill(t.rv_tend, -50, 50);
    std::vector<double> u = wrf(s.u_2), v = wrf(s.v_2);
    spec_bdyupdate_c(NX, NY, NZ, u.data(), wrf(t.ru_tend).data(), DTS, 'u', 0);
    spec_bdyupdate_c(NX, NY, NZ, v.data(), wrf(t.rv_tend).data(), DTS, 'v', 0);
    lbc.update_uv(s, t, real(DTS));
    check("update_uv u_2", s.u_2, u);
    check("update_uv v_2", s.v_2, v);
  }
  {  // spec_bdyupdate t, mu_2, muts (after advance_mu_t)
    fill(s.t_2, -20, 40); fill(s.mu_2, -2e3, 2e3); fill(s.muts, 5e4, 9e4); fill(t.t_tend, -50, 50);
    fill(t.mu_tend, -5, 5);
    std::vector<double> tt = wrf(s.t_2), mu = wrf(s.mu_2), muts = wrf(s.muts), mt = wrf(t.mu_tend);
    spec_bdyupdate_c(NX, NY, NZ, tt.data(), wrf(t.t_tend).data(), DTS, 't', 0);
    spec_bdyupdate_c(NX, NY, NZ, mu.data(), mt.data(), DTS, 'm', 1);
    spec_bdyupdate_c(NX, NY, NZ, muts.data(), mt.data(), DTS, 'm', 1);
    lbc.update_mu_t(s, t, real(DTS));
    check("update_mu_t t_2", s.t_2, tt);
    check("update_mu_t mu_2", s.mu_2, mu);
    check("update_mu_t muts", s.muts, muts);
  }
  {  // spec_bdyupdate_ph, then zero_grad_bdy w or, in a nest, spec_bdyupdate w (after advance_w)
    Field3 ph_save(g, Stagger::W, "ph_save");
    fill(ph_save, -5e3, 5e3); fill(s.ph_2, -5e3, 5e3); fill(t.ph_tend, -50, 50); fill(t.mu_tend, -5, 5);
    fill(s.muts, 5e4, 9e4); fill(s.w_2, -5, 5); fill(t.rw_tend, -50, 50);
    std::vector<double> ph = wrf(s.ph_2), w = wrf(s.w_2);
    update_ph_w_c(NX, NY, NZ, wrf(ph_save).data(), ph.data(), wrf(t.ph_tend).data(), wrf(t.mu_tend).data(),
                  wrf(s.muts).data(), c1f.data(), c2f.data(), DTS, w.data(), wrf(t.rw_tend).data());
    lbc.update_ph_w(s, ph_save, t, b, real(DTS));
    check("update_ph_w ph_2", s.ph_2, ph);
    check("update_ph_w w_2", s.w_2, w);
  }
  {  // relax_bdy_scalar, spec_bdy_scalar into moist_tend (which may hold physics tendencies)
    fill(s.q_2[qv], 0.0, 0.02); fill(t.q_tendf[qv], -1e-3, 1e-3); fill(s.mut, 5e4, 9e4);
    std::vector<double> qt = wrf(t.q_tendf[qv]);
    relax_spec_scalar_c(NX, NY, NZ, qt.data(), wrf(s.q_2[qv]).data(), wrf(s.mut).data(), c1h.data(), c2h.data(),
                        bq.packed.data(), fcx.data(), gcx.data(), DTBC, 1, 1);
    lbc.relax_scalar(qv, t.q_tendf[qv], s.q_2[qv], s, b);
    lbc.spec_scalar(qv, t.q_tendf[qv]);
    check("relax_scalar + spec_scalar", t.q_tendf[qv], qt);
  }
  {  // flow_dep_bdy (a moist species without boundary data)
    fill(s.q_2[qv], 0.0, 0.02); fill(s.ru_m, -5e5, 5e5); fill(s.rv_m, -5e5, 5e5);
    std::vector<double> q = wrf(s.q_2[qv]);
    flow_dep_bdy_c(NX, NY, NZ, q.data(), wrf(s.ru_m).data(), wrf(s.rv_m).data());
    lbc.flow_dep(s.q_2[qv], s);
    check("flow_dep q", s.q_2[qv], q);
  }
  {  // spec_bdy_final u, v, (w: nests), t, ph, mu, moist (end of the step)
    fill(s.u_2, -30, 30); fill(s.v_2, -30, 30); fill(s.w_2, -5, 5); fill(s.t_2, -20, 40); fill(s.ph_2, -5e3, 5e3);
    fill(s.mu_2, -2e3, 2e3); fill(s.q_2[qv], 0.0, 0.02);
    fill(s.muus, 5e4, 9e4); fill(s.muvs, 5e4, 9e4); fill(s.muts, 5e4, 9e4);
    std::vector<double> u = wrf(s.u_2), v = wrf(s.v_2), w = wrf(s.w_2), tt = wrf(s.t_2), ph = wrf(s.ph_2),
                        mu = wrf(s.mu_2), q = wrf(s.q_2[qv]), muts = wrf(s.muts);
    spec_bdy_final_c(NX, NY, NZ, u.data(), wrf(s.muus).data(), c1h.data(), c2h.data(), msfuy.data(),
                     bu.packed.data(), 'u', DTBC, 0);
    spec_bdy_final_c(NX, NY, NZ, v.data(), wrf(s.muvs).data(), c1h.data(), c2h.data(), msfvx.data(),
                     bv.packed.data(), 'v', DTBC, 0);
    if (nested)
      spec_bdy_final_c(NX, NY, NZ, w.data(), muts.data(), c1f.data(), c2f.data(), msfty.data(), bw.data(), 'w',
                       DTBC, 0);
    spec_bdy_final_c(NX, NY, NZ, tt.data(), muts.data(), c1h.data(), c2h.data(), msfty.data(), bt.packed.data(), 't',
                     DTBC, 0);
    spec_bdy_final_c(NX, NY, NZ, ph.data(), muts.data(), c1f.data(), c2f.data(), msfty.data(), bph.packed.data(),
                     'h', DTBC, 0);
    spec_bdy_final_c(NX, NY, NZ, mu.data(), muts.data(), c1h.data(), c2h.data(), msfty.data(), bmu.packed.data(),
                     'm', DTBC, 1);
    spec_bdy_final_c(NX, NY, NZ, q.data(), muts.data(), c1h.data(), c2h.data(), msfty.data(), bq.packed.data(), 't',
                     DTBC, 0);
    lbc.final(s, b);
    check("final u_2", s.u_2, u);
    check("final v_2", s.v_2, v);
    check("final w_2", s.w_2, w);
    check("final t_2", s.t_2, tt);
    check("final ph_2", s.ph_2, ph);
    check("final mu_2", s.mu_2, mu);
    check("final moist", s.q_2[qv], q);
  }
}

}  // namespace

int main() {
  run(false);
  run(true);
  std::printf(failures ? "LATERAL BC: %d checks FAILED\n" : "LATERAL BC: all checks identical to WRF\n", failures);
  return failures ? 1 : 0;
}
