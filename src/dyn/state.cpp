#include "dyn/state.hpp"

#include <stdexcept>

namespace cppwrf::dyn {

using S = Stagger;

State::State(const Grid& g)
    : u_1(g, S::U, "u_1"), u_2(g, S::U, "u_2"), v_1(g, S::V, "v_1"), v_2(g, S::V, "v_2"),
      w_1(g, S::W, "w_1"), w_2(g, S::W, "w_2"), t_1(g, S::M, "t_1"), t_2(g, S::M, "t_2"),
      ph_1(g, S::W, "ph_1"), ph_2(g, S::W, "ph_2"), mu_1(g, S::M, "mu_1"), mu_2(g, S::M, "mu_2"),
      p(g, S::M, "p"), al(g, S::M, "al"), alt(g, S::M, "alt"), php(g, S::M, "php"), ww(g, S::W, "ww"),
      mut(g, S::M, "mut"), muu(g, S::U, "muu"), muv(g, S::V, "muv"),
      muts(g, S::M, "muts"), muus(g, S::U, "muus"), muvs(g, S::V, "muvs"),
      ru(g, S::U, "ru"), rv(g, S::V, "rv"), rw(g, S::W, "rw"),
      ru_m(g, S::U, "ru_m"), rv_m(g, S::V, "rv_m"), ww_m(g, S::W, "ww_m") {}

Tendencies::Tendencies(const Grid& g)
    : ru_tend(g, S::U, "ru_tend"), rv_tend(g, S::V, "rv_tend"), rw_tend(g, S::W, "rw_tend"),
      t_tend(g, S::M, "t_tend"), ph_tend(g, S::W, "ph_tend"), mu_tend(g, S::M, "mu_tend"),
      ru_tendf(g, S::U, "ru_tendf"), rv_tendf(g, S::V, "rv_tendf"), rw_tendf(g, S::W, "rw_tendf"),
      t_tendf(g, S::M, "t_tendf"), ph_tendf(g, S::W, "ph_tendf") {}

int State::add_scalar(const Grid& g, const std::string& name, bool moist) {
  for (auto& n : q_names)
    if (n == name) throw std::runtime_error("duplicate scalar '" + name + "'");
  if (moist && n_moist != n_scalars()) throw std::runtime_error("moist species must precede passive tracers");
  q_names.push_back(name);
  q_1.emplace_back(g, S::M, name + "_1");
  q_2.emplace_back(g, S::M, name);
  const int idx = n_scalars() - 1;
  if (moist) {
    if (n_moist == 0) { qtot = Field3(g, S::M, "qtot"); h_diabatic = Field3(g, S::M, "h_diabatic"); }
    ++n_moist;
    if (name == "qv") P_QV = idx;
    if (name == "qc") P_QC = idx;
    if (name == "qr") P_QR = idx;
  }
  return idx;
}

void Tendencies::add_scalar(const Grid& g, const std::string& name) { q_tendf.emplace_back(g, S::M, name + "_tendf"); }

BaseState::BaseState(const Grid& g, bool flat_)
    : flat(flat_), mub(g, S::M, "mub", flat_), pb(g, S::M, "pb", flat_), alb(g, S::M, "alb", flat_),
      t_init(g, S::M, "t_init", flat_), qv_base(g, S::M, "qv_base", flat_), phb(g, S::W, "phb", flat_), p_top(g.p_top),
      znu(g.nz), znw(g.nz + 1), dnw(g.nz), dn(g.nz + 1), rdnw(g.nz), rdn(g.nz + 1), fnm(g.nz + 1),
      fnp(g.nz + 1), c1h(g.nz), c2h(g.nz), c3h(g.nz), c4h(g.nz), c1f(g.nz + 1), c2f(g.nz + 1),
      c3f(g.nz + 1), c4f(g.nz + 1), rdx(1.0 / g.dx), rdy(g.two_d() ? 0.0 : 1.0 / g.dy) {}

// WRF module_initialize_ideal.F (level setup), hybrid_opt = 0. Fortran k -> C++ k-1.
void BaseState::set_levels(const std::vector<double>& e) {
  const int nz = znu.size();
  if (int(e.size()) != nz + 1) throw std::runtime_error("set_levels: need nz+1 eta levels");
  if (e.front() != 1.0 || e.back() != 0.0) throw std::runtime_error("eta levels must go 1 -> 0");
  for (int k = 0; k < nz; ++k)
    if (!(e[k + 1] < e[k])) throw std::runtime_error("eta levels must decrease strictly");
  if (nz < 3) throw std::runtime_error("need nz >= 3");
  for (int k = 0; k <= nz; ++k) znw[k] = e[k];
  for (int k = 0; k < nz; ++k) {
    dnw[k] = znw[k + 1] - znw[k];
    rdnw[k] = 1.0 / dnw[k];
    znu[k] = 0.5 * (znw[k + 1] + znw[k]);
  }
  for (int k = 1; k < nz; ++k) {
    dn[k] = 0.5 * (dnw[k] + dnw[k - 1]);
    rdn[k] = 1.0 / dn[k];
    fnp[k] = 0.5 * dnw[k] / dn[k];
    fnm[k] = 0.5 * dnw[k - 1] / dn[k];
  }
  // Fortran dn(2), dn(3), dnw(1) -> C++ dn[1], dn[2], dnw[0]
  double cof1 = (2. * dn[1] + dn[2]) / (dn[1] + dn[2]) * dnw[0] / dn[1];
  double cof2 = dn[1] / (dn[1] + dn[2]) * dnw[0] / dn[2];
  cf1 = fnp[1] + cof1;
  cf2 = fnm[1] - cof1 - cof2;
  cf3 = cof2;
  cfn = (0.5 * dnw[nz - 1] + dn[nz - 1]) / dn[nz - 1];
  cfn1 = -0.5 * dnw[nz - 1] / dn[nz - 1];
  // hybrid_opt = 0 (pure terrain-following eta)
  for (int k = 0; k < nz; ++k) { c1h[k] = 1; c2h[k] = 0; c3h[k] = znu[k]; c4h[k] = 0; }
  for (int k = 0; k <= nz; ++k) { c1f[k] = 1; c2f[k] = 0; c3f[k] = znw[k]; c4f[k] = 0; }
}

void BaseState::upload_levels() {
  for (Array1* a : {&znu, &znw, &dnw, &dn, &rdnw, &rdn, &fnm, &fnp, &c1h, &c2h, &c3h, &c4h, &c1f, &c2f, &c3f, &c4f})
    a->upload();
}

}  // namespace cppwrf::dyn
