#pragma once
// Idealized case framework (host side, setup only).
// Initialization follows WRF dyn_em/module_initialize_ideal.F for dry cases:
//   base state from the environmental sounding (horizontally uniform, flat terrain),
//   perturbation: mu' = 0, p' = 0, alt from the EOS with the perturbed theta at pb,
//   al = alt - alb, ph' from the hydrostatic relation (TN 2.4/2.5, discrete as in WRF).
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "core/grid.hpp"
#include "dyn/state.hpp"
#include "io/config.hpp"

namespace cppwrf::cases {

// Hydrostatic sounding theta(z) [, qv(z)]: Exner pi(z) = pi_s - g/cp * int_0^z dz'/theta_rho,
// tabulated, with theta_rho = theta (1 + rvovrd qv) / (1 + qv) (moist hydrostatic balance,
// dp/dz = -g rho_d (1 + qv), as WRF get_sounding). Dry pressure p_d(z) = p(z) - W(z), with the
// vapour load W(z) = int_z^top g rho_d qv dz' (WRF: p_dry integrated down from the moist p at the
// top). Without qv, theta_rho == theta and W == 0 exactly.
class Sounding {
 public:
  Sounding(std::function<double(double)> theta, double p_surf, double zmax, double dz = 0.5,
           std::function<double(double)> qv = {});
  double theta(double z) const { return th_(z); }
  double qv(double z) const { return qv_ ? qv_(z) : 0.0; }
  bool moist() const { return bool(qv_); }
  double p(double z) const;        // full pressure at height z
  double z_of_p(double p) const;   // inverse (monotone)
  double pd(double z) const;       // dry pressure
  double vapour_load(double z) const;  // W(z) = p(z) - pd(z)
  double z_of_pd(double pd) const;
 private:
  double theta_rho(double z) const;
  std::function<double(double)> th_, qv_;
  double dz_;
  std::vector<double> pi_;  // Exner at z = n*dz
  std::vector<double> w_;   // vapour load at z = n*dz (moist only)
};

struct CaseDef {
  std::string name;
  std::function<double(double)> theta_env;             // K, function of z
  std::function<double(double)> u_env;                 // m/s, function of z
  std::function<double(double, double, double, double)> theta_pert;  // (x, y, z, exner) -> K
  std::function<double(double)> qv_env;                // kg/kg, function of z (empty: dry)
  double p_surf = 1.0e5;
  std::map<std::string, double> params;                // written to run.json ("case_params")
};

// Builds the case definition from the config ([case] section) and the grid.
CaseDef make_case(const std::string& name, const io::Config& cfg, const Grid& g);

// Eta levels: "uniform_z" (base-state heights equally spaced), "uniform_eta", or explicit list.
std::vector<double> make_eta_levels(const io::Config& cfg, const Grid& g, const Sounding& snd);

// Fills base state and the initial state (both *_1 and *_2), halos included. Base state from the
// dry sounding `snd`; with moisture (s.P_QV >= 0) the full state comes from `moist` (WRF
// module_initialize_ideal: mu_1 from the dry surface pressure, p' by moist hydrostatic
// integration from the top, theta_m = theta (1 + rvovrd qv) in t_2).
void initialize(const CaseDef& c, const Sounding& snd, const std::vector<double>& eta, const Grid& g,
                dyn::BaseState& b, dyn::State& s, const Sounding* moist = nullptr);

}  // namespace cppwrf::cases
