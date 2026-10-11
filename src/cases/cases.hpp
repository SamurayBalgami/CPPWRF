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

// Dry hydrostatic sounding theta(z): Exner pi(z) = pi_s - g/cp * int_0^z dz'/theta, tabulated.
class Sounding {
 public:
  Sounding(std::function<double(double)> theta, double p_surf, double zmax, double dz = 0.5);
  double theta(double z) const { return th_(z); }
  double p(double z) const;        // dry pressure at height z
  double z_of_p(double p) const;   // inverse (monotone)
 private:
  std::function<double(double)> th_;
  double dz_;
  std::vector<double> pi_;  // Exner at z = n*dz
};

struct CaseDef {
  std::string name;
  std::function<double(double)> theta_env;             // K, function of z
  std::function<double(double)> u_env;                 // m/s, function of z
  std::function<double(double, double, double, double)> theta_pert;  // (x, y, z, exner) -> K
  double p_surf = 1.0e5;
  std::map<std::string, double> params;                // written to run.json ("case_params")
};

// Builds the case definition from the config ([case] section) and the grid.
CaseDef make_case(const std::string& name, const io::Config& cfg, const Grid& g);

// Eta levels: "uniform_z" (base-state heights equally spaced), "uniform_eta", or explicit list.
std::vector<double> make_eta_levels(const io::Config& cfg, const Grid& g, const Sounding& snd);

// Fills base state and the initial state (both *_1 and *_2), halos included.
void initialize(const CaseDef& c, const Sounding& snd, const std::vector<double>& eta, const Grid& g,
                dyn::BaseState& b, dyn::State& s);

}  // namespace cppwrf::cases
