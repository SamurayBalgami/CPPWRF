#pragma once
// Domain diagnostics (diag.csv row) and snapshot variables (host side, output steps only).
#include <string>
#include <vector>

#include "core/backend.hpp"
#include "dyn/state.hpp"
#include "io/output.hpp"

namespace cppwrf::dyn {

struct DiagRow {
  double max_abs_u, max_abs_v, max_abs_w, min_th_pert, max_th_pert, dry_mass, theta_mass, total_energy;
  double max_abs_q = 0;  // over all scalars (moist species, tracers): NaN checks only, not in diag.csv
  bool finite() const;
};

class Diagnostics {
 public:
  explicit Diagnostics(const Grid& g) : g_(g), red_(g.nz + 1) {}
  // Deterministic device reductions. dry_mass [kg], theta_mass [kg K], total_energy [J]
  // (kinetic + c_v T + geopotential, dry). NaN/Inf in a field -> +inf in its max.
  DiagRow compute(const State& s, const BaseState& b);

 private:
  Grid g_;
  Reducer red_;
};

// Snapshot variables (u, v, w, th_pert, theta, p, p_pert, ph, mu, mu_pert, z_m) computed on the
// host from downloaded fields; `store` keeps the buffers alive for the writer.
// Extra 2D fields (e.g. rainnc from the physics) are appended under their field names.
std::vector<io::VarView> snapshot_vars(const State& s, const BaseState& b, const Grid& g,
                                       std::vector<std::vector<real>>& store,
                                       const std::vector<const Field2*>& extra_2d = {});

}  // namespace cppwrf::dyn
