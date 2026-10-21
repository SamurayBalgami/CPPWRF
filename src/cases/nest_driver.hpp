#pragma once

// ============================================================================
// CPPWRF Multi-Domain Static Nest Driver (one-way & two-way nesting)
// Ports WRF 4.8.0 frame/module_integrate.F, share/mediation_integrate.F,
// share/mediation_force_domain.F, and share/mediation_feedback_domain.F.
// ============================================================================

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "cases/domain_runner.hpp"
#include "cases/model.hpp"
#include "io/config.hpp"

namespace cppwrf::cases {

// Per-domain nesting configuration ([domains] + [d0N] overrides)
struct NestDomainSpec {
  int dom_id = 1;                  // 1-based domain ID (1 = d01, 2 = d02, ...)
  int parent_id = 0;               // 1-based parent domain ID (0 for root d01)
  std::string tag = "d01";         // "d01", "d02", ...
  int i_parent_start = 1;          // 1-based coarse cell index in parent
  int j_parent_start = 1;          // 1-based coarse cell index in parent
  int parent_grid_ratio = 1;       // horizontal refinement ratio (odd/even)
  int parent_time_step_ratio = 1;  // temporal sub-stepping ratio
  int feedback = 1;                // 0 = one-way, 1 = two-way
  int smooth_option = 2;           // 0 = none, 1 = 1-2-1, 2 = smoother-desmoother
  int spec_bdy_width = 5;
  int blend_width = 5;
  io::Config cfg;                  // domain-specific Config passed to Model::build & DomainRunner
};

// Returns true if the configuration requests multi-domain nesting ([domains] max_dom > 1)
bool has_nests(const io::Config& root_cfg);

// Parse [domains] and per-domain [d0N] sections into a vector of NestDomainSpec (index 0 = d01)
std::vector<NestDomainSpec> parse_nest_specs(const io::Config& root_cfg);

// WRF med_nest_initial(parent, child):
//   1. SINT-interpolates parent ht, mub, phb onto child
//   2. Loads raw wrfinput_d0N fields (ht_fine, save_mub, phb, p, w=0)
//   3. Blends terrain (ht, mub, phb) if save_topo_from_real == 1
//   4. Adjusts child potential temperature (t_2) and water vapor (qv) via adjust_tempqv
//   5. If feedback != 0, applies med_nest_feedback(parent, child) at t = 0
//   6. Runs start_domain(child) with press_adj = true
//   7. Runs start_domain(parent) with press_adj = false
void med_nest_initial(const io::Config& parent_cfg, Model& parent,
                      const NestDomainSpec& child_spec, Model& child);

// WRF med_nest_force(parent, child):
//   Couples parent and child in-place, fills child LateralBC slabs via SINT (bdy_interp1)
//   with cdt = parent.g.dt, uncouples child and parent in-place, and resets child.bc->dtbc = 0.
void med_nest_force(Model& parent, Model& child, const NestDomainSpec& child_spec);

// WRF med_nest_feedback(parent, child):
//   Cell-averages (copy_fcn, spec_zone = 1) + smooths (smooth_option) child interior onto parent
//   for u_2, v_2, w_2, ph_2, t_2, mu_2, all scalars q_2, and base fields phb, mub, alb, pb, ht,
//   then refreshes parent halos (update_after_feedback_em).
void med_nest_feedback(Model& parent, const Model& child, const NestDomainSpec& child_spec);

// Build and run a multi-domain nested simulation using DomainRunner per domain.
// Returns 0 on success, non-zero on error (matching src/main.cpp exit codes).
int run_nested(const io::Config& root_cfg, const std::string& out_dir);

}  // namespace cppwrf::cases
