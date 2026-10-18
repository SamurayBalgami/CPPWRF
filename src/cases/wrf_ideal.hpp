#pragma once

#include <string>

#include "cases/model.hpp"
#include "io/config.hpp"

namespace cppwrf::cases {

// WRF v4.8.0 idealized initialization, statement by statement: dyn_em/module_initialize_ideal.F
// init_domain_rk (with get_sounding / read_sounding and module_init_utilities interp_0) for
// ideal_case hill2d_x, quarter_ss, squall2d_x, squall2d_y and grav2d_x, then what start_em does
// with the resulting wrfinput on the first trip of an ideal run (t_1 = t_2, ph_1 = ph_2,
// mu_1 = mu_2, pb and alb from MUB and t_init, al from ph, p from theta_m); WRF's input_wrf
// forces hypsometric_opt = 1 for ideal runs.
// run.case = em_hill2d_x | em_quarter_ss | em_squall2d_x | em_squall2d_y | em_grav2d_x.
// Namelist equivalents: domain.nx = e_we - 1, domain.ny = e_sn - 1 (1 for the x-z cases, WRF
// e_sn = 3), domain.nz = e_vert - 1, domain.dx, domain.dy, domain.ztop, dyn.hybrid_opt,
// dyn.etac, case.sounding_file (WRF input_sounding).
bool is_wrf_ideal_case(const std::string& name);
void init_wrf_ideal(const io::Config& cfg, Model& m, const std::string& name);

}  // namespace cppwrf::cases
