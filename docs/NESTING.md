# Nesting in CPPWRF: design

Scope: WRF 4.6.1 two-way (feedback = 1) and one-way (feedback = 0) static nests with odd
parent_grid_ratio, the same vertical levels in all domains (no vertical refinement), nests that read
their own wrfinput_d0N (input_from_file = .true.), real-data cases. Not implemented: ideal-case nests,
moving nests, vertical refinement.

## WRF call sequence (frame/module_integrate.F integrate, recursive per domain)
```
loop over parent steps:
  solve_interface(parent)                     one parent step (solve_em)
  for each child:
    med_nest_force(parent, child)             couple -> interp_domain_em_part1 (copy the coarse
                                              region onto the intermediate grid) ->
                                              force_domain_em_part2 -> nest_forcedown_interp.inc:
                                              bdy_interp for every field with a nest LBC
                                              -> uncouple
    integrate(child)                          parent_time_step_ratio child steps (recursive:
                                              the child forces/feeds back its own children)
    med_nest_feedback(parent, child)          if feedback: copy_fcn (+ smooth_option) of the
                                              child interior onto the parent
```
Nest LBC (share/interp_fcn.F bdy_interp1, SINT): for the boundary zone sz = min(max(spec_zone,
relax_zone+1), spec_bdy_width) = 5 of the child,
  value    = the child's current field (coupled),
  tendency = (sint(parent field at the new parent time) - child field) / dt_parent,
in the same _bxs/_bxe/_bys/_bye slabs as wrfbdy. dtbc restarts at 0 for each parent step and
advances by dt_child (so the child boundary moves linearly to the interpolated parent state).

## Differences of a nested domain from a specified (real-data) domain in solve_em
| piece | specified | nested |
|---|---|---|
| relax coefficients fcx, gcx | 0.1/dt, 1/dt/50 ramp x spongeweight(spec_exp) | same ramp, no spongeweight (module_bc_em.F 1312-1330) |
| W | zero_grad_bdy after advance_w | relax_bdy_dry also for w (mass_weight w, mut), spec_bdytend rw_tend, spec_bdyupdate w (solve_em 1603), spec_bdy_final w (solve_em 4566) |
| moist species | qv only (have_bcs_moist .false.), others flow_dep_bdy | all species relax_bdy_scalar + spec_bdy_scalar (solve_em 2257, 4622) |
| scalars | have_bcs_scalar | all (solve_em 2804) |
| adaptive dt | - | only when dtbc ~ 0 |
| microphysics spec-zone skip | yes | same |

## Initialisation (share/mediation_integrate.F med_nest_initial, input_from_file)
med_interp_domain(parent, nest) (all 'd' fields by interp_fcn; ht_int, mub_fine, phb_fine kept),
read wrfinput_d0N, blend_terrain(ht, mub, phb) when save_topo_from_real = 1, adjust_tempqv,
med_nest_feedback(parent, nest) right away (the parent changes at t = 0 under feedback),
start_domain(nest), start_domain(parent) again. Details: docs/NEST_NUMERICS.md.

## CPPWRF design
- A domain is a cases::Model (grid, base, state, tendencies, RK3, LateralBC,
  diffusion, physics). A `NestLink` holds parent*, child*, i/j_parent_start, ratio, time ratio,
  feedback, the intermediate coarse arrays and the interpolation kernels (src/nest).
- The time loop is `advance(domain)`: rk->step(domain); for each child: force(link);
  repeat ratio times advance(child); if feedback: feedback(link).
- LateralBC has (a) W as a boundary variable with the nested relax/spec/update/final path,
  (b) any number of scalar species with boundary data (all moist species in nests), (c) the
  nested fcx/gcx, (d) `set_interval_from_device` so force() writes slabs on the device
  without host round trips, and dtbc reset per parent step.
- Verification: WRF runs N1 (one-way) and N2 (two-way), real*8. CPPWRF on the CPU in FP64 is
  bit-identical to WRF on both domains (N2 for 1 h); see docs/NEST_NUMERICS.md section 8.
