#pragma once

// Master header for WRF v4.6.1 Diffusion and Damping kernels in CPPWRF.
//

#include "dyn/diffusion/diffusion_params.hpp"
#include "dyn/diffusion/diff_metrics.hpp"
#include "dyn/diffusion/deform_div.hpp"
#include "dyn/diffusion/km_kh.hpp"
#include "dyn/diffusion/diffusion_2.hpp"
#include "dyn/diffusion/big_step_diffusion.hpp"
#include "dyn/diffusion/sixth_order.hpp"
#include "dyn/diffusion/w_damp.hpp"

// =============================================================================
// Input Halo Requirements per Launcher
// =============================================================================
//
// All halos refer to lateral (x, y) cells outside the logical tile domain:
//
// 1. compute_diff_metrics:
//    - ph, phb                     : Halo 0 (boundary/periodic edge conditions handled internally)
//
// 2. cal_deform_and_div:
//    - u                           : Halo 1 in Y (at vorticity points j-1; staggered 0..nx in X)
//    - v                           : Halo 1 in X (at vorticity points i-1; staggered 0..ny in Y)
//    - w                           : Halo 1 in X and Y (w(i-1) in defor13, w(j-1) in defor23)
//    - zx                          : Halo 1 in Y (zx(i, j-1) in defor12)
//    - zy                          : Halo 1 in X (zy(i-1, j) in defor12)
//    - rdzw                        : Halo 1 in X and Y (rdzw(i-1, j-1) in defor12)
//    - rdz                         : Halo 1 in X and Y (rdz(i-1) in defor13, rdz(j-1) in defor23)
//    - msft*, msfu*, msfv*         : Halo 0 (pointwise on respective staggered grids)
//    - u_base, v_base, 1D arrays   : Halo 0
//
// 3. calculate_km_kh:
//    - defor11, defor22, rdzw      : Halo 0 (pointwise at mass points 0..nx-1, 0..ny-1)
//    - defor12                     : Halo 1 on East & North (defor12(i+1, j+1); needs valid 0..nx, 0..ny)
//    - zx                          : Halo 1 on East (zx(i+1, j); needs valid 0..nx, 0..ny-1)
//    - zy                          : Halo 1 on North (zy(i, j+1); needs valid 0..nx-1, 0..ny)
//    - msftx, msfty                : Halo 0
//
// 4. horizontal_diffusion_2 (diff_opt = 2):
//    - thp, moist (scalars)        : Halo 1 in X and Y (stencils evaluate i-1..i+1, j-1..j+1)
//    - rho, xkhh, xkmh, xkmv       : Halo 1 in X and Y (averaged onto faces and edges)
//    - zx                          : Halo 1 on East & South (zx(i+1, j), zx(i, j-1))
//    - zy                          : Halo 1 on North & West (zy(i, j+1), zy(i-1, j))
//    - rdz, rdzw                   : Halo 1 on West & South (rdz(i-1, j), rdzw(i-1, j), rdzw(i, j-1))
//    - msft*, msfu*, msfv*         : Halo 0
//
// 5. vertical_diffusion_2 (diff_opt = 2):
//    - All fields                  : Halo 0 (purely 1D vertical columns per (i, j) point)
//
// 6. horizontal_diffusion (diff_opt = 1):
//    - field (u, v, w, or scalar m): Halo 1 in X and Y (2nd-order Laplacian stencils i+/-1, j+/-1)
//    - xkmhd                       : Halo 1 in X and Y (face-averaged i+/-1, j+/-1)
//    - mut                         : Halo 1 in X and Y (hybrid coordinate factor c1*mut + c2)
//    - msft*, msfu*, msfv*         : Halo 0
//
// 7. horizontal_diffusion_3dmp (diff_opt = 1):
//    - field, base_3d              : Halo 1 in X and Y
//    - xkmhd, mut                  : Halo 1 in X and Y
//    - msft*, msfu*, msfv*         : Halo 0
//
// 8. vertical_diffusion, vertical_diffusion_mp, vertical_diffusion_u, vertical_diffusion_v (diff_opt = 1):
//    - field (u, v, w, or scalar m): Halo 0 (purely 1D vertical columns)
//    - base, base_3d               : Halo 0
//    - mut, muu, muv               : Halo 0 (pointwise at column (i, j))
//    - alt                         : Halo 0 for w, m;
//                                    Halo 1 on West for u (alt(i-1, j) across u-staggered faces);
//                                    Halo 1 on South for v (alt(i, j-1) across v-staggered faces)
//
// 9. sixth_order_diffusion (diff_6th_opt = 1 or 2):
//    - field (u, v, w, or scalar m): Halo 3 in X and Y (7-point stencil: i-3..i+3, j-3..j+3)
//    - mut                         : Halo 1 in X and Y (face-averaged c1*mut + c2)
//    - msft*, msfu*, msfv*         : Halo 0
//
// 10. w_damp:
//    - All fields (u, v, ww, w, mut, msf*): Halo 0 (pointwise diagnostics and CFL damping)
// =============================================================================
