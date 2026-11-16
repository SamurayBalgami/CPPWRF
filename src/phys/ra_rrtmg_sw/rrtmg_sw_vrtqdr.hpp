// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "phys/ra_rrtmg_common/rrtmg_real.hpp"

#include "core/backend.hpp"
#include "phys/ra_rrtmg_common/strided.hpp"

namespace cppwrf::phys {

template <int ML>
CPPWRF_FORCEINLINE void vrtqdr_sw_direct(
    int klev,
    const rr* pref,       // (klev+1)
    const rr* prefd,      // (klev+1)
    const rr* ptra,       // (klev+1)
    const rr* ptrad,      // (klev+1)
    const rr* pdbt,       // (klev+1)
    const rr* ptdbt,      // (klev+1)
    rr zincflx,
    SD2 pfd_out,
    SD2 pfu_out,
    int iw
) {
  rr ztdn[ML + 1], prdnd[ML + 1], prup[ML + 1], prupd[ML + 1];

  // Link lowest layer with surface
  rr zreflect = rr(1.0) / (rr(1.0) - prefd[klev] * prefd[klev - 1]);
  prup[klev - 1] = pref[klev - 1] + (ptrad[klev - 1] *
      ((ptra[klev - 1] - pdbt[klev - 1]) * prefd[klev] +
       pdbt[klev - 1] * pref[klev])) * zreflect;
  prupd[klev - 1] = prefd[klev - 1] + ptrad[klev - 1] * ptrad[klev - 1] *
      prefd[klev] * zreflect;

  // Pass from bottom to top
  for (int ikx = klev - 2; ikx >= 0; --ikx) {
    int ikp = ikx + 1;
    zreflect = rr(1.0) / (rr(1.0) - prupd[ikp] * prefd[ikx]);
    prup[ikx] = pref[ikx] + (ptrad[ikx] *
        ((ptra[ikx] - pdbt[ikx]) * prupd[ikp] +
         pdbt[ikx] * prup[ikp])) * zreflect;
    prupd[ikx] = prefd[ikx] + ptrad[ikx] * ptrad[ikx] *
        prupd[ikp] * zreflect;
  }

  // Upper boundary conditions
  ztdn[0] = rr(1.0);
  prdnd[0] = rr(0.0);
  ztdn[1] = ptra[0];
  prdnd[1] = prefd[0];

  // Pass from top to bottom
  for (int jk = 1; jk < klev; ++jk) {
    int ikp = jk + 1;
    zreflect = rr(1.0) / (rr(1.0) - prefd[jk] * prdnd[jk]);
    ztdn[ikp] = ptdbt[jk] * ptra[jk] +
        (ptrad[jk] * ((ztdn[jk] - ptdbt[jk]) +
         ptdbt[jk] * pref[jk] * prdnd[jk])) * zreflect;
    prdnd[ikp] = prefd[jk] + ptrad[jk] * ptrad[jk] *
        prdnd[jk] * zreflect;
  }

  // Up- and down-welling fluxes at levels
  prup[klev] = pref[klev];
  prupd[klev] = prefd[klev];
  for (int jk = 0; jk <= klev; ++jk) {
    int ikl = klev - jk;
    zreflect = rr(1.0) / (rr(1.0) - prdnd[jk] * prupd[jk]);
    rr pfu = (ptdbt[jk] * prup[jk] +
        (ztdn[jk] - ptdbt[jk]) * prupd[jk]) * zreflect;
    rr pfd = ptdbt[jk] + (ztdn[jk] - ptdbt[jk] +
        ptdbt[jk] * prup[jk] * prdnd[jk]) * zreflect;
    pfu_out[ikl][iw] = zincflx * pfu;
    pfd_out[ikl][iw] = zincflx * pfd;
  }
}

}  // namespace cppwrf::phys
