// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "phys/ra_rrtmg_common/rrtmg_real.hpp"

#include "core/backend.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"
#include <algorithm>
#include <cmath>

namespace cppwrf::phys {

// cldprmc_sw for one (layer, g-point): the McICA subcolumn's cloud fraction (0 or 1) and paths,
// the layer's effective sizes; taucm/ssacm/asmcm hold the input optical properties (0, 1, 0 from
// the subcolumn generator) and return the cloud's, taorm the unscaled optical depth.
CPPWRF_HD inline void cldprmc_sw_elem(int inflag, int iceflag, int liqflag, int ig, rr cldf, rr ciwp,
                                      rr clwp, rr cswp, rr radice, rr radliq, rr radsno,
                                      rr fsfcm, rr& taorm, rr& taucm, rr& ssacm, rr& asmcm,
                                      const RrtmgSwTables& tables) {
  constexpr rr eps = rr(1.0e-6);
  constexpr rr cldmin = rr(1.0e-20);
  (void)eps;
  taorm = taucm;
  rr cwp = ciwp + clwp + cswp;

  if (cldf >= cldmin && (cwp >= cldmin || taucm >= cldmin)) {
    if (inflag == 0) {
      rr taucldorig_a = taucm;
      rr ffp = fsfcm;
      rr ffp1 = rr(1.0) - ffp;
      rr ffpssa = rr(1.0) - ffp * ssacm;
      rr ssacloud_a = ffp1 * ssacm / ffpssa;
      rr taucloud_a = ffpssa * taucldorig_a;

      taorm = taucldorig_a;
      ssacm = ssacloud_a;
      taucm = taucloud_a;
      asmcm = (asmcm - ffp) / ffp1;
    } else if (inflag >= 2) {
      rr extcoice = rr(0.0), ssacoice = rr(0.0), gice = rr(0.0), forwice = rr(0.0);
      rr extcosno = rr(0.0), ssacosno = rr(0.0), gsno = rr(0.0), forwsno = rr(0.0);

      if (ciwp + cswp > rr(0.0)) {
        int ib = tables.ngb[ig]; // 16 to 29
        int ib_idx = ib - 16;   // 0 to 13

        if (iceflag == 1) {
          int icx = 4;
          if (tables.wavenum2[ib_idx] > rr(1.43e04)) icx = 0;
          else if (tables.wavenum2[ib_idx] > rr(7.7e03)) icx = 1;
          else if (tables.wavenum2[ib_idx] > rr(5.3e03)) icx = 2;
          else if (tables.wavenum2[ib_idx] > rr(4.0e03)) icx = 3;
          else if (tables.wavenum2[ib_idx] >= rr(2.5e03)) icx = 4;

          extcoice = tables.abari[icx] + tables.bbari[icx] / radice;
          ssacoice = rr(1.0) - tables.cbari[icx] - tables.dbari[icx] * radice;
          gice = tables.ebari[icx] + tables.fbari[icx] * radice;
          if (gice >= rr(1.0)) gice = rr(1.0) - eps;
          forwice = gice * gice;
        } else if (iceflag == 2) {
          rr factor = (radice - rr(2.0)) / rr(3.0);
          int index = static_cast<int>(factor);
          if (index < 1) index = 1;  // WRF stops ('RADIUS OUT OF BOUNDS') below 5 um
          if (index == 43) index = 42;
          rr fint = factor - static_cast<rr>(index);
          int i0 = index - 1;
          extcoice = tables.extice2[ib_idx * 43 + i0] + fint * (tables.extice2[ib_idx * 43 + i0 + 1] - tables.extice2[ib_idx * 43 + i0]);
          ssacoice = tables.ssaice2[ib_idx * 43 + i0] + fint * (tables.ssaice2[ib_idx * 43 + i0 + 1] - tables.ssaice2[ib_idx * 43 + i0]);
          gice = tables.asyice2[ib_idx * 43 + i0] + fint * (tables.asyice2[ib_idx * 43 + i0 + 1] - tables.asyice2[ib_idx * 43 + i0]);
          forwice = gice * gice;
        } else if (iceflag >= 3) {
          rr factor = (radice - rr(2.0)) / rr(3.0);
          int index = static_cast<int>(factor);
          if (index < 1) index = 1;  // WRF stops ('RADIUS OUT OF BOUNDS') below 5 um
          if (index >= 46) index = 45;
          rr fint = factor - static_cast<rr>(index);
          int i0 = index - 1;
          extcoice = tables.extice3[ib_idx * 46 + i0] + fint * (tables.extice3[ib_idx * 46 + i0 + 1] - tables.extice3[ib_idx * 46 + i0]);
          ssacoice = tables.ssaice3[ib_idx * 46 + i0] + fint * (tables.ssaice3[ib_idx * 46 + i0 + 1] - tables.ssaice3[ib_idx * 46 + i0]);
          gice = tables.asyice3[ib_idx * 46 + i0] + fint * (tables.asyice3[ib_idx * 46 + i0 + 1] - tables.asyice3[ib_idx * 46 + i0]);
          rr fdelta = tables.fdlice3[ib_idx * 46 + i0] + fint * (tables.fdlice3[ib_idx * 46 + i0 + 1] - tables.fdlice3[ib_idx * 46 + i0]);
          forwice = fdelta + rr(0.5) / ssacoice;
          if (forwice > gice) forwice = gice;
        }

        if (cswp > rr(0.0) && iceflag == 5) {
          rr factor = (radsno - rr(2.0)) / rr(3.0);
          int index = static_cast<int>(factor);
          if (index < 1) index = 1;  // WRF stops ('RADIUS OUT OF BOUNDS') below 5 um
          if (index >= 46) index = 45;
          rr fint = factor - static_cast<rr>(index);
          int i0 = index - 1;
          extcosno = tables.extice3[ib_idx * 46 + i0] + fint * (tables.extice3[ib_idx * 46 + i0 + 1] - tables.extice3[ib_idx * 46 + i0]);
          ssacosno = tables.ssaice3[ib_idx * 46 + i0] + fint * (tables.ssaice3[ib_idx * 46 + i0 + 1] - tables.ssaice3[ib_idx * 46 + i0]);
          gsno = tables.asyice3[ib_idx * 46 + i0] + fint * (tables.asyice3[ib_idx * 46 + i0 + 1] - tables.asyice3[ib_idx * 46 + i0]);
          rr fdelta = tables.fdlice3[ib_idx * 46 + i0] + fint * (tables.fdlice3[ib_idx * 46 + i0 + 1] - tables.fdlice3[ib_idx * 46 + i0]);
          forwsno = fdelta + rr(0.5) / ssacosno;
          if (forwsno > gsno) forwsno = gsno;
        }
      }

      // Liquid water clouds
      rr extcoliq = rr(0.0), ssacoliq = rr(0.0), gliq = rr(0.0), forwliq = rr(0.0);
      if (clwp > rr(0.0) && liqflag == 1) {
        int ib = tables.ngb[ig];
        int ib_idx = ib - 16;
        int index = static_cast<int>(radliq - rr(1.5));
        if (index <= 0) index = 1;
        if (index >= 58) index = 57;
        rr fint = radliq - rr(1.5) - static_cast<rr>(index);
        int i0 = index - 1;
        extcoliq = tables.extliq1[ib_idx * 58 + i0] + fint * (tables.extliq1[ib_idx * 58 + i0 + 1] - tables.extliq1[ib_idx * 58 + i0]);
        ssacoliq = tables.ssaliq1[ib_idx * 58 + i0] + fint * (tables.ssaliq1[ib_idx * 58 + i0 + 1] - tables.ssaliq1[ib_idx * 58 + i0]);
        if (fint < rr(0.0) && ssacoliq > rr(1.0)) ssacoliq = tables.ssaliq1[ib_idx * 58 + i0];
        gliq = tables.asyliq1[ib_idx * 58 + i0] + fint * (tables.asyliq1[ib_idx * 58 + i0 + 1] - tables.asyliq1[ib_idx * 58 + i0]);
        forwliq = gliq * gliq;
      }

      rr tauliqorig = clwp * extcoliq;
      rr tauiceorig = ciwp * extcoice;
      rr ssaliq = (clwp > rr(0.0)) ? ssacoliq * (rr(1.0) - forwliq) / (rr(1.0) - forwliq * ssacoliq) : rr(0.0);
      rr tauliq = (rr(1.0) - forwliq * ssacoliq) * tauliqorig;
      rr ssaice = (ciwp > rr(0.0)) ? ssacoice * (rr(1.0) - forwice) / (rr(1.0) - forwice * ssacoice) : rr(0.0);
      rr tauice = (rr(1.0) - forwice * ssacoice) * tauiceorig;
      rr scatliq = ssaliq * tauliq;
      rr scatice = ssaice * tauice;

      if (iceflag < 5) {
        taorm = tauliqorig + tauiceorig;
        taucm = tauliq + tauice;
        if (taucm == rr(0.0)) taucm = cldmin;
        if (scatice == rr(0.0)) scatice = cldmin;
        ssacm = (scatliq + scatice) / taucm;

        if (iceflag == 3 || iceflag == 4) {
          // WRF grouping: scatliq*(gliq - forwliq)/(1 - forwliq) + scatice*((gice - forwice)/(1 - forwice))
          asmcm = (rr(1.0) / (scatliq + scatice)) *
                        (scatliq * (gliq - forwliq) / (rr(1.0) - forwliq) +
                         scatice * ((gice - forwice) / (rr(1.0) - forwice)));
        } else {
          asmcm = (scatliq * (gliq - forwliq) / (rr(1.0) - forwliq) +
                         scatice * (gice - forwice) / (rr(1.0) - forwice)) / (scatliq + scatice);
        }
      } else {
        rr tausnoorig = cswp * extcosno;
        taorm = tauliqorig + tauiceorig + tausnoorig;
        rr ssasno = (cswp > rr(0.0)) ? ssacosno * (rr(1.0) - forwsno) / (rr(1.0) - forwsno * ssacosno) : rr(0.0);
        rr tausno = (rr(1.0) - forwsno * ssacosno) * tausnoorig;
        rr scatsno = ssasno * tausno;
        taucm = tauliq + tauice + tausno;
        if (taucm == rr(0.0)) taucm = cldmin;
        if (scatice == rr(0.0)) scatice = cldmin;
        if (scatsno == rr(0.0)) scatsno = cldmin;
        ssacm = (scatliq + scatice + scatsno) / taucm;

        asmcm = (rr(1.0) / (scatliq + scatice + scatsno)) *
                      (scatliq * (gliq - forwliq) / (rr(1.0) - forwliq) +
                       scatice * ((gice - forwice) / (rr(1.0) - forwice)) +
                       scatsno * ((gsno - forwsno) / (rr(1.0) - forwsno)));
      }
    }
  }
}

}  // namespace cppwrf::phys
