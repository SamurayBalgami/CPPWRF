#pragma once

#include "core/backend.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"
#include <algorithm>
#include <cmath>

namespace cppwrf::phys {

CPPWRF_HD inline void cldprmc_sw(
    int nlayers,
    int inflag,
    int iceflag,
    int liqflag,
    const double* cldfmc,      // (nlayers, ngptsw)
    const double* ciwpmc,      // (nlayers, ngptsw)
    const double* clwpmc,      // (nlayers, ngptsw)
    const double* cswpmc,      // (nlayers, ngptsw)
    const double* reicmc,      // (nlayers)
    const double* relqmc,      // (nlayers)
    const double* resnmc,      // (nlayers)
    const double* fsfcmc,      // (nlayers, ngptsw)
    double* taormc,            // (nlayers, ngptsw)
    double* taucmc,            // (nlayers, ngptsw)
    double* ssacmc,            // (nlayers, ngptsw)
    double* asmcmc,            // (nlayers, ngptsw)
    const RrtmgSwTables& tables
) {
  constexpr int ngptsw = 112;
  constexpr double eps = 1.0e-6;
  constexpr double cldmin = 1.0e-20;

  for (int lay = 0; lay < nlayers; ++lay) {
    for (int ig = 0; ig < ngptsw; ++ig) {
      int idx = lay * ngptsw + ig;
      taormc[idx] = taucmc[idx];
    }
  }

  for (int lay = 0; lay < nlayers; ++lay) {
    double radice = reicmc[lay];
    double radliq = relqmc[lay];
    double radsno = resnmc[lay];

    for (int ig = 0; ig < ngptsw; ++ig) {
      int idx = lay * ngptsw + ig;
      double cwp = ciwpmc[idx] + clwpmc[idx] + cswpmc[idx];

      if (cldfmc[idx] >= cldmin && (cwp >= cldmin || taucmc[idx] >= cldmin)) {
        if (inflag == 0) {
          double taucldorig_a = taucmc[idx];
          double ffp = fsfcmc[idx];
          double ffp1 = 1.0 - ffp;
          double ffpssa = 1.0 - ffp * ssacmc[idx];
          double ssacloud_a = ffp1 * ssacmc[idx] / ffpssa;
          double taucloud_a = ffpssa * taucldorig_a;

          taormc[idx] = taucldorig_a;
          ssacmc[idx] = ssacloud_a;
          taucmc[idx] = taucloud_a;
          asmcmc[idx] = (asmcmc[idx] - ffp) / ffp1;
        } else if (inflag >= 2) {
          double extcoice = 0.0, ssacoice = 0.0, gice = 0.0, forwice = 0.0;
          double extcosno = 0.0, ssacosno = 0.0, gsno = 0.0, forwsno = 0.0;

          if (ciwpmc[idx] + cswpmc[idx] > 0.0) {
            int ib = tables.ngb[ig]; // 16 to 29
            int ib_idx = ib - 16;   // 0 to 13

            if (iceflag == 1) {
              int icx = 4;
              if (tables.wavenum2[ib_idx] > 1.43e04) icx = 0;
              else if (tables.wavenum2[ib_idx] > 7.7e03) icx = 1;
              else if (tables.wavenum2[ib_idx] > 5.3e03) icx = 2;
              else if (tables.wavenum2[ib_idx] > 4.0e03) icx = 3;
              else if (tables.wavenum2[ib_idx] >= 2.5e03) icx = 4;

              extcoice = tables.abari[icx] + tables.bbari[icx] / radice;
              ssacoice = 1.0 - tables.cbari[icx] - tables.dbari[icx] * radice;
              gice = tables.ebari[icx] + tables.fbari[icx] * radice;
              if (gice >= 1.0) gice = 1.0 - eps;
              forwice = gice * gice;
            } else if (iceflag == 2) {
              double factor = (radice - 2.0) / 3.0;
              int index = static_cast<int>(factor);
              if (index == 43) index = 42;
              double fint = factor - static_cast<double>(index);
              int i0 = index - 1;
              extcoice = tables.extice2[ib_idx * 43 + i0] + fint * (tables.extice2[ib_idx * 43 + i0 + 1] - tables.extice2[ib_idx * 43 + i0]);
              ssacoice = tables.ssaice2[ib_idx * 43 + i0] + fint * (tables.ssaice2[ib_idx * 43 + i0 + 1] - tables.ssaice2[ib_idx * 43 + i0]);
              gice = tables.asyice2[ib_idx * 43 + i0] + fint * (tables.asyice2[ib_idx * 43 + i0 + 1] - tables.asyice2[ib_idx * 43 + i0]);
              forwice = gice * gice;
            } else if (iceflag >= 3) {
              double factor = (radice - 2.0) / 3.0;
              int index = static_cast<int>(factor);
              if (index >= 46) index = 45;
              double fint = factor - static_cast<double>(index);
              int i0 = index - 1;
              extcoice = tables.extice3[ib_idx * 46 + i0] + fint * (tables.extice3[ib_idx * 46 + i0 + 1] - tables.extice3[ib_idx * 46 + i0]);
              ssacoice = tables.ssaice3[ib_idx * 46 + i0] + fint * (tables.ssaice3[ib_idx * 46 + i0 + 1] - tables.ssaice3[ib_idx * 46 + i0]);
              gice = tables.asyice3[ib_idx * 46 + i0] + fint * (tables.asyice3[ib_idx * 46 + i0 + 1] - tables.asyice3[ib_idx * 46 + i0]);
              double fdelta = tables.fdlice3[ib_idx * 46 + i0] + fint * (tables.fdlice3[ib_idx * 46 + i0 + 1] - tables.fdlice3[ib_idx * 46 + i0]);
              forwice = fdelta + 0.5 / ssacoice;
              if (forwice > gice) forwice = gice;
            }

            if (cswpmc[idx] > 0.0 && iceflag == 5) {
              double factor = (radsno - 2.0) / 3.0;
              int index = static_cast<int>(factor);
              if (index >= 46) index = 45;
              double fint = factor - static_cast<double>(index);
              int i0 = index - 1;
              extcosno = tables.extice3[ib_idx * 46 + i0] + fint * (tables.extice3[ib_idx * 46 + i0 + 1] - tables.extice3[ib_idx * 46 + i0]);
              ssacosno = tables.ssaice3[ib_idx * 46 + i0] + fint * (tables.ssaice3[ib_idx * 46 + i0 + 1] - tables.ssaice3[ib_idx * 46 + i0]);
              gsno = tables.asyice3[ib_idx * 46 + i0] + fint * (tables.asyice3[ib_idx * 46 + i0 + 1] - tables.asyice3[ib_idx * 46 + i0]);
              double fdelta = tables.fdlice3[ib_idx * 46 + i0] + fint * (tables.fdlice3[ib_idx * 46 + i0 + 1] - tables.fdlice3[ib_idx * 46 + i0]);
              forwsno = fdelta + 0.5 / ssacosno;
              if (forwsno > gsno) forwsno = gsno;
            }
          }

          // Liquid water clouds
          double extcoliq = 0.0, ssacoliq = 0.0, gliq = 0.0, forwliq = 0.0;
          if (clwpmc[idx] > 0.0 && liqflag == 1) {
            int ib = tables.ngb[ig];
            int ib_idx = ib - 16;
            int index = static_cast<int>(radliq - 1.5);
            if (index <= 0) index = 1;
            if (index >= 58) index = 57;
            double fint = radliq - 1.5 - static_cast<double>(index);
            int i0 = index - 1;
            extcoliq = tables.extliq1[ib_idx * 58 + i0] + fint * (tables.extliq1[ib_idx * 58 + i0 + 1] - tables.extliq1[ib_idx * 58 + i0]);
            ssacoliq = tables.ssaliq1[ib_idx * 58 + i0] + fint * (tables.ssaliq1[ib_idx * 58 + i0 + 1] - tables.ssaliq1[ib_idx * 58 + i0]);
            if (fint < 0.0 && ssacoliq > 1.0) ssacoliq = tables.ssaliq1[ib_idx * 58 + i0];
            gliq = tables.asyliq1[ib_idx * 58 + i0] + fint * (tables.asyliq1[ib_idx * 58 + i0 + 1] - tables.asyliq1[ib_idx * 58 + i0]);
            forwliq = gliq * gliq;
          }

          double tauliqorig = clwpmc[idx] * extcoliq;
          double tauiceorig = ciwpmc[idx] * extcoice;
          double ssaliq = (clwpmc[idx] > 0.0) ? ssacoliq * (1.0 - forwliq) / (1.0 - forwliq * ssacoliq) : 0.0;
          double tauliq = (1.0 - forwliq * ssacoliq) * tauliqorig;
          double ssaice = (ciwpmc[idx] > 0.0) ? ssacoice * (1.0 - forwice) / (1.0 - forwice * ssacoice) : 0.0;
          double tauice = (1.0 - forwice * ssacoice) * tauiceorig;
          double scatliq = ssaliq * tauliq;
          double scatice = ssaice * tauice;

          if (iceflag < 5) {
            taormc[idx] = tauliqorig + tauiceorig;
            taucmc[idx] = tauliq + tauice;
            if (taucmc[idx] == 0.0) taucmc[idx] = cldmin;
            if (scatice == 0.0) scatice = cldmin;
            ssacmc[idx] = (scatliq + scatice) / taucmc[idx];

            if (iceflag == 3 || iceflag == 4) {
              asmcmc[idx] = (1.0 / (scatliq + scatice)) *
                            (scatliq * (gliq - forwliq) / (1.0 - forwliq) +
                             scatice * (gice - forwice) / (1.0 - forwice));
            } else {
              asmcmc[idx] = (scatliq * (gliq - forwliq) / (1.0 - forwliq) +
                             scatice * (gice - forwice) / (1.0 - forwice)) / (scatliq + scatice);
            }
          } else {
            double tausnoorig = cswpmc[idx] * extcosno;
            taormc[idx] = tauliqorig + tauiceorig + tausnoorig;
            double ssasno = (cswpmc[idx] > 0.0) ? ssacosno * (1.0 - forwsno) / (1.0 - forwsno * ssacosno) : 0.0;
            double tausno = (1.0 - forwsno * ssacosno) * tausnoorig;
            double scatsno = ssasno * tausno;
            taucmc[idx] = tauliq + tauice + tausno;
            if (taucmc[idx] == 0.0) taucmc[idx] = cldmin;
            if (scatice == 0.0) scatice = cldmin;
            if (scatsno == 0.0) scatsno = cldmin;
            ssacmc[idx] = (scatliq + scatice + scatsno) / taucmc[idx];

            asmcmc[idx] = (1.0 / (scatliq + scatice + scatsno)) *
                          (scatliq * (gliq - forwliq) / (1.0 - forwliq) +
                           scatice * (gice - forwice) / (1.0 - forwice) +
                           scatsno * (gsno - forwsno) / (1.0 - forwsno));
          }
        }
      }
    }
  }
}

}  // namespace cppwrf::phys
