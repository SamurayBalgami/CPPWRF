#pragma once

#include "core/backend.hpp"
#include "phys/ra_rrtmg_common/strided.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_tables.hpp"

namespace cppwrf::phys {

// Computes in-cloud optical depth taucmc[ig][lay] for McICA subcolumns.
// Matches WRF module_ra_rrtmg_lw subroutine cldprmc.
CPPWRF_HD inline void rrtmg_lw_cldprmc(
    int nlayers,
    int inflag,
    int iceflag,
    int liqflag,
    SD2 cldfmc,   // [140][nlayers] McICA cloud flags (0 or 1)
    SD ciwpth,    // [nlayers] in-cloud paths: the McICA subcolumn paths are these where cldfmc = 1
    SD clwpth,    //           and 0 elsewhere (generate_stochastic_clouds)
    SD cswpth,
    SD reicmc,          // [nlayers]
    SD relqmc,          // [nlayers]
    SD resnmc,          // [nlayers]
    SD2 taucmc,         // [140][nlayers] (inout)
    const RrtmgLwTables& tables,
    int lay0, int lay1  // layers [lay0, lay1)
) {
  constexpr double cldmin = 1.0e-20;
  const int icb[16] = {0, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4};

  for (int lay = lay0; lay < (nlayers < lay1 ? nlayers : lay1); ++lay) {
    for (int ig = 0; ig < 140; ++ig) {
      const bool cl = cldfmc[ig][lay] == 1.0;
      const double ciwp = cl ? ciwpth[lay] : 0.0, clwp = cl ? clwpth[lay] : 0.0, cswp = cl ? cswpth[lay] : 0.0;
      double cwp = ciwp + clwp + cswp;
      if (cldfmc[ig][lay] >= cldmin && (cwp >= cldmin || taucmc[ig][lay] >= cldmin)) {
        if (inflag == 0) {
          // Cloud optical depth already defined in taucmc
          return;
        } else if (inflag >= 2) {
          double radice = reicmc[lay];
          double abscoice = 0.0;
          double abscosno = 0.0;

          if ((ciwp + cswp) > 0.0) {
            if (iceflag == 0) {
              abscoice = tables.absice0[0] + tables.absice0[1] / radice;
            } else if (iceflag == 1) {
              int ib = icb[tables.ngb[ig] - 1];
              abscoice = tables.absice1[ib * 2 + 0] + tables.absice1[ib * 2 + 1] / radice;
            } else if (iceflag == 2) {
              double factor = (radice - 2.0) / 3.0;
              int index = static_cast<int>(factor);
              if (index >= 43) index = 42;
              if (index < 1) index = 1;
              double fint = factor - static_cast<double>(index);
              int ib = tables.ngb[ig] - 1;
              int idx0 = index - 1;
              abscoice = tables.absice2[ib * 43 + idx0] +
                         fint * (tables.absice2[ib * 43 + idx0 + 1] - tables.absice2[ib * 43 + idx0]);
            } else if (iceflag >= 3) {
              // Fu (1996) generalized effective size
              double factor = (radice - 2.0) / 3.0;
              int index = static_cast<int>(factor);
              if (index >= 46) index = 45;
              if (index < 1) index = 1;
              double fint = factor - static_cast<double>(index);
              int ib = tables.ngb[ig] - 1;
              int idx0 = index - 1;
              abscoice = tables.absice3[ib * 46 + idx0] +
                         fint * (tables.absice3[ib * 46 + idx0 + 1] - tables.absice3[ib * 46 + idx0]);
            }
          }

          // Incorporate snow if iceflag == 5
          if (cswp > 0.0 && iceflag == 5) {
            double radsno = resnmc[lay];
            double factor = (radsno - 2.0) / 3.0;
            int index = static_cast<int>(factor);
            if (index >= 46) index = 45;
            if (index < 1) index = 1;
            double fint = factor - static_cast<double>(index);
            int ib = tables.ngb[ig] - 1;
            int idx0 = index - 1;
            abscosno = tables.absice3[ib * 46 + idx0] +
                       fint * (tables.absice3[ib * 46 + idx0 + 1] - tables.absice3[ib * 46 + idx0]);
          }

          // Liquid water droplets
          double abscoliq = 0.0;
          if (clwp > 0.0) {
            if (liqflag == 0) {
              abscoliq = tables.absliq0;
            } else if (liqflag == 1) {
              // Hu and Stamnes (1993)
              double radliq = relqmc[lay];
              int index = static_cast<int>(radliq - 1.5);
              if (index < 1) index = 1;
              if (index >= 58) index = 57;
              double fint = radliq - 1.5 - static_cast<double>(index);
              int ib = tables.ngb[ig] - 1;
              int idx0 = index - 1;
              abscoliq = tables.absliq1[ib * 58 + idx0] +
                         fint * (tables.absliq1[ib * 58 + idx0 + 1] - tables.absliq1[ib * 58 + idx0]);
            }
          }

          taucmc[ig][lay] = ciwp * abscoice +
                            clwp * abscoliq +
                            cswp * abscosno;
        }
      } else {
        taucmc[ig][lay] = 0.0;
      }
    }
  }
}

}  // namespace cppwrf::phys
