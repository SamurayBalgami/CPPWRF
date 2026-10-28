#pragma once
// Fog (cloud) water deposition, grav_settling = 1 / 2: WRF 4.8.0 phys/module_sf_fogdes.F sf_fogdes (deposition
// velocity VDFG and the deposit DFGDP of one surface point, called by the surface driver after the surface layer)
// and phys/module_bl_fogdes.F bl_fogdes (gravitational settling of qc, added to RQCBLTEN by the PBL driver after
// the PBL scheme). Every statement follows the Fortran in order; x**2.0 is x*x, other real powers are pow.
#include "core/backend.hpp"
#include "phys/driver/wrf_constants.hpp"

namespace cppwrf::phys::fogdes {

inline constexpr real gno = real(1.0);  // "original value seems too aggressive: 4.64158883361278196"
inline constexpr real gpw = real(5.) / real(3.), qcgmin = real(1.e-8);

// lu_fogdes of sf_fogdes: 0 OTHERS, 1 CROP_GRASS, 2 MIXED_CROP_GRASS_WOOD, 3 SHRUB, 4 BROAD_FOREST,
// 5 MIXED_FOREST, 6 CONIFER_FOREST_ETC (usgs: MMINLU USGS, else MODIFIED_IGBP_MODIS_NOAH)
CPPWRF_HD inline int lu_fogdes(bool usgs, int v) {
  if (usgs) {
    if ((v >= 2 && v <= 15) || (v >= 17 && v <= 18) || (v >= 20 && v <= 22)) {
      if ((v >= 2 && v <= 5) || v == 7 || v == 17 || v == 20) return 1;
      if (v == 6 || v == 9) return 2;
      if (v == 8) return 3;
      if (v == 11 || v == 13) return 4;
      if (v == 15 || v == 22) return 5;
      return 6;
    }
    return 0;
  }
  if ((v >= 1 && v <= 10) || v == 12 || v == 14 || (v >= 18 && v <= 19)) {
    if (v == 10 || v == 12) return 1;
    if (v == 14) return 2;
    if (v == 6 || v == 7) return 3;
    if (v == 2 || v == 4) return 4;
    if (v == 5 || v == 19) return 5;
    return 6;
  }
  return 0;
}

// sf_fogdes for one point: qc, rho, dz8w of level kts; dfgdp = MAX(vdfg*lwc*dtbl, 0) (FGDP, the accumulation
// of DFGDP, feeds nothing and is not output)
CPPWRF_HD inline void sf_fogdes(int grav_settling, bool usgs, int ivgtyp, real lai, real wspd, real qc, real rho,
                                real dz8w, real dtbl, real& vdfg, real& dfgdp) {
  using wrfc::g;
  using wrfc::rhowater;
  constexpr real myu = real(1.8e-5);  // air viscosity (m^2/s)
  constexpr real gpw2 = real(0.66666666666667);
  const real lwc = rho * qc;
  if (grav_settling == 2) {
    const int lu = lu_fogdes(usgs, ivgtyp);
    if (lu == 0) {
      const real dp_fog = (real(17.03) * lwc * real(1.e3) + real(9.72)) * real(1.e-6);  // Katata et al. (2008)
      vdfg = (rhowater - rho) * (dp_fog * dp_fog) * g / (real(18.0) * myu);
    } else {
      real spcfct = real(1.0), vegh = real(20.0);
      if (lu == 1) {
        spcfct = real(0.2170);
        vegh = real(3.0);
      } else if (lu == 2) {
        spcfct = (real(1.0) + real(0.2170)) / real(2.0);
        vegh = (real(20.0) + real(3.0)) / real(2.0);
      } else if (lu == 3) {
        spcfct = real(1.0);
        vegh = real(4.0);
      } else if (lu == 4) {
        spcfct = real(0.8255);
        vegh = real(20.0);
      } else if (lu == 5) {
        spcfct = (real(1.0) + real(0.8255)) / real(2.0);
        vegh = real(20.0);
      }
      const real lad = lai / vegh;
      const real ftmp1 = real(0.0164) * rpow(lad, real(-0.5000));
      const real ftmp2 = real(0.0095) * rpow(lai, real(3.0)) - real(0.05) * (lai * lai) + real(0.0916) * lai +
                         real(0.0082);
      vdfg = spcfct * (ftmp1 < ftmp2 ? ftmp1 : ftmp2) * wspd;
    }
  } else if (grav_settling == 0) {
    vdfg = real(0.0);
  } else if (grav_settling == 1) {
    if ((qc / (real(1.) + qc)) > qcgmin)
      vdfg = gno * rpow(qc / (real(1.) + qc), gpw2);
    else
      vdfg = real(0.);
  }
  const real vmax = real(0.5) * dz8w / dtbl;
  vdfg = vmax < vdfg ? vmax : vdfg;
  dfgdp = vdfg * lwc * dtbl;
  dfgdp = dfgdp > real(0.0) ? dfgdp : real(0.0);
}

// bl_fogdes for one column (0-based k = 0 .. nz - 1 mass levels; qc, dz8w, dqc via accessors)
template <typename QC, typename DZ, typename DQ>
CPPWRF_HD inline void bl_fogdes(int grav_settling, real vdfg, int nz, const QC& qc, const DZ& dz8w, DQ& dqc) {
  const real grav_settling2 = real(grav_settling) < real(1.) ? real(grav_settling) : real(1.);  // INTEGER: 1
  const real gs2 = real(int(grav_settling2));
  real gfluxp, gfluxm;
  int k = 0;
  if (qc(k) > qcgmin)
    gfluxm = gs2 * qc(k) * vdfg;
  else
    gfluxm = real(0.);
  if (real(.5) * (qc(k + 1) + qc(k)) > qcgmin)
    gfluxp = gs2 * gno * rpow(real(.5) * (qc(k + 1) + qc(k)), gpw);
  else
    gfluxp = real(0.);
  dqc(k) = dqc(k) + (gfluxp - gfluxm) / dz8w(k);
  for (k = 1; k <= nz - 2; ++k) {
    if (real(.5) * (qc(k + 1) + qc(k)) > qcgmin)
      gfluxp = gs2 * gno * rpow(real(.5) * (qc(k + 1) + qc(k)), gpw);
    else
      gfluxp = real(0.);
    if (real(.5) * (qc(k - 1) + qc(k)) > qcgmin)
      gfluxm = gs2 * gno * rpow(real(.5) * (qc(k - 1) + qc(k)), gpw);
    else
      gfluxm = real(0.);
    dqc(k) = dqc(k) + (gfluxp - gfluxm) / dz8w(k);
  }
}

}  // namespace cppwrf::phys::fogdes
