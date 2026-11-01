#pragma once
// 5-layer thermal diffusion land surface, sf_surface_physics = 1: WRF 4.8.0 phys/module_sf_slab.F SLAB1D for one
// surface point (SLAB loops over the rows; a row without land skips to label 90, which gives each point the same
// result as the per-point land test). Every statement follows the Fortran in order. The THCON loop is vectorised
// by gfortran (libmvec pow lanes over its..ite: vm::pow with the point's i of the nx-point row); the exp and pow of
// loop 70 and of the soil sub-steps are scalar libm (loops with GOTOs). TG0**4 is (TG0*TG0)*(TG0*TG0).
#include "core/backend.hpp"

namespace cppwrf::phys::slab {

inline constexpr int kSoilLayers = 5;  // check_a_mundo: num_soil_layers = 5 for SLABSCHEME
inline constexpr real difsl = real(5.e-7);   // soil diffusion constant (m^2/s)
inline constexpr real soilfac = real(1.25);  // factor to make the soil step more conservative

struct SlabParams {
  real deltsm, rovcp, xlv, svp1, svp2, svp3, svpt0, ep2, stbolt, p1000mb;  // EOMEG: num_soil_layers = 1 only
  int ifsnow;
  bool radiation;
  real zs[kSoilLayers], dzs[kSoilLayers];
};

// t1d, qv1d, p1d: level kts of t_phy, qv, p_phy; psfcpa: surface pressure (Pa); tslb: the point's 5 layers;
// i, nx: the point's index in the row and the row length (libmvec lanes)
CPPWRF_HD inline void slab1d(const SlabParams& p, real t1d, real qv1d, real p1d, real psfcpa, real xland, real tmn,
                             real gsw, real glw, real thc, real snowc, real emiss, real mavail, real flhc, real flqc,
                             real& hfx, real& qfx, real& lh, real& tsk, real& qsfc, real& chklowq, real& capg,
                             real* tslb, int i, int nx) {
  constexpr int nsl = kSoilLayers;
  const real psfc = psfcpa / real(1000.);  // cb
  const real pl = p1d / real(1000.);
  const real scr3 = t1d;
  const real thcon = vm::pow(p.p1000mb * real(0.001) / pl, p.rovcp, i, nx);
  const real thx = scr3 * thcon;
  const real qx = qv1d;
  capg = real(5.9114e7) * thc;  // num_soil_layers > 1 (else 3.298E6*THC): CAPG*DIFSL**0.5 is the thermal inertia
  if ((xland - real(1.5)) < real(0)) {  // XLD1 = 1: land
    const real tg0 = tsk;
    const int radswtch = p.radiation ? 1 : 0;
    const real tg02 = tg0 * tg0;
    const real upflux = real(radswtch) * p.stbolt * (tg02 * tg02);
    const real xinet = emiss * (glw - upflux);
    const real rnet = gsw + xinet;
    // moisture flux calculated here (overwrites the surface layer's value for land)
    real esg = p.svp1 * std::exp(p.svp2 * (tg0 - p.svpt0) / (tg0 - p.svp3));
    real qsg = p.ep2 * esg / (psfc - esg);
    real thg = tsk * rpow(real(100.) / psfc, p.rovcp);
    hfx = flhc * (thg - thx);
    qfx = flqc * (qsg - qx);
    lh = qfx * p.xlv;
    real qs = hfx + qfx * p.xlv;
    real dthgdt = real(0.);  // num_soil_layers > 1
    const int nsoil = 1 + int(soilfac * real(4) * difsl / p.dzs[0] * p.deltsm / p.dzs[0]);
    const real rnsoil = real(1.) / real(nsoil);
    real flux[nsl];
    for (int itsoil = 1; itsoil <= nsoil; ++itsoil) {  // soil sub-timestep
      for (int l = 1; l <= nsl - 1; ++l) {
        if (l == 1 && itsoil > 1) {
          const real ps1 = rpow(psfcpa / p.p1000mb, p.rovcp);
          const real ps = psfc;
          thg = tslb[0] / ps1;
          esg = p.svp1 * std::exp(p.svp2 * (tslb[0] - p.svpt0) / (tslb[0] - p.svp3));
          qsg = p.ep2 * esg / (ps - esg);
          // update the fluxes for the new ground temperature and sum them over the soil time steps
          const real hfxt = flhc * (thg - thx);
          const real qfxt = flqc * (qsg - qx);
          qs = hfxt + qfxt * p.xlv;
          hfx = hfx + hfxt;
          qfx = qfx + qfxt;
        }
        flux[0] = rnet - qs;
        flux[l] = -(difsl * capg * (tslb[l] - tslb[l - 1]) / (p.zs[l] - p.zs[l - 1]));
        const real dtsdt = -(flux[l] - flux[l - 1]) / (p.dzs[l - 1] * capg);
        tslb[l - 1] = tslb[l - 1] + dtsdt * p.deltsm * rnsoil;
        if (p.ifsnow == 1 && l == 1) {
          if (snowc > real(0.) && tslb[0] > real(273.16)) tslb[0] = real(273.16);
        }
        if (l == 1) dthgdt = dthgdt + rnsoil * dtsdt;
        if (itsoil == nsoil && l == 1) {  // average HFX and QFX over the soil time steps for the PBL
          hfx = hfx * rnsoil;
          qfx = qfx * rnsoil;
          lh = qfx * p.xlv;
        }
      }
    }
    tsk = tg0 + p.deltsm * dthgdt;
    if (p.ifsnow != 0) {  // snow cover limits the ground temperature to 0 C
      if (snowc > real(0.) && tsk > real(273.16)) tsk = real(273.16);
    }
  }
  // label 90, every point: QSFC and CHKLOWQ for the Eta PBL (QX where FLQC = 0)
  if (flqc != real(0.)) {
    qsfc = qx + qfx / flqc;
  } else {
    qsfc = qx;
  }
  chklowq = mavail;
}

}  // namespace cppwrf::phys::slab
