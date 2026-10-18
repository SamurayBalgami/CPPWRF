#pragma once

#include <algorithm>
#include <cmath>
#include "core/backend.hpp"
#include "core/precision.hpp"
#include "phys/lsm_noahmp/noahmp_var.hpp"

namespace cppwrf::phys::noahmp {

// ============================================================================
// VaporPressureSaturationMod.F90
// ============================================================================
CPPWRF_HD inline void VaporPressureSaturation(real T,
                                              real& VapPresSatWat,
                                              real& VapPresSatIce,
                                              real& VapPresSatWatD,
                                              real& VapPresSatIceD) {
  constexpr real A0 = real(6.107799961);
  constexpr real A1 = real(4.436518521e-01);
  constexpr real A2 = real(1.428945805e-02);
  constexpr real A3 = real(2.650648471e-04);
  constexpr real A4 = real(3.031240396e-06);
  constexpr real A5 = real(2.034080948e-08);
  constexpr real A6 = real(6.136820929e-11);
  constexpr real B0 = real(6.109177956);
  constexpr real B1 = real(5.034698970e-01);
  constexpr real B2 = real(1.886013408e-02);
  constexpr real B3 = real(4.176223716e-04);
  constexpr real B4 = real(5.824720280e-06);
  constexpr real B5 = real(4.838803174e-08);
  constexpr real B6 = real(1.838826904e-10);
  constexpr real C0 = real(4.438099984e-01);
  constexpr real C1 = real(2.857002636e-02);
  constexpr real C2 = real(7.938054040e-04);
  constexpr real C3 = real(1.215215065e-05);
  constexpr real C4 = real(1.036561403e-07);
  constexpr real C5 = real(3.532421810e-10);
  constexpr real C6 = real(-7.090244804e-13);
  constexpr real D0 = real(5.030305237e-01);
  constexpr real D1 = real(3.773255020e-02);
  constexpr real D2 = real(1.267995369e-03);
  constexpr real D3 = real(2.477563108e-05);
  constexpr real D4 = real(3.005693132e-07);
  constexpr real D5 = real(2.158542548e-09);
  constexpr real D6 = real(7.131097725e-12);

  VapPresSatWat  = real(100.0) * (A0 + T * (A1 + T * (A2 + T * (A3 + T * (A4 + T * (A5 + T * A6))))));
  VapPresSatIce  = real(100.0) * (B0 + T * (B1 + T * (B2 + T * (B3 + T * (B4 + T * (B5 + T * B6))))));
  VapPresSatWatD = real(100.0) * (C0 + T * (C1 + T * (C2 + T * (C3 + T * (C4 + T * (C5 + T * C6))))));
  VapPresSatIceD = real(100.0) * (D0 + T * (D1 + T * (D2 + T * (D3 + T * (D4 + T * (D5 + T * D6))))));
}

// ============================================================================
// HumiditySaturationMod.F90
// ============================================================================
CPPWRF_HD inline void HumiditySaturation(real TemperatureAir,
                                         real PressureAir,
                                         real& MixingRatioSat,
                                         real& MixingRatioSatTempD) {
  constexpr real Const1          = real(17.67);
  constexpr real TemperatureFrz  = real(273.15);
  constexpr real Const2          = real(29.65);
  constexpr real ConstLatHeatVap = real(2.501e6);
  constexpr real Const3          = Const1 * (TemperatureFrz - Const2);
  constexpr real VapPresSatFrz   = real(0.611);
  constexpr real GasConstWatVap  = real(461.0);
  constexpr real RatioGasConst   = real(0.622);

  real VapPresSatTemp = VapPresSatFrz * rexp(ConstLatHeatVap / GasConstWatVap *
                                             (real(1.0) / TemperatureFrz - real(1.0) / TemperatureAir));
  real PressureAirKpa = PressureAir * real(1.0e-3);
  MixingRatioSat = RatioGasConst * VapPresSatTemp / (PressureAirKpa - VapPresSatTemp);
  MixingRatioSat = MixingRatioSat * real(1.0e3);
  MixingRatioSatTempD = (MixingRatioSat / (real(1.0) + MixingRatioSat)) * Const3 /
                        ((TemperatureAir - Const2) * (TemperatureAir - Const2));
  MixingRatioSat = MixingRatioSat / real(1.0e3);
}

// ============================================================================
// SnowCoverGroundNiu07Mod.F90
// ============================================================================
CPPWRF_HD inline void SnowCoverGroundNiu07(NoahmpType& noahmp) {
  auto& water = noahmp.water;

  water.state.SnowCoverFrac = real(0.0);
  if (water.state.SnowDepth > real(0.0)) {
    real SnowDensBulk = water.state.SnowWaterEquiv / water.state.SnowDepth;
    real MeltFac = rpow(SnowDensBulk / real(100.0), water.param.SnowMeltFac);
    water.state.SnowCoverFrac = std::tanh(water.state.SnowDepth / (water.param.SnowCoverFac * MeltFac));
  }
}

// ============================================================================
// SnowCoverGroundAR25Mod.F90
// ============================================================================
CPPWRF_HD inline void SnowCoverGroundAR25(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& water = noahmp.water;

  water.state.SnowCoverFrac = real(0.0);
  if (water.state.SnowDepth > real(0.0)) {
    real SnowMeltFac =
        water.param.SnowCoverM1AR25 +
        std::tanh(water.param.SnowCoverM2AR25 * std::min((std::max(domain.GridSize, real(500.0)) / real(1000.0)), real(36.0)));
    real SnowCoverFac =
        water.param.SnowCoverFac1AR25 *
            std::sinh(water.param.SnowCoverFac2AR25 * std::min((std::max(domain.GridSize, real(500.0)) / real(1000.0)), real(36.0))) +
        water.param.SnowCoverFac2AR25;
    real SnowDensBulk = water.state.SnowWaterEquiv / water.state.SnowDepth;
    real MeltFac = rpow(SnowDensBulk / real(100.0), SnowMeltFac);
    water.state.SnowCoverFrac = std::tanh(water.state.SnowDepth / (SnowCoverFac * MeltFac));
  }
}

// ============================================================================
// GroundRoughnessPropertyMod.F90
// ============================================================================
CPPWRF_HD inline void GroundRoughnessProperty(NoahmpType& noahmp, bool FlagVegSfc) {
  auto& domain = noahmp.config.domain;
  auto& water = noahmp.water;
  auto& energy = noahmp.energy;

  if (domain.SurfaceType == 2) {  // Lake
    if (energy.state.TemperatureGrd <= ConstFreezePoint) {
      energy.state.RoughLenMomGrd =
          energy.param.RoughLenMomLake * (real(1.0) - water.state.SnowCoverFrac) +
          water.state.SnowCoverFrac * energy.param.RoughLenMomSnow;
    } else {
      energy.state.RoughLenMomGrd = energy.param.RoughLenMomLake;
    }
  } else {  // soil
    energy.state.RoughLenMomGrd =
        energy.param.RoughLenMomSoil * (real(1.0) - water.state.SnowCoverFrac) +
        water.state.SnowCoverFrac * energy.param.RoughLenMomSnow;
  }

  energy.state.ZeroPlaneDispGrd = water.state.SnowDepth;
  if (FlagVegSfc) {
    energy.state.RoughLenMomSfc = energy.param.RoughLenMomVeg;
    energy.state.ZeroPlaneDispSfc = real(0.65) * energy.param.HeightCanopyTop;
    if (water.state.SnowDepth > energy.state.ZeroPlaneDispSfc) {
      energy.state.ZeroPlaneDispSfc = water.state.SnowDepth;
    }
  } else {
    energy.state.RoughLenMomSfc = energy.state.RoughLenMomGrd;
    energy.state.ZeroPlaneDispSfc = energy.state.ZeroPlaneDispGrd;
  }

  if (domain.FlagUrban) {
    energy.state.RoughLenMomGrd = energy.param.RoughLenMomVeg;
    energy.state.ZeroPlaneDispGrd = real(0.65) * energy.param.HeightCanopyTop;
    energy.state.RoughLenMomSfc = energy.state.RoughLenMomGrd;
    energy.state.ZeroPlaneDispSfc = energy.state.ZeroPlaneDispGrd;
  }

  energy.state.RefHeightAboveGrd = energy.state.ZeroPlaneDispSfc + domain.RefHeightAboveSfc;
  if (energy.state.ZeroPlaneDispGrd >= energy.state.RefHeightAboveGrd) {
    energy.state.RefHeightAboveGrd = energy.state.ZeroPlaneDispGrd + domain.RefHeightAboveSfc;
  }
}

// ============================================================================
// SnowAgingBatsMod.F90
// ============================================================================
CPPWRF_HD inline void SnowAgingBats(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& water = noahmp.water;
  auto& energy = noahmp.energy;

  if (water.state.SnowWaterEquiv <= real(0.0)) {
    energy.state.SnowAgeNondim = real(0.0);
  } else {
    real SnowAgeTimeFac = domain.MainTimeStep / energy.param.SnowAgeFacBats;
    real SnowGrowVapExp =
        energy.param.SnowGrowVapFacBats * (real(1.0) / ConstFreezePoint - real(1.0) / energy.state.TemperatureGrd);
    real SnowAgeVapEff = rexp(SnowGrowVapExp);
    real SnowAgeFrzEff = rexp(std::min(real(0.0), energy.param.SnowGrowFrzFacBats * SnowGrowVapExp));
    real SnowAgeSootEff = energy.param.SnowSootFacBats;
    real SnowAgeFacTot = SnowAgeVapEff + SnowAgeFrzEff + SnowAgeSootEff;
    real SnowAgeChg = SnowAgeTimeFac * SnowAgeFacTot;
    real SnowFreshFac =
        std::max(real(0.0), water.state.SnowWaterEquiv - water.state.SnowWaterEquivPrev) /
        water.param.SnowMassFullCoverOld;
    real SnowAgeTmp = (energy.state.SnowAgeNondim + SnowAgeChg) * (real(1.0) - SnowFreshFac);
    energy.state.SnowAgeNondim = std::max(real(0.0), SnowAgeTmp);
  }

  energy.state.SnowAgeFac = energy.state.SnowAgeNondim / (energy.state.SnowAgeNondim + real(1.0));
}

// ============================================================================
// SnowAlbedoBatsMod.F90
// ============================================================================
CPPWRF_HD inline void SnowAlbedoBats(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& energy = noahmp.energy;

  for (int IndBand = 1; IndBand <= domain.NumSwRadBand; ++IndBand) {
    energy.state.AlbedoSnowDir(IndBand) = real(0.0);
    energy.state.AlbedoSnowDif(IndBand) = real(0.0);
  }

  real SolarAngleFac = energy.param.SolarZenithAdjBats;
  real SolarAngleFac1 = real(1.0) / SolarAngleFac;
  real SolarAngleFac2 = real(2.0) * SolarAngleFac;
  real ZenithAngFacTmp =
      (real(1.0) + SolarAngleFac1) / (real(1.0) + SolarAngleFac2 * domain.CosSolarZenithAngle) - SolarAngleFac1;
  real ZenithAngFac = std::max(ZenithAngFacTmp, real(0.0));
  energy.state.AlbedoSnowDif(1) =
      energy.param.FreshSnoAlbVisBats * (real(1.0) - energy.param.SnoAgeFacDifVisBats * energy.state.SnowAgeFac);
  energy.state.AlbedoSnowDif(2) =
      energy.param.FreshSnoAlbNirBats * (real(1.0) - energy.param.SnoAgeFacDifNirBats * energy.state.SnowAgeFac);
  energy.state.AlbedoSnowDir(1) =
      energy.state.AlbedoSnowDif(1) +
      energy.param.SzaFacDirVisBats * ZenithAngFac * (real(1.0) - energy.state.AlbedoSnowDif(1));
  energy.state.AlbedoSnowDir(2) =
      energy.state.AlbedoSnowDif(2) +
      energy.param.SzaFacDirNirBats * ZenithAngFac * (real(1.0) - energy.state.AlbedoSnowDif(2));
}

// ============================================================================
// SnowAlbedoClassMod.F90
// ============================================================================
CPPWRF_HD inline void SnowAlbedoClass(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& water = noahmp.water;
  auto& energy = noahmp.energy;

  for (int IndBand = 1; IndBand <= domain.NumSwRadBand; ++IndBand) {
    energy.state.AlbedoSnowDir(IndBand) = real(0.0);
    energy.state.AlbedoSnowDif(IndBand) = real(0.0);
  }

  real SnowAlbedoTmp =
      energy.param.SnowAlbRefClass +
      (energy.state.AlbedoSnowPrev - energy.param.SnowAlbRefClass) *
          rexp(real(-0.01) * domain.MainTimeStep / energy.param.SnowAgeFacClass);

  if (water.flux.SnowfallGround > real(0.0)) {
    SnowAlbedoTmp = SnowAlbedoTmp +
                    std::min(water.flux.SnowfallGround, water.param.SnowMassFullCoverOld / domain.MainTimeStep) *
                        (energy.param.SnowAlbFreshClass - SnowAlbedoTmp) /
                        (water.param.SnowMassFullCoverOld / domain.MainTimeStep);
  }

  energy.state.AlbedoSnowDif(1) = SnowAlbedoTmp;
  energy.state.AlbedoSnowDif(2) = SnowAlbedoTmp;
  energy.state.AlbedoSnowDir(1) = SnowAlbedoTmp;
  energy.state.AlbedoSnowDir(2) = SnowAlbedoTmp;

  energy.state.AlbedoSnowPrev = SnowAlbedoTmp;
}

// ============================================================================
// GroundAlbedoMod.F90
// ============================================================================
CPPWRF_HD inline void GroundAlbedo(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& water = noahmp.water;
  auto& energy = noahmp.energy;

  for (int IndSwBnd = 1; IndSwBnd <= domain.NumSwRadBand; ++IndSwBnd) {
    real AlbedoSoilAdjWet = std::max(real(0.11) - real(0.40) * water.state.SoilMoisture(1), real(0.0));

    if (domain.SurfaceType == 1) {  // soil
      energy.state.AlbedoSoilDir(IndSwBnd) =
          std::min(energy.param.AlbedoSoilSat(IndSwBnd) + AlbedoSoilAdjWet, energy.param.AlbedoSoilDry(IndSwBnd));
      energy.state.AlbedoSoilDif(IndSwBnd) = energy.state.AlbedoSoilDir(IndSwBnd);
    } else if (energy.state.TemperatureGrd > ConstFreezePoint) {  // unfrozen lake, wetland
      energy.state.AlbedoSoilDir(IndSwBnd) =
          real(0.06) / (rpow(std::max(real(0.01), domain.CosSolarZenithAngle), real(1.7)) + real(0.15));
      energy.state.AlbedoSoilDif(IndSwBnd) = real(0.06);
    } else {  // frozen lake, wetland
      energy.state.AlbedoSoilDir(IndSwBnd) = energy.param.AlbedoLakeFrz(IndSwBnd);
      energy.state.AlbedoSoilDif(IndSwBnd) = energy.state.AlbedoSoilDir(IndSwBnd);
    }

    energy.state.AlbedoGrdDir(IndSwBnd) =
        energy.state.AlbedoSoilDir(IndSwBnd) * (real(1.0) - water.state.SnowCoverFrac) +
        energy.state.AlbedoSnowDir(IndSwBnd) * water.state.SnowCoverFrac;
    energy.state.AlbedoGrdDif(IndSwBnd) =
        energy.state.AlbedoSoilDif(IndSwBnd) * (real(1.0) - water.state.SnowCoverFrac) +
        energy.state.AlbedoSnowDif(IndSwBnd) * water.state.SnowCoverFrac;
  }
}

// ============================================================================
// CanopyRadiationTwoStreamMod.F90
// ============================================================================
CPPWRF_HD inline void CanopyRadiationTwoStream(NoahmpType& noahmp, int IndSwBnd, int IndSwDif) {
  auto& nml = noahmp.config.nmlist;
  auto& domain = noahmp.config.domain;
  auto& water = noahmp.water;
  auto& energy = noahmp.energy;

  if (energy.state.VegAreaIndEff == real(0.0)) {
    energy.state.GapCanopyDir = real(1.0);
    energy.state.GapCanopyDif = real(1.0);
  } else {
    if (nml.OptCanopyRadiationTransfer == 1) {
      real VegDensity =
          -rlog(std::max(real(1.0) - energy.state.VegFrac, real(0.01))) /
          (ConstPI * (energy.param.TreeCrownRadius * energy.param.TreeCrownRadius));
      real CrownDepth = energy.param.HeightCanopyTop - energy.param.HeightCanopyBot;
      real CrownRadiusVert = real(0.5) * CrownDepth;
      real SolarAngleTmp =
          std::atan(CrownRadiusVert / energy.param.TreeCrownRadius *
                    std::tan(std::acos(std::max(real(0.01), domain.CosSolarZenithAngle))));
      energy.state.GapBtwCanopy =
          rexp(-VegDensity * ConstPI * (energy.param.TreeCrownRadius * energy.param.TreeCrownRadius) /
               std::cos(SolarAngleTmp));
      real FoliageDensity =
          energy.state.VegAreaIndEff /
          (real(1.33) * ConstPI * rpow(energy.param.TreeCrownRadius, real(3.0)) *
           (CrownRadiusVert / energy.param.TreeCrownRadius) * VegDensity);
      real VegAreaIndTmp = CrownDepth * FoliageDensity;
      energy.state.GapInCanopy =
          (real(1.0) - energy.state.GapBtwCanopy) * rexp(real(-0.5) * VegAreaIndTmp / domain.CosSolarZenithAngle);
      energy.state.GapCanopyDir =
          std::min(real(1.0) - energy.state.VegFrac, energy.state.GapBtwCanopy + energy.state.GapInCanopy);
      energy.state.GapCanopyDif = real(0.05);
    }
    if (nml.OptCanopyRadiationTransfer == 2) {
      energy.state.GapCanopyDir = real(0.0);
      energy.state.GapCanopyDif = real(0.0);
    }
    if (nml.OptCanopyRadiationTransfer == 3) {
      energy.state.GapCanopyDir = real(1.0) - energy.state.VegFrac;
      energy.state.GapCanopyDif = real(1.0) - energy.state.VegFrac;
    }
  }

  real CosSolarZenithAngleTmp = std::max(real(0.001), domain.CosSolarZenithAngle);
  real LeafOrientIndex = std::min(std::max(energy.param.CanopyOrientIndex, real(-0.4)), real(0.6));
  if (std::abs(LeafOrientIndex) <= real(0.01)) LeafOrientIndex = real(0.01);
  real Phi1 = real(0.5) - real(0.633) * LeafOrientIndex - real(0.330) * LeafOrientIndex * LeafOrientIndex;
  real Phi2 = real(0.877) * (real(1.0) - real(2.0) * Phi1);
  energy.state.VegAreaProjDir = Phi1 + Phi2 * CosSolarZenithAngleTmp;
  real OpticDepthDir = energy.state.VegAreaProjDir / CosSolarZenithAngleTmp;
  real OpticDepthDif = (real(1.0) - Phi1 / Phi2 * rlog((Phi1 + Phi2) / Phi1)) / Phi2;
  real ScatCoeffLeaf = energy.state.ReflectanceVeg(IndSwBnd) + energy.state.TransmittanceVeg(IndSwBnd);
  real Tmp0 = energy.state.VegAreaProjDir + Phi2 * CosSolarZenithAngleTmp;
  real Tmp1 = Phi1 * CosSolarZenithAngleTmp;
  real SingleScatAlb =
      real(0.5) * ScatCoeffLeaf * energy.state.VegAreaProjDir / Tmp0 *
      (real(1.0) - Tmp1 / Tmp0 * rlog((Tmp1 + Tmp0) / Tmp1));
  real UpscatCoeffLeafDir =
      (real(1.0) + OpticDepthDif * OpticDepthDir) / (ScatCoeffLeaf * OpticDepthDif * OpticDepthDir) * SingleScatAlb;
  real UpscatCoeffLeafDif =
      real(0.5) *
      (energy.state.ReflectanceVeg(IndSwBnd) + energy.state.TransmittanceVeg(IndSwBnd) +
       (energy.state.ReflectanceVeg(IndSwBnd) - energy.state.TransmittanceVeg(IndSwBnd)) *
           cppwrf::powi((real(1.0) + LeafOrientIndex) / real(2.0), 2)) /
      ScatCoeffLeaf;

  real Tmp2 = real(0.0);
  if (energy.state.TemperatureCanopy > ConstFreezePoint) {
    Tmp0 = ScatCoeffLeaf;
    Tmp1 = UpscatCoeffLeafDir;
    Tmp2 = UpscatCoeffLeafDif;
  } else {
    Tmp0 = (real(1.0) - water.state.CanopyWetFrac) * ScatCoeffLeaf +
           water.state.CanopyWetFrac * energy.param.ScatterCoeffSnow(IndSwBnd);
    Tmp1 = ((real(1.0) - water.state.CanopyWetFrac) * ScatCoeffLeaf * UpscatCoeffLeafDir +
            water.state.CanopyWetFrac * energy.param.ScatterCoeffSnow(IndSwBnd) * energy.param.UpscatterCoeffSnowDir) /
           Tmp0;
    Tmp2 = ((real(1.0) - water.state.CanopyWetFrac) * ScatCoeffLeaf * UpscatCoeffLeafDif +
            water.state.CanopyWetFrac * energy.param.ScatterCoeffSnow(IndSwBnd) * energy.param.UpscatterCoeffSnowDif) /
           Tmp0;
  }
  real ScatCoeffCan = Tmp0;
  real UpscatCoeffCanDir = Tmp1;
  real UpscatCoeffCanDif = Tmp2;

  real B = real(1.0) - ScatCoeffCan + ScatCoeffCan * UpscatCoeffCanDif;
  real C = ScatCoeffCan * UpscatCoeffCanDif;
  Tmp0 = OpticDepthDif * OpticDepthDir;
  real D = Tmp0 * ScatCoeffCan * UpscatCoeffCanDir;
  real F = Tmp0 * ScatCoeffCan * (real(1.0) - UpscatCoeffCanDir);
  Tmp1 = B * B - C * C;
  real H = std::sqrt(Tmp1) / OpticDepthDif;
  real Sigma = Tmp0 * Tmp0 - Tmp1;
  if (std::abs(Sigma) < real(1.0e-6)) Sigma = sign_f(real(1.0e-6), Sigma);
  real P1 = B + OpticDepthDif * H;
  real P2 = B - OpticDepthDif * H;
  real P3 = B + Tmp0;
  real P4 = B - Tmp0;
  real S1 = rexp(-H * energy.state.VegAreaIndEff);
  real S2 = rexp(-OpticDepthDir * energy.state.VegAreaIndEff);
  real U1 = real(0.0);
  real U2 = real(0.0);
  real U3 = real(0.0);
  if (IndSwDif == 0) {
    U1 = B - C / energy.state.AlbedoGrdDir(IndSwBnd);
    U2 = B - C * energy.state.AlbedoGrdDir(IndSwBnd);
    U3 = F + C * energy.state.AlbedoGrdDir(IndSwBnd);
  } else {
    U1 = B - C / energy.state.AlbedoGrdDif(IndSwBnd);
    U2 = B - C * energy.state.AlbedoGrdDif(IndSwBnd);
    U3 = F + C * energy.state.AlbedoGrdDif(IndSwBnd);
  }
  Tmp2 = U1 - OpticDepthDif * H;
  real Tmp3 = U1 + OpticDepthDif * H;
  real D1 = P1 * Tmp2 / S1 - P2 * Tmp3 * S1;
  real Tmp4 = U2 + OpticDepthDif * H;
  real Tmp5 = U2 - OpticDepthDif * H;
  real D2 = Tmp4 / S1 - Tmp5 * S1;
  real H1 = -D * P4 - C * F;
  real Tmp6 = D - H1 * P3 / Sigma;
  real Tmp7 = (D - C - H1 / Sigma * (U1 + Tmp0)) * S2;
  real H2 = (Tmp6 * Tmp2 / S1 - P2 * Tmp7) / D1;
  real H3 = -(Tmp6 * Tmp3 * S1 - P1 * Tmp7) / D1;
  real H4 = -F * P3 - C * D;
  real Tmp8 = H4 / Sigma;
  real Tmp9 = (U3 - Tmp8 * (U2 - Tmp0)) * S2;
  real H5 = -(Tmp8 * Tmp4 / S1 + Tmp9) / D2;
  real H6 = (Tmp8 * Tmp5 * S1 + Tmp9) / D2;
  real H7 = (C * Tmp2) / (D1 * S1);
  real H8 = (-C * Tmp3 * S1) / D1;
  real H9 = Tmp4 / (D2 * S1);
  real H10 = (-Tmp5 * S1) / D2;

  real RadSwTransDir = real(0.0);
  real RadSwTransDif = real(0.0);
  if (IndSwDif == 0) {
    RadSwTransDir = S2 * (real(1.0) - energy.state.GapCanopyDir) + energy.state.GapCanopyDir;
    RadSwTransDif = (H4 * S2 / Sigma + H5 * S1 + H6 / S1) * (real(1.0) - energy.state.GapCanopyDir);
  } else {
    RadSwTransDir = real(0.0);
    RadSwTransDif = (H9 * S1 + H10 / S1) * (real(1.0) - energy.state.GapCanopyDif) + energy.state.GapCanopyDif;
  }
  if (IndSwDif == 0) {
    energy.flux.RadSwDirTranGrdDir(IndSwBnd) = RadSwTransDir;
    energy.flux.RadSwDifTranGrdDir(IndSwBnd) = RadSwTransDif;
  } else {
    energy.flux.RadSwDirTranGrdDif(IndSwBnd) = RadSwTransDir;
    energy.flux.RadSwDifTranGrdDif(IndSwBnd) = RadSwTransDif;
  }

  real RadSwReflTot = real(0.0);
  real RadSwReflCan = real(0.0);
  real RadSwReflGrd = real(0.0);
  if (IndSwDif == 0) {
    RadSwReflTot = (H1 / Sigma + H2 + H3) * (real(1.0) - energy.state.GapCanopyDir) +
                   energy.state.AlbedoGrdDir(IndSwBnd) * energy.state.GapCanopyDir;
    RadSwReflCan = (H1 / Sigma + H2 + H3) * (real(1.0) - energy.state.GapCanopyDir);
    RadSwReflGrd = energy.state.AlbedoGrdDir(IndSwBnd) * energy.state.GapCanopyDir;
  } else {
    RadSwReflTot = (H7 + H8) * (real(1.0) - energy.state.GapCanopyDif) +
                   energy.state.AlbedoGrdDif(IndSwBnd) * energy.state.GapCanopyDif;
    RadSwReflCan = (H7 + H8) * (real(1.0) - energy.state.GapCanopyDif) +
                   energy.state.AlbedoGrdDif(IndSwBnd) * energy.state.GapCanopyDif;
    RadSwReflGrd = real(0.0);
  }
  if (IndSwDif == 0) {
    energy.state.AlbedoSfcDir(IndSwBnd) = RadSwReflTot;
    energy.flux.RadSwReflVegDir(IndSwBnd) = RadSwReflCan;
    energy.flux.RadSwReflGrdDir(IndSwBnd) = RadSwReflGrd;
  } else {
    energy.state.AlbedoSfcDif(IndSwBnd) = RadSwReflTot;
    energy.flux.RadSwReflVegDif(IndSwBnd) = RadSwReflCan;
    energy.flux.RadSwReflGrdDif(IndSwBnd) = RadSwReflGrd;
  }

  if (IndSwDif == 0) {
    energy.flux.RadSwAbsVegDir(IndSwBnd) =
        real(1.0) - energy.state.AlbedoSfcDir(IndSwBnd) -
        (real(1.0) - energy.state.AlbedoGrdDir(IndSwBnd)) * energy.flux.RadSwDirTranGrdDir(IndSwBnd) -
        (real(1.0) - energy.state.AlbedoGrdDif(IndSwBnd)) * energy.flux.RadSwDifTranGrdDir(IndSwBnd);
  } else {
    energy.flux.RadSwAbsVegDif(IndSwBnd) =
        real(1.0) - energy.state.AlbedoSfcDif(IndSwBnd) -
        (real(1.0) - energy.state.AlbedoGrdDir(IndSwBnd)) * energy.flux.RadSwDirTranGrdDif(IndSwBnd) -
        (real(1.0) - energy.state.AlbedoGrdDif(IndSwBnd)) * energy.flux.RadSwDifTranGrdDif(IndSwBnd);
  }
}

// ============================================================================
// SurfaceAlbedoMod.F90
// ============================================================================
CPPWRF_HD inline void SurfaceAlbedo(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& nml = noahmp.config.nmlist;
  auto& energy = noahmp.energy;

  real MinThr = real(1.0e-06);
  energy.state.GapBtwCanopy = real(0.0);
  energy.state.GapInCanopy = real(0.0);
  energy.state.VegAreaProjDir = real(0.0);
  energy.state.ReflectanceVeg.fill(real(0.0));
  energy.state.TransmittanceVeg.fill(real(0.0));
  energy.state.CanopySunlitFrac = real(0.0);
  for (int IndBand = 1; IndBand <= domain.NumSwRadBand; ++IndBand) {
    energy.state.AlbedoSfcDir(IndBand) = real(0.0);
    energy.state.AlbedoSfcDif(IndBand) = real(0.0);
    energy.state.AlbedoGrdDir(IndBand) = real(0.0);
    energy.state.AlbedoGrdDif(IndBand) = real(0.0);
    energy.state.AlbedoSnowDir(IndBand) = real(0.0);
    energy.state.AlbedoSnowDif(IndBand) = real(0.0);
    energy.flux.RadSwAbsVegDir(IndBand) = real(0.0);
    energy.flux.RadSwAbsVegDif(IndBand) = real(0.0);
    energy.flux.RadSwDirTranGrdDir(IndBand) = real(0.0);
    energy.flux.RadSwDirTranGrdDif(IndBand) = real(0.0);
    energy.flux.RadSwDifTranGrdDir(IndBand) = real(0.0);
    energy.flux.RadSwDifTranGrdDif(IndBand) = real(0.0);
    energy.flux.RadSwReflVegDir(IndBand) = real(0.0);
    energy.flux.RadSwReflVegDif(IndBand) = real(0.0);
    energy.flux.RadSwReflGrdDir(IndBand) = real(0.0);
    energy.flux.RadSwReflGrdDif(IndBand) = real(0.0);
  }
  energy.state.VegAreaIndEff = energy.state.LeafAreaIndEff + energy.state.StemAreaIndEff;

  if (nml.OptSnowAlbedo == 1) SnowAgingBats(noahmp);

  if (domain.CosSolarZenithAngle > real(0.0)) {
    real LeafWgt = energy.state.LeafAreaIndEff / std::max(energy.state.VegAreaIndEff, MinThr);
    real StemWgt = energy.state.StemAreaIndEff / std::max(energy.state.VegAreaIndEff, MinThr);
    for (int IndBand = 1; IndBand <= domain.NumSwRadBand; ++IndBand) {
      energy.state.ReflectanceVeg(IndBand) =
          std::max(energy.param.ReflectanceLeaf(IndBand) * LeafWgt +
                       energy.param.ReflectanceStem(IndBand) * StemWgt,
                   MinThr);
      energy.state.TransmittanceVeg(IndBand) =
          std::max(energy.param.TransmittanceLeaf(IndBand) * LeafWgt +
                       energy.param.TransmittanceStem(IndBand) * StemWgt,
                   MinThr);
    }

    if (nml.OptSnowAlbedo == 1) SnowAlbedoBats(noahmp);
    if (nml.OptSnowAlbedo == 2) SnowAlbedoClass(noahmp);

    GroundAlbedo(noahmp);

    for (int IndBand = 1; IndBand <= domain.NumSwRadBand; ++IndBand) {
      int IndDif = 0;
      CanopyRadiationTwoStream(noahmp, IndBand, IndDif);
      IndDif = 1;
      CanopyRadiationTwoStream(noahmp, IndBand, IndDif);
    }

    real LightExtDir =
        energy.state.VegAreaProjDir / domain.CosSolarZenithAngle *
        std::sqrt(real(1.0) - energy.state.ReflectanceVeg(1) - energy.state.TransmittanceVeg(1));
    energy.state.CanopySunlitFrac =
        (real(1.0) - rexp(-LightExtDir * energy.state.VegAreaIndEff)) /
        std::max(LightExtDir * energy.state.VegAreaIndEff, MinThr);
    LightExtDir = energy.state.CanopySunlitFrac;
    if (LightExtDir < real(0.01)) {
      LeafWgt = real(0.0);
    } else {
      LeafWgt = LightExtDir;
    }
    energy.state.CanopySunlitFrac = LeafWgt;
  }

  energy.state.CanopyShadeFrac = real(1.0) - energy.state.CanopySunlitFrac;
  energy.state.LeafAreaIndSunlit = energy.state.LeafAreaIndEff * energy.state.CanopySunlitFrac;
  energy.state.LeafAreaIndShade = energy.state.LeafAreaIndEff * energy.state.CanopyShadeFrac;
}

// ============================================================================
// SurfaceRadiationMod.F90
// ============================================================================
CPPWRF_HD inline void SurfaceRadiation(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& energy = noahmp.energy;

  FArray1D<real, 1, kNumRadBands> RadSwAbsCanDir;
  FArray1D<real, 1, kNumRadBands> RadSwAbsCanDif;
  RadSwAbsCanDir.fill(real(0.0));
  RadSwAbsCanDif.fill(real(0.0));

  real MinThr = real(1.0e-6);
  energy.flux.RadSwAbsGrd = real(0.0);
  energy.flux.RadSwAbsVeg = real(0.0);
  energy.flux.RadSwAbsSfc = real(0.0);
  energy.flux.RadSwReflSfc = real(0.0);
  energy.flux.RadSwReflVeg = real(0.0);
  energy.flux.RadSwReflGrd = real(0.0);
  energy.flux.RadPhotoActAbsSunlit = real(0.0);
  energy.flux.RadPhotoActAbsShade = real(0.0);

  for (int IndBand = 1; IndBand <= domain.NumSwRadBand; ++IndBand) {
    RadSwAbsCanDir(IndBand) = energy.flux.RadSwDownDir(IndBand) * energy.flux.RadSwAbsVegDir(IndBand);
    RadSwAbsCanDif(IndBand) = energy.flux.RadSwDownDif(IndBand) * energy.flux.RadSwAbsVegDif(IndBand);
    energy.flux.RadSwAbsVeg = energy.flux.RadSwAbsVeg + RadSwAbsCanDir(IndBand) + RadSwAbsCanDif(IndBand);
    energy.flux.RadSwAbsSfc = energy.flux.RadSwAbsSfc + RadSwAbsCanDir(IndBand) + RadSwAbsCanDif(IndBand);

    real RadSwTranGrdDir = energy.flux.RadSwDownDir(IndBand) * energy.flux.RadSwDirTranGrdDir(IndBand);
    real RadSwTranGrdDif = energy.flux.RadSwDownDir(IndBand) * energy.flux.RadSwDifTranGrdDir(IndBand) +
                           energy.flux.RadSwDownDif(IndBand) * energy.flux.RadSwDifTranGrdDif(IndBand);

    real RadSwAbsGrdTmp = RadSwTranGrdDir * (real(1.0) - energy.state.AlbedoGrdDir(IndBand)) +
                          RadSwTranGrdDif * (real(1.0) - energy.state.AlbedoGrdDif(IndBand));
    energy.flux.RadSwAbsGrd = energy.flux.RadSwAbsGrd + RadSwAbsGrdTmp;
    energy.flux.RadSwAbsSfc = energy.flux.RadSwAbsSfc + RadSwAbsGrdTmp;
  }

  real LeafAreaIndFrac = energy.state.LeafAreaIndEff / std::max(energy.state.VegAreaIndEff, MinThr);
  if (energy.state.CanopySunlitFrac > real(0.0)) {
    energy.flux.RadPhotoActAbsSunlit =
        (RadSwAbsCanDir(1) + energy.state.CanopySunlitFrac * RadSwAbsCanDif(1)) *
        LeafAreaIndFrac / std::max(energy.state.LeafAreaIndSunlit, MinThr);
    energy.flux.RadPhotoActAbsShade =
        (energy.state.CanopyShadeFrac * RadSwAbsCanDif(1)) *
        LeafAreaIndFrac / std::max(energy.state.LeafAreaIndShade, MinThr);
  } else {
    energy.flux.RadPhotoActAbsSunlit = real(0.0);
    energy.flux.RadPhotoActAbsShade =
        (RadSwAbsCanDir(1) + RadSwAbsCanDif(1)) *
        LeafAreaIndFrac / std::max(energy.state.LeafAreaIndShade, MinThr);
  }

  real RadSwReflSfcVis =
      energy.state.AlbedoSfcDir(1) * energy.flux.RadSwDownDir(1) +
      energy.state.AlbedoSfcDif(1) * energy.flux.RadSwDownDif(1);
  real RadSwReflSfcNir =
      energy.state.AlbedoSfcDir(2) * energy.flux.RadSwDownDir(2) +
      energy.state.AlbedoSfcDif(2) * energy.flux.RadSwDownDif(2);
  energy.flux.RadSwReflSfc = RadSwReflSfcVis + RadSwReflSfcNir;

  energy.flux.RadSwReflVeg =
      energy.flux.RadSwReflVegDir(1) * energy.flux.RadSwDownDir(1) +
      energy.flux.RadSwReflVegDif(1) * energy.flux.RadSwDownDif(1) +
      energy.flux.RadSwReflVegDir(2) * energy.flux.RadSwDownDir(2) +
      energy.flux.RadSwReflVegDif(2) * energy.flux.RadSwDownDif(2);
  energy.flux.RadSwReflGrd =
      energy.flux.RadSwReflGrdDir(1) * energy.flux.RadSwDownDir(1) +
      energy.flux.RadSwReflGrdDif(1) * energy.flux.RadSwDownDif(1) +
      energy.flux.RadSwReflGrdDir(2) * energy.flux.RadSwDownDir(2) +
      energy.flux.RadSwReflGrdDif(2) * energy.flux.RadSwDownDif(2);
}

// ============================================================================
// SurfaceEmissivityMod.F90
// ============================================================================
CPPWRF_HD inline void SurfaceEmissivity(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& water = noahmp.water;
  auto& energy = noahmp.energy;

  energy.state.EmissivityVeg =
      real(1.0) - rexp(-(energy.state.LeafAreaIndEff + energy.state.StemAreaIndEff) / real(1.0));

  if (domain.IndicatorIceSfc == 1) {
    energy.state.EmissivityGrd =
        energy.param.EmissivityIceSfc * (real(1.0) - water.state.SnowCoverFrac) +
        energy.param.EmissivitySnow * water.state.SnowCoverFrac;
  } else {
    energy.state.EmissivityGrd =
        energy.param.EmissivitySoilLake(domain.SurfaceType) * (real(1.0) - water.state.SnowCoverFrac) +
        energy.param.EmissivitySnow * water.state.SnowCoverFrac;
  }

  energy.state.EmissivitySfc =
      energy.state.VegFrac *
          (energy.state.EmissivityGrd * (real(1.0) - energy.state.EmissivityVeg) + energy.state.EmissivityVeg +
           energy.state.EmissivityVeg * (real(1.0) - energy.state.EmissivityVeg) *
               (real(1.0) - energy.state.EmissivityGrd)) +
      (real(1.0) - energy.state.VegFrac) * energy.state.EmissivityGrd;
}

// ============================================================================
// PsychrometricVariableMod.F90
// ============================================================================
CPPWRF_HD inline void PsychrometricVariable(NoahmpType& noahmp) {
  auto& forcing = noahmp.forcing;
  auto& energy = noahmp.energy;

  if (energy.state.TemperatureCanopy > ConstFreezePoint) {
    energy.state.LatHeatVapCanopy = ConstLatHeatEvap;
    energy.state.FlagFrozenCanopy = false;
  } else {
    energy.state.LatHeatVapCanopy = ConstLatHeatSublim;
    energy.state.FlagFrozenCanopy = true;
  }
  energy.state.PsychConstCanopy =
      ConstHeatCapacAir * forcing.PressureAirRefHeight / (real(0.622) * energy.state.LatHeatVapCanopy);

  if (energy.state.TemperatureGrd > ConstFreezePoint) {
    energy.state.LatHeatVapGrd = ConstLatHeatEvap;
    energy.state.FlagFrozenGround = false;
  } else {
    energy.state.LatHeatVapGrd = ConstLatHeatSublim;
    energy.state.FlagFrozenGround = true;
  }
  energy.state.PsychConstGrd =
      ConstHeatCapacAir * forcing.PressureAirRefHeight / (real(0.622) * energy.state.LatHeatVapGrd);
}

// ============================================================================
// SnowThermalPropertyMod.F90
// ============================================================================
CPPWRF_HD inline void SnowThermalProperty(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& nml = noahmp.config.nmlist;
  auto& water = noahmp.water;
  auto& energy = noahmp.energy;

  FArray1D<real, -kMaxSnowLayers + 1, 0> SnowDensBulk;
  SnowDensBulk.fill(real(0.0));

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
    water.state.SnowIceVol(LoopInd) =
        std::min(real(1.0),
                 water.state.SnowIce(LoopInd) / (domain.ThicknessSnowSoilLayer(LoopInd) * ConstDensityIce));
    water.state.SnowEffPorosity(LoopInd) = real(1.0) - water.state.SnowIceVol(LoopInd);
    water.state.SnowLiqWaterVol(LoopInd) =
        std::min(water.state.SnowEffPorosity(LoopInd),
                 water.state.SnowLiqWater(LoopInd) / (domain.ThicknessSnowSoilLayer(LoopInd) * ConstDensityWater));
  }

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
    SnowDensBulk(LoopInd) =
        (water.state.SnowIce(LoopInd) + water.state.SnowLiqWater(LoopInd)) / domain.ThicknessSnowSoilLayer(LoopInd);
    energy.state.HeatCapacVolSnow(LoopInd) =
        ConstHeatCapacIce * water.state.SnowIceVol(LoopInd) +
        ConstHeatCapacWater * water.state.SnowLiqWaterVol(LoopInd);
  }

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
    if (nml.OptSnowThermConduct == 1) {
      energy.state.ThermConductSnow(LoopInd) = real(3.2217e-6) * (SnowDensBulk(LoopInd) * SnowDensBulk(LoopInd));
    }
    if (nml.OptSnowThermConduct == 2) {
      energy.state.ThermConductSnow(LoopInd) =
          real(2e-2) + real(2.5e-6) * SnowDensBulk(LoopInd) * SnowDensBulk(LoopInd);
    }
    if (nml.OptSnowThermConduct == 3) {
      energy.state.ThermConductSnow(LoopInd) = real(0.35);
    }
    if (nml.OptSnowThermConduct == 4) {
      energy.state.ThermConductSnow(LoopInd) =
          real(2.576e-6) * (SnowDensBulk(LoopInd) * SnowDensBulk(LoopInd)) + real(0.074);
    }
    if (nml.OptSnowThermConduct == 5) {
      energy.state.ThermConductSnow(LoopInd) =
          real(2.22) * rpow(SnowDensBulk(LoopInd) / real(1000.0), real(1.88));
    }
  }
}

// ============================================================================
// SoilThermalPropertyMod.F90
// ============================================================================
CPPWRF_HD inline void SoilThermalProperty(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& water = noahmp.water;
  auto& energy = noahmp.energy;

  for (int LoopInd = 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    real SoilIceTmp = water.state.SoilMoisture(LoopInd) - water.state.SoilLiqWater(LoopInd);
    energy.state.HeatCapacVolSoil(LoopInd) =
        water.state.SoilLiqWater(LoopInd) * ConstHeatCapacWater +
        (real(1.0) - water.param.SoilMoistureSat(LoopInd)) * energy.param.SoilHeatCapacity +
        (water.param.SoilMoistureSat(LoopInd) - water.state.SoilMoisture(LoopInd)) * ConstHeatCapacAir +
        SoilIceTmp * ConstHeatCapacIce;

    real SoilSatRatio = water.state.SoilMoisture(LoopInd) / water.param.SoilMoistureSat(LoopInd);

    real ThermConductSolid =
        rpow(ConstThermConductQuartz, energy.param.SoilQuartzFrac(LoopInd)) *
        rpow(ConstThermConductSoilOth, real(1.0) - energy.param.SoilQuartzFrac(LoopInd));

    real SoilWatFrac = real(1.0);
    if (water.state.SoilMoisture(LoopInd) > real(0.0)) {
      SoilWatFrac = water.state.SoilLiqWater(LoopInd) / water.state.SoilMoisture(LoopInd);
    }
    real SoilWatFracSat = SoilWatFrac * water.param.SoilMoistureSat(LoopInd);

    real ThermConductSoilSat =
        rpow(ThermConductSolid, real(1.0) - water.param.SoilMoistureSat(LoopInd)) *
        rpow(ConstThermConductIce, water.param.SoilMoistureSat(LoopInd) - SoilWatFracSat) *
        rpow(ConstThermConductWater, SoilWatFracSat);

    real SoilGamFac = (real(1.0) - water.param.SoilMoistureSat(LoopInd)) * real(2700.0);
    real ThermConductSoilDry = (real(0.135) * SoilGamFac + real(64.7)) / (real(2700.0) - real(0.947) * SoilGamFac);

    real KerstenFac = real(0.0);
    if ((water.state.SoilLiqWater(LoopInd) + real(0.0005)) < water.state.SoilMoisture(LoopInd)) {
      KerstenFac = SoilSatRatio;
    } else {
      if (SoilSatRatio > real(0.1)) {
        KerstenFac = std::log10(SoilSatRatio) + real(1.0);
      } else {
        KerstenFac = real(0.0);
      }
    }

    energy.state.ThermConductSoil(LoopInd) =
        KerstenFac * (ThermConductSoilSat - ThermConductSoilDry) + ThermConductSoilDry;
  }
}

// ============================================================================
// GroundThermalPropertyMod.F90
// ============================================================================
CPPWRF_HD inline void GroundThermalProperty(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& water = noahmp.water;
  auto& energy = noahmp.energy;

  energy.state.HeatCapacSoilSnow.fill(real(0.0));
  energy.state.ThermConductSoilSnow.fill(real(0.0));

  SnowThermalProperty(noahmp);
  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
    energy.state.ThermConductSoilSnow(LoopInd) = energy.state.ThermConductSnow(LoopInd);
    energy.state.HeatCapacSoilSnow(LoopInd) = energy.state.HeatCapacVolSnow(LoopInd);
  }

  SoilThermalProperty(noahmp);
  for (int LoopInd = 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    energy.state.ThermConductSoilSnow(LoopInd) = energy.state.ThermConductSoil(LoopInd);
    energy.state.HeatCapacSoilSnow(LoopInd) = energy.state.HeatCapacVolSoil(LoopInd);
  }
  if (domain.FlagUrban) {
    for (int LoopInd = 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
      energy.state.ThermConductSoilSnow(LoopInd) = real(3.24);
    }
  }

  if (domain.SurfaceType == 2) {
    for (int LoopInd = 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
      if (energy.state.TemperatureSoilSnow(LoopInd) > ConstFreezePoint) {
        energy.state.HeatCapacSoilSnow(LoopInd) = ConstHeatCapacWater;
        energy.state.ThermConductSoilSnow(LoopInd) = ConstThermConductWater;
      } else {
        energy.state.HeatCapacSoilSnow(LoopInd) = ConstHeatCapacIce;
        energy.state.ThermConductSoilSnow(LoopInd) = ConstThermConductIce;
      }
    }
  }

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    energy.state.PhaseChgFacSoilSnow(LoopInd) =
        domain.MainTimeStep / (energy.state.HeatCapacSoilSnow(LoopInd) * domain.ThicknessSnowSoilLayer(LoopInd));
  }

  if (domain.NumSnowLayerNeg == 0) {
    energy.state.ThermConductSoilSnow(1) =
        (energy.state.ThermConductSoilSnow(1) * domain.ThicknessSnowSoilLayer(1) +
         real(0.35) * water.state.SnowDepth) /
        (water.state.SnowDepth + domain.ThicknessSnowSoilLayer(1));
  } else {
    energy.state.ThermConductSoilSnow(1) =
        (energy.state.ThermConductSoilSnow(1) * domain.ThicknessSnowSoilLayer(1) +
         energy.state.ThermConductSoilSnow(0) * domain.ThicknessSnowSoilLayer(0)) /
        (domain.ThicknessSnowSoilLayer(0) + domain.ThicknessSnowSoilLayer(1));
  }
}

// ============================================================================
// SoilWaterTranspirationMod.F90
// ============================================================================
CPPWRF_HD inline void SoilWaterTranspiration(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& nml = noahmp.config.nmlist;
  auto& water = noahmp.water;

  real MinThr = real(1.0e-6);
  water.state.SoilTranspFacAcc = real(0.0);

  if (domain.SurfaceType == 1) {
    for (int IndSoil = 1; IndSoil <= water.param.NumSoilLayerRoot; ++IndSoil) {
      real SoilWetFac = real(0.0);
      if (nml.OptSoilWaterTranspiration == 1) {
        SoilWetFac = (water.state.SoilLiqWater(IndSoil) - water.param.SoilMoistureWilt(IndSoil)) /
                     (water.param.SoilMoistureFieldCap(IndSoil) - water.param.SoilMoistureWilt(IndSoil));
      }
      if (nml.OptSoilWaterTranspiration == 2) {
        water.state.SoilMatPotential(IndSoil) =
            std::max(water.param.SoilMatPotentialWilt,
                     -water.param.SoilMatPotentialSat(IndSoil) *
                         rpow(std::max(real(0.01), water.state.SoilLiqWater(IndSoil)) /
                                  water.param.SoilMoistureSat(IndSoil),
                              -water.param.SoilExpCoeffB(IndSoil)));
        SoilWetFac = (real(1.0) - water.state.SoilMatPotential(IndSoil) / water.param.SoilMatPotentialWilt) /
                     (real(1.0) + water.param.SoilMatPotentialSat(IndSoil) / water.param.SoilMatPotentialWilt);
      }
      if (nml.OptSoilWaterTranspiration == 3) {
        water.state.SoilMatPotential(IndSoil) =
            std::max(water.param.SoilMatPotentialWilt,
                     -water.param.SoilMatPotentialSat(IndSoil) *
                         rpow(std::max(real(0.01), water.state.SoilLiqWater(IndSoil)) /
                                  water.param.SoilMoistureSat(IndSoil),
                              -water.param.SoilExpCoeffB(IndSoil)));
        SoilWetFac = real(1.0) -
                     rexp(real(-5.8) * rlog(water.param.SoilMatPotentialWilt / water.state.SoilMatPotential(IndSoil)));
      }
      SoilWetFac = std::min(real(1.0), std::max(real(0.0), SoilWetFac));

      water.state.SoilTranspFac(IndSoil) =
          std::max(MinThr,
                   domain.ThicknessSnowSoilLayer(IndSoil) /
                       (-domain.DepthSoilLayer(water.param.NumSoilLayerRoot)) * SoilWetFac);
      water.state.SoilTranspFacAcc = water.state.SoilTranspFacAcc + water.state.SoilTranspFac(IndSoil);
    }

    water.state.SoilTranspFacAcc = std::max(MinThr, water.state.SoilTranspFacAcc);
    for (int IndSoil = 1; IndSoil <= water.param.NumSoilLayerRoot; ++IndSoil) {
      water.state.SoilTranspFac(IndSoil) = water.state.SoilTranspFac(IndSoil) / water.state.SoilTranspFacAcc;
    }
  }
}

// ============================================================================
// ResistanceGroundEvaporationMod.F90
// ============================================================================
CPPWRF_HD inline void ResistanceGroundEvaporation(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& nml = noahmp.config.nmlist;
  auto& water = noahmp.water;
  auto& energy = noahmp.energy;

  real SoilEvapFac = std::max(real(0.0), water.state.SoilLiqWater(1) / water.param.SoilMoistureSat(1));

  if (domain.SurfaceType == 2) {  // lake point
    energy.state.ResistanceGrdEvap = real(1.0);
    energy.state.RelHumidityGrd = real(1.0);
  } else {  // soil point
    if (nml.OptGroundResistanceEvap == 1 || nml.OptGroundResistanceEvap == 4) {
      real DrySoilThickness =
          (-domain.DepthSoilLayer(1)) *
          (rexp(rpow(real(1.0) - std::min(real(1.0), water.state.SoilLiqWater(1) / water.param.SoilMoistureSat(1)),
                     energy.param.ResistanceSoilExp)) -
           real(1.0)) /
          (real(2.71828) - real(1.0));
      real VapDiffuseRed =
          real(2.2e-5) * water.param.SoilMoistureSat(1) * water.param.SoilMoistureSat(1) *
          rpow(real(1.0) - water.param.SoilMoistureWilt(1) / water.param.SoilMoistureSat(1),
               real(2.0) + real(3.0) / water.param.SoilExpCoeffB(1));
      energy.state.ResistanceGrdEvap = DrySoilThickness / VapDiffuseRed;
    } else if (nml.OptGroundResistanceEvap == 2) {
      energy.state.ResistanceGrdEvap =
          water.state.SnowCoverFrac * real(1.0) +
          (real(1.0) - water.state.SnowCoverFrac) * rexp(real(8.25) - real(4.225) * SoilEvapFac);
    } else if (nml.OptGroundResistanceEvap == 3) {
      energy.state.ResistanceGrdEvap =
          water.state.SnowCoverFrac * real(1.0) +
          (real(1.0) - water.state.SnowCoverFrac) * rexp(real(8.25) - real(6.0) * SoilEvapFac);
    }

    if (nml.OptGroundResistanceEvap == 4) {
      energy.state.ResistanceGrdEvap =
          real(1.0) / (water.state.SnowCoverFrac * (real(1.0) / energy.param.ResistanceSnowSfc) +
                       (real(1.0) - water.state.SnowCoverFrac) *
                           (real(1.0) / std::max(energy.state.ResistanceGrdEvap, real(0.001))));
    }
    if (water.state.SoilLiqWater(1) < real(0.01) && water.state.SnowDepth == real(0.0)) {
      energy.state.ResistanceGrdEvap = real(1.0e6);
    }

    real SoilMatPotentialSfc =
        -water.param.SoilMatPotentialSat(1) *
        rpow(std::max(real(0.01), water.state.SoilLiqWater(1)) / water.param.SoilMoistureSat(1),
             -water.param.SoilExpCoeffB(1));
    energy.state.RelHumidityGrd =
        water.state.SnowCoverFrac +
        (real(1.0) - water.state.SnowCoverFrac) *
            rexp(SoilMatPotentialSfc * ConstGravityAcc / (ConstGasWaterVapor * energy.state.TemperatureGrd));
  }

  if (domain.FlagUrban && water.state.SnowDepth == real(0.0)) {
    energy.state.ResistanceGrdEvap = real(1.0e6);
  }
}

// ============================================================================
// ResistanceLeafToGroundMod.F90
// ============================================================================
CPPWRF_HD inline void ResistanceLeafToGround(NoahmpType& noahmp,
                                             int IndIter,
                                             real VegAreaIndEff,
                                             real HeatSenGrdTmp) {
  auto& energy = noahmp.energy;

  real MPE = real(1.0e-6);
  energy.state.MoStabParaUndCan = real(0.0);
  energy.state.MoLengthUndCan = real(0.0);

  if (IndIter > 1) {
    real TMP1 = ConstVonKarman * (ConstGravityAcc / energy.state.TemperatureCanopyAir) * HeatSenGrdTmp /
                (energy.state.DensityAirRefHeight * ConstHeatCapacAir);
    if (std::abs(TMP1) <= MPE) TMP1 = MPE;
    energy.state.MoLengthUndCan =
        real(-1.0) * (energy.state.FrictionVelVeg * energy.state.FrictionVelVeg * energy.state.FrictionVelVeg) / TMP1;
    energy.state.MoStabParaUndCan =
        std::min((energy.state.ZeroPlaneDispSfc - energy.state.RoughLenMomGrd) / energy.state.MoLengthUndCan,
                 real(1.0));
  }
  real FHGNEW = real(0.0);
  if (energy.state.MoStabParaUndCan < real(0.0)) {
    FHGNEW = rpow(real(1.0) - real(15.0) * energy.state.MoStabParaUndCan, real(-0.25));
  } else {
    FHGNEW = real(1.0) + real(4.7) * energy.state.MoStabParaUndCan;
  }
  if (IndIter == 1) {
    energy.state.MoStabCorrShUndCan = FHGNEW;
  } else {
    energy.state.MoStabCorrShUndCan = real(0.5) * (energy.state.MoStabCorrShUndCan + FHGNEW);
  }

  energy.state.WindExtCoeffCanopy =
      rpow(energy.param.CanopyWindExtFac * VegAreaIndEff * energy.state.CanopyHeight * energy.state.MoStabCorrShUndCan,
           real(0.5));
  real TMP1 = rexp(-energy.state.WindExtCoeffCanopy * energy.state.RoughLenShVegGrd / energy.state.CanopyHeight);
  real TMP2 = rexp(-energy.state.WindExtCoeffCanopy *
                   (energy.state.RoughLenShCanopy + energy.state.ZeroPlaneDispSfc) / energy.state.CanopyHeight);
  real TMPRAH2 =
      energy.state.CanopyHeight * rexp(energy.state.WindExtCoeffCanopy) / energy.state.WindExtCoeffCanopy *
      (TMP1 - TMP2);

  real KH = std::max(ConstVonKarman * energy.state.FrictionVelVeg *
                         (energy.state.CanopyHeight - energy.state.ZeroPlaneDispSfc),
                     MPE);
  energy.state.ResistanceMomUndCan = real(0.0);
  energy.state.ResistanceShUndCan = TMPRAH2 / KH;
  energy.state.ResistanceLhUndCan = energy.state.ResistanceShUndCan;

  real TMPRB = energy.state.WindExtCoeffCanopy * real(50.0) /
               (real(1.0) - rexp(-energy.state.WindExtCoeffCanopy / real(2.0)));
  energy.state.ResistanceLeafBoundary =
      TMPRB * std::sqrt(energy.param.LeafDimLength / energy.state.WindSpdCanopyTop);
  energy.state.ResistanceLeafBoundary =
      std::min(std::max(energy.state.ResistanceLeafBoundary, real(5.0)), real(50.0));
}

}  // namespace cppwrf::phys::noahmp
