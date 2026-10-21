#pragma once

#include <algorithm>
#include <cmath>
#include "core/backend.hpp"
#include "core/precision.hpp"
#include "phys/lsm_noahmp/noahmp_energy.hpp"
#include "phys/lsm_noahmp/noahmp_rad_therm.hpp"
#include "phys/lsm_noahmp/noahmp_var.hpp"
#include "phys/lsm_noahmp/noahmp_water.hpp"

namespace cppwrf::phys::noahmp {

// ============================================================================
// GeneralInitGlacierMod.F90
// ============================================================================
CPPWRF_HD inline void GeneralInitGlacier(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    if (LoopInd == (domain.NumSnowLayerNeg + 1)) {
      domain.ThicknessSnowSoilLayer(LoopInd) = -domain.DepthSnowSoilLayer(LoopInd);
    } else {
      domain.ThicknessSnowSoilLayer(LoopInd) =
          domain.DepthSnowSoilLayer(LoopInd - 1) - domain.DepthSnowSoilLayer(LoopInd);
    }
  }
}

// ============================================================================
// PrecipitationHeatAdvectGlacierMod.F90
// ============================================================================
CPPWRF_HD inline void PrecipitationHeatAdvectGlacier(NoahmpType& noahmp) {
  auto& forcing = noahmp.forcing;
  auto& energy = noahmp.energy;
  auto& water = noahmp.water;

  real HeatPrcpAirToGrd = real(0.0);
  energy.flux.HeatPrecipAdvBareGrd = real(0.0);
  water.flux.RainfallGround = water.flux.RainfallRefHeight;
  water.flux.SnowfallGround = water.flux.SnowfallRefHeight;

  HeatPrcpAirToGrd = water.flux.RainfallGround * (ConstHeatCapacWater / real(1000.0)) *
                     (forcing.TemperatureAirRefHeight - energy.state.TemperatureGrd);
  HeatPrcpAirToGrd = HeatPrcpAirToGrd +
                     water.flux.SnowfallGround * (ConstHeatCapacIce / real(1000.0)) *
                         (forcing.TemperatureAirRefHeight - energy.state.TemperatureGrd);

  energy.flux.HeatPrecipAdvBareGrd = HeatPrcpAirToGrd;
  energy.flux.HeatPrecipAdvBareGrd = std::max(energy.flux.HeatPrecipAdvBareGrd, real(-20.0));
  energy.flux.HeatPrecipAdvBareGrd = std::min(energy.flux.HeatPrecipAdvBareGrd, real(20.0));
}

// ============================================================================
// SnowCoverGlacierMod.F90
// ============================================================================
CPPWRF_HD inline void SnowCoverGlacier(NoahmpType& noahmp) {
  auto& water = noahmp.water;
  water.state.SnowCoverFrac = real(0.0);
  if (water.state.SnowWaterEquiv > real(0.0)) {
    water.state.SnowCoverFrac = real(1.0);
  }
}

// ============================================================================
// GroundRoughnessPropertyGlacierMod.F90
// ============================================================================
CPPWRF_HD inline void GroundRoughnessPropertyGlacier(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& energy = noahmp.energy;
  auto& water = noahmp.water;

  energy.state.RoughLenMomGrd = energy.param.RoughLenMomSnow;
  energy.state.RoughLenMomSfc = energy.state.RoughLenMomGrd;

  energy.state.ZeroPlaneDispGrd = water.state.SnowDepth;
  energy.state.ZeroPlaneDispSfc = energy.state.ZeroPlaneDispGrd;

  energy.state.RefHeightAboveGrd = energy.state.ZeroPlaneDispSfc + domain.RefHeightAboveSfc;
}

// ============================================================================
// GlacierIceThermalPropertyMod.F90
// ============================================================================
CPPWRF_HD inline void GlacierIceThermalProperty(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& energy = noahmp.energy;

  for (int LoopInd1 = 1; LoopInd1 <= domain.NumSoilLayer; ++LoopInd1) {
    real DepthIceLayerMid = real(0.5) * domain.ThicknessSnowSoilLayer(LoopInd1);
    for (int LoopInd2 = 1; LoopInd2 <= LoopInd1 - 1; ++LoopInd2) {
      DepthIceLayerMid = DepthIceLayerMid + domain.ThicknessSnowSoilLayer(LoopInd2);
    }
    energy.state.HeatCapacGlaIce(LoopInd1) = real(1.0e6) * (real(0.8194) + real(0.1309) * DepthIceLayerMid);
    energy.state.ThermConductGlaIce(LoopInd1) = real(0.32333) + (real(0.10073) * DepthIceLayerMid);
  }
}

// ============================================================================
// GroundThermalPropertyGlacierMod.F90
// ============================================================================
CPPWRF_HD inline void GroundThermalPropertyGlacier(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& energy = noahmp.energy;
  auto& water = noahmp.water;

  energy.state.HeatCapacSoilSnow.fill(real(0.0));
  energy.state.ThermConductSoilSnow.fill(real(0.0));

  SnowThermalProperty(noahmp);
  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
    energy.state.ThermConductSoilSnow(LoopInd) = energy.state.ThermConductSnow(LoopInd);
    energy.state.HeatCapacSoilSnow(LoopInd) = energy.state.HeatCapacVolSnow(LoopInd);
  }

  GlacierIceThermalProperty(noahmp);
  for (int LoopInd = 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    energy.state.ThermConductSoilSnow(LoopInd) = energy.state.ThermConductGlaIce(LoopInd);
    energy.state.HeatCapacSoilSnow(LoopInd) = energy.state.HeatCapacGlaIce(LoopInd);
  }

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    energy.state.PhaseChgFacSoilSnow(LoopInd) =
        domain.MainTimeStep / (energy.state.HeatCapacSoilSnow(LoopInd) * domain.ThicknessSnowSoilLayer(LoopInd));
  }

  if (domain.NumSnowLayerNeg == 0) {
    energy.state.ThermConductSoilSnow(1) =
        (energy.state.ThermConductSoilSnow(1) * domain.ThicknessSnowSoilLayer(1) + real(0.35) * water.state.SnowDepth) /
        (water.state.SnowDepth + domain.ThicknessSnowSoilLayer(1));
  } else {
    energy.state.ThermConductSoilSnow(1) =
        (energy.state.ThermConductSoilSnow(1) * domain.ThicknessSnowSoilLayer(1) +
         energy.state.ThermConductSoilSnow(0) * domain.ThicknessSnowSoilLayer(0)) /
        (domain.ThicknessSnowSoilLayer(0) + domain.ThicknessSnowSoilLayer(1));
  }
}

// ============================================================================
// GroundAlbedoGlacierMod.F90
// ============================================================================
CPPWRF_HD inline void GroundAlbedoGlacier(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& energy = noahmp.energy;
  auto& water = noahmp.water;

  for (int IndSwBnd = 1; IndSwBnd <= domain.NumSwRadBand; ++IndSwBnd) {
    energy.state.AlbedoGrdDir(IndSwBnd) =
        energy.param.AlbedoLandIce(IndSwBnd) * (real(1.0) - water.state.SnowCoverFrac) +
        energy.state.AlbedoSnowDir(IndSwBnd) * water.state.SnowCoverFrac;
    energy.state.AlbedoGrdDif(IndSwBnd) =
        energy.param.AlbedoLandIce(IndSwBnd) * (real(1.0) - water.state.SnowCoverFrac) +
        energy.state.AlbedoSnowDif(IndSwBnd) * water.state.SnowCoverFrac;
  }
}

// ============================================================================
// SurfaceAlbedoGlacierMod.F90
// ============================================================================
CPPWRF_HD inline void SurfaceAlbedoGlacier(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& nmlist = noahmp.config.nmlist;
  auto& energy = noahmp.energy;

  for (int IndBand = 1; IndBand <= domain.NumSwRadBand; ++IndBand) {
    energy.state.AlbedoSfcDir(IndBand) = real(0.0);
    energy.state.AlbedoSfcDif(IndBand) = real(0.0);
    energy.state.AlbedoGrdDir(IndBand) = real(0.0);
    energy.state.AlbedoGrdDif(IndBand) = real(0.0);
    energy.state.AlbedoSnowDir(IndBand) = real(0.0);
    energy.state.AlbedoSnowDif(IndBand) = real(0.0);
  }

  if (nmlist.OptSnowAlbedo == 1) {
    SnowAgingBats(noahmp);
  }

  if (domain.CosSolarZenithAngle > real(0.0)) {
    if (nmlist.OptSnowAlbedo == 1) SnowAlbedoBats(noahmp);
    if (nmlist.OptSnowAlbedo == 2) SnowAlbedoClass(noahmp);

    GroundAlbedoGlacier(noahmp);

    for (int IndBand = 1; IndBand <= domain.NumSwRadBand; ++IndBand) {
      energy.state.AlbedoSfcDir(IndBand) = energy.state.AlbedoGrdDir(IndBand);
      energy.state.AlbedoSfcDif(IndBand) = energy.state.AlbedoGrdDif(IndBand);
    }
  }
}

// ============================================================================
// SurfaceRadiationGlacierMod.F90
// ============================================================================
CPPWRF_HD inline void SurfaceRadiationGlacier(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& energy = noahmp.energy;

  energy.flux.RadSwAbsGrd = real(0.0);
  energy.flux.RadSwAbsSfc = real(0.0);
  energy.flux.RadSwReflSfc = real(0.0);

  for (int IndBand = 1; IndBand <= domain.NumSwRadBand; ++IndBand) {
    real RadSwAbsGrdTmp = energy.flux.RadSwDownDir(IndBand) * (real(1.0) - energy.state.AlbedoGrdDir(IndBand)) +
                          energy.flux.RadSwDownDif(IndBand) * (real(1.0) - energy.state.AlbedoGrdDif(IndBand));
    energy.flux.RadSwAbsGrd = energy.flux.RadSwAbsGrd + RadSwAbsGrdTmp;
    energy.flux.RadSwAbsSfc = energy.flux.RadSwAbsSfc + RadSwAbsGrdTmp;

    real RadSwReflGrdTmp = energy.flux.RadSwDownDir(IndBand) * energy.state.AlbedoGrdDir(IndBand) +
                           energy.flux.RadSwDownDif(IndBand) * energy.state.AlbedoGrdDif(IndBand);
    energy.flux.RadSwReflSfc = energy.flux.RadSwReflSfc + RadSwReflGrdTmp;
  }
}

// ============================================================================
// SurfaceEmissivityGlacierMod.F90
// ============================================================================
CPPWRF_HD inline void SurfaceEmissivityGlacier(NoahmpType& noahmp) {
  auto& energy = noahmp.energy;
  auto& water = noahmp.water;

  energy.state.EmissivityGrd = energy.param.EmissivityIceSfc * (real(1.0) - water.state.SnowCoverFrac) +
                               energy.param.EmissivitySnow * water.state.SnowCoverFrac;
  energy.state.EmissivitySfc = energy.state.EmissivityGrd;
}

// ============================================================================
// ResistanceGroundEvaporationGlacierMod.F90
// ============================================================================
CPPWRF_HD inline void ResistanceGroundEvaporationGlacier(NoahmpType& noahmp) {
  auto& energy = noahmp.energy;
  energy.state.ResistanceGrdEvap = real(1.0);
  energy.state.RelHumidityGrd = real(1.0);
}

// ============================================================================
// PsychrometricVariableGlacierMod.F90
// ============================================================================
CPPWRF_HD inline void PsychrometricVariableGlacier(NoahmpType& noahmp) {
  auto& forcing = noahmp.forcing;
  auto& energy = noahmp.energy;

  energy.state.LatHeatVapGrd = ConstLatHeatSublim;
  energy.state.PsychConstGrd =
      ConstHeatCapacAir * forcing.PressureAirRefHeight / (real(0.622) * energy.state.LatHeatVapGrd);
}

// ============================================================================
// SurfaceEnergyFluxGlacierMod.F90
// ============================================================================
CPPWRF_HD inline void SurfaceEnergyFluxGlacier(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& nmlist = noahmp.config.nmlist;
  auto& forcing = noahmp.forcing;
  auto& energy = noahmp.energy;
  auto& water = noahmp.water;

  constexpr int NumIter = 5;
  FArray1D<real, 1, kMaxSoilLayers> SoilIceTmp;
  SoilIceTmp.fill(real(0.0));

  real TemperatureGrdChg = real(0.0);
  energy.state.MoStabParaBare = real(0.0);
  int MoStabParaSgn = 0;
  energy.state.MoStabCorrShBare2m = real(0.0);
  real HeatSensibleTmp = real(0.0);
  energy.state.FrictionVelBare = real(0.1);
  real LwRadCoeff = energy.state.EmissivityGrd * ConstStefanBoltzmann;
  real GrdHeatCoeff = real(2.0) * energy.state.ThermConductSoilSnow(domain.NumSnowLayerNeg + 1) /
                      domain.ThicknessSnowSoilLayer(domain.NumSnowLayerNeg + 1);

  real ExchCoeffShTmp = real(0.0);
  real ShCoeff = real(0.0);
  real LhCoeff = real(0.0);

  auto TempUnitConv = [](real TempTmp) {
    return std::min(real(50.0), std::max(real(-50.0), (TempTmp - ConstFreezePoint)));
  };

  for (int IndIter = 1; IndIter <= NumIter; ++IndIter) {
    energy.state.RoughLenShBareGrd = energy.state.RoughLenMomGrd;

    ResistanceBareGroundMOST(noahmp, IndIter, HeatSensibleTmp, MoStabParaSgn);

    ExchCoeffShTmp = real(1.0) / energy.state.ResistanceShBareGrd;

    real TempTmp = TempUnitConv(energy.state.TemperatureGrdBare);
    real VapPresSatWatTmp = real(0.0);
    real VapPresSatIceTmp = real(0.0);
    real VapPresSatWatTmpD = real(0.0);
    real VapPresSatIceTmpD = real(0.0);
    VaporPressureSaturation(TempTmp, VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD);
    if (TempTmp > real(0.0)) {
      energy.state.VapPresSatGrdBare = VapPresSatWatTmp;
      energy.state.VapPresSatGrdBareTempD = VapPresSatWatTmpD;
    } else {
      energy.state.VapPresSatGrdBare = VapPresSatIceTmp;
      energy.state.VapPresSatGrdBareTempD = VapPresSatIceTmpD;
    }

    ShCoeff = energy.state.DensityAirRefHeight * ConstHeatCapacAir / energy.state.ResistanceShBareGrd;
    if ((water.state.SnowDepth > real(0.0)) || (nmlist.OptGlacierTreatment == 1)) {
      LhCoeff = energy.state.DensityAirRefHeight * ConstHeatCapacAir / energy.state.PsychConstGrd /
                (energy.state.ResistanceGrdEvap + energy.state.ResistanceLhBareGrd);
    } else {
      LhCoeff = real(0.0);
    }
    energy.flux.RadLwNetBareGrd =
        LwRadCoeff * cppwrf::powi(energy.state.TemperatureGrdBare, 4) -
        energy.state.EmissivityGrd * forcing.RadLwDownRefHeight;
    energy.flux.HeatSensibleBareGrd =
        ShCoeff * (energy.state.TemperatureGrdBare - forcing.TemperatureAirRefHeight);
    energy.flux.HeatLatentBareGrd =
        LhCoeff * (energy.state.VapPresSatGrdBare * energy.state.RelHumidityGrd - energy.state.PressureVaporRefHeight);
    energy.flux.HeatGroundBareGrd =
        GrdHeatCoeff * (energy.state.TemperatureGrdBare - energy.state.TemperatureSoilSnow(domain.NumSnowLayerNeg + 1));
    real EnergyResTmp = energy.flux.RadSwAbsGrd - energy.flux.RadLwNetBareGrd - energy.flux.HeatSensibleBareGrd -
                        energy.flux.HeatLatentBareGrd - energy.flux.HeatGroundBareGrd + energy.flux.HeatPrecipAdvBareGrd;
    real FluxTotCoeff = real(4.0) * LwRadCoeff * cppwrf::powi(energy.state.TemperatureGrdBare, 3) + ShCoeff +
                        LhCoeff * energy.state.VapPresSatGrdBareTempD + GrdHeatCoeff;
    TemperatureGrdChg = EnergyResTmp / FluxTotCoeff;
    energy.flux.RadLwNetBareGrd =
        energy.flux.RadLwNetBareGrd +
        real(4.0) * LwRadCoeff * cppwrf::powi(energy.state.TemperatureGrdBare, 3) * TemperatureGrdChg;
    energy.flux.HeatSensibleBareGrd = energy.flux.HeatSensibleBareGrd + ShCoeff * TemperatureGrdChg;
    energy.flux.HeatLatentBareGrd =
        energy.flux.HeatLatentBareGrd + LhCoeff * energy.state.VapPresSatGrdBareTempD * TemperatureGrdChg;
    energy.flux.HeatGroundBareGrd = energy.flux.HeatGroundBareGrd + GrdHeatCoeff * TemperatureGrdChg;
    energy.state.TemperatureGrdBare = energy.state.TemperatureGrdBare + TemperatureGrdChg;

    HeatSensibleTmp = ShCoeff * (energy.state.TemperatureGrdBare - forcing.TemperatureAirRefHeight);

    TempTmp = TempUnitConv(energy.state.TemperatureGrdBare);
    VaporPressureSaturation(TempTmp, VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD);
    if (TempTmp > real(0.0)) {
      energy.state.VapPresSatGrdBare = VapPresSatWatTmp;
    } else {
      energy.state.VapPresSatGrdBare = VapPresSatIceTmp;
    }
    energy.state.SpecHumiditySfc =
        real(0.622) * (energy.state.VapPresSatGrdBare * energy.state.RelHumidityGrd) /
        (forcing.PressureAirRefHeight - real(0.378) * (energy.state.VapPresSatGrdBare * energy.state.RelHumidityGrd));
  }

  real max_soil_ice_tmp = real(-1.0e30);
  for (int k = 1; k <= domain.NumSoilLayer; ++k) {
    SoilIceTmp(k) = water.state.SoilMoisture(k) - water.state.SoilLiqWater(k);
    if (SoilIceTmp(k) > max_soil_ice_tmp) max_soil_ice_tmp = SoilIceTmp(k);
  }
  if ((nmlist.OptSnowSoilTempTime == 1) || (nmlist.OptSnowSoilTempTime == 3)) {
    if ((max_soil_ice_tmp > real(0.0) || water.state.SnowDepth > real(0.05)) &&
        (energy.state.TemperatureGrdBare > ConstFreezePoint) && (nmlist.OptGlacierTreatment == 1)) {
      energy.state.TemperatureGrdBare = ConstFreezePoint;
      real TempTmp = TempUnitConv(energy.state.TemperatureGrdBare);
      real VapPresSatWatTmp = real(0.0);
      real VapPresSatIceTmp = real(0.0);
      real VapPresSatWatTmpD = real(0.0);
      real VapPresSatIceTmpD = real(0.0);
      VaporPressureSaturation(TempTmp, VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD);
      energy.state.VapPresSatGrdBare = VapPresSatIceTmp;
      energy.state.SpecHumiditySfc =
          real(0.622) * (energy.state.VapPresSatGrdBare * energy.state.RelHumidityGrd) /
          (forcing.PressureAirRefHeight - real(0.378) * (energy.state.VapPresSatGrdBare * energy.state.RelHumidityGrd));
      energy.flux.RadLwNetBareGrd =
          LwRadCoeff * cppwrf::powi(energy.state.TemperatureGrdBare, 4) -
          energy.state.EmissivityGrd * forcing.RadLwDownRefHeight;
      energy.flux.HeatSensibleBareGrd =
          ShCoeff * (energy.state.TemperatureGrdBare - forcing.TemperatureAirRefHeight);
      energy.flux.HeatLatentBareGrd =
          LhCoeff * (energy.state.VapPresSatGrdBare * energy.state.RelHumidityGrd - energy.state.PressureVaporRefHeight);
      energy.flux.HeatGroundBareGrd =
          energy.flux.RadSwAbsGrd + energy.flux.HeatPrecipAdvBareGrd -
          (energy.flux.RadLwNetBareGrd + energy.flux.HeatSensibleBareGrd + energy.flux.HeatLatentBareGrd);
    }
  }

  energy.state.WindStressEwBare = -energy.state.DensityAirRefHeight * energy.state.ExchCoeffMomBare *
                                  energy.state.WindSpdRefHeight * forcing.WindEastwardRefHeight;
  energy.state.WindStressNsBare = -energy.state.DensityAirRefHeight * energy.state.ExchCoeffMomBare *
                                  energy.state.WindSpdRefHeight * forcing.WindNorthwardRefHeight;

  energy.state.ExchCoeffSh2mBare =
      energy.state.FrictionVelBare * ConstVonKarman /
      (rlog((real(2.0) + energy.state.RoughLenShBareGrd) / energy.state.RoughLenShBareGrd) -
       energy.state.MoStabCorrShBare2m);
  if (energy.state.ExchCoeffSh2mBare < real(1.0e-5)) {
    energy.state.TemperatureAir2mBare = energy.state.TemperatureGrdBare;
    energy.state.SpecHumidity2mBare = energy.state.SpecHumiditySfc;
  } else {
    energy.state.TemperatureAir2mBare =
        energy.state.TemperatureGrdBare -
        energy.flux.HeatSensibleBareGrd / (energy.state.DensityAirRefHeight * ConstHeatCapacAir) *
            real(1.0) / energy.state.ExchCoeffSh2mBare;
    energy.state.SpecHumidity2mBare =
        energy.state.SpecHumiditySfc -
        energy.flux.HeatLatentBareGrd / (energy.state.LatHeatVapGrd * energy.state.DensityAirRefHeight) *
            (real(1.0) / energy.state.ExchCoeffSh2mBare + energy.state.ResistanceGrdEvap);
  }

  energy.state.ExchCoeffShBare = ExchCoeffShTmp;
}

// ============================================================================
// GlacierThermalDiffusionMod.F90
// ============================================================================
CPPWRF_HD inline void GlacierThermalDiffusion(
    NoahmpType& noahmp,
    FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers>& MatLeft1,
    FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers>& MatLeft2,
    FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers>& MatLeft3,
    FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers>& MatRight) {
  auto& domain = noahmp.config.domain;
  auto& nmlist = noahmp.config.nmlist;
  auto& forcing = noahmp.forcing;
  auto& energy = noahmp.energy;

  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> DepthSnowSoilInv;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> HeatCapacPerArea;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> TempGradDepth;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> EnergyExcess;

  MatRight.fill(real(0.0));
  MatLeft1.fill(real(0.0));
  MatLeft2.fill(real(0.0));
  MatLeft3.fill(real(0.0));
  DepthSnowSoilInv.fill(real(0.0));
  HeatCapacPerArea.fill(real(0.0));
  TempGradDepth.fill(real(0.0));
  EnergyExcess.fill(real(0.0));

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    if (LoopInd == (domain.NumSnowLayerNeg + 1)) {
      HeatCapacPerArea(LoopInd) = -domain.DepthSnowSoilLayer(LoopInd) * energy.state.HeatCapacSoilSnow(LoopInd);
      real DepthSnowSoilTmp = -domain.DepthSnowSoilLayer(LoopInd + 1);
      DepthSnowSoilInv(LoopInd) = real(2.0) / DepthSnowSoilTmp;
      TempGradDepth(LoopInd) =
          real(2.0) * (energy.state.TemperatureSoilSnow(LoopInd) - energy.state.TemperatureSoilSnow(LoopInd + 1)) /
          DepthSnowSoilTmp;
      EnergyExcess(LoopInd) = energy.state.ThermConductSoilSnow(LoopInd) * TempGradDepth(LoopInd) -
                              energy.flux.HeatGroundTot - energy.flux.RadSwPenetrateGrd(LoopInd);
    } else if (LoopInd < domain.NumSoilLayer) {
      HeatCapacPerArea(LoopInd) =
          (domain.DepthSnowSoilLayer(LoopInd - 1) - domain.DepthSnowSoilLayer(LoopInd)) *
          energy.state.HeatCapacSoilSnow(LoopInd);
      real DepthSnowSoilTmp = domain.DepthSnowSoilLayer(LoopInd - 1) - domain.DepthSnowSoilLayer(LoopInd + 1);
      DepthSnowSoilInv(LoopInd) = real(2.0) / DepthSnowSoilTmp;
      TempGradDepth(LoopInd) =
          real(2.0) * (energy.state.TemperatureSoilSnow(LoopInd) - energy.state.TemperatureSoilSnow(LoopInd + 1)) /
          DepthSnowSoilTmp;
      EnergyExcess(LoopInd) =
          (energy.state.ThermConductSoilSnow(LoopInd) * TempGradDepth(LoopInd) -
           energy.state.ThermConductSoilSnow(LoopInd - 1) * TempGradDepth(LoopInd - 1)) -
          energy.flux.RadSwPenetrateGrd(LoopInd);
    } else if (LoopInd == domain.NumSoilLayer) {
      HeatCapacPerArea(LoopInd) =
          (domain.DepthSnowSoilLayer(LoopInd - 1) - domain.DepthSnowSoilLayer(LoopInd)) *
          energy.state.HeatCapacSoilSnow(LoopInd);
      if (nmlist.OptSoilTemperatureBottom == 1) {
        energy.flux.HeatFromSoilBot = real(0.0);
      }
      if (nmlist.OptSoilTemperatureBottom == 2) {
        TempGradDepth(LoopInd) =
            (energy.state.TemperatureSoilSnow(LoopInd) - forcing.TemperatureSoilBottom) /
            (real(0.5) * (domain.DepthSnowSoilLayer(LoopInd - 1) + domain.DepthSnowSoilLayer(LoopInd)) -
             energy.state.DepthSoilTempBotToSno);
        energy.flux.HeatFromSoilBot = -energy.state.ThermConductSoilSnow(LoopInd) * TempGradDepth(LoopInd);
      }
      EnergyExcess(LoopInd) =
          (-energy.flux.HeatFromSoilBot -
           energy.state.ThermConductSoilSnow(LoopInd - 1) * TempGradDepth(LoopInd - 1)) -
          energy.flux.RadSwPenetrateGrd(LoopInd);
    }
  }

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    if (LoopInd == (domain.NumSnowLayerNeg + 1)) {
      MatLeft1(LoopInd) = real(0.0);
      MatLeft3(LoopInd) =
          -energy.state.ThermConductSoilSnow(LoopInd) * DepthSnowSoilInv(LoopInd) / HeatCapacPerArea(LoopInd);
      if ((nmlist.OptSnowSoilTempTime == 1) || (nmlist.OptSnowSoilTempTime == 3)) {
        MatLeft2(LoopInd) = -MatLeft3(LoopInd);
      }
      if (nmlist.OptSnowSoilTempTime == 2) {
        MatLeft2(LoopInd) =
            -MatLeft3(LoopInd) +
            energy.state.ThermConductSoilSnow(LoopInd) /
                (real(0.5) * domain.DepthSnowSoilLayer(LoopInd) * domain.DepthSnowSoilLayer(LoopInd) *
                 energy.state.HeatCapacSoilSnow(LoopInd));
      }
    } else if (LoopInd < domain.NumSoilLayer) {
      MatLeft1(LoopInd) =
          -energy.state.ThermConductSoilSnow(LoopInd - 1) * DepthSnowSoilInv(LoopInd - 1) / HeatCapacPerArea(LoopInd);
      MatLeft3(LoopInd) =
          -energy.state.ThermConductSoilSnow(LoopInd) * DepthSnowSoilInv(LoopInd) / HeatCapacPerArea(LoopInd);
      MatLeft2(LoopInd) = -(MatLeft1(LoopInd) + MatLeft3(LoopInd));
    } else if (LoopInd == domain.NumSoilLayer) {
      MatLeft1(LoopInd) =
          -energy.state.ThermConductSoilSnow(LoopInd - 1) * DepthSnowSoilInv(LoopInd - 1) / HeatCapacPerArea(LoopInd);
      MatLeft3(LoopInd) = real(0.0);
      MatLeft2(LoopInd) = -(MatLeft1(LoopInd) + MatLeft3(LoopInd));
    }
    MatRight(LoopInd) = EnergyExcess(LoopInd) / (-HeatCapacPerArea(LoopInd));
  }
}

// ============================================================================
// GlacierTemperatureSolverMod.F90
// ============================================================================
CPPWRF_HD inline void GlacierTemperatureSolver(
    NoahmpType& noahmp,
    real TimeStep,
    FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers>& MatLeft1,
    FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers>& MatLeft2,
    FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers>& MatLeft3,
    FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers>& MatRight) {
  auto& domain = noahmp.config.domain;
  auto& energy = noahmp.energy;

  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MatRightTmp;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MatLeft3Tmp;
  MatRightTmp.fill(real(0.0));
  MatLeft3Tmp.fill(real(0.0));

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    MatRight(LoopInd) = MatRight(LoopInd) * TimeStep;
    MatLeft1(LoopInd) = MatLeft1(LoopInd) * TimeStep;
    MatLeft2(LoopInd) = real(1.0) + MatLeft2(LoopInd) * TimeStep;
    MatLeft3(LoopInd) = MatLeft3(LoopInd) * TimeStep;
  }

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    MatRightTmp(LoopInd) = MatRight(LoopInd);
    MatLeft3Tmp(LoopInd) = MatLeft3(LoopInd);
  }

  MatrixSolverTriDiagonal(MatLeft3, MatLeft1, MatLeft2, MatLeft3Tmp, MatRightTmp, MatRight,
                          domain.NumSnowLayerNeg + 1, domain.NumSoilLayer, domain.NumSnowLayerMax);

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    energy.state.TemperatureSoilSnow(LoopInd) =
        energy.state.TemperatureSoilSnow(LoopInd) + MatLeft3(LoopInd);
  }
}

// ============================================================================
// GlacierTemperatureMainMod.F90
// ============================================================================
CPPWRF_HD inline void GlacierTemperatureMain(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& energy = noahmp.energy;
  auto& water = noahmp.water;

  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MatRight;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MatLeft1;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MatLeft2;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MatLeft3;
  MatRight.fill(real(0.0));
  MatLeft1.fill(real(0.0));
  MatLeft2.fill(real(0.0));
  MatLeft3.fill(real(0.0));

  for (int k = -domain.NumSnowLayerMax + 1; k <= domain.NumSoilLayer; ++k) {
    energy.flux.RadSwPenetrateGrd(k) = real(0.0);
  }

  energy.state.DepthSoilTempBotToSno = domain.DepthSoilTempBottom - water.state.SnowDepth;

  GlacierThermalDiffusion(noahmp, MatLeft1, MatLeft2, MatLeft3, MatRight);
  GlacierTemperatureSolver(noahmp, domain.MainTimeStep, MatLeft1, MatLeft2, MatLeft3, MatRight);
}

// ============================================================================
// GlacierPhaseChangeMod.F90
// ============================================================================
CPPWRF_HD inline void GlacierPhaseChange(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& nmlist = noahmp.config.nmlist;
  auto& energy = noahmp.energy;
  auto& water = noahmp.water;

  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> EnergyRes;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> GlacierPhaseChg;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MassWatTotInit;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MassWatIceInit;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MassWatLiqInit;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MassWatIceTmp;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MassWatLiqTmp;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> EnergyResLeft;

  EnergyRes.fill(real(0.0));
  GlacierPhaseChg.fill(real(0.0));
  MassWatTotInit.fill(real(0.0));
  MassWatIceInit.fill(real(0.0));
  MassWatLiqInit.fill(real(0.0));
  MassWatIceTmp.fill(real(0.0));
  MassWatLiqTmp.fill(real(0.0));
  EnergyResLeft.fill(real(0.0));
  water.flux.MeltGroundSnow = real(0.0);
  water.state.PondSfcThinSnwMelt = real(0.0);
  real HeatLhTotPhsChg = real(0.0);

  for (int LoopInd1 = domain.NumSnowLayerNeg + 1; LoopInd1 <= 0; ++LoopInd1) {
    MassWatIceTmp(LoopInd1) = water.state.SnowIce(LoopInd1);
    MassWatLiqTmp(LoopInd1) = water.state.SnowLiqWater(LoopInd1);
  }

  for (int LoopInd1 = domain.NumSnowLayerNeg + 1; LoopInd1 <= 0; ++LoopInd1) {
    water.state.IndexPhaseChange(LoopInd1) = 0;
    EnergyRes(LoopInd1) = real(0.0);
    GlacierPhaseChg(LoopInd1) = real(0.0);
    EnergyResLeft(LoopInd1) = real(0.0);
    MassWatIceInit(LoopInd1) = MassWatIceTmp(LoopInd1);
    MassWatLiqInit(LoopInd1) = MassWatLiqTmp(LoopInd1);
    MassWatTotInit(LoopInd1) = MassWatIceTmp(LoopInd1) + MassWatLiqTmp(LoopInd1);
  }

  for (int LoopInd1 = domain.NumSnowLayerNeg + 1; LoopInd1 <= 0; ++LoopInd1) {
    if ((MassWatIceTmp(LoopInd1) > real(0.0)) && (energy.state.TemperatureSoilSnow(LoopInd1) >= ConstFreezePoint)) {
      water.state.IndexPhaseChange(LoopInd1) = 1;
    }
    if ((MassWatLiqTmp(LoopInd1) > real(0.0)) && (energy.state.TemperatureSoilSnow(LoopInd1) < ConstFreezePoint)) {
      water.state.IndexPhaseChange(LoopInd1) = 2;
    }
  }

  for (int LoopInd1 = domain.NumSnowLayerNeg + 1; LoopInd1 <= 0; ++LoopInd1) {
    if (water.state.IndexPhaseChange(LoopInd1) > 0) {
      EnergyRes(LoopInd1) = (energy.state.TemperatureSoilSnow(LoopInd1) - ConstFreezePoint) /
                            energy.state.PhaseChgFacSoilSnow(LoopInd1);
      energy.state.TemperatureSoilSnow(LoopInd1) = ConstFreezePoint;
    }
    if ((water.state.IndexPhaseChange(LoopInd1) == 1) && (EnergyRes(LoopInd1) < real(0.0))) {
      EnergyRes(LoopInd1) = real(0.0);
      water.state.IndexPhaseChange(LoopInd1) = 0;
    }
    if ((water.state.IndexPhaseChange(LoopInd1) == 2) && (EnergyRes(LoopInd1) > real(0.0))) {
      EnergyRes(LoopInd1) = real(0.0);
      water.state.IndexPhaseChange(LoopInd1) = 0;
    }
    GlacierPhaseChg(LoopInd1) = EnergyRes(LoopInd1) * domain.MainTimeStep / ConstLatHeatFusion;
  }

  if (nmlist.OptGlacierTreatment == 2) {
    if ((domain.NumSnowLayerNeg == 0) && (water.state.SnowWaterEquiv > real(0.0)) &&
        (energy.state.TemperatureSoilSnow(1) > ConstFreezePoint)) {
      EnergyRes(1) =
          (energy.state.TemperatureSoilSnow(1) - ConstFreezePoint) / energy.state.PhaseChgFacSoilSnow(1);
      energy.state.TemperatureSoilSnow(1) = ConstFreezePoint;
      GlacierPhaseChg(1) = EnergyRes(1) * domain.MainTimeStep / ConstLatHeatFusion;
      real SnowWaterPrev = water.state.SnowWaterEquiv;
      water.state.SnowWaterEquiv = std::max(real(0.0), SnowWaterPrev - GlacierPhaseChg(1));
      real SnowWaterRatio = water.state.SnowWaterEquiv / SnowWaterPrev;
      water.state.SnowDepth = std::max(real(0.0), SnowWaterRatio * water.state.SnowDepth);
      water.state.SnowDepth =
          std::min(std::max(water.state.SnowDepth, water.state.SnowWaterEquiv / real(500.0)),
                   water.state.SnowWaterEquiv / real(50.0));
      EnergyResLeft(1) = EnergyRes(1) -
                         ConstLatHeatFusion * (SnowWaterPrev - water.state.SnowWaterEquiv) / domain.MainTimeStep;
      if (EnergyResLeft(1) > real(0.0)) {
        GlacierPhaseChg(1) = EnergyResLeft(1) * domain.MainTimeStep / ConstLatHeatFusion;
        energy.state.TemperatureSoilSnow(1) =
            energy.state.TemperatureSoilSnow(1) + energy.state.PhaseChgFacSoilSnow(1) * EnergyResLeft(1);
      } else {
        GlacierPhaseChg(1) = real(0.0);
        EnergyRes(1) = real(0.0);
      }
      water.flux.MeltGroundSnow =
          std::max(real(0.0), SnowWaterPrev - water.state.SnowWaterEquiv) / domain.MainTimeStep;
      HeatLhTotPhsChg = ConstLatHeatFusion * water.flux.MeltGroundSnow;
      water.state.PondSfcThinSnwMelt = SnowWaterPrev - water.state.SnowWaterEquiv;
    }
  }

  for (int LoopInd1 = domain.NumSnowLayerNeg + 1; LoopInd1 <= 0; ++LoopInd1) {
    if ((water.state.IndexPhaseChange(LoopInd1) > 0) && (std::abs(EnergyRes(LoopInd1)) > real(0.0))) {
      EnergyResLeft(LoopInd1) = real(0.0);
      if (GlacierPhaseChg(LoopInd1) > real(0.0)) {
        MassWatIceTmp(LoopInd1) = std::max(real(0.0), MassWatIceInit(LoopInd1) - GlacierPhaseChg(LoopInd1));
        EnergyResLeft(LoopInd1) =
            EnergyRes(LoopInd1) -
            ConstLatHeatFusion * (MassWatIceInit(LoopInd1) - MassWatIceTmp(LoopInd1)) / domain.MainTimeStep;
      } else if (GlacierPhaseChg(LoopInd1) < real(0.0)) {
        MassWatIceTmp(LoopInd1) =
            std::min(MassWatTotInit(LoopInd1), MassWatIceInit(LoopInd1) - GlacierPhaseChg(LoopInd1));
        EnergyResLeft(LoopInd1) =
            EnergyRes(LoopInd1) -
            ConstLatHeatFusion * (MassWatIceInit(LoopInd1) - MassWatIceTmp(LoopInd1)) / domain.MainTimeStep;
      }
      MassWatLiqTmp(LoopInd1) = std::max(real(0.0), MassWatTotInit(LoopInd1) - MassWatIceTmp(LoopInd1));

      if (std::abs(EnergyResLeft(LoopInd1)) > real(0.0)) {
        energy.state.TemperatureSoilSnow(LoopInd1) =
            energy.state.TemperatureSoilSnow(LoopInd1) +
            energy.state.PhaseChgFacSoilSnow(LoopInd1) * EnergyResLeft(LoopInd1);
        if ((MassWatLiqTmp(LoopInd1) * MassWatIceTmp(LoopInd1)) > real(0.0)) {
          energy.state.TemperatureSoilSnow(LoopInd1) = ConstFreezePoint;
        }
      }
      HeatLhTotPhsChg =
          HeatLhTotPhsChg +
          ConstLatHeatFusion * (MassWatIceInit(LoopInd1) - MassWatIceTmp(LoopInd1)) / domain.MainTimeStep;

      water.flux.MeltGroundSnow =
          water.flux.MeltGroundSnow +
          std::max(real(0.0), (MassWatIceInit(LoopInd1) - MassWatIceTmp(LoopInd1))) / domain.MainTimeStep;
    }
  }

  if (nmlist.OptGlacierTreatment == 1) {
    for (int LoopInd1 = 1; LoopInd1 <= domain.NumSoilLayer; ++LoopInd1) {
      MassWatLiqTmp(LoopInd1) =
          water.state.SoilLiqWater(LoopInd1) * domain.ThicknessSnowSoilLayer(LoopInd1) * real(1000.0);
      MassWatIceTmp(LoopInd1) =
          (water.state.SoilMoisture(LoopInd1) - water.state.SoilLiqWater(LoopInd1)) *
          domain.ThicknessSnowSoilLayer(LoopInd1) * real(1000.0);
    }

    for (int LoopInd1 = 1; LoopInd1 <= domain.NumSoilLayer; ++LoopInd1) {
      water.state.IndexPhaseChange(LoopInd1) = 0;
      EnergyRes(LoopInd1) = real(0.0);
      GlacierPhaseChg(LoopInd1) = real(0.0);
      EnergyResLeft(LoopInd1) = real(0.0);
      MassWatIceInit(LoopInd1) = MassWatIceTmp(LoopInd1);
      MassWatLiqInit(LoopInd1) = MassWatLiqTmp(LoopInd1);
      MassWatTotInit(LoopInd1) = MassWatIceTmp(LoopInd1) + MassWatLiqTmp(LoopInd1);
    }

    for (int LoopInd1 = 1; LoopInd1 <= domain.NumSoilLayer; ++LoopInd1) {
      if ((MassWatIceTmp(LoopInd1) > real(0.0)) && (energy.state.TemperatureSoilSnow(LoopInd1) >= ConstFreezePoint)) {
        water.state.IndexPhaseChange(LoopInd1) = 1;
      }
      if ((MassWatLiqTmp(LoopInd1) > real(0.0)) && (energy.state.TemperatureSoilSnow(LoopInd1) < ConstFreezePoint)) {
        water.state.IndexPhaseChange(LoopInd1) = 2;
      }
      if ((domain.NumSnowLayerNeg == 0) && (water.state.SnowWaterEquiv > real(0.0)) && (LoopInd1 == 1)) {
        if (energy.state.TemperatureSoilSnow(LoopInd1) >= ConstFreezePoint) {
          water.state.IndexPhaseChange(LoopInd1) = 1;
        }
      }
    }

    for (int LoopInd1 = 1; LoopInd1 <= domain.NumSoilLayer; ++LoopInd1) {
      if (water.state.IndexPhaseChange(LoopInd1) > 0) {
        EnergyRes(LoopInd1) = (energy.state.TemperatureSoilSnow(LoopInd1) - ConstFreezePoint) /
                              energy.state.PhaseChgFacSoilSnow(LoopInd1);
        energy.state.TemperatureSoilSnow(LoopInd1) = ConstFreezePoint;
      }
      if ((water.state.IndexPhaseChange(LoopInd1) == 1) && (EnergyRes(LoopInd1) < real(0.0))) {
        EnergyRes(LoopInd1) = real(0.0);
        water.state.IndexPhaseChange(LoopInd1) = 0;
      }
      if ((water.state.IndexPhaseChange(LoopInd1) == 2) && (EnergyRes(LoopInd1) > real(0.0))) {
        EnergyRes(LoopInd1) = real(0.0);
        water.state.IndexPhaseChange(LoopInd1) = 0;
      }
      GlacierPhaseChg(LoopInd1) = EnergyRes(LoopInd1) * domain.MainTimeStep / ConstLatHeatFusion;
    }

    if ((domain.NumSnowLayerNeg == 0) && (water.state.SnowWaterEquiv > real(0.0)) &&
        (GlacierPhaseChg(1) > real(0.0))) {
      real SnowWaterPrev = water.state.SnowWaterEquiv;
      water.state.SnowWaterEquiv = std::max(real(0.0), SnowWaterPrev - GlacierPhaseChg(1));
      real SnowWaterRatio = water.state.SnowWaterEquiv / SnowWaterPrev;
      water.state.SnowDepth = std::max(real(0.0), SnowWaterRatio * water.state.SnowDepth);
      water.state.SnowDepth =
          std::min(std::max(water.state.SnowDepth, water.state.SnowWaterEquiv / real(500.0)),
                   water.state.SnowWaterEquiv / real(50.0));
      EnergyResLeft(1) = EnergyRes(1) -
                         ConstLatHeatFusion * (SnowWaterPrev - water.state.SnowWaterEquiv) / domain.MainTimeStep;
      if (EnergyResLeft(1) > real(0.0)) {
        GlacierPhaseChg(1) = EnergyResLeft(1) * domain.MainTimeStep / ConstLatHeatFusion;
        EnergyRes(1) = EnergyResLeft(1);
        water.state.IndexPhaseChange(1) = 1;
      } else {
        GlacierPhaseChg(1) = real(0.0);
        EnergyRes(1) = real(0.0);
        water.state.IndexPhaseChange(1) = 0;
      }
      water.flux.MeltGroundSnow =
          std::max(real(0.0), (SnowWaterPrev - water.state.SnowWaterEquiv)) / domain.MainTimeStep;
      HeatLhTotPhsChg = ConstLatHeatFusion * water.flux.MeltGroundSnow;
      water.state.PondSfcThinSnwMelt = SnowWaterPrev - water.state.SnowWaterEquiv;
    }

    for (int LoopInd1 = 1; LoopInd1 <= domain.NumSoilLayer; ++LoopInd1) {
      if ((water.state.IndexPhaseChange(LoopInd1) > 0) && (std::abs(EnergyRes(LoopInd1)) > real(0.0))) {
        EnergyResLeft(LoopInd1) = real(0.0);
        if (GlacierPhaseChg(LoopInd1) > real(0.0)) {
          MassWatIceTmp(LoopInd1) = std::max(real(0.0), MassWatIceInit(LoopInd1) - GlacierPhaseChg(LoopInd1));
          EnergyResLeft(LoopInd1) =
              EnergyRes(LoopInd1) -
              ConstLatHeatFusion * (MassWatIceInit(LoopInd1) - MassWatIceTmp(LoopInd1)) / domain.MainTimeStep;
        } else if (GlacierPhaseChg(LoopInd1) < real(0.0)) {
          MassWatIceTmp(LoopInd1) =
              std::min(MassWatTotInit(LoopInd1), MassWatIceInit(LoopInd1) - GlacierPhaseChg(LoopInd1));
          EnergyResLeft(LoopInd1) =
              EnergyRes(LoopInd1) -
              ConstLatHeatFusion * (MassWatIceInit(LoopInd1) - MassWatIceTmp(LoopInd1)) / domain.MainTimeStep;
        }
        MassWatLiqTmp(LoopInd1) = std::max(real(0.0), MassWatTotInit(LoopInd1) - MassWatIceTmp(LoopInd1));

        if (std::abs(EnergyResLeft(LoopInd1)) > real(0.0)) {
          energy.state.TemperatureSoilSnow(LoopInd1) =
              energy.state.TemperatureSoilSnow(LoopInd1) +
              energy.state.PhaseChgFacSoilSnow(LoopInd1) * EnergyResLeft(LoopInd1);
        }
        HeatLhTotPhsChg =
            HeatLhTotPhsChg +
            ConstLatHeatFusion * (MassWatIceInit(LoopInd1) - MassWatIceTmp(LoopInd1)) / domain.MainTimeStep;
      }
    }
    EnergyResLeft.fill(real(0.0));
    GlacierPhaseChg.fill(real(0.0));

    auto any_temp_gt_frz = [&]() {
      for (int k = 1; k <= domain.NumSoilLayer; ++k) {
        if (energy.state.TemperatureSoilSnow(k) > ConstFreezePoint) return true;
      }
      return false;
    };
    auto any_temp_lt_frz = [&]() {
      for (int k = 1; k <= domain.NumSoilLayer; ++k) {
        if (energy.state.TemperatureSoilSnow(k) < ConstFreezePoint) return true;
      }
      return false;
    };
    auto any_ice_gt_zero = [&]() {
      for (int k = 1; k <= domain.NumSoilLayer; ++k) {
        if (MassWatIceTmp(k) > real(0.0)) return true;
      }
      return false;
    };
    auto any_liq_gt_zero = [&]() {
      for (int k = 1; k <= domain.NumSoilLayer; ++k) {
        if (MassWatLiqTmp(k) > real(0.0)) return true;
      }
      return false;
    };

    if (any_temp_gt_frz() && any_temp_lt_frz()) {
      for (int LoopInd1 = 1; LoopInd1 <= domain.NumSoilLayer; ++LoopInd1) {
        if (energy.state.TemperatureSoilSnow(LoopInd1) > ConstFreezePoint) {
          EnergyResLeft(LoopInd1) = (energy.state.TemperatureSoilSnow(LoopInd1) - ConstFreezePoint) /
                                    energy.state.PhaseChgFacSoilSnow(LoopInd1);
          for (int LoopInd2 = 1; LoopInd2 <= domain.NumSoilLayer; ++LoopInd2) {
            if ((LoopInd1 != LoopInd2) && (energy.state.TemperatureSoilSnow(LoopInd2) < ConstFreezePoint) &&
                (EnergyResLeft(LoopInd1) > real(0.1))) {
              EnergyResLeft(LoopInd2) = (energy.state.TemperatureSoilSnow(LoopInd2) - ConstFreezePoint) /
                                        energy.state.PhaseChgFacSoilSnow(LoopInd2);
              if (std::abs(EnergyResLeft(LoopInd2)) > EnergyResLeft(LoopInd1)) {
                EnergyResLeft(LoopInd2) = EnergyResLeft(LoopInd2) + EnergyResLeft(LoopInd1);
                energy.state.TemperatureSoilSnow(LoopInd2) =
                    ConstFreezePoint + EnergyResLeft(LoopInd2) * energy.state.PhaseChgFacSoilSnow(LoopInd2);
                EnergyResLeft(LoopInd1) = real(0.0);
              } else {
                EnergyResLeft(LoopInd1) = EnergyResLeft(LoopInd1) + EnergyResLeft(LoopInd2);
                EnergyResLeft(LoopInd2) = real(0.0);
                energy.state.TemperatureSoilSnow(LoopInd2) = ConstFreezePoint;
              }
            }
          }
          energy.state.TemperatureSoilSnow(LoopInd1) =
              ConstFreezePoint + EnergyResLeft(LoopInd1) * energy.state.PhaseChgFacSoilSnow(LoopInd1);
        }
      }
    }

    if (any_temp_gt_frz() && any_temp_lt_frz()) {
      for (int LoopInd1 = 1; LoopInd1 <= domain.NumSoilLayer; ++LoopInd1) {
        if (energy.state.TemperatureSoilSnow(LoopInd1) < ConstFreezePoint) {
          EnergyResLeft(LoopInd1) = (energy.state.TemperatureSoilSnow(LoopInd1) - ConstFreezePoint) /
                                    energy.state.PhaseChgFacSoilSnow(LoopInd1);
          for (int LoopInd2 = 1; LoopInd2 <= domain.NumSoilLayer; ++LoopInd2) {
            if ((LoopInd1 != LoopInd2) && (energy.state.TemperatureSoilSnow(LoopInd2) > ConstFreezePoint) &&
                (EnergyResLeft(LoopInd1) < real(-0.1))) {
              EnergyResLeft(LoopInd2) = (energy.state.TemperatureSoilSnow(LoopInd2) - ConstFreezePoint) /
                                        energy.state.PhaseChgFacSoilSnow(LoopInd2);
              if (EnergyResLeft(LoopInd2) > std::abs(EnergyResLeft(LoopInd1))) {
                EnergyResLeft(LoopInd2) = EnergyResLeft(LoopInd2) + EnergyResLeft(LoopInd1);
                energy.state.TemperatureSoilSnow(LoopInd2) =
                    ConstFreezePoint + EnergyResLeft(LoopInd2) * energy.state.PhaseChgFacSoilSnow(LoopInd2);
                EnergyResLeft(LoopInd1) = real(0.0);
              } else {
                EnergyResLeft(LoopInd1) = EnergyResLeft(LoopInd1) + EnergyResLeft(LoopInd2);
                EnergyResLeft(LoopInd2) = real(0.0);
                energy.state.TemperatureSoilSnow(LoopInd2) = ConstFreezePoint;
              }
            }
          }
          energy.state.TemperatureSoilSnow(LoopInd1) =
              ConstFreezePoint + EnergyResLeft(LoopInd1) * energy.state.PhaseChgFacSoilSnow(LoopInd1);
        }
      }
    }

    if (any_temp_gt_frz() && any_ice_gt_zero()) {
      for (int LoopInd1 = 1; LoopInd1 <= domain.NumSoilLayer; ++LoopInd1) {
        if (energy.state.TemperatureSoilSnow(LoopInd1) > ConstFreezePoint) {
          EnergyResLeft(LoopInd1) = (energy.state.TemperatureSoilSnow(LoopInd1) - ConstFreezePoint) /
                                    energy.state.PhaseChgFacSoilSnow(LoopInd1);
          GlacierPhaseChg(LoopInd1) = EnergyResLeft(LoopInd1) * domain.MainTimeStep / ConstLatHeatFusion;
          for (int LoopInd2 = 1; LoopInd2 <= domain.NumSoilLayer; ++LoopInd2) {
            if ((LoopInd1 != LoopInd2) && (MassWatIceTmp(LoopInd2) > real(0.0)) &&
                (GlacierPhaseChg(LoopInd1) > real(0.1))) {
              if (MassWatIceTmp(LoopInd2) > GlacierPhaseChg(LoopInd1)) {
                MassWatIceTmp(LoopInd2) = MassWatIceTmp(LoopInd2) - GlacierPhaseChg(LoopInd1);
                HeatLhTotPhsChg =
                    HeatLhTotPhsChg + ConstLatHeatFusion * GlacierPhaseChg(LoopInd1) / domain.MainTimeStep;
                energy.state.TemperatureSoilSnow(LoopInd2) = ConstFreezePoint;
                GlacierPhaseChg(LoopInd1) = real(0.0);
              } else {
                GlacierPhaseChg(LoopInd1) = GlacierPhaseChg(LoopInd1) - MassWatIceTmp(LoopInd2);
                HeatLhTotPhsChg =
                    HeatLhTotPhsChg + ConstLatHeatFusion * MassWatIceTmp(LoopInd2) / domain.MainTimeStep;
                MassWatIceTmp(LoopInd2) = real(0.0);
                energy.state.TemperatureSoilSnow(LoopInd2) = ConstFreezePoint;
              }
              MassWatLiqTmp(LoopInd2) = std::max(real(0.0), MassWatTotInit(LoopInd2) - MassWatIceTmp(LoopInd2));
            }
          }
          EnergyResLeft(LoopInd1) = GlacierPhaseChg(LoopInd1) * ConstLatHeatFusion / domain.MainTimeStep;
          energy.state.TemperatureSoilSnow(LoopInd1) =
              ConstFreezePoint + EnergyResLeft(LoopInd1) * energy.state.PhaseChgFacSoilSnow(LoopInd1);
        }
      }
    }

    if (any_temp_lt_frz() && any_liq_gt_zero()) {
      for (int LoopInd1 = 1; LoopInd1 <= domain.NumSoilLayer; ++LoopInd1) {
        if (energy.state.TemperatureSoilSnow(LoopInd1) < ConstFreezePoint) {
          EnergyResLeft(LoopInd1) = (energy.state.TemperatureSoilSnow(LoopInd1) - ConstFreezePoint) /
                                    energy.state.PhaseChgFacSoilSnow(LoopInd1);
          GlacierPhaseChg(LoopInd1) = EnergyResLeft(LoopInd1) * domain.MainTimeStep / ConstLatHeatFusion;
          for (int LoopInd2 = 1; LoopInd2 <= domain.NumSoilLayer; ++LoopInd2) {
            if ((LoopInd1 != LoopInd2) && (MassWatLiqTmp(LoopInd2) > real(0.0)) &&
                (GlacierPhaseChg(LoopInd1) < real(-0.1))) {
              if (MassWatLiqTmp(LoopInd2) > std::abs(GlacierPhaseChg(LoopInd1))) {
                MassWatIceTmp(LoopInd2) = MassWatIceTmp(LoopInd2) - GlacierPhaseChg(LoopInd1);
                HeatLhTotPhsChg =
                    HeatLhTotPhsChg + ConstLatHeatFusion * GlacierPhaseChg(LoopInd1) / domain.MainTimeStep;
                energy.state.TemperatureSoilSnow(LoopInd2) = ConstFreezePoint;
                GlacierPhaseChg(LoopInd1) = real(0.0);
              } else {
                GlacierPhaseChg(LoopInd1) = GlacierPhaseChg(LoopInd1) + MassWatLiqTmp(LoopInd2);
                HeatLhTotPhsChg =
                    HeatLhTotPhsChg - ConstLatHeatFusion * MassWatLiqTmp(LoopInd2) / domain.MainTimeStep;
                MassWatIceTmp(LoopInd2) = MassWatTotInit(LoopInd2);
                energy.state.TemperatureSoilSnow(LoopInd2) = ConstFreezePoint;
              }
              MassWatLiqTmp(LoopInd2) = std::max(real(0.0), MassWatTotInit(LoopInd2) - MassWatIceTmp(LoopInd2));
            }
          }
          EnergyResLeft(LoopInd1) = GlacierPhaseChg(LoopInd1) * ConstLatHeatFusion / domain.MainTimeStep;
          energy.state.TemperatureSoilSnow(LoopInd1) =
              ConstFreezePoint + EnergyResLeft(LoopInd1) * energy.state.PhaseChgFacSoilSnow(LoopInd1);
        }
      }
    }
  }
  (void)HeatLhTotPhsChg;

  for (int LoopInd1 = domain.NumSnowLayerNeg + 1; LoopInd1 <= 0; ++LoopInd1) {
    water.state.SnowLiqWater(LoopInd1) = MassWatLiqTmp(LoopInd1);
    water.state.SnowIce(LoopInd1) = MassWatIceTmp(LoopInd1);
  }
  for (int LoopInd1 = 1; LoopInd1 <= domain.NumSoilLayer; ++LoopInd1) {
    if (nmlist.OptGlacierTreatment == 1) {
      water.state.SoilLiqWater(LoopInd1) =
          MassWatLiqTmp(LoopInd1) / (real(1000.0) * domain.ThicknessSnowSoilLayer(LoopInd1));
      water.state.SoilLiqWater(LoopInd1) =
          std::max(real(0.0), std::min(real(1.0), water.state.SoilLiqWater(LoopInd1)));
    } else if (nmlist.OptGlacierTreatment == 2) {
      water.state.SoilLiqWater(LoopInd1) = real(0.0);
    }
    water.state.SoilMoisture(LoopInd1) = real(1.0);
  }
}

// ============================================================================
// EnergyMainGlacierMod.F90
// ============================================================================
CPPWRF_HD inline void EnergyMainGlacier(NoahmpType& noahmp) {
  auto& nmlist = noahmp.config.nmlist;
  auto& forcing = noahmp.forcing;
  auto& energy = noahmp.energy;
  auto& water = noahmp.water;

  SnowCoverGlacier(noahmp);
  GroundRoughnessPropertyGlacier(noahmp);
  GroundThermalPropertyGlacier(noahmp);
  SurfaceAlbedoGlacier(noahmp);
  SurfaceRadiationGlacier(noahmp);
  SurfaceEmissivityGlacier(noahmp);
  ResistanceGroundEvaporationGlacier(noahmp);
  PsychrometricVariableGlacier(noahmp);

  energy.state.TemperatureGrdBare = energy.state.TemperatureGrd;
  energy.state.ExchCoeffMomBare = energy.state.ExchCoeffMomSfc;
  energy.state.ExchCoeffShBare = energy.state.ExchCoeffShSfc;
  SurfaceEnergyFluxGlacier(noahmp);

  energy.state.WindStressEwSfc = energy.state.WindStressEwBare;
  energy.state.WindStressNsSfc = energy.state.WindStressNsBare;
  energy.flux.RadLwNetSfc = energy.flux.RadLwNetBareGrd;
  energy.flux.HeatSensibleSfc = energy.flux.HeatSensibleBareGrd;
  energy.flux.HeatLatentGrd = energy.flux.HeatLatentBareGrd;
  energy.flux.HeatGroundTot = energy.flux.HeatGroundBareGrd;
  energy.state.TemperatureGrd = energy.state.TemperatureGrdBare;
  energy.state.TemperatureAir2m = energy.state.TemperatureAir2mBare;
  energy.flux.HeatPrecipAdvSfc = energy.flux.HeatPrecipAdvBareGrd;
  energy.state.TemperatureSfc = energy.state.TemperatureGrd;
  energy.state.ExchCoeffMomSfc = energy.state.ExchCoeffMomBare;
  energy.state.ExchCoeffShSfc = energy.state.ExchCoeffShBare;
  energy.state.SpecHumiditySfcMean = energy.state.SpecHumiditySfc;
  energy.state.SpecHumidity2m = energy.state.SpecHumidity2mBare;
  energy.state.RoughLenMomSfcToAtm = energy.state.RoughLenMomGrd;

  energy.flux.RadLwEmitSfc = forcing.RadLwDownRefHeight + energy.flux.RadLwNetSfc;
  if (energy.flux.RadLwEmitSfc <= real(0.0)) raise_fatal(kFatalNoahmpLongwave);

  energy.state.TemperatureRadSfc = rpow(
      (energy.flux.RadLwEmitSfc - (real(1.0) - energy.state.EmissivitySfc) * forcing.RadLwDownRefHeight) /
          (energy.state.EmissivitySfc * ConstStefanBoltzmann),
      real(0.25));

  GlacierTemperatureMain(noahmp);

  if (nmlist.OptSnowSoilTempTime == 2) {
    if ((water.state.SnowDepth > real(0.05)) && (energy.state.TemperatureGrd > ConstFreezePoint)) {
      energy.state.TemperatureGrdBare = ConstFreezePoint;
      energy.state.TemperatureGrd = energy.state.TemperatureGrdBare;
      energy.state.TemperatureSfc = energy.state.TemperatureGrdBare;
    }
  }

  GlacierPhaseChange(noahmp);

  if (forcing.RadSwDownRefHeight > real(0.0)) {
    energy.state.AlbedoSfc = energy.flux.RadSwReflSfc / forcing.RadSwDownRefHeight;
  } else {
    energy.state.AlbedoSfc = undefined_real;
  }
}

// ============================================================================
// SnowpackHydrologyGlacierMod.F90
// ============================================================================
CPPWRF_HD inline void SnowpackHydrologyGlacier(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& nmlist = noahmp.config.nmlist;
  auto& energy = noahmp.energy;
  auto& water = noahmp.water;

  water.state.SnowEffPorosity.fill(real(0.0));
  water.flux.SnowBotOutflow = real(0.0);
  real InflowSnowLayer = real(0.0);
  water.flux.OutflowSnowLayer.fill(real(0.0));

  if (water.state.SnowWaterEquiv == real(0.0)) {
    if (nmlist.OptGlacierTreatment == 1) {
      water.state.SoilIce(1) =
          water.state.SoilIce(1) +
          (water.flux.FrostSnowSfcIce - water.flux.SublimSnowSfcIce) * domain.MainTimeStep /
              (domain.ThicknessSnowSoilLayer(1) * real(1000.0));
    } else if (nmlist.OptGlacierTreatment == 2) {
      energy.flux.HeatSensibleSfc =
          energy.flux.HeatSensibleSfc -
          (water.flux.FrostSnowSfcIce - water.flux.SublimSnowSfcIce) * ConstLatHeatSublim;
      water.flux.FrostSnowSfcIce = real(0.0);
      water.flux.SublimSnowSfcIce = real(0.0);
    }
  }

  if ((domain.NumSnowLayerNeg == 0) && (water.state.SnowWaterEquiv > real(0.0))) {
    if (nmlist.OptGlacierTreatment == 1) {
      real SnowWaterTmp = water.state.SnowWaterEquiv;
      water.state.SnowWaterEquiv =
          water.state.SnowWaterEquiv - water.flux.SublimSnowSfcIce * domain.MainTimeStep +
          water.flux.FrostSnowSfcIce * domain.MainTimeStep;
      real SnowWaterRatio = water.state.SnowWaterEquiv / SnowWaterTmp;
      water.state.SnowDepth = std::max(real(0.0), SnowWaterRatio * water.state.SnowDepth);
      water.state.SnowDepth =
          std::min(std::max(water.state.SnowDepth, water.state.SnowWaterEquiv / real(500.0)),
                   water.state.SnowWaterEquiv / real(50.0));
    } else if (nmlist.OptGlacierTreatment == 2) {
      energy.flux.HeatSensibleSfc =
          energy.flux.HeatSensibleSfc -
          (water.flux.FrostSnowSfcIce - water.flux.SublimSnowSfcIce) * ConstLatHeatSublim;
      water.flux.FrostSnowSfcIce = real(0.0);
      water.flux.SublimSnowSfcIce = real(0.0);
    }
    if (water.state.SnowWaterEquiv < real(0.0)) {
      water.state.SoilIce(1) =
          water.state.SoilIce(1) + water.state.SnowWaterEquiv / (domain.ThicknessSnowSoilLayer(1) * real(1000.0));
      water.state.SnowWaterEquiv = real(0.0);
      water.state.SnowDepth = real(0.0);
    }
    if (water.state.SoilIce(1) < real(0.0)) {
      water.state.SoilLiqWater(1) = water.state.SoilLiqWater(1) + water.state.SoilIce(1);
      water.state.SoilIce(1) = real(0.0);
    }
  }

  if ((water.state.SnowDepth <= real(1.0e-8)) || (water.state.SnowWaterEquiv <= real(1.0e-6))) {
    water.state.SnowDepth = real(0.0);
    water.state.SnowWaterEquiv = real(0.0);
  }

  if (domain.NumSnowLayerNeg < 0) {
    real SnowWaterTmp =
        water.state.SnowIce(domain.NumSnowLayerNeg + 1) + water.state.SnowLiqWater(domain.NumSnowLayerNeg + 1);
    real SnowIceTmp = water.state.SnowIce(domain.NumSnowLayerNeg + 1) -
                      water.flux.SublimSnowSfcIce * domain.MainTimeStep +
                      water.flux.FrostSnowSfcIce * domain.MainTimeStep;
    water.state.SnowIce(domain.NumSnowLayerNeg + 1) = SnowIceTmp;
    if ((SnowIceTmp < real(1.0e-6)) && (domain.NumSnowLayerNeg < 0)) {
      SnowLayerCombine(noahmp);
    }
    if ((SnowIceTmp >= real(1.0e-6)) && (domain.NumSnowLayerNeg < 0)) {
      domain.ThicknessSnowSoilLayer(domain.NumSnowLayerNeg + 1) =
          domain.ThicknessSnowSoilLayer(domain.NumSnowLayerNeg + 1) *
          (water.state.SnowIce(domain.NumSnowLayerNeg + 1) + water.state.SnowLiqWater(domain.NumSnowLayerNeg + 1)) /
          SnowWaterTmp;
    }
    if (domain.NumSnowLayerNeg < 0) {
      water.state.SnowLiqWater(domain.NumSnowLayerNeg + 1) =
          water.state.SnowLiqWater(domain.NumSnowLayerNeg + 1) +
          water.flux.RainfallGround * domain.MainTimeStep;
      water.state.SnowLiqWater(domain.NumSnowLayerNeg + 1) =
          std::max(real(0.0), water.state.SnowLiqWater(domain.NumSnowLayerNeg + 1));
    }
  }

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
    water.state.SnowIceVol(LoopInd) =
        std::min(real(1.0),
                 water.state.SnowIce(LoopInd) / (domain.ThicknessSnowSoilLayer(LoopInd) * ConstDensityIce));
    water.state.SnowEffPorosity(LoopInd) = real(1.0) - water.state.SnowIceVol(LoopInd);
  }

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
    water.state.SnowLiqWater(LoopInd) = water.state.SnowLiqWater(LoopInd) + InflowSnowLayer;
    water.state.SnowLiqWaterVol(LoopInd) =
        water.state.SnowLiqWater(LoopInd) / (domain.ThicknessSnowSoilLayer(LoopInd) * ConstDensityWater);
    water.flux.OutflowSnowLayer(LoopInd) = std::max(
        real(0.0),
        (water.state.SnowLiqWaterVol(LoopInd) - water.param.SnowLiqHoldCap * water.state.SnowEffPorosity(LoopInd)) *
            domain.ThicknessSnowSoilLayer(LoopInd));
    if (LoopInd == 0) {
      water.flux.OutflowSnowLayer(LoopInd) =
          std::max((water.state.SnowLiqWaterVol(LoopInd) - water.state.SnowEffPorosity(LoopInd)) *
                       domain.ThicknessSnowSoilLayer(LoopInd),
                   water.param.SnowLiqReleaseFac * domain.MainTimeStep * water.flux.OutflowSnowLayer(LoopInd));
    }
    water.flux.OutflowSnowLayer(LoopInd) = water.flux.OutflowSnowLayer(LoopInd) * ConstDensityWater;
    water.state.SnowLiqWater(LoopInd) = water.state.SnowLiqWater(LoopInd) - water.flux.OutflowSnowLayer(LoopInd);
    if ((water.state.SnowLiqWater(LoopInd) / (water.state.SnowIce(LoopInd) + water.state.SnowLiqWater(LoopInd))) >
        water.param.SnowLiqFracMax) {
      water.flux.OutflowSnowLayer(LoopInd) =
          water.flux.OutflowSnowLayer(LoopInd) +
          (water.state.SnowLiqWater(LoopInd) -
           water.param.SnowLiqFracMax / (real(1.0) - water.param.SnowLiqFracMax) * water.state.SnowIce(LoopInd));
      water.state.SnowLiqWater(LoopInd) =
          water.param.SnowLiqFracMax / (real(1.0) - water.param.SnowLiqFracMax) * water.state.SnowIce(LoopInd);
    }
    InflowSnowLayer = water.flux.OutflowSnowLayer(LoopInd);
    water.state.SnowLiqWaterVol(LoopInd) =
        water.state.SnowLiqWater(LoopInd) / (domain.ThicknessSnowSoilLayer(LoopInd) * ConstDensityWater);
  }

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
    domain.ThicknessSnowSoilLayer(LoopInd) =
        std::max(domain.ThicknessSnowSoilLayer(LoopInd),
                 water.state.SnowLiqWater(LoopInd) / ConstDensityWater +
                     water.state.SnowIce(LoopInd) / ConstDensityIce);
  }

  water.flux.SnowBotOutflow = water.flux.OutflowSnowLayer(0) / domain.MainTimeStep;
  for (int k = -domain.NumSnowLayerMax + 1; k <= 0; ++k) {
    water.flux.OutflowSnowLayer(k) = water.flux.OutflowSnowLayer(k) / domain.MainTimeStep;
  }
}

// ============================================================================
// SnowWaterMainGlacierMod.F90
// ============================================================================
CPPWRF_HD inline void SnowWaterMainGlacier(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& nmlist = noahmp.config.nmlist;
  auto& energy = noahmp.energy;
  auto& water = noahmp.water;

  water.flux.GlacierExcessFlow = real(0.0);
  water.state.PondSfcThinSnwComb = real(0.0);
  water.state.PondSfcThinSnwTrans = real(0.0);

  SnowfallAfterCanopyIntercept(noahmp);

  if (domain.NumSnowLayerNeg < 0 && nmlist.OptSnowCompaction == 1) SnowpackCompaction(noahmp);
  if (domain.NumSnowLayerNeg < 0 && nmlist.OptSnowCompaction == 2) SnowpackCompactionAR24(noahmp);

  if (domain.NumSnowLayerNeg < 0) SnowLayerCombine(noahmp);
  if (domain.NumSnowLayerNeg < 0) SnowLayerDivide(noahmp);

  SnowpackHydrologyGlacier(noahmp);

  for (int LoopInd = -domain.NumSnowLayerMax + 1; LoopInd <= domain.NumSnowLayerNeg; ++LoopInd) {
    water.state.SnowIce(LoopInd) = real(0.0);
    water.state.SnowLiqWater(LoopInd) = real(0.0);
    energy.state.TemperatureSoilSnow(LoopInd) = real(0.0);
    domain.ThicknessSnowSoilLayer(LoopInd) = real(0.0);
    domain.DepthSnowSoilLayer(LoopInd) = real(0.0);
  }

  if (water.state.SnowWaterEquiv > water.param.SnoWatEqvMaxGlacier) {
    real SnowDensBulk = water.state.SnowIce(0) / domain.ThicknessSnowSoilLayer(0);
    water.flux.GlacierExcessFlow = water.state.SnowWaterEquiv - water.param.SnoWatEqvMaxGlacier;
    water.state.SnowIce(0) = water.state.SnowIce(0) - water.flux.GlacierExcessFlow;
    domain.ThicknessSnowSoilLayer(0) =
        domain.ThicknessSnowSoilLayer(0) - water.flux.GlacierExcessFlow / SnowDensBulk;
    water.flux.GlacierExcessFlow = water.flux.GlacierExcessFlow / domain.MainTimeStep;
  }

  if (domain.NumSnowLayerNeg < 0) {
    water.state.SnowWaterEquiv = real(0.0);
    for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
      water.state.SnowWaterEquiv =
          water.state.SnowWaterEquiv + water.state.SnowIce(LoopInd) + water.state.SnowLiqWater(LoopInd);
    }
  }

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
    domain.ThicknessSnowSoilLayer(LoopInd) = -domain.ThicknessSnowSoilLayer(LoopInd);
  }

  domain.ThicknessSnowSoilLayer(1) = domain.DepthSoilLayer(1);
  for (int LoopInd = 2; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    domain.ThicknessSnowSoilLayer(LoopInd) =
        domain.DepthSoilLayer(LoopInd) - domain.DepthSoilLayer(LoopInd - 1);
  }

  domain.DepthSnowSoilLayer(domain.NumSnowLayerNeg + 1) =
      domain.ThicknessSnowSoilLayer(domain.NumSnowLayerNeg + 1);
  for (int LoopInd = domain.NumSnowLayerNeg + 2; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    domain.DepthSnowSoilLayer(LoopInd) =
        domain.DepthSnowSoilLayer(LoopInd - 1) + domain.ThicknessSnowSoilLayer(LoopInd);
  }

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    domain.ThicknessSnowSoilLayer(LoopInd) = -domain.ThicknessSnowSoilLayer(LoopInd);
  }

  if (domain.NumSnowLayerNeg < 0) {
    water.state.SnowDepth = real(0.0);
    for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
      water.state.SnowDepth = water.state.SnowDepth + domain.ThicknessSnowSoilLayer(LoopInd);
    }
  }

  if ((water.state.SnowDepth <= real(1.0e-6)) || (water.state.SnowWaterEquiv <= real(1.0e-6))) {
    water.state.SnowDepth = real(0.0);
    water.state.SnowWaterEquiv = real(0.0);
  }
}

// ============================================================================
// WaterMainGlacierMod.F90
// ============================================================================
CPPWRF_HD inline void WaterMainGlacier(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& nmlist = noahmp.config.nmlist;
  auto& energy = noahmp.energy;
  auto& water = noahmp.water;

  FArray1D<real, 1, kMaxSoilLayers> SoilIceTmp;
  FArray1D<real, 1, kMaxSoilLayers> SoilLiqWaterTmp;
  SoilIceTmp.fill(real(0.0));
  SoilLiqWaterTmp.fill(real(0.0));
  water.flux.GlacierExcessFlow = real(0.0);
  water.flux.RunoffSubsurface = real(0.0);
  water.flux.RunoffSurface = real(0.0);
  water.flux.SnowDepthIncr = real(0.0);

  for (int k = 1; k <= domain.NumSoilLayer; ++k) {
    water.state.SoilIce(k) = std::max(real(0.0), water.state.SoilMoisture(k) - water.state.SoilLiqWater(k));
    SoilIceTmp(k) = water.state.SoilIce(k);
    SoilLiqWaterTmp(k) = water.state.SoilLiqWater(k);
  }
  water.state.SnowWaterEquivPrev = water.state.SnowWaterEquiv;

  water.flux.VaporizeGrd = std::max(energy.flux.HeatLatentGrd / energy.state.LatHeatVapGrd, real(0.0));
  water.flux.CondenseVapGrd =
      std::abs(std::min(energy.flux.HeatLatentGrd / energy.state.LatHeatVapGrd, real(0.0)));
  water.flux.EvapGroundNet = water.flux.VaporizeGrd - water.flux.CondenseVapGrd;

  water.flux.SnowDepthIncr = water.flux.SnowfallGround / water.state.SnowfallDensity;
  water.flux.SublimSnowSfcIce = water.flux.VaporizeGrd;
  water.flux.FrostSnowSfcIce = water.flux.CondenseVapGrd;

  SnowWaterMainGlacier(noahmp);

  water.flux.SoilSfcInflow =
      (water.state.PondSfcThinSnwMelt + water.state.PondSfcThinSnwComb + water.state.PondSfcThinSnwTrans) /
      domain.MainTimeStep * real(0.001);
  if (domain.NumSnowLayerNeg == 0) {
    water.flux.SoilSfcInflow =
        water.flux.SoilSfcInflow + (water.flux.SnowBotOutflow + water.flux.RainfallGround) * real(0.001);
  } else {
    water.flux.SoilSfcInflow = water.flux.SoilSfcInflow + water.flux.SnowBotOutflow * real(0.001);
  }

  water.flux.RunoffSurface = water.flux.SoilSfcInflow * real(1000.0);

  real WatReplaceSublim = real(0.0);
  if (nmlist.OptGlacierTreatment == 1) {
    WatReplaceSublim = real(0.0);
    for (int LoopInd = 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
      WatReplaceSublim =
          WatReplaceSublim +
          domain.ThicknessSnowSoilLayer(LoopInd) *
              (water.state.SoilIce(LoopInd) - SoilIceTmp(LoopInd) +
               water.state.SoilLiqWater(LoopInd) - SoilLiqWaterTmp(LoopInd));
    }
    WatReplaceSublim = WatReplaceSublim * real(1000.0) / domain.MainTimeStep;
    for (int LoopInd = 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
      water.state.SoilIce(LoopInd) = std::min(real(1.0), SoilIceTmp(LoopInd));
    }
  } else if (nmlist.OptGlacierTreatment == 2) {
    for (int LoopInd = 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
      water.state.SoilIce(LoopInd) = real(1.0);
    }
  }
  for (int LoopInd = 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    water.state.SoilLiqWater(LoopInd) = real(1.0) - water.state.SoilIce(LoopInd);
  }

  if (nmlist.OptGlacierTreatment == 1) {
    water.flux.RunoffSubsurface = water.flux.GlacierExcessFlow + WatReplaceSublim;
  } else if (nmlist.OptGlacierTreatment == 2) {
    water.flux.RunoffSubsurface = water.flux.GlacierExcessFlow;
    water.flux.VaporizeGrd = water.flux.SublimSnowSfcIce;
    water.flux.CondenseVapGrd = water.flux.FrostSnowSfcIce;
  }

  if (nmlist.OptGlacierTreatment == 2) {
    water.flux.EvapGroundNet = water.flux.VaporizeGrd - water.flux.CondenseVapGrd;
    energy.flux.HeatLatentGrd = water.flux.EvapGroundNet * energy.state.LatHeatVapGrd;
  }
}

// ============================================================================
// BalanceErrorCheckGlacierMod.F90
// ============================================================================
CPPWRF_HD inline void BalanceWaterInitGlacier(NoahmpType& noahmp) {
  auto& water = noahmp.water;
  water.state.WaterStorageTotBeg = water.state.SnowWaterEquiv;
}

CPPWRF_HD inline void BalanceWaterCheckGlacier(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& water = noahmp.water;

  water.state.WaterStorageTotEnd = water.state.SnowWaterEquiv;
  water.state.WaterBalanceError =
      water.state.WaterStorageTotEnd - water.state.WaterStorageTotBeg -
      (water.flux.PrecipTotRefHeight - water.flux.EvapGroundNet - water.flux.RunoffSurface -
       water.flux.RunoffSubsurface) *
          domain.MainTimeStep;
  if (std::abs(water.state.WaterBalanceError) > real(0.1)) raise_fatal(kFatalNoahmpWater);
}

CPPWRF_HD inline void BalanceEnergyCheckGlacier(NoahmpType& noahmp) {
  auto& forcing = noahmp.forcing;
  auto& energy = noahmp.energy;

  energy.state.RadSwBalanceError =
      forcing.RadSwDownRefHeight - (energy.flux.RadSwAbsSfc + energy.flux.RadSwReflSfc);
  energy.state.EnergyBalanceError =
      energy.flux.RadSwAbsGrd + energy.flux.HeatPrecipAdvSfc -
      (energy.flux.RadLwNetSfc + energy.flux.HeatSensibleSfc + energy.flux.HeatLatentGrd +
       energy.flux.HeatGroundTot);
  if (std::abs(energy.state.RadSwBalanceError) > real(0.01)) raise_fatal(kFatalNoahmpSolar);
  if (std::abs(energy.state.EnergyBalanceError) > real(0.01)) raise_fatal(kFatalNoahmpGlacierEnergy);
}

// ============================================================================
// NoahmpMainGlacierMod.F90
// ============================================================================
CPPWRF_HD inline void NoahmpMainGlacier(NoahmpType& noahmp) {
  ProcessAtmosForcing(noahmp);
  GeneralInitGlacier(noahmp);
  BalanceWaterInitGlacier(noahmp);
  PrecipitationHeatAdvectGlacier(noahmp);
  EnergyMainGlacier(noahmp);
  WaterMainGlacier(noahmp);
  BalanceWaterCheckGlacier(noahmp);
  BalanceEnergyCheckGlacier(noahmp);
}

}  // namespace cppwrf::phys::noahmp
