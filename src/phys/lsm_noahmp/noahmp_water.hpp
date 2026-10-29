#pragma once

#include <algorithm>
#include <cmath>
#include "core/backend.hpp"
#include "core/precision.hpp"
#include "phys/lsm_noahmp/noahmp_energy.hpp"
#include "phys/lsm_noahmp/noahmp_rad_therm.hpp"
#include "phys/lsm_noahmp/noahmp_var.hpp"

namespace cppwrf::phys::noahmp {

// ============================================================================
// AtmosForcingMod.F90
// ============================================================================
CPPWRF_HD inline void ProcessAtmosForcing(NoahmpType& noahmp) {
  constexpr int LoopNum = 10;
  real PrecipFrozenTot  = real(0.0);
  real RadDirFrac       = real(0.0);
  real RadVisFrac       = real(0.0);
  real VapPresSat       = real(0.0);
  real LatHeatVap       = real(0.0);
  real PsychConst       = real(0.0);
  real TemperatureDegC  = real(0.0);
  real TemperatureWetBulb = real(0.0);

  const real& CosSolarZenithAngle     = noahmp.config.domain.CosSolarZenithAngle;
  const int&  OptRainSnowPartition    = noahmp.config.nmlist.OptRainSnowPartition;
  const real& PressureAirRefHeight    = noahmp.forcing.PressureAirRefHeight;
  const real& TemperatureAirRefHeight = noahmp.forcing.TemperatureAirRefHeight;
  const real& SpecHumidityRefHeight   = noahmp.forcing.SpecHumidityRefHeight;
  const real& PrecipConvRefHeight     = noahmp.forcing.PrecipConvRefHeight;
  const real& PrecipNonConvRefHeight  = noahmp.forcing.PrecipNonConvRefHeight;
  const real& PrecipShConvRefHeight   = noahmp.forcing.PrecipShConvRefHeight;
  const real& PrecipSnowRefHeight     = noahmp.forcing.PrecipSnowRefHeight;
  const real& PrecipGraupelRefHeight  = noahmp.forcing.PrecipGraupelRefHeight;
  const real& PrecipHailRefHeight     = noahmp.forcing.PrecipHailRefHeight;
  real&       RadSwDownRefHeight      = noahmp.forcing.RadSwDownRefHeight;
  const real& WindEastwardRefHeight   = noahmp.forcing.WindEastwardRefHeight;
  const real& WindNorthwardRefHeight  = noahmp.forcing.WindNorthwardRefHeight;
  real&       RadSwVisFrac            = noahmp.forcing.RadSwVisFrac;
  real&       RadSwDirFrac            = noahmp.forcing.RadSwDirFrac;
  const real& SnowfallDensityMax      = noahmp.water.param.SnowfallDensityMax;
  real&       TemperaturePotRefHeight = noahmp.energy.state.TemperaturePotRefHeight;
  real&       PressureVaporRefHeight  = noahmp.energy.state.PressureVaporRefHeight;
  real&       DensityAirRefHeight     = noahmp.energy.state.DensityAirRefHeight;
  real&       WindSpdRefHeight        = noahmp.energy.state.WindSpdRefHeight;
  auto&       RadSwDownDir            = noahmp.energy.flux.RadSwDownDir;
  auto&       RadSwDownDif            = noahmp.energy.flux.RadSwDownDif;
  real&       RainfallRefHeight       = noahmp.water.flux.RainfallRefHeight;
  real&       SnowfallRefHeight       = noahmp.water.flux.SnowfallRefHeight;
  real&       PrecipTotRefHeight      = noahmp.water.flux.PrecipTotRefHeight;
  real&       PrecipConvTotRefHeight  = noahmp.water.flux.PrecipConvTotRefHeight;
  real&       PrecipLargeSclRefHeight = noahmp.water.flux.PrecipLargeSclRefHeight;
  real&       PrecipAreaFrac          = noahmp.water.state.PrecipAreaFrac;
  real&       FrozenPrecipFrac        = noahmp.water.state.FrozenPrecipFrac;
  real&       SnowfallDensity         = noahmp.water.state.SnowfallDensity;

  // surface air variables
  TemperaturePotRefHeight = TemperatureAirRefHeight *
                            rpow(PressureAirRefHeight / PressureAirRefHeight, ConstGasDryAir / ConstHeatCapacAir);
  PressureVaporRefHeight  = SpecHumidityRefHeight * PressureAirRefHeight / (real(0.622) + real(0.378) * SpecHumidityRefHeight);
  DensityAirRefHeight     = (PressureAirRefHeight - real(0.378) * PressureVaporRefHeight) /
                            (ConstGasDryAir * TemperatureAirRefHeight);

  // downward solar radiation
  RadDirFrac = real(0.7);
  RadVisFrac = real(0.5);

  if (RadSwDirFrac >= real(0.0) && RadSwDirFrac <= real(1.0)) {
    RadDirFrac = RadSwDirFrac;
  } else {
    RadSwDirFrac = RadDirFrac;
  }

  if (RadSwVisFrac >= real(0.0) && RadSwVisFrac <= real(1.0)) {
    RadVisFrac = RadSwVisFrac;
  } else {
    RadSwVisFrac = RadVisFrac;
  }

  if (CosSolarZenithAngle <= real(0.0)) RadSwDownRefHeight = real(0.0);
  RadSwDownDir(1) = RadSwDownRefHeight * RadDirFrac * RadVisFrac;
  RadSwDownDir(2) = RadSwDownRefHeight * RadDirFrac * (real(1.0) - RadVisFrac);
  RadSwDownDif(1) = RadSwDownRefHeight * (real(1.0) - RadDirFrac) * RadVisFrac;
  RadSwDownDif(2) = RadSwDownRefHeight * (real(1.0) - RadDirFrac) * (real(1.0) - RadVisFrac);

  // precipitation
  PrecipTotRefHeight = PrecipConvRefHeight + PrecipNonConvRefHeight + PrecipShConvRefHeight;
  if (OptRainSnowPartition == 4) {
    PrecipConvTotRefHeight  = PrecipConvRefHeight + PrecipShConvRefHeight;
    PrecipLargeSclRefHeight = PrecipNonConvRefHeight;
  } else {
    PrecipConvTotRefHeight  = real(0.10) * PrecipTotRefHeight;
    PrecipLargeSclRefHeight = real(0.90) * PrecipTotRefHeight;
  }

  // fractional area that receives precipitation (see, Niu et al. 2005)
  PrecipAreaFrac = real(0.0);
  if ((PrecipConvTotRefHeight + PrecipLargeSclRefHeight) > real(0.0)) {
    PrecipAreaFrac = (PrecipConvTotRefHeight + PrecipLargeSclRefHeight) /
                     (real(10.0) * PrecipConvTotRefHeight + PrecipLargeSclRefHeight);
  }

  // partition precipitation into rain and snow
  if (OptRainSnowPartition == 1) {
    if (TemperatureAirRefHeight > (ConstFreezePoint + real(2.5))) {
      FrozenPrecipFrac = real(0.0);
    } else {
      if (TemperatureAirRefHeight <= (ConstFreezePoint + real(0.5))) {
        FrozenPrecipFrac = real(1.0);
      } else if (TemperatureAirRefHeight <= (ConstFreezePoint + real(2.0))) {
        FrozenPrecipFrac = real(1.0) - (real(-54.632) + real(0.2) * TemperatureAirRefHeight);
      } else {
        FrozenPrecipFrac = real(0.6);
      }
    }
  }

  if (OptRainSnowPartition == 2) {
    if (TemperatureAirRefHeight >= (ConstFreezePoint + real(2.2))) {
      FrozenPrecipFrac = real(0.0);
    } else {
      FrozenPrecipFrac = real(1.0);
    }
  }

  if (OptRainSnowPartition == 3) {
    if (TemperatureAirRefHeight >= ConstFreezePoint) {
      FrozenPrecipFrac = real(0.0);
    } else {
      FrozenPrecipFrac = real(1.0);
    }
  }

  SnowfallDensity = std::min(SnowfallDensityMax,
                             real(67.92) + real(51.25) * rexp((TemperatureAirRefHeight - ConstFreezePoint) / real(2.59)));
  if (OptRainSnowPartition == 4) {
    PrecipFrozenTot = PrecipSnowRefHeight + PrecipGraupelRefHeight + PrecipHailRefHeight;
    if ((PrecipNonConvRefHeight > real(0.0)) && (PrecipFrozenTot > real(0.0))) {
      FrozenPrecipFrac = std::min(real(1.0), PrecipFrozenTot / PrecipNonConvRefHeight);
      FrozenPrecipFrac = std::max(real(0.0), FrozenPrecipFrac);
      SnowfallDensity  = SnowfallDensity     * (PrecipSnowRefHeight / PrecipFrozenTot) +
                         ConstDensityGraupel * (PrecipGraupelRefHeight / PrecipFrozenTot) +
                         ConstDensityHail    * (PrecipHailRefHeight / PrecipFrozenTot);
    } else {
      FrozenPrecipFrac = real(0.0);
    }
  }

  if (OptRainSnowPartition == 5) {
    if (TemperatureAirRefHeight >= (ConstFreezePoint + real(8.0))) {
      FrozenPrecipFrac = real(0.0);
    } else {
      TemperatureDegC = std::min(real(50.0), std::max(real(-50.0), (TemperatureAirRefHeight - ConstFreezePoint)));
      if (TemperatureAirRefHeight > ConstFreezePoint) {
        LatHeatVap = ConstLatHeatEvap;
      } else {
        LatHeatVap = ConstLatHeatSublim;
      }
      PsychConst         = ConstHeatCapacAir * PressureAirRefHeight / (real(0.622) * LatHeatVap);
      TemperatureWetBulb = TemperatureDegC - real(5.0);
      for (int LoopInd = 1; LoopInd <= LoopNum; ++LoopInd) {
        VapPresSat         = real(610.8) * rexp((real(17.27) * TemperatureWetBulb) / (real(237.3) + TemperatureWetBulb));
        TemperatureWetBulb = TemperatureWetBulb - (VapPresSat - PressureVaporRefHeight) / PsychConst;
      }
      FrozenPrecipFrac = real(1.0) / (real(1.0) + real(6.99e-5) * rexp(real(2.0) * (TemperatureWetBulb + real(3.97))));
    }
  }

  // rain-snow partitioning
  RainfallRefHeight = PrecipTotRefHeight * (real(1.0) - FrozenPrecipFrac);
  SnowfallRefHeight = PrecipTotRefHeight * FrozenPrecipFrac;

  // wind speed at reference height for turbulence calculation
  WindSpdRefHeight = std::max(std::sqrt(rpow(WindEastwardRefHeight, real(2.0)) + rpow(WindNorthwardRefHeight, real(2.0))), real(1.0));
}

// ============================================================================
// GeneralInitMod.F90
// ============================================================================
CPPWRF_HD inline void GeneralInit(NoahmpType& noahmp) {
  const int&  IsModisLanduse         = noahmp.config.domain.IsModisLanduse;
  const int&  VegType                = noahmp.config.domain.VegType;
  const int&  NumSoilLayer           = noahmp.config.domain.NumSoilLayer;
  const auto& DepthSoilLayer         = noahmp.config.domain.DepthSoilLayer;
  const int&  NumSoilLayerRoot       = noahmp.water.param.NumSoilLayerRoot;
  const int&  NumSnowLayerNeg        = noahmp.config.domain.NumSnowLayerNeg;
  const auto& DepthSnowSoilLayer     = noahmp.config.domain.DepthSnowSoilLayer;
  const auto& TemperatureSoilSnow    = noahmp.energy.state.TemperatureSoilSnow;
  bool&       FlagCropland           = noahmp.config.domain.FlagCropland;
  bool&       FlagWetland            = noahmp.config.domain.FlagWetland;
  auto&       ThicknessSnowSoilLayer = noahmp.config.domain.ThicknessSnowSoilLayer;
  real&       TemperatureRootZone    = noahmp.energy.state.TemperatureRootZone;

  for (int LoopInd = NumSnowLayerNeg + 1; LoopInd <= NumSoilLayer; ++LoopInd) {
    if (LoopInd == NumSnowLayerNeg + 1) {
      ThicknessSnowSoilLayer(LoopInd) = -DepthSnowSoilLayer(LoopInd);
    } else {
      ThicknessSnowSoilLayer(LoopInd) = DepthSnowSoilLayer(LoopInd - 1) - DepthSnowSoilLayer(LoopInd);
    }
  }

  TemperatureRootZone = real(0.0);
  for (int LoopInd = 1; LoopInd <= NumSoilLayerRoot; ++LoopInd) {
    TemperatureRootZone = TemperatureRootZone +
                          TemperatureSoilSnow(LoopInd) * ThicknessSnowSoilLayer(LoopInd) / (-DepthSoilLayer(NumSoilLayerRoot));
  }

  FlagCropland = false;
  FlagWetland  = false;
  if (IsModisLanduse == 0) {
    if ((VegType >= 3) && (VegType <= 6))   FlagCropland = true;
    if ((VegType >= 17) && (VegType <= 18)) FlagWetland  = true;
  } else if (IsModisLanduse == 1) {
    if ((VegType == 12) || (VegType == 14)) FlagCropland = true;
    if (VegType == 11)                      FlagWetland  = true;
  }
}

// ============================================================================
// PhenologyMainMod.F90
// ============================================================================
CPPWRF_HD inline void PhenologyMain(NoahmpType& noahmp) {
  int  IntpMonth1 = 0, IntpMonth2 = 0;
  real ThicknessCanBury = real(0.0);
  real SnowDepthVegBury = real(0.0);
  real DayCurrent       = real(0.0);
  real IntpWgt1         = real(0.0);
  real IntpWgt2         = real(0.0);
  real MonthCurrent     = real(0.0);

  const int&  OptDynamicVeg          = noahmp.config.nmlist.OptDynamicVeg;
  const int&  OptCropModel           = noahmp.config.nmlist.OptCropModel;
  const int&  VegType                = noahmp.config.domain.VegType;
  const int&  CropType               = noahmp.config.domain.CropType;
  const int&  IndexIcePoint          = noahmp.config.domain.IndexIcePoint;
  const int&  IndexBarrenPoint       = noahmp.config.domain.IndexBarrenPoint;
  const int&  IndexWaterPoint        = noahmp.config.domain.IndexWaterPoint;
  const bool& FlagUrban              = noahmp.config.domain.FlagUrban;
  bool&       FlagDynamicVeg         = noahmp.config.domain.FlagDynamicVeg;
  bool&       FlagDynamicCrop        = noahmp.config.domain.FlagDynamicCrop;
  const real& Latitude               = noahmp.config.domain.Latitude;
  const int&  NumDayInYear           = noahmp.config.domain.NumDayInYear;
  const real& DayJulianInYear        = noahmp.config.domain.DayJulianInYear;
  const real& HeightCanopyTop        = noahmp.energy.param.HeightCanopyTop;
  const real& HeightCanopyBot        = noahmp.energy.param.HeightCanopyBot;
  const auto& LeafAreaIndexMon       = noahmp.energy.param.LeafAreaIndexMon;
  const auto& StemAreaIndexMon       = noahmp.energy.param.StemAreaIndexMon;
  const real& VegFracAnnMax          = noahmp.energy.param.VegFracAnnMax;
  const real& VegFracGreen           = noahmp.energy.param.VegFracGreen;
  const real& TemperatureMinPhotosyn = noahmp.biochem.param.TemperatureMinPhotosyn;
  const int&  PlantGrowStage         = noahmp.biochem.state.PlantGrowStage;
  const real& SnowDepth              = noahmp.water.state.SnowDepth;
  const real& TemperatureCanopy      = noahmp.energy.state.TemperatureCanopy;
  real&       LeafAreaIndex          = noahmp.energy.state.LeafAreaIndex;
  real&       StemAreaIndex          = noahmp.energy.state.StemAreaIndex;
  real&       LeafAreaIndEff         = noahmp.energy.state.LeafAreaIndEff;
  real&       StemAreaIndEff         = noahmp.energy.state.StemAreaIndEff;
  real&       VegFrac                = noahmp.energy.state.VegFrac;
  real&       CanopyFracSnowBury     = noahmp.energy.state.CanopyFracSnowBury;
  real&       IndexGrowSeason        = noahmp.biochem.state.IndexGrowSeason;

  if (CropType == 0) {
    if ((OptDynamicVeg == 1) || (OptDynamicVeg == 3) || (OptDynamicVeg == 4)) {
      if (Latitude >= real(0.0)) {
        DayCurrent = DayJulianInYear;
      } else {
        DayCurrent = std::fmod(DayJulianInYear + (real(0.5) * static_cast<real>(NumDayInYear)),
                               static_cast<real>(static_cast<float>(NumDayInYear)));
      }
      MonthCurrent = real(12.0) * DayCurrent / static_cast<real>(static_cast<float>(NumDayInYear));
      IntpMonth1   = static_cast<int>(MonthCurrent + real(0.5));
      IntpMonth2   = IntpMonth1 + 1;
      IntpWgt1     = (static_cast<real>(IntpMonth1) + real(0.5)) - MonthCurrent;
      IntpWgt2     = real(1.0) - IntpWgt1;
      if (IntpMonth1 < 1)  IntpMonth1 = 12;
      if (IntpMonth2 > 12) IntpMonth2 = 1;
      LeafAreaIndex = IntpWgt1 * LeafAreaIndexMon(IntpMonth1) + IntpWgt2 * LeafAreaIndexMon(IntpMonth2);
      StemAreaIndex = IntpWgt1 * StemAreaIndexMon(IntpMonth1) + IntpWgt2 * StemAreaIndexMon(IntpMonth2);
    }

    if ((OptDynamicVeg == 7) || (OptDynamicVeg == 8) || (OptDynamicVeg == 9)) {
      StemAreaIndex = std::max(real(0.05), real(0.1) * LeafAreaIndex);
      if (LeafAreaIndex < real(0.05)) StemAreaIndex = real(0.0);
    }
    if (StemAreaIndex < real(0.05)) StemAreaIndex = real(0.0);
    if ((LeafAreaIndex < real(0.05)) || (StemAreaIndex == real(0.0))) LeafAreaIndex = real(0.0);

    if ((VegType == IndexWaterPoint) || (VegType == IndexBarrenPoint) ||
        (VegType == IndexIcePoint) || FlagUrban) {
      LeafAreaIndex = real(0.0);
      StemAreaIndex = real(0.0);
    }
  }

  ThicknessCanBury   = std::min(std::max(SnowDepth - HeightCanopyBot, real(0.0)), (HeightCanopyTop - HeightCanopyBot));
  CanopyFracSnowBury = ThicknessCanBury / std::max(real(1.0e-06), (HeightCanopyTop - HeightCanopyBot));
  if ((HeightCanopyTop > real(0.0)) && (HeightCanopyTop <= real(1.0))) {
    SnowDepthVegBury   = HeightCanopyTop * rexp(-std::min(SnowDepth, real(10.0)) / real(0.2));
    CanopyFracSnowBury = std::min(SnowDepth, SnowDepthVegBury) / SnowDepthVegBury;
  }

  LeafAreaIndEff = LeafAreaIndex * (real(1.0) - CanopyFracSnowBury);
  StemAreaIndEff = StemAreaIndex * (real(1.0) - CanopyFracSnowBury);
  if ((StemAreaIndEff < real(0.05)) && (CropType == 0)) StemAreaIndEff = real(0.0);
  if (((LeafAreaIndEff < real(0.05)) || (StemAreaIndEff == real(0.0))) && (CropType == 0)) {
    LeafAreaIndEff = real(0.0);
  }

  if (((TemperatureCanopy > TemperatureMinPhotosyn) && (CropType == 0)) ||
      ((PlantGrowStage > 2) && (PlantGrowStage < 7) && (CropType > 0))) {
    IndexGrowSeason = real(1.0);
  } else {
    IndexGrowSeason = real(0.0);
  }

  if ((OptDynamicVeg == 1) || (OptDynamicVeg == 6) || (OptDynamicVeg == 7)) {
    VegFrac = VegFracGreen;
  } else if ((OptDynamicVeg == 2) || (OptDynamicVeg == 3) || (OptDynamicVeg == 8)) {
    VegFrac = real(1.0) - rexp(real(-0.52) * (LeafAreaIndex + StemAreaIndex));
  } else if ((OptDynamicVeg == 4) || (OptDynamicVeg == 5) || (OptDynamicVeg == 9)) {
    VegFrac = VegFracAnnMax;
  }

  if ((OptCropModel > 0) && (CropType > 0)) {
    VegFrac = VegFracAnnMax;
  }

  if (VegFrac <= real(0.05)) VegFrac = real(0.05);
  if (FlagUrban || (VegType == IndexBarrenPoint)) VegFrac = real(0.0);
  if ((LeafAreaIndEff + StemAreaIndEff) == real(0.0)) VegFrac = real(0.0);

  FlagDynamicCrop = false;
  FlagDynamicVeg  = false;
  if ((OptDynamicVeg == 2) || (OptDynamicVeg == 5) || (OptDynamicVeg == 6)) {
    FlagDynamicVeg = true;
  }
  if ((OptCropModel > 0) && (CropType > 0)) {
    FlagDynamicCrop = true;
    FlagDynamicVeg  = false;
  }
}

// ============================================================================
// IrrigationPrepareMod.F90 (for OptIrrigation == 0 / default)
// ============================================================================
CPPWRF_HD inline void IrrigationPrepare(NoahmpType& noahmp) {
  const bool& FlagSoilProcess         = noahmp.config.domain.FlagSoilProcess;
  const int&  OptIrrigationMethod     = noahmp.config.nmlist.OptIrrigationMethod;
  const real& IrriFracThreshold       = noahmp.water.param.IrriFracThreshold;
  const real& IrriStopPrecipThr       = noahmp.water.param.IrriStopPrecipThr;
  const real& IrrigationFracGrid      = noahmp.water.state.IrrigationFracGrid;
  real&       IrrigationAmtSprinkler  = noahmp.water.state.IrrigationAmtSprinkler;
  real&       IrrigationAmtFlood      = noahmp.water.state.IrrigationAmtFlood;
  real&       IrrigationAmtMicro      = noahmp.water.state.IrrigationAmtMicro;
  const real& RainfallRefHeight       = noahmp.water.flux.RainfallRefHeight;
  real&       IrrigationFracSprinkler = noahmp.water.state.IrrigationFracSprinkler;
  real&       IrrigationFracMicro     = noahmp.water.state.IrrigationFracMicro;
  real&       IrrigationFracFlood     = noahmp.water.state.IrrigationFracFlood;

  if ((OptIrrigationMethod == 0) && (IrrigationFracSprinkler == real(0.0)) && (IrrigationFracMicro == real(0.0)) &&
      (IrrigationFracFlood == real(0.0)) && (IrrigationFracGrid >= IrriFracThreshold)) {
    IrrigationFracSprinkler = real(1.0);
  }

  if (OptIrrigationMethod == 1) {
    IrrigationFracSprinkler = real(1.0);
    IrrigationFracMicro     = real(0.0);
    IrrigationFracFlood     = real(0.0);
  } else if (OptIrrigationMethod == 2) {
    IrrigationFracSprinkler = real(0.0);
    IrrigationFracMicro     = real(1.0);
    IrrigationFracFlood     = real(0.0);
  } else if (OptIrrigationMethod == 3) {
    IrrigationFracSprinkler = real(0.0);
    IrrigationFracMicro     = real(0.0);
    IrrigationFracFlood     = real(1.0);
  }

  if (FlagSoilProcess) {
    if ((RainfallRefHeight >= (IrriStopPrecipThr / real(3600.0))) || (IrrigationFracGrid < IrriFracThreshold)) {
      IrrigationAmtSprinkler = real(0.0);
      IrrigationAmtMicro     = real(0.0);
      IrrigationAmtFlood     = real(0.0);
    }
  }
}

// ============================================================================
// CanopyWaterInterceptMod.F90
// ============================================================================
CPPWRF_HD inline void CanopyWaterIntercept(NoahmpType& noahmp) {
  real IceDripFacTemp = real(0.0);
  real IceDripFacWind = real(0.0);
  real CanopySnowDrip = real(0.0);

  const int&  SurfaceType            = noahmp.config.domain.SurfaceType;
  const real& MainTimeStep           = noahmp.config.domain.MainTimeStep;
  const real& WindEastwardRefHeight  = noahmp.forcing.WindEastwardRefHeight;
  const real& WindNorthwardRefHeight = noahmp.forcing.WindNorthwardRefHeight;
  const real& LeafAreaIndEff         = noahmp.energy.state.LeafAreaIndEff;
  const real& StemAreaIndEff         = noahmp.energy.state.StemAreaIndEff;
  const real& VegFrac                = noahmp.energy.state.VegFrac;
  const real& TemperatureCanopy      = noahmp.energy.state.TemperatureCanopy;
  const real& TemperatureGrd         = noahmp.energy.state.TemperatureGrd;
  const real& CanopyLiqHoldCap       = noahmp.water.param.CanopyLiqHoldCap;
  const real& RainfallRefHeight      = noahmp.water.flux.RainfallRefHeight;
  const real& SnowfallRefHeight      = noahmp.water.flux.SnowfallRefHeight;
  const real& SnowfallDensity        = noahmp.water.state.SnowfallDensity;
  const real& PrecipAreaFrac         = noahmp.water.state.PrecipAreaFrac;
  real&       CanopyLiqWater         = noahmp.water.state.CanopyLiqWater;
  real&       CanopyIce              = noahmp.water.state.CanopyIce;
  real&       CanopyWetFrac          = noahmp.water.state.CanopyWetFrac;
  real&       CanopyTotalWater       = noahmp.water.state.CanopyTotalWater;
  real&       CanopyIceMax           = noahmp.water.state.CanopyIceMax;
  real&       CanopyLiqWaterMax      = noahmp.water.state.CanopyLiqWaterMax;
  real&       InterceptCanopyRain    = noahmp.water.flux.InterceptCanopyRain;
  real&       DripCanopyRain         = noahmp.water.flux.DripCanopyRain;
  real&       ThroughfallRain        = noahmp.water.flux.ThroughfallRain;
  real&       InterceptCanopySnow    = noahmp.water.flux.InterceptCanopySnow;
  real&       DripCanopySnow         = noahmp.water.flux.DripCanopySnow;
  real&       ThroughfallSnow        = noahmp.water.flux.ThroughfallSnow;
  real&       RainfallGround         = noahmp.water.flux.RainfallGround;
  real&       SnowfallGround         = noahmp.water.flux.SnowfallGround;
  real&       SnowDepthIncr          = noahmp.water.flux.SnowDepthIncr;

  InterceptCanopyRain = real(0.0);
  DripCanopyRain      = real(0.0);
  ThroughfallRain     = real(0.0);
  InterceptCanopySnow = real(0.0);
  DripCanopySnow      = real(0.0);
  ThroughfallSnow     = real(0.0);
  RainfallGround      = real(0.0);
  SnowfallGround      = real(0.0);
  SnowDepthIncr       = real(0.0);

  // canopy liquid water
  CanopyLiqWaterMax = VegFrac * CanopyLiqHoldCap * (LeafAreaIndEff + StemAreaIndEff);

  if ((LeafAreaIndEff + StemAreaIndEff) > real(0.0)) {
    InterceptCanopyRain = VegFrac * RainfallRefHeight * PrecipAreaFrac;
    InterceptCanopyRain = std::min(InterceptCanopyRain,
                                   (CanopyLiqWaterMax - CanopyLiqWater) / MainTimeStep *
                                   (real(1.0) - rexp(-RainfallRefHeight * MainTimeStep / CanopyLiqWaterMax)));
    InterceptCanopyRain = std::max(InterceptCanopyRain, real(0.0));
    DripCanopyRain      = VegFrac * RainfallRefHeight - InterceptCanopyRain;
    ThroughfallRain     = (real(1.0) - VegFrac) * RainfallRefHeight;
    CanopyLiqWater      = std::max(real(0.0), CanopyLiqWater + InterceptCanopyRain * MainTimeStep);
  } else {
    InterceptCanopyRain = real(0.0);
    DripCanopyRain      = real(0.0);
    ThroughfallRain     = RainfallRefHeight;
    if (CanopyLiqWater > real(0.0)) {
      DripCanopyRain = DripCanopyRain + CanopyLiqWater / MainTimeStep;
      CanopyLiqWater = real(0.0);
    }
  }

  // canopy ice
  CanopyIceMax = VegFrac * real(6.6) * (real(0.27) + real(46.0) / SnowfallDensity) * (LeafAreaIndEff + StemAreaIndEff);

  if ((LeafAreaIndEff + StemAreaIndEff) > real(0.0)) {
    InterceptCanopySnow = VegFrac * SnowfallRefHeight * PrecipAreaFrac;
    InterceptCanopySnow = std::min(InterceptCanopySnow,
                                   (CanopyIceMax - CanopyIce) / MainTimeStep *
                                   (real(1.0) - rexp(-SnowfallRefHeight * MainTimeStep / CanopyIceMax)));
    InterceptCanopySnow = std::max(InterceptCanopySnow, real(0.0));
    IceDripFacTemp      = std::max(real(0.0), (TemperatureCanopy - real(270.15)) / real(1.87e5));
    IceDripFacWind      = std::sqrt(rpow(WindEastwardRefHeight, real(2.0)) + rpow(WindNorthwardRefHeight, real(2.0))) / real(1.56e5);
    CanopySnowDrip      = std::max(real(0.0), CanopyIce) * (IceDripFacWind + IceDripFacTemp);
    CanopySnowDrip      = std::min(CanopyIce / MainTimeStep + InterceptCanopySnow, CanopySnowDrip);
    DripCanopySnow      = (VegFrac * SnowfallRefHeight - InterceptCanopySnow) + CanopySnowDrip;
    ThroughfallSnow     = (real(1.0) - VegFrac) * SnowfallRefHeight;
    CanopyIce           = std::max(real(0.0), CanopyIce + (InterceptCanopySnow - CanopySnowDrip) * MainTimeStep);
  } else {
    InterceptCanopySnow = real(0.0);
    DripCanopySnow      = real(0.0);
    ThroughfallSnow     = SnowfallRefHeight;
    if (CanopyIce > real(0.0)) {
      DripCanopySnow = DripCanopySnow + CanopyIce / MainTimeStep;
      CanopyIce      = real(0.0);
    }
  }

  if (CanopyIce > real(0.0)) {
    CanopyWetFrac = std::max(real(0.0), CanopyIce) / std::max(CanopyIceMax, real(1.0e-06));
  } else {
    CanopyWetFrac = std::max(real(0.0), CanopyLiqWater) / std::max(CanopyLiqWaterMax, real(1.0e-06));
  }
  CanopyWetFrac    = rpow(std::min(CanopyWetFrac, real(1.0)), real(0.667));
  CanopyTotalWater = CanopyLiqWater + CanopyIce;

  RainfallGround = DripCanopyRain + ThroughfallRain;
  SnowfallGround = DripCanopySnow + ThroughfallSnow;
  SnowDepthIncr  = SnowfallGround / SnowfallDensity;
  if ((SurfaceType == 2) && (TemperatureGrd > ConstFreezePoint)) {
    SnowfallGround = real(0.0);
    SnowDepthIncr  = real(0.0);
  }
}

// ============================================================================
// PrecipitationHeatAdvectMod.F90
// ============================================================================
CPPWRF_HD inline void PrecipitationHeatAdvect(NoahmpType& noahmp) {
  real HeatPrcpAirToCan = real(0.0);
  real HeatPrcpCanToGrd = real(0.0);
  real HeatPrcpAirToGrd = real(0.0);

  const real& TemperatureAirRefHeight = noahmp.forcing.TemperatureAirRefHeight;
  const real& TemperatureCanopy       = noahmp.energy.state.TemperatureCanopy;
  const real& TemperatureGrd          = noahmp.energy.state.TemperatureGrd;
  const real& VegFrac                 = noahmp.energy.state.VegFrac;
  const real& RainfallRefHeight       = noahmp.water.flux.RainfallRefHeight;
  const real& SnowfallRefHeight       = noahmp.water.flux.SnowfallRefHeight;
  const real& DripCanopyRain          = noahmp.water.flux.DripCanopyRain;
  const real& ThroughfallRain         = noahmp.water.flux.ThroughfallRain;
  const real& DripCanopySnow          = noahmp.water.flux.DripCanopySnow;
  const real& ThroughfallSnow         = noahmp.water.flux.ThroughfallSnow;
  real&       HeatPrecipAdvCanopy     = noahmp.energy.flux.HeatPrecipAdvCanopy;
  real&       HeatPrecipAdvVegGrd     = noahmp.energy.flux.HeatPrecipAdvVegGrd;
  real&       HeatPrecipAdvBareGrd    = noahmp.energy.flux.HeatPrecipAdvBareGrd;

  HeatPrecipAdvCanopy  = real(0.0);
  HeatPrecipAdvVegGrd  = real(0.0);
  HeatPrecipAdvBareGrd = real(0.0);

  HeatPrcpAirToCan = VegFrac * RainfallRefHeight * (ConstHeatCapacWater / real(1000.0)) * (TemperatureAirRefHeight - TemperatureCanopy);
  HeatPrcpCanToGrd = DripCanopyRain * (ConstHeatCapacWater / real(1000.0)) * (TemperatureCanopy - TemperatureGrd);
  HeatPrcpAirToGrd = ThroughfallRain * (ConstHeatCapacWater / real(1000.0)) * (TemperatureAirRefHeight - TemperatureGrd);

  HeatPrcpAirToCan = HeatPrcpAirToCan +
                     VegFrac * SnowfallRefHeight * (ConstHeatCapacIce / real(1000.0)) * (TemperatureAirRefHeight - TemperatureCanopy);
  HeatPrcpCanToGrd = HeatPrcpCanToGrd +
                     DripCanopySnow * (ConstHeatCapacIce / real(1000.0)) * (TemperatureCanopy - TemperatureGrd);
  HeatPrcpAirToGrd = HeatPrcpAirToGrd +
                     ThroughfallSnow * (ConstHeatCapacIce / real(1000.0)) * (TemperatureAirRefHeight - TemperatureGrd);

  HeatPrecipAdvCanopy  = HeatPrcpAirToCan - HeatPrcpCanToGrd;
  HeatPrecipAdvVegGrd  = HeatPrcpCanToGrd;
  HeatPrecipAdvBareGrd = HeatPrcpAirToGrd;

  if ((VegFrac > real(0.0)) && (VegFrac < real(1.0))) {
    HeatPrecipAdvVegGrd  = HeatPrecipAdvVegGrd / VegFrac;
    HeatPrecipAdvBareGrd = HeatPrecipAdvBareGrd / (real(1.0) - VegFrac);
  } else if (VegFrac <= real(0.0)) {
    HeatPrecipAdvBareGrd = HeatPrecipAdvVegGrd + HeatPrecipAdvBareGrd;
    HeatPrecipAdvVegGrd  = real(0.0);
    HeatPrecipAdvCanopy  = real(0.0);
  } else if (VegFrac >= real(1.0)) {
    HeatPrecipAdvBareGrd = real(0.0);
  }

  HeatPrecipAdvCanopy  = std::max(HeatPrecipAdvCanopy,  real(-20.0));
  HeatPrecipAdvCanopy  = std::min(HeatPrecipAdvCanopy,  real( 20.0));
  HeatPrecipAdvVegGrd  = std::max(HeatPrecipAdvVegGrd,  real(-20.0));
  HeatPrecipAdvVegGrd  = std::min(HeatPrecipAdvVegGrd,  real( 20.0));
  HeatPrecipAdvBareGrd = std::max(HeatPrecipAdvBareGrd, real(-20.0));
  HeatPrecipAdvBareGrd = std::min(HeatPrecipAdvBareGrd, real( 20.0));
}

// ============================================================================
// CanopyHydrologyMod.F90
// ============================================================================
CPPWRF_HD inline void CanopyHydrology(NoahmpType& noahmp) {
  const real& MainTimeStep      = noahmp.config.domain.MainTimeStep;
  const real& HeatLatentCanopy  = noahmp.energy.flux.HeatLatentCanopy;
  const real& HeatLatentTransp  = noahmp.energy.flux.HeatLatentTransp;
  const real& LeafAreaIndEff    = noahmp.energy.state.LeafAreaIndEff;
  const real& StemAreaIndEff    = noahmp.energy.state.StemAreaIndEff;
  const bool& FlagFrozenCanopy  = noahmp.energy.state.FlagFrozenCanopy;
  const real& VegFrac           = noahmp.energy.state.VegFrac;
  const real& SnowfallDensity   = noahmp.water.state.SnowfallDensity;
  const real& CanopyLiqHoldCap  = noahmp.water.param.CanopyLiqHoldCap;
  real&       CanopyLiqWater    = noahmp.water.state.CanopyLiqWater;
  real&       CanopyIce         = noahmp.water.state.CanopyIce;
  real&       TemperatureCanopy = noahmp.energy.state.TemperatureCanopy;
  real&       CanopyTotalWater  = noahmp.water.state.CanopyTotalWater;
  real&       CanopyWetFrac     = noahmp.water.state.CanopyWetFrac;
  real&       CanopyIceMax      = noahmp.water.state.CanopyIceMax;
  real&       CanopyLiqWaterMax = noahmp.water.state.CanopyLiqWaterMax;
  real&       EvapCanopyNet     = noahmp.water.flux.EvapCanopyNet;
  real&       Transpiration     = noahmp.water.flux.Transpiration;
  real&       EvapCanopyLiq     = noahmp.water.flux.EvapCanopyLiq;
  real&       DewCanopyLiq      = noahmp.water.flux.DewCanopyLiq;
  real&       FrostCanopyIce    = noahmp.water.flux.FrostCanopyIce;
  real&       SublimCanopyIce   = noahmp.water.flux.SublimCanopyIce;
  real&       MeltCanopyIce     = noahmp.water.flux.MeltCanopyIce;
  real&       FreezeCanopyLiq   = noahmp.water.flux.FreezeCanopyLiq;

  EvapCanopyNet     = real(0.0);
  Transpiration     = real(0.0);
  EvapCanopyLiq     = real(0.0);
  DewCanopyLiq      = real(0.0);
  FrostCanopyIce    = real(0.0);
  SublimCanopyIce   = real(0.0);
  MeltCanopyIce     = real(0.0);
  FreezeCanopyLiq   = real(0.0);
  CanopyLiqWaterMax = real(0.0);
  CanopyIceMax      = real(0.0);
  CanopyWetFrac     = real(0.0);
  CanopyTotalWater  = real(0.0);

  CanopyLiqWaterMax = VegFrac * CanopyLiqHoldCap * (LeafAreaIndEff + StemAreaIndEff);

  // max(x, 0.0): gfortran's maxsd(x, 0.0) = x > 0 ? x : 0 (+0.0 for x = -0.0; std::max keeps -0.0)
  if (!FlagFrozenCanopy) {
    Transpiration   = HeatLatentTransp / ConstLatHeatEvap > real(0.0) ? HeatLatentTransp / ConstLatHeatEvap : real(0.0);
    EvapCanopyLiq   = HeatLatentCanopy / ConstLatHeatEvap > real(0.0) ? HeatLatentCanopy / ConstLatHeatEvap : real(0.0);
    DewCanopyLiq    = std::abs(std::min(HeatLatentCanopy / ConstLatHeatEvap, real(0.0)));
    SublimCanopyIce = real(0.0);
    FrostCanopyIce  = real(0.0);
  } else {
    Transpiration   = HeatLatentTransp / ConstLatHeatSublim > real(0.0) ? HeatLatentTransp / ConstLatHeatSublim : real(0.0);
    EvapCanopyLiq   = real(0.0);
    DewCanopyLiq    = real(0.0);
    SublimCanopyIce = HeatLatentCanopy / ConstLatHeatSublim > real(0.0) ? HeatLatentCanopy / ConstLatHeatSublim : real(0.0);
    FrostCanopyIce  = std::abs(std::min(HeatLatentCanopy / ConstLatHeatSublim, real(0.0)));
  }

  EvapCanopyLiq  = std::min(CanopyLiqWater / MainTimeStep, EvapCanopyLiq);
  CanopyLiqWater = std::max(real(0.0), CanopyLiqWater + (DewCanopyLiq - EvapCanopyLiq) * MainTimeStep);
  if (CanopyLiqWater <= real(1.0e-06)) CanopyLiqWater = real(0.0);

  CanopyIceMax = VegFrac * real(6.6) * (real(0.27) + real(46.0) / SnowfallDensity) * (LeafAreaIndEff + StemAreaIndEff);

  SublimCanopyIce = std::min(CanopyIce / MainTimeStep, SublimCanopyIce);
  CanopyIce       = std::max(real(0.0), CanopyIce + (FrostCanopyIce - SublimCanopyIce) * MainTimeStep);
  if (CanopyIce <= real(1.0e-6)) CanopyIce = real(0.0);

  if ((CanopyIce > real(0.0)) && (CanopyIce >= CanopyLiqWater)) {
    CanopyWetFrac = std::max(real(0.0), CanopyIce) / std::max(CanopyIceMax, real(1.0e-06));
  } else {
    CanopyWetFrac = std::max(real(0.0), CanopyLiqWater) / std::max(CanopyLiqWaterMax, real(1.0e-06));
  }
  CanopyWetFrac    = rpow(std::min(CanopyWetFrac, real(1.0)), real(0.667));
  CanopyTotalWater = CanopyLiqWater + CanopyIce;

  if ((CanopyIce > real(1.0e-6)) && (TemperatureCanopy > ConstFreezePoint)) {
    MeltCanopyIce     = std::min(CanopyIce / MainTimeStep,
                                 (TemperatureCanopy - ConstFreezePoint) * ConstHeatCapacIce *
                                 CanopyIce / ConstDensityIce / (MainTimeStep * ConstLatHeatFusion));
    CanopyIce         = std::max(real(0.0), CanopyIce - MeltCanopyIce * MainTimeStep);
    CanopyLiqWater    = std::max(real(0.0), CanopyTotalWater - CanopyIce);
    TemperatureCanopy = CanopyWetFrac * ConstFreezePoint + (real(1.0) - CanopyWetFrac) * TemperatureCanopy;
  }

  if ((CanopyLiqWater > real(1.0e-6)) && (TemperatureCanopy < ConstFreezePoint)) {
    FreezeCanopyLiq   = std::min(CanopyLiqWater / MainTimeStep,
                                 (ConstFreezePoint - TemperatureCanopy) * ConstHeatCapacWater *
                                 CanopyLiqWater / ConstDensityWater / (MainTimeStep * ConstLatHeatFusion));
    CanopyLiqWater    = std::max(real(0.0), CanopyLiqWater - FreezeCanopyLiq * MainTimeStep);
    CanopyIce         = std::max(real(0.0), CanopyTotalWater - CanopyLiqWater);
    TemperatureCanopy = CanopyWetFrac * ConstFreezePoint + (real(1.0) - CanopyWetFrac) * TemperatureCanopy;
  }

  CanopyTotalWater = CanopyLiqWater + CanopyIce;
  EvapCanopyNet    = EvapCanopyLiq + SublimCanopyIce - DewCanopyLiq - FrostCanopyIce;
}

// ============================================================================
// SnowLayerWaterComboMod.F90
// ============================================================================
CPPWRF_HD inline void SnowLayerWaterCombo(real& ThickLayer1,
                                          real& LiqLayer1,
                                          real& IceLayer1,
                                          real& TempLayer1,
                                          real  ThickLayer2,
                                          real  LiqLayer2,
                                          real  IceLayer2,
                                          real  TempLayer2) {
  real ThickLayerComb = ThickLayer1 + ThickLayer2;
  real IceLayerComb   = IceLayer1 + IceLayer2;
  real LiqLayerComb   = LiqLayer1 + LiqLayer2;
  real EnthLayer1     = (ConstHeatCapacIce * IceLayer1 + ConstHeatCapacWater * LiqLayer1) *
                        (TempLayer1 - ConstFreezePoint) + ConstLatHeatFusion * LiqLayer1;
  real EnthLayer2     = (ConstHeatCapacIce * IceLayer2 + ConstHeatCapacWater * LiqLayer2) *
                        (TempLayer2 - ConstFreezePoint) + ConstLatHeatFusion * LiqLayer2;
  real EnthLayerComb  = EnthLayer1 + EnthLayer2;
  real TempLayerComb  = real(0.0);

  if (EnthLayerComb < real(0.0)) {
    TempLayerComb = ConstFreezePoint + EnthLayerComb /
                    (ConstHeatCapacIce * IceLayerComb + ConstHeatCapacWater * LiqLayerComb);
  } else if (EnthLayerComb <= (ConstLatHeatFusion * LiqLayerComb)) {
    TempLayerComb = ConstFreezePoint;
  } else {
    TempLayerComb = ConstFreezePoint + (EnthLayerComb - ConstLatHeatFusion * LiqLayerComb) /
                    (ConstHeatCapacIce * IceLayerComb + ConstHeatCapacWater * LiqLayerComb);
  }

  ThickLayer1 = ThickLayerComb;
  IceLayer1   = IceLayerComb;
  LiqLayer1   = LiqLayerComb;
  TempLayer1  = TempLayerComb;
}

// ============================================================================
// SnowfallBelowCanopyMod.F90
// ============================================================================
CPPWRF_HD inline void SnowfallAfterCanopyIntercept(NoahmpType& noahmp) {
  int IndNewSnowLayer = 0;

  const real& MainTimeStep            = noahmp.config.domain.MainTimeStep;
  const real& SnowfallGround          = noahmp.water.flux.SnowfallGround;
  const real& SnowDepthIncr           = noahmp.water.flux.SnowDepthIncr;
  const real& TemperatureAirRefHeight = noahmp.forcing.TemperatureAirRefHeight;
  int&        NumSnowLayerNeg         = noahmp.config.domain.NumSnowLayerNeg;
  auto&       ThicknessSnowSoilLayer  = noahmp.config.domain.ThicknessSnowSoilLayer;
  real&       SnowDepth               = noahmp.water.state.SnowDepth;
  real&       SnowWaterEquiv          = noahmp.water.state.SnowWaterEquiv;
  auto&       SnowIce                 = noahmp.water.state.SnowIce;
  auto&       SnowLiqWater            = noahmp.water.state.SnowLiqWater;
  auto&       TemperatureSoilSnow     = noahmp.energy.state.TemperatureSoilSnow;

  if ((NumSnowLayerNeg == 0) && (SnowfallGround > real(0.0))) {
    SnowDepth      = SnowDepth + SnowDepthIncr * MainTimeStep;
    SnowWaterEquiv = SnowWaterEquiv + SnowfallGround * MainTimeStep;
  }

  if ((NumSnowLayerNeg == 0) && (SnowDepth >= real(0.025))) {
    NumSnowLayerNeg           = -1;
    IndNewSnowLayer           = 1;
    ThicknessSnowSoilLayer(0) = SnowDepth;
    SnowDepth                 = real(0.0);
    TemperatureSoilSnow(0)    = std::min(real(273.16), TemperatureAirRefHeight);
    SnowIce(0)                = SnowWaterEquiv;
    SnowLiqWater(0)           = real(0.0);
  }

  if ((NumSnowLayerNeg < 0) && (IndNewSnowLayer == 0) && (SnowfallGround > real(0.0))) {
    SnowIce(NumSnowLayerNeg + 1)                = SnowIce(NumSnowLayerNeg + 1) + SnowfallGround * MainTimeStep;
    ThicknessSnowSoilLayer(NumSnowLayerNeg + 1) = ThicknessSnowSoilLayer(NumSnowLayerNeg + 1) +
                                                  SnowDepthIncr * MainTimeStep;
  }
}

// ============================================================================
// SnowpackCompactionMod.F90
// ============================================================================
CPPWRF_HD inline void SnowpackCompaction(NoahmpType& noahmp) {
  real SnowBurden           = real(0.0);
  real SnowCompactAgeExpFac = real(0.0);
  real TempDiff             = real(0.0);
  real SnowVoid             = real(0.0);
  real SnowWatTotTmp        = real(0.0);
  real SnowIceDens          = real(0.0);

  const real& MainTimeStep           = noahmp.config.domain.MainTimeStep;
  const auto& TemperatureSoilSnow    = noahmp.energy.state.TemperatureSoilSnow;
  const auto& SnowIce                = noahmp.water.state.SnowIce;
  const auto& SnowLiqWater           = noahmp.water.state.SnowLiqWater;
  const auto& IndexPhaseChange       = noahmp.water.state.IndexPhaseChange;
  const auto& SnowIceFracPrev        = noahmp.water.state.SnowIceFracPrev;
  const real& SnowCompactBurdenFac   = noahmp.water.param.SnowCompactBurdenFac;
  const real& SnowCompactAgingFac1   = noahmp.water.param.SnowCompactAgingFac1;
  const real& SnowCompactAgingFac2   = noahmp.water.param.SnowCompactAgingFac2;
  const real& SnowCompactAgingFac3   = noahmp.water.param.SnowCompactAgingFac3;
  const real& SnowCompactAgingMax    = noahmp.water.param.SnowCompactAgingMax;
  const real& SnowViscosityCoeff     = noahmp.water.param.SnowViscosityCoeff;
  const int&  NumSnowLayerNeg        = noahmp.config.domain.NumSnowLayerNeg;
  auto&       ThicknessSnowSoilLayer = noahmp.config.domain.ThicknessSnowSoilLayer;
  auto&       CompactionSnowAging    = noahmp.water.flux.CompactionSnowAging;
  auto&       CompactionSnowBurden   = noahmp.water.flux.CompactionSnowBurden;
  auto&       CompactionSnowMelt     = noahmp.water.flux.CompactionSnowMelt;
  auto&       CompactionSnowTot      = noahmp.water.flux.CompactionSnowTot;
  auto&       SnowIceFrac            = noahmp.water.state.SnowIceFrac;

  CompactionSnowAging.fill(real(0.0));
  CompactionSnowBurden.fill(real(0.0));
  CompactionSnowMelt.fill(real(0.0));
  CompactionSnowTot.fill(real(0.0));
  SnowIceFrac.fill(real(0.0));

  SnowBurden = real(0.0);
  for (int LoopInd = NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
    SnowWatTotTmp        = SnowIce(LoopInd) + SnowLiqWater(LoopInd);
    SnowIceFrac(LoopInd) = SnowIce(LoopInd) / SnowWatTotTmp;
    SnowVoid             = real(1.0) - (SnowIce(LoopInd) / ConstDensityIce + SnowLiqWater(LoopInd) / ConstDensityWater) /
                                       ThicknessSnowSoilLayer(LoopInd);

    if ((SnowVoid > real(0.001)) && (SnowIce(LoopInd) > real(0.1))) {
      SnowIceDens = SnowIce(LoopInd) / ThicknessSnowSoilLayer(LoopInd);
      TempDiff    = std::max(real(0.0), ConstFreezePoint - TemperatureSoilSnow(LoopInd));

      SnowCompactAgeExpFac         = rexp(-SnowCompactAgingFac2 * TempDiff);
      CompactionSnowAging(LoopInd) = -SnowCompactAgingFac1 * SnowCompactAgeExpFac;
      if (SnowIceDens > SnowCompactAgingMax) {
        CompactionSnowAging(LoopInd) = CompactionSnowAging(LoopInd) *
                                       rexp(real(-46.0e-3) * (SnowIceDens - SnowCompactAgingMax));
      }
      if (SnowLiqWater(LoopInd) > (real(0.01) * ThicknessSnowSoilLayer(LoopInd))) {
        CompactionSnowAging(LoopInd) = CompactionSnowAging(LoopInd) * SnowCompactAgingFac3;
      }

      CompactionSnowBurden(LoopInd) = -(SnowBurden + real(0.5) * SnowWatTotTmp) *
                                      rexp(real(-0.08) * TempDiff - SnowCompactBurdenFac * SnowIceDens) / SnowViscosityCoeff;

      if (IndexPhaseChange(LoopInd) == 1) {
        CompactionSnowMelt(LoopInd) = std::max(real(0.0),
                                               (SnowIceFracPrev(LoopInd) - SnowIceFrac(LoopInd)) /
                                               std::max(real(1.0e-6), SnowIceFracPrev(LoopInd)));
        CompactionSnowMelt(LoopInd) = -CompactionSnowMelt(LoopInd) / MainTimeStep;
      } else {
        CompactionSnowMelt(LoopInd) = real(0.0);
      }

      CompactionSnowTot(LoopInd) = (CompactionSnowAging(LoopInd) + CompactionSnowBurden(LoopInd) +
                                    CompactionSnowMelt(LoopInd)) * MainTimeStep;
      CompactionSnowTot(LoopInd) = std::max(real(-0.5), CompactionSnowTot(LoopInd));

      ThicknessSnowSoilLayer(LoopInd) = ThicknessSnowSoilLayer(LoopInd) * (real(1.0) + CompactionSnowTot(LoopInd));
      ThicknessSnowSoilLayer(LoopInd) = std::max(ThicknessSnowSoilLayer(LoopInd),
                                                 SnowIce(LoopInd) / ConstDensityIce + SnowLiqWater(LoopInd) / ConstDensityWater);

      ThicknessSnowSoilLayer(LoopInd) = std::min(std::max(ThicknessSnowSoilLayer(LoopInd),
                                                          (SnowIce(LoopInd) + SnowLiqWater(LoopInd)) / real(500.0)),
                                                 (SnowIce(LoopInd) + SnowLiqWater(LoopInd)) / real(50.0));
    }

    SnowBurden = SnowBurden + SnowWatTotTmp;
  }
}

// ============================================================================
// SnowpackCompactionAR24Mod.F90
// ============================================================================
CPPWRF_HD inline void SnowpackCompactionAR24(NoahmpType& noahmp) {
  real SnowBurden           = real(0.0);
  real SnowCompactAgeExpFac = real(0.0);
  real TempDiff             = real(0.0);
  real SnowVoid             = real(0.0);
  real SnowWatTotTmp        = real(0.0);
  real SnowIceDens          = real(0.0);
  real SnowCompactBurdenFac = real(0.0);

  const real& MainTimeStep            = noahmp.config.domain.MainTimeStep;
  const real& TemperatureAirRefHeight = noahmp.forcing.TemperatureAirRefHeight;
  const real& PressureAirRefHeight    = noahmp.forcing.PressureAirRefHeight;
  const auto& TemperatureSoilSnow     = noahmp.energy.state.TemperatureSoilSnow;
  const auto& SnowIce                 = noahmp.water.state.SnowIce;
  const auto& SnowLiqWater            = noahmp.water.state.SnowLiqWater;
  const auto& IndexPhaseChange        = noahmp.water.state.IndexPhaseChange;
  const auto& SnowIceFracPrev         = noahmp.water.state.SnowIceFracPrev;
  const real& SnowCompactAgingFac1    = noahmp.water.param.SnowCompactAgingFac1;
  const real& SnowCompactAgingFac2    = noahmp.water.param.SnowCompactAgingFac2;
  const real& SnowCompactAgingFac3    = noahmp.water.param.SnowCompactAgingFac3;
  const real& SnowCompactm            = noahmp.water.param.SnowCompactmAR24;
  const real& SnowCompactb            = noahmp.water.param.SnowCompactbAR24;
  const real& SnowCompactPSFC1        = noahmp.water.param.SnowCompactP1AR24;
  const real& SnowCompactPSFC2        = noahmp.water.param.SnowCompactP2AR24;
  const real& SnowCompactPSFC3        = noahmp.water.param.SnowCompactP3AR24;
  const real& BurdenFacUpper          = noahmp.water.param.BurdenFacUpAR24;
  const real& SnowCompactAgingMax     = noahmp.water.param.SnowCompactAgingMax;
  const real& SnowViscosityCoeff      = noahmp.water.param.SnowViscosityCoeff;
  const int&  NumSnowLayerNeg         = noahmp.config.domain.NumSnowLayerNeg;
  auto&       ThicknessSnowSoilLayer  = noahmp.config.domain.ThicknessSnowSoilLayer;
  auto&       CompactionSnowAging     = noahmp.water.flux.CompactionSnowAging;
  auto&       CompactionSnowBurden    = noahmp.water.flux.CompactionSnowBurden;
  auto&       CompactionSnowMelt      = noahmp.water.flux.CompactionSnowMelt;
  auto&       CompactionSnowTot       = noahmp.water.flux.CompactionSnowTot;
  auto&       SnowIceFrac             = noahmp.water.state.SnowIceFrac;

  CompactionSnowAging.fill(real(0.0));
  CompactionSnowBurden.fill(real(0.0));
  CompactionSnowMelt.fill(real(0.0));
  CompactionSnowTot.fill(real(0.0));
  SnowIceFrac.fill(real(0.0));

  SnowBurden = real(0.0);
  SnowCompactBurdenFac = SnowCompactm * TemperatureAirRefHeight + SnowCompactb;
  if (PressureAirRefHeight >= real(85000.0)) {
    SnowCompactBurdenFac = std::max(SnowCompactBurdenFac, SnowCompactPSFC1);
  }
  if ((PressureAirRefHeight >= real(80000.0)) && (PressureAirRefHeight < real(85000.0))) {
    SnowCompactBurdenFac = std::max(SnowCompactBurdenFac, SnowCompactPSFC2);
  }
  if (PressureAirRefHeight < real(80000.0)) {
    SnowCompactBurdenFac = std::max(SnowCompactBurdenFac, SnowCompactPSFC3);
  }
  SnowCompactBurdenFac = std::min(SnowCompactBurdenFac, BurdenFacUpper);

  for (int LoopInd = NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
    SnowWatTotTmp        = SnowIce(LoopInd) + SnowLiqWater(LoopInd);
    SnowIceFrac(LoopInd) = SnowIce(LoopInd) / SnowWatTotTmp;
    SnowVoid             = real(1.0) - (SnowIce(LoopInd) / ConstDensityIce + SnowLiqWater(LoopInd) / ConstDensityWater) /
                                       ThicknessSnowSoilLayer(LoopInd);

    if ((SnowVoid > real(0.001)) && (SnowIce(LoopInd) > real(0.1))) {
      SnowIceDens = SnowIce(LoopInd) / ThicknessSnowSoilLayer(LoopInd);
      TempDiff    = std::max(real(0.0), ConstFreezePoint - TemperatureSoilSnow(LoopInd));

      SnowCompactAgeExpFac         = rexp(-SnowCompactAgingFac2 * TempDiff);
      CompactionSnowAging(LoopInd) = -SnowCompactAgingFac1 * SnowCompactAgeExpFac;
      if (SnowIceDens > SnowCompactAgingMax) {
        CompactionSnowAging(LoopInd) = CompactionSnowAging(LoopInd) *
                                       rexp(real(-46.0e-3) * (SnowIceDens - SnowCompactAgingMax));
      }
      if (SnowLiqWater(LoopInd) > (real(0.01) * ThicknessSnowSoilLayer(LoopInd))) {
        CompactionSnowAging(LoopInd) = CompactionSnowAging(LoopInd) * SnowCompactAgingFac3;
      }

      CompactionSnowBurden(LoopInd) = -(SnowBurden + real(0.5) * SnowWatTotTmp) *
                                      rexp(real(-0.08) * TempDiff - SnowCompactBurdenFac * SnowIceDens) / SnowViscosityCoeff;

      if (IndexPhaseChange(LoopInd) == 1) {
        CompactionSnowMelt(LoopInd) = std::max(real(0.0),
                                               (SnowIceFracPrev(LoopInd) - SnowIceFrac(LoopInd)) /
                                               std::max(real(1.0e-6), SnowIceFracPrev(LoopInd)));
        CompactionSnowMelt(LoopInd) = -CompactionSnowMelt(LoopInd) / MainTimeStep;
      } else {
        CompactionSnowMelt(LoopInd) = real(0.0);
      }

      CompactionSnowTot(LoopInd) = (CompactionSnowAging(LoopInd) + CompactionSnowBurden(LoopInd) +
                                    CompactionSnowMelt(LoopInd)) * MainTimeStep;
      CompactionSnowTot(LoopInd) = std::max(real(-0.5), CompactionSnowTot(LoopInd));

      ThicknessSnowSoilLayer(LoopInd) = ThicknessSnowSoilLayer(LoopInd) * (real(1.0) + CompactionSnowTot(LoopInd));
      ThicknessSnowSoilLayer(LoopInd) = std::max(ThicknessSnowSoilLayer(LoopInd),
                                                 SnowIce(LoopInd) / ConstDensityIce + SnowLiqWater(LoopInd) / ConstDensityWater);

      ThicknessSnowSoilLayer(LoopInd) = std::min(std::max(ThicknessSnowSoilLayer(LoopInd),
                                                          (SnowIce(LoopInd) + SnowLiqWater(LoopInd)) / real(500.0)),
                                                 (SnowIce(LoopInd) + SnowLiqWater(LoopInd)) / real(50.0));
    }

    SnowBurden = SnowBurden + SnowWatTotTmp;
  }
}

// ============================================================================
// SnowLayerCombineMod.F90
// ============================================================================
CPPWRF_HD inline void SnowLayerCombine(NoahmpType& noahmp) {
  int  NumSnowLayerOld = 0;
  int  IndLayer        = 0;
  int  IndNeighbor     = 0;
  real SnowIceTmp      = real(0.0);
  real SnowLiqTmp      = real(0.0);
  const real SnowThickMin[3] = {real(0.025), real(0.025), real(0.1)};

  int&       NumSnowLayerNeg        = noahmp.config.domain.NumSnowLayerNeg;
  auto&      ThicknessSnowSoilLayer = noahmp.config.domain.ThicknessSnowSoilLayer;
  auto&      TemperatureSoilSnow    = noahmp.energy.state.TemperatureSoilSnow;
  real&      SnowDepth              = noahmp.water.state.SnowDepth;
  real&      SnowWaterEquiv         = noahmp.water.state.SnowWaterEquiv;
  auto&      SnowIce                = noahmp.water.state.SnowIce;
  auto&      SnowLiqWater           = noahmp.water.state.SnowLiqWater;
  auto&      SoilLiqWater           = noahmp.water.state.SoilLiqWater;
  auto&      SoilIce                = noahmp.water.state.SoilIce;
  real&      PondSfcThinSnwComb     = noahmp.water.state.PondSfcThinSnwComb;
  real&      PondSfcThinSnwTrans    = noahmp.water.state.PondSfcThinSnwTrans;

  NumSnowLayerOld = NumSnowLayerNeg;

  for (int J = NumSnowLayerOld + 1; J <= 0; ++J) {
    if (SnowIce(J) <= real(0.1)) {
      if (J != 0) {
        SnowLiqWater(J + 1)           = SnowLiqWater(J + 1) + SnowLiqWater(J);
        SnowIce(J + 1)                = SnowIce(J + 1) + SnowIce(J);
        ThicknessSnowSoilLayer(J + 1) = ThicknessSnowSoilLayer(J + 1) + ThicknessSnowSoilLayer(J);
      } else {
        if (NumSnowLayerNeg < -1) {
          SnowLiqWater(J - 1)           = SnowLiqWater(J - 1) + SnowLiqWater(J);
          SnowIce(J - 1)                = SnowIce(J - 1) + SnowIce(J);
          ThicknessSnowSoilLayer(J - 1) = ThicknessSnowSoilLayer(J - 1) + ThicknessSnowSoilLayer(J);
        } else {
          if (SnowIce(J) >= real(0.0)) {
            PondSfcThinSnwComb = SnowLiqWater(J);
            SnowWaterEquiv     = SnowIce(J);
            SnowDepth          = ThicknessSnowSoilLayer(J);
          } else {
            PondSfcThinSnwComb = SnowLiqWater(J) + SnowIce(J);
            if (PondSfcThinSnwComb < real(0.0)) {
              SoilIce(1)         = SoilIce(1) + PondSfcThinSnwComb / (ThicknessSnowSoilLayer(1) * real(1000.0));
              PondSfcThinSnwComb = real(0.0);
            }
            SnowWaterEquiv = real(0.0);
            SnowDepth      = real(0.0);
          }
          SnowLiqWater(J)           = real(0.0);
          SnowIce(J)                = real(0.0);
          ThicknessSnowSoilLayer(J) = real(0.0);
        }
      }

      if ((J > NumSnowLayerNeg + 1) && (NumSnowLayerNeg < -1)) {
        for (int I = J; I >= NumSnowLayerNeg + 2; --I) {
          TemperatureSoilSnow(I)    = TemperatureSoilSnow(I - 1);
          SnowLiqWater(I)           = SnowLiqWater(I - 1);
          SnowIce(I)                = SnowIce(I - 1);
          ThicknessSnowSoilLayer(I) = ThicknessSnowSoilLayer(I - 1);
        }
      }
      NumSnowLayerNeg = NumSnowLayerNeg + 1;
    }
  }

  if (SoilIce(1) < real(0.0)) {
    SoilLiqWater(1) = SoilLiqWater(1) + SoilIce(1);
    SoilIce(1)      = real(0.0);
  }

  if (NumSnowLayerNeg == 0) return;

  SnowWaterEquiv = real(0.0);
  SnowDepth      = real(0.0);
  SnowIceTmp     = real(0.0);
  SnowLiqTmp     = real(0.0);

  for (int J = NumSnowLayerNeg + 1; J <= 0; ++J) {
    SnowWaterEquiv = SnowWaterEquiv + SnowIce(J) + SnowLiqWater(J);
    SnowDepth      = SnowDepth + ThicknessSnowSoilLayer(J);
    SnowIceTmp     = SnowIceTmp + SnowIce(J);
    SnowLiqTmp     = SnowLiqTmp + SnowLiqWater(J);
  }

  if ((SnowDepth < real(0.025)) && (NumSnowLayerNeg < 0)) {
    NumSnowLayerNeg     = 0;
    SnowWaterEquiv      = SnowIceTmp;
    PondSfcThinSnwTrans = SnowLiqTmp;
    if (SnowWaterEquiv <= real(0.0)) SnowDepth = real(0.0);
  }

  if (NumSnowLayerNeg < -1) {
    NumSnowLayerOld = NumSnowLayerNeg;
    IndLayer        = 1;
    for (int I = NumSnowLayerOld + 1; I <= 0; ++I) {
      if (ThicknessSnowSoilLayer(I) < SnowThickMin[IndLayer - 1]) {
        if (I == NumSnowLayerNeg + 1) {
          IndNeighbor = I + 1;
        } else if (I == 0) {
          IndNeighbor = I - 1;
        } else {
          IndNeighbor = I + 1;
          if ((ThicknessSnowSoilLayer(I - 1) + ThicknessSnowSoilLayer(I)) <
              (ThicknessSnowSoilLayer(I + 1) + ThicknessSnowSoilLayer(I))) {
            IndNeighbor = I - 1;
          }
        }
        int J = 0, L = 0;
        if (IndNeighbor > I) {
          J = IndNeighbor;
          L = I;
        } else {
          J = I;
          L = IndNeighbor;
        }

        SnowLayerWaterCombo(ThicknessSnowSoilLayer(J), SnowLiqWater(J), SnowIce(J), TemperatureSoilSnow(J),
                            ThicknessSnowSoilLayer(L), SnowLiqWater(L), SnowIce(L), TemperatureSoilSnow(L));

        if ((J - 1) > (NumSnowLayerNeg + 1)) {
          for (int K = J - 1; K >= NumSnowLayerNeg + 2; --K) {
            TemperatureSoilSnow(K)    = TemperatureSoilSnow(K - 1);
            SnowIce(K)                = SnowIce(K - 1);
            SnowLiqWater(K)           = SnowLiqWater(K - 1);
            ThicknessSnowSoilLayer(K) = ThicknessSnowSoilLayer(K - 1);
          }
        }
        NumSnowLayerNeg = NumSnowLayerNeg + 1;
        if (NumSnowLayerNeg >= -1) break;
      } else {
        IndLayer = IndLayer + 1;
      }
    }
  }
}

// ============================================================================
// SnowLayerDivideMod.F90
// ============================================================================
CPPWRF_HD inline void SnowLayerDivide(NoahmpType& noahmp) {
  int  NumSnowLayerTmp  = 0;
  real SnowThickCombTmp = real(0.0);
  real SnowIceExtra     = real(0.0);
  real SnowLiqExtra     = real(0.0);
  real SnowFracExtra    = real(0.0);
  real SnowTempGrad     = real(0.0);

  FArray1D<real, 1, kMaxSnowLayers> SnowThickTmp;
  FArray1D<real, 1, kMaxSnowLayers> SnowIceTmp;
  FArray1D<real, 1, kMaxSnowLayers> SnowLiqTmp;
  FArray1D<real, 1, kMaxSnowLayers> TemperatureSnowTmp;

  const int& NumSnowLayerMax        = noahmp.config.domain.NumSnowLayerMax;
  int&       NumSnowLayerNeg        = noahmp.config.domain.NumSnowLayerNeg;
  auto&      ThicknessSnowSoilLayer = noahmp.config.domain.ThicknessSnowSoilLayer;
  auto&      TemperatureSoilSnow    = noahmp.energy.state.TemperatureSoilSnow;
  auto&      SnowIce                = noahmp.water.state.SnowIce;
  auto&      SnowLiqWater           = noahmp.water.state.SnowLiqWater;

  SnowIceTmp.fill(real(0.0));
  SnowLiqTmp.fill(real(0.0));
  TemperatureSnowTmp.fill(real(0.0));
  SnowThickTmp.fill(real(0.0));

  for (int LoopInd = 1; LoopInd <= NumSnowLayerMax; ++LoopInd) {
    if (LoopInd <= std::abs(NumSnowLayerNeg)) {
      SnowThickTmp(LoopInd)       = ThicknessSnowSoilLayer(LoopInd + NumSnowLayerNeg);
      SnowIceTmp(LoopInd)         = SnowIce(LoopInd + NumSnowLayerNeg);
      SnowLiqTmp(LoopInd)         = SnowLiqWater(LoopInd + NumSnowLayerNeg);
      TemperatureSnowTmp(LoopInd) = TemperatureSoilSnow(LoopInd + NumSnowLayerNeg);
    }
  }

  NumSnowLayerTmp = std::abs(NumSnowLayerNeg);

  if (NumSnowLayerTmp == 1) {
    if (SnowThickTmp(1) > real(0.05)) {
      NumSnowLayerTmp       = 2;
      SnowThickTmp(1)       = SnowThickTmp(1) / real(2.0);
      SnowIceTmp(1)         = SnowIceTmp(1) / real(2.0);
      SnowLiqTmp(1)         = SnowLiqTmp(1) / real(2.0);
      SnowThickTmp(2)       = SnowThickTmp(1);
      SnowIceTmp(2)         = SnowIceTmp(1);
      SnowLiqTmp(2)         = SnowLiqTmp(1);
      TemperatureSnowTmp(2) = TemperatureSnowTmp(1);
    }
  }

  if (NumSnowLayerTmp > 1) {
    if (SnowThickTmp(1) > real(0.05)) {
      SnowThickCombTmp = SnowThickTmp(1) - real(0.05);
      SnowFracExtra    = SnowThickCombTmp / SnowThickTmp(1);
      SnowIceExtra     = SnowFracExtra * SnowIceTmp(1);
      SnowLiqExtra     = SnowFracExtra * SnowLiqTmp(1);

      SnowFracExtra   = real(0.05) / SnowThickTmp(1);
      SnowIceTmp(1)   = SnowFracExtra * SnowIceTmp(1);
      SnowLiqTmp(1)   = SnowFracExtra * SnowLiqTmp(1);
      SnowThickTmp(1) = real(0.05);

      SnowLayerWaterCombo(SnowThickTmp(2), SnowLiqTmp(2), SnowIceTmp(2), TemperatureSnowTmp(2),
                          SnowThickCombTmp, SnowLiqExtra, SnowIceExtra, TemperatureSnowTmp(1));

      if ((NumSnowLayerTmp <= 2) && (SnowThickTmp(2) > real(0.20))) {
        NumSnowLayerTmp       = 3;
        SnowTempGrad          = (TemperatureSnowTmp(1) - TemperatureSnowTmp(2)) /
                                ((SnowThickTmp(1) + SnowThickTmp(2)) / real(2.0));
        SnowThickTmp(2)       = SnowThickTmp(2) / real(2.0);
        SnowIceTmp(2)         = SnowIceTmp(2) / real(2.0);
        SnowLiqTmp(2)         = SnowLiqTmp(2) / real(2.0);
        SnowThickTmp(3)       = SnowThickTmp(2);
        SnowIceTmp(3)         = SnowIceTmp(2);
        SnowLiqTmp(3)         = SnowLiqTmp(2);
        TemperatureSnowTmp(3) = TemperatureSnowTmp(2) - SnowTempGrad * SnowThickTmp(2) / real(2.0);
        if (TemperatureSnowTmp(3) >= ConstFreezePoint) {
          TemperatureSnowTmp(3) = TemperatureSnowTmp(2);
        } else {
          TemperatureSnowTmp(2) = TemperatureSnowTmp(2) + SnowTempGrad * SnowThickTmp(2) / real(2.0);
        }
      }
    }
  }

  if (NumSnowLayerTmp > 2) {
    if (SnowThickTmp(2) > real(0.2)) {
      SnowThickCombTmp = SnowThickTmp(2) - real(0.2);
      SnowFracExtra    = SnowThickCombTmp / SnowThickTmp(2);
      SnowIceExtra     = SnowFracExtra * SnowIceTmp(2);
      SnowLiqExtra     = SnowFracExtra * SnowLiqTmp(2);

      SnowFracExtra   = real(0.2) / SnowThickTmp(2);
      SnowIceTmp(2)   = SnowFracExtra * SnowIceTmp(2);
      SnowLiqTmp(2)   = SnowFracExtra * SnowLiqTmp(2);
      SnowThickTmp(2) = real(0.2);

      SnowLayerWaterCombo(SnowThickTmp(3), SnowLiqTmp(3), SnowIceTmp(3), TemperatureSnowTmp(3),
                          SnowThickCombTmp, SnowLiqExtra, SnowIceExtra, TemperatureSnowTmp(2));
    }
  }

  NumSnowLayerNeg = -NumSnowLayerTmp;

  for (int LoopInd = NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
    ThicknessSnowSoilLayer(LoopInd) = SnowThickTmp(LoopInd - NumSnowLayerNeg);
    SnowIce(LoopInd)                = SnowIceTmp(LoopInd - NumSnowLayerNeg);
    SnowLiqWater(LoopInd)           = SnowLiqTmp(LoopInd - NumSnowLayerNeg);
    TemperatureSoilSnow(LoopInd)    = TemperatureSnowTmp(LoopInd - NumSnowLayerNeg);
  }
}

// ============================================================================
// SnowpackHydrologyMod.F90
// ============================================================================
CPPWRF_HD inline void SnowpackHydrology(NoahmpType& noahmp) {
  real InflowSnowLayer = real(0.0);
  real SnowIceTmp      = real(0.0);
  real SnowWaterRatio  = real(0.0);
  real SnowWaterTmp    = real(0.0);

  const real& MainTimeStep           = noahmp.config.domain.MainTimeStep;
  const real& FrostSnowSfcIce        = noahmp.water.flux.FrostSnowSfcIce;
  const real& SublimSnowSfcIce       = noahmp.water.flux.SublimSnowSfcIce;
  const real& RainfallGround         = noahmp.water.flux.RainfallGround;
  const real& SnowLiqFracMax         = noahmp.water.param.SnowLiqFracMax;
  const real& SnowLiqHoldCap         = noahmp.water.param.SnowLiqHoldCap;
  const real& SnowLiqReleaseFac      = noahmp.water.param.SnowLiqReleaseFac;
  int&        NumSnowLayerNeg        = noahmp.config.domain.NumSnowLayerNeg;
  auto&       ThicknessSnowSoilLayer = noahmp.config.domain.ThicknessSnowSoilLayer;
  real&       SnowDepth              = noahmp.water.state.SnowDepth;
  real&       SnowWaterEquiv         = noahmp.water.state.SnowWaterEquiv;
  auto&       SnowIce                = noahmp.water.state.SnowIce;
  auto&       SnowLiqWater           = noahmp.water.state.SnowLiqWater;
  auto&       SoilLiqWater           = noahmp.water.state.SoilLiqWater;
  auto&       SoilIce                = noahmp.water.state.SoilIce;
  auto&       SnowIceVol             = noahmp.water.state.SnowIceVol;
  auto&       SnowLiqWaterVol        = noahmp.water.state.SnowLiqWaterVol;
  auto&       SnowEffPorosity        = noahmp.water.state.SnowEffPorosity;
  real&       SnowBotOutflow         = noahmp.water.flux.SnowBotOutflow;
  auto&       OutflowSnowLayer       = noahmp.water.flux.OutflowSnowLayer;

  SnowEffPorosity.fill(real(0.0));
  SnowBotOutflow  = real(0.0);
  InflowSnowLayer = real(0.0);
  OutflowSnowLayer.fill(real(0.0));

  if (SnowWaterEquiv == real(0.0)) {
    SoilIce(1) = SoilIce(1) + (FrostSnowSfcIce - SublimSnowSfcIce) * MainTimeStep /
                              (ThicknessSnowSoilLayer(1) * real(1000.0));
    if (SoilIce(1) < real(0.0)) {
      SoilLiqWater(1) = SoilLiqWater(1) + SoilIce(1);
      SoilIce(1)      = real(0.0);
    }
  }

  if ((NumSnowLayerNeg == 0) && (SnowWaterEquiv > real(0.0))) {
    SnowWaterTmp   = SnowWaterEquiv;
    SnowWaterEquiv = SnowWaterEquiv - SublimSnowSfcIce * MainTimeStep + FrostSnowSfcIce * MainTimeStep;
    SnowWaterRatio = SnowWaterEquiv / SnowWaterTmp;
    SnowDepth      = std::max(real(0.0), SnowWaterRatio * SnowDepth);
    SnowDepth      = std::min(std::max(SnowDepth, SnowWaterEquiv / real(500.0)), SnowWaterEquiv / real(50.0));
    if (SnowWaterEquiv < real(0.0)) {
      SoilIce(1)     = SoilIce(1) + SnowWaterEquiv / (ThicknessSnowSoilLayer(1) * real(1000.0));
      SnowWaterEquiv = real(0.0);
      SnowDepth      = real(0.0);
    }
    if (SoilIce(1) < real(0.0)) {
      SoilLiqWater(1) = SoilLiqWater(1) + SoilIce(1);
      SoilIce(1)      = real(0.0);
    }
  }

  if ((SnowDepth <= real(1.0e-8)) || (SnowWaterEquiv <= real(1.0e-6))) {
    SnowDepth      = real(0.0);
    SnowWaterEquiv = real(0.0);
  }

  if (NumSnowLayerNeg < 0) {
    SnowWaterTmp               = SnowIce(NumSnowLayerNeg + 1) + SnowLiqWater(NumSnowLayerNeg + 1);
    SnowIceTmp                 = SnowIce(NumSnowLayerNeg + 1) - SublimSnowSfcIce * MainTimeStep + FrostSnowSfcIce * MainTimeStep;
    SnowIce(NumSnowLayerNeg + 1) = SnowIceTmp;
    if ((SnowIceTmp < real(1.0e-6)) && (NumSnowLayerNeg < 0)) SnowLayerCombine(noahmp);
    if ((SnowIceTmp >= real(1.0e-6)) && (NumSnowLayerNeg < 0)) {
      ThicknessSnowSoilLayer(NumSnowLayerNeg + 1) = ThicknessSnowSoilLayer(NumSnowLayerNeg + 1) *
          (SnowIce(NumSnowLayerNeg + 1) + SnowLiqWater(NumSnowLayerNeg + 1)) / SnowWaterTmp;
    }
    if (NumSnowLayerNeg < 0) {
      SnowLiqWater(NumSnowLayerNeg + 1) = SnowLiqWater(NumSnowLayerNeg + 1) + RainfallGround * MainTimeStep;
      SnowLiqWater(NumSnowLayerNeg + 1) = std::max(real(0.0), SnowLiqWater(NumSnowLayerNeg + 1));
    }
  }

  for (int LoopInd = NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
    SnowIceVol(LoopInd)      = std::min(real(1.0), SnowIce(LoopInd) / (ThicknessSnowSoilLayer(LoopInd) * ConstDensityIce));
    SnowEffPorosity(LoopInd) = real(1.0) - SnowIceVol(LoopInd);
  }

  for (int LoopInd = NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
    SnowLiqWater(LoopInd)     = SnowLiqWater(LoopInd) + InflowSnowLayer;
    SnowLiqWaterVol(LoopInd)  = SnowLiqWater(LoopInd) / (ThicknessSnowSoilLayer(LoopInd) * ConstDensityWater);
    OutflowSnowLayer(LoopInd) = std::max(real(0.0),
                                         (SnowLiqWaterVol(LoopInd) - SnowLiqHoldCap * SnowEffPorosity(LoopInd)) *
                                         ThicknessSnowSoilLayer(LoopInd));
    if (LoopInd == 0) {
      OutflowSnowLayer(LoopInd) = std::max((SnowLiqWaterVol(LoopInd) - SnowEffPorosity(LoopInd)) * ThicknessSnowSoilLayer(LoopInd),
                                           SnowLiqReleaseFac * MainTimeStep * OutflowSnowLayer(LoopInd));
    }
    OutflowSnowLayer(LoopInd) = OutflowSnowLayer(LoopInd) * ConstDensityWater;
    SnowLiqWater(LoopInd)     = SnowLiqWater(LoopInd) - OutflowSnowLayer(LoopInd);
    if ((SnowLiqWater(LoopInd) / (SnowIce(LoopInd) + SnowLiqWater(LoopInd))) > SnowLiqFracMax) {
      OutflowSnowLayer(LoopInd) = OutflowSnowLayer(LoopInd) +
                                  (SnowLiqWater(LoopInd) - SnowLiqFracMax / (real(1.0) - SnowLiqFracMax) * SnowIce(LoopInd));
      SnowLiqWater(LoopInd)     = SnowLiqFracMax / (real(1.0) - SnowLiqFracMax) * SnowIce(LoopInd);
    }
    InflowSnowLayer          = OutflowSnowLayer(LoopInd);
    SnowLiqWaterVol(LoopInd) = SnowLiqWater(LoopInd) / (ThicknessSnowSoilLayer(LoopInd) * ConstDensityWater);
  }

  for (int LoopInd = NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
    ThicknessSnowSoilLayer(LoopInd) = std::max(ThicknessSnowSoilLayer(LoopInd),
                                               SnowLiqWater(LoopInd) / ConstDensityWater + SnowIce(LoopInd) / ConstDensityIce);
  }

  SnowBotOutflow = OutflowSnowLayer(0) / MainTimeStep;
  for (int LoopInd = -noahmp.config.domain.NumSnowLayerMax + 1; LoopInd <= 0; ++LoopInd) {
    OutflowSnowLayer(LoopInd) = OutflowSnowLayer(LoopInd) / MainTimeStep;
  }
}

// ============================================================================
// SnowWaterMainMod.F90
// ============================================================================
CPPWRF_HD inline void SnowWaterMain(NoahmpType& noahmp) {
  real SnowDensBulk            = real(0.0);

  const int&  NumSnowLayerMax        = noahmp.config.domain.NumSnowLayerMax;
  const int&  NumSoilLayer           = noahmp.config.domain.NumSoilLayer;
  const real& MainTimeStep           = noahmp.config.domain.MainTimeStep;
  const auto& DepthSoilLayer         = noahmp.config.domain.DepthSoilLayer;
  const int&  OptSnowCompaction      = noahmp.config.nmlist.OptSnowCompaction;
  const real& SnoWatEqvMaxGlacier    = noahmp.water.param.SnoWatEqvMaxGlacier;
  auto&       ThicknessSnowSoilLayer = noahmp.config.domain.ThicknessSnowSoilLayer;
  auto&       DepthSnowSoilLayer     = noahmp.config.domain.DepthSnowSoilLayer;
  int&        NumSnowLayerNeg        = noahmp.config.domain.NumSnowLayerNeg;
  auto&       TemperatureSoilSnow    = noahmp.energy.state.TemperatureSoilSnow;
  real&       SnowDepth              = noahmp.water.state.SnowDepth;
  real&       SnowWaterEquiv         = noahmp.water.state.SnowWaterEquiv;
  auto&       SnowIce                = noahmp.water.state.SnowIce;
  auto&       SnowLiqWater           = noahmp.water.state.SnowLiqWater;
  real&       GlacierExcessFlow      = noahmp.water.flux.GlacierExcessFlow;
  real&       PondSfcThinSnwComb     = noahmp.water.state.PondSfcThinSnwComb;
  real&       PondSfcThinSnwTrans    = noahmp.water.state.PondSfcThinSnwTrans;

  GlacierExcessFlow       = real(0.0);
  PondSfcThinSnwComb      = real(0.0);
  PondSfcThinSnwTrans     = real(0.0);

  SnowfallAfterCanopyIntercept(noahmp);

  if (NumSnowLayerNeg < 0 && OptSnowCompaction == 1) SnowpackCompaction(noahmp);
  if (NumSnowLayerNeg < 0 && OptSnowCompaction == 2) SnowpackCompactionAR24(noahmp);

  if (NumSnowLayerNeg < 0) SnowLayerCombine(noahmp);

  if (NumSnowLayerNeg < 0) SnowLayerDivide(noahmp);

  SnowpackHydrology(noahmp);

  for (int LoopInd = -NumSnowLayerMax + 1; LoopInd <= NumSnowLayerNeg; ++LoopInd) {
    SnowIce(LoopInd)                = real(0.0);
    SnowLiqWater(LoopInd)           = real(0.0);
    TemperatureSoilSnow(LoopInd)    = real(0.0);
    ThicknessSnowSoilLayer(LoopInd) = real(0.0);
    DepthSnowSoilLayer(LoopInd)     = real(0.0);
  }

  if (SnowWaterEquiv > SnoWatEqvMaxGlacier) {
    SnowDensBulk              = SnowIce(0) / ThicknessSnowSoilLayer(0);
    GlacierExcessFlow         = SnowWaterEquiv - SnoWatEqvMaxGlacier;
    SnowIce(0)                = SnowIce(0) - GlacierExcessFlow;
    ThicknessSnowSoilLayer(0) = ThicknessSnowSoilLayer(0) - GlacierExcessFlow / SnowDensBulk;
    GlacierExcessFlow         = GlacierExcessFlow / MainTimeStep;
  }

  if (NumSnowLayerNeg < 0) {
    SnowWaterEquiv = real(0.0);
    for (int LoopInd = NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
      SnowWaterEquiv = SnowWaterEquiv + SnowIce(LoopInd) + SnowLiqWater(LoopInd);
    }
  }

  for (int LoopInd = NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
    ThicknessSnowSoilLayer(LoopInd) = -ThicknessSnowSoilLayer(LoopInd);
  }

  ThicknessSnowSoilLayer(1) = DepthSoilLayer(1);
  for (int LoopInd = 2; LoopInd <= NumSoilLayer; ++LoopInd) {
    ThicknessSnowSoilLayer(LoopInd) = DepthSoilLayer(LoopInd) - DepthSoilLayer(LoopInd - 1);
  }

  DepthSnowSoilLayer(NumSnowLayerNeg + 1) = ThicknessSnowSoilLayer(NumSnowLayerNeg + 1);
  for (int LoopInd = NumSnowLayerNeg + 2; LoopInd <= NumSoilLayer; ++LoopInd) {
    DepthSnowSoilLayer(LoopInd) = DepthSnowSoilLayer(LoopInd - 1) + ThicknessSnowSoilLayer(LoopInd);
  }

  for (int LoopInd = NumSnowLayerNeg + 1; LoopInd <= NumSoilLayer; ++LoopInd) {
    ThicknessSnowSoilLayer(LoopInd) = -ThicknessSnowSoilLayer(LoopInd);
  }

  if (NumSnowLayerNeg < 0) {
    SnowDepth = real(0.0);
    for (int LoopInd = NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
      SnowDepth = SnowDepth + ThicknessSnowSoilLayer(LoopInd);
    }
  }

  if ((SnowDepth <= real(1.0e-6)) || (SnowWaterEquiv <= real(1.0e-6))) {
    SnowDepth      = real(0.0);
    SnowWaterEquiv = real(0.0);
  }
}

// ============================================================================
// SoilHydraulicPropertyMod.F90
// ============================================================================
CPPWRF_HD inline void SoilDiffusivityConductivityOpt1(NoahmpType& noahmp,
                                                      real& SoilWatDiffusivity,
                                                      real& SoilWatConductivity,
                                                      real  SoilMoisture,
                                                      real  SoilImpervFrac,
                                                      int   IndLayer) {
  const auto& SoilMoistureSat        = noahmp.water.param.SoilMoistureSat;
  const auto& SoilExpCoeffB          = noahmp.water.param.SoilExpCoeffB;
  const auto& SoilWatDiffusivitySat  = noahmp.water.param.SoilWatDiffusivitySat;
  const auto& SoilWatConductivitySat = noahmp.water.param.SoilWatConductivitySat;

  real SoilPreFac = std::max(real(0.01), SoilMoisture / SoilMoistureSat(IndLayer));
  real SoilExpTmp = SoilExpCoeffB(IndLayer) + real(2.0);
  SoilWatDiffusivity = SoilWatDiffusivitySat(IndLayer) * rpow(SoilPreFac, SoilExpTmp);
  SoilWatDiffusivity = SoilWatDiffusivity * (real(1.0) - SoilImpervFrac);

  SoilExpTmp          = real(2.0) * SoilExpCoeffB(IndLayer) + real(3.0);
  SoilWatConductivity = SoilWatConductivitySat(IndLayer) * rpow(SoilPreFac, SoilExpTmp);
  SoilWatConductivity = SoilWatConductivity * (real(1.0) - SoilImpervFrac);
}

CPPWRF_HD inline void SoilDiffusivityConductivityOpt2(NoahmpType& noahmp,
                                                      real& SoilWatDiffusivity,
                                                      real& SoilWatConductivity,
                                                      real  SoilMoisture,
                                                      real  SoilIce,
                                                      int   IndLayer) {
  const auto& SoilMoistureSat        = noahmp.water.param.SoilMoistureSat;
  const auto& SoilExpCoeffB          = noahmp.water.param.SoilExpCoeffB;
  const auto& SoilWatDiffusivitySat  = noahmp.water.param.SoilWatDiffusivitySat;
  const auto& SoilWatConductivitySat = noahmp.water.param.SoilWatConductivitySat;

  real SoilPreFac1 = real(0.05) / SoilMoistureSat(IndLayer);
  real SoilPreFac2 = std::max(real(0.01), SoilMoisture / SoilMoistureSat(IndLayer));
  SoilPreFac1      = std::min(SoilPreFac1, SoilPreFac2);

  real SoilExpTmp = SoilExpCoeffB(IndLayer) + real(2.0);
  SoilWatDiffusivity = SoilWatDiffusivitySat(IndLayer) * rpow(SoilPreFac2, SoilExpTmp);
  if (SoilIce > real(0.0)) {
    real SoilIceWgt    = real(1.0) / (real(1.0) + rpow(real(500.0) * SoilIce, real(3.0)));
    SoilWatDiffusivity = SoilIceWgt * SoilWatDiffusivity +
                         (real(1.0) - SoilIceWgt) * SoilWatDiffusivitySat(IndLayer) * rpow(SoilPreFac1, SoilExpTmp);
  }

  SoilExpTmp          = real(2.0) * SoilExpCoeffB(IndLayer) + real(3.0);
  SoilWatConductivity = SoilWatConductivitySat(IndLayer) * rpow(SoilPreFac2, SoilExpTmp);
}

// ============================================================================
// RunoffSurfaceFreeDrainMod.F90 (OptRunoffSurface == 3, default)
// ============================================================================
CPPWRF_HD inline void RunoffSurfaceFreeDrain(NoahmpType& noahmp, real TimeStep) {
  constexpr int FrzSoilFac = 3;
  real FracVoidRem         = real(0.0);
  real SoilWatHoldMaxRem   = real(0.0);
  real WaterInSfc          = real(0.0);
  real TimeStepDay         = real(0.0);
  real SoilWatHoldMaxAcc   = real(0.0);
  real SoilIceWatTmp       = real(0.0);
  real SoilImpervFrac      = real(0.0);
  real IndAcc              = real(0.0);
  real SoilIceCoeff        = real(0.0);
  real SoilWatDiffusivity  = real(0.0);
  real SoilWatConductivity = real(0.0);
  real SoilWatHoldCap      = real(0.0);
  real InfilRateMax        = real(0.0);
  FArray1D<real, 1, kMaxSoilLayers> SoilWatMaxHold;

  const int&  NumSoilLayer        = noahmp.config.domain.NumSoilLayer;
  const auto& DepthSoilLayer      = noahmp.config.domain.DepthSoilLayer;
  const auto& SoilLiqWater        = noahmp.water.state.SoilLiqWater;
  const auto& SoilIce             = noahmp.water.state.SoilIce;
  const real& SoilIceMax          = noahmp.water.state.SoilIceMax;
  const real& SoilSfcInflowMean   = noahmp.water.flux.SoilSfcInflowMean;
  const auto& SoilMoistureSat     = noahmp.water.param.SoilMoistureSat;
  const auto& SoilMoistureWilt    = noahmp.water.param.SoilMoistureWilt;
  const real& SoilInfilMaxCoeff   = noahmp.water.param.SoilInfilMaxCoeff;
  const real& SoilImpervFracCoeff = noahmp.water.param.SoilImpervFracCoeff;
  real&       RunoffSurface       = noahmp.water.flux.RunoffSurface;
  real&       InfilRateSfc        = noahmp.water.flux.InfilRateSfc;

  SoilWatMaxHold.fill(real(0.0));

  if (SoilSfcInflowMean > real(0.0)) {
    TimeStepDay    = TimeStep / real(86400.0);
    SoilWatHoldCap = SoilMoistureSat(1) - SoilMoistureWilt(1);

    SoilWatMaxHold(1) = -DepthSoilLayer(1) * SoilWatHoldCap;
    SoilIceWatTmp     = -DepthSoilLayer(1) * SoilIce(1);
    SoilWatMaxHold(1) = SoilWatMaxHold(1) * (real(1.0) - (SoilLiqWater(1) + SoilIce(1) - SoilMoistureWilt(1)) / SoilWatHoldCap);
    SoilWatHoldMaxAcc = SoilWatMaxHold(1);
    for (int LoopInd3 = 2; LoopInd3 <= NumSoilLayer; ++LoopInd3) {
      SoilIceWatTmp            = SoilIceWatTmp + (DepthSoilLayer(LoopInd3 - 1) - DepthSoilLayer(LoopInd3)) * SoilIce(LoopInd3);
      SoilWatMaxHold(LoopInd3) = (DepthSoilLayer(LoopInd3 - 1) - DepthSoilLayer(LoopInd3)) * SoilWatHoldCap;
      SoilWatMaxHold(LoopInd3) = SoilWatMaxHold(LoopInd3) * (real(1.0) - (SoilLiqWater(LoopInd3) + SoilIce(LoopInd3) -
                                                             SoilMoistureWilt(LoopInd3)) / SoilWatHoldCap);
      SoilWatHoldMaxAcc        = SoilWatHoldMaxAcc + SoilWatMaxHold(LoopInd3);
    }
    FracVoidRem       = real(1.0) - rexp(real(-1.0) * SoilInfilMaxCoeff * TimeStepDay);
    SoilWatHoldMaxRem = SoilWatHoldMaxAcc * FracVoidRem;
    WaterInSfc        = std::max(real(0.0), SoilSfcInflowMean * TimeStep);
    InfilRateMax      = (WaterInSfc * (SoilWatHoldMaxRem / (WaterInSfc + SoilWatHoldMaxRem))) / TimeStep;

    SoilImpervFrac = real(1.0);
    if (SoilIceWatTmp > real(1.0e-2)) {
      SoilIceCoeff   = static_cast<real>(FrzSoilFac) * SoilImpervFracCoeff / SoilIceWatTmp;
      IndAcc         = real(1.0);
      int IndSoilFrz = FrzSoilFac - 1;
      for (int LoopInd1 = 1; LoopInd1 <= IndSoilFrz; ++LoopInd1) {
        int LoopInd3 = 1;
        for (int LoopInd2 = LoopInd1 + 1; LoopInd2 <= IndSoilFrz; ++LoopInd2) {
          LoopInd3 = LoopInd3 * LoopInd2;
        }
        IndAcc = IndAcc + cppwrf::powi(SoilIceCoeff, FrzSoilFac - LoopInd1) / static_cast<real>(static_cast<float>(LoopInd3));
      }
      SoilImpervFrac = real(1.0) - rexp(-SoilIceCoeff) * IndAcc;
    }

    InfilRateMax = InfilRateMax * SoilImpervFrac;

    SoilDiffusivityConductivityOpt2(noahmp, SoilWatDiffusivity, SoilWatConductivity, SoilLiqWater(1), SoilIceMax, 1);

    InfilRateMax = std::max(InfilRateMax, SoilWatConductivity);
    InfilRateMax = std::min(InfilRateMax, WaterInSfc / TimeStep);

    RunoffSurface = std::max(real(0.0), SoilSfcInflowMean - InfilRateMax);
    InfilRateSfc  = SoilSfcInflowMean - RunoffSurface;
  }
}

// ============================================================================
// RunoffSurfaceTopModelGrdMod.F90 (OptRunoffSurface == 1)
// ============================================================================
CPPWRF_HD inline void RunoffSurfaceTopModelGrd(NoahmpType& noahmp) {
  const real& SoilSfcInflowMean = noahmp.water.flux.SoilSfcInflowMean;
  real&       RunoffDecayFac    = noahmp.water.param.RunoffDecayFac;
  const real& SoilSfcSatFracMax = noahmp.water.param.SoilSfcSatFracMax;
  const auto& SoilExpCoeffB     = noahmp.water.param.SoilExpCoeffB;
  const auto& SoilImpervFrac    = noahmp.water.state.SoilImpervFrac;
  const real& WaterTableDepth   = noahmp.water.state.WaterTableDepth;
  real&       SoilSaturateFrac  = noahmp.water.state.SoilSaturateFrac;
  real&       RunoffSurface     = noahmp.water.flux.RunoffSurface;
  real&       InfilRateSfc      = noahmp.water.flux.InfilRateSfc;

  RunoffDecayFac   = SoilExpCoeffB(1) / real(3.0);
  SoilSaturateFrac = SoilSfcSatFracMax * rexp(real(-0.5) * RunoffDecayFac * WaterTableDepth);

  if (SoilSfcInflowMean > real(0.0)) {
    RunoffSurface = SoilSfcInflowMean * ((real(1.0) - SoilImpervFrac(1)) * SoilSaturateFrac + SoilImpervFrac(1));
    InfilRateSfc  = SoilSfcInflowMean - RunoffSurface;
  }
}

// ============================================================================
// RunoffSurfaceTopModelEquiMod.F90 (OptRunoffSurface == 2)
// ============================================================================
CPPWRF_HD inline void RunoffSurfaceTopModelEqui(NoahmpType& noahmp) {
  const real& SoilSfcInflowMean = noahmp.water.flux.SoilSfcInflowMean;
  real&       RunoffDecayFac    = noahmp.water.param.RunoffDecayFac;
  const real& SoilSfcSatFracMax = noahmp.water.param.SoilSfcSatFracMax;
  const auto& SoilImpervFrac    = noahmp.water.state.SoilImpervFrac;
  const real& WaterTableDepth   = noahmp.water.state.WaterTableDepth;
  real&       SoilSaturateFrac  = noahmp.water.state.SoilSaturateFrac;
  real&       RunoffSurface     = noahmp.water.flux.RunoffSurface;
  real&       InfilRateSfc      = noahmp.water.flux.InfilRateSfc;

  RunoffDecayFac   = real(2.0);
  SoilSaturateFrac = SoilSfcSatFracMax * rexp(real(-0.5) * RunoffDecayFac * WaterTableDepth);

  if (SoilSfcInflowMean > real(0.0)) {
    RunoffSurface = SoilSfcInflowMean * ((real(1.0) - SoilImpervFrac(1)) * SoilSaturateFrac + SoilImpervFrac(1));
    InfilRateSfc  = SoilSfcInflowMean - RunoffSurface;
  }
}

// ============================================================================
// RunoffSurfaceBatsMod.F90 (OptRunoffSurface == 4)
// ============================================================================
CPPWRF_HD inline void RunoffSurfaceBATS(NoahmpType& noahmp) {
  real SoilMoistureTmp = real(0.0);
  real SoilDepthTmp    = real(0.0);

  const int&  NumSoilLayer           = noahmp.config.domain.NumSoilLayer;
  const auto& ThicknessSnowSoilLayer = noahmp.config.domain.ThicknessSnowSoilLayer;
  const auto& SoilMoisture           = noahmp.water.state.SoilMoisture;
  const auto& SoilImpervFrac         = noahmp.water.state.SoilImpervFrac;
  const real& SoilSfcInflowMean      = noahmp.water.flux.SoilSfcInflowMean;
  const auto& SoilMoistureSat        = noahmp.water.param.SoilMoistureSat;
  real&       SoilSaturateFrac       = noahmp.water.state.SoilSaturateFrac;
  real&       RunoffSurface          = noahmp.water.flux.RunoffSurface;
  real&       InfilRateSfc           = noahmp.water.flux.InfilRateSfc;

  for (int LoopInd = 1; LoopInd <= NumSoilLayer; ++LoopInd) {
    SoilDepthTmp    = SoilDepthTmp + ThicknessSnowSoilLayer(LoopInd);
    SoilMoistureTmp = SoilMoistureTmp +
                      SoilMoisture(LoopInd) / SoilMoistureSat(LoopInd) * ThicknessSnowSoilLayer(LoopInd);
    if (SoilDepthTmp >= real(2.0)) break;
  }
  SoilMoistureTmp  = SoilMoistureTmp / SoilDepthTmp;
  SoilSaturateFrac = rpow(std::max(real(0.01), SoilMoistureTmp), real(4.0));

  if (SoilSfcInflowMean > real(0.0)) {
    RunoffSurface = SoilSfcInflowMean * ((real(1.0) - SoilImpervFrac(1)) * SoilSaturateFrac + SoilImpervFrac(1));
    InfilRateSfc  = SoilSfcInflowMean - RunoffSurface;
  }
}

// ============================================================================
// WaterTableEquilibriumMod.F90 (OptRunoffSubsurface == 2)
// ============================================================================
CPPWRF_HD inline void WaterTableEquilibrium(NoahmpType& noahmp) {
  constexpr int NumSoilFineLy = 100;
  real WatDeficitCoarse = real(0.0);
  real WatDeficitFine   = real(0.0);
  real ThickSoilFineLy  = real(0.0);
  real TmpVar           = real(0.0);
  FArray1D<real, 1, NumSoilFineLy> DepthSoilFineLy;
  DepthSoilFineLy.fill(real(0.0));

  const int&  NumSoilLayer           = noahmp.config.domain.NumSoilLayer;
  const auto& DepthSoilLayer         = noahmp.config.domain.DepthSoilLayer;
  const auto& ThicknessSnowSoilLayer = noahmp.config.domain.ThicknessSnowSoilLayer;
  const auto& SoilLiqWater           = noahmp.water.state.SoilLiqWater;
  const auto& SoilMoistureSat        = noahmp.water.param.SoilMoistureSat;
  const auto& SoilMatPotentialSat    = noahmp.water.param.SoilMatPotentialSat;
  const auto& SoilExpCoeffB          = noahmp.water.param.SoilExpCoeffB;
  real&       WaterTableDepth        = noahmp.water.state.WaterTableDepth;

  for (int IndSoil = 1; IndSoil <= NumSoilLayer; ++IndSoil) {
    WatDeficitCoarse = WatDeficitCoarse + (SoilMoistureSat(1) - SoilLiqWater(IndSoil)) *
                                          ThicknessSnowSoilLayer(IndSoil);
  }

  ThickSoilFineLy = real(3.0) * (-DepthSoilLayer(NumSoilLayer)) / static_cast<real>(NumSoilFineLy);
  for (int IndSoil = 1; IndSoil <= NumSoilFineLy; ++IndSoil) {
    DepthSoilFineLy(IndSoil) = static_cast<real>(static_cast<float>(IndSoil)) * ThickSoilFineLy;
  }

  WaterTableDepth = real(-3.0) * DepthSoilLayer(NumSoilLayer) - real(0.001);

  WatDeficitFine = real(0.0);
  for (int IndSoil = 1; IndSoil <= NumSoilFineLy; ++IndSoil) {
    TmpVar         = real(1.0) + (WaterTableDepth - DepthSoilFineLy(IndSoil)) / SoilMatPotentialSat(1);
    WatDeficitFine = WatDeficitFine + SoilMoistureSat(1) *
                                      (real(1.0) - rpow(TmpVar, real(-1.0) / SoilExpCoeffB(1))) * ThickSoilFineLy;
    if (std::abs(WatDeficitFine - WatDeficitCoarse) <= real(0.01)) {
      WaterTableDepth = DepthSoilFineLy(IndSoil);
      break;
    }
  }
}

// ============================================================================
// RunoffSubSurfaceEquiWaterTableMod.F90 (OptRunoffSubsurface == 2)
// ============================================================================
CPPWRF_HD inline void RunoffSubSurfaceEquiWaterTable(NoahmpType& noahmp) {
  const real& SoilImpervFracMax = noahmp.water.state.SoilImpervFracMax;
  const real& GridTopoIndex     = noahmp.water.param.GridTopoIndex;
  real&       RunoffDecayFac    = noahmp.water.param.RunoffDecayFac;
  real&       BaseflowCoeff     = noahmp.water.param.BaseflowCoeff;
  const real& WaterTableDepth   = noahmp.water.state.WaterTableDepth;
  real&       RunoffSubsurface  = noahmp.water.flux.RunoffSubsurface;

  RunoffDecayFac = real(2.0);
  BaseflowCoeff  = real(4.0);

  WaterTableEquilibrium(noahmp);

  RunoffSubsurface = (real(1.0) - SoilImpervFracMax) * BaseflowCoeff *
                     rexp(-GridTopoIndex) * rexp(-RunoffDecayFac * WaterTableDepth);
}

// ============================================================================
// GroundWaterTopModelMod.F90 & RunoffSubSurfaceGroundWaterMod.F90 (OptRunoffSubsurface == 1)
// ============================================================================
CPPWRF_HD inline void GroundWaterTopModel(NoahmpType& noahmp) {
  int  IndUnsatSoil      = 0;
  double SatDegUnsatSoil = 0.0;
  real SoilMatPotFrz     = real(0.0);
  real AquiferWatConduct = real(0.0);
  real WaterHeadTbl      = real(0.0);
  real WaterHead         = real(0.0);
  real WaterFillPore     = real(0.0);
  real WatConductAcc     = real(0.0);
  real SoilMoistureMin   = real(0.0);
  real WaterExcessSat    = real(0.0);
  FArray1D<real, 1, kMaxSoilLayers> ThicknessSoil;
  FArray1D<real, 1, kMaxSoilLayers> DepthSoilMid;
  FArray1D<real, 1, kMaxSoilLayers> SoilLiqTmp;
  FArray1D<real, 1, kMaxSoilLayers> SoilEffPorosity;
  FArray1D<real, 1, kMaxSoilLayers> SoilWatConductTmp;
  FArray1D<real, 1, kMaxSoilLayers> SoilMoisture;

  const int&  NumSoilLayer           = noahmp.config.domain.NumSoilLayer;
  const real& SoilTimeStep           = noahmp.config.domain.SoilTimeStep;
  const auto& DepthSoilLayer         = noahmp.config.domain.DepthSoilLayer;
  const real& SoilImpervFracMax      = noahmp.water.state.SoilImpervFracMax;
  const auto& SoilIce                = noahmp.water.state.SoilIce;
  const auto& SoilWatConductivity    = noahmp.water.state.SoilWatConductivity;
  const auto& SoilMoistureSat        = noahmp.water.param.SoilMoistureSat;
  const real& GridTopoIndex          = noahmp.water.param.GridTopoIndex;
  const auto& SoilMatPotentialSat    = noahmp.water.param.SoilMatPotentialSat;
  const auto& SoilExpCoeffB          = noahmp.water.param.SoilExpCoeffB;
  const real& SpecYieldGw            = noahmp.water.param.SpecYieldGw;
  const real& MicroPoreContent       = noahmp.water.param.MicroPoreContent;
  const auto& SoilWatConductivitySat = noahmp.water.param.SoilWatConductivitySat;
  auto&       SoilLiqWater           = noahmp.water.state.SoilLiqWater;
  real&       WaterTableDepth        = noahmp.water.state.WaterTableDepth;
  real&       WaterStorageAquifer    = noahmp.water.state.WaterStorageAquifer;
  real&       WaterStorageSoilAqf    = noahmp.water.state.WaterStorageSoilAqf;
  real&       RunoffDecayFac         = noahmp.water.param.RunoffDecayFac;
  real&       BaseflowCoeff          = noahmp.water.param.BaseflowCoeff;
  real&       RechargeGw             = noahmp.water.flux.RechargeGw;
  real&       DischargeGw            = noahmp.water.flux.DischargeGw;

  DepthSoilMid.fill(real(0.0));
  ThicknessSoil.fill(real(0.0));
  SoilLiqTmp.fill(real(0.0));
  SoilEffPorosity.fill(real(0.0));
  SoilWatConductTmp.fill(real(0.0));
  SoilMoisture.fill(real(0.0));
  DischargeGw = real(0.0);
  RechargeGw  = real(0.0);

  ThicknessSoil(1) = -DepthSoilLayer(1) * real(1.0e3);
  for (int LoopInd = 2; LoopInd <= NumSoilLayer; ++LoopInd) {
    ThicknessSoil(LoopInd) = real(1.0e3) * (DepthSoilLayer(LoopInd - 1) - DepthSoilLayer(LoopInd));
  }

  DepthSoilMid(1) = -DepthSoilLayer(1) / real(2.0);
  for (int LoopInd = 2; LoopInd <= NumSoilLayer; ++LoopInd) {
    DepthSoilMid(LoopInd) = -DepthSoilLayer(LoopInd - 1) +
                            real(0.5) * (DepthSoilLayer(LoopInd - 1) - DepthSoilLayer(LoopInd));
  }

  for (int LoopInd = 1; LoopInd <= NumSoilLayer; ++LoopInd) {
    SoilMoisture(LoopInd)      = SoilLiqWater(LoopInd) + SoilIce(LoopInd);
    SoilLiqTmp(LoopInd)        = SoilLiqWater(LoopInd) * ThicknessSoil(LoopInd);
    SoilEffPorosity(LoopInd)   = std::max(real(0.01), SoilMoistureSat(LoopInd) - SoilIce(LoopInd));
    SoilWatConductTmp(LoopInd) = real(1.0e3) * SoilWatConductivity(LoopInd);
  }

  IndUnsatSoil = NumSoilLayer;
  for (int LoopInd = 2; LoopInd <= NumSoilLayer; ++LoopInd) {
    if (WaterTableDepth <= -DepthSoilLayer(LoopInd)) {
      IndUnsatSoil = LoopInd - 1;
      break;
    }
  }

  RunoffDecayFac = SoilExpCoeffB(IndUnsatSoil) / real(3.0);
  BaseflowCoeff  = SoilWatConductTmp(IndUnsatSoil) * real(1.0e3) * rexp(real(3.0));
  DischargeGw    = (real(1.0) - SoilImpervFracMax) * BaseflowCoeff * rexp(-GridTopoIndex) *
                   rexp(-RunoffDecayFac * WaterTableDepth);

  SatDegUnsatSoil = static_cast<double>(std::min(real(1.0), SoilMoisture(IndUnsatSoil) / SoilMoistureSat(IndUnsatSoil)));
  SatDegUnsatSoil = std::max(SatDegUnsatSoil, 0.01);
  if (SatDegUnsatSoil < 0.01) SatDegUnsatSoil = 0.01;
  SoilMatPotFrz = -SoilMatPotentialSat(IndUnsatSoil) * real(1000.0) *
                  static_cast<real>(rpow(SatDegUnsatSoil, static_cast<double>(-SoilExpCoeffB(IndUnsatSoil))));
  SoilMatPotFrz = std::max(real(-120000.0), MicroPoreContent * SoilMatPotFrz);

  AquiferWatConduct = real(2.0) * (SoilWatConductTmp(IndUnsatSoil) * SoilWatConductivitySat(IndUnsatSoil) * real(1.0e3)) /
                      (SoilWatConductTmp(IndUnsatSoil) + SoilWatConductivitySat(IndUnsatSoil) * real(1.0e3));
  WaterHeadTbl      = -WaterTableDepth * real(1.0e3);
  WaterHead         = SoilMatPotFrz - DepthSoilMid(IndUnsatSoil) * real(1.0e3);
  RechargeGw        = -AquiferWatConduct * (WaterHeadTbl - WaterHead) /
                      ((WaterTableDepth - DepthSoilMid(IndUnsatSoil)) * real(1.0e3));
  RechargeGw        = std::max(real(-10.0) / SoilTimeStep, std::min(real(10.0) / SoilTimeStep, RechargeGw));

  WaterStorageSoilAqf = WaterStorageSoilAqf + (RechargeGw - DischargeGw) * SoilTimeStep;
  if (IndUnsatSoil == NumSoilLayer) {
    WaterStorageAquifer      = WaterStorageAquifer + (RechargeGw - DischargeGw) * SoilTimeStep;
    WaterStorageSoilAqf      = WaterStorageAquifer;
    WaterTableDepth          = (-DepthSoilLayer(NumSoilLayer) + real(25.0)) -
                               WaterStorageAquifer / real(1000.0) / SpecYieldGw;
    SoilLiqTmp(NumSoilLayer) = SoilLiqTmp(NumSoilLayer) - RechargeGw * SoilTimeStep;
    SoilLiqTmp(NumSoilLayer) = SoilLiqTmp(NumSoilLayer) + std::max(real(0.0), (WaterStorageAquifer - real(5000.0)));
    WaterStorageAquifer      = std::min(WaterStorageAquifer, real(5000.0));
  } else {
    if (IndUnsatSoil == NumSoilLayer - 1) {
      WaterTableDepth = -DepthSoilLayer(NumSoilLayer) -
                        (WaterStorageSoilAqf - SpecYieldGw * real(1000.0) * real(25.0)) /
                        (SoilEffPorosity(NumSoilLayer)) / real(1000.0);
    } else {
      WaterFillPore = real(0.0);
      for (int LoopInd = IndUnsatSoil + 2; LoopInd <= NumSoilLayer; ++LoopInd) {
        WaterFillPore = WaterFillPore + SoilEffPorosity(LoopInd) * ThicknessSoil(LoopInd);
      }
      WaterTableDepth = -DepthSoilLayer(IndUnsatSoil + 1) -
                        (WaterStorageSoilAqf - SpecYieldGw * real(1000.0) * real(25.0) - WaterFillPore) /
                        (SoilEffPorosity(IndUnsatSoil + 1)) / real(1000.0);
    }
    WatConductAcc = real(0.0);
    for (int LoopInd = 1; LoopInd <= NumSoilLayer; ++LoopInd) {
      WatConductAcc = WatConductAcc + SoilWatConductTmp(LoopInd) * ThicknessSoil(LoopInd);
    }
    for (int LoopInd = 1; LoopInd <= NumSoilLayer; ++LoopInd) {
      SoilLiqTmp(LoopInd) = SoilLiqTmp(LoopInd) - DischargeGw * SoilTimeStep *
                                                  SoilWatConductTmp(LoopInd) * ThicknessSoil(LoopInd) / WatConductAcc;
    }
  }
  WaterTableDepth = std::max(real(1.5), WaterTableDepth);

  SoilMoistureMin = real(0.01);
  for (int LoopInd = 1; LoopInd <= NumSoilLayer - 1; ++LoopInd) {
    if (SoilLiqTmp(LoopInd) < real(0.0)) {
      WaterExcessSat = SoilMoistureMin - SoilLiqTmp(LoopInd);
    } else {
      WaterExcessSat = real(0.0);
    }
    SoilLiqTmp(LoopInd)     = SoilLiqTmp(LoopInd) + WaterExcessSat;
    SoilLiqTmp(LoopInd + 1) = SoilLiqTmp(LoopInd + 1) - WaterExcessSat;
  }
  int LoopInd = NumSoilLayer;
  if (SoilLiqTmp(LoopInd) < SoilMoistureMin) {
    WaterExcessSat = SoilMoistureMin - SoilLiqTmp(LoopInd);
  } else {
    WaterExcessSat = real(0.0);
  }
  SoilLiqTmp(LoopInd) = SoilLiqTmp(LoopInd) + WaterExcessSat;
  WaterStorageAquifer = WaterStorageAquifer - WaterExcessSat;
  WaterStorageSoilAqf = WaterStorageSoilAqf - WaterExcessSat;

  for (int LoopInd_i = 1; LoopInd_i <= NumSoilLayer; ++LoopInd_i) {
    SoilLiqWater(LoopInd_i) = SoilLiqTmp(LoopInd_i) / ThicknessSoil(LoopInd_i);
  }
}

CPPWRF_HD inline void RunoffSubSurfaceGroundWater(NoahmpType& noahmp) {
  GroundWaterTopModel(noahmp);
  noahmp.water.flux.RunoffSubsurface = noahmp.water.flux.DischargeGw;
}

// ============================================================================
// RunoffSubSurfaceDrainageMod.F90
// ============================================================================
CPPWRF_HD inline void RunoffSubSurfaceDrainage(NoahmpType& noahmp) {
  noahmp.water.flux.RunoffSubsurface = noahmp.water.flux.RunoffSubsurface + noahmp.water.flux.DrainSoilBot;
}

// ============================================================================
// SoilWaterDiffusionRichardsMod.F90
// ============================================================================
template <typename Arr>
CPPWRF_HD inline void SoilWaterDiffusionRichards(NoahmpType& noahmp,
                                                 Arr& MatLeft1,
                                                 Arr& MatLeft2,
                                                 Arr& MatLeft3,
                                                 Arr& MatRight) {
  real DepthSnowSoilTmp = real(0.0);
  real SoilMoistTmpToWT = real(0.0);
  real SoilMoistBotTmp  = real(0.0);

  FArray1D<real, 1, kMaxSoilLayers> DepthSnowSoilInv;
  FArray1D<real, 1, kMaxSoilLayers> SoilThickTmp;
  FArray1D<real, 1, kMaxSoilLayers> SoilWaterGrad;
  FArray1D<real, 1, kMaxSoilLayers> WaterExcess;
  FArray1D<real, 1, kMaxSoilLayers> SoilMoistureTmp;

  const int&  NumSoilLayer              = noahmp.config.domain.NumSoilLayer;
  const auto& DepthSoilLayer            = noahmp.config.domain.DepthSoilLayer;
  const int&  OptSoilPermeabilityFrozen = noahmp.config.nmlist.OptSoilPermeabilityFrozen;
  const int&  OptRunoffSubsurface       = noahmp.config.nmlist.OptRunoffSubsurface;
  const real& SoilDrainSlope            = noahmp.water.param.SoilDrainSlope;
  const real& InfilRateSfc              = noahmp.water.flux.InfilRateSfc;
  const real& EvapSoilSfcLiqMean        = noahmp.water.flux.EvapSoilSfcLiqMean;
  const auto& TranspWatLossSoilMean     = noahmp.water.flux.TranspWatLossSoilMean;
  const auto& SoilLiqWater              = noahmp.water.state.SoilLiqWater;
  const auto& SoilMoisture              = noahmp.water.state.SoilMoisture;
  const real& WaterTableDepth           = noahmp.water.state.WaterTableDepth;
  const auto& SoilImpervFrac            = noahmp.water.state.SoilImpervFrac;
  const real& SoilImpervFracMax         = noahmp.water.state.SoilImpervFracMax;
  const real& SoilIceMax                = noahmp.water.state.SoilIceMax;
  const real& SoilMoistureToWT          = noahmp.water.state.SoilMoistureToWT;
  auto&       SoilWatConductivity       = noahmp.water.state.SoilWatConductivity;
  auto&       SoilWatDiffusivity        = noahmp.water.state.SoilWatDiffusivity;
  real&       DrainSoilBot              = noahmp.water.flux.DrainSoilBot;

  MatRight.fill(real(0.0));
  MatLeft1.fill(real(0.0));
  MatLeft2.fill(real(0.0));
  MatLeft3.fill(real(0.0));
  DepthSnowSoilInv.fill(real(0.0));
  SoilThickTmp.fill(real(0.0));
  SoilWaterGrad.fill(real(0.0));
  WaterExcess.fill(real(0.0));
  SoilMoistureTmp.fill(real(0.0));

  if (OptSoilPermeabilityFrozen == 1) {
    for (int LoopInd = 1; LoopInd <= NumSoilLayer; ++LoopInd) {
      SoilDiffusivityConductivityOpt1(noahmp, SoilWatDiffusivity(LoopInd), SoilWatConductivity(LoopInd),
                                      SoilMoisture(LoopInd), SoilImpervFrac(LoopInd), LoopInd);
      SoilMoistureTmp(LoopInd) = SoilMoisture(LoopInd);
    }
    if (OptRunoffSubsurface == 5) SoilMoistTmpToWT = SoilMoistureToWT;
  }

  if (OptSoilPermeabilityFrozen == 2) {
    for (int LoopInd = 1; LoopInd <= NumSoilLayer; ++LoopInd) {
      SoilDiffusivityConductivityOpt2(noahmp, SoilWatDiffusivity(LoopInd), SoilWatConductivity(LoopInd),
                                      SoilLiqWater(LoopInd), SoilIceMax, LoopInd);
      SoilMoistureTmp(LoopInd) = SoilLiqWater(LoopInd);
    }
    if (OptRunoffSubsurface == 5) {
      SoilMoistTmpToWT = SoilMoistureToWT * SoilLiqWater(NumSoilLayer) / SoilMoisture(NumSoilLayer);
    }
  }

  for (int LoopInd = 1; LoopInd <= NumSoilLayer; ++LoopInd) {
    if (LoopInd == 1) {
      SoilThickTmp(LoopInd)     = -DepthSoilLayer(LoopInd);
      DepthSnowSoilTmp          = -DepthSoilLayer(LoopInd + 1);
      DepthSnowSoilInv(LoopInd) = real(2.0) / DepthSnowSoilTmp;
      SoilWaterGrad(LoopInd)    = real(2.0) * (SoilMoistureTmp(LoopInd) - SoilMoistureTmp(LoopInd + 1)) / DepthSnowSoilTmp;
      WaterExcess(LoopInd)      = SoilWatDiffusivity(LoopInd) * SoilWaterGrad(LoopInd) + SoilWatConductivity(LoopInd) -
                                  InfilRateSfc + TranspWatLossSoilMean(LoopInd) + EvapSoilSfcLiqMean;
    } else if (LoopInd < NumSoilLayer) {
      SoilThickTmp(LoopInd)     = (DepthSoilLayer(LoopInd - 1) - DepthSoilLayer(LoopInd));
      DepthSnowSoilTmp          = (DepthSoilLayer(LoopInd - 1) - DepthSoilLayer(LoopInd + 1));
      DepthSnowSoilInv(LoopInd) = real(2.0) / DepthSnowSoilTmp;
      SoilWaterGrad(LoopInd)    = real(2.0) * (SoilMoistureTmp(LoopInd) - SoilMoistureTmp(LoopInd + 1)) / DepthSnowSoilTmp;
      WaterExcess(LoopInd)      = SoilWatDiffusivity(LoopInd) * SoilWaterGrad(LoopInd) + SoilWatConductivity(LoopInd) -
                                  SoilWatDiffusivity(LoopInd - 1) * SoilWaterGrad(LoopInd - 1) - SoilWatConductivity(LoopInd - 1) +
                                  TranspWatLossSoilMean(LoopInd);
    } else {
      SoilThickTmp(LoopInd) = (DepthSoilLayer(LoopInd - 1) - DepthSoilLayer(LoopInd));
      if ((OptRunoffSubsurface == 1) || (OptRunoffSubsurface == 2)) {
        DrainSoilBot = real(0.0);
      }
      if ((OptRunoffSubsurface == 3) || (OptRunoffSubsurface == 6) ||
          (OptRunoffSubsurface == 7) || (OptRunoffSubsurface == 8)) {
        DrainSoilBot = SoilDrainSlope * SoilWatConductivity(LoopInd);
      }
      if (OptRunoffSubsurface == 4) {
        DrainSoilBot = (real(1.0) - SoilImpervFracMax) * SoilWatConductivity(LoopInd);
      }
      if (OptRunoffSubsurface == 5) {
        DepthSnowSoilTmp = real(2.0) * SoilThickTmp(LoopInd);
        if (WaterTableDepth < (DepthSoilLayer(NumSoilLayer) - SoilThickTmp(NumSoilLayer))) {
          SoilMoistBotTmp = SoilMoistureTmp(LoopInd) - (SoilMoistureTmp(LoopInd) - SoilMoistTmpToWT) *
                            SoilThickTmp(LoopInd) * real(2.0) / (SoilThickTmp(LoopInd) + DepthSoilLayer(LoopInd) - WaterTableDepth);
        } else {
          SoilMoistBotTmp = SoilMoistTmpToWT;
        }
        SoilWaterGrad(LoopInd) = real(2.0) * (SoilMoistureTmp(LoopInd) - SoilMoistBotTmp) / DepthSnowSoilTmp;
        DrainSoilBot           = SoilWatDiffusivity(LoopInd) * SoilWaterGrad(LoopInd) + SoilWatConductivity(LoopInd);
      }
      WaterExcess(LoopInd) = -(SoilWatDiffusivity(LoopInd - 1) * SoilWaterGrad(LoopInd - 1)) - SoilWatConductivity(LoopInd - 1) +
                             TranspWatLossSoilMean(LoopInd) + DrainSoilBot;
    }
  }

  for (int LoopInd = 1; LoopInd <= NumSoilLayer; ++LoopInd) {
    if (LoopInd == 1) {
      MatLeft1(LoopInd) = real(0.0);
      MatLeft2(LoopInd) = SoilWatDiffusivity(LoopInd) * DepthSnowSoilInv(LoopInd) / SoilThickTmp(LoopInd);
      MatLeft3(LoopInd) = -MatLeft2(LoopInd);
    } else if (LoopInd < NumSoilLayer) {
      MatLeft1(LoopInd) = -SoilWatDiffusivity(LoopInd - 1) * DepthSnowSoilInv(LoopInd - 1) / SoilThickTmp(LoopInd);
      MatLeft3(LoopInd) = -SoilWatDiffusivity(LoopInd) * DepthSnowSoilInv(LoopInd) / SoilThickTmp(LoopInd);
      MatLeft2(LoopInd) = -(MatLeft1(LoopInd) + MatLeft3(LoopInd));
    } else {
      MatLeft1(LoopInd) = -SoilWatDiffusivity(LoopInd - 1) * DepthSnowSoilInv(LoopInd - 1) / SoilThickTmp(LoopInd);
      MatLeft3(LoopInd) = real(0.0);
      MatLeft2(LoopInd) = -(MatLeft1(LoopInd) + MatLeft3(LoopInd));
    }
    MatRight(LoopInd) = WaterExcess(LoopInd) / (-SoilThickTmp(LoopInd));
  }
}

// ============================================================================
// SoilMoistureSolverMod.F90
// ============================================================================
template <typename Arr>
CPPWRF_HD inline void SoilMoistureSolver(NoahmpType& noahmp,
                                         real TimeStep,
                                         Arr& MatLeft1,
                                         Arr& MatLeft2,
                                         Arr& MatLeft3,
                                         Arr& MatRight) {
  real WatDefiTmp = real(0.0);
  FArray1D<real, 1, kMaxSoilLayers> MatRightTmp;
  FArray1D<real, 1, kMaxSoilLayers> MatLeft3Tmp;

  const int&  NumSoilLayer           = noahmp.config.domain.NumSoilLayer;
  const auto& DepthSoilLayer         = noahmp.config.domain.DepthSoilLayer;
  const auto& ThicknessSnowSoilLayer = noahmp.config.domain.ThicknessSnowSoilLayer;
  const int&  OptRunoffSubsurface    = noahmp.config.nmlist.OptRunoffSubsurface;
  const auto& SoilMoistureSat        = noahmp.water.param.SoilMoistureSat;
  const real& WaterTableDepth        = noahmp.water.state.WaterTableDepth;
  const auto& SoilIce                = noahmp.water.state.SoilIce;
  auto&       SoilLiqWater           = noahmp.water.state.SoilLiqWater;
  auto&       SoilMoisture           = noahmp.water.state.SoilMoisture;
  real&       SoilMoistureToWT       = noahmp.water.state.SoilMoistureToWT;
  real&       RechargeGwDeepWT       = noahmp.water.state.RechargeGwDeepWT;
  real&       DrainSoilBot           = noahmp.water.flux.DrainSoilBot;
  auto&       SoilEffPorosity        = noahmp.water.state.SoilEffPorosity;
  real&       SoilSaturationExcess   = noahmp.water.state.SoilSaturationExcess;

  MatRightTmp.fill(real(0.0));
  MatLeft3Tmp.fill(real(0.0));
  SoilSaturationExcess = real(0.0);
  SoilEffPorosity.fill(real(0.0));

  for (int LoopInd = 1; LoopInd <= NumSoilLayer; ++LoopInd) {
    MatRight(LoopInd) = MatRight(LoopInd) * TimeStep;
    MatLeft1(LoopInd) = MatLeft1(LoopInd) * TimeStep;
    MatLeft2(LoopInd) = real(1.0) + MatLeft2(LoopInd) * TimeStep;
    MatLeft3(LoopInd) = MatLeft3(LoopInd) * TimeStep;
  }

  for (int LoopInd = 1; LoopInd <= NumSoilLayer; ++LoopInd) {
    MatRightTmp(LoopInd) = MatRight(LoopInd);
    MatLeft3Tmp(LoopInd) = MatLeft3(LoopInd);
  }

  MatrixSolverTriDiagonal(MatLeft3, MatLeft1, MatLeft2, MatLeft3Tmp, MatRightTmp, MatRight, 1, NumSoilLayer, 0);

  for (int LoopInd = 1; LoopInd <= NumSoilLayer; ++LoopInd) {
    SoilLiqWater(LoopInd) = SoilLiqWater(LoopInd) + MatLeft3(LoopInd);
  }

  if (OptRunoffSubsurface == 5) {
    if (WaterTableDepth < (DepthSoilLayer(NumSoilLayer) - ThicknessSnowSoilLayer(NumSoilLayer))) {
      RechargeGwDeepWT = RechargeGwDeepWT + TimeStep * DrainSoilBot;
    } else {
      SoilMoistureToWT           = SoilMoistureToWT + TimeStep * DrainSoilBot / ThicknessSnowSoilLayer(NumSoilLayer);
      SoilSaturationExcess       = std::max((SoilMoistureToWT - SoilMoistureSat(NumSoilLayer)), real(0.0)) *
                                   ThicknessSnowSoilLayer(NumSoilLayer);
      WatDefiTmp                 = std::max((real(1.0e-4) - SoilMoistureToWT), real(0.0)) * ThicknessSnowSoilLayer(NumSoilLayer);
      SoilMoistureToWT           = std::max(std::min(SoilMoistureToWT, SoilMoistureSat(NumSoilLayer)), real(1.0e-4));
      SoilLiqWater(NumSoilLayer) = SoilLiqWater(NumSoilLayer) +
                                   SoilSaturationExcess / ThicknessSnowSoilLayer(NumSoilLayer);
      DrainSoilBot               = DrainSoilBot - SoilSaturationExcess / TimeStep;
      RechargeGwDeepWT           = RechargeGwDeepWT - WatDefiTmp;
    }
  }

  for (int LoopInd = NumSoilLayer; LoopInd >= 2; --LoopInd) {
    SoilEffPorosity(LoopInd)  = std::max(real(1.0e-4), (SoilMoistureSat(LoopInd) - SoilIce(LoopInd)));
    SoilSaturationExcess      = std::max((SoilLiqWater(LoopInd) - SoilEffPorosity(LoopInd)), real(0.0)) *
                                ThicknessSnowSoilLayer(LoopInd);
    SoilLiqWater(LoopInd)     = std::min(SoilEffPorosity(LoopInd), SoilLiqWater(LoopInd));
    SoilLiqWater(LoopInd - 1) = SoilLiqWater(LoopInd - 1) + SoilSaturationExcess / ThicknessSnowSoilLayer(LoopInd - 1);
  }

  SoilEffPorosity(1)   = std::max(real(1.0e-4), (SoilMoistureSat(1) - SoilIce(1)));
  SoilSaturationExcess = std::max((SoilLiqWater(1) - SoilEffPorosity(1)), real(0.0)) * ThicknessSnowSoilLayer(1);
  SoilLiqWater(1)      = std::min(SoilEffPorosity(1), SoilLiqWater(1));

  if (SoilSaturationExcess > real(0.0)) {
    SoilLiqWater(2) = SoilLiqWater(2) + SoilSaturationExcess / ThicknessSnowSoilLayer(2);
    for (int LoopInd = 2; LoopInd <= NumSoilLayer - 1; ++LoopInd) {
      SoilEffPorosity(LoopInd)  = std::max(real(1.0e-4), (SoilMoistureSat(LoopInd) - SoilIce(LoopInd)));
      SoilSaturationExcess      = std::max((SoilLiqWater(LoopInd) - SoilEffPorosity(LoopInd)), real(0.0)) *
                                  ThicknessSnowSoilLayer(LoopInd);
      SoilLiqWater(LoopInd)     = std::min(SoilEffPorosity(LoopInd), SoilLiqWater(LoopInd));
      SoilLiqWater(LoopInd + 1) = SoilLiqWater(LoopInd + 1) + SoilSaturationExcess / ThicknessSnowSoilLayer(LoopInd + 1);
    }
    SoilEffPorosity(NumSoilLayer) = std::max(real(1.0e-4), (SoilMoistureSat(NumSoilLayer) - SoilIce(NumSoilLayer)));
    SoilSaturationExcess          = std::max((SoilLiqWater(NumSoilLayer) - SoilEffPorosity(NumSoilLayer)), real(0.0)) *
                                    ThicknessSnowSoilLayer(NumSoilLayer);
    SoilLiqWater(NumSoilLayer)    = std::min(SoilEffPorosity(NumSoilLayer), SoilLiqWater(NumSoilLayer));
  }

  for (int LoopInd = 1; LoopInd <= NumSoilLayer; ++LoopInd) {
    SoilMoisture(LoopInd) = SoilLiqWater(LoopInd) + SoilIce(LoopInd);
  }
}

// ============================================================================
// SoilWaterMainMod.F90
// ============================================================================
CPPWRF_HD inline void SoilWaterMain(NoahmpType& noahmp) {
  int  NumIterSoilWat    = 0;
  real TimeStepFine      = real(0.0);
  real SoilSatExcAcc     = real(0.0);
  real SoilWatConductAcc = real(0.0);
  real WaterRemove       = real(0.0);
  real SoilWatRem        = real(0.0);
  real SoilWaterMin      = real(0.0);
  real DrainSoilBotAcc   = real(0.0);
  real RunoffSurfaceAcc  = real(0.0);
  constexpr real SoilImpPara = real(4.0);

  FArray1D<real, 1, kMaxSoilLayers> MatRight;
  FArray1D<real, 1, kMaxSoilLayers> MatLeft1;
  FArray1D<real, 1, kMaxSoilLayers> MatLeft2;
  FArray1D<real, 1, kMaxSoilLayers> MatLeft3;
  FArray1D<real, 1, kMaxSoilLayers> SoilLiqTmp;

  const int&  NumSoilLayer           = noahmp.config.domain.NumSoilLayer;
  const real& SoilTimeStep           = noahmp.config.domain.SoilTimeStep;
  const auto& ThicknessSnowSoilLayer = noahmp.config.domain.ThicknessSnowSoilLayer;
  const bool& FlagUrban              = noahmp.config.domain.FlagUrban;
  const int&  OptRunoffSurface       = noahmp.config.nmlist.OptRunoffSurface;
  const int&  OptRunoffSubsurface    = noahmp.config.nmlist.OptRunoffSubsurface;
  const auto& SoilIce                = noahmp.water.state.SoilIce;
  const real& SoilSfcInflowMean      = noahmp.water.flux.SoilSfcInflowMean;
  const auto& SoilMoistureSat        = noahmp.water.param.SoilMoistureSat;
  auto&       SoilLiqWater           = noahmp.water.state.SoilLiqWater;
  auto&       SoilMoisture           = noahmp.water.state.SoilMoisture;
  real&       RechargeGwDeepWT       = noahmp.water.state.RechargeGwDeepWT;
  real&       DrainSoilBot           = noahmp.water.flux.DrainSoilBot;
  real&       RunoffSurface          = noahmp.water.flux.RunoffSurface;
  real&       RunoffSubsurface       = noahmp.water.flux.RunoffSubsurface;
  real&       InfilRateSfc           = noahmp.water.flux.InfilRateSfc;
  real&       TileDrain              = noahmp.water.flux.TileDrain;
  real&       SoilImpervFracMax      = noahmp.water.state.SoilImpervFracMax;
  const auto& SoilWatConductivity    = noahmp.water.state.SoilWatConductivity;
  auto&       SoilEffPorosity        = noahmp.water.state.SoilEffPorosity;
  auto&       SoilImpervFrac         = noahmp.water.state.SoilImpervFrac;
  auto&       SoilIceFrac            = noahmp.water.state.SoilIceFrac;
  const real& SoilSaturationExcess   = noahmp.water.state.SoilSaturationExcess;
  real&       SoilIceMax             = noahmp.water.state.SoilIceMax;
  real&       SoilLiqWaterMin        = noahmp.water.state.SoilLiqWaterMin;

  MatRight.fill(real(0.0));
  MatLeft1.fill(real(0.0));
  MatLeft2.fill(real(0.0));
  MatLeft3.fill(real(0.0));
  SoilLiqTmp.fill(real(0.0));
  RunoffSurface    = real(0.0);
  RunoffSubsurface = real(0.0);
  InfilRateSfc     = real(0.0);
  SoilSatExcAcc    = real(0.0);

  for (int LoopInd1 = 1; LoopInd1 <= NumSoilLayer; ++LoopInd1) {
    SoilEffPorosity(LoopInd1) = std::max(real(1.0e-4), (SoilMoistureSat(LoopInd1) - SoilIce(LoopInd1)));
    SoilSatExcAcc             = SoilSatExcAcc + std::max(real(0.0), SoilLiqWater(LoopInd1) - SoilEffPorosity(LoopInd1)) *
                                                ThicknessSnowSoilLayer(LoopInd1);
    SoilLiqWater(LoopInd1)    = std::min(SoilEffPorosity(LoopInd1), SoilLiqWater(LoopInd1));
  }

  for (int LoopInd1 = 1; LoopInd1 <= NumSoilLayer; ++LoopInd1) {
    SoilIceFrac(LoopInd1)    = std::min(real(1.0), SoilIce(LoopInd1) / SoilMoistureSat(LoopInd1));
    SoilImpervFrac(LoopInd1) = std::max(real(0.0),
                                        cppwrf::vm::exp(-SoilImpPara * (real(1.0) - SoilIceFrac(LoopInd1)),
                                                        LoopInd1 - 1, NumSoilLayer) -
                                            rexp(-SoilImpPara)) /
                               (real(1.0) - rexp(-SoilImpPara));
  }

  SoilIceMax        = real(0.0);
  SoilImpervFracMax = real(0.0);
  SoilLiqWaterMin   = SoilMoistureSat(1);
  for (int LoopInd1 = 1; LoopInd1 <= NumSoilLayer; ++LoopInd1) {
    if (SoilIce(LoopInd1) > SoilIceMax)               SoilIceMax        = SoilIce(LoopInd1);
    if (SoilImpervFrac(LoopInd1) > SoilImpervFracMax) SoilImpervFracMax = SoilImpervFrac(LoopInd1);
    if (SoilLiqWater(LoopInd1) < SoilLiqWaterMin)     SoilLiqWaterMin   = SoilLiqWater(LoopInd1);
  }

  if (OptRunoffSubsurface == 2) RunoffSubSurfaceEquiWaterTable(noahmp);

  if (FlagUrban) SoilImpervFrac(1) = real(0.95);

  if (OptRunoffSurface == 1) RunoffSurfaceTopModelGrd(noahmp);
  if (OptRunoffSurface == 2) RunoffSurfaceTopModelEqui(noahmp);
  if (OptRunoffSurface == 3) RunoffSurfaceFreeDrain(noahmp, SoilTimeStep);
  if (OptRunoffSurface == 4) RunoffSurfaceBATS(noahmp);

  NumIterSoilWat = 3;
  if ((InfilRateSfc * SoilTimeStep) > (ThicknessSnowSoilLayer(1) * SoilMoistureSat(1))) {
    NumIterSoilWat = NumIterSoilWat * 2;
  }
  TimeStepFine = SoilTimeStep / static_cast<real>(NumIterSoilWat);

  DrainSoilBotAcc  = real(0.0);
  RunoffSurfaceAcc = real(0.0);

  for (int IndIter = 1; IndIter <= NumIterSoilWat; ++IndIter) {
    if (SoilSfcInflowMean > real(0.0)) {
      if (OptRunoffSurface == 3) RunoffSurfaceFreeDrain(noahmp, TimeStepFine);
    }
    SoilWaterDiffusionRichards(noahmp, MatLeft1, MatLeft2, MatLeft3, MatRight);
    SoilMoistureSolver(noahmp, TimeStepFine, MatLeft1, MatLeft2, MatLeft3, MatRight);
    SoilSatExcAcc    = SoilSatExcAcc + SoilSaturationExcess;
    DrainSoilBotAcc  = DrainSoilBotAcc + DrainSoilBot;
    RunoffSurfaceAcc = RunoffSurfaceAcc + RunoffSurface;
  }

  DrainSoilBot  = DrainSoilBotAcc / static_cast<real>(NumIterSoilWat);
  RunoffSurface = RunoffSurfaceAcc / static_cast<real>(NumIterSoilWat);
  RunoffSurface = RunoffSurface * real(1000.0) + SoilSatExcAcc * real(1000.0) / SoilTimeStep;
  DrainSoilBot  = DrainSoilBot * real(1000.0);

  if (OptRunoffSubsurface == 2) {
    SoilWatConductAcc = real(0.0);
    for (int LoopInd1 = 1; LoopInd1 <= NumSoilLayer; ++LoopInd1) {
      SoilWatConductAcc = SoilWatConductAcc + SoilWatConductivity(LoopInd1) * ThicknessSnowSoilLayer(LoopInd1);
    }
    for (int LoopInd1 = 1; LoopInd1 <= NumSoilLayer; ++LoopInd1) {
      WaterRemove          = RunoffSubsurface * SoilTimeStep *
                             (SoilWatConductivity(LoopInd1) * ThicknessSnowSoilLayer(LoopInd1)) / SoilWatConductAcc;
      SoilLiqWater(LoopInd1) = SoilLiqWater(LoopInd1) - WaterRemove / (ThicknessSnowSoilLayer(LoopInd1) * real(1000.0));
    }
  }

  if (OptRunoffSubsurface != 1) {
    for (int LoopInd2 = 1; LoopInd2 <= NumSoilLayer; ++LoopInd2) {
      SoilLiqTmp(LoopInd2) = SoilLiqWater(LoopInd2) * ThicknessSnowSoilLayer(LoopInd2) * real(1000.0);
    }

    SoilWaterMin = real(0.01);
    for (int LoopInd2 = 1; LoopInd2 <= NumSoilLayer - 1; ++LoopInd2) {
      if (SoilLiqTmp(LoopInd2) < real(0.0)) {
        SoilWatRem = SoilWaterMin - SoilLiqTmp(LoopInd2);
      } else {
        SoilWatRem = real(0.0);
      }
      SoilLiqTmp(LoopInd2)     = SoilLiqTmp(LoopInd2) + SoilWatRem;
      SoilLiqTmp(LoopInd2 + 1) = SoilLiqTmp(LoopInd2 + 1) - SoilWatRem;
    }
    int LoopInd2 = NumSoilLayer;
    if (SoilLiqTmp(LoopInd2) < SoilWaterMin) {
      SoilWatRem = SoilWaterMin - SoilLiqTmp(LoopInd2);
    } else {
      SoilWatRem = real(0.0);
    }
    SoilLiqTmp(LoopInd2) = SoilLiqTmp(LoopInd2) + SoilWatRem;
    RunoffSubsurface     = RunoffSubsurface - SoilWatRem / SoilTimeStep;

    if (OptRunoffSubsurface == 5) RechargeGwDeepWT = RechargeGwDeepWT - SoilWatRem * real(1.0e-3);

    for (int LoopInd2_i = 1; LoopInd2_i <= NumSoilLayer; ++LoopInd2_i) {
      SoilLiqWater(LoopInd2_i) = SoilLiqTmp(LoopInd2_i) / (ThicknessSnowSoilLayer(LoopInd2_i) * real(1000.0));
    }
  }

  if (OptRunoffSubsurface == 1) RunoffSubSurfaceGroundWater(noahmp);

  if ((OptRunoffSubsurface == 3) || (OptRunoffSubsurface == 4) || (OptRunoffSubsurface == 6) ||
      (OptRunoffSubsurface == 7) || (OptRunoffSubsurface == 8)) {
    RunoffSubSurfaceDrainage(noahmp);
  }

  for (int LoopInd2 = 1; LoopInd2 <= NumSoilLayer; ++LoopInd2) {
    SoilMoisture(LoopInd2) = SoilLiqWater(LoopInd2) + SoilIce(LoopInd2);
  }

  RunoffSurface    = RunoffSurface * SoilTimeStep;
  RunoffSubsurface = RunoffSubsurface * SoilTimeStep;
  TileDrain        = TileDrain * SoilTimeStep;
}

// ============================================================================
// WaterMainMod.F90
// ============================================================================
CPPWRF_HD inline void WaterMain(NoahmpType& noahmp) {
  const real& MainTimeStep           = noahmp.config.domain.MainTimeStep;
  const real& SoilTimeStep           = noahmp.config.domain.SoilTimeStep;
  const int&  SurfaceType            = noahmp.config.domain.SurfaceType;
  const bool& FlagUrban              = noahmp.config.domain.FlagUrban;
  const bool& FlagSoilProcess        = noahmp.config.domain.FlagSoilProcess;
  const int&  NumSoilTimeStep        = noahmp.config.domain.NumSoilTimeStep;
  real&       VaporizeGrd            = noahmp.water.flux.VaporizeGrd;
  real&       CondenseVapGrd         = noahmp.water.flux.CondenseVapGrd;
  const real& RainfallGround         = noahmp.water.flux.RainfallGround;
  const auto& SoilTranspFac          = noahmp.water.state.SoilTranspFac;
  const real& WaterStorageLakeMax    = noahmp.water.param.WaterStorageLakeMax;
  const int&  NumSoilLayerRoot       = noahmp.water.param.NumSoilLayerRoot;
  const bool& FlagFrozenGround       = noahmp.energy.state.FlagFrozenGround;
  const real& LatHeatVapGrd          = noahmp.energy.state.LatHeatVapGrd;
  const real& DensityAirRefHeight    = noahmp.energy.state.DensityAirRefHeight;
  const real& ExchCoeffShSfc         = noahmp.energy.state.ExchCoeffShSfc;
  const real& SpecHumidityRefHeight  = noahmp.forcing.SpecHumidityRefHeight;
  const real& HeatLatentGrd          = noahmp.energy.flux.HeatLatentGrd;
  int&        NumSnowLayerNeg        = noahmp.config.domain.NumSnowLayerNeg;
  auto&       ThicknessSnowSoilLayer = noahmp.config.domain.ThicknessSnowSoilLayer;
  real&       SnowWaterEquiv         = noahmp.water.state.SnowWaterEquiv;
  real&       SnowWaterEquivPrev     = noahmp.water.state.SnowWaterEquivPrev;
  auto&       SoilLiqWater           = noahmp.water.state.SoilLiqWater;
  auto&       SoilIce                = noahmp.water.state.SoilIce;
  auto&       SoilMoisture           = noahmp.water.state.SoilMoisture;
  real&       WaterStorageLake       = noahmp.water.state.WaterStorageLake;
  real&       PondSfcThinSnwMelt     = noahmp.water.state.PondSfcThinSnwMelt;
  real&       SoilSfcInflow          = noahmp.water.flux.SoilSfcInflow;
  real&       EvapSoilSfcLiq         = noahmp.water.flux.EvapSoilSfcLiq;
  real&       DewSoilSfcLiq          = noahmp.water.flux.DewSoilSfcLiq;
  real&       FrostSnowSfcIce        = noahmp.water.flux.FrostSnowSfcIce;
  real&       SublimSnowSfcIce       = noahmp.water.flux.SublimSnowSfcIce;
  auto&       TranspWatLossSoil      = noahmp.water.flux.TranspWatLossSoil;
  real&       GlacierExcessFlow      = noahmp.water.flux.GlacierExcessFlow;
  real&       GlacierExcessFlowAcc   = noahmp.water.flux.GlacierExcessFlowAcc;
  real&       SoilSfcInflowAcc       = noahmp.water.flux.SoilSfcInflowAcc;
  real&       EvapSoilSfcLiqAcc      = noahmp.water.flux.EvapSoilSfcLiqAcc;
  auto&       TranspWatLossSoilAcc   = noahmp.water.flux.TranspWatLossSoilAcc;
  real&       SpecHumidity2mBare     = noahmp.energy.state.SpecHumidity2mBare;
  real&       SpecHumiditySfc        = noahmp.energy.state.SpecHumiditySfc;
  real&       EvapGroundNet          = noahmp.water.flux.EvapGroundNet;
  real&       Transpiration          = noahmp.water.flux.Transpiration;
  real&       EvapCanopyNet          = noahmp.water.flux.EvapCanopyNet;
  real&       RunoffSurface          = noahmp.water.flux.RunoffSurface;
  real&       RunoffSubsurface       = noahmp.water.flux.RunoffSubsurface;
  real&       TileDrain              = noahmp.water.flux.TileDrain;
  real&       SnowBotOutflow         = noahmp.water.flux.SnowBotOutflow;
  real&       WaterToAtmosTotal      = noahmp.water.flux.WaterToAtmosTotal;
  real&       SoilSfcInflowMean      = noahmp.water.flux.SoilSfcInflowMean;
  real&       EvapSoilSfcLiqMean     = noahmp.water.flux.EvapSoilSfcLiqMean;
  auto&       TranspWatLossSoilMean  = noahmp.water.flux.TranspWatLossSoilMean;
  real&       PondSfcThinSnwComb     = noahmp.water.state.PondSfcThinSnwComb;
  real&       PondSfcThinSnwTrans    = noahmp.water.state.PondSfcThinSnwTrans;

  TranspWatLossSoil.fill(real(0.0));
  SoilSfcInflow    = real(0.0);
  RunoffSurface    = real(0.0);
  RunoffSubsurface = real(0.0);
  TileDrain        = real(0.0);

  for (int LoopInd = 1; LoopInd <= noahmp.config.domain.NumSoilLayer; ++LoopInd) {
    SoilIce(LoopInd) = std::max(real(0.0), SoilMoisture(LoopInd) - SoilLiqWater(LoopInd));
  }
  SnowWaterEquivPrev = SnowWaterEquiv;
  VaporizeGrd        = std::max(HeatLatentGrd / LatHeatVapGrd, real(0.0));
  CondenseVapGrd     = std::abs(std::min(HeatLatentGrd / LatHeatVapGrd, real(0.0)));
  EvapGroundNet      = VaporizeGrd - CondenseVapGrd;

  CanopyHydrology(noahmp);

  SublimSnowSfcIce = real(0.0);
  if (SnowWaterEquiv > real(0.0)) {
    SublimSnowSfcIce = std::min(VaporizeGrd, SnowWaterEquiv / MainTimeStep);
  }
  EvapSoilSfcLiq = VaporizeGrd - SublimSnowSfcIce;

  FrostSnowSfcIce = real(0.0);
  if (SnowWaterEquiv > real(0.0)) {
    FrostSnowSfcIce = CondenseVapGrd;
  }
  DewSoilSfcLiq = CondenseVapGrd - FrostSnowSfcIce;

  SnowWaterMain(noahmp);

  GlacierExcessFlowAcc = GlacierExcessFlowAcc + GlacierExcessFlow * MainTimeStep;

  if (FlagFrozenGround) {
    SoilIce(1)     = SoilIce(1) + (DewSoilSfcLiq - EvapSoilSfcLiq) * MainTimeStep /
                                  (ThicknessSnowSoilLayer(1) * real(1000.0));
    DewSoilSfcLiq  = real(0.0);
    EvapSoilSfcLiq = real(0.0);
    if (SoilIce(1) < real(0.0)) {
      SoilLiqWater(1) = SoilLiqWater(1) + SoilIce(1);
      SoilIce(1)      = real(0.0);
    }
    SoilMoisture(1) = SoilLiqWater(1) + SoilIce(1);
  }
  EvapSoilSfcLiq = EvapSoilSfcLiq * real(0.001);

  for (int LoopInd = 1; LoopInd <= NumSoilLayerRoot; ++LoopInd) {
    TranspWatLossSoil(LoopInd) = Transpiration * SoilTranspFac(LoopInd) * real(0.001);
  }

  SoilSfcInflow = (PondSfcThinSnwMelt + PondSfcThinSnwComb + PondSfcThinSnwTrans) /
                  MainTimeStep * real(0.001);
  if (NumSnowLayerNeg == 0) {
    SoilSfcInflow = SoilSfcInflow + (SnowBotOutflow + DewSoilSfcLiq + RainfallGround) * real(0.001);
  } else {
    SoilSfcInflow = SoilSfcInflow + (SnowBotOutflow + DewSoilSfcLiq) * real(0.001);
  }

  SoilSfcInflowAcc  = SoilSfcInflowAcc + SoilSfcInflow;
  EvapSoilSfcLiqAcc = EvapSoilSfcLiqAcc + EvapSoilSfcLiq;
  for (int LoopInd = 1; LoopInd <= noahmp.config.domain.NumSoilLayer; ++LoopInd) {
    TranspWatLossSoilAcc(LoopInd) = TranspWatLossSoilAcc(LoopInd) + TranspWatLossSoil(LoopInd);
  }

  if (FlagSoilProcess) {
    SoilSfcInflowMean  = SoilSfcInflowAcc / static_cast<real>(NumSoilTimeStep);
    EvapSoilSfcLiqMean = EvapSoilSfcLiqAcc / static_cast<real>(NumSoilTimeStep);
    for (int LoopInd = 1; LoopInd <= noahmp.config.domain.NumSoilLayer; ++LoopInd) {
      TranspWatLossSoilMean(LoopInd) = TranspWatLossSoilAcc(LoopInd) / static_cast<real>(NumSoilTimeStep);
    }

    if (SurfaceType == 2) {
      RunoffSurface = real(0.0);
      if (WaterStorageLake >= WaterStorageLakeMax) {
        RunoffSurface = SoilSfcInflowMean * real(1000.0) * SoilTimeStep;
      }
      WaterStorageLake = WaterStorageLake + (SoilSfcInflowMean - EvapSoilSfcLiqMean) * real(1000.0) * SoilTimeStep - RunoffSurface;
    } else {
      SoilWaterMain(noahmp);
    }

    RunoffSubsurface = RunoffSubsurface + GlacierExcessFlowAcc;
  }

  WaterToAtmosTotal = Transpiration + EvapCanopyNet + EvapGroundNet;
  if (FlagUrban) {
    SpecHumiditySfc    = WaterToAtmosTotal / (DensityAirRefHeight * ExchCoeffShSfc) + SpecHumidityRefHeight;
    SpecHumidity2mBare = SpecHumiditySfc;
  }
}

// ============================================================================
// BalanceErrorCheckMod.F90
// ============================================================================
CPPWRF_HD inline void BalanceWaterInit(NoahmpType& noahmp) {
  const int&  NumSoilLayer           = noahmp.config.domain.NumSoilLayer;
  const int&  SurfaceType            = noahmp.config.domain.SurfaceType;
  const auto& ThicknessSnowSoilLayer = noahmp.config.domain.ThicknessSnowSoilLayer;
  const real& CanopyLiqWater         = noahmp.water.state.CanopyLiqWater;
  const real& CanopyIce              = noahmp.water.state.CanopyIce;
  const real& SnowWaterEquiv         = noahmp.water.state.SnowWaterEquiv;
  const auto& SoilMoisture           = noahmp.water.state.SoilMoisture;
  const real& WaterStorageAquifer    = noahmp.water.state.WaterStorageAquifer;
  const real& WaterStorageWetland    = noahmp.water.state.WaterStorageWetland;
  real&       WaterStorageTotBeg     = noahmp.water.state.WaterStorageTotBeg;

  if (SurfaceType == 1) {
    WaterStorageTotBeg = CanopyLiqWater + CanopyIce + SnowWaterEquiv + WaterStorageAquifer + WaterStorageWetland;
    for (int LoopInd = 1; LoopInd <= NumSoilLayer; ++LoopInd) {
      WaterStorageTotBeg = WaterStorageTotBeg + SoilMoisture(LoopInd) * ThicknessSnowSoilLayer(LoopInd) * real(1000.0);
    }
  }
}

CPPWRF_HD inline void BalanceWaterCheck(NoahmpType& noahmp) {
  const int&  NumSoilLayer            = noahmp.config.domain.NumSoilLayer;
  const int&  SurfaceType             = noahmp.config.domain.SurfaceType;
  const auto& ThicknessSnowSoilLayer  = noahmp.config.domain.ThicknessSnowSoilLayer;
  const real& MainTimeStep            = noahmp.config.domain.MainTimeStep;
  const bool& FlagCropland            = noahmp.config.domain.FlagCropland;
  const bool& FlagSoilProcess         = noahmp.config.domain.FlagSoilProcess;
  const real& IrriFracThreshold       = noahmp.water.param.IrriFracThreshold;
  const real& IrrigationFracGrid      = noahmp.water.state.IrrigationFracGrid;
  const real& CanopyLiqWater          = noahmp.water.state.CanopyLiqWater;
  const real& CanopyIce               = noahmp.water.state.CanopyIce;
  const real& SnowWaterEquiv          = noahmp.water.state.SnowWaterEquiv;
  const auto& SoilMoisture            = noahmp.water.state.SoilMoisture;
  const real& WaterStorageAquifer     = noahmp.water.state.WaterStorageAquifer;
  const real& WaterStorageWetland     = noahmp.water.state.WaterStorageWetland;
  const real& WaterStorageTotBeg      = noahmp.water.state.WaterStorageTotBeg;
  real&       PrecipTotRefHeight      = noahmp.water.flux.PrecipTotRefHeight;
  const real& EvapCanopyNet           = noahmp.water.flux.EvapCanopyNet;
  const real& Transpiration           = noahmp.water.flux.Transpiration;
  const real& EvapGroundNet           = noahmp.water.flux.EvapGroundNet;
  const real& RunoffSurface           = noahmp.water.flux.RunoffSurface;
  const real& RunoffSubsurface        = noahmp.water.flux.RunoffSubsurface;
  const real& TileDrain               = noahmp.water.flux.TileDrain;
  const real& IrrigationRateSprinkler = noahmp.water.flux.IrrigationRateSprinkler;
  const real& IrrigationRateMicro     = noahmp.water.flux.IrrigationRateMicro;
  const real& IrrigationRateFlood     = noahmp.water.flux.IrrigationRateFlood;
  real&       SfcWaterTotChgAcc       = noahmp.water.flux.SfcWaterTotChgAcc;
  real&       PrecipTotAcc            = noahmp.water.flux.PrecipTotAcc;
  real&       EvapCanopyNetAcc        = noahmp.water.flux.EvapCanopyNetAcc;
  real&       TranspirationAcc        = noahmp.water.flux.TranspirationAcc;
  real&       EvapGroundNetAcc        = noahmp.water.flux.EvapGroundNetAcc;
  real&       WaterStorageTotEnd      = noahmp.water.state.WaterStorageTotEnd;
  real&       WaterBalanceError       = noahmp.water.state.WaterBalanceError;

  if (FlagCropland && (IrrigationFracGrid >= IrriFracThreshold)) {
    PrecipTotRefHeight = PrecipTotRefHeight + IrrigationRateSprinkler * real(1000.0) / MainTimeStep;
  }

  if (SurfaceType == 1) {
    WaterStorageTotEnd = CanopyLiqWater + CanopyIce + SnowWaterEquiv + WaterStorageAquifer + WaterStorageWetland;
    for (int LoopInd = 1; LoopInd <= NumSoilLayer; ++LoopInd) {
      WaterStorageTotEnd = WaterStorageTotEnd + SoilMoisture(LoopInd) * ThicknessSnowSoilLayer(LoopInd) * real(1000.0);
    }
    SfcWaterTotChgAcc = SfcWaterTotChgAcc + (WaterStorageTotEnd - WaterStorageTotBeg);
    PrecipTotAcc      = PrecipTotAcc      + PrecipTotRefHeight * MainTimeStep;
    EvapCanopyNetAcc  = EvapCanopyNetAcc  + EvapCanopyNet      * MainTimeStep;
    TranspirationAcc  = TranspirationAcc  + Transpiration      * MainTimeStep;
    EvapGroundNetAcc  = EvapGroundNetAcc  + EvapGroundNet      * MainTimeStep;

    if (FlagSoilProcess) {
      WaterBalanceError = SfcWaterTotChgAcc -
                          (PrecipTotAcc + IrrigationRateMicro * real(1000.0) + IrrigationRateFlood * real(1000.0) -
                           EvapCanopyNetAcc - TranspirationAcc - EvapGroundNetAcc - RunoffSurface - RunoffSubsurface -
                           TileDrain);
      if (std::abs(WaterBalanceError) > real(0.1)) raise_fatal(kFatalNoahmpWater);
    }
  } else {
    WaterBalanceError = real(0.0);
  }
}

CPPWRF_HD inline void BalanceEnergyCheck(NoahmpType& noahmp) {
  const real& RadSwDownRefHeight   = noahmp.forcing.RadSwDownRefHeight;
  const real& RadSwAbsSfc          = noahmp.energy.flux.RadSwAbsSfc;
  const real& RadSwReflSfc         = noahmp.energy.flux.RadSwReflSfc;
  const real& RadLwNetSfc          = noahmp.energy.flux.RadLwNetSfc;
  const real& HeatSensibleSfc      = noahmp.energy.flux.HeatSensibleSfc;
  const real& HeatLatentCanopy     = noahmp.energy.flux.HeatLatentCanopy;
  const real& HeatLatentGrd        = noahmp.energy.flux.HeatLatentGrd;
  const real& HeatLatentTransp     = noahmp.energy.flux.HeatLatentTransp;
  const real& HeatGroundTot        = noahmp.energy.flux.HeatGroundTot;
  const real& RadSwAbsVeg          = noahmp.energy.flux.RadSwAbsVeg;
  const real& RadSwAbsGrd          = noahmp.energy.flux.RadSwAbsGrd;
  const real& HeatPrecipAdvSfc     = noahmp.energy.flux.HeatPrecipAdvSfc;
  const real& HeatLatentIrriEvap   = noahmp.energy.flux.HeatLatentIrriEvap;
  const real& HeatCanStorageChg    = noahmp.energy.flux.HeatCanStorageChg;
  real&       EnergyBalanceError   = noahmp.energy.state.EnergyBalanceError;
  real&       RadSwBalanceError    = noahmp.energy.state.RadSwBalanceError;

  RadSwBalanceError  = RadSwDownRefHeight - (RadSwAbsSfc + RadSwReflSfc);
  EnergyBalanceError = RadSwAbsVeg + RadSwAbsGrd + HeatPrecipAdvSfc -
                       (RadLwNetSfc + HeatSensibleSfc + HeatLatentCanopy + HeatLatentGrd +
                        HeatLatentTransp + HeatGroundTot + HeatLatentIrriEvap + HeatCanStorageChg);
  // BalanceErrorCheckMod: the run stops (the SNICAR check needs opt_alb = 3, not ported)
  if (std::abs(RadSwBalanceError) > real(0.01)) raise_fatal(kFatalNoahmpSolar);
  if (std::abs(EnergyBalanceError) > real(0.01)) raise_fatal(kFatalNoahmpEnergy);
}

// ============================================================================
// NoahmpMainMod.F90
// ============================================================================
CPPWRF_HD inline void NoahmpMain(NoahmpType& noahmp) {
  ProcessAtmosForcing(noahmp);
  GeneralInit(noahmp);
  BalanceWaterInit(noahmp);
  PhenologyMain(noahmp);
  IrrigationPrepare(noahmp);
  CanopyWaterIntercept(noahmp);
  PrecipitationHeatAdvect(noahmp);
  EnergyMain(noahmp);
  WaterMain(noahmp);
  BalanceWaterCheck(noahmp);
  BalanceEnergyCheck(noahmp);
}

}  // namespace cppwrf::phys::noahmp
