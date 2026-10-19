#pragma once

#include <algorithm>
#include <cmath>
#include "core/backend.hpp"
#include "core/precision.hpp"
#include "phys/lsm_noahmp/noahmp_rad_therm.hpp"
#include "phys/lsm_noahmp/noahmp_var.hpp"

namespace cppwrf::phys::noahmp {

// ============================================================================
// ResistanceCanopyStomataBallBerryMod.F90
// ============================================================================
CPPWRF_HD inline void ResistanceCanopyStomataBallBerry(NoahmpType& noahmp, int IndexShade) {
  auto& forcing = noahmp.forcing;
  auto& water = noahmp.water;
  auto& energy = noahmp.energy;
  auto& biochem = noahmp.biochem;

  constexpr int NumIter = 3;
  real MPE = real(1.0e-6);

  auto F1 = [](real AB, real BC) -> real {
    return rpow(AB, (BC - real(25.0)) / real(10.0));
  };
  auto F2 = [](real AB) -> real {
    return real(1.0) + rexp((real(-2.2e05) + real(710.0) * (AB + real(273.16))) /
                            (real(8.314) * (AB + real(273.16))));
  };

  real CF = forcing.PressureAirRefHeight / (real(8.314) * forcing.TemperatureAirRefHeight) * real(1.0e06);
  real ResistanceStomataTmp = real(1.0) / energy.param.ConductanceLeafMin * CF;
  real PhotosynLeafTmp = real(0.0);
  real RadPhotoActAbsTmp = real(0.0);
  if (IndexShade == 0) {
    RadPhotoActAbsTmp = energy.flux.RadPhotoActAbsSunlit / std::max(energy.state.VegFrac, real(1.0e-6));
  }
  if (IndexShade == 1) {
    RadPhotoActAbsTmp = energy.flux.RadPhotoActAbsShade / std::max(energy.state.VegFrac, real(1.0e-6));
  }

  if (RadPhotoActAbsTmp > real(0.0)) {
    real NitrogenFoliageFac =
        std::min(biochem.state.NitrogenConcFoliage / std::max(MPE, biochem.param.NitrogenConcFoliageMax),
                 real(1.0));
    real TC = energy.state.TemperatureCanopy - ConstFreezePoint;
    real PPF = real(4.6) * RadPhotoActAbsTmp;
    real J = PPF * biochem.param.QuantumEfficiency25C;
    real KC = energy.param.Co2MmConst25C * F1(energy.param.Co2MmConstQ10, TC);
    real KO = energy.param.O2MmConst25C * F1(energy.param.O2MmConstQ10, TC);
    real AWC = KC * (real(1.0) + energy.state.PressureAtmosO2 / KO);
    real CP = real(0.5) * KC / KO * energy.state.PressureAtmosO2 * real(0.21);
    real CarboxylRateMax = biochem.param.CarboxylRateMax25C / F2(TC) * NitrogenFoliageFac *
                           water.state.SoilTranspFacAcc * F1(biochem.param.CarboxylRateMaxQ10, TC);
    real CI = real(0.7) * energy.state.PressureAtmosCO2 * biochem.param.PhotosynPathC3 +
              real(0.4) * energy.state.PressureAtmosCO2 * (real(1.0) - biochem.param.PhotosynPathC3);
    real RLB = energy.state.ResistanceLeafBoundary / CF;
    real CEA = std::max(real(0.25) * energy.state.VapPresSatCanopy * biochem.param.PhotosynPathC3 +
                            real(0.40) * energy.state.VapPresSatCanopy * (real(1.0) - biochem.param.PhotosynPathC3),
                        std::min(energy.state.PressureVaporCanAir, energy.state.VapPresSatCanopy));

    for (int IndIter = 1; IndIter <= NumIter; ++IndIter) {
      real WJ = std::max(CI - CP, real(0.0)) * J / (CI + real(2.0) * CP) * biochem.param.PhotosynPathC3 +
                J * (real(1.0) - biochem.param.PhotosynPathC3);
      real WC = std::max(CI - CP, real(0.0)) * CarboxylRateMax / (CI + AWC) * biochem.param.PhotosynPathC3 +
                CarboxylRateMax * (real(1.0) - biochem.param.PhotosynPathC3);
      real WE = real(0.5) * CarboxylRateMax * biochem.param.PhotosynPathC3 +
                real(4000.0) * CarboxylRateMax * CI / forcing.PressureAirRefHeight *
                    (real(1.0) - biochem.param.PhotosynPathC3);
      PhotosynLeafTmp = std::min({WJ, WC, WE}) * biochem.state.IndexGrowSeason;
      real CS = std::max(energy.state.PressureAtmosCO2 -
                             real(1.37) * RLB * forcing.PressureAirRefHeight * PhotosynLeafTmp,
                         MPE);
      real A = biochem.param.SlopeConductToPhotosyn * PhotosynLeafTmp * forcing.PressureAirRefHeight * CEA /
                   (CS * energy.state.VapPresSatCanopy) +
               energy.param.ConductanceLeafMin;
      real B = (biochem.param.SlopeConductToPhotosyn * PhotosynLeafTmp * forcing.PressureAirRefHeight / CS +
                energy.param.ConductanceLeafMin) *
                   RLB -
               real(1.0);
      real C = -RLB;
      real Q = real(0.0);
      if (B >= real(0.0)) {
        Q = real(-0.5) * (B + std::sqrt(B * B - real(4.0) * A * C));
      } else {
        Q = real(-0.5) * (B - std::sqrt(B * B - real(4.0) * A * C));
      }
      real R1 = Q / A;
      real R2 = C / Q;
      ResistanceStomataTmp = std::max(R1, R2);
      CI = std::max(CS - PhotosynLeafTmp * forcing.PressureAirRefHeight * real(1.65) * ResistanceStomataTmp,
                    real(0.0));
    }

    ResistanceStomataTmp = ResistanceStomataTmp * CF;
  }

  if (IndexShade == 0) {
    energy.state.ResistanceStomataSunlit = ResistanceStomataTmp;
    biochem.flux.PhotosynLeafSunlit = PhotosynLeafTmp;
  }
  if (IndexShade == 1) {
    energy.state.ResistanceStomataShade = ResistanceStomataTmp;
    biochem.flux.PhotosynLeafShade = PhotosynLeafTmp;
  }
}

// ============================================================================
// ResistanceCanopyStomataJarvisMod.F90
// ============================================================================
CPPWRF_HD inline void ResistanceCanopyStomataJarvis(NoahmpType& noahmp, int IndexShade) {
  auto& forcing = noahmp.forcing;
  auto& water = noahmp.water;
  auto& energy = noahmp.energy;
  auto& biochem = noahmp.biochem;

  real RadPhotoActAbsTmp = real(0.0);
  if (IndexShade == 0) {
    RadPhotoActAbsTmp = energy.flux.RadPhotoActAbsSunlit / std::max(energy.state.VegFrac, real(1.0e-6));
  }
  if (IndexShade == 1) {
    RadPhotoActAbsTmp = energy.flux.RadPhotoActAbsShade / std::max(energy.state.VegFrac, real(1.0e-6));
  }

  real SpecHumidityTmp = real(0.622) * energy.state.PressureVaporCanAir /
                         (forcing.PressureAirRefHeight - real(0.378) * energy.state.PressureVaporCanAir);
  real MixingRatioTmp = SpecHumidityTmp / (real(1.0) - SpecHumidityTmp);
  real MixingRatioSat = real(0.0);
  real MixingRatioSatTempD = real(0.0);
  HumiditySaturation(energy.state.TemperatureCanopy, forcing.PressureAirRefHeight,
                     MixingRatioSat, MixingRatioSatTempD);

  real RadFac = real(2.0) * RadPhotoActAbsTmp / energy.param.RadiationStressFac;
  real ResistanceSolar = (RadFac + energy.param.ResistanceStomataMin / energy.param.ResistanceStomataMax) /
                         (real(1.0) + RadFac);
  ResistanceSolar = std::max(ResistanceSolar, real(0.0001));

  real ResistanceTemp =
      real(1.0) - real(0.0016) * rpow(energy.param.AirTempOptimTransp - energy.state.TemperatureCanopy, real(2.0));
  ResistanceTemp = std::max(ResistanceTemp, real(0.0001));

  real ResistanceVapDef =
      real(1.0) / (real(1.0) + energy.param.VaporPresDeficitFac * std::max(real(0.0), MixingRatioSat - MixingRatioTmp));
  ResistanceVapDef = std::max(ResistanceVapDef, real(0.01));

  real ResistanceStomataTmp =
      energy.param.ResistanceStomataMin /
      (ResistanceSolar * ResistanceTemp * ResistanceVapDef * water.state.SoilTranspFacAcc);
  real PhotosynLeafTmp = real(-999.99);

  if (IndexShade == 0) {
    energy.state.ResistanceStomataSunlit = ResistanceStomataTmp;
    biochem.flux.PhotosynLeafSunlit = PhotosynLeafTmp;
  }
  if (IndexShade == 1) {
    energy.state.ResistanceStomataShade = ResistanceStomataTmp;
    biochem.flux.PhotosynLeafShade = PhotosynLeafTmp;
  }
}

// ============================================================================
// ResistanceAboveCanopyMostMod.F90
// ============================================================================
CPPWRF_HD inline void ResistanceAboveCanopyMOST(NoahmpType& noahmp,
                                                int IterationInd,
                                                real HeatSensibleTmp,
                                                int& MoStabParaSgn) {
  auto& forcing = noahmp.forcing;
  auto& energy = noahmp.energy;

  real MPE = real(1.0e-6);
  real MOZOLD = energy.state.MoStabParaAbvCan;

  // gfortran SLP-vectorizes these 4 log calls into two 2-lane _ZGVbN2v_log calls
  real TMPCM  = cppwrf::vm::log((energy.state.RefHeightAboveGrd - energy.state.ZeroPlaneDispSfc) / energy.state.RoughLenMomSfc, 0, 2);
  real TMPCH  = cppwrf::vm::log((energy.state.RefHeightAboveGrd - energy.state.ZeroPlaneDispSfc) / energy.state.RoughLenShCanopy, 1, 2);
  real TMPCM2 = cppwrf::vm::log((real(2.0) + energy.state.RoughLenMomSfc) / energy.state.RoughLenMomSfc, 0, 2);
  real TMPCH2 = cppwrf::vm::log((real(2.0) + energy.state.RoughLenShCanopy) / energy.state.RoughLenShCanopy, 1, 2);

  if (IterationInd == 1) {
    energy.state.FrictionVelVeg = real(0.0);
    energy.state.MoStabParaAbvCan = real(0.0);
    energy.state.MoLengthAbvCan = real(0.0);
    energy.state.MoStabParaVeg2m = real(0.0);
  } else {
    real TVIR = (real(1.0) + real(0.61) * forcing.SpecHumidityRefHeight) * forcing.TemperatureAirRefHeight;
    real TMP1 = ConstVonKarman * (ConstGravityAcc / TVIR) * HeatSensibleTmp /
                (energy.state.DensityAirRefHeight * ConstHeatCapacAir);
    if (std::abs(TMP1) <= MPE) TMP1 = MPE;
    energy.state.MoLengthAbvCan = real(-1.0) * cppwrf::powi(energy.state.FrictionVelVeg, 3) / TMP1;
    energy.state.MoStabParaAbvCan =
        std::min((energy.state.RefHeightAboveGrd - energy.state.ZeroPlaneDispSfc) / energy.state.MoLengthAbvCan,
                 real(1.0));
    energy.state.MoStabParaVeg2m =
        std::min((real(2.0) + energy.state.RoughLenShCanopy) / energy.state.MoLengthAbvCan, real(1.0));
  }

  if (MOZOLD * energy.state.MoStabParaAbvCan < real(0.0)) MoStabParaSgn = MoStabParaSgn + 1;
  if (MoStabParaSgn >= 2) {
    energy.state.MoStabParaAbvCan = real(0.0);
    energy.state.MoStabCorrMomAbvCan = real(0.0);
    energy.state.MoStabCorrShAbvCan = real(0.0);
    energy.state.MoStabParaVeg2m = real(0.0);
    energy.state.MoStabCorrMomVeg2m = real(0.0);
    energy.state.MoStabCorrShVeg2m = real(0.0);
  }

  real FMNEW = real(0.0);
  real FHNEW = real(0.0);
  real FM2NEW = real(0.0);
  real FH2NEW = real(0.0);
  if (energy.state.MoStabParaAbvCan < real(0.0)) {
    real TMP1 = rpow(real(1.0) - real(16.0) * energy.state.MoStabParaAbvCan, real(0.25));
    real TMP2 = rlog((real(1.0) + TMP1 * TMP1) / real(2.0));
    real TMP3 = rlog((real(1.0) + TMP1) / real(2.0));
    FMNEW = real(2.0) * TMP3 + TMP2 - real(2.0) * std::atan(TMP1) + real(1.5707963);
    FHNEW = real(2.0) * TMP2;

    real TMP12 = rpow(real(1.0) - real(16.0) * energy.state.MoStabParaVeg2m, real(0.25));
    real TMP22 = rlog((real(1.0) + TMP12 * TMP12) / real(2.0));
    real TMP32 = rlog((real(1.0) + TMP12) / real(2.0));
    FM2NEW = real(2.0) * TMP32 + TMP22 - real(2.0) * std::atan(TMP12) + real(1.5707963);
    FH2NEW = real(2.0) * TMP22;
  } else {
    FMNEW = real(-5.0) * energy.state.MoStabParaAbvCan;
    FHNEW = FMNEW;
    FM2NEW = real(-5.0) * energy.state.MoStabParaVeg2m;
    FH2NEW = FM2NEW;
  }

  if (IterationInd == 1) {
    energy.state.MoStabCorrMomAbvCan = FMNEW;
    energy.state.MoStabCorrShAbvCan = FHNEW;
    energy.state.MoStabCorrMomVeg2m = FM2NEW;
    energy.state.MoStabCorrShVeg2m = FH2NEW;
  } else {
    energy.state.MoStabCorrMomAbvCan = real(0.5) * (energy.state.MoStabCorrMomAbvCan + FMNEW);
    energy.state.MoStabCorrShAbvCan = real(0.5) * (energy.state.MoStabCorrShAbvCan + FHNEW);
    energy.state.MoStabCorrMomVeg2m = real(0.5) * (energy.state.MoStabCorrMomVeg2m + FM2NEW);
    energy.state.MoStabCorrShVeg2m = real(0.5) * (energy.state.MoStabCorrShVeg2m + FH2NEW);
  }

  energy.state.MoStabCorrShAbvCan = std::min(energy.state.MoStabCorrShAbvCan, real(0.9) * TMPCH);
  energy.state.MoStabCorrMomAbvCan = std::min(energy.state.MoStabCorrMomAbvCan, real(0.9) * TMPCM);
  energy.state.MoStabCorrShVeg2m = std::min(energy.state.MoStabCorrShVeg2m, real(0.9) * TMPCH2);
  energy.state.MoStabCorrMomVeg2m = std::min(energy.state.MoStabCorrMomVeg2m, real(0.9) * TMPCM2);
  real CMFM = TMPCM - energy.state.MoStabCorrMomAbvCan;
  real CHFH = TMPCH - energy.state.MoStabCorrShAbvCan;
  real CM2FM2 = TMPCM2 - energy.state.MoStabCorrMomVeg2m;
  real CH2FH2 = TMPCH2 - energy.state.MoStabCorrShVeg2m;
  if (std::abs(CMFM) <= MPE) CMFM = MPE;
  if (std::abs(CHFH) <= MPE) CHFH = MPE;
  if (std::abs(CM2FM2) <= MPE) CM2FM2 = MPE;
  if (std::abs(CH2FH2) <= MPE) CH2FH2 = MPE;
  energy.state.ExchCoeffMomAbvCan = ConstVonKarman * ConstVonKarman / (CMFM * CMFM);
  energy.state.ExchCoeffShAbvCan = ConstVonKarman * ConstVonKarman / (CMFM * CHFH);

  energy.state.FrictionVelVeg = energy.state.WindSpdRefHeight * std::sqrt(energy.state.ExchCoeffMomAbvCan);
  energy.state.ExchCoeffSh2mVegMo = ConstVonKarman * energy.state.FrictionVelVeg / CH2FH2;

  energy.state.ResistanceMomAbvCan =
      std::max(real(1.0), real(1.0) / (energy.state.ExchCoeffMomAbvCan * energy.state.WindSpdRefHeight));
  energy.state.ResistanceShAbvCan =
      std::max(real(1.0), real(1.0) / (energy.state.ExchCoeffShAbvCan * energy.state.WindSpdRefHeight));
  energy.state.ResistanceLhAbvCan = energy.state.ResistanceShAbvCan;
}

// ============================================================================
// ResistanceAboveCanopyChen97Mod.F90
// ============================================================================
CPPWRF_HD inline void ResistanceAboveCanopyChen97(NoahmpType& noahmp, int IterationInd) {
  auto& energy = noahmp.energy;

  constexpr real WWST   = real(1.2);
  constexpr real WWST2  = WWST * WWST;
  constexpr real VKRM   = real(0.40);
  constexpr real EXCM   = real(0.001);
  constexpr real BETA   = real(1.0) / real(270.0);
  constexpr real BTG    = BETA * ConstGravityAcc;
  constexpr real ELFC   = VKRM * BTG;
  constexpr real WOLD   = real(0.15);
  constexpr real WNEW   = real(1.0) - WOLD;
  constexpr real PIHF   = real(3.14159265) / real(2.0);
  constexpr real EPSU2  = real(1.0e-4);
  constexpr real EPSUST = real(0.07);
  constexpr real ZTMIN  = real(-5.0);
  constexpr real ZTMAX  = real(1.0);
  constexpr real HPBL   = real(1000.0);
  constexpr real SQVISC = real(258.2);
  constexpr real RIC    = real(0.183);
  constexpr real FHNEU  = real(0.8);
  constexpr real RFC    = real(0.191);
  constexpr real RFAC   = RIC / (FHNEU * RFC * RFC);

  auto PSLMU = [](real ZZ) -> real { return real(-0.96) * rlog(real(1.0) - real(4.5) * ZZ); };
  auto PSLMS = [&](real ZZ) -> real { return ZZ / RFC - real(2.076) * (real(1.0) - real(1.0) / (ZZ + real(1.0))); };
  auto PSLHU = [](real ZZ) -> real { return real(-0.96) * rlog(real(1.0) - real(4.5) * ZZ); };
  auto PSLHS = [&](real ZZ) -> real { return ZZ * RFAC - real(2.076) * (real(1.0) - rexp(real(-1.2) * ZZ)); };
  auto PSPMU = [&](real XX) -> real {
    return real(-2.0) * rlog((XX + real(1.0)) * real(0.5)) -
           rlog((XX * XX + real(1.0)) * real(0.5)) + real(2.0) * std::atan(XX) - PIHF;
  };
  auto PSPMS = [](real YY) -> real { return real(5.0) * YY; };
  auto PSPHU = [](real XX) -> real { return real(-2.0) * rlog((XX * XX + real(1.0)) * real(0.5)); };
  auto PSPHS = [](real YY) -> real { return real(5.0) * YY; };

  int ILECH = 0;
  real ZILFC = -energy.param.ZilitinkevichCoeff * VKRM * SQVISC;
  real ZU = energy.state.RoughLenMomSfc;
  real RDZ = real(1.0) / energy.state.RefHeightAboveGrd;
  real CXCH = EXCM * RDZ;
  real DTHV = energy.state.TemperaturePotRefHeight - energy.state.TemperatureCanopyAir;

  real DU2 = std::max(energy.state.WindSpdRefHeight * energy.state.WindSpdRefHeight, EPSU2);
  real BTGH = BTG * HPBL;
  if (IterationInd == 1) {
    if ((BTGH * energy.state.ExchCoeffShAbvCan * DTHV) != real(0.0)) {
      energy.state.FrictionVelVertVeg =
          WWST2 * rpow(std::abs(BTGH * energy.state.ExchCoeffShAbvCan * DTHV), real(2.0) / real(3.0));
    } else {
      energy.state.FrictionVelVertVeg = real(0.0);
    }
    energy.state.FrictionVelVeg =
        std::max(std::sqrt(energy.state.ExchCoeffMomAbvCan * std::sqrt(DU2 + energy.state.FrictionVelVertVeg)),
                 EPSUST);
    energy.state.MoStabParaAbvCan =
        ELFC * energy.state.ExchCoeffShAbvCan * DTHV / cppwrf::powi(energy.state.FrictionVelVeg, 3);
  }

  real ZT = std::max(real(1.0e-6),
                     rexp(ZILFC * std::sqrt(energy.state.FrictionVelVeg * energy.state.RoughLenMomSfc)) *
                         energy.state.RoughLenMomSfc);
  real ZSLU = energy.state.RefHeightAboveGrd + ZU;
  real ZSLT = energy.state.RefHeightAboveGrd + ZT;
  real RLOGU = rlog(ZSLU / ZU);
  real RLOGT = rlog(ZSLT / ZT);

  real ZETALT = std::max(ZSLT * energy.state.MoStabParaAbvCan, ZTMIN);
  energy.state.MoStabParaAbvCan = ZETALT / ZSLT;
  real ZETALU = ZSLU * energy.state.MoStabParaAbvCan;
  real ZETAU = ZU * energy.state.MoStabParaAbvCan;
  real ZETAT = ZT * energy.state.MoStabParaAbvCan;
  real PSMZ = real(0.0);
  real SIMM = real(0.0);
  real PSHZ = real(0.0);
  real SIMH = real(0.0);
  if (ILECH == 0) {
    if (energy.state.MoStabParaAbvCan < real(0.0)) {
      real XLU4 = real(1.0) - real(16.0) * ZETALU;
      real XLT4 = real(1.0) - real(16.0) * ZETALT;
      real XU4 = real(1.0) - real(16.0) * ZETAU;
      real XT4 = real(1.0) - real(16.0) * ZETAT;
      real XLU = std::sqrt(std::sqrt(XLU4));
      real XLT = std::sqrt(std::sqrt(XLT4));
      real XU = std::sqrt(std::sqrt(XU4));
      real XT = std::sqrt(std::sqrt(XT4));
      PSMZ = PSPMU(XU);
      SIMM = PSPMU(XLU) - PSMZ + RLOGU;
      PSHZ = PSPHU(XT);
      SIMH = PSPHU(XLT) - PSHZ + RLOGT;
    } else {
      ZETALU = std::min(ZETALU, ZTMAX);
      ZETALT = std::min(ZETALT, ZTMAX);
      ZETAU = std::min(ZETAU, ZTMAX / (ZSLU / ZU));
      ZETAT = std::min(ZETAT, ZTMAX / (ZSLT / ZT));
      PSMZ = PSPMS(ZETAU);
      SIMM = PSPMS(ZETALU) - PSMZ + RLOGU;
      PSHZ = PSPHS(ZETAT);
      SIMH = PSPHS(ZETALT) - PSHZ + RLOGT;
    }
  } else {
    if (energy.state.MoStabParaAbvCan < real(0.0)) {
      PSMZ = PSLMU(ZETAU);
      SIMM = PSLMU(ZETALU) - PSMZ + RLOGU;
      PSHZ = PSLHU(ZETAT);
      SIMH = PSLHU(ZETALT) - PSHZ + RLOGT;
    } else {
      ZETALU = std::min(ZETALU, ZTMAX);
      ZETALT = std::min(ZETALT, ZTMAX);
      PSMZ = PSLMS(ZETAU);
      SIMM = PSLMS(ZETALU) - PSMZ + RLOGU;
      PSHZ = PSLHS(ZETAT);
      SIMH = PSLHS(ZETALT) - PSHZ + RLOGT;
    }
  }

  energy.state.FrictionVelVeg =
      std::max(std::sqrt(energy.state.ExchCoeffMomAbvCan * std::sqrt(DU2 + energy.state.FrictionVelVertVeg)),
               EPSUST);

  ZT = std::max(real(1.0e-6),
                rexp(ZILFC * std::sqrt(energy.state.FrictionVelVeg * energy.state.RoughLenMomSfc)) *
                    energy.state.RoughLenMomSfc);
  ZSLT = energy.state.RefHeightAboveGrd + ZT;
  RLOGT = rlog(ZSLT / ZT);
  real USTARK = energy.state.FrictionVelVeg * VKRM;

  if (SIMM < real(1.0e-6)) SIMM = real(1.0e-6);
  energy.state.ExchCoeffMomAbvCan = std::max(USTARK / SIMM, CXCH);
  if (SIMH < real(1.0e-6)) SIMH = real(1.0e-6);
  energy.state.ExchCoeffShAbvCan = std::max(USTARK / SIMH, CXCH);

  if ((BTGH * energy.state.ExchCoeffShAbvCan * DTHV) != real(0.0)) {
    energy.state.FrictionVelVertVeg =
        WWST2 * rpow(std::abs(BTGH * energy.state.ExchCoeffShAbvCan * DTHV), real(2.0) / real(3.0));
  } else {
    energy.state.FrictionVelVertVeg = real(0.0);
  }

  real RLMN = ELFC * energy.state.ExchCoeffShAbvCan * DTHV / cppwrf::powi(energy.state.FrictionVelVeg, 3);
  real RLMA = energy.state.MoStabParaAbvCan * WOLD + RLMN * WNEW;
  energy.state.MoStabParaAbvCan = RLMA;

  energy.state.ExchCoeffShAbvCan = energy.state.ExchCoeffShAbvCan / energy.state.WindSpdRefHeight;
  energy.state.ExchCoeffMomAbvCan = energy.state.ExchCoeffMomAbvCan / energy.state.WindSpdRefHeight;

  energy.state.ResistanceMomAbvCan =
      std::max(real(1.0), real(1.0) / (energy.state.ExchCoeffMomAbvCan * energy.state.WindSpdRefHeight));
  energy.state.ResistanceShAbvCan =
      std::max(real(1.0), real(1.0) / (energy.state.ExchCoeffShAbvCan * energy.state.WindSpdRefHeight));
  energy.state.ResistanceLhAbvCan = energy.state.ResistanceShAbvCan;
}

// ============================================================================
// ResistanceBareGroundMostMod.F90
// ============================================================================
CPPWRF_HD inline void ResistanceBareGroundMOST(NoahmpType& noahmp,
                                               int IndIter,
                                               real HeatSensibleTmp,
                                               int& MoStabParaSgn) {
  auto& forcing = noahmp.forcing;
  auto& energy = noahmp.energy;

  real MPE = real(1.0e-6);
  real MOZOLD = energy.state.MoStabParaBare;

  // gfortran SLP-vectorizes these 4 log calls into two 2-lane _ZGVbN2v_log calls
  real TMPCM  = cppwrf::vm::log((energy.state.RefHeightAboveGrd - energy.state.ZeroPlaneDispGrd) / energy.state.RoughLenMomGrd, 0, 2);
  real TMPCH  = cppwrf::vm::log((energy.state.RefHeightAboveGrd - energy.state.ZeroPlaneDispGrd) / energy.state.RoughLenShBareGrd, 1, 2);
  real TMPCM2 = cppwrf::vm::log((real(2.0) + energy.state.RoughLenMomGrd) / energy.state.RoughLenMomGrd, 0, 2);
  real TMPCH2 = cppwrf::vm::log((real(2.0) + energy.state.RoughLenShBareGrd) / energy.state.RoughLenShBareGrd, 1, 2);

  if (IndIter == 1) {
    energy.state.FrictionVelBare = real(0.0);
    energy.state.MoStabParaBare = real(0.0);
    energy.state.MoLengthBare = real(0.0);
    energy.state.MoStabParaBare2m = real(0.0);
  } else {
    real TVIR = (real(1.0) + real(0.61) * forcing.SpecHumidityRefHeight) * forcing.TemperatureAirRefHeight;
    real TMP1 = ConstVonKarman * (ConstGravityAcc / TVIR) * HeatSensibleTmp /
                (energy.state.DensityAirRefHeight * ConstHeatCapacAir);
    if (std::abs(TMP1) <= MPE) TMP1 = MPE;
    energy.state.MoLengthBare = real(-1.0) * cppwrf::powi(energy.state.FrictionVelBare, 3) / TMP1;
    energy.state.MoStabParaBare =
        std::min((energy.state.RefHeightAboveGrd - energy.state.ZeroPlaneDispGrd) / energy.state.MoLengthBare,
                 real(1.0));
    energy.state.MoStabParaBare2m =
        std::min((real(2.0) + energy.state.RoughLenShBareGrd) / energy.state.MoLengthBare, real(1.0));
  }

  if (MOZOLD * energy.state.MoStabParaBare < real(0.0)) MoStabParaSgn = MoStabParaSgn + 1;
  if (MoStabParaSgn >= 2) {
    energy.state.MoStabParaBare = real(0.0);
    energy.state.MoStabCorrMomBare = real(0.0);
    energy.state.MoStabCorrShBare = real(0.0);
    energy.state.MoStabParaBare2m = real(0.0);
    energy.state.MoStabCorrMomBare2m = real(0.0);
    energy.state.MoStabCorrShBare2m = real(0.0);
  }

  real FMNEW = real(0.0);
  real FHNEW = real(0.0);
  real FM2NEW = real(0.0);
  real FH2NEW = real(0.0);
  if (energy.state.MoStabParaBare < real(0.0)) {
    real TMP1 = rpow(real(1.0) - real(16.0) * energy.state.MoStabParaBare, real(0.25));
    real TMP2 = rlog((real(1.0) + TMP1 * TMP1) / real(2.0));
    real TMP3 = rlog((real(1.0) + TMP1) / real(2.0));
    FMNEW = real(2.0) * TMP3 + TMP2 - real(2.0) * std::atan(TMP1) + real(1.5707963);
    FHNEW = real(2.0) * TMP2;

    real TMP12 = rpow(real(1.0) - real(16.0) * energy.state.MoStabParaBare2m, real(0.25));
    real TMP22 = rlog((real(1.0) + TMP12 * TMP12) / real(2.0));
    real TMP32 = rlog((real(1.0) + TMP12) / real(2.0));
    FM2NEW = real(2.0) * TMP32 + TMP22 - real(2.0) * std::atan(TMP12) + real(1.5707963);
    FH2NEW = real(2.0) * TMP22;
  } else {
    FMNEW = real(-5.0) * energy.state.MoStabParaBare;
    FHNEW = FMNEW;
    FM2NEW = real(-5.0) * energy.state.MoStabParaBare2m;
    FH2NEW = FM2NEW;
  }

  if (IndIter == 1) {
    energy.state.MoStabCorrMomBare = FMNEW;
    energy.state.MoStabCorrShBare = FHNEW;
    energy.state.MoStabCorrMomBare2m = FM2NEW;
    energy.state.MoStabCorrShBare2m = FH2NEW;
  } else {
    energy.state.MoStabCorrMomBare = real(0.5) * (energy.state.MoStabCorrMomBare + FMNEW);
    energy.state.MoStabCorrShBare = real(0.5) * (energy.state.MoStabCorrShBare + FHNEW);
    energy.state.MoStabCorrMomBare2m = real(0.5) * (energy.state.MoStabCorrMomBare2m + FM2NEW);
    energy.state.MoStabCorrShBare2m = real(0.5) * (energy.state.MoStabCorrShBare2m + FH2NEW);
  }

  energy.state.MoStabCorrShBare = std::min(energy.state.MoStabCorrShBare, real(0.9) * TMPCH);
  energy.state.MoStabCorrMomBare = std::min(energy.state.MoStabCorrMomBare, real(0.9) * TMPCM);
  energy.state.MoStabCorrShBare2m = std::min(energy.state.MoStabCorrShBare2m, real(0.9) * TMPCH2);
  energy.state.MoStabCorrMomBare2m = std::min(energy.state.MoStabCorrMomBare2m, real(0.9) * TMPCM2);
  real CMFM = TMPCM - energy.state.MoStabCorrMomBare;
  real CHFH = TMPCH - energy.state.MoStabCorrShBare;
  real CM2FM2 = TMPCM2 - energy.state.MoStabCorrMomBare2m;
  real CH2FH2 = TMPCH2 - energy.state.MoStabCorrShBare2m;
  if (std::abs(CMFM) <= MPE) CMFM = MPE;
  if (std::abs(CHFH) <= MPE) CHFH = MPE;
  if (std::abs(CM2FM2) <= MPE) CM2FM2 = MPE;
  if (std::abs(CH2FH2) <= MPE) CH2FH2 = MPE;
  energy.state.ExchCoeffMomBare = ConstVonKarman * ConstVonKarman / (CMFM * CMFM);
  energy.state.ExchCoeffShBare = ConstVonKarman * ConstVonKarman / (CMFM * CHFH);

  energy.state.FrictionVelBare = energy.state.WindSpdRefHeight * std::sqrt(energy.state.ExchCoeffMomBare);
  energy.state.ExchCoeffSh2mBareMo = ConstVonKarman * energy.state.FrictionVelBare / CH2FH2;

  energy.state.ResistanceMomBareGrd =
      std::max(real(1.0), real(1.0) / (energy.state.ExchCoeffMomBare * energy.state.WindSpdRefHeight));
  energy.state.ResistanceShBareGrd =
      std::max(real(1.0), real(1.0) / (energy.state.ExchCoeffShBare * energy.state.WindSpdRefHeight));
  energy.state.ResistanceLhBareGrd = energy.state.ResistanceShBareGrd;
}

// ============================================================================
// ResistanceBareGroundChen97Mod.F90
// ============================================================================
CPPWRF_HD inline void ResistanceBareGroundChen97(NoahmpType& noahmp, int IndIter) {
  auto& water = noahmp.water;
  auto& energy = noahmp.energy;

  constexpr real WWST   = real(1.2);
  constexpr real WWST2  = WWST * WWST;
  constexpr real VKRM   = real(0.40);
  constexpr real EXCM   = real(0.001);
  constexpr real BETA   = real(1.0) / real(270.0);
  constexpr real BTG    = BETA * ConstGravityAcc;
  constexpr real ELFC   = VKRM * BTG;
  constexpr real WOLD   = real(0.15);
  constexpr real WNEW   = real(1.0) - WOLD;
  constexpr real PIHF   = real(3.14159265) / real(2.0);
  constexpr real EPSU2  = real(1.0e-4);
  constexpr real EPSUST = real(0.07);
  constexpr real ZTMIN  = real(-5.0);
  constexpr real ZTMAX  = real(1.0);
  constexpr real HPBL   = real(1000.0);
  constexpr real SQVISC = real(258.2);
  constexpr real RIC    = real(0.183);
  constexpr real FHNEU  = real(0.8);
  constexpr real RFC    = real(0.191);
  constexpr real RFAC   = RIC / (FHNEU * RFC * RFC);

  auto PSLMU = [](real ZZ) -> real { return real(-0.96) * rlog(real(1.0) - real(4.5) * ZZ); };
  auto PSLMS = [&](real ZZ) -> real { return ZZ / RFC - real(2.076) * (real(1.0) - real(1.0) / (ZZ + real(1.0))); };
  auto PSLHU = [](real ZZ) -> real { return real(-0.96) * rlog(real(1.0) - real(4.5) * ZZ); };
  auto PSLHS = [&](real ZZ) -> real { return ZZ * RFAC - real(2.076) * (real(1.0) - rexp(real(-1.2) * ZZ)); };
  auto PSPMU = [&](real XX) -> real {
    return real(-2.0) * rlog((XX + real(1.0)) * real(0.5)) -
           rlog((XX * XX + real(1.0)) * real(0.5)) + real(2.0) * std::atan(XX) - PIHF;
  };
  auto PSPMS = [](real YY) -> real { return real(5.0) * YY; };
  auto PSPHU = [](real XX) -> real { return real(-2.0) * rlog((XX * XX + real(1.0)) * real(0.5)); };
  auto PSPHS = [](real YY) -> real { return real(5.0) * YY; };

  int ILECH = 0;
  real ZILFC = -energy.param.ZilitinkevichCoeff * VKRM * SQVISC;
  real ZU = energy.state.RoughLenMomGrd;
  real RDZ = real(1.0) / energy.state.RefHeightAboveGrd;
  real CXCH = EXCM * RDZ;
  real DTHV = energy.state.TemperaturePotRefHeight - energy.state.TemperatureGrdBare;

  real DU2 = std::max(energy.state.WindSpdRefHeight * energy.state.WindSpdRefHeight, EPSU2);
  real BTGH = BTG * HPBL;
  if (IndIter == 1) {
    if ((BTGH * energy.state.ExchCoeffShBare * DTHV) != real(0.0)) {
      energy.state.FrictionVelVertBare =
          WWST2 * rpow(std::abs(BTGH * energy.state.ExchCoeffShBare * DTHV), real(2.0) / real(3.0));
    } else {
      energy.state.FrictionVelVertBare = real(0.0);
    }
    energy.state.FrictionVelBare =
        std::max(std::sqrt(energy.state.ExchCoeffMomBare * std::sqrt(DU2 + energy.state.FrictionVelVertBare)),
                 EPSUST);
    energy.state.MoStabParaBare =
        ELFC * energy.state.ExchCoeffShBare * DTHV / cppwrf::powi(energy.state.FrictionVelBare, 3);
  }

  real ZT = std::max(real(1.0e-6),
                     rexp(ZILFC * std::sqrt(energy.state.FrictionVelBare * energy.state.RoughLenMomGrd)) *
                         energy.state.RoughLenMomGrd);
  real ZSLU = energy.state.RefHeightAboveGrd + ZU;
  real ZSLT = energy.state.RefHeightAboveGrd + ZT;
  real RLOGU = rlog(ZSLU / ZU);
  real RLOGT = rlog(ZSLT / ZT);

  real ZETALT = std::max(ZSLT * energy.state.MoStabParaBare, ZTMIN);
  energy.state.MoStabParaBare = ZETALT / ZSLT;
  real ZETALU = ZSLU * energy.state.MoStabParaBare;
  real ZETAU = ZU * energy.state.MoStabParaBare;
  real ZETAT = ZT * energy.state.MoStabParaBare;
  real PSMZ = real(0.0);
  real SIMM = real(0.0);
  real PSHZ = real(0.0);
  real SIMH = real(0.0);
  if (ILECH == 0) {
    if (energy.state.MoStabParaBare < real(0.0)) {
      real XLU4 = real(1.0) - real(16.0) * ZETALU;
      real XLT4 = real(1.0) - real(16.0) * ZETALT;
      real XU4 = real(1.0) - real(16.0) * ZETAU;
      real XT4 = real(1.0) - real(16.0) * ZETAT;
      real XLU = std::sqrt(std::sqrt(XLU4));
      real XLT = std::sqrt(std::sqrt(XLT4));
      real XU = std::sqrt(std::sqrt(XU4));
      real XT = std::sqrt(std::sqrt(XT4));
      PSMZ = PSPMU(XU);
      SIMM = PSPMU(XLU) - PSMZ + RLOGU;
      PSHZ = PSPHU(XT);
      SIMH = PSPHU(XLT) - PSHZ + RLOGT;
    } else {
      ZETALU = std::min(ZETALU, ZTMAX);
      ZETALT = std::min(ZETALT, ZTMAX);
      ZETAU = std::min(ZETAU, ZTMAX / (ZSLU / ZU));
      ZETAT = std::min(ZETAT, ZTMAX / (ZSLT / ZT));
      PSMZ = PSPMS(ZETAU);
      SIMM = PSPMS(ZETALU) - PSMZ + RLOGU;
      PSHZ = PSPHS(ZETAT);
      SIMH = PSPHS(ZETALT) - PSHZ + RLOGT;
    }
  } else {
    if (energy.state.MoStabParaBare < real(0.0)) {
      PSMZ = PSLMU(ZETAU);
      SIMM = PSLMU(ZETALU) - PSMZ + RLOGU;
      PSHZ = PSLHU(ZETAT);
      SIMH = PSLHU(ZETALT) - PSHZ + RLOGT;
    } else {
      ZETALU = std::min(ZETALU, ZTMAX);
      ZETALT = std::min(ZETALT, ZTMAX);
      PSMZ = PSLMS(ZETAU);
      SIMM = PSLMS(ZETALU) - PSMZ + RLOGU;
      PSHZ = PSLHS(ZETAT);
      SIMH = PSLHS(ZETALT) - PSHZ + RLOGT;
    }
  }

  energy.state.FrictionVelBare =
      std::max(std::sqrt(energy.state.ExchCoeffMomBare * std::sqrt(DU2 + energy.state.FrictionVelVertBare)),
               EPSUST);

  ZT = std::max(real(1.0e-6),
                rexp(ZILFC * std::sqrt(energy.state.FrictionVelBare * energy.state.RoughLenMomGrd)) *
                    energy.state.RoughLenMomGrd);
  ZSLT = energy.state.RefHeightAboveGrd + ZT;
  RLOGT = rlog(ZSLT / ZT);
  real USTARK = energy.state.FrictionVelBare * VKRM;

  if (SIMM < real(1.0e-6)) SIMM = real(1.0e-6);
  energy.state.ExchCoeffMomBare = std::max(USTARK / SIMM, CXCH);
  if (SIMH < real(1.0e-6)) SIMH = real(1.0e-6);
  energy.state.ExchCoeffShBare = std::max(USTARK / SIMH, CXCH);

  if (BTGH * energy.state.ExchCoeffShBare * DTHV != real(0.0)) {
    energy.state.FrictionVelVertBare =
        WWST2 * rpow(std::abs(BTGH * energy.state.ExchCoeffShBare * DTHV), real(2.0) / real(3.0));
  } else {
    energy.state.FrictionVelVertBare = real(0.0);
  }

  real RLMN = ELFC * energy.state.ExchCoeffShBare * DTHV / cppwrf::powi(energy.state.FrictionVelBare, 3);
  real RLMA = energy.state.MoStabParaBare * WOLD + RLMN * WNEW;
  energy.state.MoStabParaBare = RLMA;

  energy.state.ExchCoeffShBare = energy.state.ExchCoeffShBare / energy.state.WindSpdRefHeight;
  energy.state.ExchCoeffMomBare = energy.state.ExchCoeffMomBare / energy.state.WindSpdRefHeight;
  if (water.state.SnowDepth > real(0.0)) {
    energy.state.ExchCoeffMomBare = std::min(real(0.01), energy.state.ExchCoeffMomBare);
    energy.state.ExchCoeffShBare = std::min(real(0.01), energy.state.ExchCoeffShBare);
  }

  energy.state.ResistanceMomBareGrd =
      std::max(real(1.0), real(1.0) / (energy.state.ExchCoeffMomBare * energy.state.WindSpdRefHeight));
  energy.state.ResistanceShBareGrd =
      std::max(real(1.0), real(1.0) / (energy.state.ExchCoeffShBare * energy.state.WindSpdRefHeight));
  energy.state.ResistanceLhBareGrd = energy.state.ResistanceShBareGrd;
}

// ============================================================================
// SurfaceEnergyFluxVegetatedMod.F90
// ============================================================================
CPPWRF_HD inline void SurfaceEnergyFluxVegetated(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& nml = noahmp.config.nmlist;
  auto& forcing = noahmp.forcing;
  auto& water = noahmp.water;
  auto& energy = noahmp.energy;

  constexpr int NumIterC = 20;
  constexpr int NumIterG = 5;

  auto TempUnitConv = [](real TempTmp) -> real {
    return std::min(real(50.0), std::max(real(-50.0), TempTmp - ConstFreezePoint));
  };

  int LastIter = 0;
  energy.state.FrictionVelVeg = real(0.1);
  real TemperatureCanChg = real(0.0);
  real TemperatureGrdChg = real(0.0);
  energy.state.MoStabParaAbvCan = real(0.0);
  int MoStabParaSgn = 0;
  energy.state.MoStabCorrShVeg2m = real(0.0);
  real ShGrdTmp = real(0.0);
  real ShCanTmp = real(0.0);

  real VegAreaIndTmp = std::min(real(6.0), energy.state.VegAreaIndEff);
  real LeafAreaIndSunEff = std::min(real(6.0), energy.state.LeafAreaIndSunlit);
  real LeafAreaIndShdEff = std::min(real(6.0), energy.state.LeafAreaIndShade);

  real TempTmp = TempUnitConv(energy.state.TemperatureGrdVeg);
  real VapPresSatWatTmp = real(0.0);
  real VapPresSatIceTmp = real(0.0);
  real VapPresSatWatTmpD = real(0.0);
  real VapPresSatIceTmpD = real(0.0);
  VaporPressureSaturation(TempTmp, VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD);
  if (TempTmp > real(0.0)) {
    energy.state.VapPresSatGrdVeg = VapPresSatWatTmp;
  } else {
    energy.state.VapPresSatGrdVeg = VapPresSatIceTmp;
  }

  energy.state.CanopyHeight = energy.param.HeightCanopyTop;
  energy.state.WindSpdCanopyTop =
      energy.state.WindSpdRefHeight *
      rlog((energy.state.CanopyHeight - energy.state.ZeroPlaneDispSfc + energy.state.RoughLenMomSfc) /
           energy.state.RoughLenMomSfc) /
      rlog(energy.state.RefHeightAboveGrd / energy.state.RoughLenMomSfc);

  real LwCoeffAir =
      -energy.state.EmissivityVeg *
          (real(1.0) + (real(1.0) - energy.state.EmissivityVeg) * (real(1.0) - energy.state.EmissivityGrd)) *
          forcing.RadLwDownRefHeight -
      energy.state.EmissivityVeg * energy.state.EmissivityGrd * ConstStefanBoltzmann *
          cppwrf::powi(energy.state.TemperatureGrdVeg, 4);
  real LwCoeffCan =
      (real(2.0) - energy.state.EmissivityVeg * (real(1.0) - energy.state.EmissivityGrd)) *
      energy.state.EmissivityVeg * ConstStefanBoltzmann;

  real ExchCoeffShAbvCanTmp = real(0.0);
  real ExchCoeffShLeafTmp = real(0.0);
  for (int IndIter = 1; IndIter <= NumIterC; ++IndIter) {
    if (IndIter == 1) {
      energy.state.RoughLenShCanopy = energy.state.RoughLenMomSfc;
      energy.state.RoughLenShVegGrd = energy.state.RoughLenMomGrd;
    } else {
      energy.state.RoughLenShCanopy = energy.state.RoughLenMomSfc;
      energy.state.RoughLenShVegGrd = energy.state.RoughLenMomGrd;
    }

    if (nml.OptSurfaceDrag == 1) ResistanceAboveCanopyMOST(noahmp, IndIter, ShCanTmp, MoStabParaSgn);
    if (nml.OptSurfaceDrag == 2) ResistanceAboveCanopyChen97(noahmp, IndIter);

    ResistanceLeafToGround(noahmp, IndIter, VegAreaIndTmp, ShGrdTmp);

    TempTmp = TempUnitConv(energy.state.TemperatureCanopy);
    VaporPressureSaturation(TempTmp, VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD);
    if (TempTmp > real(0.0)) {
      energy.state.VapPresSatCanopy = VapPresSatWatTmp;
      energy.state.VapPresSatCanTempD = VapPresSatWatTmpD;
    } else {
      energy.state.VapPresSatCanopy = VapPresSatIceTmp;
      energy.state.VapPresSatCanTempD = VapPresSatIceTmpD;
    }

    if (IndIter == 1) {
      if (nml.OptStomataResistance == 1) {
        ResistanceCanopyStomataBallBerry(noahmp, 0);
        ResistanceCanopyStomataBallBerry(noahmp, 1);
      }
      if (nml.OptStomataResistance == 2) {
        ResistanceCanopyStomataJarvis(noahmp, 0);
        ResistanceCanopyStomataJarvis(noahmp, 1);
      }
    }

    ExchCoeffShAbvCanTmp = real(1.0) / energy.state.ResistanceShAbvCan;
    ExchCoeffShLeafTmp = real(2.0) * VegAreaIndTmp / energy.state.ResistanceLeafBoundary;
    real GrdHeatCoeff = real(1.0) / energy.state.ResistanceShUndCan;
    real ExchCoeffTot = ExchCoeffShAbvCanTmp + ExchCoeffShLeafTmp + GrdHeatCoeff;
    real TempShGhTmp =
        (forcing.TemperatureAirRefHeight * ExchCoeffShAbvCanTmp + energy.state.TemperatureGrdVeg * GrdHeatCoeff) /
        ExchCoeffTot;
    real ExchCoeffShFrac = ExchCoeffShLeafTmp / ExchCoeffTot;
    real ShCoeff =
        (real(1.0) - ExchCoeffShFrac) * energy.state.DensityAirRefHeight * ConstHeatCapacAir * ExchCoeffShLeafTmp;

    energy.state.ExchCoeffLhAbvCan = real(1.0) / energy.state.ResistanceLhAbvCan;
    energy.state.ExchCoeffLhEvap = water.state.CanopyWetFrac * VegAreaIndTmp / energy.state.ResistanceLeafBoundary;
    energy.state.ExchCoeffLhTransp =
        (real(1.0) - water.state.CanopyWetFrac) *
        (LeafAreaIndSunEff / (energy.state.ResistanceLeafBoundary + energy.state.ResistanceStomataSunlit) +
         LeafAreaIndShdEff / (energy.state.ResistanceLeafBoundary + energy.state.ResistanceStomataShade));
    energy.state.ExchCoeffLhUndCan = real(1.0) / (energy.state.ResistanceLhUndCan + energy.state.ResistanceGrdEvap);
    ExchCoeffTot = energy.state.ExchCoeffLhAbvCan + energy.state.ExchCoeffLhEvap +
                   energy.state.ExchCoeffLhTransp + energy.state.ExchCoeffLhUndCan;
    real VapPresLhTot =
        (energy.state.PressureVaporRefHeight * energy.state.ExchCoeffLhAbvCan +
         energy.state.VapPresSatGrdVeg * energy.state.ExchCoeffLhUndCan) /
        ExchCoeffTot;
    real ExchCoeffEtFrac = (energy.state.ExchCoeffLhEvap + energy.state.ExchCoeffLhTransp) / ExchCoeffTot;
    real LhCoeff =
        (real(1.0) - ExchCoeffEtFrac) * energy.state.ExchCoeffLhEvap * energy.state.DensityAirRefHeight *
        ConstHeatCapacAir / energy.state.PsychConstCanopy;
    real TranspHeatCoeff =
        (real(1.0) - ExchCoeffEtFrac) * energy.state.ExchCoeffLhTransp * energy.state.DensityAirRefHeight *
        ConstHeatCapacAir / energy.state.PsychConstCanopy;

    energy.state.TemperatureCanopyAir = TempShGhTmp + ExchCoeffShFrac * energy.state.TemperatureCanopy;
    energy.state.PressureVaporCanAir = VapPresLhTot + ExchCoeffEtFrac * energy.state.VapPresSatCanopy;
    energy.flux.RadLwNetCanopy =
        energy.state.VegFrac * (LwCoeffAir + LwCoeffCan * cppwrf::powi(energy.state.TemperatureCanopy, 4));
    energy.flux.HeatSensibleCanopy =
        energy.state.VegFrac * energy.state.DensityAirRefHeight * ConstHeatCapacAir * ExchCoeffShLeafTmp *
        (energy.state.TemperatureCanopy - energy.state.TemperatureCanopyAir);
    energy.flux.HeatLatentCanEvap =
        energy.state.VegFrac * energy.state.DensityAirRefHeight * ConstHeatCapacAir * energy.state.ExchCoeffLhEvap *
        (energy.state.VapPresSatCanopy - energy.state.PressureVaporCanAir) / energy.state.PsychConstCanopy;
    energy.flux.HeatLatentCanTransp =
        energy.state.VegFrac * energy.state.DensityAirRefHeight * ConstHeatCapacAir * energy.state.ExchCoeffLhTransp *
        (energy.state.VapPresSatCanopy - energy.state.PressureVaporCanAir) / energy.state.PsychConstCanopy;
    if (energy.state.TemperatureCanopy > ConstFreezePoint) {
      energy.flux.HeatLatentCanEvap =
          std::min(water.state.CanopyLiqWater * energy.state.LatHeatVapCanopy / domain.MainTimeStep,
                   energy.flux.HeatLatentCanEvap);
    } else {
      energy.flux.HeatLatentCanEvap =
          std::min(water.state.CanopyIce * energy.state.LatHeatVapCanopy / domain.MainTimeStep,
                   energy.flux.HeatLatentCanEvap);
    }

    real HeatCapacCan =
        energy.param.HeatCapacCanFac * VegAreaIndTmp * ConstHeatCapacWater +
        water.state.CanopyLiqWater * ConstHeatCapacWater / ConstDensityWater +
        water.state.CanopyIce * ConstHeatCapacIce / ConstDensityIce;
    real EnergyResTmp =
        energy.flux.RadSwAbsVeg - energy.flux.RadLwNetCanopy - energy.flux.HeatSensibleCanopy -
        energy.flux.HeatLatentCanEvap - energy.flux.HeatLatentCanTransp + energy.flux.HeatPrecipAdvCanopy;
    real FluxTotCoeff =
        energy.state.VegFrac *
        (real(4.0) * LwCoeffCan * cppwrf::powi(energy.state.TemperatureCanopy, 3) + ShCoeff +
         (LhCoeff + TranspHeatCoeff) * energy.state.VapPresSatCanTempD + HeatCapacCan / domain.MainTimeStep);
    TemperatureCanChg = EnergyResTmp / FluxTotCoeff;

    energy.flux.RadLwNetCanopy =
        energy.flux.RadLwNetCanopy +
        energy.state.VegFrac * real(4.0) * LwCoeffCan * cppwrf::powi(energy.state.TemperatureCanopy, 3) *
            TemperatureCanChg;
    energy.flux.HeatSensibleCanopy =
        energy.flux.HeatSensibleCanopy + energy.state.VegFrac * ShCoeff * TemperatureCanChg;
    energy.flux.HeatLatentCanEvap =
        energy.flux.HeatLatentCanEvap +
        energy.state.VegFrac * LhCoeff * energy.state.VapPresSatCanTempD * TemperatureCanChg;
    energy.flux.HeatLatentCanTransp =
        energy.flux.HeatLatentCanTransp +
        energy.state.VegFrac * TranspHeatCoeff * energy.state.VapPresSatCanTempD * TemperatureCanChg;
    energy.flux.HeatCanStorageChg =
        energy.state.VegFrac * HeatCapacCan / domain.MainTimeStep * TemperatureCanChg;
    energy.state.TemperatureCanopy = energy.state.TemperatureCanopy + TemperatureCanChg;

    ShCanTmp = energy.state.DensityAirRefHeight * ConstHeatCapacAir *
               (energy.state.TemperatureCanopyAir - forcing.TemperatureAirRefHeight) /
               energy.state.ResistanceShAbvCan;
    ShGrdTmp = energy.state.DensityAirRefHeight * ConstHeatCapacAir *
               (energy.state.TemperatureGrdVeg - energy.state.TemperatureCanopyAir) /
               energy.state.ResistanceShUndCan;

    energy.state.SpecHumiditySfc =
        (real(0.622) * energy.state.PressureVaporCanAir) /
        (forcing.PressureAirRefHeight - real(0.378) * energy.state.PressureVaporCanAir);
    if (LastIter == 1) {
      break;
    }
    if (IndIter >= 5 && std::abs(TemperatureCanChg) <= real(0.01) && LastIter == 0) {
      LastIter = 1;
    }
  }

  LwCoeffAir = -energy.state.EmissivityGrd * (real(1.0) - energy.state.EmissivityVeg) * forcing.RadLwDownRefHeight -
               energy.state.EmissivityGrd * energy.state.EmissivityVeg * ConstStefanBoltzmann *
                   cppwrf::powi(energy.state.TemperatureCanopy, 4);
  LwCoeffCan = energy.state.EmissivityGrd * ConstStefanBoltzmann;
  real ShCoeff = energy.state.DensityAirRefHeight * ConstHeatCapacAir / energy.state.ResistanceShUndCan;
  real LhCoeff = energy.state.DensityAirRefHeight * ConstHeatCapacAir /
                 (energy.state.PsychConstGrd * (energy.state.ResistanceLhUndCan + energy.state.ResistanceGrdEvap));
  real GrdHeatCoeff =
      real(2.0) * energy.state.ThermConductSoilSnow(domain.NumSnowLayerNeg + 1) /
      domain.ThicknessSnowSoilLayer(domain.NumSnowLayerNeg + 1);

  for (int IndIter = 1; IndIter <= NumIterG; ++IndIter) {
    TempTmp = TempUnitConv(energy.state.TemperatureGrdVeg);
    VaporPressureSaturation(TempTmp, VapPresSatWatTmp, VapPresSatIceTmp, VapPresSatWatTmpD, VapPresSatIceTmpD);
    if (TempTmp > real(0.0)) {
      energy.state.VapPresSatGrdVeg = VapPresSatWatTmp;
      energy.state.VapPresSatGrdVegTempD = VapPresSatWatTmpD;
    } else {
      energy.state.VapPresSatGrdVeg = VapPresSatIceTmp;
      energy.state.VapPresSatGrdVegTempD = VapPresSatIceTmpD;
    }
    energy.flux.RadLwNetVegGrd = LwCoeffCan * cppwrf::powi(energy.state.TemperatureGrdVeg, 4) + LwCoeffAir;
    energy.flux.HeatSensibleVegGrd = ShCoeff * (energy.state.TemperatureGrdVeg - energy.state.TemperatureCanopyAir);
    energy.flux.HeatLatentVegGrd =
        LhCoeff * (energy.state.VapPresSatGrdVeg * energy.state.RelHumidityGrd - energy.state.PressureVaporCanAir);
    energy.flux.HeatGroundVegGrd =
        GrdHeatCoeff *
        (energy.state.TemperatureGrdVeg - energy.state.TemperatureSoilSnow(domain.NumSnowLayerNeg + 1));
    real EnergyResTmp =
        energy.flux.RadSwAbsGrd - energy.flux.RadLwNetVegGrd - energy.flux.HeatSensibleVegGrd -
        energy.flux.HeatLatentVegGrd - energy.flux.HeatGroundVegGrd + energy.flux.HeatPrecipAdvVegGrd;
    real FluxTotCoeff =
        real(4.0) * LwCoeffCan * cppwrf::powi(energy.state.TemperatureGrdVeg, 3) + ShCoeff +
        LhCoeff * energy.state.VapPresSatGrdVegTempD + GrdHeatCoeff;
    TemperatureGrdChg = EnergyResTmp / FluxTotCoeff;
    energy.flux.RadLwNetVegGrd =
        energy.flux.RadLwNetVegGrd +
        real(4.0) * LwCoeffCan * cppwrf::powi(energy.state.TemperatureGrdVeg, 3) * TemperatureGrdChg;
    energy.flux.HeatSensibleVegGrd = energy.flux.HeatSensibleVegGrd + ShCoeff * TemperatureGrdChg;
    energy.flux.HeatLatentVegGrd =
        energy.flux.HeatLatentVegGrd + LhCoeff * energy.state.VapPresSatGrdVegTempD * TemperatureGrdChg;
    energy.flux.HeatGroundVegGrd = energy.flux.HeatGroundVegGrd + GrdHeatCoeff * TemperatureGrdChg;
    energy.state.TemperatureGrdVeg = energy.state.TemperatureGrdVeg + TemperatureGrdChg;
  }

  if (nml.OptSnowSoilTempTime == 1 || nml.OptSnowSoilTempTime == 3) {
    if (water.state.SnowDepth > real(0.05) && energy.state.TemperatureGrdVeg > ConstFreezePoint) {
      if (nml.OptSnowSoilTempTime == 1) {
        energy.state.TemperatureGrdVeg = ConstFreezePoint;
      }
      if (nml.OptSnowSoilTempTime == 3) {
        energy.state.TemperatureGrdVeg =
            (real(1.0) - water.state.SnowCoverFrac) * energy.state.TemperatureGrdVeg +
            water.state.SnowCoverFrac * ConstFreezePoint;
      }

      energy.flux.RadLwNetVegGrd =
          LwCoeffCan * cppwrf::powi(energy.state.TemperatureGrdVeg, 4) -
          energy.state.EmissivityGrd * (real(1.0) - energy.state.EmissivityVeg) * forcing.RadLwDownRefHeight -
          energy.state.EmissivityGrd * energy.state.EmissivityVeg * ConstStefanBoltzmann *
              cppwrf::powi(energy.state.TemperatureCanopy, 4);
      energy.flux.HeatSensibleVegGrd = ShCoeff * (energy.state.TemperatureGrdVeg - energy.state.TemperatureCanopyAir);
      energy.flux.HeatLatentVegGrd =
          LhCoeff * (energy.state.VapPresSatGrdVeg * energy.state.RelHumidityGrd - energy.state.PressureVaporCanAir);
      energy.flux.HeatGroundVegGrd =
          energy.flux.RadSwAbsGrd + energy.flux.HeatPrecipAdvVegGrd -
          (energy.flux.RadLwNetVegGrd + energy.flux.HeatSensibleVegGrd + energy.flux.HeatLatentVegGrd);
    }
  }

  energy.state.WindStressEwVeg =
      -energy.state.DensityAirRefHeight * energy.state.ExchCoeffMomAbvCan * energy.state.WindSpdRefHeight *
      forcing.WindEastwardRefHeight;
  energy.state.WindStressNsVeg =
      -energy.state.DensityAirRefHeight * energy.state.ExchCoeffMomAbvCan * energy.state.WindSpdRefHeight *
      forcing.WindNorthwardRefHeight;

  if (nml.OptSurfaceDrag == 1 || nml.OptSurfaceDrag == 2) {
    energy.state.ExchCoeffSh2mVeg =
        energy.state.FrictionVelVeg * ConstVonKarman /
        (rlog((real(2.0) + energy.state.RoughLenShCanopy) / energy.state.RoughLenShCanopy) -
         energy.state.MoStabCorrShVeg2m);
    if (energy.state.ExchCoeffSh2mVeg < real(1.0e-5)) {
      energy.state.TemperatureAir2mVeg = energy.state.TemperatureCanopyAir;
      energy.state.SpecHumidity2mVeg = energy.state.SpecHumiditySfc;
    } else {
      energy.state.TemperatureAir2mVeg =
          energy.state.TemperatureCanopyAir -
          (energy.flux.HeatSensibleVegGrd + energy.flux.HeatSensibleCanopy / energy.state.VegFrac) /
              (energy.state.DensityAirRefHeight * ConstHeatCapacAir) * real(1.0) / energy.state.ExchCoeffSh2mVeg;
      energy.state.SpecHumidity2mVeg =
          energy.state.SpecHumiditySfc -
          ((energy.flux.HeatLatentCanEvap + energy.flux.HeatLatentCanTransp) / energy.state.VegFrac +
           energy.flux.HeatLatentVegGrd) /
              (energy.state.LatHeatVapCanopy * energy.state.DensityAirRefHeight) * real(1.0) /
              energy.state.ExchCoeffSh2mVeg;
    }
  }

  energy.state.ExchCoeffShAbvCan = ExchCoeffShAbvCanTmp;
  energy.state.ExchCoeffShLeaf = ExchCoeffShLeafTmp;
  energy.state.ExchCoeffShUndCan = real(1.0) / energy.state.ResistanceShUndCan;
}

// ============================================================================
// SurfaceEnergyFluxBareGroundMod.F90
// ============================================================================
CPPWRF_HD inline void SurfaceEnergyFluxBareGround(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& nml = noahmp.config.nmlist;
  auto& forcing = noahmp.forcing;
  auto& water = noahmp.water;
  auto& energy = noahmp.energy;

  constexpr int NumIter = 5;

  auto TempUnitConv = [](real TempTmp) -> real {
    return std::min(real(50.0), std::max(real(-50.0), TempTmp - ConstFreezePoint));
  };

  real TemperatureGrdChg = real(0.0);
  energy.state.MoStabParaBare = real(0.0);
  int MoStabParaSgn = 0;
  energy.state.MoStabCorrShBare2m = real(0.0);
  energy.state.FrictionVelBare = real(0.1);
  real HeatSensibleTmp = real(0.0);
  real LwRadCoeff = energy.state.EmissivityGrd * ConstStefanBoltzmann;
  real GrdHeatCoeff =
      real(2.0) * energy.state.ThermConductSoilSnow(domain.NumSnowLayerNeg + 1) /
      domain.ThicknessSnowSoilLayer(domain.NumSnowLayerNeg + 1);

  real ExchCoeffShTmp = real(0.0);
  real ShCoeff = real(0.0);
  real LhCoeff = real(0.0);
  for (int IndIter = 1; IndIter <= NumIter; ++IndIter) {
    if (IndIter == 1) {
      energy.state.RoughLenShBareGrd = energy.state.RoughLenMomGrd;
    } else {
      energy.state.RoughLenShBareGrd = energy.state.RoughLenMomGrd;
    }

    if (nml.OptSurfaceDrag == 1) ResistanceBareGroundMOST(noahmp, IndIter, HeatSensibleTmp, MoStabParaSgn);
    if (nml.OptSurfaceDrag == 2) ResistanceBareGroundChen97(noahmp, IndIter);

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
    LhCoeff = energy.state.DensityAirRefHeight * ConstHeatCapacAir / energy.state.PsychConstGrd /
              (energy.state.ResistanceGrdEvap + energy.state.ResistanceLhBareGrd);
    energy.flux.RadLwNetBareGrd =
        LwRadCoeff * cppwrf::powi(energy.state.TemperatureGrdBare, 4) -
        energy.state.EmissivityGrd * forcing.RadLwDownRefHeight;
    energy.flux.HeatSensibleBareGrd = ShCoeff * (energy.state.TemperatureGrdBare - forcing.TemperatureAirRefHeight);
    energy.flux.HeatLatentBareGrd =
        LhCoeff * (energy.state.VapPresSatGrdBare * energy.state.RelHumidityGrd - energy.state.PressureVaporRefHeight);
    energy.flux.HeatGroundBareGrd =
        GrdHeatCoeff *
        (energy.state.TemperatureGrdBare - energy.state.TemperatureSoilSnow(domain.NumSnowLayerNeg + 1));
    real EnergyResTmp =
        energy.flux.RadSwAbsGrd - energy.flux.RadLwNetBareGrd - energy.flux.HeatSensibleBareGrd -
        energy.flux.HeatLatentBareGrd - energy.flux.HeatGroundBareGrd + energy.flux.HeatPrecipAdvBareGrd;
    real FluxTotCoeff =
        real(4.0) * LwRadCoeff * cppwrf::powi(energy.state.TemperatureGrdBare, 3) + ShCoeff +
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

  if (nml.OptSnowSoilTempTime == 1 || nml.OptSnowSoilTempTime == 3) {
    if (water.state.SnowDepth > real(0.05) && energy.state.TemperatureGrdBare > ConstFreezePoint) {
      if (nml.OptSnowSoilTempTime == 1) {
        energy.state.TemperatureGrdBare = ConstFreezePoint;
      }
      if (nml.OptSnowSoilTempTime == 3) {
        energy.state.TemperatureGrdBare =
            (real(1.0) - water.state.SnowCoverFrac) * energy.state.TemperatureGrdBare +
            water.state.SnowCoverFrac * ConstFreezePoint;
      }

      energy.flux.RadLwNetBareGrd =
          LwRadCoeff * cppwrf::powi(energy.state.TemperatureGrdBare, 4) -
          energy.state.EmissivityGrd * forcing.RadLwDownRefHeight;
      energy.flux.HeatSensibleBareGrd = ShCoeff * (energy.state.TemperatureGrdBare - forcing.TemperatureAirRefHeight);
      energy.flux.HeatLatentBareGrd =
          LhCoeff *
          (energy.state.VapPresSatGrdBare * energy.state.RelHumidityGrd - energy.state.PressureVaporRefHeight);
      energy.flux.HeatGroundBareGrd =
          energy.flux.RadSwAbsGrd + energy.flux.HeatPrecipAdvBareGrd -
          (energy.flux.RadLwNetBareGrd + energy.flux.HeatSensibleBareGrd + energy.flux.HeatLatentBareGrd);
    }
  }

  energy.state.WindStressEwBare =
      -energy.state.DensityAirRefHeight * energy.state.ExchCoeffMomBare * energy.state.WindSpdRefHeight *
      forcing.WindEastwardRefHeight;
  energy.state.WindStressNsBare =
      -energy.state.DensityAirRefHeight * energy.state.ExchCoeffMomBare * energy.state.WindSpdRefHeight *
      forcing.WindNorthwardRefHeight;

  if (nml.OptSurfaceDrag == 1 || nml.OptSurfaceDrag == 2) {
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
          energy.flux.HeatSensibleBareGrd / (energy.state.DensityAirRefHeight * ConstHeatCapacAir) * real(1.0) /
              energy.state.ExchCoeffSh2mBare;
      energy.state.SpecHumidity2mBare =
          energy.state.SpecHumiditySfc -
          energy.flux.HeatLatentBareGrd / (energy.state.LatHeatVapGrd * energy.state.DensityAirRefHeight) *
              (real(1.0) / energy.state.ExchCoeffSh2mBare + energy.state.ResistanceGrdEvap);
    }
    if (domain.FlagUrban) energy.state.SpecHumidity2mBare = energy.state.SpecHumiditySfc;
  }

  energy.state.ExchCoeffShBare = ExchCoeffShTmp;
}

// ============================================================================
// SoilSnowThermalDiffusionMod.F90
// ============================================================================
template <typename Arr1, typename Arr2, typename Arr3, typename Arr4>
CPPWRF_HD inline void SoilSnowThermalDiffusion(NoahmpType& noahmp,
                                               Arr1& MatLeft1,
                                               Arr2& MatLeft2,
                                               Arr3& MatLeft3,
                                               Arr4& MatRight) {
  auto& domain = noahmp.config.domain;
  auto& nml = noahmp.config.nmlist;
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
      EnergyExcess(LoopInd) =
          energy.state.ThermConductSoilSnow(LoopInd) * TempGradDepth(LoopInd) -
          energy.flux.HeatGroundTotMean - energy.flux.RadSwPenetrateGrd(LoopInd);
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
      if (nml.OptSoilTemperatureBottom == 1) {
        energy.flux.HeatFromSoilBot = real(0.0);
      }
      if (nml.OptSoilTemperatureBottom == 2) {
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
      if (nml.OptSnowSoilTempTime == 1 || nml.OptSnowSoilTempTime == 3) {
        MatLeft2(LoopInd) = -MatLeft3(LoopInd);
      }
      if (nml.OptSnowSoilTempTime == 2) {
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
// SoilSnowTemperatureSolverMod.F90
// ============================================================================
template <typename Arr1, typename Arr2, typename Arr3, typename Arr4>
CPPWRF_HD inline void SoilSnowTemperatureSolver(NoahmpType& noahmp,
                                                real TimeStep,
                                                Arr1& MatLeft1,
                                                Arr2& MatLeft2,
                                                Arr3& MatLeft3,
                                                Arr4& MatRight) {
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
    energy.state.TemperatureSoilSnow(LoopInd) = energy.state.TemperatureSoilSnow(LoopInd) + MatLeft3(LoopInd);
  }
}

// ============================================================================
// SoilSnowTemperatureMainMod.F90
// ============================================================================
CPPWRF_HD inline void SoilSnowTemperatureMain(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& water = noahmp.water;
  auto& energy = noahmp.energy;

  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MatRight;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MatLeft1;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MatLeft2;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MatLeft3;
  MatRight.fill(real(0.0));
  MatLeft1.fill(real(0.0));
  MatLeft2.fill(real(0.0));
  MatLeft3.fill(real(0.0));

  for (int IndLoop = -domain.NumSnowLayerMax + 1; IndLoop <= domain.NumSoilLayer; ++IndLoop) {
    energy.flux.RadSwPenetrateGrd(IndLoop) = real(0.0);
  }

  energy.state.DepthSoilTempBotToSno = domain.DepthSoilTempBottom - water.state.SnowDepth;

  SoilSnowThermalDiffusion(noahmp, MatLeft1, MatLeft2, MatLeft3, MatRight);
  SoilSnowTemperatureSolver(noahmp, domain.SoilTimeStep, MatLeft1, MatLeft2, MatLeft3, MatRight);

  energy.flux.HeatFromSoilBot = energy.flux.HeatFromSoilBot * domain.SoilTimeStep;
}

// ============================================================================
// SoilWaterSupercoolNiu06Mod.F90
// ============================================================================
CPPWRF_HD inline void SoilWaterSupercoolNiu06(NoahmpType& noahmp,
                                              int IndSoil,
                                              real& SoilWatSupercool,
                                              real SoilTemperature) {
  auto& water = noahmp.water;

  real SoilWatPotFrz = ConstLatHeatFusion * (ConstFreezePoint - SoilTemperature) / (ConstGravityAcc * SoilTemperature);
  SoilWatSupercool = water.param.SoilMoistureSat(IndSoil) *
                     rpow(SoilWatPotFrz / water.param.SoilMatPotentialSat(IndSoil),
                          real(-1.0) / water.param.SoilExpCoeffB(IndSoil));
}

// ============================================================================
// SoilWaterSupercoolKoren99Mod.F90
// ============================================================================
CPPWRF_HD inline void SoilWaterSupercoolKoren99(NoahmpType& noahmp,
                                                int IndSoil,
                                                real& SoilWatSupercool,
                                                real SoilTemperature,
                                                real SoilMoisture,
                                                real SoilLiqWater) {
  auto& water = noahmp.water;

  constexpr real CK = real(8.0);
  constexpr real SoilExpBMax = real(5.5);
  constexpr real ErrorThr = real(0.005);

  real SoilExpB = water.param.SoilExpCoeffB(IndSoil);
  if (water.param.SoilExpCoeffB(IndSoil) > SoilExpBMax) SoilExpB = SoilExpBMax;
  int NumIter = 0;
  int IndCnt = 0;

  if (SoilTemperature > (ConstFreezePoint - real(1.0e-3))) {
    SoilWatSupercool = SoilMoisture;
  } else {
    if (CK != real(0.0)) {
      real SoilIce = SoilMoisture - SoilLiqWater;
      if (SoilIce > (SoilMoisture - real(0.02))) SoilIce = SoilMoisture - real(0.02);
      if (SoilIce < real(0.0)) SoilIce = real(0.0);
      while (NumIter < 10 && IndCnt == 0) {
        NumIter = NumIter + 1;
        real DF = rlog((water.param.SoilMatPotentialSat(IndSoil) * ConstGravityAcc / ConstLatHeatFusion) *
                       rpow(real(1.0) + CK * SoilIce, real(2.0)) *
                       rpow(water.param.SoilMoistureSat(IndSoil) / (SoilMoisture - SoilIce), SoilExpB)) -
                  rlog(-(SoilTemperature - ConstFreezePoint) / SoilTemperature);
        real Denom = real(2.0) * CK / (real(1.0) + CK * SoilIce) + SoilExpB / (SoilMoisture - SoilIce);
        real SoilIceTmp = SoilIce - DF / Denom;
        if (SoilIceTmp > (SoilMoisture - real(0.02))) SoilIceTmp = SoilMoisture - real(0.02);
        if (SoilIceTmp < real(0.0)) SoilIceTmp = real(0.0);
        real SoilIceChg = std::abs(SoilIceTmp - SoilIce);
        SoilIce = SoilIceTmp;
        if (SoilIceChg <= ErrorThr) {
          IndCnt = IndCnt + 1;
        }
      }
      SoilWatSupercool = SoilMoisture - SoilIce;
    }

    if (IndCnt == 0) {
      real FlerFac = rpow((ConstLatHeatFusion / (ConstGravityAcc * (-water.param.SoilMatPotentialSat(IndSoil)))) *
                              ((SoilTemperature - ConstFreezePoint) / SoilTemperature),
                          real(-1.0) / SoilExpB) *
                     water.param.SoilMoistureSat(IndSoil);
      if (FlerFac < real(0.02)) FlerFac = real(0.02);
      SoilWatSupercool = std::min(FlerFac, SoilMoisture);
    }
  }
}

// ============================================================================
// SoilSnowWaterPhaseChangeMod.F90
// ============================================================================
CPPWRF_HD inline void SoilSnowWaterPhaseChange(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& nml = noahmp.config.nmlist;
  auto& water = noahmp.water;
  auto& energy = noahmp.energy;

  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> EnergyRes;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> WaterPhaseChg;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MassWatTotInit;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MassWatIceInit;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MassWatLiqInit;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MassWatIceTmp;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> MassWatLiqTmp;
  EnergyRes.fill(real(0.0));
  WaterPhaseChg.fill(real(0.0));
  MassWatTotInit.fill(real(0.0));
  MassWatIceInit.fill(real(0.0));
  MassWatLiqInit.fill(real(0.0));
  MassWatIceTmp.fill(real(0.0));
  MassWatLiqTmp.fill(real(0.0));
  water.flux.MeltGroundSnow = real(0.0);
  water.state.PondSfcThinSnwMelt = real(0.0);
  real HeatLhTotPhsChg = real(0.0);

  for (int LoopInd = -domain.NumSnowLayerMax + 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    water.state.SoilSupercoolWater(LoopInd) = real(0.0);
  }

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
    MassWatIceTmp(LoopInd) = water.state.SnowIce(LoopInd);
    MassWatLiqTmp(LoopInd) = water.state.SnowLiqWater(LoopInd);
  }

  for (int LoopInd = 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    MassWatLiqTmp(LoopInd) = water.state.SoilLiqWater(LoopInd) * domain.ThicknessSnowSoilLayer(LoopInd) * real(1000.0);
    MassWatIceTmp(LoopInd) = (water.state.SoilMoisture(LoopInd) - water.state.SoilLiqWater(LoopInd)) *
                             domain.ThicknessSnowSoilLayer(LoopInd) * real(1000.0);
  }

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    water.state.IndexPhaseChange(LoopInd) = 0;
    EnergyRes(LoopInd) = real(0.0);
    WaterPhaseChg(LoopInd) = real(0.0);
    MassWatIceInit(LoopInd) = MassWatIceTmp(LoopInd);
    MassWatLiqInit(LoopInd) = MassWatLiqTmp(LoopInd);
    MassWatTotInit(LoopInd) = MassWatIceTmp(LoopInd) + MassWatLiqTmp(LoopInd);
  }

  if (domain.SurfaceType == 1) {
    for (int LoopInd = 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
      if (nml.OptSoilSupercoolWater == 1) {
        if (energy.state.TemperatureSoilSnow(LoopInd) < ConstFreezePoint) {
          SoilWaterSupercoolNiu06(noahmp, LoopInd, water.state.SoilSupercoolWater(LoopInd),
                                  energy.state.TemperatureSoilSnow(LoopInd));
          water.state.SoilSupercoolWater(LoopInd) =
              water.state.SoilSupercoolWater(LoopInd) * domain.ThicknessSnowSoilLayer(LoopInd) * real(1000.0);
        }
      }
      if (nml.OptSoilSupercoolWater == 2) {
        if (energy.state.TemperatureSoilSnow(LoopInd) < ConstFreezePoint) {
          SoilWaterSupercoolKoren99(noahmp, LoopInd, water.state.SoilSupercoolWater(LoopInd),
                                    energy.state.TemperatureSoilSnow(LoopInd),
                                    water.state.SoilMoisture(LoopInd), water.state.SoilLiqWater(LoopInd));
          water.state.SoilSupercoolWater(LoopInd) =
              water.state.SoilSupercoolWater(LoopInd) * domain.ThicknessSnowSoilLayer(LoopInd) * real(1000.0);
        }
      }
    }
  }

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    if (MassWatIceTmp(LoopInd) > real(0.0) && energy.state.TemperatureSoilSnow(LoopInd) >= ConstFreezePoint) {
      water.state.IndexPhaseChange(LoopInd) = 1;
    }
    if (MassWatLiqTmp(LoopInd) > water.state.SoilSupercoolWater(LoopInd) &&
        energy.state.TemperatureSoilSnow(LoopInd) < ConstFreezePoint) {
      water.state.IndexPhaseChange(LoopInd) = 2;
    }
    if (domain.NumSnowLayerNeg == 0 && water.state.SnowWaterEquiv > real(0.0) && LoopInd == 1) {
      if (energy.state.TemperatureSoilSnow(LoopInd) >= ConstFreezePoint) {
        water.state.IndexPhaseChange(LoopInd) = 1;
      }
    }
  }

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    if (water.state.IndexPhaseChange(LoopInd) > 0) {
      EnergyRes(LoopInd) =
          (energy.state.TemperatureSoilSnow(LoopInd) - ConstFreezePoint) / energy.state.PhaseChgFacSoilSnow(LoopInd);
      energy.state.TemperatureSoilSnow(LoopInd) = ConstFreezePoint;
    }
    if (water.state.IndexPhaseChange(LoopInd) == 1 && EnergyRes(LoopInd) < real(0.0)) {
      EnergyRes(LoopInd) = real(0.0);
      water.state.IndexPhaseChange(LoopInd) = 0;
    }
    if (water.state.IndexPhaseChange(LoopInd) == 2 && EnergyRes(LoopInd) > real(0.0)) {
      EnergyRes(LoopInd) = real(0.0);
      water.state.IndexPhaseChange(LoopInd) = 0;
    }
    WaterPhaseChg(LoopInd) = EnergyRes(LoopInd) * domain.MainTimeStep / ConstLatHeatFusion;
  }

  if (domain.NumSnowLayerNeg == 0 && water.state.SnowWaterEquiv > real(0.0) && WaterPhaseChg(1) > real(0.0)) {
    real SnowWaterPrev = water.state.SnowWaterEquiv;
    water.state.SnowWaterEquiv = std::max(real(0.0), SnowWaterPrev - WaterPhaseChg(1));
    real SnowWaterRatio = water.state.SnowWaterEquiv / SnowWaterPrev;
    water.state.SnowDepth = std::max(real(0.0), SnowWaterRatio * water.state.SnowDepth);
    water.state.SnowDepth =
        std::min(std::max(water.state.SnowDepth, water.state.SnowWaterEquiv / real(500.0)),
                 water.state.SnowWaterEquiv / real(50.0));
    real EnergyResLeft =
        EnergyRes(1) - ConstLatHeatFusion * (SnowWaterPrev - water.state.SnowWaterEquiv) / domain.MainTimeStep;
    if (EnergyResLeft > real(0.0)) {
      WaterPhaseChg(1) = EnergyResLeft * domain.MainTimeStep / ConstLatHeatFusion;
      EnergyRes(1) = EnergyResLeft;
    } else {
      WaterPhaseChg(1) = real(0.0);
      EnergyRes(1) = real(0.0);
    }
    water.flux.MeltGroundSnow = std::max(real(0.0), SnowWaterPrev - water.state.SnowWaterEquiv) / domain.MainTimeStep;
    HeatLhTotPhsChg = ConstLatHeatFusion * water.flux.MeltGroundSnow;
    water.state.PondSfcThinSnwMelt = SnowWaterPrev - water.state.SnowWaterEquiv;
  }

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    if (water.state.IndexPhaseChange(LoopInd) > 0 && std::abs(EnergyRes(LoopInd)) > real(0.0)) {
      real EnergyResLeft = real(0.0);
      if (WaterPhaseChg(LoopInd) > real(0.0)) {
        MassWatIceTmp(LoopInd) = std::max(real(0.0), MassWatIceInit(LoopInd) - WaterPhaseChg(LoopInd));
        EnergyResLeft = EnergyRes(LoopInd) -
                        ConstLatHeatFusion * (MassWatIceInit(LoopInd) - MassWatIceTmp(LoopInd)) / domain.MainTimeStep;
      } else if (WaterPhaseChg(LoopInd) < real(0.0)) {
        if (LoopInd <= 0) {
          MassWatIceTmp(LoopInd) =
              std::min(MassWatTotInit(LoopInd), MassWatIceInit(LoopInd) - WaterPhaseChg(LoopInd));
        } else {
          if (MassWatTotInit(LoopInd) < water.state.SoilSupercoolWater(LoopInd)) {
            MassWatIceTmp(LoopInd) = real(0.0);
          } else {
            MassWatIceTmp(LoopInd) =
                std::min(MassWatTotInit(LoopInd) - water.state.SoilSupercoolWater(LoopInd),
                         MassWatIceInit(LoopInd) - WaterPhaseChg(LoopInd));
            MassWatIceTmp(LoopInd) = std::max(MassWatIceTmp(LoopInd), real(0.0));
          }
        }
        EnergyResLeft = EnergyRes(LoopInd) -
                        ConstLatHeatFusion * (MassWatIceInit(LoopInd) - MassWatIceTmp(LoopInd)) / domain.MainTimeStep;
      }
      MassWatLiqTmp(LoopInd) = std::max(real(0.0), MassWatTotInit(LoopInd) - MassWatIceTmp(LoopInd));

      if (std::abs(EnergyResLeft) > real(0.0)) {
        energy.state.TemperatureSoilSnow(LoopInd) =
            energy.state.TemperatureSoilSnow(LoopInd) + energy.state.PhaseChgFacSoilSnow(LoopInd) * EnergyResLeft;
        if (LoopInd <= 0) {
          if ((MassWatLiqTmp(LoopInd) * MassWatIceTmp(LoopInd)) > real(0.0)) {
            energy.state.TemperatureSoilSnow(LoopInd) = ConstFreezePoint;
          }
          if (MassWatIceTmp(LoopInd) == real(0.0)) {
            energy.state.TemperatureSoilSnow(LoopInd) = ConstFreezePoint;
            EnergyRes(LoopInd + 1) = EnergyRes(LoopInd + 1) + EnergyResLeft;
            WaterPhaseChg(LoopInd + 1) = EnergyRes(LoopInd + 1) * domain.MainTimeStep / ConstLatHeatFusion;
          }
        }
      }
      HeatLhTotPhsChg =
          HeatLhTotPhsChg +
          ConstLatHeatFusion * (MassWatIceInit(LoopInd) - MassWatIceTmp(LoopInd)) / domain.MainTimeStep;
      (void)HeatLhTotPhsChg;
      if (LoopInd < 1) {
        water.flux.MeltGroundSnow =
            water.flux.MeltGroundSnow +
            std::max(real(0.0), MassWatIceInit(LoopInd) - MassWatIceTmp(LoopInd)) / domain.MainTimeStep;
      }
    }
  }

  for (int LoopInd = domain.NumSnowLayerNeg + 1; LoopInd <= 0; ++LoopInd) {
    water.state.SnowLiqWater(LoopInd) = MassWatLiqTmp(LoopInd);
    water.state.SnowIce(LoopInd) = MassWatIceTmp(LoopInd);
  }
  for (int LoopInd = 1; LoopInd <= domain.NumSoilLayer; ++LoopInd) {
    water.state.SoilLiqWater(LoopInd) =
        MassWatLiqTmp(LoopInd) / (real(1000.0) * domain.ThicknessSnowSoilLayer(LoopInd));
    water.state.SoilMoisture(LoopInd) =
        (MassWatLiqTmp(LoopInd) + MassWatIceTmp(LoopInd)) / (real(1000.0) * domain.ThicknessSnowSoilLayer(LoopInd));
  }
}

// ============================================================================
// EnergyMainMod.F90
// ============================================================================
CPPWRF_HD inline void EnergyMain(NoahmpType& noahmp) {
  auto& domain = noahmp.config.domain;
  auto& nml = noahmp.config.nmlist;
  auto& forcing = noahmp.forcing;
  auto& water = noahmp.water;
  auto& energy = noahmp.energy;
  auto& biochem = noahmp.biochem;

  energy.state.WindStressEwVeg = real(0.0);
  energy.state.WindStressNsVeg = real(0.0);
  energy.flux.RadLwNetCanopy = real(0.0);
  energy.flux.HeatSensibleCanopy = real(0.0);
  energy.flux.RadLwNetVegGrd = real(0.0);
  energy.flux.HeatSensibleVegGrd = real(0.0);
  energy.flux.HeatLatentVegGrd = real(0.0);
  energy.flux.HeatLatentCanEvap = real(0.0);
  energy.flux.HeatLatentCanTransp = real(0.0);
  energy.flux.HeatGroundVegGrd = real(0.0);
  biochem.flux.PhotosynLeafSunlit = real(0.0);
  biochem.flux.PhotosynLeafShade = real(0.0);
  energy.state.TemperatureAir2mVeg = real(0.0);
  energy.state.SpecHumidity2mVeg = real(0.0);
  energy.state.ExchCoeffShAbvCan = real(0.0);
  energy.state.ExchCoeffShLeaf = real(0.0);
  energy.state.ExchCoeffShUndCan = real(0.0);
  energy.state.ExchCoeffSh2mVeg = real(0.0);
  energy.flux.HeatPrecipAdvSfc = real(0.0);
  energy.flux.HeatCanStorageChg = real(0.0);

  energy.state.VegAreaIndEff = energy.state.LeafAreaIndEff + energy.state.StemAreaIndEff;
  bool FlagVegSfc = false;
  if (energy.state.VegAreaIndEff > real(0.0)) FlagVegSfc = true;

  if (nml.OptSnowCoverGround == 1) SnowCoverGroundNiu07(noahmp);
  if (nml.OptSnowCoverGround == 2) SnowCoverGroundAR25(noahmp);

  GroundRoughnessProperty(noahmp, FlagVegSfc);
  GroundThermalProperty(noahmp);
  SurfaceAlbedo(noahmp);
  SurfaceRadiation(noahmp);
  SurfaceEmissivity(noahmp);
  SoilWaterTranspiration(noahmp);
  ResistanceGroundEvaporation(noahmp);
  PsychrometricVariable(noahmp);

  if (FlagVegSfc && energy.state.VegFrac > real(0.0)) {
    energy.state.TemperatureGrdVeg = energy.state.TemperatureGrd;
    energy.state.ExchCoeffMomAbvCan = energy.state.ExchCoeffMomSfc;
    energy.state.ExchCoeffShAbvCan = energy.state.ExchCoeffShSfc;
    SurfaceEnergyFluxVegetated(noahmp);
  }

  energy.state.TemperatureGrdBare = energy.state.TemperatureGrd;
  energy.state.ExchCoeffMomBare = energy.state.ExchCoeffMomSfc;
  energy.state.ExchCoeffShBare = energy.state.ExchCoeffShSfc;
  SurfaceEnergyFluxBareGround(noahmp);

  if (FlagVegSfc && energy.state.VegFrac > real(0.0)) {
    energy.state.WindStressEwSfc =
        energy.state.VegFrac * energy.state.WindStressEwVeg +
        (real(1.0) - energy.state.VegFrac) * energy.state.WindStressEwBare;
    energy.state.WindStressNsSfc =
        energy.state.VegFrac * energy.state.WindStressNsVeg +
        (real(1.0) - energy.state.VegFrac) * energy.state.WindStressNsBare;
    energy.flux.RadLwNetSfc =
        energy.state.VegFrac * energy.flux.RadLwNetVegGrd +
        (real(1.0) - energy.state.VegFrac) * energy.flux.RadLwNetBareGrd + energy.flux.RadLwNetCanopy;
    energy.flux.HeatSensibleSfc =
        energy.state.VegFrac * energy.flux.HeatSensibleVegGrd +
        (real(1.0) - energy.state.VegFrac) * energy.flux.HeatSensibleBareGrd + energy.flux.HeatSensibleCanopy;
    energy.flux.HeatLatentGrd =
        energy.state.VegFrac * energy.flux.HeatLatentVegGrd +
        (real(1.0) - energy.state.VegFrac) * energy.flux.HeatLatentBareGrd;
    energy.flux.HeatGroundTot =
        energy.state.VegFrac * energy.flux.HeatGroundVegGrd +
        (real(1.0) - energy.state.VegFrac) * energy.flux.HeatGroundBareGrd;
    energy.flux.HeatLatentCanopy = energy.flux.HeatLatentCanEvap;
    energy.flux.HeatLatentTransp = energy.flux.HeatLatentCanTransp;
    energy.flux.HeatPrecipAdvSfc =
        energy.state.VegFrac * energy.flux.HeatPrecipAdvVegGrd +
        (real(1.0) - energy.state.VegFrac) * energy.flux.HeatPrecipAdvBareGrd + energy.flux.HeatPrecipAdvCanopy;
    energy.state.TemperatureGrd =
        energy.state.VegFrac * energy.state.TemperatureGrdVeg +
        (real(1.0) - energy.state.VegFrac) * energy.state.TemperatureGrdBare;
    energy.state.TemperatureAir2m =
        energy.state.VegFrac * energy.state.TemperatureAir2mVeg +
        (real(1.0) - energy.state.VegFrac) * energy.state.TemperatureAir2mBare;
    energy.state.TemperatureSfc =
        energy.state.VegFrac * energy.state.TemperatureCanopy +
        (real(1.0) - energy.state.VegFrac) * energy.state.TemperatureGrdBare;
    energy.state.ExchCoeffMomSfc =
        energy.state.VegFrac * energy.state.ExchCoeffMomAbvCan +
        (real(1.0) - energy.state.VegFrac) * energy.state.ExchCoeffMomBare;
    energy.state.ExchCoeffShSfc =
        energy.state.VegFrac * energy.state.ExchCoeffShAbvCan +
        (real(1.0) - energy.state.VegFrac) * energy.state.ExchCoeffShBare;
    energy.state.SpecHumidity2m =
        energy.state.VegFrac * energy.state.SpecHumidity2mVeg +
        (real(1.0) - energy.state.VegFrac) * energy.state.SpecHumidity2mBare;
    energy.state.SpecHumiditySfcMean =
        energy.state.VegFrac *
            (energy.state.PressureVaporCanAir * real(0.622) /
             (forcing.PressureAirRefHeight - real(0.378) * energy.state.PressureVaporCanAir)) +
        (real(1.0) - energy.state.VegFrac) * energy.state.SpecHumiditySfc;
    energy.state.RoughLenMomSfcToAtm = energy.state.RoughLenMomSfc;
  } else {
    energy.state.WindStressEwSfc = energy.state.WindStressEwBare;
    energy.state.WindStressNsSfc = energy.state.WindStressNsBare;
    energy.flux.RadLwNetSfc = energy.flux.RadLwNetBareGrd;
    energy.flux.HeatSensibleSfc = energy.flux.HeatSensibleBareGrd;
    energy.flux.HeatLatentGrd = energy.flux.HeatLatentBareGrd;
    energy.flux.HeatGroundTot = energy.flux.HeatGroundBareGrd;
    energy.state.TemperatureGrd = energy.state.TemperatureGrdBare;
    energy.state.TemperatureAir2m = energy.state.TemperatureAir2mBare;
    energy.flux.HeatLatentCanopy = real(0.0);
    energy.flux.HeatLatentTransp = real(0.0);
    energy.flux.HeatPrecipAdvSfc = energy.flux.HeatPrecipAdvBareGrd;
    energy.state.TemperatureSfc = energy.state.TemperatureGrd;
    energy.state.ExchCoeffMomSfc = energy.state.ExchCoeffMomBare;
    energy.state.ExchCoeffShSfc = energy.state.ExchCoeffShBare;
    energy.state.SpecHumiditySfcMean = energy.state.SpecHumiditySfc;
    energy.state.SpecHumidity2m = energy.state.SpecHumidity2mBare;
    energy.state.ResistanceStomataSunlit = real(0.0);
    energy.state.ResistanceStomataShade = real(0.0);
    energy.state.TemperatureGrdVeg = energy.state.TemperatureGrdBare;
    energy.state.ExchCoeffShAbvCan = energy.state.ExchCoeffShBare;
    energy.state.RoughLenMomSfcToAtm = energy.state.RoughLenMomGrd;
  }

  energy.flux.RadLwEmitSfc = forcing.RadLwDownRefHeight + energy.flux.RadLwNetSfc;
  energy.state.TemperatureRadSfc =
      rpow((energy.flux.RadLwEmitSfc - (real(1.0) - energy.state.EmissivitySfc) * forcing.RadLwDownRefHeight) /
               (energy.state.EmissivitySfc * ConstStefanBoltzmann),
           real(0.25));

  energy.flux.RadPhotoActAbsCan =
      energy.flux.RadPhotoActAbsSunlit * energy.state.LeafAreaIndSunlit +
      energy.flux.RadPhotoActAbsShade * energy.state.LeafAreaIndShade;
  biochem.flux.PhotosynTotal =
      biochem.flux.PhotosynLeafSunlit * energy.state.LeafAreaIndSunlit +
      biochem.flux.PhotosynLeafShade * energy.state.LeafAreaIndShade;

  energy.flux.HeatFromSoilBot = real(0.0);
  energy.flux.HeatGroundTotAcc = energy.flux.HeatGroundTotAcc + energy.flux.HeatGroundTot;
  if (domain.FlagSoilProcess) {
    energy.flux.HeatGroundTotMean = energy.flux.HeatGroundTotAcc / domain.NumSoilTimeStep;
    SoilSnowTemperatureMain(noahmp);
  }

  if (nml.OptSnowSoilTempTime == 2) {
    if (water.state.SnowDepth > real(0.05) && energy.state.TemperatureGrd > ConstFreezePoint) {
      energy.state.TemperatureGrdVeg = ConstFreezePoint;
      energy.state.TemperatureGrdBare = ConstFreezePoint;
      if (FlagVegSfc && energy.state.VegFrac > real(0.0)) {
        energy.state.TemperatureGrd =
            energy.state.VegFrac * energy.state.TemperatureGrdVeg +
            (real(1.0) - energy.state.VegFrac) * energy.state.TemperatureGrdBare;
        energy.state.TemperatureSfc =
            energy.state.VegFrac * energy.state.TemperatureCanopy +
            (real(1.0) - energy.state.VegFrac) * energy.state.TemperatureGrdBare;
      } else {
        energy.state.TemperatureGrd = energy.state.TemperatureGrdBare;
        energy.state.TemperatureSfc = energy.state.TemperatureGrdBare;
      }
    }
  }

  SoilSnowWaterPhaseChange(noahmp);

  if (domain.FlagCropland && water.state.IrrigationFracGrid >= water.param.IrriFracThreshold) {
    energy.flux.HeatSensibleSfc = energy.flux.HeatSensibleSfc - energy.flux.HeatLatentIrriEvap;
  }

  if (forcing.RadSwDownRefHeight > real(0.0)) {
    energy.state.AlbedoSfc = energy.flux.RadSwReflSfc / forcing.RadSwDownRefHeight;
  } else {
    energy.state.AlbedoSfc = undefined_real;
  }
}

}  // namespace cppwrf::phys::noahmp
