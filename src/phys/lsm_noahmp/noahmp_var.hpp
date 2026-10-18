#pragma once

#include <algorithm>
#include <cmath>
#include "core/backend.hpp"
#include "core/precision.hpp"
#include "phys/lsm_noahmp/noahmp_tables.hpp"

namespace cppwrf::phys::noahmp {

CPPWRF_HD inline float rexp(float a) {
#if defined(__CUDA_ARCH__)
  return expf(a);
#else
  return std::exp(a);
#endif
}
CPPWRF_HD inline double rexp(double a) {
#if defined(__CUDA_ARCH__)
  return ::exp(a);
#else
  return std::exp(a);
#endif
}

CPPWRF_HD inline float rpow(float a, float b) {
#if defined(__CUDA_ARCH__)
  return powf(a, b);
#else
  return std::pow(a, b);
#endif
}
CPPWRF_HD inline double rpow(double a, double b) {
#if defined(__CUDA_ARCH__)
  return ::pow(a, b);
#else
  return std::pow(a, b);
#endif
}

CPPWRF_HD inline float rlog(float a) {
#if defined(__CUDA_ARCH__)
  return logf(a);
#else
  return std::log(a);
#endif
}
CPPWRF_HD inline double rlog(double a) {
#if defined(__CUDA_ARCH__)
  return ::log(a);
#else
  return std::log(a);
#endif
}

// Physical constants from ConstantDefineMod.F90 (in working precision `real`)
inline constexpr real ConstGravityAcc          = real(9.80616);    // acceleration due to gravity [m/s2]
inline constexpr real ConstStefanBoltzmann     = real(5.67e-08);   // Stefan-Boltzmann constant [w/m2/k4]
inline constexpr real ConstBoltzmann           = real(5.67e-08);   // alias for Stefan-Boltzmann constant [w/m2/k4]
inline constexpr real ConstVonKarman           = real(0.40);       // Von Karman constant
inline constexpr real ConstFreezePoint         = real(273.16);     // freezing/melting point [K]
inline constexpr real ConstLatHeatSublim       = real(2.8440e06);  // latent heat of sublimation [j/kg]
inline constexpr real ConstLatHeatEvap         = real(2.5104e06);  // latent heat of vaporization [j/kg]
inline constexpr real ConstLatHeatVapor        = real(2.5104e06);  // alias for latent heat of vaporization [j/kg]
inline constexpr real ConstLatHeatFusion       = real(0.3336e06);  // latent heat of fusion [j/kg]
inline constexpr real ConstHeatCapacWater      = real(4.188e06);   // specific heat capacity of water [j/m3/k]
inline constexpr real ConstHeatCapacIce        = real(2.094e06);   // specific heat capacity of ice [j/m3/k]
inline constexpr real ConstHeatCapacAir        = real(1004.64);    // heat capacity dry air at const pres [j/kg/k]
inline constexpr real ConstThermConductWater   = real(0.57);       // thermal conductivity of liquid water [w/m/k]
inline constexpr real ConstThermalCondWat      = real(0.57);       // alias
inline constexpr real ConstThermConductIce     = real(2.2);        // thermal conductivity of ice [w/m/k]
inline constexpr real ConstThermalCondIce      = real(2.2);        // alias
inline constexpr real ConstThermConductAir     = real(0.023);      // thermal conductivity of air [w/m/k]
inline constexpr real ConstThermalCondAir      = real(0.023);      // alias
inline constexpr real ConstThermConductQuartz  = real(7.7);        // thermal conductivity for quartz [w/m/k]
inline constexpr real ConstThermConductSoilOth = real(2.0);        // thermal conductivity for other soil [w/m/k]
inline constexpr real ConstGasDryAir           = real(287.04);     // gas constant for dry air [j/kg/k]
inline constexpr real ConstGasWaterVapor       = real(461.269);    // gas constant for water vapor [j/kg/k]
inline constexpr real ConstDensityWater        = real(1000.0);     // density of liquid water [kg/m3]
inline constexpr real ConstDensityIce          = real(917.0);      // density of ice [kg/m3]
inline constexpr real ConstPI                  = real(3.14159265); // pi (exact Fortran literal in ConstantDefineMod.F90)
inline constexpr real ConstPi                  = real(3.14159265); // alias
inline constexpr real ConstDegToRad            = ConstPI / real(180.0); // degree to radian
inline constexpr real ConstDensityGraupel      = real(500.0);      // bulk density of graupel [kg/m3]
inline constexpr real ConstDensityHail         = real(917.0);      // bulk density of hail [kg/m3]

inline constexpr int  undefined_int            = -9999;
inline constexpr real undefined_real           = real(-9999.0);

// Fortran NINT(x)
CPPWRF_HD inline int nint_f(real x) {
  return static_cast<int>(x + (x >= real(0.0) ? real(0.5) : real(-0.5)));
}

// Fortran SIGN(a, b)
CPPWRF_HD inline real sign_f(real a, real b) {
  return std::copysign(std::abs(a), b);
}

// 1D fixed-bounds array with Fortran-style indexing [Low..High] inclusive
template <typename T, int Low, int High>
struct FArray1D {
  T data_[High - Low + 1];

  CPPWRF_HD constexpr T& operator()(int i) { return data_[i - Low]; }
  CPPWRF_HD constexpr const T& operator()(int i) const { return data_[i - Low]; }

  CPPWRF_HD constexpr void fill(T val) {
    for (int i = 0; i < High - Low + 1; ++i) data_[i] = val;
  }
};

// 2D fixed-bounds array with Fortran-style indexing [Low1..High1, Low2..High2] inclusive
template <typename T, int Low1, int High1, int Low2, int High2>
struct FArray2D {
  static constexpr int N1 = High1 - Low1 + 1;
  static constexpr int N2 = High2 - Low2 + 1;
  T data_[N1 * N2];

  CPPWRF_HD constexpr T& operator()(int i, int j) { return data_[(i - Low1) * N2 + (j - Low2)]; }
  CPPWRF_HD constexpr const T& operator()(int i, int j) const { return data_[(i - Low1) * N2 + (j - Low2)]; }

  CPPWRF_HD constexpr void fill(T val) {
    for (int k = 0; k < N1 * N2; ++k) data_[k] = val;
  }
};

// Tridiagonal matrix solver from MatrixSolverTriDiagonalMod.F90
template <typename ArrP, typename ArrA, typename ArrB, typename ArrC, typename ArrD, typename ArrDelta>
CPPWRF_HD inline void MatrixSolverTriDiagonal(ArrP& P,
                                              const ArrA& A,
                                              const ArrB& B,
                                              ArrC& C,
                                              const ArrD& D,
                                              ArrDelta& Delta,
                                              int IndTopLayer,
                                              int NumSoilLayer,
                                              int NumSnowLayerMax) {
  (void)NumSnowLayerMax;
  // initialize matrix coef for the lowest soil layer
  C(NumSoilLayer) = real(0.0);
  P(IndTopLayer) = -C(IndTopLayer) / B(IndTopLayer);

  // solve the coefs for the 1st soil layer
  Delta(IndTopLayer) = D(IndTopLayer) / B(IndTopLayer);

  // solve the coefs for soil layers 2 to NumSoilLayer
  for (int K = IndTopLayer + 1; K <= NumSoilLayer; ++K) {
    P(K) = -C(K) * (real(1.0) / (B(K) + A(K) * P(K - 1)));
    Delta(K) = (D(K) - A(K) * Delta(K - 1)) * (real(1.0) / (B(K) + A(K) * P(K - 1)));
  }

  // set p to delta for lowest soil layer
  P(NumSoilLayer) = Delta(NumSoilLayer);

  // adjust p for soil layers 2 to NumSoilLayer
  for (int K = IndTopLayer + 1; K <= NumSoilLayer; ++K) {
    int KK = NumSoilLayer - K + (IndTopLayer - 1) + 1;
    P(KK) = P(KK) * P(KK + 1) + Delta(KK);
  }
}


// === ConfigVarType ===
struct NamelistType {
  int OptDynamicVeg              = undefined_int;
  int OptRainSnowPartition       = undefined_int;
  int OptSoilWaterTranspiration  = undefined_int;
  int OptGroundResistanceEvap    = undefined_int;
  int OptSurfaceDrag             = undefined_int;
  int OptStomataResistance       = undefined_int;
  int OptSnowAlbedo              = undefined_int;
  int OptCanopyRadiationTransfer = undefined_int;
  int OptSnowSoilTempTime        = undefined_int;
  int OptSnowThermConduct        = undefined_int;
  int OptSoilTemperatureBottom   = undefined_int;
  int OptSoilSupercoolWater      = undefined_int;
  int OptRunoffSurface           = undefined_int;
  int OptRunoffSubsurface        = undefined_int;
  int OptSoilPermeabilityFrozen  = undefined_int;
  int OptDynVicInfiltration      = undefined_int;
  int OptTileDrainage            = undefined_int;
  int OptIrrigation              = undefined_int;
  int OptIrrigationMethod        = undefined_int;
  int OptCropModel               = undefined_int;
  int OptSoilProperty            = undefined_int;
  int OptPedotransfer            = undefined_int;
  int OptGlacierTreatment        = undefined_int;
  int OptSnowCompaction          = undefined_int;
  int OptSnowCoverGround         = undefined_int;
  int OptWetlandModel            = undefined_int;
  int OptSnicarSnowShape         = undefined_int;
  int OptSnicarRTSolver          = undefined_int;
  int OptSnicarBandNum           = undefined_int;
  int OptSnicarSolarSpec         = undefined_int;
  int OptSnicarSnwOptic          = undefined_int;
  int OptSnicarDustOptic         = undefined_int;
  bool FlagSnicarSnowBCIntmix     = false;
  bool FlagSnicarSnowDustIntmix   = false;
  bool FlagSnicarUseAerosol       = false;
  bool FlagSnicarUseOC            = false;
  bool FlagSnicarAerosolReadTable = false;
};

struct DomainType {
  int IsModisLanduse              = 1;  // 1 if MODIFIED_IGBP_MODIS_NOAH, 0 if USGS
  bool FlagUrban                  = false;
  bool FlagCropland               = false;
  bool FlagWetland                = false;
  bool FlagDynamicCrop            = false;
  bool FlagDynamicVeg             = false;
  bool FlagSoilProcess            = false;
  int GridIndexI                  = undefined_int;
  int GridIndexJ                  = undefined_int;
  int VegType                     = undefined_int;
  int CropType                    = undefined_int;
  int NumSoilLayer                = undefined_int;
  int NumSnowLayerMax             = undefined_int;
  int NumSnowLayerNeg             = undefined_int;
  int SurfaceType                 = undefined_int;
  int NumSwRadBand                = undefined_int;
  int SoilColor                   = undefined_int;
  int IndicatorIceSfc             = undefined_int;
  int IndexWaterPoint             = undefined_int;
  int IndexBarrenPoint            = undefined_int;
  int IndexIcePoint               = undefined_int;
  int IndexCropPoint              = undefined_int;
  int IndexEBLForest              = undefined_int;
  int NumCropGrowStage            = undefined_int;
  int NumDayInYear                = undefined_int;
  int RunoffSlopeType             = undefined_int;
  int NumSoilTimeStep             = undefined_int;
  int NumTempSnwAgeSnicar         = undefined_int;
  int NumTempGradSnwAgeSnicar     = undefined_int;
  int NumDensitySnwAgeSnicar      = undefined_int;
  int NumSnicarRadBand            = undefined_int;
  int NumRadiusSnwMieSnicar       = undefined_int;
  real MainTimeStep               = undefined_real;
  real SoilTimeStep               = undefined_real;
  real GridSize                   = undefined_real;
  real DayJulianInYear            = undefined_real;
  real CosSolarZenithAngle        = undefined_real;
  real RefHeightAboveSfc          = undefined_real;
  real ThicknessAtmosBotLayer     = undefined_real;
  real Latitude                   = undefined_real;
  real DepthSoilTempBottom        = undefined_real;

  FArray1D<int, 1, kMaxSoilLayers> SoilType;
  FArray1D<real, 1, kMaxSoilLayers> DepthSoilLayer;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> ThicknessSnowSoilLayer;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> DepthSnowSoilLayer;
  FArray1D<real, 1, kMaxSoilLayers> ThicknessSoilLayer;
};

struct ConfigType {
  NamelistType nmlist;
  DomainType domain;
};

// === ForcingVarType ===
struct ForcingType {
  real SpecHumidityRefHeight   = undefined_real;
  real TemperatureAirRefHeight = undefined_real;
  real WindEastwardRefHeight   = undefined_real;
  real WindNorthwardRefHeight  = undefined_real;
  real RadSwDownRefHeight      = undefined_real;
  real RadLwDownRefHeight      = undefined_real;
  real PressureAirRefHeight    = undefined_real;
  real PressureAirSurface      = undefined_real;
  real PrecipConvRefHeight     = undefined_real;
  real PrecipNonConvRefHeight  = undefined_real;
  real PrecipShConvRefHeight   = undefined_real;
  real PrecipSnowRefHeight     = undefined_real;
  real PrecipGraupelRefHeight  = undefined_real;
  real PrecipHailRefHeight     = undefined_real;
  real TemperatureSoilBottom   = undefined_real;
  real DepBChydropho           = undefined_real;
  real DepBChydrophi           = undefined_real;
  real DepOChydropho           = undefined_real;
  real DepOChydrophi           = undefined_real;
  real DepDust1                = undefined_real;
  real DepDust2                = undefined_real;
  real DepDust3                = undefined_real;
  real DepDust4                = undefined_real;
  real DepDust5                = undefined_real;
  real RadSwVisFrac            = undefined_real;
  real RadSwDirFrac            = undefined_real;
};

// === EnergyVarType ===
struct EnergyFluxType {
  real HeatLatentCanopy        = undefined_real;
  real HeatLatentTransp        = undefined_real;
  real HeatLatentGrd           = undefined_real;
  real HeatLatentIrriEvap      = real(0.0);
  real HeatPrecipAdvCanopy     = undefined_real;
  real HeatPrecipAdvVegGrd     = undefined_real;
  real HeatPrecipAdvBareGrd    = undefined_real;
  real HeatPrecipAdvSfc        = undefined_real;
  real HeatSensibleCanopy      = undefined_real;
  real HeatLatentCanEvap       = undefined_real;
  real HeatSensibleVegGrd      = undefined_real;
  real HeatSensibleSfc         = undefined_real;
  real HeatLatentVegGrd        = undefined_real;
  real HeatLatentCanTransp     = undefined_real;
  real HeatGroundVegGrd        = undefined_real;
  real HeatSensibleBareGrd     = undefined_real;
  real HeatLatentBareGrd       = undefined_real;
  real HeatGroundBareGrd       = undefined_real;
  real HeatGroundTot           = undefined_real;
  real HeatGroundTotMean       = undefined_real;
  real HeatFromSoilBot         = undefined_real;
  real HeatCanStorageChg       = undefined_real;
  real HeatGroundTotAcc        = undefined_real;
  real RadPhotoActAbsSunlit    = undefined_real;
  real RadPhotoActAbsShade     = undefined_real;
  real RadSwAbsVeg             = undefined_real;
  real RadSwAbsGrd             = undefined_real;
  real RadSwAbsSfc             = undefined_real;
  real RadSwReflSfc            = undefined_real;
  real RadSwReflVeg            = undefined_real;
  real RadSwReflGrd            = undefined_real;
  real RadLwNetCanopy          = undefined_real;
  real RadLwNetSfc             = undefined_real;
  real RadPhotoActAbsCan       = undefined_real;
  real RadLwEmitSfc            = undefined_real;
  real RadLwNetVegGrd          = undefined_real;
  real RadLwNetBareGrd         = undefined_real;

  FArray1D<real, 1, kNumRadBands> RadSwAbsVegDir;
  FArray1D<real, 1, kNumRadBands> RadSwAbsVegDif;
  FArray1D<real, 1, kNumRadBands> RadSwDirTranGrdDir;
  FArray1D<real, 1, kNumRadBands> RadSwDirTranGrdDif;
  FArray1D<real, 1, kNumRadBands> RadSwDifTranGrdDir;
  FArray1D<real, 1, kNumRadBands> RadSwDifTranGrdDif;
  FArray1D<real, 1, kNumRadBands> RadSwReflVegDir;
  FArray1D<real, 1, kNumRadBands> RadSwReflVegDif;
  FArray1D<real, 1, kNumRadBands> RadSwReflGrdDir;
  FArray1D<real, 1, kNumRadBands> RadSwReflGrdDif;
  FArray1D<real, 1, kNumRadBands> RadSwDownDir;
  FArray1D<real, 1, kNumRadBands> RadSwDownDif;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> RadSwPenetrateGrd;
  FArray1D<real, -kMaxSnowLayers + 1, 1> RadSwAbsSnowSoilLayer;
  FArray2D<real, -kMaxSnowLayers + 1, 1, 1, kNumRadBands> FracRadSwAbsSnowDir;
  FArray2D<real, -kMaxSnowLayers + 1, 1, 1, kNumRadBands> FracRadSwAbsSnowDif;
};

struct EnergyStateType {
  bool FlagFrozenCanopy          = false;
  bool FlagFrozenGround          = false;
  real LeafAreaIndEff            = undefined_real;
  real StemAreaIndEff            = undefined_real;
  real LeafAreaIndex             = undefined_real;
  real StemAreaIndex             = undefined_real;
  real VegAreaIndEff             = undefined_real;
  real VegFrac                   = undefined_real;
  real TemperatureGrd            = undefined_real;
  real TemperatureCanopy         = undefined_real;
  real TemperatureSfc            = undefined_real;
  real TemperatureRootZone       = undefined_real;
  real PressureVaporRefHeight    = undefined_real;
  real SnowAgeFac                = undefined_real;
  real SnowAgeNondim             = undefined_real;
  real AlbedoSnowPrev            = undefined_real;
  real VegAreaProjDir            = undefined_real;
  real GapBtwCanopy              = undefined_real;
  real GapInCanopy               = undefined_real;
  real GapCanopyDif              = undefined_real;
  real GapCanopyDir              = undefined_real;
  real CanopySunlitFrac          = undefined_real;
  real CanopyShadeFrac           = undefined_real;
  real LeafAreaIndSunlit         = undefined_real;
  real LeafAreaIndShade          = undefined_real;
  real VapPresSatCanopy          = undefined_real;
  real VapPresSatGrdVeg          = undefined_real;
  real VapPresSatGrdBare         = undefined_real;
  real VapPresSatCanTempD        = undefined_real;
  real VapPresSatGrdVegTempD     = undefined_real;
  real VapPresSatGrdBareTempD    = undefined_real;
  real PressureVaporCanAir       = undefined_real;
  real PressureAtmosCO2          = undefined_real;
  real PressureAtmosO2           = undefined_real;
  real ResistanceStomataSunlit   = undefined_real;
  real ResistanceStomataShade    = undefined_real;
  real DensityAirRefHeight       = undefined_real;
  real TemperatureCanopyAir      = undefined_real;
  real ZeroPlaneDispSfc          = undefined_real;
  real ZeroPlaneDispGrd          = undefined_real;
  real RoughLenMomGrd            = undefined_real;
  real RoughLenMomSfc            = undefined_real;
  real RoughLenShCanopy          = undefined_real;
  real RoughLenShVegGrd          = undefined_real;
  real RoughLenShBareGrd         = undefined_real;
  real CanopyHeight              = undefined_real;
  real WindSpdCanopyTop          = undefined_real;
  real FrictionVelVeg            = undefined_real;
  real FrictionVelBare           = undefined_real;
  real WindExtCoeffCanopy        = undefined_real;
  real MoStabParaUndCan          = undefined_real;
  real MoStabParaAbvCan          = undefined_real;
  real MoStabParaBare            = undefined_real;
  real MoStabParaVeg2m           = undefined_real;
  real MoStabParaBare2m          = undefined_real;
  real MoLengthUndCan            = undefined_real;
  real MoLengthAbvCan            = undefined_real;
  real MoLengthBare              = undefined_real;
  real MoStabCorrShUndCan        = undefined_real;
  real MoStabCorrMomAbvCan       = undefined_real;
  real MoStabCorrShAbvCan        = undefined_real;
  real MoStabCorrMomVeg2m        = undefined_real;
  real MoStabCorrShVeg2m         = undefined_real;
  real MoStabCorrShBare          = undefined_real;
  real MoStabCorrMomBare         = undefined_real;
  real MoStabCorrMomBare2m       = undefined_real;
  real MoStabCorrShBare2m        = undefined_real;
  real ExchCoeffMomSfc           = undefined_real;
  real ExchCoeffMomAbvCan        = undefined_real;
  real ExchCoeffMomBare          = undefined_real;
  real ExchCoeffShSfc            = undefined_real;
  real ExchCoeffShAbvCan         = undefined_real;
  real ExchCoeffShBare           = undefined_real;
  real ExchCoeffSh2mVegMo        = undefined_real;
  real ExchCoeffSh2mBareMo       = undefined_real;
  real ExchCoeffSh2mVeg          = undefined_real;
  real ExchCoeffLhAbvCan         = undefined_real;
  real ExchCoeffLhTransp         = undefined_real;
  real ExchCoeffLhEvap           = undefined_real;
  real ExchCoeffLhUndCan         = undefined_real;
  real ResistanceMomUndCan       = undefined_real;
  real ResistanceShUndCan        = undefined_real;
  real ResistanceLhUndCan        = undefined_real;
  real ResistanceMomAbvCan       = undefined_real;
  real ResistanceShAbvCan        = undefined_real;
  real ResistanceLhAbvCan        = undefined_real;
  real ResistanceMomBareGrd      = undefined_real;
  real ResistanceShBareGrd       = undefined_real;
  real ResistanceLhBareGrd       = undefined_real;
  real ResistanceLeafBoundary    = undefined_real;
  real TemperaturePotRefHeight   = undefined_real;
  real WindSpdRefHeight          = undefined_real;
  real FrictionVelVertVeg        = undefined_real;
  real FrictionVelVertBare       = undefined_real;
  real EmissivityVeg             = undefined_real;
  real EmissivityGrd             = undefined_real;
  real ResistanceGrdEvap         = undefined_real;
  real PsychConstCanopy          = undefined_real;
  real LatHeatVapCanopy          = undefined_real;
  real PsychConstGrd             = undefined_real;
  real LatHeatVapGrd             = undefined_real;
  real RelHumidityGrd            = undefined_real;
  real SpecHumiditySfc           = undefined_real;
  real SpecHumiditySfcMean       = undefined_real;
  real SpecHumidity2mVeg         = undefined_real;
  real SpecHumidity2mBare        = undefined_real;
  real SpecHumidity2m            = undefined_real;
  real TemperatureGrdVeg         = undefined_real;
  real TemperatureGrdBare        = undefined_real;
  real WindStressEwVeg           = undefined_real;
  real WindStressNsVeg           = undefined_real;
  real WindStressEwBare          = undefined_real;
  real WindStressNsBare          = undefined_real;
  real WindStressEwSfc           = undefined_real;
  real WindStressNsSfc           = undefined_real;
  real TemperatureAir2mVeg       = undefined_real;
  real TemperatureAir2mBare      = undefined_real;
  real TemperatureAir2m          = undefined_real;
  real ExchCoeffShLeaf           = undefined_real;
  real ExchCoeffShUndCan         = undefined_real;
  real ExchCoeffSh2mBare         = undefined_real;
  real RefHeightAboveGrd         = undefined_real;
  real CanopyFracSnowBury        = undefined_real;
  real DepthSoilTempBotToSno     = undefined_real;
  real RoughLenMomSfcToAtm       = undefined_real;
  real TemperatureRadSfc         = undefined_real;
  real EmissivitySfc             = undefined_real;
  real AlbedoSfc                 = undefined_real;
  real EnergyBalanceError        = undefined_real;
  real RadSwBalanceError         = undefined_real;

  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> TemperatureSoilSnow;
  FArray1D<real, -kMaxSnowLayers + 1, 0> HeatCapacVolSnow;
  FArray1D<real, -kMaxSnowLayers + 1, 0> ThermConductSnow;
  FArray1D<real, 1, kMaxSoilLayers> HeatCapacVolSoil;
  FArray1D<real, 1, kMaxSoilLayers> ThermConductSoil;
  FArray1D<real, 1, kMaxSoilLayers> HeatCapacGlaIce;
  FArray1D<real, 1, kMaxSoilLayers> ThermConductGlaIce;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> ThermConductSoilSnow;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> HeatCapacSoilSnow;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> PhaseChgFacSoilSnow;
  FArray1D<real, 1, kNumRadBands> AlbedoSnowDir;
  FArray1D<real, 1, kNumRadBands> AlbedoSnowDif;
  FArray1D<real, 1, kNumRadBands> AlbedoSoilDir;
  FArray1D<real, 1, kNumRadBands> AlbedoSoilDif;
  FArray1D<real, 1, kNumRadBands> AlbedoGrdDir;
  FArray1D<real, 1, kNumRadBands> AlbedoGrdDif;
  FArray1D<real, 1, kNumRadBands> ReflectanceVeg;
  FArray1D<real, 1, kNumRadBands> TransmittanceVeg;
  FArray1D<real, 1, kNumRadBands> AlbedoSfcDir;
  FArray1D<real, 1, kNumRadBands> AlbedoSfcDif;
};

struct EnergyParamType {
  real TreeCrownRadius       = undefined_real;
  real HeightCanopyTop       = undefined_real;
  real HeightCanopyBot       = undefined_real;
  real RoughLenMomVeg        = undefined_real;
  real TreeDensity           = undefined_real;
  real CanopyOrientIndex     = undefined_real;
  real UpscatterCoeffSnowDir = undefined_real;
  real UpscatterCoeffSnowDif = undefined_real;
  real SoilHeatCapacity      = undefined_real;
  real SnowAgeFacBats        = undefined_real;
  real SnowGrowVapFacBats    = undefined_real;
  real SnowSootFacBats       = undefined_real;
  real SnowGrowFrzFacBats    = undefined_real;
  real SolarZenithAdjBats    = undefined_real;
  real FreshSnoAlbVisBats    = undefined_real;
  real FreshSnoAlbNirBats    = undefined_real;
  real SnoAgeFacDifVisBats   = undefined_real;
  real SnoAgeFacDifNirBats   = undefined_real;
  real SzaFacDirVisBats      = undefined_real;
  real SzaFacDirNirBats      = undefined_real;
  real SnowAlbRefClass       = undefined_real;
  real SnowAgeFacClass       = undefined_real;
  real SnowAlbFreshClass     = undefined_real;
  real ConductanceLeafMin    = undefined_real;
  real Co2MmConst25C         = undefined_real;
  real O2MmConst25C          = undefined_real;
  real Co2MmConstQ10         = undefined_real;
  real O2MmConstQ10          = undefined_real;
  real RadiationStressFac    = undefined_real;
  real ResistanceStomataMin  = undefined_real;
  real ResistanceStomataMax  = undefined_real;
  real AirTempOptimTransp    = undefined_real;
  real VaporPresDeficitFac   = undefined_real;
  real LeafDimLength         = undefined_real;
  real ZilitinkevichCoeff    = undefined_real;
  real EmissivitySnow        = undefined_real;
  real CanopyWindExtFac      = undefined_real;
  real RoughLenMomSnow       = undefined_real;
  real RoughLenMomSoil       = undefined_real;
  real RoughLenMomLake       = undefined_real;
  real EmissivityIceSfc      = undefined_real;
  real ResistanceSoilExp     = undefined_real;
  real ResistanceSnowSfc     = undefined_real;
  real VegFracGreen          = undefined_real;
  real VegFracAnnMax         = undefined_real;
  real HeatCapacCanFac       = undefined_real;

  FArray1D<real, 1, 12> LeafAreaIndexMon;
  FArray1D<real, 1, 12> StemAreaIndexMon;
  FArray1D<real, 1, kMaxSoilLayers> SoilQuartzFrac;
  FArray1D<real, 1, kNumRadBands> AlbedoSoilSat;
  FArray1D<real, 1, kNumRadBands> AlbedoSoilDry;
  FArray1D<real, 1, kNumRadBands> AlbedoLakeFrz;
  FArray1D<real, 1, kNumRadBands> ScatterCoeffSnow;
  FArray1D<real, 1, kNumRadBands> ReflectanceLeaf;
  FArray1D<real, 1, kNumRadBands> ReflectanceStem;
  FArray1D<real, 1, kNumRadBands> TransmittanceLeaf;
  FArray1D<real, 1, kNumRadBands> TransmittanceStem;
  FArray1D<real, 1, 2> EmissivitySoilLake;
  FArray1D<real, 1, kNumRadBands> AlbedoLandIce;
};

struct EnergyType {
  EnergyFluxType flux;
  EnergyStateType state;
  EnergyParamType param;
};

// === WaterVarType ===
struct WaterFluxType {
  real RainfallRefHeight       = undefined_real;
  real SnowfallRefHeight       = undefined_real;
  real PrecipTotRefHeight      = undefined_real;
  real PrecipConvTotRefHeight  = undefined_real;
  real PrecipLargeSclRefHeight = undefined_real;
  real EvapCanopyNet           = undefined_real;
  real Transpiration           = undefined_real;
  real EvapCanopyLiq           = undefined_real;
  real DewCanopyLiq            = undefined_real;
  real FrostCanopyIce          = undefined_real;
  real SublimCanopyIce         = undefined_real;
  real MeltCanopyIce           = undefined_real;
  real FreezeCanopyLiq         = undefined_real;
  real SnowfallGround          = undefined_real;
  real SnowDepthIncr           = undefined_real;
  real FrostSnowSfcIce         = undefined_real;
  real SublimSnowSfcIce        = undefined_real;
  real RainfallGround          = undefined_real;
  real SnowBotOutflow          = undefined_real;
  real GlacierExcessFlow       = undefined_real;
  real IrrigationRateFlood     = real(0.0);
  real IrrigationRateMicro     = real(0.0);
  real IrrigationRateSprinkler = real(0.0);
  real IrriEvapLossSprinkler   = real(0.0);
  real SoilSfcInflow           = undefined_real;
  real RunoffSurface           = undefined_real;
  real RunoffSubsurface        = undefined_real;
  real InfilRateSfc            = undefined_real;
  real EvapSoilSfcLiq          = undefined_real;
  real DrainSoilBot            = undefined_real;
  real TileDrain               = real(0.0);
  real RechargeGw              = undefined_real;
  real DischargeGw             = undefined_real;
  real VaporizeGrd             = undefined_real;
  real CondenseVapGrd          = undefined_real;
  real DewSoilSfcLiq           = undefined_real;
  real EvapIrriSprinkler       = real(0.0);
  real InterceptCanopyRain     = undefined_real;
  real DripCanopyRain          = undefined_real;
  real ThroughfallRain         = undefined_real;
  real InterceptCanopySnow     = undefined_real;
  real DripCanopySnow          = undefined_real;
  real ThroughfallSnow         = undefined_real;
  real EvapGroundNet           = undefined_real;
  real MeltGroundSnow          = undefined_real;
  real WaterToAtmosTotal       = undefined_real;
  real EvapSoilSfcLiqAcc       = undefined_real;
  real SoilSfcInflowAcc        = undefined_real;
  real SfcWaterTotChgAcc       = undefined_real;
  real PrecipTotAcc            = undefined_real;
  real EvapCanopyNetAcc        = undefined_real;
  real TranspirationAcc        = undefined_real;
  real EvapGroundNetAcc        = undefined_real;
  real GlacierExcessFlowAcc    = undefined_real;
  real EvapSoilSfcLiqMean      = undefined_real;
  real SoilSfcInflowMean       = undefined_real;

  FArray1D<real, 1, kMaxSoilLayers> TranspWatLossSoil;
  FArray1D<real, 1, kMaxSoilLayers> TranspWatLossSoilAcc;
  FArray1D<real, 1, kMaxSoilLayers> TranspWatLossSoilMean;
  FArray1D<real, -kMaxSnowLayers + 1, 0> CompactionSnowAging;
  FArray1D<real, -kMaxSnowLayers + 1, 0> CompactionSnowBurden;
  FArray1D<real, -kMaxSnowLayers + 1, 0> CompactionSnowMelt;
  FArray1D<real, -kMaxSnowLayers + 1, 0> CompactionSnowTot;
  FArray1D<real, -kMaxSnowLayers + 1, 0> OutflowSnowLayer;
  FArray1D<real, -kMaxSnowLayers + 1, 0> SnowFreezeRate;
};

struct WaterStateType {
  int IrrigationCntSprinkler     = undefined_int;
  int IrrigationCntMicro         = undefined_int;
  int IrrigationCntFlood         = undefined_int;
  real CanopyTotalWater          = undefined_real;
  real CanopyWetFrac             = undefined_real;
  real SnowfallDensity           = undefined_real;
  real CanopyLiqWater            = undefined_real;
  real CanopyIce                 = undefined_real;
  real CanopyIceMax              = undefined_real;
  real CanopyLiqWaterMax         = undefined_real;
  real SnowDepth                 = undefined_real;
  real SnowWaterEquiv            = undefined_real;
  real SnowWaterEquivPrev        = undefined_real;
  real PondSfcThinSnwMelt        = undefined_real;
  real PondSfcThinSnwComb        = undefined_real;
  real PondSfcThinSnwTrans       = undefined_real;
  real IrrigationFracFlood       = undefined_real;
  real IrrigationAmtFlood        = undefined_real;
  real IrrigationFracMicro       = undefined_real;
  real IrrigationAmtMicro        = undefined_real;
  real IrrigationFracSprinkler   = undefined_real;
  real IrrigationAmtSprinkler    = undefined_real;
  real WaterTableDepth           = undefined_real;
  real SoilIceMax                = undefined_real;
  real SoilLiqWaterMin           = undefined_real;
  real SoilSaturateFrac          = undefined_real;
  real SoilImpervFracMax         = undefined_real;
  real SoilMoistureToWT          = undefined_real;
  real RechargeGwDeepWT          = undefined_real;
  real RechargeGwShallowWT       = undefined_real;
  real SoilSaturationExcess      = undefined_real;
  real WaterTableHydro           = undefined_real;
  real TileDrainFrac             = undefined_real;
  real WaterStorageAquifer       = undefined_real;
  real WaterStorageSoilAqf       = undefined_real;
  real WaterStorageLake          = undefined_real;
  real WaterStorageWetland       = undefined_real;
  real WaterHeadSfc              = undefined_real;
  real IrrigationFracGrid        = undefined_real;
  real PrecipAreaFrac            = undefined_real;
  real SnowCoverFrac             = undefined_real;
  real SoilTranspFacAcc          = undefined_real;
  real FrozenPrecipFrac          = undefined_real;
  real SoilWaterRootZone         = undefined_real;
  real SoilWaterStress           = undefined_real;
  real WaterStorageTotBeg        = undefined_real;
  real WaterBalanceError         = undefined_real;
  real WaterStorageTotEnd        = undefined_real;
  real SnowRadiusFresh           = undefined_real;

  FArray1D<int, -kMaxSnowLayers + 1, kMaxSoilLayers> IndexPhaseChange;
  FArray1D<real, -kMaxSnowLayers + 1, 0> SnowIce;
  FArray1D<real, -kMaxSnowLayers + 1, 0> SnowLiqWater;
  FArray1D<real, -kMaxSnowLayers + 1, 0> SnowIceFracPrev;
  FArray1D<real, -kMaxSnowLayers + 1, 0> SnowIceFrac;
  FArray1D<real, 1, kMaxSoilLayers> SoilIceFrac;
  FArray1D<real, -kMaxSnowLayers + 1, 0> SnowEffPorosity;
  FArray1D<real, 1, kMaxSoilLayers> SoilLiqWater;
  FArray1D<real, 1, kMaxSoilLayers> SoilIce;
  FArray1D<real, 1, kMaxSoilLayers> SoilMoisture;
  FArray1D<real, 1, kMaxSoilLayers> SoilImpervFrac;
  FArray1D<real, 1, kMaxSoilLayers> SoilWatConductivity;
  FArray1D<real, 1, kMaxSoilLayers> SoilWatDiffusivity;
  FArray1D<real, 1, kMaxSoilLayers> SoilEffPorosity;
  FArray1D<real, 1, kMaxSoilLayers> SoilMoistureEqui;
  FArray1D<real, 1, kMaxSoilLayers> SoilTranspFac;
  FArray1D<real, -kMaxSnowLayers + 1, 0> SnowIceVol;
  FArray1D<real, -kMaxSnowLayers + 1, 0> SnowLiqWaterVol;
  FArray1D<real, -kMaxSnowLayers + 1, kMaxSoilLayers> SoilSupercoolWater;
  FArray1D<real, 1, kMaxSoilLayers> SoilMatPotential;
};

struct WaterParamType {
  int DrainSoilLayerInd        = undefined_int;
  int TileDrainTubeDepth       = undefined_int;
  int NumSoilLayerRoot         = undefined_int;
  int IrriStopDayBfHarvest     = undefined_int;
  real CanopyLiqHoldCap        = undefined_real;
  real SnowCompactBurdenFac    = undefined_real;
  real SnowCompactAgingFac1    = undefined_real;
  real SnowCompactAgingFac2    = undefined_real;
  real SnowCompactAgingFac3    = undefined_real;
  real SnowCompactAgingMax     = undefined_real;
  real SnowViscosityCoeff      = undefined_real;
  real SnowCompactmAR24        = undefined_real;
  real SnowCompactbAR24        = undefined_real;
  real SnowCompactP1AR24       = undefined_real;
  real SnowCompactP2AR24       = undefined_real;
  real SnowCompactP3AR24       = undefined_real;
  real BurdenFacUpAR24         = undefined_real;
  real SnowCoverM1AR25         = undefined_real;
  real SnowCoverM2AR25         = undefined_real;
  real SnowCoverFac1AR25       = undefined_real;
  real SnowCoverFac2AR25       = undefined_real;
  real SnowLiqFracMax          = undefined_real;
  real SnowLiqHoldCap          = undefined_real;
  real SnowLiqReleaseFac       = undefined_real;
  real IrriFloodRateFac        = undefined_real;
  real IrriMicroRate           = undefined_real;
  real SoilInfilMaxCoeff       = undefined_real;
  real SoilImpervFracCoeff     = undefined_real;
  real InfilFacVic             = undefined_real;
  real TensionWatDistrInfl     = undefined_real;
  real TensionWatDistrShp      = undefined_real;
  real FreeWatDistrShp         = undefined_real;
  real InfilHeteroDynVic       = undefined_real;
  real InfilCapillaryDynVic    = undefined_real;
  real InfilFacDynVic          = undefined_real;
  real SoilDrainSlope          = undefined_real;
  real TileDrainCoeffSp        = undefined_real;
  real DrainFacSoilWat         = undefined_real;
  real TileDrainCoeff          = undefined_real;
  real DrainDepthToImperv      = undefined_real;
  real LateralWatCondFac       = undefined_real;
  real TileDrainDepth          = undefined_real;
  real DrainTubeDist           = undefined_real;
  real DrainTubeRadius         = undefined_real;
  real DrainWatDepToImperv     = undefined_real;
  real RunoffDecayFac          = undefined_real;
  real BaseflowCoeff           = undefined_real;
  real GridTopoIndex           = undefined_real;
  real SoilSfcSatFracMax       = undefined_real;
  real SpecYieldGw             = undefined_real;
  real MicroPoreContent        = undefined_real;
  real WaterStorageLakeMax     = undefined_real;
  real SnoWatEqvMaxGlacier     = undefined_real;
  real SoilConductivityRef     = undefined_real;
  real SoilInfilFacRef         = undefined_real;
  real GroundFrzCoeff          = undefined_real;
  real IrriTriggerLaiMin       = undefined_real;
  real SoilWatDeficitAllow     = undefined_real;
  real IrriFloodLossFrac       = undefined_real;
  real IrriSprinklerRate       = undefined_real;
  real IrriFracThreshold       = undefined_real;
  real IrriStopPrecipThr       = undefined_real;
  real SnowfallDensityMax      = undefined_real;
  real SnowMassFullCoverOld    = undefined_real;
  real SoilMatPotentialWilt    = undefined_real;
  real SnowMeltFac             = undefined_real;
  real SnowCoverFac            = undefined_real;
  real WetlandCapMax           = undefined_real;

  FArray1D<real, 1, kMaxSoilLayers> SoilMoistureSat;
  FArray1D<real, 1, kMaxSoilLayers> SoilMoistureWilt;
  FArray1D<real, 1, kMaxSoilLayers> SoilMoistureFieldCap;
  FArray1D<real, 1, kMaxSoilLayers> SoilMoistureDry;
  FArray1D<real, 1, kMaxSoilLayers> SoilWatDiffusivitySat;
  FArray1D<real, 1, kMaxSoilLayers> SoilWatConductivitySat;
  FArray1D<real, 1, kMaxSoilLayers> SoilExpCoeffB;
  FArray1D<real, 1, kMaxSoilLayers> SoilMatPotentialSat;
};

struct WaterType {
  WaterFluxType flux;
  WaterStateType state;
  WaterParamType param;
};

// === BiochemVarType ===
struct BiochemFluxType {
  real PhotosynTotal           = undefined_real;
  real PhotosynLeafSunlit      = undefined_real;
  real PhotosynLeafShade       = undefined_real;
  real PhotosynCrop            = undefined_real;
  real GrossPriProduction      = undefined_real;
  real NetEcoExchange          = undefined_real;
  real NetPriProductionTot     = undefined_real;
  real NetPriProductionLeaf    = undefined_real;
  real NetPriProductionRoot    = undefined_real;
  real NetPriProductionWood    = undefined_real;
  real NetPriProductionStem    = undefined_real;
  real NetPriProductionGrain   = undefined_real;
  real RespirationPlantTot     = undefined_real;
  real RespirationSoilOrg      = undefined_real;
  real CarbonToAtmos           = undefined_real;
  real GrowthRespLeaf          = undefined_real;
  real GrowthRespRoot          = undefined_real;
  real GrowthRespWood          = undefined_real;
  real GrowthRespStem          = undefined_real;
  real GrowthRespGrain         = undefined_real;
  real LeafMassMaxChg          = undefined_real;
  real StemMassMaxChg          = undefined_real;
  real CarbonDecayToStable     = undefined_real;
  real RespirationLeaf         = undefined_real;
  real RespirationStem         = undefined_real;
  real RespirationWood         = undefined_real;
  real RespirationLeafMaint    = undefined_real;
  real RespirationRoot         = undefined_real;
  real RespirationSoil         = undefined_real;
  real RespirationGrain        = undefined_real;
  real ConvRootToGrain         = undefined_real;
  real ConvStemToGrain         = undefined_real;
  real ConvLeafToGrain         = undefined_real;
  real TurnoverLeaf            = undefined_real;
  real TurnoverStem            = undefined_real;
  real TurnoverWood            = undefined_real;
  real TurnoverRoot            = undefined_real;
  real TurnoverGrain           = undefined_real;
  real DeathLeaf               = undefined_real;
  real DeathStem               = undefined_real;
  real CarbonAssim             = undefined_real;
  real CarbohydrAssim          = undefined_real;
};

struct BiochemStateType {
  int PlantGrowStage            = undefined_int;
  int IndexPlanting             = undefined_int;
  int IndexHarvest              = undefined_int;
  real IndexGrowSeason          = undefined_real;
  real NitrogenConcFoliage      = undefined_real;
  real LeafMass                 = undefined_real;
  real RootMass                 = undefined_real;
  real StemMass                 = undefined_real;
  real WoodMass                 = undefined_real;
  real GrainMass                = undefined_real;
  real CarbonMassDeepSoil       = undefined_real;
  real CarbonMassShallowSoil    = undefined_real;
  real CarbonMassSoilTot        = undefined_real;
  real CarbonMassLiveTot        = undefined_real;
  real LeafAreaPerMass          = undefined_real;
  real StemAreaPerMass          = undefined_real;
  real LeafMassMin              = undefined_real;
  real StemMassMin              = undefined_real;
  real CarbonFracToLeaf         = undefined_real;
  real CarbonFracToRoot         = undefined_real;
  real CarbonFracToWood         = undefined_real;
  real CarbonFracToStem         = undefined_real;
  real WoodCarbonFrac           = undefined_real;
  real CarbonFracToWoodRoot     = undefined_real;
  real MicroRespFactorSoilWater = undefined_real;
  real MicroRespFactorSoilTemp  = undefined_real;
  real RespFacNitrogenFoliage   = undefined_real;
  real RespFacTemperature       = undefined_real;
  real RespReductionFac         = undefined_real;
  real GrowDegreeDay            = undefined_real;
};

struct BiochemParamType {
  int DatePlanting             = undefined_int;
  int DateHarvest              = undefined_int;
  real QuantumEfficiency25C    = undefined_real;
  real CarboxylRateMax25C      = undefined_real;
  real CarboxylRateMaxQ10      = undefined_real;
  real PhotosynPathC3          = undefined_real;
  real SlopeConductToPhotosyn  = undefined_real;
  real TemperatureMinPhotosyn  = undefined_real;
  real LeafAreaPerMass1side    = undefined_real;
  real NitrogenConcFoliageMax  = undefined_real;
  real WoodToRootRatio         = undefined_real;
  real WoodPoolIndex           = undefined_real;
  real TurnoverCoeffLeafVeg    = undefined_real;
  real LeafDeathWaterCoeffVeg  = undefined_real;
  real LeafDeathTempCoeffVeg   = undefined_real;
  real MicroRespCoeff          = undefined_real;
  real RespMaintQ10            = undefined_real;
  real RespMaintLeaf25C        = undefined_real;
  real RespMaintStem25C        = undefined_real;
  real RespMaintRoot25C        = undefined_real;
  real RespMaintGrain25C       = undefined_real;
  real GrowthRespFrac          = undefined_real;
  real TemperaureLeafFreeze    = undefined_real;
  real LeafAreaPerBiomass      = undefined_real;
  real TempBaseGrowDegDay      = undefined_real;
  real TempMaxGrowDegDay       = undefined_real;
  real GrowDegDayEmerg         = undefined_real;
  real GrowDegDayInitVeg       = undefined_real;
  real GrowDegDayPostVeg       = undefined_real;
  real GrowDegDayInitReprod    = undefined_real;
  real GrowDegDayMature        = undefined_real;
  real PhotosynRadFrac         = undefined_real;
  real TempMinCarbonAssim      = undefined_real;
  real TempMaxCarbonAssim      = undefined_real;
  real TempMaxCarbonAssimMax   = undefined_real;
  real CarbonAssimRefMax       = undefined_real;
  real LightExtCoeff           = undefined_real;
  real LightUseEfficiency      = undefined_real;
  real CarbonAssimReducFac     = undefined_real;
  real StemAreaIndexMin        = undefined_real;
  real WoodAllocFac            = undefined_real;
  real WaterStressCoeff        = undefined_real;
  real LeafAreaIndexMin        = undefined_real;
  real TurnoverCoeffRootVeg    = undefined_real;
  real WoodRespCoeff           = undefined_real;

  FArray1D<real, 1, kNumCropStages> LeafDeathTempCoeffCrop;
  FArray1D<real, 1, kNumCropStages> LeafDeathWaterCoeffCrop;
  FArray1D<real, 1, kNumCropStages> CarbohydrLeafToGrain;
  FArray1D<real, 1, kNumCropStages> CarbohydrStemToGrain;
  FArray1D<real, 1, kNumCropStages> CarbohydrRootToGrain;
  FArray1D<real, 1, kNumCropStages> CarbohydrFracToLeaf;
  FArray1D<real, 1, kNumCropStages> CarbohydrFracToStem;
  FArray1D<real, 1, kNumCropStages> CarbohydrFracToRoot;
  FArray1D<real, 1, kNumCropStages> CarbohydrFracToGrain;
  FArray1D<real, 1, kNumCropStages> TurnoverCoeffLeafCrop;
  FArray1D<real, 1, kNumCropStages> TurnoverCoeffStemCrop;
  FArray1D<real, 1, kNumCropStages> TurnoverCoeffRootCrop;
};

struct BiochemType {
  BiochemFluxType flux;
  BiochemStateType state;
  BiochemParamType param;
};

struct NoahmpType {
  ConfigType config;
  ForcingType forcing;
  EnergyType energy;
  WaterType water;
  BiochemType biochem;

  CPPWRF_HD void InitDefault() {
    config = ConfigType{};
    config.domain.SoilType.fill(undefined_int);
    config.domain.DepthSoilLayer.fill(undefined_real);
    config.domain.ThicknessSoilLayer.fill(undefined_real);
    config.domain.ThicknessSnowSoilLayer.fill(undefined_real);
    config.domain.DepthSnowSoilLayer.fill(undefined_real);

    forcing = ForcingType{};

    energy = EnergyType{};
    energy.state.TemperatureSoilSnow.fill(undefined_real);
    energy.state.ThermConductSoilSnow.fill(undefined_real);
    energy.state.HeatCapacSoilSnow.fill(undefined_real);
    energy.state.PhaseChgFacSoilSnow.fill(undefined_real);
    energy.state.HeatCapacVolSnow.fill(undefined_real);
    energy.state.ThermConductSnow.fill(undefined_real);
    energy.state.HeatCapacVolSoil.fill(undefined_real);
    energy.state.ThermConductSoil.fill(undefined_real);
    energy.state.HeatCapacGlaIce.fill(undefined_real);
    energy.state.ThermConductGlaIce.fill(undefined_real);
    energy.state.AlbedoSnowDir.fill(undefined_real);
    energy.state.AlbedoSnowDif.fill(undefined_real);
    energy.state.AlbedoSoilDir.fill(real(0.0));
    energy.state.AlbedoSoilDif.fill(real(0.0));
    energy.state.AlbedoGrdDir.fill(undefined_real);
    energy.state.AlbedoGrdDif.fill(undefined_real);
    energy.state.ReflectanceVeg.fill(undefined_real);
    energy.state.TransmittanceVeg.fill(undefined_real);
    energy.state.AlbedoSfcDir.fill(undefined_real);
    energy.state.AlbedoSfcDif.fill(undefined_real);

    energy.flux.RadSwAbsVegDir.fill(undefined_real);
    energy.flux.RadSwAbsVegDif.fill(undefined_real);
    energy.flux.RadSwDirTranGrdDir.fill(undefined_real);
    energy.flux.RadSwDirTranGrdDif.fill(undefined_real);
    energy.flux.RadSwDifTranGrdDir.fill(undefined_real);
    energy.flux.RadSwDifTranGrdDif.fill(undefined_real);
    energy.flux.RadSwReflVegDir.fill(undefined_real);
    energy.flux.RadSwReflVegDif.fill(undefined_real);
    energy.flux.RadSwReflGrdDir.fill(undefined_real);
    energy.flux.RadSwReflGrdDif.fill(undefined_real);
    energy.flux.RadSwDownDir.fill(undefined_real);
    energy.flux.RadSwDownDif.fill(undefined_real);
    energy.flux.RadSwPenetrateGrd.fill(undefined_real);

    energy.param.LeafAreaIndexMon.fill(undefined_real);
    energy.param.StemAreaIndexMon.fill(undefined_real);
    energy.param.SoilQuartzFrac.fill(undefined_real);
    energy.param.AlbedoSoilSat.fill(undefined_real);
    energy.param.AlbedoSoilDry.fill(undefined_real);
    energy.param.AlbedoLakeFrz.fill(undefined_real);
    energy.param.ScatterCoeffSnow.fill(undefined_real);
    energy.param.ReflectanceLeaf.fill(undefined_real);
    energy.param.ReflectanceStem.fill(undefined_real);
    energy.param.TransmittanceLeaf.fill(undefined_real);
    energy.param.TransmittanceStem.fill(undefined_real);
    energy.param.EmissivitySoilLake.fill(undefined_real);
    energy.param.AlbedoLandIce.fill(undefined_real);

    water = WaterType{};
    water.state.IndexPhaseChange.fill(undefined_int);
    water.state.SoilSupercoolWater.fill(undefined_real);
    water.state.SnowIce.fill(undefined_real);
    water.state.SnowLiqWater.fill(undefined_real);
    water.state.SnowIceVol.fill(undefined_real);
    water.state.SnowLiqWaterVol.fill(undefined_real);
    water.state.SnowIceFracPrev.fill(undefined_real);
    water.state.SnowIceFrac.fill(undefined_real);
    water.state.SoilIceFrac.fill(undefined_real);
    water.state.SnowEffPorosity.fill(undefined_real);
    water.state.SoilLiqWater.fill(undefined_real);
    water.state.SoilIce.fill(undefined_real);
    water.state.SoilMoisture.fill(undefined_real);
    water.state.SoilImpervFrac.fill(undefined_real);
    water.state.SoilWatConductivity.fill(undefined_real);
    water.state.SoilWatDiffusivity.fill(undefined_real);
    water.state.SoilEffPorosity.fill(undefined_real);
    water.state.SoilMoistureEqui.fill(undefined_real);
    water.state.SoilTranspFac.fill(undefined_real);
    water.state.SoilMatPotential.fill(undefined_real);

    water.flux.CompactionSnowAging.fill(undefined_real);
    water.flux.CompactionSnowBurden.fill(undefined_real);
    water.flux.CompactionSnowMelt.fill(undefined_real);
    water.flux.CompactionSnowTot.fill(undefined_real);
    water.flux.TranspWatLossSoil.fill(undefined_real);
    water.flux.TranspWatLossSoilAcc.fill(undefined_real);
    water.flux.TranspWatLossSoilMean.fill(undefined_real);
    water.flux.OutflowSnowLayer.fill(undefined_real);
    water.flux.SnowFreezeRate.fill(undefined_real);

    water.param.SoilMoistureSat.fill(undefined_real);
    water.param.SoilMoistureWilt.fill(undefined_real);
    water.param.SoilMoistureFieldCap.fill(undefined_real);
    water.param.SoilMoistureDry.fill(undefined_real);
    water.param.SoilWatDiffusivitySat.fill(undefined_real);
    water.param.SoilWatConductivitySat.fill(undefined_real);
    water.param.SoilExpCoeffB.fill(undefined_real);
    water.param.SoilMatPotentialSat.fill(undefined_real);

    biochem = BiochemType{};
    biochem.param.LeafDeathTempCoeffCrop.fill(undefined_real);
    biochem.param.LeafDeathWaterCoeffCrop.fill(undefined_real);
    biochem.param.CarbohydrLeafToGrain.fill(undefined_real);
    biochem.param.CarbohydrStemToGrain.fill(undefined_real);
    biochem.param.CarbohydrRootToGrain.fill(undefined_real);
    biochem.param.CarbohydrFracToLeaf.fill(undefined_real);
    biochem.param.CarbohydrFracToStem.fill(undefined_real);
    biochem.param.CarbohydrFracToRoot.fill(undefined_real);
    biochem.param.CarbohydrFracToGrain.fill(undefined_real);
    biochem.param.TurnoverCoeffLeafCrop.fill(undefined_real);
    biochem.param.TurnoverCoeffStemCrop.fill(undefined_real);
    biochem.param.TurnoverCoeffRootCrop.fill(undefined_real);
  }
};

using NoahmpVarType = NoahmpType;

}  // namespace cppwrf::phys::noahmp
