#pragma once

#include <string>
#include "core/precision.hpp"

namespace cppwrf::phys::noahmp {

static constexpr int kMaxVegCats = 65;     // MVT = 61 in NoahmpReadTableMod (USGS/MODIS + LCZ 51..61)
static constexpr int kMaxSoilCats = 30;    // MAX_SOILTYP = 30 in NoahmpReadTableMod
static constexpr int kMaxSlopeCats = 9;    // NUM_SLOPE = 9
static constexpr int kMaxSoilColors = 8;   // MSC = 8
static constexpr int kNumRadBands = 2;     // MBAND = 2 (1=vis, 2=nir)
static constexpr int kNumMonths = 12;      // 12 months
static constexpr int kNumCropCats = 5;     // NCROP = 5
static constexpr int kNumCropStages = 8;   // NSTAGE = 8
static constexpr int kMaxSnowLayers = 3;   // NSNOW = 3 (-2..0)
static constexpr int kMaxSoilLayers = 4;   // NSOIL = 4 (1..4)

// Namelist options for Noah-MP (WRF Registry.EM_COMMON defaults)
struct NoahmpOptions {
  int dveg = 4;          // dynamic_veg_option (4 = table LAI, max vegfrac)
  int opt_crs = 1;       // canopy_stomatal_resistance_option (1 = Ball-Berry)
  int opt_btr = 1;       // btr_option (1 = Noah)
  int opt_runsrf = 3;    // surface_runoff_option (3 = Schaake96)
  int opt_runsub = 3;    // subsurface_runoff_option (3 = Schaake96)
  int opt_sfc = 1;       // surface_drag_option (1 = M-O)
  int opt_frz = 1;       // supercooled_water_option (1 = NY06)
  int opt_inf = 1;       // frozen_soil_option (1 = NY06)
  int opt_rad = 3;       // radiative_transfer_option (3 = gap=1-Fveg)
  int opt_alb = 2;       // snow_albedo_option (2 = CLASS)
  int opt_snf = 1;       // pcp_partition_option (1 = Jordan91)
  int opt_tbot = 2;      // tbot_option (2 = Noah)
  int opt_stc = 1;       // temp_time_scheme_option (1 = semi-implicit)
  int opt_gla = 1;       // glacier_option (1 = phase change)
  int opt_rsf = 1;       // surface_resistance_option (1 = Sakaguchi and Zeng 2009)
  int opt_soil = 1;      // soil_data_option (1 = input dominant soil texture)
  int opt_pedo = 1;      // pedotransfer_option (1 = Saxton and Rawls 2006)
  int opt_crop = 0;      // crop_option (0 = none)
  int opt_irr = 0;       // irrigation_option (0 = none)
  int opt_irrm = 0;      // irrigation_method (0 = method_isp)
  int opt_infdv = 0;     // dvic_infiltration_option (0 = default)
  int opt_tdrn = 0;      // tile_drainage_option (0 = none)
  int opt_tksno = 1;     // snow_thermal_cond (1 = Yen 1965)
  int opt_compact = 2;   // snow_compaction (2 = Abolafia-Rosenzweig 2024)
  int opt_scf = 1;       // snowcover_fraction (1 = Niu and Yang 2007)
  int opt_wetland = 0;   // wetland_option (0 = off)
  int sf_urban_physics = 0;
  int iz0tlnd = 0;
  double soiltstep = 0.0;
};

// Throws std::runtime_error if any option in `opt` is not supported by this port.
void validate_noahmp_options(const NoahmpOptions& opt);

// Parameter tables read from NoahmpTable.TBL (1-based indexing in index 1..N; index 0 unused,
// or 0-based where noted; here we store 1-based arrays [kMaxVegCats + 1] so table(ivgtyp)
// matches Fortran NoahmpIO%*_TABLE(IVGTYP) directly).
struct NoahmpTables {
  // Vegetation category indicators
  int nveg = 0;
  int isurban = -9999;
  int iswater = -9999;
  int isbarren = -9999;
  int isice = -9999;
  int iscrop = -9999;
  int eblforest = -9999;
  int natural = -9999;
  int urbtype_beg = -9999;
  int lcz_1 = -9999;
  int lcz_2 = -9999;
  int lcz_3 = -9999;
  int lcz_4 = -9999;
  int lcz_5 = -9999;
  int lcz_6 = -9999;
  int lcz_7 = -9999;
  int lcz_8 = -9999;
  int lcz_9 = -9999;
  int lcz_10 = -9999;
  int lcz_11 = -9999;

  // Vegetation parameters (1-indexed: [1..kMaxVegCats])
  real ch2op[kMaxVegCats + 1];
  real dleaf[kMaxVegCats + 1];
  real z0mvt[kMaxVegCats + 1];
  real hvt[kMaxVegCats + 1];
  real hvb[kMaxVegCats + 1];
  real den[kMaxVegCats + 1];
  real rc[kMaxVegCats + 1];
  real mfsno[kMaxVegCats + 1];
  real scffac[kMaxVegCats + 1];
  real cbiom[kMaxVegCats + 1];
  real xl[kMaxVegCats + 1];
  real cwpvt[kMaxVegCats + 1];
  real c3psn[kMaxVegCats + 1];
  real kc25[kMaxVegCats + 1];
  real akc[kMaxVegCats + 1];
  real ko25[kMaxVegCats + 1];
  real ako[kMaxVegCats + 1];
  real avcmx[kMaxVegCats + 1];
  real aqe[kMaxVegCats + 1];
  real ltovrc[kMaxVegCats + 1];
  real dilefc[kMaxVegCats + 1];
  real dilefw[kMaxVegCats + 1];
  real rmf25[kMaxVegCats + 1];
  real sla[kMaxVegCats + 1];
  real fragr[kMaxVegCats + 1];
  real tmin[kMaxVegCats + 1];
  real vcmx25[kMaxVegCats + 1];
  real tdlef[kMaxVegCats + 1];
  real bp[kMaxVegCats + 1];
  real mp[kMaxVegCats + 1];
  real qe25[kMaxVegCats + 1];
  real rms25[kMaxVegCats + 1];
  real rmr25[kMaxVegCats + 1];
  real arm[kMaxVegCats + 1];
  real folnmx[kMaxVegCats + 1];
  real wdpool[kMaxVegCats + 1];
  real wrrat[kMaxVegCats + 1];
  real mrp[kMaxVegCats + 1];
  real nroot[kMaxVegCats + 1];
  real rgl[kMaxVegCats + 1];
  real rs[kMaxVegCats + 1];
  real hs[kMaxVegCats + 1];
  real topt[kMaxVegCats + 1];
  real rsmax[kMaxVegCats + 1];
  real rtovrc[kMaxVegCats + 1];
  real rswoodc[kMaxVegCats + 1];
  real bf[kMaxVegCats + 1];
  real wstrc[kMaxVegCats + 1];
  real laimin[kMaxVegCats + 1];
  real xsamin[kMaxVegCats + 1];

  // Monthly SAI and LAI (1-indexed: [1..kMaxVegCats][1..12])
  real saim[kMaxVegCats + 1][kNumMonths + 1];
  real laim[kMaxVegCats + 1][kNumMonths + 1];

  // Band-dependent leaf/stem reflectance & transmittance (1-indexed: [1..kMaxVegCats][1..2])
  real rhol[kMaxVegCats + 1][kNumRadBands + 1];
  real rhos[kMaxVegCats + 1][kNumRadBands + 1];
  real taul[kMaxVegCats + 1][kNumRadBands + 1];
  real taus[kMaxVegCats + 1][kNumRadBands + 1];

  // Soil parameters (1-indexed: [1..kMaxSoilCats])
  int slcats = -9999;
  real bexp[kMaxSoilCats + 1];
  real smcdry[kMaxSoilCats + 1];
  real smcmax[kMaxSoilCats + 1];
  real smcref[kMaxSoilCats + 1];
  real psisat[kMaxSoilCats + 1];
  real dksat[kMaxSoilCats + 1];
  real dwsat[kMaxSoilCats + 1];
  real smcwlt[kMaxSoilCats + 1];
  real quartz[kMaxSoilCats + 1];
  real bvic[kMaxSoilCats + 1];
  real axaj[kMaxSoilCats + 1];
  real bxaj[kMaxSoilCats + 1];
  real xxaj[kMaxSoilCats + 1];
  real bdvic[kMaxSoilCats + 1];
  real gdvic[kMaxSoilCats + 1];
  real bbvic[kMaxSoilCats + 1];

  // General parameters
  real slope[kMaxSlopeCats + 1];
  real csoil = real(-9999.0);
  real refdk = real(-9999.0);
  real refkdt = real(-9999.0);
  real frzk = real(-9999.0);
  real zbot = real(-9999.0);
  real czil = real(-9999.0);

  // Radiation parameters (1-indexed)
  real albsat[kMaxSoilColors + 1][kNumRadBands + 1];
  real albdry[kMaxSoilColors + 1][kNumRadBands + 1];
  real albice[kNumRadBands + 1];
  real alblak[kNumRadBands + 1];
  real omegas[kNumRadBands + 1];
  real betads = real(-9999.0);
  real betais = real(-9999.0);
  real eg[3];  // [1..2]
  real eice = real(-9999.0);

  // Global parameters
  real co2 = real(-9999.0);
  real o2 = real(-9999.0);
  real timean = real(-9999.0);
  real fsatmx = real(-9999.0);
  real z0sno = real(-9999.0);
  real ssi = real(-9999.0);
  real snow_ret_fac = real(-9999.0);
  real snow_emis = real(-9999.0);
  real swemx = real(-9999.0);
  real tau0 = real(-9999.0);
  real grain_growth = real(-9999.0);
  real extra_growth = real(-9999.0);
  real dirt_soot = real(-9999.0);
  real bats_cosz = real(-9999.0);
  real bats_vis_new = real(-9999.0);
  real bats_nir_new = real(-9999.0);
  real bats_vis_age = real(-9999.0);
  real bats_nir_age = real(-9999.0);
  real bats_vis_dir = real(-9999.0);
  real bats_nir_dir = real(-9999.0);
  real rsurf_snow = real(-9999.0);
  real rsurf_exp = real(-9999.0);
  real c2_snowcompact = real(-9999.0);
  real c3_snowcompact = real(-9999.0);
  real c4_snowcompact = real(-9999.0);
  real c5_snowcompact = real(-9999.0);
  real dm_snowcompact = real(-9999.0);
  real eta0_snowcompact = real(-9999.0);
  real snowcompactm_ar24 = real(-9999.0);
  real snowcompactb_ar24 = real(-9999.0);
  real snowcompact_p1_ar24 = real(-9999.0);
  real snowcompact_p2_ar24 = real(-9999.0);
  real snowcompact_p3_ar24 = real(-9999.0);
  real snowcompact_up_ar24 = real(-9999.0);
  real scfm1_ar25 = real(-9999.0);
  real scfm2_ar25 = real(-9999.0);
  real scfac1_ar25 = real(-9999.0);
  real scfac2_ar25 = real(-9999.0);
  real snliqmaxfrac = real(-9999.0);
  real swemaxgla = real(-9999.0);
  real wslmax = real(-9999.0);
  real rous = real(-9999.0);
  real cmic = real(-9999.0);
  real snowden_max = real(-9999.0);
  real class_alb_ref = real(-9999.0);
  real class_sno_age = real(-9999.0);
  real class_alb_new = real(-9999.0);
  real psiwlt = real(-9999.0);
  real z0soil = real(-9999.0);
  real z0lake = real(-9999.0);
  real wcap = real(-9999.0);

  // Irrigation parameters
  int irr_har = -9999;
  real irr_frac = real(-9999.0);
  real irr_lai = real(-9999.0);
  real irr_mad = real(-9999.0);
  real filoss = real(-9999.0);
  real sprir_rate = real(-9999.0);
  real micir_rate = real(-9999.0);
  real firtfac = real(-9999.0);
  real ir_rain = real(-9999.0);

  // Crop parameters (1-indexed: [1..kNumCropCats], stages [1..kNumCropStages])
  int default_crop = -9999;
  int pltday[kNumCropCats + 1];
  int hsday[kNumCropCats + 1];
  real plantpop[kNumCropCats + 1];
  real irri[kNumCropCats + 1];
  real gddtbase[kNumCropCats + 1];
  real gddtcut[kNumCropCats + 1];
  real gdds1[kNumCropCats + 1];
  real gdds2[kNumCropCats + 1];
  real gdds3[kNumCropCats + 1];
  real gdds4[kNumCropCats + 1];
  real gdds5[kNumCropCats + 1];
  real c3psni[kNumCropCats + 1];
  real kc25i[kNumCropCats + 1];
  real akci[kNumCropCats + 1];
  real ko25i[kNumCropCats + 1];
  real akoi[kNumCropCats + 1];
  real avcmxi[kNumCropCats + 1];
  real vcmx25i[kNumCropCats + 1];
  real bpi[kNumCropCats + 1];
  real mpi[kNumCropCats + 1];
  real folnmxi[kNumCropCats + 1];
  real qe25i[kNumCropCats + 1];
  real aref[kNumCropCats + 1];
  real psnrf[kNumCropCats + 1];
  real i2par[kNumCropCats + 1];
  real tassim0[kNumCropCats + 1];
  real tassim1[kNumCropCats + 1];
  real tassim2[kNumCropCats + 1];
  real k_crop[kNumCropCats + 1];
  real epsi[kNumCropCats + 1];
  real q10mr[kNumCropCats + 1];
  real lefreez[kNumCropCats + 1];
  real fra_gr[kNumCropCats + 1];
  real lfmr25[kNumCropCats + 1];
  real stmr25[kNumCropCats + 1];
  real rtmr25[kNumCropCats + 1];
  real grainmr25[kNumCropCats + 1];
  real bio2lai[kNumCropCats + 1];
  real dile_fc[kNumCropCats + 1][kNumCropStages + 1];
  real dile_fw[kNumCropCats + 1][kNumCropStages + 1];
  real lf_ovrc[kNumCropCats + 1][kNumCropStages + 1];
  real st_ovrc[kNumCropCats + 1][kNumCropStages + 1];
  real rt_ovrc[kNumCropCats + 1][kNumCropStages + 1];
  real lfpt[kNumCropCats + 1][kNumCropStages + 1];
  real stpt[kNumCropCats + 1][kNumCropStages + 1];
  real rtpt[kNumCropCats + 1][kNumCropStages + 1];
  real grainpt[kNumCropCats + 1][kNumCropStages + 1];
  real lfct[kNumCropCats + 1][kNumCropStages + 1];
  real stct[kNumCropCats + 1][kNumCropStages + 1];
  real rtct[kNumCropCats + 1][kNumCropStages + 1];

  // Tile drainage parameters (1-indexed: [1..kMaxSoilCats])
  int drain_layer_opt = -9999;
  int td_depth[kMaxSoilCats + 1];
  real tdsmc_fac[kMaxSoilCats + 1];
  real td_dc[kMaxSoilCats + 1];
  real td_dcoef[kMaxSoilCats + 1];
  real td_d[kMaxSoilCats + 1];
  real td_adepth[kMaxSoilCats + 1];
  real td_radi[kMaxSoilCats + 1];
  real td_spac[kMaxSoilCats + 1];
  real td_ddrain[kMaxSoilCats + 1];
  real klat_fac[kMaxSoilCats + 1];

  // Optional SR2006 pedotransfer parameters
  real sr2006_theta_1500t_a = real(-9999.0);
  real sr2006_theta_1500t_b = real(-9999.0);
  real sr2006_theta_1500t_c = real(-9999.0);
  real sr2006_theta_1500t_d = real(-9999.0);
  real sr2006_theta_1500t_e = real(-9999.0);
  real sr2006_theta_1500t_f = real(-9999.0);
  real sr2006_theta_1500t_g = real(-9999.0);
  real sr2006_theta_1500_a = real(-9999.0);
  real sr2006_theta_1500_b = real(-9999.0);
  real sr2006_theta_33t_a = real(-9999.0);
  real sr2006_theta_33t_b = real(-9999.0);
  real sr2006_theta_33t_c = real(-9999.0);
  real sr2006_theta_33t_d = real(-9999.0);
  real sr2006_theta_33t_e = real(-9999.0);
  real sr2006_theta_33t_f = real(-9999.0);
  real sr2006_theta_33t_g = real(-9999.0);
  real sr2006_theta_33_a = real(-9999.0);
  real sr2006_theta_33_b = real(-9999.0);
  real sr2006_theta_33_c = real(-9999.0);
  real sr2006_theta_s33t_a = real(-9999.0);
  real sr2006_theta_s33t_b = real(-9999.0);
  real sr2006_theta_s33t_c = real(-9999.0);
  real sr2006_theta_s33t_d = real(-9999.0);
  real sr2006_theta_s33t_e = real(-9999.0);
  real sr2006_theta_s33t_f = real(-9999.0);
  real sr2006_theta_s33t_g = real(-9999.0);
  real sr2006_theta_s33_a = real(-9999.0);
  real sr2006_theta_s33_b = real(-9999.0);
  real sr2006_psi_et_a = real(-9999.0);
  real sr2006_psi_et_b = real(-9999.0);
  real sr2006_psi_et_c = real(-9999.0);
  real sr2006_psi_et_d = real(-9999.0);
  real sr2006_psi_et_e = real(-9999.0);
  real sr2006_psi_et_f = real(-9999.0);
  real sr2006_psi_et_g = real(-9999.0);
  real sr2006_psi_e_a = real(-9999.0);
  real sr2006_psi_e_b = real(-9999.0);
  real sr2006_psi_e_c = real(-9999.0);
  real sr2006_smcmax_a = real(-9999.0);
  real sr2006_smcmax_b = real(-9999.0);

  // SNICAR parameters
  real depbchydropho = real(-9999.0);
  real depbchydrophi = real(-9999.0);
  real depochydropho = real(-9999.0);
  real depochydrophi = real(-9999.0);
  real depdust1 = real(-9999.0);
  real depdust2 = real(-9999.0);
  real depdust3 = real(-9999.0);
  real depdust4 = real(-9999.0);
  real depdust5 = real(-9999.0);
  real snowradiusmin = real(-9999.0);
  real freshsnowradiusmax = real(-9999.0);
  real snowradiusrefrz = real(-9999.0);
  real scaveffmeltscale = real(-9999.0);
  real scaveffmeltbcphi = real(-9999.0);
  real scaveffmeltbcpho = real(-9999.0);
  real scaveffmeltocphi = real(-9999.0);
  real scaveffmeltocpho = real(-9999.0);
  real scaveffmeltdust1 = real(-9999.0);
  real scaveffmeltdust2 = real(-9999.0);
  real scaveffmeltdust3 = real(-9999.0);
  real scaveffmeltdust4 = real(-9999.0);
  real scaveffmeltdust5 = real(-9999.0);
  real snowradiusmax = real(-9999.0);
  real snowwetagec1brun89 = real(-9999.0);
  real snowwetagec2brun89 = real(-9999.0);
  real snowagescalefac = real(-9999.0);
};

// Reads NoahmpTable.TBL from `table_dir` (or direct file path if it ends with .TBL).
bool read_noahmp_tables(const std::string& table_dir,
                        const std::string& veg_dataset,
                        NoahmpTables& tables);

}  // namespace cppwrf::phys::noahmp
