#include "wrf_io/global_atts.hpp"

#include <cctype>
#include <cstdio>
#include <cstdlib>

namespace cppwrf::wrf_io {

void julgmt(const std::string& date, int& julyr, int& julday, double& gmt) {
  const int y = std::atoi(date.substr(0, 4).c_str()), m = std::atoi(date.substr(5, 2).c_str()),
            dd = std::atoi(date.substr(8, 2).c_str()), h = std::atoi(date.substr(11, 2).c_str()),
            mi = std::atoi(date.substr(14, 2).c_str());
  static const int mdays[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  const bool leap = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
  julday = dd;
  for (int k = 1; k < m; ++k) julday += mdays[k - 1] + (k == 2 && leap ? 1 : 0);
  julyr = y;
  gmt = double(h) + double(mi) / 60.;
}

std::string add_seconds(const std::string& date, long seconds) {
  int y = std::atoi(date.substr(0, 4).c_str()), m = std::atoi(date.substr(5, 2).c_str()),
      d = std::atoi(date.substr(8, 2).c_str());
  const int h = std::atoi(date.substr(11, 2).c_str()), mi = std::atoi(date.substr(14, 2).c_str()),
            s = std::atoi(date.substr(17, 2).c_str());
  long t = ((long(h) * 60 + mi) * 60 + s) + seconds;
  long days = t / 86400;
  t %= 86400;
  if (t < 0) { t += 86400; --days; }
  static const int mdays[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  auto ml = [&](int yy, int mm) {
    const bool leap = (yy % 4 == 0 && yy % 100 != 0) || yy % 400 == 0;
    return mdays[mm - 1] + (mm == 2 && leap ? 1 : 0);
  };
  d += int(days);
  while (d > ml(y, m)) { d -= ml(y, m); if (++m > 12) { m = 1; ++y; } }
  while (d < 1) { if (--m < 1) { m = 12; --y; } d += ml(y, m); }
  char b[32];
  std::snprintf(b, sizeof b, "%04d-%02d-%02d_%02ld:%02ld:%02ld", y, m, d, t / 3600, (t / 60) % 60, t % 60);
  return b;
}

void write_global_atts(NcOut& nc, const WrfConfig& c, const DomainAtts& d, Stream s) {
  const bool input = s == Stream::input, history = s == Stream::history, restart = s == Stream::restart;
  auto r = [&](const std::string& n, double v) {  // wrf_put_dom_ti_real: WRF's REAL kind
    if (d.real_bytes == 8) nc.att_d(n, v); else nc.att_f(n, float(v));
  };
  auto i = [&](const std::string& n, int v) { nc.att(n, v); };
  auto ci = [&](const std::string& n, const char* opt) { nc.att(n, c.integer(opt)); };
  auto cr = [&](const std::string& n, const char* opt) { r(n, c.num(opt)); };

  nc.att("TITLE", d.title);
  nc.att("START_DATE", d.start_date);
  if (input) nc.att("SIMULATION_START_DATE", d.start_date);
  else if (restart || history) nc.att("SIMULATION_START_DATE", d.simulation_start_date);
  if (restart) {
    i("FLAG_RESTART", 1);
    i("MAX_WRF_ALARMS", 55);
    for (const auto& [n, v] : d.restart_alarms) i(n, v);
  }
  i("WEST-EAST_GRID_DIMENSION", d.ide);
  i("SOUTH-NORTH_GRID_DIMENSION", d.jde);
  i("BOTTOM-TOP_GRID_DIMENSION", d.kde);
  cr("DX", "dx"); cr("DY", "dy");
  ci("AERCU_OPT", "aercu_opt"); cr("AERCU_FCT", "aercu_fct");
  ci("IDEAL_CASE", "ideal_case"); ci("DIFF_6TH_SLOPEOPT", "diff_6th_slopeopt");
  ci("AUTO_LEVELS_OPT", "auto_levels_opt"); cr("DIFF_6TH_THRESH", "diff_6th_thresh");
  cr("DZBOT", "dzbot"); cr("DZSTRETCH_S", "dzstretch_s"); cr("DZSTRETCH_U", "dzstretch_u");
  if (history) {
    ci("SKEBS_ON", "skebs_on");
    if (c.integer("perturb_bdy") != 0) ci("PERTURB_BDY", "perturb_bdy");
    ci("USE_Q_DIABATIC", "use_q_diabatic");
  }
  nc.att("GRIDTYPE", std::string("C"));
  ci("DIFF_OPT", "diff_opt"); ci("KM_OPT", "km_opt"); ci("DAMP_OPT", "damp_opt");
  cr("DAMPCOEF", "dampcoef"); cr("KHDIF", "khdif"); cr("KVDIF", "kvdif");
  ci("MP_PHYSICS", "mp_physics"); ci("RA_LW_PHYSICS", "ra_lw_physics"); ci("RA_SW_PHYSICS", "ra_sw_physics");
  ci("SF_SFCLAY_PHYSICS", "sf_sfclay_physics"); ci("SF_SURFACE_PHYSICS", "sf_surface_physics");
  ci("BL_PBL_PHYSICS", "bl_pbl_physics"); ci("CU_PHYSICS", "cu_physics"); ci("SF_LAKE_PHYSICS", "sf_lake_physics");
  ci("SURFACE_INPUT_SOURCE", "surface_input_source"); ci("SST_UPDATE", "sst_update"); ci("GHG_INPUT", "ghg_input");
  ci("GRID_FDDA", "grid_fdda"); ci("GFDDA_INTERVAL_M", "auxinput10_interval_m"); ci("GFDDA_END_H", "auxinput10_end_h");
  ci("GRID_SFDDA", "grid_sfdda"); ci("SGFDDA_INTERVAL_M", "auxinput9_interval_m"); ci("SGFDDA_END_H", "auxinput9_end_h");
  ci("HYPSOMETRIC_OPT", "hypsometric_opt"); ci("USE_THETA_M", "use_theta_m");
  if (input) { ci("USE_MAXW_LEVEL", "use_maxw_level"); ci("USE_TROP_LEVEL", "use_trop_level"); }
  ci("GWD_OPT", "gwd_opt"); ci("SF_URBAN_PHYSICS", "sf_urban_physics");
  ci("SF_SURFACE_MOSAIC", "sf_surface_mosaic"); ci("SF_OCEAN_PHYSICS", "sf_ocean_physics");
  i("SLUCM_DISTRIBUTED_DRAG", c.logical("slucm_distributed_drag") ? 1 : 0);
  ci("DISTRIBUTED_AHE_OPT", "distributed_ahe_opt");
  if (history) {
    ci("SHCU_PHYSICS", "shcu_physics"); ci("MFSHCONV", "mfshconv"); ci("FEEDBACK", "feedback");
    ci("SMOOTH_OPTION", "smooth_option"); cr("SWRAD_SCAT", "swrad_scat"); ci("W_DAMPING", "w_damping");
    r("DT", d.dt);
    cr("RADT", "radt"); cr("BLDT", "bldt"); cr("CUDT", "cudt");
    ci("AER_OPT", "aer_opt"); ci("SWINT_OPT", "swint_opt"); ci("AER_TYPE", "aer_type");
    ci("AER_AOD550_OPT", "aer_aod550_opt"); ci("AER_ANGEXP_OPT", "aer_angexp_opt");
    ci("AER_SSA_OPT", "aer_ssa_opt"); ci("AER_ASY_OPT", "aer_asy_opt");
    cr("AER_AOD550_VAL", "aer_aod550_val"); cr("AER_ANGEXP_VAL", "aer_angexp_val");
    cr("AER_SSA_VAL", "aer_ssa_val"); cr("AER_ASY_VAL", "aer_asy_val");
    ci("MOIST_ADV_OPT", "moist_adv_opt"); ci("SCALAR_ADV_OPT", "scalar_adv_opt"); ci("TKE_ADV_OPT", "tke_adv_opt");
    ci("DIFF_6TH_OPT", "diff_6th_opt"); cr("DIFF_6TH_FACTOR", "diff_6th_factor");
    if (c.integer("topo_wind") == 1) ci("TOPO_WIND", "topo_wind");
    ci("OBS_NUDGE_OPT", "obs_nudge_opt");
    cr("BUCKET_MM", "bucket_mm"); cr("BUCKET_J", "bucket_j"); cr("PREC_ACC_DT", "prec_acc_dt");
    ci("ISFTCFLX", "isftcflx"); ci("ISHALLOW", "ishallow"); ci("ISFFLX", "isfflx");
    ci("ICLOUD", "icloud"); ci("ICLOUD_CU", "icloud_cu"); ci("TRACER_PBLMIX", "tracer_pblmix");
    ci("SCALAR_PBLMIX", "scalar_pblmix"); ci("YSU_TOPDOWN_PBLMIX", "ysu_topdown_pblmix");
    ci("GRAV_SETTLING", "grav_settling"); ci("CLDOVRLP", "cldovrlp"); ci("IDCOR", "idcor");
    if (c.integer("sf_surface_physics") == 4)
      for (const char* o : {"opt_sfc", "dveg", "opt_crs", "opt_btr", "opt_runsub", "opt_runsrf", "opt_frz", "opt_inf",
                            "opt_rad", "opt_alb", "opt_snf", "opt_tbot", "opt_stc", "opt_gla", "opt_rsf", "opt_soil",
                            "opt_pedo", "opt_crop", "opt_irr", "opt_irrm", "opt_compact", "opt_tksno", "opt_scf",
                            "opt_wetland", "opt_tdrn", "opt_infdv"}) {
        std::string n = o;
        for (char& ch : n) ch = char(std::toupper(static_cast<unsigned char>(ch)));
        ci(n, o);
      }
    ci("DFI_OPT", "dfi_opt");
    // serial run; the ideal build keeps the namelist's nproc_x = nproc_y = -1
    i("NTASKS_X", d.ideal ? -1 : 1); i("NTASKS_Y", d.ideal ? -1 : 1); i("NTASKS_TOTAL", 1);
  }
  if (input || history || restart)
    nc.att("SIMULATION_INITIALIZATION_TYPE", std::string(d.ideal ? "IDEALIZED DATA" : "REAL-DATA CASE"));
  const int ide = d.ide, jde = d.jde, kde = d.kde;
  i("WEST-EAST_PATCH_START_UNSTAG", 1); i("WEST-EAST_PATCH_END_UNSTAG", ide - 1);
  i("WEST-EAST_PATCH_START_STAG", 1); i("WEST-EAST_PATCH_END_STAG", ide);
  i("SOUTH-NORTH_PATCH_START_UNSTAG", 1); i("SOUTH-NORTH_PATCH_END_UNSTAG", jde - 1);
  i("SOUTH-NORTH_PATCH_START_STAG", 1); i("SOUTH-NORTH_PATCH_END_STAG", jde);
  i("BOTTOM-TOP_PATCH_START_UNSTAG", 1); i("BOTTOM-TOP_PATCH_END_UNSTAG", kde - 1);
  i("BOTTOM-TOP_PATCH_START_STAG", 1); i("BOTTOM-TOP_PATCH_END_STAG", kde);
  i("GRID_ID", d.grid_id); i("PARENT_ID", d.parent_id); i("I_PARENT_START", d.i_parent_start);
  i("J_PARENT_START", d.j_parent_start); i("PARENT_GRID_RATIO", d.parent_grid_ratio);
  r("DT", d.dt);  // already defined for history: the value is rewritten in place
  for (const char* a : {"CEN_LAT", "CEN_LON", "TRUELAT1", "TRUELAT2", "MOAD_CEN_LAT", "STAND_LON", "POLE_LAT", "POLE_LON"}) {
    auto it = d.proj.find(a);
    r(a, it == d.proj.end() ? 0.0 : it->second);
  }
  if (s != Stream::boundary) {
    int jy, jd;
    double gm;
    julgmt(d.gmt_date.empty() ? d.start_date : d.gmt_date, jy, jd, gm);  // config_flags%gmt/julyr/julday: this domain's start
    if (d.ideal) jy = 0;  // output_wrf: nl_get_julyr, which only real.exe's input sets
    r("GMT", gm); i("JULYR", jy); i("JULDAY", jd);
  }
  const int mp = d.map_proj;
  i("MAP_PROJ", mp);
  nc.att("MAP_PROJ_CHAR", std::string(mp == 0 ? "Cartesian" : mp == 1 ? "Lambert Conformal" : mp == 2 ? "Polar Stereographic"
                                     : mp == 3 ? "Mercator" : mp == 6 ? "Cylindrical Equidistant" : "unknown"));
  std::string mminlu = d.mminlu.substr(0, 4) == "USGS" ? std::string("USGS") : d.mminlu;
  nc.att("MMINLU", d.ideal ? std::string() : mminlu.empty() || mminlu[0] == ' ' ? std::string("    ") : mminlu);
  i("NUM_LAND_CAT", d.num_land_cat); i("ISWATER", d.iswater); i("ISLAKE", d.islake);
  i("ISICE", d.isice); i("ISURBAN", d.isurban); i("ISOILWATER", d.isoilwater);
  ci("HYBRID_OPT", "hybrid_opt"); r("ETAC", c.integer("hybrid_opt") == 0 ? 0.0 : c.num("etac"));
}

}  // namespace cppwrf::wrf_io
