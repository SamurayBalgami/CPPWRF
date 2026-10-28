#include "init_real/real_domain.hpp"
#include "wrf_io/wrf_schema.hpp"

#include <cstdio>
#include <iostream>
#include <stdexcept>

namespace cppwrf::init_real {

// ---------------------------------------------------------------- options and dates
RealOptions RealOptions::from_namelist(const Namelist& nl, const WrfSchema& sc, int dom) {
  RealOptions o;
  o.dom = dom;
  // WRF's namelist read: a max_domains option without a value for this domain takes the Registry default,
  // other options their first value (wrf_schema entries)
  WrfConfig c(nl, sc.defaults(), dom, &sc.entries());
  o.max_dom = c.integer("max_dom");
  auto date = [&](const char* p) {
    char b[32];
    const std::string q(p);
    std::snprintf(b, sizeof b, "%04d-%02d-%02d_%02d:%02d:%02d", c.integer(q + "_year"), c.integer(q + "_month"),
                  c.integer(q + "_day"), c.integer(q + "_hour"), c.integer(q + "_minute"), c.integer(q + "_second"));
    return std::string(b);
  };
  o.start_date = date("start");
  o.end_date = date("end");
  o.interval_seconds = c.integer("interval_seconds");
  o.e_we = c.integer("e_we");
  o.e_sn = c.integer("e_sn");
  o.e_vert = c.integer("e_vert");
  o.num_metgrid_levels = c.integer("num_metgrid_levels");
  o.num_metgrid_soil_levels = c.integer("num_metgrid_soil_levels");
  // check_nml_consistency (module_check_a_mundo.F:263-278): a nest's dx / dy are dx(1) / dy(1) over the product of
  // the parent grid ratios up to domain 1 (get_moad_factor), whatever the namelist gives for the nest
  int factor = 1;
  for (int id = dom, n = 0; id > 1 && n < o.max_dom; id = c.integer("parent_id", id), ++n)
    factor *= c.integer("parent_grid_ratio", id);
  o.dx = c.num("dx", 1) / double(factor);
  o.dy = c.num("dy", 1) / double(factor);
  o.p_top_requested = c.num("p_top_requested");
  for (const std::string& v : nl.values("domains", "eta_levels")) o.eta_levels.push_back(std::strtod(v.c_str(), nullptr));
  o.auto_levels_opt = c.integer("auto_levels_opt");
  o.max_dz = c.num("max_dz");
  o.dzbot = c.num("dzbot");
  o.dzstretch_s = c.num("dzstretch_s");
  o.dzstretch_u = c.num("dzstretch_u");
  o.interp_type = c.integer("interp_type");
  o.lagrange_order = c.integer("lagrange_order");
  o.linear_interp = c.integer("linear_interp");
  o.extrap_type = c.integer("extrap_type");
  o.t_extrap_type = c.integer("t_extrap_type");
  o.lowest_lev_from_sfc = c.logical("lowest_lev_from_sfc");
  o.use_levels_below_ground = c.logical("use_levels_below_ground");
  o.use_surface = c.logical("use_surface");
  o.zap_close_levels = c.num("zap_close_levels");
  o.force_sfc_in_vinterp = c.integer("force_sfc_in_vinterp");
  // check_a_mundo (use_wps_input = 1): without the surface level, neither force_sfc_in_vinterp nor
  // lowest_lev_from_sfc
  if (!o.use_surface && o.force_sfc_in_vinterp > 0) o.force_sfc_in_vinterp = 0;
  if (!o.use_surface && o.lowest_lev_from_sfc) o.lowest_lev_from_sfc = false;
  o.interp_theta = c.logical("interp_theta");
  o.sfcp_to_sfcp = c.logical("sfcp_to_sfcp");
  o.adjust_heights = c.logical("adjust_heights");
  o.smooth_cg_topo = c.logical("smooth_cg_topo");
  o.rebalance = c.integer("rebalance");
  o.hypsometric_opt = c.integer("hypsometric_opt");
  o.rh2qv_wrt_liquid = c.logical("rh2qv_wrt_liquid");
  o.rh2qv_method = c.integer("rh2qv_method");
  o.use_sh_qv = c.logical("use_sh_qv");
  o.qv_max_p_safe = c.num("qv_max_p_safe");
  o.qv_max_flag = c.num("qv_max_flag");
  o.qv_max_value = c.num("qv_max_value");
  o.qv_min_p_safe = c.num("qv_min_p_safe");
  o.qv_min_flag = c.num("qv_min_flag");
  o.qv_min_value = c.num("qv_min_value");
  o.use_maxw_level = c.integer("use_maxw_level");
  o.use_trop_level = c.integer("use_trop_level");
  o.hybrid_opt = c.integer("hybrid_opt");
  o.use_theta_m = c.integer("use_theta_m");
  o.etac = c.num("etac");
  o.base_pres = c.num("base_pres");
  o.base_temp = c.num("base_temp");
  o.base_lapse = c.num("base_lapse");
  o.iso_temp = c.num("iso_temp");
  o.base_pres_strat = c.num("base_pres_strat");
  o.base_lapse_strat = c.num("base_lapse_strat");
  // setup_physics_suite, then set_physics_rconfigs: num_soil_layers from domain 1's sf_surface_physics
  const std::string suite = c.get("physics_suite");
  o.sf_surface_physics = wrf_io::physics_suite_option(suite, "sf_surface_physics", c.integer("sf_surface_physics"));
  o.mp_physics = wrf_io::physics_suite_option(suite, "mp_physics", c.integer("mp_physics"));
  const int lsm1 = wrf_io::physics_suite_option(suite, "sf_surface_physics", c.integer("sf_surface_physics", 1));
  if (lsm1 == 0 || lsm1 == 1) o.num_soil_layers = 5;
  else if (lsm1 == 2 || lsm1 == 4) o.num_soil_layers = 4;
  else throw std::runtime_error("real: sf_surface_physics = " + std::to_string(lsm1) +
                                " is not supported (WRF: unknown or not ported land surface scheme)");
  if (o.mp_physics < 0) throw std::runtime_error("real: mp_physics = -1 needs a physics_suite (CONUS or tropical)");
  o.surface_input_source = c.integer("surface_input_source");
  o.sf_urban_physics = c.integer("sf_urban_physics");
  o.fractional_seaice = c.integer("fractional_seaice");
  o.default_soiltype = c.integer("default_soiltype");
  o.sst_update = c.integer("sst_update");
  o.grid_fdda = c.integer("grid_fdda");
  o.qna_update = c.integer("qna_update");
  o.all_ic_times = c.logical("all_ic_times");
  o.polar = c.logical("polar");
  o.seaice_threshold = c.num("seaice_threshold");
  o.use_tavg_for_tsk = c.logical("use_tavg_for_tsk");
  o.spec_bdy_width = c.integer("spec_bdy_width");
  o.input_from_file = c.logical("input_from_file");
  o.parent_id = c.integer("parent_id");
  o.i_parent_start = c.integer("i_parent_start");
  o.j_parent_start = c.integer("j_parent_start");
  o.parent_grid_ratio = c.integer("parent_grid_ratio");
  if (o.e_we <= 0 || o.e_sn <= 0 || o.e_vert <= 0 || o.num_metgrid_levels <= 0)
    throw std::runtime_error("real: e_we, e_sn, e_vert and num_metgrid_levels are required in &domains");
  // options of real.exe that are not ported
  for (const char* k : {"use_maxw_level", "use_trop_level", "madwrf_opt", "madwrf_cldinit"})
    if (c.integer(k) != 0) throw std::runtime_error(std::string("real: ") + k + " is not ported");
  for (const char* k : {"insert_init_cloud", "all_ic_times", "multi_bdy_files", "slucm_distributed_drag"})
    if (c.logical(k)) throw std::runtime_error(std::string("real: ") + k + " = .true. is not ported");
  return o;
}

static bool leap(int y) { return (y % 4 == 0 && y % 100 != 0) || y % 400 == 0; }

void get_julgmt(const std::string& date, int& julyr, int& julday, double& gmt) {
  int y = std::atoi(date.substr(0, 4).c_str()), m = std::atoi(date.substr(5, 2).c_str()),
      dd = std::atoi(date.substr(8, 2).c_str()), h = std::atoi(date.substr(11, 2).c_str()),
      mi = std::atoi(date.substr(14, 2).c_str()), sec = std::atoi(date.substr(17, 2).c_str());
  // module_date_time.F: DATA MMD (implicitly SAVE) gets MMD(2) = 29 at the first leap-year date and keeps it
  static int mmd[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  if (leap(y)) mmd[1] = 29;
  julday = dd;
  for (int k = 1; k < m; ++k) julday += mmd[k - 1];
  julyr = y;
  gmt = double(h) + double(mi) / 60. + double(sec) / 3600.;
}

std::string add_seconds(const std::string& date, long seconds) {
  int y = std::atoi(date.substr(0, 4).c_str()), m = std::atoi(date.substr(5, 2).c_str()),
      d = std::atoi(date.substr(8, 2).c_str()), h = std::atoi(date.substr(11, 2).c_str()),
      mi = std::atoi(date.substr(14, 2).c_str()), s = std::atoi(date.substr(17, 2).c_str());
  long t = ((long(h) * 60 + mi) * 60 + s) + seconds;
  long days = t / 86400;
  t %= 86400;
  if (t < 0) { t += 86400; --days; }
  static const int mdays[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  d += int(days);
  auto ml = [&](int yy, int mm) { return mdays[mm - 1] + (mm == 2 && leap(yy) ? 1 : 0); };
  while (d > ml(y, m)) { d -= ml(y, m); if (++m > 12) { m = 1; ++y; } }
  while (d < 1) { if (--m < 1) { m = 12; --y; } d += ml(y, m); }
  char b[32];
  std::snprintf(b, sizeof b, "%04d-%02d-%02d_%02ld:%02ld:%02ld", y, m, d, t / 3600, (t / 60) % 60, t % 60);
  return b;
}

// ---------------------------------------------------------------- domain
template <typename T>
RealDomain<T>::RealDomain(const RealOptions& opt) : o(opt) {
  d = Dims{1, o.e_we, 1, o.e_sn, 1, o.e_vert};
  nlev = o.num_metgrid_levels;
  const int ide = d.ide, jde = d.jde, kde = d.kde;
  for (F1<T>* a : {&znw, &znu, &dnw, &rdnw, &dn, &rdn, &fnm, &fnp, &c1f, &c2f, &c3f, &c4f, &c1h, &c2h, &c3h, &c4h})
    *a = F1<T>(1, kde);
  for (F3<T>* a : {&u_1, &u_2, &v_1, &v_2, &t_2, &t_init, &ph_2, &ph0, &phb, &p, &pb, &al, &alb, &alt, &p_hyd, &php, &qv,
                   &th_phy_m_t0})
    *a = F3<T>(1, ide, 1, kde, 1, jde);
  for (F2<T>* a : {&mu_2, &mub, &mu0, &psfc, &ht, &tsk, &tmn, &xlat, &xlong, &xlat_u, &xlong_u, &xlat_v, &xlong_v, &clat,
                   &clong, &landmask, &xland, &lu_index, &vegcat, &soilcat, &snow, &snowh, &snowc, &canwat, &vegfra,
                   &shdmax, &shdmin, &shdavg, &albbck, &snoalb, &lai, &u10, &v10, &t2, &th2, &q2, &sst, &xice,
                   &lakemask, &msft, &msfu, &msfv, &msftx, &msfty, &msfux, &msfuy, &msfvx, &msfvy, &msfvx_inv, &f, &e,
                   &sina, &cosa, &tsk_save, &var_sso, &var2d, &oc12d, &oa1, &oa2, &oa3, &oa4, &ol1, &ol2, &ol3, &ol4,
                   &water_depth, &pdhs, &toposoil, &sct_dom_gc, &scb_dom_gc, &frc_urb2d, &psfc_gc, &pslv_gc, &ht_gc,
                   &intq_gc, &tsk_gc, &tmn_gc, &xlat_gc, &xlong_gc})
    *a = F2<T>(1, ide, 1, jde);
  ivgtyp = F2<int>(1, ide, 1, jde);
  isltyp = F2<int>(1, ide, 1, jde);
  for (F3<T>* a : {&p_gc, &t_gc, &rh_gc, &u_gc, &v_gc, &ght_gc, &qv_gc, &pd_gc}) *a = F3<T>(1, ide, 1, nlev, 1, jde);
  const int ns = std::max(1, o.num_soil_layers);
  tslb = F3<T>(1, ide, 1, ns, 1, jde);
  smois = F3<T>(1, ide, 1, ns, 1, jde);
  sh2o = F3<T>(1, ide, 1, ns, 1, jde);
  zs = F1<T>(1, ns);
  dzs = F1<T>(1, ns);
}

namespace {
// met_em variable -> F2 (mass or staggered), F3 (levels / categories in the middle index)
template <typename T>
void put2(const MetEm& m, const char* name, F2<T>& f) {
  const MetVar& v = m.var(name);
  if (v.shape.size() != 2) throw std::runtime_error(std::string("met_em: ") + name + " is not 2-D");
  const std::size_t ny = v.shape[0], nx = v.shape[1];
  for (std::size_t j = 0; j < ny; ++j)
    for (std::size_t i = 0; i < nx; ++i) f(int(i) + 1, int(j) + 1) = T(v.v[j * nx + i]);
}
template <typename T>
void put3(const MetEm& m, const char* name, F3<T>& f, int k0 = 1) {
  const MetVar& v = m.var(name);
  if (v.shape.size() != 3) throw std::runtime_error(std::string("met_em: ") + name + " is not 3-D");
  const std::size_t nk = v.shape[0], ny = v.shape[1], nx = v.shape[2];
  for (std::size_t k = 0; k < nk; ++k)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t i = 0; i < nx; ++i) f(int(i) + 1, int(k) + k0, int(j) + 1) = T(v.v[(k * ny + j) * nx + i]);
}
}  // namespace

template <typename T>
void RealDomain<T>::read_input(const MetEm& m) {
  const int ide = d.ide, jde = d.jde;
  flags.clear();
  for (const auto& [name, v] : m.int_atts())
    if (name.rfind("FLAG_", 0) == 0 && !v.empty()) flags[name.substr(5)] = v[0];
  if (m.att_int("NUM_METGRID_LEVELS", nlev) != nlev)
    throw std::runtime_error("real: num_metgrid_levels in the namelist does not match " + m.path);
  if (m.att_int("WEST-EAST_GRID_DIMENSION", ide) != ide || m.att_int("SOUTH-NORTH_GRID_DIMENSION", jde) != jde)
    throw std::runtime_error("real: e_we / e_sn do not match " + m.path);
  iswater = m.att_int("ISWATER", 16);
  islake = m.att_int("ISLAKE", -1);
  isice = m.att_int("ISICE", 24);
  isurban = m.att_int("ISURBAN", 1);
  isoilwater = m.att_int("ISOILWATER", 14);
  mminlu = m.att_str("MMINLU", "USGS");
  while (!mminlu.empty() && (mminlu.back() == ' ' || mminlu.back() == '\0')) mminlu.pop_back();
  // 3-D metgrid fields
  put3(m, "PRES", p_gc); put3(m, "TT", t_gc); put3(m, "RH", rh_gc); put3(m, "UU", u_gc); put3(m, "VV", v_gc);
  put3(m, "GHT", ght_gc);
  if (m.has("QV")) put3(m, "QV", qv_gc);
  // 2-D fields (Registry i1 names)
  put2(m, "PSFC", psfc); put2(m, "PMSL", pslv_gc); put2(m, "HGT_M", ht_gc); put2(m, "SKINTEMP", tsk_gc);
  put2(m, "SOILTEMP", tmn_gc); put2(m, "XLAT_M", xlat_gc); put2(m, "XLONG_M", xlong_gc);
  put2(m, "XLAT_U", xlat_u); put2(m, "XLONG_U", xlong_u); put2(m, "XLAT_V", xlat_v); put2(m, "XLONG_V", xlong_v);
  put2(m, "CLAT", clat); put2(m, "CLONG", clong);
  put2(m, "MAPFAC_M", msft); put2(m, "MAPFAC_U", msfu); put2(m, "MAPFAC_V", msfv);
  put2(m, "MAPFAC_MX", msftx); put2(m, "MAPFAC_MY", msfty); put2(m, "MAPFAC_UX", msfux); put2(m, "MAPFAC_UY", msfuy);
  put2(m, "MAPFAC_VX", msfvx); put2(m, "MAPFAC_VY", msfvy);
  put2(m, "F", f); put2(m, "E", e); put2(m, "SINALPHA", sina); put2(m, "COSALPHA", cosa);
  put2(m, "LANDMASK", landmask); put2(m, "LU_INDEX", lu_index);
  if (m.has("SOILHGT")) put2(m, "SOILHGT", toposoil);
  if (m.has("SNOW")) put2(m, "SNOW", snow);
  if (m.has("SNOWH")) put2(m, "SNOWH", snowh);
  if (m.has("SEAICE")) put2(m, "SEAICE", xice);
  if (m.has("SST")) put2(m, "SST", sst);  // Registry SST: input stream i1 (used with FLAG_SST)
  // optional input-stream fields WRF passes from met_em to wrfinput, not carried by cppwrf_real
  for (const char* n : {"UOCE", "VOCE", "TOPOSLPX", "TOPOSLPY", "EROD", "ICEDEPTH", "ALBSI", "SNOWSI", "LANDUSEF2",
                        "SOILT1", "o3_gfs_du"})
    if (m.has(n)) throw std::runtime_error(std::string("real: met_em field ") + n + " is not ported");
  if (m.has("SNOALB")) put2(m, "SNOALB", snoalb);
  if (m.has("SCT_DOM")) put2(m, "SCT_DOM", sct_dom_gc);
  if (m.has("SCB_DOM")) put2(m, "SCB_DOM", scb_dom_gc);
  if (m.has("VAR")) put2(m, "VAR", var2d);
  if (m.has("CON")) put2(m, "CON", oc12d);
  for (auto [n, fp] : {std::pair<const char*, F2<T>*>{"OA1", &oa1}, {"OA2", &oa2}, {"OA3", &oa3}, {"OA4", &oa4},
                       {"OL1", &ol1}, {"OL2", &ol2}, {"OL3", &ol3}, {"OL4", &ol4}})
    if (m.has(n)) put2(m, n, *fp);
  // Registry dataname "var_sso" (lower case; NetCDF names are case-sensitive, so WRF does not read VAR_SSO)
  if (m.has("var_sso")) put2(m, "var_sso", var_sso);
  // optional_input: soil levels (module_optional_input.F optional_lsm_levels, FLAG_SOIL_LAYERS path)
  num_st_levels_input = num_sm_levels_input = num_sw_levels_input = 0;
  if (flags.count("SOIL_LEVELS") && flags["SOIL_LEVELS"] == 1)
    throw std::runtime_error("real: soil-level (RUC) met_em input is not ported");
  if (flag("UM_SOIL") == 1) throw std::runtime_error("real: UM soil input (FLAG_UM_SOIL) is not ported");
  if (m.has("SW")) throw std::runtime_error("real: met_em soil water (SW) input is not ported");
  if (o.sf_surface_physics == 4 && mminlu == "NLCD40")  // module_initialize_real: fatal with Noah-MP
    throw std::runtime_error("real: Noah-MP does not support the NLCD40 land use");
  if (flag("SOIL_LAYERS") == 1)
    num_st_levels_input = num_sm_levels_input = num_sw_levels_input = o.num_metgrid_soil_levels;
  const int sfs = o.sf_surface_physics;
  if (sfs == 1 || sfs == 2 || sfs == 3 || sfs == 4 || sfs == 5 || sfs == 7 || sfs == 8 || sfs == 88) {
    const int nalloc = o.num_soil_layers * 3;
    if (st_input.empty()) {
      st_input = F3<T>(1, ide, 1, nalloc, 1, jde);
      sm_input = F3<T>(1, ide, 1, nalloc, 1, jde);
      sw_input = F3<T>(1, ide, 1, nalloc, 1, jde);
    }
    if (flag("SOIL_LAYERS") != 1) throw std::runtime_error("real: legacy soil input (no FLAG_SOIL_LAYERS) is not ported");
    const int n = num_st_levels_input;
    F3<T> sl(1, ide, 1, n, 1, jde), st(1, ide, 1, n, 1, jde), sm(1, ide, 1, n, 1, jde);
    put3(m, "SOIL_LAYERS", sl);
    put3(m, "ST", st);
    put3(m, "SM", sm);
    st_levels_input.assign(std::size_t(n), -1);
    sm_levels_input.assign(std::size_t(n), -1);
    sw_levels_input.assign(std::size_t(n), -1);
    int level_above = 0;
    for (int k = 1; k <= n; ++k) {
      // INTEGER = (INTEGER + REAL) / 2 (real division, truncated)
      T v = (T(level_above) + sl(1, n + 1 - k, 1)) / T(2);
      st_levels_input[std::size_t(k - 1)] = sm_levels_input[std::size_t(k - 1)] = sw_levels_input[std::size_t(k - 1)] = int(v);
      level_above = int(sl(1, n + 1 - k, 1));
    }
    for (int j = 1; j <= jde - 1; ++j)
      for (int k = 1; k <= n; ++k)
        for (int i = 1; i <= ide - 1; ++i) {
          st_input(i, k + 1, j) = st(i, n + 1 - k, j);
          sm_input(i, k + 1, j) = sm(i, n + 1 - k, j);
          sw_input(i, k + 1, j) = T(0);  // grid%sw: not in the metgrid output
        }
    if (n > nalloc) throw std::runtime_error("real: soil_levels_too_few");
  }
  // monthly and category fields
  auto put_cat = [&](const char* name, F3<T>& fld, int& ncat) {
    if (!m.has(name)) { ncat = 0; return; }
    ncat = int(m.var(name).shape[0]);
    fld = F3<T>(1, ide, 1, ncat, 1, jde);
    put3(m, name, fld);
  };
  int n12 = 0;
  put_cat("GREENFRAC", greenfrac, n12);
  put_cat("ALBEDO12M", albedo12m, n12);
  put_cat("LAI12M", lai12m, n12);
  put_cat("LANDUSEF", landusef, num_veg_cat);
  put_cat("SOILCTOP", soilctop, num_soil_top_cat);
  put_cat("SOILCBOT", soilcbot, num_soil_bot_cat);
}

// ---------------------------------------------------------------- init_domain_rk
template <typename T>
void RealDomain<T>::init_domain(int loop, const std::string& current_date) {
  using C = Const<T>;
  const int ids = d.ids, ide = d.ide, jds = d.jds, jde = d.jde, kds = d.kds, kde = d.kde;
  const int its = d.its(), ite = d.ite(), jts = d.jts(), jte = d.jte(), kts = d.kts(), kte = d.kte();
  const T g = C::g, r_d = C::r_d, cp = C::cp, t0 = C::t0, p1000mb = C::p1000mb, cvpm = C::cvpm, cpovcv = C::cpovcv,
          rvovrd = C::rvovrd;
  const int sfs = o.sf_surface_physics;
  if (sfs != 0 && sfs != 2 && sfs != 4)
    throw std::runtime_error("real: only sf_surface_physics = 0, 2 (Noah) and 4 (Noah-MP) are ported");
  const bool lsm = (sfs == 2 || sfs == 4);  // LSMSCHEME, NOAHMPSCHEME: module_initialize_real treats them alike
  SkipCtx sk;
  sk.width = o.spec_bdy_width;
  sk.hold_ups = (loop == 1) || o.grid_fdda != 0 || o.sst_update == 1 || o.qna_update == 1 || o.all_ic_times || o.polar;
  auto skip = [&](int i, int j) { return skip_middle(d, i, j, sk.width, sk.hold_ups); };
  const int i1 = std::min(ide - 1, ite), j1 = std::min(jde - 1, jte);

  // use_maxw_level / use_trop_level = 0: the max-wind and tropopause levels are not used
  for (int j = jts; j <= j1; ++j)
    for (int i = its; i <= i1; ++i) lakemask(i, j) = (int(lu_index(i, j)) != islake) ? T(0) : T(1);
  if (flag("BATHYMETRY") == 0)
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) water_depth(i, j) = T(-4.0);
  else
    throw std::runtime_error("real: bathymetry input is not ported");
  step_number = 0;
  // const_module_initialize
  p00 = T(o.base_pres); t00 = T(o.base_temp); a = T(o.base_lapse); tiso = T(o.iso_temp);
  p_strat = T(o.base_pres_strat); a_strat = T(o.base_lapse_strat);
  if (flag("EXCLUDED_MIDDLE") != 0) throw std::runtime_error("real: excluded-middle met_em input is not ported");
  bool any_valid = false;
  int i_valid = 0, j_valid = 0;
  for (int j = jts; j <= jte && !any_valid; ++j)
    for (int i = its; i <= ite; ++i) {
      if (skip(i, j)) continue;
      any_valid = true; i_valid = i; j_valid = j;
      break;
    }
  if (flag("ICEFRAC") == 1 || flag("ICEPCT") == 1) throw std::runtime_error("real: ICEFRAC/ICEPCT input is not ported");
  if (flag("SNOW") == 0 && flag("SNOWH") == 0) {
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        snow(i, j) = T(0.); snowh(i, j) = T(0.);
      }
  } else if (flag("SNOW") == 0 && flag("SNOWH") == 1) {
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        snow(i, j) = snowh(i, j) * T(1000.) / T(5.);
      }
  } else if (flag("SNOW") == 1 && flag("SNOWH") == 0) {
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        snowh(i, j) = snow(i, j) / T(1000.) * T(5.);
      }
  }
  if (flag("MF_XY") == 1) {
    if (msfvx(its, jts) == T(0)) throw std::runtime_error("real: msfvx = 0 (global domain without polar)");
    for (int j = jts; j <= std::min(jde, jte); ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        msfvx_inv(i, j) = T(1.) / msfvx(i, j);
      }
  } else {
    throw std::runtime_error("real: met_em without FLAG_MF_XY is not supported");
  }
  const bool we_have_tavgsfc = flag("TAVGSFC") == 1;
  if (o.use_tavg_for_tsk) throw std::runtime_error("real: use_tavg_for_tsk is not ported");
  if (flag("METGRID") != 1) throw std::runtime_error("real: input must come from metgrid (FLAG_METGRID = 1)");

  // ================================ vertical interpolation part
  if (flag("PINTERP") == 1) throw std::runtime_error("real: P_INTERP input is not ported");
  for (int j = jts; j <= std::min(jte, jde - 1); ++j)
    for (int i = its; i <= std::min(ite, ide - 1); ++i) {
      if (skip(i, j)) continue;
      tsk(i, j) = tsk_gc(i, j); tmn(i, j) = tmn_gc(i, j); xlat(i, j) = xlat_gc(i, j); xlong(i, j) = xlong_gc(i, j);
      ht(i, j) = ht_gc(i, j);
    }
  if (o.smooth_cg_topo) throw std::runtime_error("real: smooth_cg_topo is not ported");
  if (flag("PSFC") == 1)
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        psfc_gc(i, j) = psfc(i, j);
        p_gc(i, 1, j) = psfc(i, j);
      }
  if (flag("SOILHGT") == 1)
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        ght_gc(i, 1, j) = toposoil(i, j);
        ht_gc(i, j) = toposoil(i, j);
      }
  if (flag("PTHETA") == 1 || flag("PRHO") == 1) throw std::runtime_error("real: PTHETA/PRHO input is not ported");
  bool sfcp_to_sfcp = o.sfcp_to_sfcp, use_surface = o.use_surface, use_sh_qv = o.use_sh_qv;
  QvLimits ql{o.rh2qv_wrt_liquid, o.qv_max_p_safe, o.qv_max_flag, o.qv_max_value,
              o.qv_min_p_safe, o.qv_min_flag, o.qv_min_value};
  if (o.rh2qv_method != 1) throw std::runtime_error("real: only rh2qv_method = 1 is ported");
  if (any_valid) {
    int k = (p_gc(i_valid, nlev, j_valid) < p_gc(i_valid, 2, j_valid)) ? 2 : nlev;
    if (t_gc(i_valid, 1, j_valid) == T(-1.E30) || rh_gc(i_valid, 1, j_valid) == T(-1.E30) ||
        u_gc(i_valid, 1, j_valid) == T(-1.E30) || v_gc(i_valid, 1, j_valid) == T(-1.E30))
      throw std::runtime_error("real: missing surface fields (use_surface fallback) are not ported");
    (void)k;
    if (flag("QV") != 1 && flag("SH") != 1) {
      use_sh_qv = false;
      Dims dm{ids, ide, jds, jde, 1, nlev};
      rh_to_mxrat1(rh_gc, t_gc, p_gc, qv_gc, ql, nlev, dm, sk);
    } else {
      throw std::runtime_error("real: QV / SH metgrid input is not ported");
    }
    if (ght_gc(i_valid, nlev / 2, j_valid) < T(1)) throw std::runtime_error("real: GHT reconstruction is not ported");
    if (p_gc(i_valid, nlev / 2, j_valid) < p_gc(i_valid, nlev / 2 + 1, j_valid)) sfcp_to_sfcp = true;
  }
  for (int j = jts; j <= std::min(jde - 1, jte); ++j)
    for (int i = its; i <= std::min(ide, ite); ++i) {
      if (skip(i, j)) continue;
      u10(i, j) = u_gc(i, 1, j);
    }
  for (int j = jts; j <= std::min(jde, jte); ++j)
    for (int i = its; i <= std::min(ide - 1, ite); ++i) {
      if (skip(i, j)) continue;
      v10(i, j) = v_gc(i, 1, j);
    }
  for (int j = jts; j <= j1; ++j)
    for (int i = its; i <= i1; ++i) {
      if (skip(i, j)) continue;
      t2(i, j) = t_gc(i, 1, j);
    }
  {
    Dims dm{ids, ide, jds, jde, 1, nlev};
    T pt = find_p_top(p_gc, nlev, dm, sk);
    if (loop == 1) {
      if (T(o.p_top_requested) < pt)
        throw std::runtime_error("real: p_top_requested < p_top possible from data");
      p_top = T(o.p_top_requested);
      p_top_save = p_top;
    } else {
      if (pt > p_top_save) throw std::runtime_error("real: p_top > previous value");
      p_top = p_top_save;
    }
  }
  // monthly_interp_to_date / monthly_min_max / monthly_avg
  {
    int middle[14];
    for (int l = 1; l <= 12; ++l) {
      char b[32];
      std::snprintf(b, sizeof b, "%s-%02d-15_00:00:00", current_date.substr(0, 4).c_str(), l);
      int jy, jd; double gm;
      get_julgmt(b, jy, jd, gm);
      middle[l] = jy * 1000 + jd;
    }
    middle[0] = middle[1] - 31;
    middle[13] = middle[12] + 31;
    int ty, td; double gm;
    get_julgmt(current_date, ty, td, gm);
    const int target_date = ty * 1000 + td;
    auto interp = [&](const F3<T>& fin, F2<T>& fout) {
      if (fin.empty()) return;
      for (int l = 0; l <= 12; ++l) {
        if (middle[l] < target_date && middle[l + 1] >= target_date) {
          for (int j = jts; j <= j1; ++j)
            for (int i = its; i <= i1; ++i) {
              if (skip(i, j)) continue;
              int m1, m2;
              if (l == 0 || l == 12) { m1 = 12; m2 = 1; } else { m1 = l; m2 = l + 1; }
              fout(i, j) = (fin(i, m2, j) * T(target_date - middle[l]) + fin(i, m1, j) * T(middle[l + 1] - target_date)) /
                           T(middle[l + 1] - middle[l]);
            }
          break;
        }
      }
    };
    interp(greenfrac, vegfra);
    interp(albedo12m, albbck);
    interp(lai12m, lai);
    if (!greenfrac.empty())
      for (int j = jts; j <= j1; ++j)
        for (int i = its; i <= i1; ++i) {
          T mn = greenfrac(i, 1, j), mx = greenfrac(i, 1, j);
          for (int l = 2; l <= 12; ++l) {
            if (greenfrac(i, l, j) < mn) mn = greenfrac(i, l, j);
            if (greenfrac(i, l, j) > mx) mx = greenfrac(i, l, j);
          }
          shdmin(i, j) = mn; shdmax(i, j) = mx;
          T s = T(0);
          for (int l = 1; l <= 12; ++l) s = s + greenfrac(i, l, j);
          shdavg(i, j) = s / T(12);
        }
  }
  for (int j = jts; j <= j1; ++j)
    for (int i = its; i <= i1; ++i) {
      if (skip(i, j)) continue;
      vegfra(i, j) = vegfra(i, j) * T(100.); shdmax(i, j) = shdmax(i, j) * T(100.);
      shdmin(i, j) = shdmin(i, j) * T(100.); shdavg(i, j) = shdavg(i, j) * T(100.);
    }
  for (int j = jts; j <= j1; ++j)
    for (int i = its; i <= i1; ++i) {
      if (skip(i, j)) continue;
      albbck(i, j) = albbck(i, j) / T(100.);
      snoalb(i, j) = snoalb(i, j) / T(100.);
      if (landmask(i, j) < T(0.5)) { albbck(i, j) = T(0.08); snoalb(i, j) = T(0.08); }
    }
  {
    Dims dm{ids, ide, jds, jde, 1, nlev};
    if (flag("PSFC") == 1 && flag("SOILHGT") == 1 && flag("SLP") == 1 && !sfcp_to_sfcp) {
      sfcprs3(ght_gc, p_gc, ht, pslv_gc, psfc, nlev, dm, sk);
    } else {
      throw std::runtime_error("real: only the sfcprs3 surface pressure path is ported");
    }
    integ_moist(qv_gc, p_gc, pd_gc, t_gc, ght_gc, intq_gc, nlev, dm, sk);
  }
  p_dts(mu0, intq_gc, psfc, p_top, d, sk);
  p_dhs(pdhs, ht, p00, t00, a, d, sk);
  if (znw(1) != T(1.0)) {
    EtaOpts eo;
    bool have_eta = !o.eta_levels.empty() && o.eta_levels[0] != -1.0;
    if (have_eta) eo.eta_levels = o.eta_levels;
    eo.auto_levels_opt = o.auto_levels_opt;
    eo.max_dz = o.max_dz; eo.dzbot = o.dzbot; eo.dzstretch_s = o.dzstretch_s; eo.dzstretch_u = o.dzstretch_u;
    compute_eta(znw, eo, p_top, kde);
  }
  compute_vcoord_1d_coeffs(T(o.etac), znw, o.hybrid_opt, p_top, kde, znu, c1f, c2f, c3f, c4f, c1h, c2h, c3h, c4h);
  if (o.interp_theta) throw std::runtime_error("real: interp_theta is not ported");
  if (flag("SLP") == 1) {
    p_dry(mu0, znw, p_top, pb, true, c3f, c3h, c4f, c4h, d, sk);
    VInterpOpts vo;
    vo.interp_type = 2; vo.lagrange_order = o.lagrange_order; vo.lowest_lev_from_sfc = false;
    vo.use_levels_below_ground = true; vo.use_surface = true; vo.zap_close_levels = o.zap_close_levels;
    vo.force_sfc_in_vinterp = 0; vo.extrap_type = 1;
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        psfc_gc(i, j) = pd_gc(i, 1, j);
        pd_gc(i, 1, j) = pslv_gc(i, j) - (p_gc(i, 1, j) - pd_gc(i, 1, j));
        ht_gc(i, j) = ght_gc(i, 1, j);
        ght_gc(i, 1, j) = T(0.);
      }
    vert_interp(ght_gc, pd_gc, ph0, pb, nlev, 'Z', vo, d, sk);
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        pd_gc(i, 1, j) = psfc_gc(i, j);
        ght_gc(i, 1, j) = ht_gc(i, j);
      }
  }
  p_dry(mu0, znw, p_top, pb, false, c3f, c3h, c4f, c4h, d, sk);
  VInterpOpts vo;
  vo.interp_type = o.interp_type; vo.lagrange_order = o.lagrange_order; vo.lowest_lev_from_sfc = o.lowest_lev_from_sfc;
  vo.use_levels_below_ground = o.use_levels_below_ground; vo.use_surface = use_surface;
  vo.zap_close_levels = o.zap_close_levels; vo.force_sfc_in_vinterp = o.force_sfc_in_vinterp;
  vo.extrap_type = o.extrap_type;
  vert_interp(rh_gc, pd_gc, u_1, pb, nlev, 'Q', vo, d, sk);
  if (use_sh_qv) throw std::runtime_error("real: use_sh_qv is not ported");
  {
    VInterpOpts vt = vo;
    vt.interp_type = 2;
    vt.extrap_type = o.t_extrap_type;
    vert_interp(t_gc, pd_gc, t_2, pb, nlev, 'T', vt, d, sk);
    vt.interp_type = 1;
    vert_interp(p_gc, pd_gc, p, pb, nlev, 'T', vt, d, sk);
  }
  v_1 = t_2;
  rh_to_mxrat1(u_1, v_1, p, qv, ql, kte - 1, d, sk);
  t_to_theta(t_2, p, p00, kte - 1, d, sk);
  for (const char* q : {"QR", "QC", "QI", "QS", "QG", "QH", "QNI", "QNC", "QNR", "QNS", "QNG", "QNH", "SPECCLDL",
                        "SPECCLDF"})
    if (flag(q) == 1) throw std::runtime_error(std::string("real: hydrometeor input (FLAG_") + q + ") is not ported");
  vert_interp(u_gc, pd_gc, u_2, pb, nlev, 'U', vo, d, sk);
  vert_interp(v_gc, pd_gc, v_2, pb, nlev, 'V', vo, d, sk);
  // ================================ end of the vertical interpolation part

  if (islake < 0) throw std::runtime_error("real: input without inland lake information is not ported");
  lakeflag = 1;
  if (we_have_tavgsfc) throw std::runtime_error("real: TAVGSFC with inland lakes is not ported");
  for (int j = jts; j <= j1; ++j)
    for (int i = its; i <= i1; ++i) {
      if (skip(i, j)) continue;
      landusef(i, iswater, j) = landusef(i, iswater, j) + landusef(i, islake, j);
      landusef(i, islake, j) = T(0.);
    }
  if (o.surface_input_source == 3)
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        if (lu_index(i, j) == T(islake)) lu_index(i, j) = T(iswater);
      }
  for (int j = jts; j <= j1; ++j)
    for (int i = its; i <= i1; ++i) {
      if (skip(i, j)) continue;
      tsk_save(i, j) = tsk(i, j);
    }
  const int flag_sst = flag("SST");
  for (int j = jts; j <= j1; ++j)
    for (int i = its; i <= i1; ++i) {
      if (skip(i, j)) continue;
      if (landmask(i, j) < T(0.5) && flag_sst == 1 && sst(i, j) > T(170.) && sst(i, j) < T(400.)) tsk(i, j) = sst(i, j);
    }
  for (int j = jts; j <= j1; ++j)
    for (int i = its; i <= i1; ++i) {
      if (skip(i, j)) continue;
      snowc(i, j) = snow(i, j) >= T(10.) ? T(1.) : T(0.0);
    }
  ifndsnowh = flag("SNOWH");
  ifndsoilw = num_sw_levels_input >= 1 ? 1 : 0;
  if (lsm) {
    if (num_st_levels_input < 2) throw std::runtime_error("real: not enough soil temperature data for Noah LSM");
    // process_soil_real: fix_bottom_level_for_temp has no action for LSMSCHEME
    // adjust_soil_temp_new
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        if (landmask(i, j) > T(0.5)) tmn(i, j) = tmn(i, j) - T(0.0065) * ht(i, j);
      }
    if (flag("SOILHGT") == 1)
      for (int j = jts; j <= j1; ++j)
        for (int i = its; i <= i1; ++i) {
          if (skip(i, j)) continue;
          T vmin = toposoil(i, j), vmax = toposoil(i, j);
          T dmin = ht(i, j) - toposoil(i, j), dmax = ht(i, j) - toposoil(i, j);
          if (vmin < T(-1000)) continue;
          if (vmax > T(10000)) continue;
          if (dmin < T(-3000) || dmax > T(3000)) continue;
          if (landmask(i, j) > T(0.5)) {
            tsk(i, j) = tsk(i, j) - T(0.0065) * (ht(i, j) - toposoil(i, j));
            for (int k = 2; k <= num_st_levels_input + 1; ++k)
              st_input(i, k, j) = st_input(i, k, j) - T(0.0065) * (ht(i, j) - toposoil(i, j));
          }
        }
    if (o.num_soil_layers > 1) {
      // init_soil_depth_2
      if (o.num_soil_layers != 4) throw std::runtime_error("real: the Noah LSM uses num_soil_layers = 4");
      const T dz4[4] = {T(0.1), T(0.3), T(0.6), T(1.0)};
      for (int l = 1; l <= 4; ++l) dzs(l) = dz4[l - 1];
      zs(1) = T(.5) * dzs(1);
      for (int l = 2; l <= o.num_soil_layers; ++l) zs(l) = zs(l - 1) + T(.5) * dzs(l - 1) + T(.5) * dzs(l);
      // init_soil_2_real (flag_soil_layers = 1)
      const int nst = num_st_levels_input, nsm = num_sm_levels_input, nsw = num_sw_levels_input;
      if (nst * nsm < 1) throw std::runtime_error("real: no input soil level data");
      auto sort_levels = [&](std::vector<int>& lev, F3<T>& in, int n) {
        for (int lout = 1; lout <= n - 1; ++lout)
          for (int lin = lout + 1; lin <= n; ++lin)
            if (lev[std::size_t(lout - 1)] > lev[std::size_t(lin - 1)]) {
              T temp = T(lev[std::size_t(lout - 1)]);
              lev[std::size_t(lout - 1)] = lev[std::size_t(lin - 1)];
              lev[std::size_t(lin - 1)] = int(std::lround(double(temp)));
              for (int j = jts; j <= j1; ++j)
                for (int i = its; i <= i1; ++i) {
                  if (skip(i, j)) continue;
                  T t2v = in(i, lout + 1, j);
                  in(i, lout + 1, j) = in(i, lin + 1, j);
                  in(i, lin + 1, j) = t2v;
                }
            }
      };
      sort_levels(st_levels_input, st_input, nst);
      for (int j = jts; j <= j1; ++j)
        for (int i = its; i <= i1; ++i) {
          if (skip(i, j)) continue;
          st_input(i, 1, j) = tsk(i, j);
          st_input(i, nst + 2, j) = tmn(i, j);
        }
      sort_levels(sm_levels_input, sm_input, nsm);
      for (int j = jts; j <= j1; ++j)
        for (int i = its; i <= i1; ++i) {
          if (skip(i, j)) continue;
          sm_input(i, 1, j) = sm_input(i, 2, j);
          sm_input(i, nsm + 2, j) = sm_input(i, nsm + 1, j);
        }
      sort_levels(sw_levels_input, sw_input, nsw);
      if (nsw > 1)
        for (int j = jts; j <= j1; ++j)
          for (int i = its; i <= i1; ++i) {
            if (skip(i, j)) continue;
            sw_input(i, 1, j) = sw_input(i, 2, j);
            sw_input(i, nsw + 2, j) = sw_input(i, nsw + 1, j);
          }
      std::vector<T> zhave(std::size_t(std::max({nst, nsm, nsw}) + 2 + 1));
      auto interp_layers = [&](const std::vector<int>& lev, const F3<T>& in, F3<T>& out, int n) {
        zhave[1] = T(0.);
        for (int l = 1; l <= n; ++l) zhave[std::size_t(l + 1)] = T(lev[std::size_t(l - 1)]) / T(100.);
        zhave[std::size_t(n + 2)] = T(300.) / T(100.);
        for (int lwant = 1; lwant <= o.num_soil_layers; ++lwant)
          for (int lhave = 1; lhave <= n + 2 - 1; ++lhave)
            if (zs(lwant) >= zhave[std::size_t(lhave)] && zs(lwant) <= zhave[std::size_t(lhave + 1)]) {
              for (int j = jts; j <= j1; ++j)
                for (int i = its; i <= i1; ++i) {
                  if (skip(i, j)) continue;
                  out(i, lwant, j) = (in(i, lhave, j) * (zhave[std::size_t(lhave + 1)] - zs(lwant)) +
                                      in(i, lhave + 1, j) * (zs(lwant) - zhave[std::size_t(lhave)])) /
                                     (zhave[std::size_t(lhave + 1)] - zhave[std::size_t(lhave)]);
                }
              break;
            }
      };
      interp_layers(st_levels_input, st_input, tslb, nst);
      interp_layers(sm_levels_input, sm_input, smois, nsm);
      if (nsw > 1) interp_layers(sw_levels_input, sw_input, sh2o, nsw);
      for (int j = jts; j <= j1; ++j)
        for (int i = its; i <= i1; ++i) {
          if (skip(i, j)) continue;
          if (landmask(i, j) < T(0.5))
            for (int l = 1; l <= o.num_soil_layers; ++l) {
              tslb(i, l, j) = flag_sst == 1 ? sst(i, j) : tsk(i, j);
              smois(i, l, j) = T(1.0);
              sh2o(i, l, j) = T(1.0);
            }
        }
    }
  }
  if (any_valid) {
    if (o.surface_input_source != 3) throw std::runtime_error("real: only surface_input_source = 3 is ported");
    if (sct_dom_gc(i_valid, j_valid) > T(0.5)) {
      for (int j = jts; j <= j1; ++j)
        for (int i = its; i <= i1; ++i) {
          if (skip(i, j)) continue;
          isltyp(i, j) = int(std::lround(double(sct_dom_gc(i, j))));
          soilcat(i, j) = T(isltyp(i, j));
        }
    } else {
      throw std::runtime_error("real: surface_input_source = 3 needs dominant soil data");
    }
    if (lu_index(i_valid, j_valid) > T(0.5)) {
      for (int j = jts; j <= j1; ++j)
        for (int i = its; i <= i1; ++i) {
          if (skip(i, j)) continue;
          ivgtyp(i, j) = int(std::lround(double(lu_index(i, j))));
          vegcat(i, j) = T(ivgtyp(i, j));
        }
    } else {
      throw std::runtime_error("real: surface_input_source = 3 needs dominant land use data");
    }
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        if (landmask(i, j) > T(0.5) && isltyp(i, j) == isoilwater) isltyp(i, j) = o.default_soiltype;
        else if (landmask(i, j) < T(0.5) && isltyp(i, j) != isoilwater) isltyp(i, j) = isoilwater;
      }
    if (o.sf_urban_physics != 0) throw std::runtime_error("real: urban physics input is not ported");
  }
  // adjust_for_seaice_pre (LSMSCHEME: no sea ice on land points; no-op for sf_surface_physics = 0)
  if (lsm)
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        if (landmask(i, j) > T(0.5)) xice(i, j) = T(0.);
      }
  for (int j = jts; j <= j1; ++j)
    for (int i = its; i <= i1; ++i) {
      if (skip(i, j)) continue;
      lu_index(i, j) = T(ivgtyp(i, j));
      if (int(lu_index(i, j)) != iswater) { landmask(i, j) = T(1); xland(i, j) = T(1); }
      else { landmask(i, j) = T(0); xland(i, j) = T(2); }
    }
  if (lsm)  // fix_tsk_tmn
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        if (landmask(i, j) < T(0.5) && flag_sst == 1 && sst(i, j) > T(170.) && sst(i, j) < T(400.)) {
          tmn(i, j) = sst(i, j);
          tsk(i, j) = sst(i, j);
        } else if (landmask(i, j) < T(0.5)) {
          tmn(i, j) = tsk(i, j);
        }
      }
  if (loop != 1) {
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        if (tsk(i, j) < T(170) || tsk(i, j) > T(400.)) tsk(i, j) = t_2(i, 1, j);
      }
  } else {
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        if (tsk(i, j) < T(170) || tsk(i, j) > T(400.)) {
          if (tmn(i, j) > T(170.) && tmn(i, j) < T(400.)) tsk(i, j) = tmn(i, j);
          else if (sst(i, j) > T(170.) && sst(i, j) < T(400.)) tsk(i, j) = sst(i, j);
          else throw std::runtime_error("real: tsk unreasonable");
        }
      }
  }
  for (int j = jts; j <= j1; ++j)
    for (int i = its; i <= i1; ++i) {
      if (skip(i, j)) continue;
      if ((tmn(i, j) < T(170.) || tmn(i, j) > T(400.)) && landmask(i, j) > T(0.5)) {
        if (tsk(i, j) > T(170.) && tsk(i, j) < T(400.)) tmn(i, j) = tsk(i, j);
        else if (sst(i, j) > T(170.) && sst(i, j) < T(400.)) tmn(i, j) = sst(i, j);
        else throw std::runtime_error("real: tmn unreasonable");
      }
    }
  // account_for_zero_soil_moisture (only at the start time)
  if (lsm && current_date == o.start_date && flag("SOIL_LAYERS") == 1)
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        if (landmask(i, j) > T(0.5) && tslb(i, 1, j) > T(170) && tslb(i, 1, j) < T(400) && smois(i, 1, j) < T(0.005))
          for (int l = 1; l <= o.num_soil_layers; ++l) smois(i, l, j) = T(0.005);
      }
  // soil temperature checks (sf_surface_physics = 0: tslb stays 0 and WRF resets smois(1:4) to 0.3 on land)
  if (loop != 1) {
    for (int j = jts; j <= j1; ++j)
      for (int ns = 1; ns <= o.num_soil_layers; ++ns)
        for (int i = its; i <= i1; ++i) {
          if (skip(i, j)) continue;
          if (tslb(i, ns, j) < T(170) || tslb(i, ns, j) > T(400.)) {
            tslb(i, ns, j) = t_2(i, 1, j);
            smois(i, ns, j) = T(0.3);
          }
        }
  } else {
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        if ((tslb(i, 1, j) < T(170.) || tslb(i, 1, j) > T(400.)) && landmask(i, j) > T(0.5)) {
          if (!lsm)
            for (int l = 1; l <= 4 && l <= o.num_soil_layers; ++l) smois(i, l, j) = T(0.3);
          if (tsk(i, j) > T(170.) && tsk(i, j) < T(400.) && tmn(i, j) > T(170.) && tmn(i, j) < T(400.)) {
            if (lsm)  // fake_soil_temp
              for (int ns = 1; ns <= o.num_soil_layers; ++ns)
                tslb(i, ns, j) = (tsk(i, j) * (T(3.0) - zs(ns)) + tmn(i, j) * (T(0.0) - zs(ns))) / (T(3.0) - T(0.0));
          } else {
            throw std::runtime_error("real: tslb unreasonable");
          }
        }
      }
  }
  // adjust_for_seaice_post (LSMSCHEME; no-op for sf_surface_physics = 0)
  if (lsm) {
    const T xice_threshold = o.fractional_seaice == 0 ? T(0.5) : T(0.02);
    const T seaice_threshold = T(o.seaice_threshold);
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        if (xice(i, j) > T(200.0)) xice(i, j) = T(0.);
      }
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        if ((tsk(i, j) < T(170) || tsk(i, j) > T(400)) && tsk_save(i, j) > T(170) && tsk_save(i, j) < T(400))
          tsk(i, j) = tsk_save(i, j);
        if ((tsk(i, j) < T(170) || tsk(i, j) > T(400)) && (tsk_save(i, j) < T(170) || tsk_save(i, j) > T(400))) {
          throw std::runtime_error("real: TSK is unrealistic, problems for seaice post");
        } else if (xice(i, j) >= xice_threshold || (landmask(i, j) < T(0.5) && tsk(i, j) < seaice_threshold)) {
          if (o.fractional_seaice == 0) xice(i, j) = T(1.0);
          if (landmask(i, j) < T(0.5)) tmn(i, j) = T(271.4);
          vegcat(i, j) = T(isice); ivgtyp(i, j) = isice; lu_index(i, j) = T(isice);
          landmask(i, j) = T(1.); xland(i, j) = T(1.); vegfra(i, j) = T(0.);
          for (int l = 1; l <= num_veg_cat; ++l) landusef(i, l, j) = T(0.);
          landusef(i, ivgtyp(i, j), j) = T(1.);
          tsk_save(i, j) = tsk(i, j);
          isltyp(i, j) = 16;
          soilcat(i, j) = T(isltyp(i, j));
          for (int l = 1; l <= num_soil_top_cat; ++l) soilctop(i, l, j) = T(0);
          for (int l = 1; l <= num_soil_bot_cat; ++l) soilcbot(i, l, j) = T(0);
          soilctop(i, isltyp(i, j), j) = T(1.);
          soilcbot(i, isltyp(i, j), j) = T(1.);
          const T total_depth = T(3.);
          for (int l = 1; l <= o.num_soil_layers; ++l) {
            T mid = (total_depth / T(o.num_soil_layers)) / T(2.) + T(l - 1) * (total_depth / T(o.num_soil_layers));
            tslb(i, l, j) = ((total_depth - mid) * tsk(i, j) + mid * tmn(i, j)) / total_depth;
          }
          for (int l = 1; l <= o.num_soil_layers; ++l) { smois(i, l, j) = T(1.0); sh2o(i, l, j) = T(0.0); }
        } else if (xice(i, j) < xice_threshold) {
          xice(i, j) = T(0.);
        }
      }
  }
  for (int j = jts; j <= j1; ++j)
    for (int i = its; i <= i1; ++i) {
      if (skip(i, j)) continue;
      bool bad = (landmask(i, j) < T(0.5) && (ivgtyp(i, j) != iswater || isltyp(i, j) != 14)) ||
                 (landmask(i, j) > T(0.5) && (ivgtyp(i, j) == iswater || isltyp(i, j) == 14));
      if (!bad) continue;
      if (tslb(i, 1, j) > T(1.)) {
        ivgtyp(i, j) = 5; isltyp(i, j) = 8; landmask(i, j) = T(1); xland(i, j) = T(1);
      } else if (sst(i, j) > T(1.)) {
        ivgtyp(i, j) = iswater; isltyp(i, j) = 14; landmask(i, j) = T(0); xland(i, j) = T(2);
      } else {
        throw std::runtime_error("real: the landmask and soil/veg categories do not match");
      }
    }
  for (int j = jts; j <= j1; ++j)
    for (int i = its; i <= i1; ++i) {
      if (skip(i, j)) continue;
      if (flag_sst != 1) sst(i, j) = tsk(i, j);
    }
  for (int j = jts; j <= j1; ++j)
    for (int i = its; i <= i1; ++i) {
      if (skip(i, j)) continue;
      if (ivgtyp(i, j) == isice) snoalb(i, j) = T(0.75);
    }
  // znw checks (full levels, top-down order) are satisfied by compute_eta
  for (int k = 1; k <= kde - 1; ++k) {
    dnw(k) = znw(k + 1) - znw(k);
    rdnw(k) = T(1.) / dnw(k);
    znu(k) = T(0.5) * (znw(k + 1) + znw(k));
  }
  for (int k = 2; k <= kde - 1; ++k) {
    dn(k) = T(0.5) * (dnw(k) + dnw(k - 1));
    rdn(k) = T(1.) / dn(k);
    fnp(k) = T(.5) * dnw(k) / dn(k);
    fnm(k) = T(.5) * dnw(k - 1) / dn(k);
  }
  {
    T cof1 = (T(2.) * dn(2) + dn(3)) / (dn(2) + dn(3)) * dnw(1) / dn(2);
    T cof2 = dn(2) / (dn(2) + dn(3)) * dnw(1) / dn(3);
    cf1 = fnp(2) + cof1;
    cf2 = fnm(2) - cof1 - cof2;
    cf3 = cof2;
    cfn = (T(.5) * dnw(kde - 1) + dn(kde - 1)) / dn(kde - 1);
    cfn1 = T(-.5) * dnw(kde - 1) / dn(kde - 1);
  }
  rdx = T(1.) / T(o.dx);
  rdy = T(1.) / T(o.dy);
  for (int j = jts; j <= jte; ++j)
    for (int i = its; i <= ite; ++i) {
      if (skip(i, j)) continue;
      ph0(i, 1, j) = ht(i, j) * g;
      ph_2(i, 1, j) = T(0.);
    }
  for (int j = jts; j <= j1; ++j)
    for (int i = its; i <= i1; ++i) {
      if (skip(i, j)) continue;
      T p_surf = p00 * std::exp(-t00 / a + std::pow(sq(t00 / a) - T(2.) * g * ht(i, j) / a / r_d, T(0.5)));
      for (int k = 1; k <= kte - 1; ++k) {
        php(i, k, j) = c3f(k) * (p_surf - p_top) + c4f(k) + p_top;
        pb(i, k, j) = c3h(k) * (p_surf - p_top) + c4h(k) + p_top;
        T temp = std::max(tiso, t00 + a * std::log(pb(i, k, j) / p00));
        if (pb(i, k, j) < p_strat) temp = tiso + a_strat * std::log(pb(i, k, j) / p_strat);
        t_init(i, k, j) = temp * std::pow(p00 / pb(i, k, j), r_d / cp) - t0;
        alb(i, k, j) = (r_d / p1000mb) * (t_init(i, k, j) + t0) * std::pow(pb(i, k, j) / p1000mb, cvpm);
      }
      php(i, kte, j) = p_top;
      mub(i, j) = p_surf - p_top;
      phb(i, 1, j) = ht(i, j) * g;
      if (o.hypsometric_opt == 1) {
        for (int kk = 2; kk <= kte; ++kk) {
          int k = kk - 1;
          phb(i, kk, j) = phb(i, kk - 1, j) - dnw(kk - 1) * (c1h(k) * mub(i, j) + c2h(k)) * alb(i, kk - 1, j);
        }
      } else if (o.hypsometric_opt == 2) {
        for (int k = 2; k <= kte; ++k) {
          T pfu = c3f(k) * mub(i, j) + c4f(k) + p_top;
          T pfd = c3f(k - 1) * mub(i, j) + c4f(k - 1) + p_top;
          T phm = c3h(k - 1) * mub(i, j) + c4h(k - 1) + p_top;
          phb(i, k, j) = phb(i, k - 1, j) + alb(i, k - 1, j) * phm * std::log(pfd / pfu);
        }
      } else {
        throw std::runtime_error("real: hypsometric_opt should be 1 or 2");
      }
    }
  if (ite == ide) {
    int i = ide;
    for (int j = jts; j <= std::min(jde - 1, jte); ++j) {
      mub(i, j) = mub(i - 1, j);
      mu_2(i, j) = mu_2(i - 1, j);
      for (int k = 1; k <= kte - 1; ++k) {
        pb(i, k, j) = pb(i - 1, k, j); t_init(i, k, j) = t_init(i - 1, k, j); alb(i, k, j) = alb(i - 1, k, j);
      }
      for (int k = 1; k <= kte; ++k) phb(i, k, j) = phb(i - 1, k, j);
    }
  }
  if (jte == jde) {
    int j = jde;
    for (int i = its; i <= ite; ++i) {
      mub(i, j) = mub(i, j - 1);
      mu_2(i, j) = mu_2(i, j - 1);
      for (int k = 1; k <= kte - 1; ++k) {
        pb(i, k, j) = pb(i, k, j - 1); t_init(i, k, j) = t_init(i, k, j - 1); alb(i, k, j) = alb(i, k, j - 1);
      }
      for (int k = 1; k <= kte; ++k) phb(i, k, j) = phb(i, k, j - 1);
    }
  }
  for (int j = jts; j <= j1; ++j)
    for (int i = its; i <= i1; ++i) {
      if (skip(i, j)) continue;
      mu_2(i, j) = mu0(i, j) - mub(i, j);
    }
  if (ite == ide)
    for (int j = jts; j <= std::min(jde - 1, jte); ++j) mu_2(ide, j) = mu_2(ide - 1, j);
  if (jte == jde)
    for (int i = its; i <= ite; ++i) mu_2(i, jde) = mu_2(i, jde - 1);
  // hydrostatic balance (rebalance = 0, adjust_heights = .false.: one pass)
  if (o.adjust_heights) throw std::runtime_error("real: adjust_heights is not ported");
  for (int j = jts; j <= j1; ++j)
    for (int i = its; i <= i1; ++i) {
      if (skip(i, j)) continue;
      for (int k = 1; k <= kde - 1; ++k) t_2(i, k, j) = t_2(i, k, j) - t0;
      {
        int kk = kte - 1, k = kk + 1;
        T qtot = T(0.) + qv(i, kk, j);
        T qvf2 = T(1.) / (T(1.) + qtot);
        T qvf1 = qtot * qvf2;
        p(i, kk, j) = -T(0.5) * ((c1f(k) * mu_2(i, j)) + qvf1 * (c1f(k) * mub(i, j) + c2f(k))) / rdnw(kk) / qvf2;
        T qvf = T(1.) + rvovrd * qv(i, kk, j);
        alt(i, kk, j) = (r_d / p1000mb) * (t_2(i, kk, j) + t0) * qvf *
                        std::pow((p(i, kk, j) + pb(i, kk, j)) / p1000mb, cvpm);
        al(i, kk, j) = alt(i, kk, j) - alb(i, kk, j);
        p_hyd(i, kk, j) = p(i, kk, j) + pb(i, kk, j);
        for (kk = kte - 2; kk >= 1; --kk) {
          k = kk + 1;
          qtot = T(0.) + T(0.5) * (qv(i, kk, j) + qv(i, kk + 1, j));
          qvf2 = T(1.) / (T(1.) + qtot);
          qvf1 = qtot * qvf2;
          p(i, kk, j) = p(i, kk + 1, j) - ((c1f(k) * mu_2(i, j)) + qvf1 * (c1f(k) * mub(i, j) + c2f(k))) / qvf2 / rdn(kk + 1);
          qvf = T(1.) + rvovrd * qv(i, kk, j);
          alt(i, kk, j) = (r_d / p1000mb) * (t_2(i, kk, j) + t0) * qvf *
                          std::pow((p(i, kk, j) + pb(i, kk, j)) / p1000mb, cvpm);
          al(i, kk, j) = alt(i, kk, j) - alb(i, kk, j);
          p_hyd(i, kk, j) = p(i, kk, j) + pb(i, kk, j);
        }
      }
      if (o.hypsometric_opt == 1) {
        for (int kk = 2; kk <= kte; ++kk) {
          int k = kk - 1;
          ph_2(i, kk, j) = ph_2(i, kk - 1, j) -
                           dnw(kk - 1) * (((c1h(k) * mub(i, j) + c2h(k)) + (c1h(k) * mu_2(i, j))) * al(i, kk - 1, j) +
                                          (c1h(k) * mu_2(i, j)) * alb(i, kk - 1, j));
          ph0(i, kk, j) = ph_2(i, kk, j) + phb(i, kk, j);
        }
      } else {
        ph_2(i, 1, j) = phb(i, 1, j);
        for (int k = 2; k <= kte; ++k) {
          T pfu = c3f(k) * mu0(i, j) + c4f(k) + p_top;
          T pfd = c3f(k - 1) * mu0(i, j) + c4f(k - 1) + p_top;
          T phm = c3h(k - 1) * mu0(i, j) + c4h(k - 1) + p_top;
          ph_2(i, k, j) = ph_2(i, k - 1, j) + alt(i, k - 1, j) * phm * std::log(pfd / pfu);
        }
        for (int k = 1; k <= kte; ++k) ph_2(i, k, j) = ph_2(i, k, j) - phb(i, k, j);
      }
      if (o.hypsometric_opt == 1) {
        for (int k = kts; k <= kte - 1; ++k)
          al(i, k, j) = -T(1.) / ((c1h(k) * mub(i, j) + c2h(k)) + (c1h(k) * mu_2(i, j))) *
                        (alb(i, k, j) * (c1h(k) * mu_2(i, j)) + rdnw(k) * (ph_2(i, k + 1, j) - ph_2(i, k, j)));
      } else {
        for (int k = kts; k <= kte - 1; ++k) {
          T pfu = c3f(k + 1) * (mub(i, j) + mu_2(i, j)) + c4f(k + 1) + p_top;
          T pfd = c3f(k) * (mub(i, j) + mu_2(i, j)) + c4f(k) + p_top;
          T phm = c3h(k) * (mub(i, j) + mu_2(i, j)) + c4h(k) + p_top;
          al(i, k, j) = (ph_2(i, k + 1, j) - ph_2(i, k, j) + phb(i, k + 1, j) - phb(i, k, j)) / phm / std::log(pfd / pfu) -
                        alb(i, k, j);
        }
      }
      for (int k = kts; k <= kte - 1; ++k) {
        T qvf = T(1.) + rvovrd * qv(i, k, j);
        p(i, k, j) = p1000mb * std::pow((r_d * (t0 + t_2(i, k, j)) * qvf) / (p1000mb * (al(i, k, j) + alb(i, k, j))), cpovcv) -
                     pb(i, k, j);
        p_hyd(i, k, j) = p(i, k, j) + pb(i, k, j);
        alt(i, k, j) = al(i, k, j) + alb(i, k, j);
      }
    }
  for (int j = 1; j <= jde; ++j)
    for (int k = 1; k <= kde; ++k)
      for (int i = 1; i <= ide; ++i) v_1(i, k, j) = t_2(i, k, j) + t0;
  theta_to_t(v_1, p_hyd, p00, d, sk);
  if (!use_sh_qv) rh_to_mxrat1(u_1, v_1, p_hyd, qv, ql, kte - 1, d, sk);
  if (o.rebalance != 0) throw std::runtime_error("real: rebalance is not ported");
  for (int j = jts; j <= std::min(jde - 1, jte); ++j)
    for (int k = kts; k <= kte - 1; ++k)
      for (int i = its; i <= std::min(ide, ite); ++i) {
        T qvf = T(1.) + rvovrd * qv(i, k, j);
        p(i, k, j) = p1000mb * std::pow((r_d * (t0 + t_2(i, k, j)) * qvf) / (p1000mb * (al(i, k, j) + alb(i, k, j))), cpovcv) -
                     pb(i, k, j);
        p_hyd(i, k, j) = p(i, k, j) + pb(i, k, j);
      }
  for (int j = jts; j <= j1; ++j)
    for (int i = its; i <= i1; ++i) {
      if (skip(i, j)) continue;
      th2(i, j) = t2(i, j) * std::pow(p00 / (p(i, 1, j) + pb(i, 1, j)), r_d / cp);
    }
  if (flag("QV") != 1)
    for (int j = jts; j <= j1; ++j)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        q2(i, j) = qv_gc(i, 1, j);
      }
  // Thompson (mp_physics = 8): number concentrations are derived only where qc, qi, qr > 0 (no hydrometeor
  // input from metgrid here, so nothing to do)
  got_var_sso = false;
  for (int j = jts; j <= j1; ++j)
    for (int i = its; i <= i1; ++i) {
      if (skip(i, j)) continue;
      if (var_sso(i, j) != T(0)) got_var_sso = true;
    }
  for (int j = jts; j <= j1; ++j)
    for (int k = kts; k <= kte; ++k)
      for (int i = its; i <= i1; ++i) {
        if (skip(i, j)) continue;
        th_phy_m_t0(i, k, j) = t_2(i, k, j);
      }
  if (o.use_theta_m == 1)
    for (int j = jts; j <= j1; ++j)
      for (int k = kts; k <= kte; ++k)
        for (int i = its; i <= std::min(ide, ite); ++i) {
          if (skip(i, j)) continue;
          t_2(i, k, j) = (t_2(i, k, j) + t0) * (T(1.) + (C::r_v / r_d) * qv(i, k, j)) - t0;
        }
}

template class RealDomain<float>;
template class RealDomain<double>;

}  // namespace cppwrf::init_real
