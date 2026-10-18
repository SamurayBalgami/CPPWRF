#pragma once
// WRF real.exe for one domain: the met_em input (input_auxinput1 / optional_input), init_domain_rk
// (dyn_em/module_initialize_real.F) and the boundary arrays of assemble_output (main/real_em.F,
// share/module_bc.F). T is WRF's default REAL kind (float, or double for the -r8 build).
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "init_real/farray.hpp"
#include "init_real/met_em.hpp"
#include "init_real/wrf_io_names.hpp"
#include "init_real/real_core.hpp"

namespace cppwrf::init_real {

// namelist.input options used by real.exe (WRF Registry defaults)
struct RealOptions {
  // time_control
  std::string start_date, end_date;  // yyyy-mm-dd_hh:mm:ss
  int interval_seconds = 21600;
  // domains
  int e_we = 0, e_sn = 0, e_vert = 0, num_metgrid_levels = 0, num_metgrid_soil_levels = 4;
  double dx = 0, dy = 0, p_top_requested = 5000;
  std::vector<double> eta_levels;
  int auto_levels_opt = 2;
  double max_dz = 1000., dzbot = 50., dzstretch_s = 1.3, dzstretch_u = 1.1;
  int interp_type = 2, lagrange_order = 2, linear_interp = 1, extrap_type = 2, t_extrap_type = 2;
  bool lowest_lev_from_sfc = false, use_levels_below_ground = true, use_surface = true;
  double zap_close_levels = 500;
  int force_sfc_in_vinterp = 1;
  bool interp_theta = false, sfcp_to_sfcp = false, adjust_heights = false, smooth_cg_topo = false;
  int rebalance = 0, hypsometric_opt = 2;
  bool rh2qv_wrt_liquid = true, use_sh_qv = false;
  int rh2qv_method = 1;
  double qv_max_p_safe = 10000, qv_max_flag = 1.E-5, qv_max_value = 3.E-6;
  double qv_min_p_safe = 110000, qv_min_flag = 1.E-6, qv_min_value = 1.E-6;
  int use_maxw_level = 0, use_trop_level = 0;
  // dynamics
  int hybrid_opt = 2, use_theta_m = 1;
  double etac = 0.2, base_pres = 100000., base_temp = 290., base_lapse = 50., iso_temp = 200.,
         base_pres_strat = 0., base_lapse_strat = -11.;
  // physics
  int sf_surface_physics = 0, num_soil_layers = 5, surface_input_source = 3, sf_urban_physics = 0;
  int mp_physics = 0, fractional_seaice = 0, default_soiltype = 8;
  double seaice_threshold = 100;
  bool use_tavg_for_tsk = false;
  // bdy_control
  int spec_bdy_width = 5;

  int dom = 1, max_dom = 1, parent_id = 0, i_parent_start = 1, j_parent_start = 1, parent_grid_ratio = 1;
  static RealOptions from_namelist(const Namelist& nl, int dom = 1);
};

std::string add_seconds(const std::string& date, long seconds);  // geth_newdate
void get_julgmt(const std::string& date, int& julyr, int& julday, double& gmt);

template <typename T>
class RealDomain {
 public:
  RealDomain(const RealOptions& o);

  // input_auxinput1 + optional_input for one time
  void read_input(const MetEm& m);
  // init_domain_rk; loop = internal_time_loop (1 = initial time)
  void init_domain(int loop, const std::string& current_date);

  const RealOptions o;
  Dims d{};
  int nlev = 0;  // num_metgrid_levels
  int iswater = 17, islake = -1, isice = 15, isurban = 13, isoilwater = 14;
  std::string mminlu;
  std::map<std::string, int> flags;  // FLAG_* of the input
  int num_veg_cat = 0, num_soil_top_cat = 0, num_soil_bot_cat = 0;

  // ---- constants of the base state
  T p00{}, t00{}, a{}, tiso{}, p_strat{}, a_strat{};
  T p_top{}, p_top_save{};
  T cf1{}, cf2{}, cf3{}, cfn{}, cfn1{}, rdx{}, rdy{};
  int lakeflag = 0, ifndsnowh = 0, ifndsoilw = 0, step_number = 0;
  bool got_var_sso = false;

  // ---- vertical coordinate (1..kde)
  F1<T> znw, znu, dnw, rdnw, dn, rdn, fnm, fnp, c1f, c2f, c3f, c4f, c1h, c2h, c3h, c4h;
  // ---- 3D model fields (i, k, j)
  F3<T> u_1, u_2, v_1, v_2, t_2, t_init, ph_2, ph0, phb, p, pb, al, alb, alt, p_hyd, php, qv, th_phy_m_t0;
  // ---- 2D model fields
  F2<T> mu_2, mub, mu0, psfc, ht, tsk, tmn, xlat, xlong, xlat_u, xlong_u, xlat_v, xlong_v, clat, clong, landmask,
      xland, lu_index, vegcat, soilcat, snow, snowh, snowc, canwat, vegfra, shdmax, shdmin, shdavg, albbck, snoalb,
      lai, u10, v10, t2, th2, q2, sst, xice, lakemask, msft, msfu, msfv, msftx, msfty, msfux, msfuy, msfvx, msfvy,
      msfvx_inv, f, e, sina, cosa, tsk_save, var_sso, var2d, oc12d, oa1, oa2, oa3, oa4, ol1, ol2, ol3, ol4,
      water_depth, pdhs, toposoil, sct_dom_gc, scb_dom_gc, frc_urb2d;
  F2<int> ivgtyp, isltyp;
  // ---- metgrid fields (i, 1..nlev, j)
  F3<T> p_gc, t_gc, rh_gc, u_gc, v_gc, ght_gc, qv_gc, pd_gc;
  F2<T> psfc_gc, pslv_gc, ht_gc, intq_gc, tsk_gc, tmn_gc, xlat_gc, xlong_gc;
  F3<T> greenfrac, albedo12m, lai12m, landusef, soilctop, soilcbot;  // (i, month/category, j)
  // ---- soil
  F3<T> tslb, smois, sh2o;  // (i, 1..num_soil_layers, j)
  F1<T> zs, dzs;
  // soil input of optional_input: (i, 1..num_levels+1, j), level 1 = surface (skin)
  F3<T> st_input, sm_input, sw_input;
  std::vector<int> st_levels_input, sm_levels_input, sw_levels_input;
  int num_st_levels_input = 0, num_sm_levels_input = 0, num_sw_levels_input = 0;

 private:
  int flag(const std::string& n) const {
    auto it = flags.find(n);
    return it == flags.end() ? 0 : it->second;
  }
};

}  // namespace cppwrf::init_real
