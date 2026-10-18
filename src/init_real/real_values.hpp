#pragma once
// Values of the wrfinput variables real.exe sets (by WRF output name); every other input-stream
// variable keeps WRF's initial value 0.
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "init_real/met_em.hpp"
#include "init_real/nc_order.hpp"
#include "init_real/wrf_io_names.hpp"
#include "init_real/real_domain.hpp"

namespace cppwrf::init_real {

using ValueFn = std::function<std::vector<double>(const OutVar&)>;

template <typename T>
std::map<std::string, ValueFn> wrfinput_values(const RealDomain<T>& dm, const MetEm& m) {
  (void)m;
  std::map<std::string, ValueFn> v;
  auto f3 = [](const F3<T>* f) {
    return ValueFn([f](const OutVar& o) { return nc_order(*f, int(o.dimlens[2]), int(o.dimlens[0]), int(o.dimlens[1])); });
  };
  auto f2 = [](const F2<T>* f) {
    return ValueFn([f](const OutVar& o) { return nc_order(*f, int(o.dimlens[1]), int(o.dimlens[0])); });
  };
  auto f1 = [](const F1<T>* f) { return ValueFn([f](const OutVar& o) { return nc_order(*f, int(o.dimlens[0])); }); };
  auto sc = [](double x) { return ValueFn([x](const OutVar&) { return std::vector<double>{x}; }); };
  auto i2 = [](const F2<int>* f) {
    return ValueFn([f](const OutVar& o) {
      const int ni = int(o.dimlens[1]), nj = int(o.dimlens[0]);
      std::vector<double> r(static_cast<std::size_t>(ni) * nj);
      for (int j = 1; j <= nj; ++j)
        for (int i = 1; i <= ni; ++i) r[std::size_t(j - 1) * ni + (i - 1)] = (*f)(i, j);
      return r;
    });
  };
  // 3-D fields
  for (auto [n, f] : {std::pair<const char*, const F3<T>*>{"U", &dm.u_2}, {"V", &dm.v_2}, {"PH", &dm.ph_2},
                      {"PHB", &dm.phb}, {"T", &dm.th_phy_m_t0}, {"THM", &dm.t_2}, {"T_INIT", &dm.t_init},
                      {"P", &dm.p}, {"AL", &dm.al}, {"ALB", &dm.alb}, {"PB", &dm.pb}, {"P_HYD", &dm.p_hyd},
                      {"QVAPOR", &dm.qv}, {"TSLB", &dm.tslb}, {"SMOIS", &dm.smois}, {"SH2O", &dm.sh2o},
                      {"LANDUSEF", &dm.landusef}, {"SOILCTOP", &dm.soilctop}, {"SOILCBOT", &dm.soilcbot}})
    v[n] = f3(f);
  // 2-D fields
  for (auto [n, f] : {std::pair<const char*, const F2<T>*>{"MU", &dm.mu_2}, {"MUB", &dm.mub}, {"Q2", &dm.q2},
                      {"T2", &dm.t2}, {"TH2", &dm.th2}, {"PSFC", &dm.psfc}, {"U10", &dm.u10}, {"V10", &dm.v10},
                      {"XLAT", &dm.xlat}, {"XLONG", &dm.xlong}, {"LU_INDEX", &dm.lu_index}, {"SEAICE", &dm.xice},
                      {"VEGFRA", &dm.vegfra}, {"SNOW", &dm.snow}, {"SNOWH", &dm.snowh}, {"CANWAT", &dm.canwat},
                      {"WATER_DEPTH", &dm.water_depth}, {"LAI", &dm.lai}, {"VAR", &dm.var2d},
                      {"VAR_SSO", &dm.var_sso}, {"MAPFAC_M", &dm.msft}, {"MAPFAC_U", &dm.msfu},
                      {"MAPFAC_V", &dm.msfv}, {"MAPFAC_MX", &dm.msftx}, {"MAPFAC_MY", &dm.msfty},
                      {"MAPFAC_UX", &dm.msfux}, {"MAPFAC_UY", &dm.msfuy}, {"MAPFAC_VX", &dm.msfvx},
                      {"MAPFAC_VY", &dm.msfvy}, {"MF_VX_INV", &dm.msfvx_inv}, {"F", &dm.f}, {"E", &dm.e},
                      {"SINALPHA", &dm.sina}, {"COSALPHA", &dm.cosa}, {"HGT", &dm.ht}, {"TSK", &dm.tsk},
                      {"XLAT_U", &dm.xlat_u}, {"XLONG_U", &dm.xlong_u}, {"XLAT_V", &dm.xlat_v},
                      {"XLONG_V", &dm.xlong_v}, {"CLAT", &dm.clat}, {"ALBBCK", &dm.albbck}, {"TMN", &dm.tmn},
                      {"XLAND", &dm.xland}, {"SHDMAX", &dm.shdmax}, {"SHDMIN", &dm.shdmin}, {"SHDAVG", &dm.shdavg},
                      {"SNOALB", &dm.snoalb}, {"SST", &dm.sst}, {"SNOWC", &dm.snowc}, {"LAKEMASK", &dm.lakemask},
                      {"LANDMASK", &dm.landmask}, {"OA1", &dm.oa1}, {"OA2", &dm.oa2}, {"OA3", &dm.oa3},
                      {"OA4", &dm.oa4}, {"OL1", &dm.ol1}, {"OL2", &dm.ol2}, {"OL3", &dm.ol3}, {"OL4", &dm.ol4},
                      {"CON", &dm.oc12d}})
    v[n] = f2(f);
  v["IVGTYP"] = i2(&dm.ivgtyp);
  v["ISLTYP"] = i2(&dm.isltyp);
  // vertical coordinate
  for (auto [n, f] : {std::pair<const char*, const F1<T>*>{"ZNU", &dm.znu}, {"ZNW", &dm.znw}, {"FNM", &dm.fnm},
                      {"FNP", &dm.fnp}, {"RDNW", &dm.rdnw}, {"RDN", &dm.rdn}, {"DNW", &dm.dnw}, {"DN", &dm.dn},
                      {"C1H", &dm.c1h}, {"C2H", &dm.c2h}, {"C3H", &dm.c3h}, {"C4H", &dm.c4h}, {"C1F", &dm.c1f},
                      {"C2F", &dm.c2f}, {"C3F", &dm.c3f}, {"C4F", &dm.c4f}, {"ZS", &dm.zs}, {"DZS", &dm.dzs}})
    v[n] = f1(f);
  // scalars
  v["CFN"] = sc(double(dm.cfn)); v["CFN1"] = sc(double(dm.cfn1)); v["CF1"] = sc(double(dm.cf1));
  v["CF2"] = sc(double(dm.cf2)); v["CF3"] = sc(double(dm.cf3)); v["RDX"] = sc(double(dm.rdx));
  v["RDY"] = sc(double(dm.rdy)); v["P_TOP"] = sc(double(dm.p_top)); v["T00"] = sc(double(dm.t00));
  v["P00"] = sc(double(dm.p00)); v["TLP"] = sc(double(dm.a)); v["TISO"] = sc(double(dm.tiso));
  v["TLP_STRAT"] = sc(double(dm.a_strat)); v["P_STRAT"] = sc(double(dm.p_strat));
  v["FNDSNOWH"] = sc(dm.ifndsnowh); v["FNDSOILW"] = sc(dm.ifndsoilw); v["LAKEFLAG"] = sc(dm.lakeflag);
  v["SAVE_TOPO_FROM_REAL"] = sc(1); v["GOT_VAR_SSO"] = sc(dm.got_var_sso ? 1 : 0);
  v["STEP_NUMBER"] = sc(dm.step_number);
  return v;
}

}  // namespace cppwrf::init_real
