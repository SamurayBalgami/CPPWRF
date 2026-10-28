#include "phys/ra_rrtmg_sw/rrtmg_sw.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw.hpp"

#include <algorithm>
#include <cmath>

namespace cppwrf::phys {

namespace {
constexpr double rrsw_scon = 1.36822e+03;
}  // namespace

RrtmgSwCallConfig rrtmg_sw_call_config(int num_levels, const RrtmgSwConfig& config) {
  RrtmgSwCallConfig c{};
  c.num_levels = num_levels;
  c.nlayers = num_levels + 1;  // one buffer layer to the TOA (1e-5 mb)
  c.o3input = config.o3input;
  c.co2 = config.co2; c.n2o = config.n2o; c.ch4 = config.ch4;
  c.g = config.g;
  c.heatfac = config.heatfac;
  c.adjflux = 1.0 * (config.solcon / rrsw_scon);  // inatm_sw: adjflux = adjes * (scon / rrsw_scon)
  c.icld = config.cldovrlp; c.idcor = config.idcor;
  c.juldat = int(config.julian);  // juldat = julian (INTEGER)
  return c;
}

// WRF rrtmg_sw (cldprmc_sw, inatm_sw, setcoef_sw, spcvmc_sw, fluxes, heating rates) for one column whose
// McICA subcolumn cloud flags are in ws.cldfmc (the core check has no clouds: the paths are zero).
// Level arrays: nlayers + 1 (0 = surface); heating rates K/day.
void rrtmg_sw_core(int nlayers, double heatfac, const double* play, const double* plev, const double* tlay,
                   const double* tlev, const double* h2ovmr, const double* o3vmr, const double* co2vmr,
                   const double* ch4vmr, const double* n2ovmr, const double* o2vmr, double asdir, double asdif,
                   double aldir, double aldif, double coszen, double adjes, double scon, const double* reicmcl,
                   const double* relqmcl, const double* resnmcl, double* swuflx, double* swdflx, double* swhr,
                   double* swuflxc, double* swdflxc, double* swhrc, double* swdkdir, double* swdkdif,
                   double* sibvisdir, double* sibvisdif, double* sibnirdir, double* sibnirdif, RrtmgSwWorkspace& ws,
                   const RrtmgSwTables& tables, int inflag, int iceflag, int liqflag) {
  SwLayout lay(nlayers);
  ensure_workspace(ws, lay);
  const SwViews v = lay.views(ws.d.data(), ws.i.data(), 0, 1);
  v.pz[0] = plev[0];
  v.tz[0] = tlev[0];
  for (int l = 0; l < nlayers; ++l) {
    v.pavel[l] = play[l];
    v.tavel[l] = tlay[l];
    v.pz[l + 1] = plev[l + 1];
    v.tz[l + 1] = tlev[l + 1];
    v.wkl[0][l] = h2ovmr[l]; v.wkl[1][l] = co2vmr[l]; v.wkl[2][l] = o3vmr[l]; v.wkl[3][l] = n2ovmr[l];
    v.wkl[4][l] = 0.0; v.wkl[5][l] = ch4vmr[l]; v.wkl[6][l] = o2vmr[l];
    v.reicmc[l] = reicmcl[l]; v.relqmc[l] = relqmcl[l]; v.resnmc[l] = resnmcl[l];
    v.ciwpth[l] = 0.0; v.clwpth[l] = 0.0; v.cswpth[l] = 0.0;
    for (int g = 0; g < 112; ++g) v.cldfmc[l][g] = ws.cldfmc[l * 112 + g];
  }
  v.iscal[2] = inflag; v.iscal[3] = iceflag; v.iscal[4] = liqflag;
  rrtmg_sw_core_prep(nlayers, v, tables, coszen, asdir, asdif, aldir, aldif);
  RrtmgSwCallConfig c{};
  c.num_levels = nlayers - 1; c.nlayers = nlayers; c.heatfac = heatfac; c.adjflux = adjes * (scon / rrsw_scon);
  for (int g = 0; g < 112; ++g) rrtmg_sw_gpoint(c, g, v, tables);
  for (int lev = 0; lev <= nlayers; ++lev) rrtmg_sw_level(lev, v, tables);
  for (int i = 0; i < nlayers; ++i) {
    const double zdpgcp = heatfac / v.pdp[i];
    swhrc[i] = ((v.cd[i + 1] - v.cu[i + 1]) - (v.cd[i] - v.cu[i])) * zdpgcp;
    swhr[i] = ((v.fd[i + 1] - v.fu[i + 1]) - (v.fd[i] - v.fu[i])) * zdpgcp;
  }
  swhrc[nlayers - 1] = 0.;
  swhr[nlayers - 1] = 0.;
  for (int i = 0; i <= nlayers; ++i) {
    swuflxc[i] = v.cu[i];
    swdflxc[i] = v.cd[i];
    swuflx[i] = v.fu[i];
    swdflx[i] = v.fd[i];
    swdkdir[i] = v.fddir[i];
    swdkdif[i] = swdflx[i] - v.fddir[i];
    sibvisdir[i] = v.uvfddir[i];
    sibvisdif[i] = v.uvfd[i] - v.uvfddir[i];
    sibnirdir[i] = v.nifddir[i];
    sibnirdif[i] = v.nifd[i] - v.nifddir[i];
  }
}

void rrtmg_sw_column(
    int num_levels,
    const RrtmgSwConfig& config,
    const RrtmgSwTables& tables,
    double coszen,
    double albedo,
    double xland,
    double xice,
    double snow,
    double xlat,
    const double* p8w,
    const double* t8w,
    const double* p3d,
    const double* pi3d,
    const double* dz8w,
    const double* t3d,
    const double* rho3d,
    const double* qv3d,
    const double* qc3d,
    const double* qr3d,
    const double* qi3d,
    const double* qs3d,
    const double* qg3d,
    const double* cldfra3d,
    const double* o33d,
    SwColumnOut& out,
    double* rthratensw,
    double* rthratenswc,
    RrtmgSwWorkspace& ws,
    const double* re_cloud,
    const double* re_ice,
    const double* re_snow
) {
  (void)rho3d; (void)qr3d; (void)qg3d;
  const RrtmgSwCallConfig c = rrtmg_sw_call_config(num_levels, config);
  SwLayout lay(c.nlayers);
  ensure_workspace(ws, lay);
  const SwViews v = lay.views(ws.d.data(), ws.i.data(), 0, 1);
  std::vector<double> zero;
  if (!cldfra3d || config.icloud == 0) zero.assign(num_levels, 0.0);
  SwColumnIn in{};
  in.coszen = coszen; in.albedo = albedo; in.xland = xland; in.xice = xice; in.snow = snow; in.xlat = xlat;
  in.dz8w = sd_const(dz8w);
  in.p8w = sd_const(p8w); in.t8w = sd_const(t8w); in.p3d = sd_const(p3d); in.pi3d = sd_const(pi3d); in.t3d = sd_const(t3d); in.qv3d = sd_const(qv3d);
  in.qc3d = sd_const(qc3d); in.qi3d = sd_const(qi3d); in.qs3d = sd_const(qs3d);
  in.cldfra3d = sd_const(zero.empty() ? cldfra3d : zero.data());
  in.o33d = sd_const(o33d ? o33d : qv3d);
  in.re_cloud = sd_const(re_cloud); in.re_ice = sd_const(re_ice); in.re_snow = sd_const(re_snow);
  in.has_qc = in.has_qi = in.has_qs = true;  // absent species are passed as zero arrays
  rrtmg_sw_prep(c, in, v, tables, rrtmg_lw_consts());
  for (int g = 0; g < 112; ++g) rrtmg_sw_gpoint(c, g, v, tables);
  for (int lev = 0; lev <= c.nlayers; ++lev) rrtmg_sw_level(lev, v, tables);
  out = rrtmg_sw_finish(c, v, sd_const(pi3d), SD{rthratensw, 1}, SD{rthratenswc, 1});
}

} // namespace cppwrf::phys
