#pragma once
// WRF 4.8.0 phys/module_ra_rrtm.F (RRTM longwave, ra_lw_physics = 1) for one column: RRTMLWRAD's column set-up,
// RRTM (INIRAD / O3DATA, MM5ATM, SETCOEF, GASABS with TAUGB1..16, RTRN) and RRTMLWRAD's tendencies, statement by
// statement. Host code: the g-point arrays of a column (TAUG, PFRAC, ITR, RTRN's BBU, ABSS, ...) are 140 x
// NLAYERS each, kept in a workspace per thread.
#include <string>
#include <vector>

#include "phys/ra_rrtm_lw/rrtm_lw_tables.hpp"

namespace cppwrf::phys {

// rrtminit: NLAYERS = kme + nint(p_top*0.01/deltap) - 1, the model layers plus the 4 hPa buffer layers up to 0 hPa
int rrtm_lw_nlayers(int nz, double p_top);

// the column arrays of RRTMLWRAD, RRTM, MM5ATM and RTRN with Fortran's bounds (index 0 unused for arrays from 1)
struct RrtmLwWorkspace {
  RrtmLwWorkspace(int nz, int nlayers);
  int nz, nlayers;
  // kts:kte(+1)
  std::vector<double> pw1d, tw1d, tten1d, tten1dc, cldfra1d, dz1d, p1d, t1d, qv1d, qc1d, qr1d, qi1d, qs1d, o3prof,
      cldfrc, o3, n2o, ch4, clwp, ciwp, plwp, piwp;
  // 1:NLAYERS
  std::vector<double> pavel, tavel, cldfrac, taucloud, coldry, colh2o, colco2, colo3, coln2o, colch4, colo2,
      co2mult, fac00, fac01, fac10, fac11, forfac, selffac, selffrac, o3prof2, pzr, varint, tlayfrac, abscld, odcld,
      efclfrac;
  std::vector<int> icldlyr, jp, jt, jt1, indself, indlay;
  // 0:NLAYERS
  std::vector<double> pz, tz, totdflux, totuflux, htr, htrc, totuclfl, totdclfl, fnet, fnetc, tlevfrac;
  std::vector<int> indlev;
  // (NGPT, NLAYERS), (35, NLAYERS), (MAXXSEC, NLAYERS), (NBANDS, 0:NLAYERS), (NBANDS, NLAYERS), NGPT*NLAYERS
  std::vector<double> taug, pfrac, wkl, wx, odclr, plvl, play, bbu, bbutot, abss, atot;
  std::vector<int> itr;
};

struct RrtmLwColumn {
  // model levels k = 0 (lowest) .. nz-1, bottom-up as WRF's 3-D arrays; p8w, t8w: nz+1 w levels
  const double *p8w, *t8w, *p, *t, *dz, *pi, *qv;
  // the hydrometeors the radiation sees (RRTMLWRAD reads them with icloud /= 0): a null qc / qr / qs is a species
  // the scheme lacks (F_Q* false); qi is read whenever F_QI is passed (a scheme without ice: moist slot 1, 0)
  const double *qc, *qr, *qi, *qs, *cldfra;
  bool has_qi;  // F_QI: without it (and not warm_rain) cloud and rain below 273.15 K are ice and snow
  double tsk, emiss;
  double co2, n2o, ch4;  // read_CAMgases / the fixed values
  double r_d, g;         // radiation_driver's R=R_d, G=G
  int icloud;
  bool warm_rain;
};

// RRTMLWRAD for one column: glw (GLW), olr (RLWTOA), rth / rthc: RTHRATEN / RTHRATENC (= RTHRATENLWC) as
// 0 + TTEN1D(NK)/pi3D (the driver zeroes them first), model levels bottom-up. message: MM5ATM's wrf_message when
// TBOUND exceeds the table limit (empty otherwise)
void rrtm_lw_column(const RrtmLwTables& t, const RrtmLwColumn& c, RrtmLwWorkspace& ws, double& glw, double& olr,
                    double* rth, double* rthc, std::string& message);

// every column of a domain on the host (OpenMP): 3-D arrays with the column index fastest (idx = k*ncol + col,
// p8w / t8w nz+1 levels), 2-D ones by column; a null qc / qr / qi / qs as in RrtmLwColumn. Outputs rth, rthc
// (k*ncol + col), glw, olr. The CUDA path's radiation columns (radiation.cpp is device code without OpenMP).
struct RrtmLwBatchInput {
  int ncol, nz;
  const double *p8w, *t8w, *p, *t, *dz, *pi, *qv, *qc, *qr, *qi, *qs, *cldfra, *tsk, *emiss;
  bool has_qi;
  double co2, n2o, ch4, r_d, g;
  int icloud;
  bool warm_rain;
};
void rrtm_lw_batch(const RrtmLwTables& t, const RrtmLwBatchInput& in, int nlayers, double* rth, double* rthc,
                   double* glw, double* olr);

}  // namespace cppwrf::phys
