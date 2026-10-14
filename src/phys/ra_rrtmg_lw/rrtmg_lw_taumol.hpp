#pragma once

#include "core/backend.hpp"
#include "phys/ra_rrtmg_common/strided.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_tables.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_setcoef.hpp"
#include <algorithm>
#include <cmath>

namespace cppwrf::phys {

CPPWRF_HD inline void rrtmg_lw_taugb01(
    int nlayers,
    SD pavel,
    SD2 wx,
    SD coldry,
    int laytrop,
    SI jp,
    SI jt,
    SI jt1,
    SD colh2o,
    SD colco2,
    SD colo3,
    SD coln2o,
    SD colco,
    SD colch4,
    SD colo2,
    SD colbrd,
    SD fac00,
    SD fac01,
    SD fac10,
    SD fac11,
    SD rat_h2oco2,
    SD rat_h2oco2_1,
    SD rat_h2oo3,
    SD rat_h2oo3_1,
    SD rat_h2on2o,
    SD rat_h2on2o_1,
    SD rat_h2och4,
    SD rat_h2och4_1,
    SD rat_n2oco2,
    SD rat_n2oco2_1,
    SD rat_o3co2,
    SD rat_o3co2_1,
    SD selffac,
    SD selffrac,
    SI indself,
    SD forfac,
    SD forfrac,
    SI indfor,
    SD minorfrac,
    SD scaleminor,
    SD scaleminorn2,
    SI indminor,
    SD2 fracs,
    SD2 taug,
    const RrtmgLwTables& tables,
    int lay0, int lay1  // layers [lay0, lay1) of the column
) {
  constexpr double oneminus = 1.0 - 1.0e-6;
  const auto& b = tables.b01;
  int lay, ind0, ind1, inds, indf, indm, ig;
  double pp, corradj, scalen2, tauself, taufor, taun2;
  for (int lay = lay0; lay < (laytrop < lay1 ? laytrop : lay1); ++lay) {
  ind0 = ((jp[lay]-1)*5+(jt[lay]-1))*1 + 1;
  ind1 = (jp[lay]*5+(jt1[lay]-1))*1 + 1;
  inds = indself[lay];
  indf = indfor[lay];
  indm = indminor[lay];
  pp = pavel[lay];
  corradj =  1.;
  if (pp  <  250.) {
  corradj = 1. - 0.15 * (250.-pp) / 154.4;
  }
  scalen2 = colbrd[lay] * scaleminorn2[lay];
  for (int ig = 1; ig <= 10; ++ig) {
  tauself = selffac[lay] * (b.selfref[(ig - 1) * 10 + (inds - 1)] + selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds+1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)]));
  taufor =  forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] -  b.forref[(ig - 1) * 4 + (indf - 1)]));
  taun2 = scalen2*(b.ka_mn2[(ig - 1) * 19 + (indm - 1)] + minorfrac[lay] * (b.ka_mn2[(ig - 1) * 19 + (indm+1 - 1)] - b.ka_mn2[(ig - 1) * 19 + (indm - 1)]));
  taug[lay][ig - 1] = corradj * (colh2o[lay] * (fac00[lay] * b.ka[(ig - 1) * 65 + (ind0 - 1)] + fac10[lay] * b.ka[(ig - 1) * 65 + (ind0+1 - 1)] + fac01[lay] * b.ka[(ig - 1) * 65 + (ind1 - 1)] + fac11[lay] * b.ka[(ig - 1) * 65 + (ind1+1 - 1)]) + tauself + taufor + taun2);
  fracs[lay][ig - 1] = b.fracrefa[ig - 1];
  }
  }
  for (int lay = (laytrop > lay0 ? laytrop : lay0); lay < (nlayers < lay1 ? nlayers : lay1); ++lay) {
  ind0 = ((jp[lay]-13)*5+(jt[lay]-1))*1 + 1;
  ind1 = ((jp[lay]-12)*5+(jt1[lay]-1))*1 + 1;
  indf = indfor[lay];
  indm = indminor[lay];
  pp = pavel[lay];
  corradj =  1. - 0.15 * (pp / 95.6);
  scalen2 = colbrd[lay] * scaleminorn2[lay];
  for (int ig = 1; ig <= 10; ++ig) {
  taufor = forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  taun2 = scalen2*(b.kb_mn2[(ig - 1) * 19 + (indm - 1)] + minorfrac[lay] * (b.kb_mn2[(ig - 1) * 19 + (indm+1 - 1)] - b.kb_mn2[(ig - 1) * 19 + (indm - 1)]));
  taug[lay][ig - 1] = corradj * (colh2o[lay] * (fac00[lay] * b.kb[(ig - 1) * 235 + (ind0 - 1)] + fac10[lay] * b.kb[(ig - 1) * 235 + (ind0+1 - 1)] + fac01[lay] * b.kb[(ig - 1) * 235 + (ind1 - 1)] + fac11[lay] * b.kb[(ig - 1) * 235 + (ind1+1 - 1)]) + taufor + taun2);
  fracs[lay][ig - 1] = b.fracrefb[ig - 1];
  }
  }
}

CPPWRF_HD inline void rrtmg_lw_taugb02(
    int nlayers,
    SD pavel,
    SD2 wx,
    SD coldry,
    int laytrop,
    SI jp,
    SI jt,
    SI jt1,
    SD colh2o,
    SD colco2,
    SD colo3,
    SD coln2o,
    SD colco,
    SD colch4,
    SD colo2,
    SD colbrd,
    SD fac00,
    SD fac01,
    SD fac10,
    SD fac11,
    SD rat_h2oco2,
    SD rat_h2oco2_1,
    SD rat_h2oo3,
    SD rat_h2oo3_1,
    SD rat_h2on2o,
    SD rat_h2on2o_1,
    SD rat_h2och4,
    SD rat_h2och4_1,
    SD rat_n2oco2,
    SD rat_n2oco2_1,
    SD rat_o3co2,
    SD rat_o3co2_1,
    SD selffac,
    SD selffrac,
    SI indself,
    SD forfac,
    SD forfrac,
    SI indfor,
    SD minorfrac,
    SD scaleminor,
    SD scaleminorn2,
    SI indminor,
    SD2 fracs,
    SD2 taug,
    const RrtmgLwTables& tables,
    int lay0, int lay1  // layers [lay0, lay1) of the column
) {
  constexpr double oneminus = 1.0 - 1.0e-6;
  const auto& b = tables.b02;
  int lay, ind0, ind1, inds, indf, ig;
  double pp, corradj, tauself, taufor;
  for (int lay = lay0; lay < (laytrop < lay1 ? laytrop : lay1); ++lay) {
  ind0 = ((jp[lay]-1)*5+(jt[lay]-1))*1 + 1;
  ind1 = (jp[lay]*5+(jt1[lay]-1))*1 + 1;
  inds = indself[lay];
  indf = indfor[lay];
  pp = pavel[lay];
  corradj = 1. - .05 * (pp - 100.) / 900.;
  for (int ig = 1; ig <= 12; ++ig) {
  tauself = selffac[lay] * (b.selfref[(ig - 1) * 10 + (inds - 1)] + selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds+1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)]));
  taufor =  forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  taug[lay][10+ig - 1] = corradj * (colh2o[lay] * (fac00[lay] * b.ka[(ig - 1) * 65 + (ind0 - 1)] + fac10[lay] * b.ka[(ig - 1) * 65 + (ind0+1 - 1)] + fac01[lay] * b.ka[(ig - 1) * 65 + (ind1 - 1)] + fac11[lay] * b.ka[(ig - 1) * 65 + (ind1+1 - 1)]) + tauself + taufor);
  fracs[lay][10+ig - 1] = b.fracrefa[ig - 1];
  }
  }
  for (int lay = (laytrop > lay0 ? laytrop : lay0); lay < (nlayers < lay1 ? nlayers : lay1); ++lay) {
  ind0 = ((jp[lay]-13)*5+(jt[lay]-1))*1 + 1;
  ind1 = ((jp[lay]-12)*5+(jt1[lay]-1))*1 + 1;
  indf = indfor[lay];
  for (int ig = 1; ig <= 12; ++ig) {
  taufor =  forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  taug[lay][10+ig - 1] = colh2o[lay] * (fac00[lay] * b.kb[(ig - 1) * 235 + (ind0 - 1)] + fac10[lay] * b.kb[(ig - 1) * 235 + (ind0+1 - 1)] + fac01[lay] * b.kb[(ig - 1) * 235 + (ind1 - 1)] + fac11[lay] * b.kb[(ig - 1) * 235 + (ind1+1 - 1)]) + taufor;
  fracs[lay][10+ig - 1] = b.fracrefb[ig - 1];
  }
  }
}

CPPWRF_HD inline void rrtmg_lw_taugb03(
    int nlayers,
    SD pavel,
    SD2 wx,
    SD coldry,
    int laytrop,
    SI jp,
    SI jt,
    SI jt1,
    SD colh2o,
    SD colco2,
    SD colo3,
    SD coln2o,
    SD colco,
    SD colch4,
    SD colo2,
    SD colbrd,
    SD fac00,
    SD fac01,
    SD fac10,
    SD fac11,
    SD rat_h2oco2,
    SD rat_h2oco2_1,
    SD rat_h2oo3,
    SD rat_h2oo3_1,
    SD rat_h2on2o,
    SD rat_h2on2o_1,
    SD rat_h2och4,
    SD rat_h2och4_1,
    SD rat_n2oco2,
    SD rat_n2oco2_1,
    SD rat_o3co2,
    SD rat_o3co2_1,
    SD selffac,
    SD selffrac,
    SI indself,
    SD forfac,
    SD forfrac,
    SI indfor,
    SD minorfrac,
    SD scaleminor,
    SD scaleminorn2,
    SI indminor,
    SD2 fracs,
    SD2 taug,
    const RrtmgLwTables& tables,
    int lay0, int lay1  // layers [lay0, lay1) of the column
) {
  constexpr double oneminus = 1.0 - 1.0e-6;
  const auto& b = tables.b03;
  int lay, ind0, ind1, inds, indf, indm, ig;
  int js, js1, jmn2o, jpl;
  double speccomb, specparm, specmult, fs;
  double speccomb1, specparm1, specmult1, fs1;
  double speccomb_mn2o, specparm_mn2o, specmult_mn2o, fmn2o, fmn2omf, chi_n2o, ratn2o, adjfac, adjcoln2o;
  double speccomb_planck, specparm_planck, specmult_planck, fpl;
  double p, p4, fk0, fk1, fk2;
  double fac000, fac100, fac200, fac010, fac110, fac210;
  double fac001, fac101, fac201, fac011, fac111, fac211;
  double tauself, taufor, n2om1, n2om2, absn2o;
  double refrat_planck_a, refrat_planck_b, refrat_m_a, refrat_m_b;
  double tau_major, tau_major1;
  refrat_planck_a = get_chi_mls(tables, 1, 9)/get_chi_mls(tables, 2, 9);
  refrat_planck_b = get_chi_mls(tables, 1, 13)/get_chi_mls(tables, 2, 13);
  refrat_m_a = get_chi_mls(tables, 1, 3)/get_chi_mls(tables, 2, 3);
  refrat_m_b = get_chi_mls(tables, 1, 13)/get_chi_mls(tables, 2, 13);
  for (int lay = lay0; lay < (laytrop < lay1 ? laytrop : lay1); ++lay) {
  speccomb = colh2o[lay] + rat_h2oco2[lay]*colco2[lay];
  specparm = colh2o[lay]/speccomb;
  if (specparm  >=  oneminus) specparm = oneminus;
  specmult = 8.*(specparm);
  js = 1 + static_cast<int>(specmult);
  fs = std::fmod(specmult, 1.0);
  speccomb1 = colh2o[lay] + rat_h2oco2_1[lay]*colco2[lay];
  specparm1 = colh2o[lay]/speccomb1;
  if (specparm1  >=  oneminus) specparm1 = oneminus;
  specmult1 = 8.*(specparm1);
  js1 = 1 + static_cast<int>(specmult1);
  fs1 = std::fmod(specmult1, 1.0);
  speccomb_mn2o = colh2o[lay] + refrat_m_a*colco2[lay];
  specparm_mn2o = colh2o[lay]/speccomb_mn2o;
  if (specparm_mn2o  >=  oneminus) specparm_mn2o = oneminus;
  specmult_mn2o = 8.*specparm_mn2o;
  jmn2o = 1 + static_cast<int>(specmult_mn2o);
  fmn2o = std::fmod(specmult_mn2o, 1.0);
  fmn2omf = minorfrac[lay]*fmn2o;
  chi_n2o = coln2o[lay]/coldry[lay];
  ratn2o = 1.e20*chi_n2o/get_chi_mls(tables, 4, jp[lay]+1);
  if (ratn2o  >  1.5) {
  adjfac = 0.5+std::pow((ratn2o-0.5), 0.65);
  adjcoln2o = adjfac*get_chi_mls(tables, 4, jp[lay]+1)*coldry[lay]*1.e-20;
  } else {
  adjcoln2o = coln2o[lay];
  }
  speccomb_planck = colh2o[lay]+refrat_planck_a*colco2[lay];
  specparm_planck = colh2o[lay]/speccomb_planck;
  if (specparm_planck  >=  oneminus) specparm_planck=oneminus;
  specmult_planck = 8.*specparm_planck;
  jpl= 1 + static_cast<int>(specmult_planck);
  fpl = std::fmod(specmult_planck, 1.0);
  ind0 = ((jp[lay]-1)*5+(jt[lay]-1))*9 + js;
  ind1 = (jp[lay]*5+(jt1[lay]-1))*9 + js1;
  inds = indself[lay];
  indf = indfor[lay];
  indm = indminor[lay];
  if (specparm  <  0.125) {
  p = fs - 1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac000 = fk0*fac00[lay];
  fac100 = fk1*fac00[lay];
  fac200 = fk2*fac00[lay];
  fac010 = fk0*fac10[lay];
  fac110 = fk1*fac10[lay];
  fac210 = fk2*fac10[lay];
  } else if (specparm  >  0.875) {
  p = -fs;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac000 = fk0*fac00[lay];
  fac100 = fk1*fac00[lay];
  fac200 = fk2*fac00[lay];
  fac010 = fk0*fac10[lay];
  fac110 = fk1*fac10[lay];
  fac210 = fk2*fac10[lay];
  } else {
  fac000 = (1. - fs) * fac00[lay];
  fac010 = (1. - fs) * fac10[lay];
  fac100 = fs * fac00[lay];
  fac110 = fs * fac10[lay];
  }
  if (specparm1  <  0.125) {
  p = fs1 - 1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac001 = fk0*fac01[lay];
  fac101 = fk1*fac01[lay];
  fac201 = fk2*fac01[lay];
  fac011 = fk0*fac11[lay];
  fac111 = fk1*fac11[lay];
  fac211 = fk2*fac11[lay];
  } else if (specparm1  >  0.875) {
  p = -fs1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac001 = fk0*fac01[lay];
  fac101 = fk1*fac01[lay];
  fac201 = fk2*fac01[lay];
  fac011 = fk0*fac11[lay];
  fac111 = fk1*fac11[lay];
  fac211 = fk2*fac11[lay];
  } else {
  fac001 = (1. - fs1) * fac01[lay];
  fac011 = (1. - fs1) * fac11[lay];
  fac101 = fs1 * fac01[lay];
  fac111 = fs1 * fac11[lay];
  }
  for (int ig = 1; ig <= 16; ++ig) {
  tauself = selffac[lay]* (b.selfref[(ig - 1) * 10 + (inds - 1)] + selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds+1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)]));
  taufor = forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  n2om1 = b.ka_mn2o[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmn2o - 1)] + fmn2o * (b.ka_mn2o[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmn2o+1 - 1)] - b.ka_mn2o[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmn2o - 1)]);
  n2om2 = b.ka_mn2o[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmn2o - 1)] + fmn2o * (b.ka_mn2o[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmn2o+1 - 1)] - b.ka_mn2o[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmn2o - 1)]);
  absn2o = n2om1 + minorfrac[lay] * (n2om2 - n2om1);
  if (specparm  <  0.125) {
  tau_major = speccomb * (fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac200 * b.ka[(ig - 1) * 585 + (ind0+2 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)] + fac210 * b.ka[(ig - 1) * 585 + (ind0+11 - 1)]);
  } else if (specparm  >  0.875) {
  tau_major = speccomb * (fac200 * b.ka[(ig - 1) * 585 + (ind0-1 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac000 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac210 * b.ka[(ig - 1) * 585 + (ind0+8 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)]);
  } else {
  tau_major = speccomb * (fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)]);
  }
  if (specparm1  <  0.125) {
  tau_major1 = speccomb1 * (fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac201 * b.ka[(ig - 1) * 585 + (ind1+2 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)] + fac211 * b.ka[(ig - 1) * 585 + (ind1+11 - 1)]);
  } else if (specparm1  >  0.875) {
  tau_major1 = speccomb1 * (fac201 * b.ka[(ig - 1) * 585 + (ind1-1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac001 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac211 * b.ka[(ig - 1) * 585 + (ind1+8 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)]);
  } else {
  tau_major1 = speccomb1 * (fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)]);
  }
  taug[lay][22+ig - 1] = tau_major + tau_major1 + tauself + taufor + adjcoln2o*absn2o;
  fracs[lay][22+ig - 1] = b.fracrefa[(jpl - 1) * 16 + (ig - 1)] + fpl * (b.fracrefa[(jpl+1 - 1) * 16 + (ig - 1)]-b.fracrefa[(jpl - 1) * 16 + (ig - 1)]);
  }
  }
  for (int lay = (laytrop > lay0 ? laytrop : lay0); lay < (nlayers < lay1 ? nlayers : lay1); ++lay) {
  speccomb = colh2o[lay] + rat_h2oco2[lay]*colco2[lay];
  specparm = colh2o[lay]/speccomb;
  if (specparm  >=  oneminus) specparm = oneminus;
  specmult = 4.*(specparm);
  js = 1 + static_cast<int>(specmult);
  fs = std::fmod(specmult, 1.0);
  speccomb1 = colh2o[lay] + rat_h2oco2_1[lay]*colco2[lay];
  specparm1 = colh2o[lay]/speccomb1;
  if (specparm1  >=  oneminus) specparm1 = oneminus;
  specmult1 = 4.*(specparm1);
  js1 = 1 + static_cast<int>(specmult1);
  fs1 = std::fmod(specmult1, 1.0);
  fac000 = (1. - fs) * fac00[lay];
  fac010 = (1. - fs) * fac10[lay];
  fac100 = fs * fac00[lay];
  fac110 = fs * fac10[lay];
  fac001 = (1. - fs1) * fac01[lay];
  fac011 = (1. - fs1) * fac11[lay];
  fac101 = fs1 * fac01[lay];
  fac111 = fs1 * fac11[lay];
  speccomb_mn2o = colh2o[lay] + refrat_m_b*colco2[lay];
  specparm_mn2o = colh2o[lay]/speccomb_mn2o;
  if (specparm_mn2o  >=  oneminus) specparm_mn2o = oneminus;
  specmult_mn2o = 4.*specparm_mn2o;
  jmn2o = 1 + static_cast<int>(specmult_mn2o);
  fmn2o = std::fmod(specmult_mn2o, 1.0);
  fmn2omf = minorfrac[lay]*fmn2o;
  chi_n2o = coln2o[lay]/coldry[lay];
  ratn2o = 1.e20*chi_n2o/get_chi_mls(tables, 4, jp[lay]+1);
  if (ratn2o  >  1.5) {
  adjfac = 0.5+std::pow((ratn2o-0.5), 0.65);
  adjcoln2o = adjfac*get_chi_mls(tables, 4, jp[lay]+1)*coldry[lay]*1.e-20;
  } else {
  adjcoln2o = coln2o[lay];
  }
  speccomb_planck = colh2o[lay]+refrat_planck_b*colco2[lay];
  specparm_planck = colh2o[lay]/speccomb_planck;
  if (specparm_planck  >=  oneminus) specparm_planck=oneminus;
  specmult_planck = 4.*specparm_planck;
  jpl= 1 + static_cast<int>(specmult_planck);
  fpl = std::fmod(specmult_planck, 1.0);
  ind0 = ((jp[lay]-13)*5+(jt[lay]-1))*5 + js;
  ind1 = ((jp[lay]-12)*5+(jt1[lay]-1))*5 + js1;
  indf = indfor[lay];
  indm = indminor[lay];
  for (int ig = 1; ig <= 16; ++ig) {
  taufor = forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  n2om1 = b.kb_mn2o[(ig - 1) * (19 * 5) + (indm - 1) * 5 + (jmn2o - 1)] + fmn2o * (b.kb_mn2o[(ig - 1) * (19 * 5) + (indm - 1) * 5 + (jmn2o+1 - 1)]-b.kb_mn2o[(ig - 1) * (19 * 5) + (indm - 1) * 5 + (jmn2o - 1)]);
  n2om2 = b.kb_mn2o[(ig - 1) * (19 * 5) + (indm+1 - 1) * 5 + (jmn2o - 1)] + fmn2o * (b.kb_mn2o[(ig - 1) * (19 * 5) + (indm+1 - 1) * 5 + (jmn2o+1 - 1)]-b.kb_mn2o[(ig - 1) * (19 * 5) + (indm+1 - 1) * 5 + (jmn2o - 1)]);
  absn2o = n2om1 + minorfrac[lay] * (n2om2 - n2om1);
  taug[lay][22+ig - 1] = speccomb * (fac000 * b.kb[(ig - 1) * 1175 + (ind0 - 1)] + fac100 * b.kb[(ig - 1) * 1175 + (ind0+1 - 1)] + fac010 * b.kb[(ig - 1) * 1175 + (ind0+5 - 1)] + fac110 * b.kb[(ig - 1) * 1175 + (ind0+6 - 1)]) + speccomb1 * (fac001 * b.kb[(ig - 1) * 1175 + (ind1 - 1)] + fac101 * b.kb[(ig - 1) * 1175 + (ind1+1 - 1)] + fac011 * b.kb[(ig - 1) * 1175 + (ind1+5 - 1)] + fac111 * b.kb[(ig - 1) * 1175 + (ind1+6 - 1)]) + taufor + adjcoln2o*absn2o;
  fracs[lay][22+ig - 1] = b.fracrefb[(jpl - 1) * 16 + (ig - 1)] + fpl * (b.fracrefb[(jpl+1 - 1) * 16 + (ig - 1)]-b.fracrefb[(jpl - 1) * 16 + (ig - 1)]);
  }
  }
}

CPPWRF_HD inline void rrtmg_lw_taugb04(
    int nlayers,
    SD pavel,
    SD2 wx,
    SD coldry,
    int laytrop,
    SI jp,
    SI jt,
    SI jt1,
    SD colh2o,
    SD colco2,
    SD colo3,
    SD coln2o,
    SD colco,
    SD colch4,
    SD colo2,
    SD colbrd,
    SD fac00,
    SD fac01,
    SD fac10,
    SD fac11,
    SD rat_h2oco2,
    SD rat_h2oco2_1,
    SD rat_h2oo3,
    SD rat_h2oo3_1,
    SD rat_h2on2o,
    SD rat_h2on2o_1,
    SD rat_h2och4,
    SD rat_h2och4_1,
    SD rat_n2oco2,
    SD rat_n2oco2_1,
    SD rat_o3co2,
    SD rat_o3co2_1,
    SD selffac,
    SD selffrac,
    SI indself,
    SD forfac,
    SD forfrac,
    SI indfor,
    SD minorfrac,
    SD scaleminor,
    SD scaleminorn2,
    SI indminor,
    SD2 fracs,
    SD2 taug,
    const RrtmgLwTables& tables,
    int lay0, int lay1  // layers [lay0, lay1) of the column
) {
  constexpr double oneminus = 1.0 - 1.0e-6;
  const auto& b = tables.b04;
  int lay, ind0, ind1, inds, indf, ig;
  int js, js1, jpl;
  double speccomb, specparm, specmult, fs;
  double speccomb1, specparm1, specmult1, fs1;
  double speccomb_planck, specparm_planck, specmult_planck, fpl;
  double p, p4, fk0, fk1, fk2;
  double fac000, fac100, fac200, fac010, fac110, fac210;
  double fac001, fac101, fac201, fac011, fac111, fac211;
  double tauself, taufor;
  double refrat_planck_a, refrat_planck_b;
  double tau_major, tau_major1;
  refrat_planck_a = get_chi_mls(tables, 1, 11)/get_chi_mls(tables, 2, 11);
  refrat_planck_b = get_chi_mls(tables, 3, 13)/get_chi_mls(tables, 2, 13);
  for (int lay = lay0; lay < (laytrop < lay1 ? laytrop : lay1); ++lay) {
  speccomb = colh2o[lay] + rat_h2oco2[lay]*colco2[lay];
  specparm = colh2o[lay]/speccomb;
  if (specparm  >=  oneminus) specparm = oneminus;
  specmult = 8.*(specparm);
  js = 1 + static_cast<int>(specmult);
  fs = std::fmod(specmult, 1.0);
  speccomb1 = colh2o[lay] + rat_h2oco2_1[lay]*colco2[lay];
  specparm1 = colh2o[lay]/speccomb1;
  if (specparm1  >=  oneminus) specparm1 = oneminus;
  specmult1 = 8.*(specparm1);
  js1 = 1 + static_cast<int>(specmult1);
  fs1 = std::fmod(specmult1, 1.0);
  speccomb_planck = colh2o[lay]+refrat_planck_a*colco2[lay];
  specparm_planck = colh2o[lay]/speccomb_planck;
  if (specparm_planck  >=  oneminus) specparm_planck=oneminus;
  specmult_planck = 8.*specparm_planck;
  jpl= 1 + static_cast<int>(specmult_planck);
  fpl = std::fmod(specmult_planck, 1.0);
  ind0 = ((jp[lay]-1)*5+(jt[lay]-1))*9 + js;
  ind1 = (jp[lay]*5+(jt1[lay]-1))*9 + js1;
  inds = indself[lay];
  indf = indfor[lay];
  if (specparm  <  0.125) {
  p = fs - 1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac000 = fk0*fac00[lay];
  fac100 = fk1*fac00[lay];
  fac200 = fk2*fac00[lay];
  fac010 = fk0*fac10[lay];
  fac110 = fk1*fac10[lay];
  fac210 = fk2*fac10[lay];
  } else if (specparm  >  0.875) {
  p = -fs;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac000 = fk0*fac00[lay];
  fac100 = fk1*fac00[lay];
  fac200 = fk2*fac00[lay];
  fac010 = fk0*fac10[lay];
  fac110 = fk1*fac10[lay];
  fac210 = fk2*fac10[lay];
  } else {
  fac000 = (1. - fs) * fac00[lay];
  fac010 = (1. - fs) * fac10[lay];
  fac100 = fs * fac00[lay];
  fac110 = fs * fac10[lay];
  }
  if (specparm1  <  0.125) {
  p = fs1 - 1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac001 = fk0*fac01[lay];
  fac101 = fk1*fac01[lay];
  fac201 = fk2*fac01[lay];
  fac011 = fk0*fac11[lay];
  fac111 = fk1*fac11[lay];
  fac211 = fk2*fac11[lay];
  } else if (specparm1  >  0.875) {
  p = -fs1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac001 = fk0*fac01[lay];
  fac101 = fk1*fac01[lay];
  fac201 = fk2*fac01[lay];
  fac011 = fk0*fac11[lay];
  fac111 = fk1*fac11[lay];
  fac211 = fk2*fac11[lay];
  } else {
  fac001 = (1. - fs1) * fac01[lay];
  fac011 = (1. - fs1) * fac11[lay];
  fac101 = fs1 * fac01[lay];
  fac111 = fs1 * fac11[lay];
  }
  for (int ig = 1; ig <= 14; ++ig) {
  tauself = selffac[lay]* (b.selfref[(ig - 1) * 10 + (inds - 1)] + selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds+1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)]));
  taufor =  forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  if (specparm  <  0.125) {
  tau_major = speccomb * (fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac200 * b.ka[(ig - 1) * 585 + (ind0+2 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)] + fac210 * b.ka[(ig - 1) * 585 + (ind0+11 - 1)]);
  } else if (specparm  >  0.875) {
  tau_major = speccomb * (fac200 * b.ka[(ig - 1) * 585 + (ind0-1 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac000 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac210 * b.ka[(ig - 1) * 585 + (ind0+8 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)]);
  } else {
  tau_major = speccomb * (fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)]);
  }
  if (specparm1  <  0.125) {
  tau_major1 = speccomb1 * (fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac201 * b.ka[(ig - 1) * 585 + (ind1+2 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)] + fac211 * b.ka[(ig - 1) * 585 + (ind1+11 - 1)]);
  } else if (specparm1  >  0.875) {
  tau_major1 = speccomb1 * (fac201 * b.ka[(ig - 1) * 585 + (ind1-1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac001 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac211 * b.ka[(ig - 1) * 585 + (ind1+8 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)]);
  } else {
  tau_major1 = speccomb1 * (fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)]);
  }
  taug[lay][38+ig - 1] = tau_major + tau_major1 + tauself + taufor;
  fracs[lay][38+ig - 1] = b.fracrefa[(jpl - 1) * 14 + (ig - 1)] + fpl * (b.fracrefa[(jpl+1 - 1) * 14 + (ig - 1)]-b.fracrefa[(jpl - 1) * 14 + (ig - 1)]);
  }
  }
  for (int lay = (laytrop > lay0 ? laytrop : lay0); lay < (nlayers < lay1 ? nlayers : lay1); ++lay) {
  speccomb = colo3[lay] + rat_o3co2[lay]*colco2[lay];
  specparm = colo3[lay]/speccomb;
  if (specparm  >=  oneminus) specparm = oneminus;
  specmult = 4.*(specparm);
  js = 1 + static_cast<int>(specmult);
  fs = std::fmod(specmult, 1.0);
  speccomb1 = colo3[lay] + rat_o3co2_1[lay]*colco2[lay];
  specparm1 = colo3[lay]/speccomb1;
  if (specparm1  >=  oneminus) specparm1 = oneminus;
  specmult1 = 4.*(specparm1);
  js1 = 1 + static_cast<int>(specmult1);
  fs1 = std::fmod(specmult1, 1.0);
  fac000 = (1. - fs) * fac00[lay];
  fac010 = (1. - fs) * fac10[lay];
  fac100 = fs * fac00[lay];
  fac110 = fs * fac10[lay];
  fac001 = (1. - fs1) * fac01[lay];
  fac011 = (1. - fs1) * fac11[lay];
  fac101 = fs1 * fac01[lay];
  fac111 = fs1 * fac11[lay];
  speccomb_planck = colo3[lay]+refrat_planck_b*colco2[lay];
  specparm_planck = colo3[lay]/speccomb_planck;
  if (specparm_planck  >=  oneminus) specparm_planck=oneminus;
  specmult_planck = 4.*specparm_planck;
  jpl= 1 + static_cast<int>(specmult_planck);
  fpl = std::fmod(specmult_planck, 1.0);
  ind0 = ((jp[lay]-13)*5+(jt[lay]-1))*5 + js;
  ind1 = ((jp[lay]-12)*5+(jt1[lay]-1))*5 + js1;
  for (int ig = 1; ig <= 14; ++ig) {
  taug[lay][38+ig - 1] =  speccomb * (fac000 * b.kb[(ig - 1) * 1175 + (ind0 - 1)] + fac100 * b.kb[(ig - 1) * 1175 + (ind0+1 - 1)] + fac010 * b.kb[(ig - 1) * 1175 + (ind0+5 - 1)] + fac110 * b.kb[(ig - 1) * 1175 + (ind0+6 - 1)]) + speccomb1 * (fac001 * b.kb[(ig - 1) * 1175 + (ind1 - 1)] + fac101 * b.kb[(ig - 1) * 1175 + (ind1+1 - 1)] + fac011 * b.kb[(ig - 1) * 1175 + (ind1+5 - 1)] + fac111 * b.kb[(ig - 1) * 1175 + (ind1+6 - 1)]);
  fracs[lay][38+ig - 1] = b.fracrefb[(jpl - 1) * 14 + (ig - 1)] + fpl * (b.fracrefb[(jpl+1 - 1) * 14 + (ig - 1)]-b.fracrefb[(jpl - 1) * 14 + (ig - 1)]);
  }
  taug[lay][38+8 - 1]=taug[lay][38+8 - 1]*0.92;
  taug[lay][38+9 - 1]=taug[lay][38+9 - 1]*0.88;
  taug[lay][38+10 - 1]=taug[lay][38+10 - 1]*1.07;
  taug[lay][38+11 - 1]=taug[lay][38+11 - 1]*1.1;
  taug[lay][38+12 - 1]=taug[lay][38+12 - 1]*0.99;
  taug[lay][38+13 - 1]=taug[lay][38+13 - 1]*0.88;
  taug[lay][38+14 - 1]=taug[lay][38+14 - 1]*0.943;
  }
}

CPPWRF_HD inline void rrtmg_lw_taugb05(
    int nlayers,
    SD pavel,
    SD2 wx,
    SD coldry,
    int laytrop,
    SI jp,
    SI jt,
    SI jt1,
    SD colh2o,
    SD colco2,
    SD colo3,
    SD coln2o,
    SD colco,
    SD colch4,
    SD colo2,
    SD colbrd,
    SD fac00,
    SD fac01,
    SD fac10,
    SD fac11,
    SD rat_h2oco2,
    SD rat_h2oco2_1,
    SD rat_h2oo3,
    SD rat_h2oo3_1,
    SD rat_h2on2o,
    SD rat_h2on2o_1,
    SD rat_h2och4,
    SD rat_h2och4_1,
    SD rat_n2oco2,
    SD rat_n2oco2_1,
    SD rat_o3co2,
    SD rat_o3co2_1,
    SD selffac,
    SD selffrac,
    SI indself,
    SD forfac,
    SD forfrac,
    SI indfor,
    SD minorfrac,
    SD scaleminor,
    SD scaleminorn2,
    SI indminor,
    SD2 fracs,
    SD2 taug,
    const RrtmgLwTables& tables,
    int lay0, int lay1  // layers [lay0, lay1) of the column
) {
  constexpr double oneminus = 1.0 - 1.0e-6;
  const auto& b = tables.b05;
  int lay, ind0, ind1, inds, indf, indm, ig;
  int js, js1, jmo3, jpl;
  double speccomb, specparm, specmult, fs;
  double speccomb1, specparm1, specmult1, fs1;
  double speccomb_mo3, specparm_mo3, specmult_mo3, fmo3;
  double speccomb_planck, specparm_planck, specmult_planck, fpl;
  double p, p4, fk0, fk1, fk2;
  double fac000, fac100, fac200, fac010, fac110, fac210;
  double fac001, fac101, fac201, fac011, fac111, fac211;
  double tauself, taufor, o3m1, o3m2, abso3;
  double refrat_planck_a, refrat_planck_b, refrat_m_a;
  double tau_major, tau_major1;
  refrat_planck_a = get_chi_mls(tables, 1, 5)/get_chi_mls(tables, 2, 5);
  refrat_planck_b = get_chi_mls(tables, 3, 43)/get_chi_mls(tables, 2, 43);
  refrat_m_a = get_chi_mls(tables, 1, 7)/get_chi_mls(tables, 2, 7);
  for (int lay = lay0; lay < (laytrop < lay1 ? laytrop : lay1); ++lay) {
  speccomb = colh2o[lay] + rat_h2oco2[lay]*colco2[lay];
  specparm = colh2o[lay]/speccomb;
  if (specparm  >=  oneminus) specparm = oneminus;
  specmult = 8.*(specparm);
  js = 1 + static_cast<int>(specmult);
  fs = std::fmod(specmult, 1.0);
  speccomb1 = colh2o[lay] + rat_h2oco2_1[lay]*colco2[lay];
  specparm1 = colh2o[lay]/speccomb1;
  if (specparm1  >=  oneminus) specparm1 = oneminus;
  specmult1 = 8.*(specparm1);
  js1 = 1 + static_cast<int>(specmult1);
  fs1 = std::fmod(specmult1, 1.0);
  speccomb_mo3 = colh2o[lay] + refrat_m_a*colco2[lay];
  specparm_mo3 = colh2o[lay]/speccomb_mo3;
  if (specparm_mo3  >=  oneminus) specparm_mo3 = oneminus;
  specmult_mo3 = 8.*specparm_mo3;
  jmo3 = 1 + static_cast<int>(specmult_mo3);
  fmo3 = std::fmod(specmult_mo3, 1.0);
  speccomb_planck = colh2o[lay]+refrat_planck_a*colco2[lay];
  specparm_planck = colh2o[lay]/speccomb_planck;
  if (specparm_planck  >=  oneminus) specparm_planck=oneminus;
  specmult_planck = 8.*specparm_planck;
  jpl= 1 + static_cast<int>(specmult_planck);
  fpl = std::fmod(specmult_planck, 1.0);
  ind0 = ((jp[lay]-1)*5+(jt[lay]-1))*9 + js;
  ind1 = (jp[lay]*5+(jt1[lay]-1))*9 + js1;
  inds = indself[lay];
  indf = indfor[lay];
  indm = indminor[lay];
  if (specparm  <  0.125) {
  p = fs - 1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac000 = fk0*fac00[lay];
  fac100 = fk1*fac00[lay];
  fac200 = fk2*fac00[lay];
  fac010 = fk0*fac10[lay];
  fac110 = fk1*fac10[lay];
  fac210 = fk2*fac10[lay];
  } else if (specparm  >  0.875) {
  p = -fs;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac000 = fk0*fac00[lay];
  fac100 = fk1*fac00[lay];
  fac200 = fk2*fac00[lay];
  fac010 = fk0*fac10[lay];
  fac110 = fk1*fac10[lay];
  fac210 = fk2*fac10[lay];
  } else {
  fac000 = (1. - fs) * fac00[lay];
  fac010 = (1. - fs) * fac10[lay];
  fac100 = fs * fac00[lay];
  fac110 = fs * fac10[lay];
  }
  if (specparm1  <  0.125) {
  p = fs1 - 1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac001 = fk0*fac01[lay];
  fac101 = fk1*fac01[lay];
  fac201 = fk2*fac01[lay];
  fac011 = fk0*fac11[lay];
  fac111 = fk1*fac11[lay];
  fac211 = fk2*fac11[lay];
  } else if (specparm1  >  0.875) {
  p = -fs1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac001 = fk0*fac01[lay];
  fac101 = fk1*fac01[lay];
  fac201 = fk2*fac01[lay];
  fac011 = fk0*fac11[lay];
  fac111 = fk1*fac11[lay];
  fac211 = fk2*fac11[lay];
  } else {
  fac001 = (1. - fs1) * fac01[lay];
  fac011 = (1. - fs1) * fac11[lay];
  fac101 = fs1 * fac01[lay];
  fac111 = fs1 * fac11[lay];
  }
  for (int ig = 1; ig <= 16; ++ig) {
  tauself = selffac[lay] * (b.selfref[(ig - 1) * 10 + (inds - 1)] + selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds+1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)]));
  taufor =  forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  o3m1 = b.ka_mo3[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmo3 - 1)] + fmo3 * (b.ka_mo3[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmo3+1 - 1)]-b.ka_mo3[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmo3 - 1)]);
  o3m2 = b.ka_mo3[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmo3 - 1)] + fmo3 * (b.ka_mo3[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmo3+1 - 1)]-b.ka_mo3[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmo3 - 1)]);
  abso3 = o3m1 + minorfrac[lay]*(o3m2-o3m1);
  if (specparm  <  0.125) {
  tau_major = speccomb * (fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac200 * b.ka[(ig - 1) * 585 + (ind0+2 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)] + fac210 * b.ka[(ig - 1) * 585 + (ind0+11 - 1)]);
  } else if (specparm  >  0.875) {
  tau_major = speccomb * (fac200 * b.ka[(ig - 1) * 585 + (ind0-1 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac000 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac210 * b.ka[(ig - 1) * 585 + (ind0+8 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)]);
  } else {
  tau_major = speccomb * (fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)]);
  }
  if (specparm1  <  0.125) {
  tau_major1 = speccomb1 * (fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac201 * b.ka[(ig - 1) * 585 + (ind1+2 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)] + fac211 * b.ka[(ig - 1) * 585 + (ind1+11 - 1)]);
  } else if (specparm1  >  0.875) {
  tau_major1 = speccomb1 * (fac201 * b.ka[(ig - 1) * 585 + (ind1-1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac001 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac211 * b.ka[(ig - 1) * 585 + (ind1+8 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)]);
  } else {
  tau_major1 = speccomb1 * (fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)]);
  }
  taug[lay][52+ig - 1] = tau_major + tau_major1 + tauself + taufor + abso3*colo3[lay] + wx[1 - 1][lay] * b.ccl4[ig - 1];
  fracs[lay][52+ig - 1] = b.fracrefa[(jpl - 1) * 16 + (ig - 1)] + fpl * (b.fracrefa[(jpl+1 - 1) * 16 + (ig - 1)]-b.fracrefa[(jpl - 1) * 16 + (ig - 1)]);
  }
  }
  for (int lay = (laytrop > lay0 ? laytrop : lay0); lay < (nlayers < lay1 ? nlayers : lay1); ++lay) {
  speccomb = colo3[lay] + rat_o3co2[lay]*colco2[lay];
  specparm = colo3[lay]/speccomb;
  if (specparm  >=  oneminus) specparm = oneminus;
  specmult = 4.*(specparm);
  js = 1 + static_cast<int>(specmult);
  fs = std::fmod(specmult, 1.0);
  speccomb1 = colo3[lay] + rat_o3co2_1[lay]*colco2[lay];
  specparm1 = colo3[lay]/speccomb1;
  if (specparm1  >=  oneminus) specparm1 = oneminus;
  specmult1 = 4.*(specparm1);
  js1 = 1 + static_cast<int>(specmult1);
  fs1 = std::fmod(specmult1, 1.0);
  fac000 = (1. - fs) * fac00[lay];
  fac010 = (1. - fs) * fac10[lay];
  fac100 = fs * fac00[lay];
  fac110 = fs * fac10[lay];
  fac001 = (1. - fs1) * fac01[lay];
  fac011 = (1. - fs1) * fac11[lay];
  fac101 = fs1 * fac01[lay];
  fac111 = fs1 * fac11[lay];
  speccomb_planck = colo3[lay]+refrat_planck_b*colco2[lay];
  specparm_planck = colo3[lay]/speccomb_planck;
  if (specparm_planck  >=  oneminus) specparm_planck=oneminus;
  specmult_planck = 4.*specparm_planck;
  jpl= 1 + static_cast<int>(specmult_planck);
  fpl = std::fmod(specmult_planck, 1.0);
  ind0 = ((jp[lay]-13)*5+(jt[lay]-1))*5 + js;
  ind1 = ((jp[lay]-12)*5+(jt1[lay]-1))*5 + js1;
  for (int ig = 1; ig <= 16; ++ig) {
  taug[lay][52+ig - 1] = speccomb * (fac000 * b.kb[(ig - 1) * 1175 + (ind0 - 1)] + fac100 * b.kb[(ig - 1) * 1175 + (ind0+1 - 1)] + fac010 * b.kb[(ig - 1) * 1175 + (ind0+5 - 1)] + fac110 * b.kb[(ig - 1) * 1175 + (ind0+6 - 1)]) + speccomb1 * (fac001 * b.kb[(ig - 1) * 1175 + (ind1 - 1)] + fac101 * b.kb[(ig - 1) * 1175 + (ind1+1 - 1)] + fac011 * b.kb[(ig - 1) * 1175 + (ind1+5 - 1)] + fac111 * b.kb[(ig - 1) * 1175 + (ind1+6 - 1)]) + wx[1 - 1][lay] * b.ccl4[ig - 1];
  fracs[lay][52+ig - 1] = b.fracrefb[(jpl - 1) * 16 + (ig - 1)] + fpl * (b.fracrefb[(jpl+1 - 1) * 16 + (ig - 1)]-b.fracrefb[(jpl - 1) * 16 + (ig - 1)]);
  }
  }
}

CPPWRF_HD inline void rrtmg_lw_taugb06(
    int nlayers,
    SD pavel,
    SD2 wx,
    SD coldry,
    int laytrop,
    SI jp,
    SI jt,
    SI jt1,
    SD colh2o,
    SD colco2,
    SD colo3,
    SD coln2o,
    SD colco,
    SD colch4,
    SD colo2,
    SD colbrd,
    SD fac00,
    SD fac01,
    SD fac10,
    SD fac11,
    SD rat_h2oco2,
    SD rat_h2oco2_1,
    SD rat_h2oo3,
    SD rat_h2oo3_1,
    SD rat_h2on2o,
    SD rat_h2on2o_1,
    SD rat_h2och4,
    SD rat_h2och4_1,
    SD rat_n2oco2,
    SD rat_n2oco2_1,
    SD rat_o3co2,
    SD rat_o3co2_1,
    SD selffac,
    SD selffrac,
    SI indself,
    SD forfac,
    SD forfrac,
    SI indfor,
    SD minorfrac,
    SD scaleminor,
    SD scaleminorn2,
    SI indminor,
    SD2 fracs,
    SD2 taug,
    const RrtmgLwTables& tables,
    int lay0, int lay1  // layers [lay0, lay1) of the column
) {
  constexpr double oneminus = 1.0 - 1.0e-6;
  const auto& b = tables.b06;
  int lay, ind0, ind1, inds, indf, indm, ig;
  double chi_co2, ratco2, adjfac, adjcolco2;
  double tauself, taufor, absco2;
  for (int lay = lay0; lay < (laytrop < lay1 ? laytrop : lay1); ++lay) {
  chi_co2 = colco2[lay]/(coldry[lay]);
  ratco2 = 1.e20*chi_co2/get_chi_mls(tables, 2, jp[lay]+1);
  if (ratco2  >  3.0) {
  adjfac = 2.0+std::pow((ratco2-2.0), 0.77);
  adjcolco2 = adjfac*get_chi_mls(tables, 2, jp[lay]+1)*coldry[lay]*1.e-20;
  } else {
  adjcolco2 = colco2[lay];
  }
  ind0 = ((jp[lay]-1)*5+(jt[lay]-1))*1 + 1;
  ind1 = (jp[lay]*5+(jt1[lay]-1))*1 + 1;
  inds = indself[lay];
  indf = indfor[lay];
  indm = indminor[lay];
  for (int ig = 1; ig <= 8; ++ig) {
  tauself = selffac[lay] * (b.selfref[(ig - 1) * 10 + (inds - 1)] + selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds+1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)]));
  taufor =  forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  absco2 =  (b.ka_mco2[(ig - 1) * 19 + (indm - 1)] + minorfrac[lay] * (b.ka_mco2[(ig - 1) * 19 + (indm+1 - 1)] - b.ka_mco2[(ig - 1) * 19 + (indm - 1)]));
  taug[lay][68+ig - 1] = colh2o[lay] * (fac00[lay] * b.ka[(ig - 1) * 65 + (ind0 - 1)] + fac10[lay] * b.ka[(ig - 1) * 65 + (ind0+1 - 1)] + fac01[lay] * b.ka[(ig - 1) * 65 + (ind1 - 1)] + fac11[lay] * b.ka[(ig - 1) * 65 + (ind1+1 - 1)]) + tauself + taufor + adjcolco2 * absco2 + wx[2 - 1][lay] * b.cfc11adj[ig - 1] + wx[3 - 1][lay] * b.cfc12[ig - 1];
  fracs[lay][68+ig - 1] = b.fracrefa[ig - 1];
  }
  }
  for (int lay = (laytrop > lay0 ? laytrop : lay0); lay < (nlayers < lay1 ? nlayers : lay1); ++lay) {
  for (int ig = 1; ig <= 8; ++ig) {
  taug[lay][68+ig - 1] = 0.0 + wx[2 - 1][lay] * b.cfc11adj[ig - 1] + wx[3 - 1][lay] * b.cfc12[ig - 1];
  fracs[lay][68+ig - 1] = b.fracrefa[ig - 1];
  }
  }
}

CPPWRF_HD inline void rrtmg_lw_taugb07(
    int nlayers,
    SD pavel,
    SD2 wx,
    SD coldry,
    int laytrop,
    SI jp,
    SI jt,
    SI jt1,
    SD colh2o,
    SD colco2,
    SD colo3,
    SD coln2o,
    SD colco,
    SD colch4,
    SD colo2,
    SD colbrd,
    SD fac00,
    SD fac01,
    SD fac10,
    SD fac11,
    SD rat_h2oco2,
    SD rat_h2oco2_1,
    SD rat_h2oo3,
    SD rat_h2oo3_1,
    SD rat_h2on2o,
    SD rat_h2on2o_1,
    SD rat_h2och4,
    SD rat_h2och4_1,
    SD rat_n2oco2,
    SD rat_n2oco2_1,
    SD rat_o3co2,
    SD rat_o3co2_1,
    SD selffac,
    SD selffrac,
    SI indself,
    SD forfac,
    SD forfrac,
    SI indfor,
    SD minorfrac,
    SD scaleminor,
    SD scaleminorn2,
    SI indminor,
    SD2 fracs,
    SD2 taug,
    const RrtmgLwTables& tables,
    int lay0, int lay1  // layers [lay0, lay1) of the column
) {
  constexpr double oneminus = 1.0 - 1.0e-6;
  const auto& b = tables.b07;
  int lay, ind0, ind1, inds, indf, indm, ig;
  int js, js1, jmco2, jpl;
  double speccomb, specparm, specmult, fs;
  double speccomb1, specparm1, specmult1, fs1;
  double speccomb_mco2, specparm_mco2, specmult_mco2, fmco2;
  double speccomb_planck, specparm_planck, specmult_planck, fpl;
  double p, p4, fk0, fk1, fk2;
  double fac000, fac100, fac200, fac010, fac110, fac210;
  double fac001, fac101, fac201, fac011, fac111, fac211;
  double tauself, taufor, co2m1, co2m2, absco2;
  double chi_co2, ratco2, adjfac, adjcolco2;
  double refrat_planck_a, refrat_m_a;
  double tau_major, tau_major1;
  refrat_planck_a = get_chi_mls(tables, 1, 3)/get_chi_mls(tables, 3, 3);
  refrat_m_a = get_chi_mls(tables, 1, 3)/get_chi_mls(tables, 3, 3);
  for (int lay = lay0; lay < (laytrop < lay1 ? laytrop : lay1); ++lay) {
  speccomb = colh2o[lay] + rat_h2oo3[lay]*colo3[lay];
  specparm = colh2o[lay]/speccomb;
  if (specparm  >=  oneminus) specparm = oneminus;
  specmult = 8.*(specparm);
  js = 1 + static_cast<int>(specmult);
  fs = std::fmod(specmult, 1.0);
  speccomb1 = colh2o[lay] + rat_h2oo3_1[lay]*colo3[lay];
  specparm1 = colh2o[lay]/speccomb1;
  if (specparm1  >=  oneminus) specparm1 = oneminus;
  specmult1 = 8.*(specparm1);
  js1 = 1 + static_cast<int>(specmult1);
  fs1 = std::fmod(specmult1, 1.0);
  speccomb_mco2 = colh2o[lay] + refrat_m_a*colo3[lay];
  specparm_mco2 = colh2o[lay]/speccomb_mco2;
  if (specparm_mco2  >=  oneminus) specparm_mco2 = oneminus;
  specmult_mco2 = 8.*specparm_mco2;
  jmco2 = 1 + static_cast<int>(specmult_mco2);
  fmco2 = std::fmod(specmult_mco2, 1.0);
  chi_co2 = colco2[lay]/(coldry[lay]);
  ratco2 = 1.e20*chi_co2/get_chi_mls(tables, 2, jp[lay]+1);
  if (ratco2  >  3.0) {
  adjfac = 3.0+std::pow((ratco2-3.0), 0.79);
  adjcolco2 = adjfac*get_chi_mls(tables, 2, jp[lay]+1)*coldry[lay]*1.e-20;
  } else {
  adjcolco2 = colco2[lay];
  }
  speccomb_planck = colh2o[lay]+refrat_planck_a*colo3[lay];
  specparm_planck = colh2o[lay]/speccomb_planck;
  if (specparm_planck  >=  oneminus) specparm_planck=oneminus;
  specmult_planck = 8.*specparm_planck;
  jpl= 1 + static_cast<int>(specmult_planck);
  fpl = std::fmod(specmult_planck, 1.0);
  ind0 = ((jp[lay]-1)*5+(jt[lay]-1))*9 + js;
  ind1 = (jp[lay]*5+(jt1[lay]-1))*9 + js1;
  inds = indself[lay];
  indf = indfor[lay];
  indm = indminor[lay];
  if (specparm  <  0.125) {
  p = fs - 1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac000 = fk0*fac00[lay];
  fac100 = fk1*fac00[lay];
  fac200 = fk2*fac00[lay];
  fac010 = fk0*fac10[lay];
  fac110 = fk1*fac10[lay];
  fac210 = fk2*fac10[lay];
  } else if (specparm  >  0.875) {
  p = -fs;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac000 = fk0*fac00[lay];
  fac100 = fk1*fac00[lay];
  fac200 = fk2*fac00[lay];
  fac010 = fk0*fac10[lay];
  fac110 = fk1*fac10[lay];
  fac210 = fk2*fac10[lay];
  } else {
  fac000 = (1. - fs) * fac00[lay];
  fac010 = (1. - fs) * fac10[lay];
  fac100 = fs * fac00[lay];
  fac110 = fs * fac10[lay];
  }
  if (specparm1  <  0.125) {
  p = fs1 - 1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac001 = fk0*fac01[lay];
  fac101 = fk1*fac01[lay];
  fac201 = fk2*fac01[lay];
  fac011 = fk0*fac11[lay];
  fac111 = fk1*fac11[lay];
  fac211 = fk2*fac11[lay];
  } else if (specparm1  >  0.875) {
  p = -fs1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac001 = fk0*fac01[lay];
  fac101 = fk1*fac01[lay];
  fac201 = fk2*fac01[lay];
  fac011 = fk0*fac11[lay];
  fac111 = fk1*fac11[lay];
  fac211 = fk2*fac11[lay];
  } else {
  fac001 = (1. - fs1) * fac01[lay];
  fac011 = (1. - fs1) * fac11[lay];
  fac101 = fs1 * fac01[lay];
  fac111 = fs1 * fac11[lay];
  }
  for (int ig = 1; ig <= 12; ++ig) {
  tauself = selffac[lay]* (b.selfref[(ig - 1) * 10 + (inds - 1)] + selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds+1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)]));
  taufor = forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  co2m1 = b.ka_mco2[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmco2 - 1)] + fmco2 * (b.ka_mco2[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmco2+1 - 1)] - b.ka_mco2[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmco2 - 1)]);
  co2m2 = b.ka_mco2[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmco2 - 1)] + fmco2 * (b.ka_mco2[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmco2+1 - 1)] - b.ka_mco2[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmco2 - 1)]);
  absco2 = co2m1 + minorfrac[lay] * (co2m2 - co2m1);
  if (specparm  <  0.125) {
  tau_major = speccomb * (fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac200 * b.ka[(ig - 1) * 585 + (ind0+2 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)] + fac210 * b.ka[(ig - 1) * 585 + (ind0+11 - 1)]);
  } else if (specparm  >  0.875) {
  tau_major = speccomb * (fac200 * b.ka[(ig - 1) * 585 + (ind0-1 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac000 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac210 * b.ka[(ig - 1) * 585 + (ind0+8 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)]);
  } else {
  tau_major = speccomb * (fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)]);
  }
  if (specparm1  <  0.125) {
  tau_major1 = speccomb1 * (fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac201 * b.ka[(ig - 1) * 585 + (ind1+2 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)] + fac211 * b.ka[(ig - 1) * 585 + (ind1+11 - 1)]);
  } else if (specparm1  >  0.875) {
  tau_major1 = speccomb1 * (fac201 * b.ka[(ig - 1) * 585 + (ind1-1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac001 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac211 * b.ka[(ig - 1) * 585 + (ind1+8 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)]);
  } else {
  tau_major1 = speccomb1 * (fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)]);
  }
  taug[lay][76+ig - 1] = tau_major + tau_major1 + tauself + taufor + adjcolco2*absco2;
  fracs[lay][76+ig - 1] = b.fracrefa[(jpl - 1) * 12 + (ig - 1)] + fpl * (b.fracrefa[(jpl+1 - 1) * 12 + (ig - 1)]-b.fracrefa[(jpl - 1) * 12 + (ig - 1)]);
  }
  }
  for (int lay = (laytrop > lay0 ? laytrop : lay0); lay < (nlayers < lay1 ? nlayers : lay1); ++lay) {
  chi_co2 = colco2[lay]/(coldry[lay]);
  ratco2 = 1.e20*chi_co2/get_chi_mls(tables, 2, jp[lay]+1);
  if (ratco2  >  3.0) {
  adjfac = 2.0+std::pow((ratco2-2.0), 0.79);
  adjcolco2 = adjfac*get_chi_mls(tables, 2, jp[lay]+1)*coldry[lay]*1.e-20;
  } else {
  adjcolco2 = colco2[lay];
  }
  ind0 = ((jp[lay]-13)*5+(jt[lay]-1))*1 + 1;
  ind1 = ((jp[lay]-12)*5+(jt1[lay]-1))*1 + 1;
  indm = indminor[lay];
  for (int ig = 1; ig <= 12; ++ig) {
  absco2 = b.kb_mco2[(ig - 1) * 19 + (indm - 1)] + minorfrac[lay] * (b.kb_mco2[(ig - 1) * 19 + (indm+1 - 1)] - b.kb_mco2[(ig - 1) * 19 + (indm - 1)]);
  taug[lay][76+ig - 1] = colo3[lay] * (fac00[lay] * b.kb[(ig - 1) * 235 + (ind0 - 1)] + fac10[lay] * b.kb[(ig - 1) * 235 + (ind0+1 - 1)] + fac01[lay] * b.kb[(ig - 1) * 235 + (ind1 - 1)] + fac11[lay] * b.kb[(ig - 1) * 235 + (ind1+1 - 1)]) + adjcolco2 * absco2;
  fracs[lay][76+ig - 1] = b.fracrefb[ig - 1];
  }
  taug[lay][76+6 - 1]=taug[lay][76+6 - 1]*0.92;
  taug[lay][76+7 - 1]=taug[lay][76+7 - 1]*0.88;
  taug[lay][76+8 - 1]=taug[lay][76+8 - 1]*1.07;
  taug[lay][76+9 - 1]=taug[lay][76+9 - 1]*1.1;
  taug[lay][76+10 - 1]=taug[lay][76+10 - 1]*0.99;
  taug[lay][76+11 - 1]=taug[lay][76+11 - 1]*0.855;
  }
}

CPPWRF_HD inline void rrtmg_lw_taugb08(
    int nlayers,
    SD pavel,
    SD2 wx,
    SD coldry,
    int laytrop,
    SI jp,
    SI jt,
    SI jt1,
    SD colh2o,
    SD colco2,
    SD colo3,
    SD coln2o,
    SD colco,
    SD colch4,
    SD colo2,
    SD colbrd,
    SD fac00,
    SD fac01,
    SD fac10,
    SD fac11,
    SD rat_h2oco2,
    SD rat_h2oco2_1,
    SD rat_h2oo3,
    SD rat_h2oo3_1,
    SD rat_h2on2o,
    SD rat_h2on2o_1,
    SD rat_h2och4,
    SD rat_h2och4_1,
    SD rat_n2oco2,
    SD rat_n2oco2_1,
    SD rat_o3co2,
    SD rat_o3co2_1,
    SD selffac,
    SD selffrac,
    SI indself,
    SD forfac,
    SD forfrac,
    SI indfor,
    SD minorfrac,
    SD scaleminor,
    SD scaleminorn2,
    SI indminor,
    SD2 fracs,
    SD2 taug,
    const RrtmgLwTables& tables,
    int lay0, int lay1  // layers [lay0, lay1) of the column
) {
  constexpr double oneminus = 1.0 - 1.0e-6;
  const auto& b = tables.b08;
  int lay, ind0, ind1, inds, indf, indm, ig;
  double tauself, taufor, absco2, abso3, absn2o;
  double chi_co2, ratco2, adjfac, adjcolco2;
  for (int lay = lay0; lay < (laytrop < lay1 ? laytrop : lay1); ++lay) {
  chi_co2 = colco2[lay]/(coldry[lay]);
  ratco2 = 1.e20*chi_co2/get_chi_mls(tables, 2, jp[lay]+1);
  if (ratco2  >  3.0) {
  adjfac = 2.0+std::pow((ratco2-2.0), 0.65);
  adjcolco2 = adjfac*get_chi_mls(tables, 2, jp[lay]+1)*coldry[lay]*1.e-20;
  } else {
  adjcolco2 = colco2[lay];
  }
  ind0 = ((jp[lay]-1)*5+(jt[lay]-1))*1 + 1;
  ind1 = (jp[lay]*5+(jt1[lay]-1))*1 + 1;
  inds = indself[lay];
  indf = indfor[lay];
  indm = indminor[lay];
  for (int ig = 1; ig <= 8; ++ig) {
  tauself = selffac[lay] * (b.selfref[(ig - 1) * 10 + (inds - 1)] + selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds+1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)]));
  taufor = forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  absco2 =  (b.ka_mco2[(ig - 1) * 19 + (indm - 1)] + minorfrac[lay] * (b.ka_mco2[(ig - 1) * 19 + (indm+1 - 1)] - b.ka_mco2[(ig - 1) * 19 + (indm - 1)]));
  abso3 =  (b.ka_mo3[(ig - 1) * 19 + (indm - 1)] + minorfrac[lay] * (b.ka_mo3[(ig - 1) * 19 + (indm+1 - 1)] - b.ka_mo3[(ig - 1) * 19 + (indm - 1)]));
  absn2o =  (b.ka_mn2o[(ig - 1) * 19 + (indm - 1)] + minorfrac[lay] * (b.ka_mn2o[(ig - 1) * 19 + (indm+1 - 1)] - b.ka_mn2o[(ig - 1) * 19 + (indm - 1)]));
  taug[lay][88+ig - 1] = colh2o[lay] * (fac00[lay] * b.ka[(ig - 1) * 65 + (ind0 - 1)] + fac10[lay] * b.ka[(ig - 1) * 65 + (ind0+1 - 1)] + fac01[lay] * b.ka[(ig - 1) * 65 + (ind1 - 1)] + fac11[lay] * b.ka[(ig - 1) * 65 + (ind1+1 - 1)]) + tauself + taufor + adjcolco2*absco2 + colo3[lay] * abso3 + coln2o[lay] * absn2o + wx[3 - 1][lay] * b.cfc12[ig - 1] + wx[4 - 1][lay] * b.cfc22adj[ig - 1];
  fracs[lay][88+ig - 1] = b.fracrefa[ig - 1];
  }
  }
  for (int lay = (laytrop > lay0 ? laytrop : lay0); lay < (nlayers < lay1 ? nlayers : lay1); ++lay) {
  chi_co2 = colco2[lay]/coldry[lay];
  ratco2 = 1.e20*chi_co2/get_chi_mls(tables, 2, jp[lay]+1);
  if (ratco2  >  3.0) {
  adjfac = 2.0+std::pow((ratco2-2.0), 0.65);
  adjcolco2 = adjfac*get_chi_mls(tables, 2, jp[lay]+1) * coldry[lay]*1.e-20;
  } else {
  adjcolco2 = colco2[lay];
  }
  ind0 = ((jp[lay]-13)*5+(jt[lay]-1))*1 + 1;
  ind1 = ((jp[lay]-12)*5+(jt1[lay]-1))*1 + 1;
  indm = indminor[lay];
  for (int ig = 1; ig <= 8; ++ig) {
  absco2 =  (b.kb_mco2[(ig - 1) * 19 + (indm - 1)] + minorfrac[lay] * (b.kb_mco2[(ig - 1) * 19 + (indm+1 - 1)] - b.kb_mco2[(ig - 1) * 19 + (indm - 1)]));
  absn2o =  (b.kb_mn2o[(ig - 1) * 19 + (indm - 1)] + minorfrac[lay] * (b.kb_mn2o[(ig - 1) * 19 + (indm+1 - 1)] - b.kb_mn2o[(ig - 1) * 19 + (indm - 1)]));
  taug[lay][88+ig - 1] = colo3[lay] * (fac00[lay] * b.kb[(ig - 1) * 235 + (ind0 - 1)] + fac10[lay] * b.kb[(ig - 1) * 235 + (ind0+1 - 1)] + fac01[lay] * b.kb[(ig - 1) * 235 + (ind1 - 1)] + fac11[lay] * b.kb[(ig - 1) * 235 + (ind1+1 - 1)]) + adjcolco2*absco2 + coln2o[lay]*absn2o + wx[3 - 1][lay] * b.cfc12[ig - 1] + wx[4 - 1][lay] * b.cfc22adj[ig - 1];
  fracs[lay][88+ig - 1] = b.fracrefb[ig - 1];
  }
  }
}

CPPWRF_HD inline void rrtmg_lw_taugb09(
    int nlayers,
    SD pavel,
    SD2 wx,
    SD coldry,
    int laytrop,
    SI jp,
    SI jt,
    SI jt1,
    SD colh2o,
    SD colco2,
    SD colo3,
    SD coln2o,
    SD colco,
    SD colch4,
    SD colo2,
    SD colbrd,
    SD fac00,
    SD fac01,
    SD fac10,
    SD fac11,
    SD rat_h2oco2,
    SD rat_h2oco2_1,
    SD rat_h2oo3,
    SD rat_h2oo3_1,
    SD rat_h2on2o,
    SD rat_h2on2o_1,
    SD rat_h2och4,
    SD rat_h2och4_1,
    SD rat_n2oco2,
    SD rat_n2oco2_1,
    SD rat_o3co2,
    SD rat_o3co2_1,
    SD selffac,
    SD selffrac,
    SI indself,
    SD forfac,
    SD forfrac,
    SI indfor,
    SD minorfrac,
    SD scaleminor,
    SD scaleminorn2,
    SI indminor,
    SD2 fracs,
    SD2 taug,
    const RrtmgLwTables& tables,
    int lay0, int lay1  // layers [lay0, lay1) of the column
) {
  constexpr double oneminus = 1.0 - 1.0e-6;
  const auto& b = tables.b09;
  int lay, ind0, ind1, inds, indf, indm, ig;
  int js, js1, jmn2o, jpl;
  double speccomb, specparm, specmult, fs;
  double speccomb1, specparm1, specmult1, fs1;
  double speccomb_mn2o, specparm_mn2o, specmult_mn2o, fmn2o;
  double speccomb_planck, specparm_planck, specmult_planck, fpl;
  double p, p4, fk0, fk1, fk2;
  double fac000, fac100, fac200, fac010, fac110, fac210;
  double fac001, fac101, fac201, fac011, fac111, fac211;
  double tauself, taufor, n2om1, n2om2, absn2o;
  double chi_n2o, ratn2o, adjfac, adjcoln2o;
  double refrat_planck_a, refrat_m_a;
  double tau_major, tau_major1;
  refrat_planck_a = get_chi_mls(tables, 1, 9)/get_chi_mls(tables, 6, 9);
  refrat_m_a = get_chi_mls(tables, 1, 3)/get_chi_mls(tables, 6, 3);
  for (int lay = lay0; lay < (laytrop < lay1 ? laytrop : lay1); ++lay) {
  speccomb = colh2o[lay] + rat_h2och4[lay]*colch4[lay];
  specparm = colh2o[lay]/speccomb;
  if (specparm  >=  oneminus) specparm = oneminus;
  specmult = 8.*(specparm);
  js = 1 + static_cast<int>(specmult);
  fs = std::fmod(specmult, 1.0);
  speccomb1 = colh2o[lay] + rat_h2och4_1[lay]*colch4[lay];
  specparm1 = colh2o[lay]/speccomb1;
  if (specparm1  >=  oneminus) specparm1 = oneminus;
  specmult1 = 8.*(specparm1);
  js1 = 1 + static_cast<int>(specmult1);
  fs1 = std::fmod(specmult1, 1.0);
  speccomb_mn2o = colh2o[lay] + refrat_m_a*colch4[lay];
  specparm_mn2o = colh2o[lay]/speccomb_mn2o;
  if (specparm_mn2o  >=  oneminus) specparm_mn2o = oneminus;
  specmult_mn2o = 8.*specparm_mn2o;
  jmn2o = 1 + static_cast<int>(specmult_mn2o);
  fmn2o = std::fmod(specmult_mn2o, 1.0);
  chi_n2o = coln2o[lay]/(coldry[lay]);
  ratn2o = 1.e20*chi_n2o/get_chi_mls(tables, 4, jp[lay]+1);
  if (ratn2o  >  1.5) {
  adjfac = 0.5+std::pow((ratn2o-0.5), 0.65);
  adjcoln2o = adjfac*get_chi_mls(tables, 4, jp[lay]+1)*coldry[lay]*1.e-20;
  } else {
  adjcoln2o = coln2o[lay];
  }
  speccomb_planck = colh2o[lay]+refrat_planck_a*colch4[lay];
  specparm_planck = colh2o[lay]/speccomb_planck;
  if (specparm_planck  >=  oneminus) specparm_planck=oneminus;
  specmult_planck = 8.*specparm_planck;
  jpl= 1 + static_cast<int>(specmult_planck);
  fpl = std::fmod(specmult_planck, 1.0);
  ind0 = ((jp[lay]-1)*5+(jt[lay]-1))*9 + js;
  ind1 = (jp[lay]*5+(jt1[lay]-1))*9 + js1;
  inds = indself[lay];
  indf = indfor[lay];
  indm = indminor[lay];
  if (specparm  <  0.125) {
  p = fs - 1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac000 = fk0*fac00[lay];
  fac100 = fk1*fac00[lay];
  fac200 = fk2*fac00[lay];
  fac010 = fk0*fac10[lay];
  fac110 = fk1*fac10[lay];
  fac210 = fk2*fac10[lay];
  } else if (specparm  >  0.875) {
  p = -fs;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac000 = fk0*fac00[lay];
  fac100 = fk1*fac00[lay];
  fac200 = fk2*fac00[lay];
  fac010 = fk0*fac10[lay];
  fac110 = fk1*fac10[lay];
  fac210 = fk2*fac10[lay];
  } else {
  fac000 = (1. - fs) * fac00[lay];
  fac010 = (1. - fs) * fac10[lay];
  fac100 = fs * fac00[lay];
  fac110 = fs * fac10[lay];
  }
  if (specparm1  <  0.125) {
  p = fs1 - 1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac001 = fk0*fac01[lay];
  fac101 = fk1*fac01[lay];
  fac201 = fk2*fac01[lay];
  fac011 = fk0*fac11[lay];
  fac111 = fk1*fac11[lay];
  fac211 = fk2*fac11[lay];
  } else if (specparm1  >  0.875) {
  p = -fs1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac001 = fk0*fac01[lay];
  fac101 = fk1*fac01[lay];
  fac201 = fk2*fac01[lay];
  fac011 = fk0*fac11[lay];
  fac111 = fk1*fac11[lay];
  fac211 = fk2*fac11[lay];
  } else {
  fac001 = (1. - fs1) * fac01[lay];
  fac011 = (1. - fs1) * fac11[lay];
  fac101 = fs1 * fac01[lay];
  fac111 = fs1 * fac11[lay];
  }
  for (int ig = 1; ig <= 12; ++ig) {
  tauself = selffac[lay]* (b.selfref[(ig - 1) * 10 + (inds - 1)] + selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds+1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)]));
  taufor = forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  n2om1 = b.ka_mn2o[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmn2o - 1)] + fmn2o * (b.ka_mn2o[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmn2o+1 - 1)] - b.ka_mn2o[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmn2o - 1)]);
  n2om2 = b.ka_mn2o[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmn2o - 1)] + fmn2o * (b.ka_mn2o[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmn2o+1 - 1)] - b.ka_mn2o[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmn2o - 1)]);
  absn2o = n2om1 + minorfrac[lay] * (n2om2 - n2om1);
  if (specparm  <  0.125) {
  tau_major = speccomb * (fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac200 * b.ka[(ig - 1) * 585 + (ind0+2 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)] + fac210 * b.ka[(ig - 1) * 585 + (ind0+11 - 1)]);
  } else if (specparm  >  0.875) {
  tau_major = speccomb * (fac200 * b.ka[(ig - 1) * 585 + (ind0-1 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac000 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac210 * b.ka[(ig - 1) * 585 + (ind0+8 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)]);
  } else {
  tau_major = speccomb * (fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)]);
  }
  if (specparm1  <  0.125) {
  tau_major1 = speccomb1 * (fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac201 * b.ka[(ig - 1) * 585 + (ind1+2 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)] + fac211 * b.ka[(ig - 1) * 585 + (ind1+11 - 1)]);
  } else if (specparm1  >  0.875) {
  tau_major1 = speccomb1 * (fac201 * b.ka[(ig - 1) * 585 + (ind1-1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac001 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac211 * b.ka[(ig - 1) * 585 + (ind1+8 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)]);
  } else {
  tau_major1 = speccomb1 * (fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)]);
  }
  taug[lay][96+ig - 1] = tau_major + tau_major1 + tauself + taufor + adjcoln2o*absn2o;
  fracs[lay][96+ig - 1] = b.fracrefa[(jpl - 1) * 12 + (ig - 1)] + fpl * (b.fracrefa[(jpl+1 - 1) * 12 + (ig - 1)]-b.fracrefa[(jpl - 1) * 12 + (ig - 1)]);
  }
  }
  for (int lay = (laytrop > lay0 ? laytrop : lay0); lay < (nlayers < lay1 ? nlayers : lay1); ++lay) {
  chi_n2o = coln2o[lay]/(coldry[lay]);
  ratn2o = 1.e20*chi_n2o/get_chi_mls(tables, 4, jp[lay]+1);
  if (ratn2o  >  1.5) {
  adjfac = 0.5+std::pow((ratn2o-0.5), 0.65);
  adjcoln2o = adjfac*get_chi_mls(tables, 4, jp[lay]+1)*coldry[lay]*1.e-20;
  } else {
  adjcoln2o = coln2o[lay];
  }
  ind0 = ((jp[lay]-13)*5+(jt[lay]-1))*1 + 1;
  ind1 = ((jp[lay]-12)*5+(jt1[lay]-1))*1 + 1;
  indm = indminor[lay];
  for (int ig = 1; ig <= 12; ++ig) {
  absn2o = b.kb_mn2o[(ig - 1) * 19 + (indm - 1)] + minorfrac[lay] * (b.kb_mn2o[(ig - 1) * 19 + (indm+1 - 1)] - b.kb_mn2o[(ig - 1) * 19 + (indm - 1)]);
  taug[lay][96+ig - 1] = colch4[lay] * (fac00[lay] * b.kb[(ig - 1) * 235 + (ind0 - 1)] + fac10[lay] * b.kb[(ig - 1) * 235 + (ind0+1 - 1)] + fac01[lay] * b.kb[(ig - 1) * 235 + (ind1 - 1)] + fac11[lay] * b.kb[(ig - 1) * 235 + (ind1+1 - 1)]) + adjcoln2o*absn2o;
  fracs[lay][96+ig - 1] = b.fracrefb[ig - 1];
  }
  }
}

CPPWRF_HD inline void rrtmg_lw_taugb10(
    int nlayers,
    SD pavel,
    SD2 wx,
    SD coldry,
    int laytrop,
    SI jp,
    SI jt,
    SI jt1,
    SD colh2o,
    SD colco2,
    SD colo3,
    SD coln2o,
    SD colco,
    SD colch4,
    SD colo2,
    SD colbrd,
    SD fac00,
    SD fac01,
    SD fac10,
    SD fac11,
    SD rat_h2oco2,
    SD rat_h2oco2_1,
    SD rat_h2oo3,
    SD rat_h2oo3_1,
    SD rat_h2on2o,
    SD rat_h2on2o_1,
    SD rat_h2och4,
    SD rat_h2och4_1,
    SD rat_n2oco2,
    SD rat_n2oco2_1,
    SD rat_o3co2,
    SD rat_o3co2_1,
    SD selffac,
    SD selffrac,
    SI indself,
    SD forfac,
    SD forfrac,
    SI indfor,
    SD minorfrac,
    SD scaleminor,
    SD scaleminorn2,
    SI indminor,
    SD2 fracs,
    SD2 taug,
    const RrtmgLwTables& tables,
    int lay0, int lay1  // layers [lay0, lay1) of the column
) {
  constexpr double oneminus = 1.0 - 1.0e-6;
  const auto& b = tables.b10;
  int lay, ind0, ind1, inds, indf, ig;
  double tauself, taufor;
  for (int lay = lay0; lay < (laytrop < lay1 ? laytrop : lay1); ++lay) {
  ind0 = ((jp[lay]-1)*5+(jt[lay]-1))*1 + 1;
  ind1 = (jp[lay]*5+(jt1[lay]-1))*1 + 1;
  inds = indself[lay];
  indf = indfor[lay];
  for (int ig = 1; ig <= 6; ++ig) {
  tauself = selffac[lay] * (b.selfref[(ig - 1) * 10 + (inds - 1)] + selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds+1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)]));
  taufor = forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  taug[lay][108+ig - 1] = colh2o[lay] * (fac00[lay] * b.ka[(ig - 1) * 65 + (ind0 - 1)] + fac10[lay] * b.ka[(ig - 1) * 65 + (ind0+1 - 1)] + fac01[lay] * b.ka[(ig - 1) * 65 + (ind1 - 1)] + fac11[lay] * b.ka[(ig - 1) * 65 + (ind1+1 - 1)]) + tauself + taufor;
  fracs[lay][108+ig - 1] = b.fracrefa[ig - 1];
  }
  }
  for (int lay = (laytrop > lay0 ? laytrop : lay0); lay < (nlayers < lay1 ? nlayers : lay1); ++lay) {
  ind0 = ((jp[lay]-13)*5+(jt[lay]-1))*1 + 1;
  ind1 = ((jp[lay]-12)*5+(jt1[lay]-1))*1 + 1;
  indf = indfor[lay];
  for (int ig = 1; ig <= 6; ++ig) {
  taufor = forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  taug[lay][108+ig - 1] = colh2o[lay] * (fac00[lay] * b.kb[(ig - 1) * 235 + (ind0 - 1)] + fac10[lay] * b.kb[(ig - 1) * 235 + (ind0+1 - 1)] + fac01[lay] * b.kb[(ig - 1) * 235 + (ind1 - 1)] + fac11[lay] * b.kb[(ig - 1) * 235 + (ind1+1 - 1)]) + taufor;
  fracs[lay][108+ig - 1] = b.fracrefb[ig - 1];
  }
  }
}

CPPWRF_HD inline void rrtmg_lw_taugb11(
    int nlayers,
    SD pavel,
    SD2 wx,
    SD coldry,
    int laytrop,
    SI jp,
    SI jt,
    SI jt1,
    SD colh2o,
    SD colco2,
    SD colo3,
    SD coln2o,
    SD colco,
    SD colch4,
    SD colo2,
    SD colbrd,
    SD fac00,
    SD fac01,
    SD fac10,
    SD fac11,
    SD rat_h2oco2,
    SD rat_h2oco2_1,
    SD rat_h2oo3,
    SD rat_h2oo3_1,
    SD rat_h2on2o,
    SD rat_h2on2o_1,
    SD rat_h2och4,
    SD rat_h2och4_1,
    SD rat_n2oco2,
    SD rat_n2oco2_1,
    SD rat_o3co2,
    SD rat_o3co2_1,
    SD selffac,
    SD selffrac,
    SI indself,
    SD forfac,
    SD forfrac,
    SI indfor,
    SD minorfrac,
    SD scaleminor,
    SD scaleminorn2,
    SI indminor,
    SD2 fracs,
    SD2 taug,
    const RrtmgLwTables& tables,
    int lay0, int lay1  // layers [lay0, lay1) of the column
) {
  constexpr double oneminus = 1.0 - 1.0e-6;
  const auto& b = tables.b11;
  int lay, ind0, ind1, inds, indf, indm, ig;
  double scaleo2, tauself, taufor, tauo2;
  for (int lay = lay0; lay < (laytrop < lay1 ? laytrop : lay1); ++lay) {
  ind0 = ((jp[lay]-1)*5+(jt[lay]-1))*1 + 1;
  ind1 = (jp[lay]*5+(jt1[lay]-1))*1 + 1;
  inds = indself[lay];
  indf = indfor[lay];
  indm = indminor[lay];
  scaleo2 = colo2[lay]*scaleminor[lay];
  for (int ig = 1; ig <= 8; ++ig) {
  tauself = selffac[lay] * (b.selfref[(ig - 1) * 10 + (inds - 1)] + selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds+1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)]));
  taufor = forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  tauo2 =  scaleo2 * (b.ka_mo2[(ig - 1) * 19 + (indm - 1)] + minorfrac[lay] * (b.ka_mo2[(ig - 1) * 19 + (indm+1 - 1)] - b.ka_mo2[(ig - 1) * 19 + (indm - 1)]));
  taug[lay][114+ig - 1] = colh2o[lay] * (fac00[lay] * b.ka[(ig - 1) * 65 + (ind0 - 1)] + fac10[lay] * b.ka[(ig - 1) * 65 + (ind0+1 - 1)] + fac01[lay] * b.ka[(ig - 1) * 65 + (ind1 - 1)] + fac11[lay] * b.ka[(ig - 1) * 65 + (ind1+1 - 1)]) + tauself + taufor + tauo2;
  fracs[lay][114+ig - 1] = b.fracrefa[ig - 1];
  }
  }
  for (int lay = (laytrop > lay0 ? laytrop : lay0); lay < (nlayers < lay1 ? nlayers : lay1); ++lay) {
  ind0 = ((jp[lay]-13)*5+(jt[lay]-1))*1 + 1;
  ind1 = ((jp[lay]-12)*5+(jt1[lay]-1))*1 + 1;
  indf = indfor[lay];
  indm = indminor[lay];
  scaleo2 = colo2[lay]*scaleminor[lay];
  for (int ig = 1; ig <= 8; ++ig) {
  taufor = forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  tauo2 =  scaleo2 * (b.kb_mo2[(ig - 1) * 19 + (indm - 1)] + minorfrac[lay] * (b.kb_mo2[(ig - 1) * 19 + (indm+1 - 1)] - b.kb_mo2[(ig - 1) * 19 + (indm - 1)]));
  taug[lay][114+ig - 1] = colh2o[lay] * (fac00[lay] * b.kb[(ig - 1) * 235 + (ind0 - 1)] + fac10[lay] * b.kb[(ig - 1) * 235 + (ind0+1 - 1)] + fac01[lay] * b.kb[(ig - 1) * 235 + (ind1 - 1)] + fac11[lay] * b.kb[(ig - 1) * 235 + (ind1+1 - 1)]) + taufor + tauo2;
  fracs[lay][114+ig - 1] = b.fracrefb[ig - 1];
  }
  }
}

CPPWRF_HD inline void rrtmg_lw_taugb12(
    int nlayers,
    SD pavel,
    SD2 wx,
    SD coldry,
    int laytrop,
    SI jp,
    SI jt,
    SI jt1,
    SD colh2o,
    SD colco2,
    SD colo3,
    SD coln2o,
    SD colco,
    SD colch4,
    SD colo2,
    SD colbrd,
    SD fac00,
    SD fac01,
    SD fac10,
    SD fac11,
    SD rat_h2oco2,
    SD rat_h2oco2_1,
    SD rat_h2oo3,
    SD rat_h2oo3_1,
    SD rat_h2on2o,
    SD rat_h2on2o_1,
    SD rat_h2och4,
    SD rat_h2och4_1,
    SD rat_n2oco2,
    SD rat_n2oco2_1,
    SD rat_o3co2,
    SD rat_o3co2_1,
    SD selffac,
    SD selffrac,
    SI indself,
    SD forfac,
    SD forfrac,
    SI indfor,
    SD minorfrac,
    SD scaleminor,
    SD scaleminorn2,
    SI indminor,
    SD2 fracs,
    SD2 taug,
    const RrtmgLwTables& tables,
    int lay0, int lay1  // layers [lay0, lay1) of the column
) {
  constexpr double oneminus = 1.0 - 1.0e-6;
  const auto& b = tables.b12;
  int lay, ind0, ind1, inds, indf, ig;
  int js, js1, jpl;
  double speccomb, specparm, specmult, fs;
  double speccomb1, specparm1, specmult1, fs1;
  double speccomb_planck, specparm_planck, specmult_planck, fpl;
  double p, p4, fk0, fk1, fk2;
  double fac000, fac100, fac200, fac010, fac110, fac210;
  double fac001, fac101, fac201, fac011, fac111, fac211;
  double tauself, taufor;
  double refrat_planck_a;
  double tau_major, tau_major1;
  refrat_planck_a = get_chi_mls(tables, 1, 10)/get_chi_mls(tables, 2, 10);
  for (int lay = lay0; lay < (laytrop < lay1 ? laytrop : lay1); ++lay) {
  speccomb = colh2o[lay] + rat_h2oco2[lay]*colco2[lay];
  specparm = colh2o[lay]/speccomb;
  if (specparm  >=  oneminus) specparm = oneminus;
  specmult = 8.*(specparm);
  js = 1 + static_cast<int>(specmult);
  fs = std::fmod(specmult, 1.0);
  speccomb1 = colh2o[lay] + rat_h2oco2_1[lay]*colco2[lay];
  specparm1 = colh2o[lay]/speccomb1;
  if (specparm1  >=  oneminus) specparm1 = oneminus;
  specmult1 = 8.*(specparm1);
  js1 = 1 + static_cast<int>(specmult1);
  fs1 = std::fmod(specmult1, 1.0);
  speccomb_planck = colh2o[lay]+refrat_planck_a*colco2[lay];
  specparm_planck = colh2o[lay]/speccomb_planck;
  if (specparm_planck  >=  oneminus) specparm_planck=oneminus;
  specmult_planck = 8.*specparm_planck;
  jpl= 1 + static_cast<int>(specmult_planck);
  fpl = std::fmod(specmult_planck, 1.0);
  ind0 = ((jp[lay]-1)*5+(jt[lay]-1))*9 + js;
  ind1 = (jp[lay]*5+(jt1[lay]-1))*9 + js1;
  inds = indself[lay];
  indf = indfor[lay];
  if (specparm  <  0.125) {
  p = fs - 1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac000 = fk0*fac00[lay];
  fac100 = fk1*fac00[lay];
  fac200 = fk2*fac00[lay];
  fac010 = fk0*fac10[lay];
  fac110 = fk1*fac10[lay];
  fac210 = fk2*fac10[lay];
  } else if (specparm  >  0.875) {
  p = -fs;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac000 = fk0*fac00[lay];
  fac100 = fk1*fac00[lay];
  fac200 = fk2*fac00[lay];
  fac010 = fk0*fac10[lay];
  fac110 = fk1*fac10[lay];
  fac210 = fk2*fac10[lay];
  } else {
  fac000 = (1. - fs) * fac00[lay];
  fac010 = (1. - fs) * fac10[lay];
  fac100 = fs * fac00[lay];
  fac110 = fs * fac10[lay];
  }
  if (specparm1  <  0.125) {
  p = fs1 - 1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac001 = fk0*fac01[lay];
  fac101 = fk1*fac01[lay];
  fac201 = fk2*fac01[lay];
  fac011 = fk0*fac11[lay];
  fac111 = fk1*fac11[lay];
  fac211 = fk2*fac11[lay];
  } else if (specparm1  >  0.875) {
  p = -fs1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac001 = fk0*fac01[lay];
  fac101 = fk1*fac01[lay];
  fac201 = fk2*fac01[lay];
  fac011 = fk0*fac11[lay];
  fac111 = fk1*fac11[lay];
  fac211 = fk2*fac11[lay];
  } else {
  fac001 = (1. - fs1) * fac01[lay];
  fac011 = (1. - fs1) * fac11[lay];
  fac101 = fs1 * fac01[lay];
  fac111 = fs1 * fac11[lay];
  }
  for (int ig = 1; ig <= 8; ++ig) {
  tauself = selffac[lay]* (b.selfref[(ig - 1) * 10 + (inds - 1)] + selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds+1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)]));
  taufor = forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  if (specparm  <  0.125) {
  tau_major = speccomb * (fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac200 * b.ka[(ig - 1) * 585 + (ind0+2 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)] + fac210 * b.ka[(ig - 1) * 585 + (ind0+11 - 1)]);
  } else if (specparm  >  0.875) {
  tau_major = speccomb * (fac200 * b.ka[(ig - 1) * 585 + (ind0-1 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac000 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac210 * b.ka[(ig - 1) * 585 + (ind0+8 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)]);
  } else {
  tau_major = speccomb * (fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)]);
  }
  if (specparm1  <  0.125) {
  tau_major1 = speccomb1 * (fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac201 * b.ka[(ig - 1) * 585 + (ind1+2 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)] + fac211 * b.ka[(ig - 1) * 585 + (ind1+11 - 1)]);
  } else if (specparm1  >  0.875) {
  tau_major1 = speccomb1 * (fac201 * b.ka[(ig - 1) * 585 + (ind1-1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac001 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac211 * b.ka[(ig - 1) * 585 + (ind1+8 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)]);
  } else {
  tau_major1 = speccomb1 * (fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)]);
  }
  taug[lay][122+ig - 1] = tau_major + tau_major1 + tauself + taufor;
  fracs[lay][122+ig - 1] = b.fracrefa[(jpl - 1) * 8 + (ig - 1)] + fpl * (b.fracrefa[(jpl+1 - 1) * 8 + (ig - 1)]-b.fracrefa[(jpl - 1) * 8 + (ig - 1)]);
  }
  }
  for (int lay = (laytrop > lay0 ? laytrop : lay0); lay < (nlayers < lay1 ? nlayers : lay1); ++lay) {
  for (int ig = 1; ig <= 8; ++ig) {
  taug[lay][122+ig - 1] = 0.0;
  fracs[lay][122+ig - 1] = 0.0;
  }
  }
}

CPPWRF_HD inline void rrtmg_lw_taugb13(
    int nlayers,
    SD pavel,
    SD2 wx,
    SD coldry,
    int laytrop,
    SI jp,
    SI jt,
    SI jt1,
    SD colh2o,
    SD colco2,
    SD colo3,
    SD coln2o,
    SD colco,
    SD colch4,
    SD colo2,
    SD colbrd,
    SD fac00,
    SD fac01,
    SD fac10,
    SD fac11,
    SD rat_h2oco2,
    SD rat_h2oco2_1,
    SD rat_h2oo3,
    SD rat_h2oo3_1,
    SD rat_h2on2o,
    SD rat_h2on2o_1,
    SD rat_h2och4,
    SD rat_h2och4_1,
    SD rat_n2oco2,
    SD rat_n2oco2_1,
    SD rat_o3co2,
    SD rat_o3co2_1,
    SD selffac,
    SD selffrac,
    SI indself,
    SD forfac,
    SD forfrac,
    SI indfor,
    SD minorfrac,
    SD scaleminor,
    SD scaleminorn2,
    SI indminor,
    SD2 fracs,
    SD2 taug,
    const RrtmgLwTables& tables,
    int lay0, int lay1  // layers [lay0, lay1) of the column
) {
  constexpr double oneminus = 1.0 - 1.0e-6;
  const auto& b = tables.b13;
  int lay, ind0, ind1, inds, indf, indm, ig;
  int js, js1, jmco2, jmco, jpl;
  double speccomb, specparm, specmult, fs;
  double speccomb1, specparm1, specmult1, fs1;
  double speccomb_mco2, specparm_mco2, specmult_mco2, fmco2;
  double speccomb_mco, specparm_mco, specmult_mco, fmco;
  double speccomb_planck, specparm_planck, specmult_planck, fpl;
  double p, p4, fk0, fk1, fk2;
  double fac000, fac100, fac200, fac010, fac110, fac210;
  double fac001, fac101, fac201, fac011, fac111, fac211;
  double tauself, taufor, co2m1, co2m2, absco2;
  double com1, com2, absco, abso3;
  double chi_co2, ratco2, adjfac, adjcolco2;
  double refrat_planck_a, refrat_m_a, refrat_m_a3;
  double tau_major, tau_major1;
  refrat_planck_a = get_chi_mls(tables, 1, 5)/get_chi_mls(tables, 4, 5);
  refrat_m_a = get_chi_mls(tables, 1, 1)/get_chi_mls(tables, 4, 1);
  refrat_m_a3 = get_chi_mls(tables, 1, 3)/get_chi_mls(tables, 4, 3);
  for (int lay = lay0; lay < (laytrop < lay1 ? laytrop : lay1); ++lay) {
  speccomb = colh2o[lay] + rat_h2on2o[lay]*coln2o[lay];
  specparm = colh2o[lay]/speccomb;
  if (specparm  >=  oneminus) specparm = oneminus;
  specmult = 8.*(specparm);
  js = 1 + static_cast<int>(specmult);
  fs = std::fmod(specmult, 1.0);
  speccomb1 = colh2o[lay] + rat_h2on2o_1[lay]*coln2o[lay];
  specparm1 = colh2o[lay]/speccomb1;
  if (specparm1  >=  oneminus) specparm1 = oneminus;
  specmult1 = 8.*(specparm1);
  js1 = 1 + static_cast<int>(specmult1);
  fs1 = std::fmod(specmult1, 1.0);
  speccomb_mco2 = colh2o[lay] + refrat_m_a*coln2o[lay];
  specparm_mco2 = colh2o[lay]/speccomb_mco2;
  if (specparm_mco2  >=  oneminus) specparm_mco2 = oneminus;
  specmult_mco2 = 8.*specparm_mco2;
  jmco2 = 1 + static_cast<int>(specmult_mco2);
  fmco2 = std::fmod(specmult_mco2, 1.0);
  chi_co2 = colco2[lay]/(coldry[lay]);
  ratco2 = 1.e20*chi_co2/3.55e-4;
  if (ratco2  >  3.0) {
  adjfac = 2.0+std::pow((ratco2-2.0), 0.68);
  adjcolco2 = adjfac*3.55e-4*coldry[lay]*1.e-20;
  } else {
  adjcolco2 = colco2[lay];
  }
  speccomb_mco = colh2o[lay] + refrat_m_a3*coln2o[lay];
  specparm_mco = colh2o[lay]/speccomb_mco;
  if (specparm_mco  >=  oneminus) specparm_mco = oneminus;
  specmult_mco = 8.*specparm_mco;
  jmco = 1 + static_cast<int>(specmult_mco);
  fmco = std::fmod(specmult_mco, 1.0);
  speccomb_planck = colh2o[lay]+refrat_planck_a*coln2o[lay];
  specparm_planck = colh2o[lay]/speccomb_planck;
  if (specparm_planck  >=  oneminus) specparm_planck=oneminus;
  specmult_planck = 8.*specparm_planck;
  jpl= 1 + static_cast<int>(specmult_planck);
  fpl = std::fmod(specmult_planck, 1.0);
  ind0 = ((jp[lay]-1)*5+(jt[lay]-1))*9 + js;
  ind1 = (jp[lay]*5+(jt1[lay]-1))*9 + js1;
  inds = indself[lay];
  indf = indfor[lay];
  indm = indminor[lay];
  if (specparm  <  0.125) {
  p = fs - 1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac000 = fk0*fac00[lay];
  fac100 = fk1*fac00[lay];
  fac200 = fk2*fac00[lay];
  fac010 = fk0*fac10[lay];
  fac110 = fk1*fac10[lay];
  fac210 = fk2*fac10[lay];
  } else if (specparm  >  0.875) {
  p = -fs;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac000 = fk0*fac00[lay];
  fac100 = fk1*fac00[lay];
  fac200 = fk2*fac00[lay];
  fac010 = fk0*fac10[lay];
  fac110 = fk1*fac10[lay];
  fac210 = fk2*fac10[lay];
  } else {
  fac000 = (1. - fs) * fac00[lay];
  fac010 = (1. - fs) * fac10[lay];
  fac100 = fs * fac00[lay];
  fac110 = fs * fac10[lay];
  }
  if (specparm1  <  0.125) {
  p = fs1 - 1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac001 = fk0*fac01[lay];
  fac101 = fk1*fac01[lay];
  fac201 = fk2*fac01[lay];
  fac011 = fk0*fac11[lay];
  fac111 = fk1*fac11[lay];
  fac211 = fk2*fac11[lay];
  } else if (specparm1  >  0.875) {
  p = -fs1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac001 = fk0*fac01[lay];
  fac101 = fk1*fac01[lay];
  fac201 = fk2*fac01[lay];
  fac011 = fk0*fac11[lay];
  fac111 = fk1*fac11[lay];
  fac211 = fk2*fac11[lay];
  } else {
  fac001 = (1. - fs1) * fac01[lay];
  fac011 = (1. - fs1) * fac11[lay];
  fac101 = fs1 * fac01[lay];
  fac111 = fs1 * fac11[lay];
  }
  for (int ig = 1; ig <= 4; ++ig) {
  tauself = selffac[lay]* (b.selfref[(ig - 1) * 10 + (inds - 1)] + selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds+1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)]));
  taufor = forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  co2m1 = b.ka_mco2[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmco2 - 1)] + fmco2 * (b.ka_mco2[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmco2+1 - 1)] - b.ka_mco2[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmco2 - 1)]);
  co2m2 = b.ka_mco2[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmco2 - 1)] + fmco2 * (b.ka_mco2[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmco2+1 - 1)] - b.ka_mco2[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmco2 - 1)]);
  absco2 = co2m1 + minorfrac[lay] * (co2m2 - co2m1);
  com1 = b.ka_mco[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmco - 1)] + fmco * (b.ka_mco[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmco+1 - 1)] - b.ka_mco[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmco - 1)]);
  com2 = b.ka_mco[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmco - 1)] + fmco * (b.ka_mco[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmco+1 - 1)] - b.ka_mco[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmco - 1)]);
  absco = com1 + minorfrac[lay] * (com2 - com1);
  if (specparm  <  0.125) {
  tau_major = speccomb * (fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac200 * b.ka[(ig - 1) * 585 + (ind0+2 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)] + fac210 * b.ka[(ig - 1) * 585 + (ind0+11 - 1)]);
  } else if (specparm  >  0.875) {
  tau_major = speccomb * (fac200 * b.ka[(ig - 1) * 585 + (ind0-1 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac000 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac210 * b.ka[(ig - 1) * 585 + (ind0+8 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)]);
  } else {
  tau_major = speccomb * (fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)]);
  }
  if (specparm1  <  0.125) {
  tau_major1 = speccomb1 * (fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac201 * b.ka[(ig - 1) * 585 + (ind1+2 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)] + fac211 * b.ka[(ig - 1) * 585 + (ind1+11 - 1)]);
  } else if (specparm1  >  0.875) {
  tau_major1 = speccomb1 * (fac201 * b.ka[(ig - 1) * 585 + (ind1-1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac001 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac211 * b.ka[(ig - 1) * 585 + (ind1+8 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)]);
  } else {
  tau_major1 = speccomb1 * (fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)]);
  }
  taug[lay][130+ig - 1] = tau_major + tau_major1 + tauself + taufor + adjcolco2*absco2 + colco[lay]*absco;
  fracs[lay][130+ig - 1] = b.fracrefa[(jpl - 1) * 4 + (ig - 1)] + fpl * (b.fracrefa[(jpl+1 - 1) * 4 + (ig - 1)]-b.fracrefa[(jpl - 1) * 4 + (ig - 1)]);
  }
  }
  for (int lay = (laytrop > lay0 ? laytrop : lay0); lay < (nlayers < lay1 ? nlayers : lay1); ++lay) {
  indm = indminor[lay];
  for (int ig = 1; ig <= 4; ++ig) {
  abso3 = b.kb_mo3[(ig - 1) * 19 + (indm - 1)] + minorfrac[lay] * (b.kb_mo3[(ig - 1) * 19 + (indm+1 - 1)] - b.kb_mo3[(ig - 1) * 19 + (indm - 1)]);
  taug[lay][130+ig - 1] = colo3[lay]*abso3;
  fracs[lay][130+ig - 1] =  b.fracrefb[ig - 1];
  }
  }
}

CPPWRF_HD inline void rrtmg_lw_taugb14(
    int nlayers,
    SD pavel,
    SD2 wx,
    SD coldry,
    int laytrop,
    SI jp,
    SI jt,
    SI jt1,
    SD colh2o,
    SD colco2,
    SD colo3,
    SD coln2o,
    SD colco,
    SD colch4,
    SD colo2,
    SD colbrd,
    SD fac00,
    SD fac01,
    SD fac10,
    SD fac11,
    SD rat_h2oco2,
    SD rat_h2oco2_1,
    SD rat_h2oo3,
    SD rat_h2oo3_1,
    SD rat_h2on2o,
    SD rat_h2on2o_1,
    SD rat_h2och4,
    SD rat_h2och4_1,
    SD rat_n2oco2,
    SD rat_n2oco2_1,
    SD rat_o3co2,
    SD rat_o3co2_1,
    SD selffac,
    SD selffrac,
    SI indself,
    SD forfac,
    SD forfrac,
    SI indfor,
    SD minorfrac,
    SD scaleminor,
    SD scaleminorn2,
    SI indminor,
    SD2 fracs,
    SD2 taug,
    const RrtmgLwTables& tables,
    int lay0, int lay1  // layers [lay0, lay1) of the column
) {
  constexpr double oneminus = 1.0 - 1.0e-6;
  const auto& b = tables.b14;
  int lay, ind0, ind1, inds, indf, ig;
  double tauself, taufor;
  for (int lay = lay0; lay < (laytrop < lay1 ? laytrop : lay1); ++lay) {
  ind0 = ((jp[lay]-1)*5+(jt[lay]-1))*1 + 1;
  ind1 = (jp[lay]*5+(jt1[lay]-1))*1 + 1;
  inds = indself[lay];
  indf = indfor[lay];
  for (int ig = 1; ig <= 2; ++ig) {
  tauself = selffac[lay] * (b.selfref[(ig - 1) * 10 + (inds - 1)] + selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds+1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)]));
  taufor =  forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  taug[lay][134+ig - 1] = colco2[lay] * (fac00[lay] * b.ka[(ig - 1) * 65 + (ind0 - 1)] + fac10[lay] * b.ka[(ig - 1) * 65 + (ind0+1 - 1)] + fac01[lay] * b.ka[(ig - 1) * 65 + (ind1 - 1)] + fac11[lay] * b.ka[(ig - 1) * 65 + (ind1+1 - 1)]) + tauself + taufor;
  fracs[lay][134+ig - 1] = b.fracrefa[ig - 1];
  }
  }
  for (int lay = (laytrop > lay0 ? laytrop : lay0); lay < (nlayers < lay1 ? nlayers : lay1); ++lay) {
  ind0 = ((jp[lay]-13)*5+(jt[lay]-1))*1 + 1;
  ind1 = ((jp[lay]-12)*5+(jt1[lay]-1))*1 + 1;
  for (int ig = 1; ig <= 2; ++ig) {
  taug[lay][134+ig - 1] = colco2[lay] * (fac00[lay] * b.kb[(ig - 1) * 235 + (ind0 - 1)] + fac10[lay] * b.kb[(ig - 1) * 235 + (ind0+1 - 1)] + fac01[lay] * b.kb[(ig - 1) * 235 + (ind1 - 1)] + fac11[lay] * b.kb[(ig - 1) * 235 + (ind1+1 - 1)]);
  fracs[lay][134+ig - 1] = b.fracrefb[ig - 1];
  }
  }
}

CPPWRF_HD inline void rrtmg_lw_taugb15(
    int nlayers,
    SD pavel,
    SD2 wx,
    SD coldry,
    int laytrop,
    SI jp,
    SI jt,
    SI jt1,
    SD colh2o,
    SD colco2,
    SD colo3,
    SD coln2o,
    SD colco,
    SD colch4,
    SD colo2,
    SD colbrd,
    SD fac00,
    SD fac01,
    SD fac10,
    SD fac11,
    SD rat_h2oco2,
    SD rat_h2oco2_1,
    SD rat_h2oo3,
    SD rat_h2oo3_1,
    SD rat_h2on2o,
    SD rat_h2on2o_1,
    SD rat_h2och4,
    SD rat_h2och4_1,
    SD rat_n2oco2,
    SD rat_n2oco2_1,
    SD rat_o3co2,
    SD rat_o3co2_1,
    SD selffac,
    SD selffrac,
    SI indself,
    SD forfac,
    SD forfrac,
    SI indfor,
    SD minorfrac,
    SD scaleminor,
    SD scaleminorn2,
    SI indminor,
    SD2 fracs,
    SD2 taug,
    const RrtmgLwTables& tables,
    int lay0, int lay1  // layers [lay0, lay1) of the column
) {
  constexpr double oneminus = 1.0 - 1.0e-6;
  const auto& b = tables.b15;
  int lay, ind0, ind1, inds, indf, indm, ig;
  int js, js1, jmn2, jpl;
  double speccomb, specparm, specmult, fs;
  double speccomb1, specparm1, specmult1, fs1;
  double speccomb_mn2, specparm_mn2, specmult_mn2, fmn2;
  double speccomb_planck, specparm_planck, specmult_planck, fpl;
  double p, p4, fk0, fk1, fk2;
  double fac000, fac100, fac200, fac010, fac110, fac210;
  double fac001, fac101, fac201, fac011, fac111, fac211;
  double scalen2, tauself, taufor, n2m1, n2m2, taun2;
  double refrat_planck_a, refrat_m_a;
  double tau_major, tau_major1;
  refrat_planck_a = get_chi_mls(tables, 4, 1)/get_chi_mls(tables, 2, 1);
  refrat_m_a = get_chi_mls(tables, 4, 1)/get_chi_mls(tables, 2, 1);
  for (int lay = lay0; lay < (laytrop < lay1 ? laytrop : lay1); ++lay) {
  speccomb = coln2o[lay] + rat_n2oco2[lay]*colco2[lay];
  specparm = coln2o[lay]/speccomb;
  if (specparm  >=  oneminus) specparm = oneminus;
  specmult = 8.*(specparm);
  js = 1 + static_cast<int>(specmult);
  fs = std::fmod(specmult, 1.0);
  speccomb1 = coln2o[lay] + rat_n2oco2_1[lay]*colco2[lay];
  specparm1 = coln2o[lay]/speccomb1;
  if (specparm1  >=  oneminus) specparm1 = oneminus;
  specmult1 = 8.*(specparm1);
  js1 = 1 + static_cast<int>(specmult1);
  fs1 = std::fmod(specmult1, 1.0);
  speccomb_mn2 = coln2o[lay] + refrat_m_a*colco2[lay];
  specparm_mn2 = coln2o[lay]/speccomb_mn2;
  if (specparm_mn2  >=  oneminus) specparm_mn2 = oneminus;
  specmult_mn2 = 8.*specparm_mn2;
  jmn2 = 1 + static_cast<int>(specmult_mn2);
  fmn2 = std::fmod(specmult_mn2, 1.0);
  speccomb_planck = coln2o[lay]+refrat_planck_a*colco2[lay];
  specparm_planck = coln2o[lay]/speccomb_planck;
  if (specparm_planck  >=  oneminus) specparm_planck=oneminus;
  specmult_planck = 8.*specparm_planck;
  jpl= 1 + static_cast<int>(specmult_planck);
  fpl = std::fmod(specmult_planck, 1.0);
  ind0 = ((jp[lay]-1)*5+(jt[lay]-1))*9 + js;
  ind1 = (jp[lay]*5+(jt1[lay]-1))*9 + js1;
  inds = indself[lay];
  indf = indfor[lay];
  indm = indminor[lay];
  scalen2 = colbrd[lay]*scaleminor[lay];
  if (specparm  <  0.125) {
  p = fs - 1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac000 = fk0*fac00[lay];
  fac100 = fk1*fac00[lay];
  fac200 = fk2*fac00[lay];
  fac010 = fk0*fac10[lay];
  fac110 = fk1*fac10[lay];
  fac210 = fk2*fac10[lay];
  } else if (specparm  >  0.875) {
  p = -fs;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac000 = fk0*fac00[lay];
  fac100 = fk1*fac00[lay];
  fac200 = fk2*fac00[lay];
  fac010 = fk0*fac10[lay];
  fac110 = fk1*fac10[lay];
  fac210 = fk2*fac10[lay];
  } else {
  fac000 = (1. - fs) * fac00[lay];
  fac010 = (1. - fs) * fac10[lay];
  fac100 = fs * fac00[lay];
  fac110 = fs * fac10[lay];
  }
  if (specparm1  <  0.125) {
  p = fs1 - 1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac001 = fk0*fac01[lay];
  fac101 = fk1*fac01[lay];
  fac201 = fk2*fac01[lay];
  fac011 = fk0*fac11[lay];
  fac111 = fk1*fac11[lay];
  fac211 = fk2*fac11[lay];
  } else if (specparm1  >  0.875) {
  p = -fs1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac001 = fk0*fac01[lay];
  fac101 = fk1*fac01[lay];
  fac201 = fk2*fac01[lay];
  fac011 = fk0*fac11[lay];
  fac111 = fk1*fac11[lay];
  fac211 = fk2*fac11[lay];
  } else {
  fac001 = (1. - fs1) * fac01[lay];
  fac011 = (1. - fs1) * fac11[lay];
  fac101 = fs1 * fac01[lay];
  fac111 = fs1 * fac11[lay];
  }
  for (int ig = 1; ig <= 2; ++ig) {
  tauself = selffac[lay]* (b.selfref[(ig - 1) * 10 + (inds - 1)] + selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds+1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)]));
  taufor =  forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  n2m1 = b.ka_mn2[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmn2 - 1)] + fmn2 * (b.ka_mn2[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmn2+1 - 1)] - b.ka_mn2[(ig - 1) * (19 * 9) + (indm - 1) * 9 + (jmn2 - 1)]);
  n2m2 = b.ka_mn2[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmn2 - 1)] + fmn2 * (b.ka_mn2[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmn2+1 - 1)] - b.ka_mn2[(ig - 1) * (19 * 9) + (indm+1 - 1) * 9 + (jmn2 - 1)]);
  taun2 = scalen2 * (n2m1 + minorfrac[lay] * (n2m2 - n2m1));
  if (specparm  <  0.125) {
  tau_major = speccomb * (fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac200 * b.ka[(ig - 1) * 585 + (ind0+2 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)] + fac210 * b.ka[(ig - 1) * 585 + (ind0+11 - 1)]);
  } else if (specparm  >  0.875) {
  tau_major = speccomb * (fac200 * b.ka[(ig - 1) * 585 + (ind0-1 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac000 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac210 * b.ka[(ig - 1) * 585 + (ind0+8 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)]);
  } else {
  tau_major = speccomb * (fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)]);
  }
  if (specparm1  <  0.125) {
  tau_major1 = speccomb1 * (fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac201 * b.ka[(ig - 1) * 585 + (ind1+2 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)] + fac211 * b.ka[(ig - 1) * 585 + (ind1+11 - 1)]);
  } else if (specparm1  >  0.875) {
  tau_major1 = speccomb1 * (fac201 * b.ka[(ig - 1) * 585 + (ind1-1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac001 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac211 * b.ka[(ig - 1) * 585 + (ind1+8 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)]);
  } else {
  tau_major1 = speccomb1 * (fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)]);
  }
  taug[lay][136+ig - 1] = tau_major + tau_major1 + tauself + taufor + taun2;
  fracs[lay][136+ig - 1] = b.fracrefa[(jpl - 1) * 2 + (ig - 1)] + fpl * (b.fracrefa[(jpl+1 - 1) * 2 + (ig - 1)]-b.fracrefa[(jpl - 1) * 2 + (ig - 1)]);
  }
  }
  for (int lay = (laytrop > lay0 ? laytrop : lay0); lay < (nlayers < lay1 ? nlayers : lay1); ++lay) {
  for (int ig = 1; ig <= 2; ++ig) {
  taug[lay][136+ig - 1] = 0.0;
  fracs[lay][136+ig - 1] = 0.0;
  }
  }
}

CPPWRF_HD inline void rrtmg_lw_taugb16(
    int nlayers,
    SD pavel,
    SD2 wx,
    SD coldry,
    int laytrop,
    SI jp,
    SI jt,
    SI jt1,
    SD colh2o,
    SD colco2,
    SD colo3,
    SD coln2o,
    SD colco,
    SD colch4,
    SD colo2,
    SD colbrd,
    SD fac00,
    SD fac01,
    SD fac10,
    SD fac11,
    SD rat_h2oco2,
    SD rat_h2oco2_1,
    SD rat_h2oo3,
    SD rat_h2oo3_1,
    SD rat_h2on2o,
    SD rat_h2on2o_1,
    SD rat_h2och4,
    SD rat_h2och4_1,
    SD rat_n2oco2,
    SD rat_n2oco2_1,
    SD rat_o3co2,
    SD rat_o3co2_1,
    SD selffac,
    SD selffrac,
    SI indself,
    SD forfac,
    SD forfrac,
    SI indfor,
    SD minorfrac,
    SD scaleminor,
    SD scaleminorn2,
    SI indminor,
    SD2 fracs,
    SD2 taug,
    const RrtmgLwTables& tables,
    int lay0, int lay1  // layers [lay0, lay1) of the column
) {
  constexpr double oneminus = 1.0 - 1.0e-6;
  const auto& b = tables.b16;
  int lay, ind0, ind1, inds, indf, ig;
  int js, js1, jpl;
  double speccomb, specparm, specmult, fs;
  double speccomb1, specparm1, specmult1, fs1;
  double speccomb_planck, specparm_planck, specmult_planck, fpl;
  double p, p4, fk0, fk1, fk2;
  double fac000, fac100, fac200, fac010, fac110, fac210;
  double fac001, fac101, fac201, fac011, fac111, fac211;
  double tauself, taufor;
  double refrat_planck_a;
  double tau_major, tau_major1;
  refrat_planck_a = get_chi_mls(tables, 1, 6)/get_chi_mls(tables, 6, 6);
  for (int lay = lay0; lay < (laytrop < lay1 ? laytrop : lay1); ++lay) {
  speccomb = colh2o[lay] + rat_h2och4[lay]*colch4[lay];
  specparm = colh2o[lay]/speccomb;
  if (specparm  >=  oneminus) specparm = oneminus;
  specmult = 8.*(specparm);
  js = 1 + static_cast<int>(specmult);
  fs = std::fmod(specmult, 1.0);
  speccomb1 = colh2o[lay] + rat_h2och4_1[lay]*colch4[lay];
  specparm1 = colh2o[lay]/speccomb1;
  if (specparm1  >=  oneminus) specparm1 = oneminus;
  specmult1 = 8.*(specparm1);
  js1 = 1 + static_cast<int>(specmult1);
  fs1 = std::fmod(specmult1, 1.0);
  speccomb_planck = colh2o[lay]+refrat_planck_a*colch4[lay];
  specparm_planck = colh2o[lay]/speccomb_planck;
  if (specparm_planck  >=  oneminus) specparm_planck=oneminus;
  specmult_planck = 8.*specparm_planck;
  jpl= 1 + static_cast<int>(specmult_planck);
  fpl = std::fmod(specmult_planck, 1.0);
  ind0 = ((jp[lay]-1)*5+(jt[lay]-1))*9 + js;
  ind1 = (jp[lay]*5+(jt1[lay]-1))*9 + js1;
  inds = indself[lay];
  indf = indfor[lay];
  if (specparm  <  0.125) {
  p = fs - 1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac000 = fk0*fac00[lay];
  fac100 = fk1*fac00[lay];
  fac200 = fk2*fac00[lay];
  fac010 = fk0*fac10[lay];
  fac110 = fk1*fac10[lay];
  fac210 = fk2*fac10[lay];
  } else if (specparm  >  0.875) {
  p = -fs;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac000 = fk0*fac00[lay];
  fac100 = fk1*fac00[lay];
  fac200 = fk2*fac00[lay];
  fac010 = fk0*fac10[lay];
  fac110 = fk1*fac10[lay];
  fac210 = fk2*fac10[lay];
  } else {
  fac000 = (1. - fs) * fac00[lay];
  fac010 = (1. - fs) * fac10[lay];
  fac100 = fs * fac00[lay];
  fac110 = fs * fac10[lay];
  }
  if (specparm1  <  0.125) {
  p = fs1 - 1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac001 = fk0*fac01[lay];
  fac101 = fk1*fac01[lay];
  fac201 = fk2*fac01[lay];
  fac011 = fk0*fac11[lay];
  fac111 = fk1*fac11[lay];
  fac211 = fk2*fac11[lay];
  } else if (specparm1  >  0.875) {
  p = -fs1;
  p4 = ((p * p) * (p * p));  // gfortran p**4
  fk0 = p4;
  fk1 = 1 - p - 2.0*p4;
  fk2 = p + p4;
  fac001 = fk0*fac01[lay];
  fac101 = fk1*fac01[lay];
  fac201 = fk2*fac01[lay];
  fac011 = fk0*fac11[lay];
  fac111 = fk1*fac11[lay];
  fac211 = fk2*fac11[lay];
  } else {
  fac001 = (1. - fs1) * fac01[lay];
  fac011 = (1. - fs1) * fac11[lay];
  fac101 = fs1 * fac01[lay];
  fac111 = fs1 * fac11[lay];
  }
  for (int ig = 1; ig <= 2; ++ig) {
  tauself = selffac[lay]* (b.selfref[(ig - 1) * 10 + (inds - 1)] + selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds+1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)]));
  taufor =  forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] + forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf+1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
  if (specparm  <  0.125) {
  tau_major = speccomb * (fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac200 * b.ka[(ig - 1) * 585 + (ind0+2 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)] + fac210 * b.ka[(ig - 1) * 585 + (ind0+11 - 1)]);
  } else if (specparm  >  0.875) {
  tau_major = speccomb * (fac200 * b.ka[(ig - 1) * 585 + (ind0-1 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac000 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac210 * b.ka[(ig - 1) * 585 + (ind0+8 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)]);
  } else {
  tau_major = speccomb * (fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] + fac100 * b.ka[(ig - 1) * 585 + (ind0+1 - 1)] + fac010 * b.ka[(ig - 1) * 585 + (ind0+9 - 1)] + fac110 * b.ka[(ig - 1) * 585 + (ind0+10 - 1)]);
  }
  if (specparm1  <  0.125) {
  tau_major1 = speccomb1 * (fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac201 * b.ka[(ig - 1) * 585 + (ind1+2 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)] + fac211 * b.ka[(ig - 1) * 585 + (ind1+11 - 1)]);
  } else if (specparm1  >  0.875) {
  tau_major1 = speccomb1 * (fac201 * b.ka[(ig - 1) * 585 + (ind1-1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac001 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac211 * b.ka[(ig - 1) * 585 + (ind1+8 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)]);
  } else {
  tau_major1 = speccomb1 * (fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] + fac101 * b.ka[(ig - 1) * 585 + (ind1+1 - 1)] + fac011 * b.ka[(ig - 1) * 585 + (ind1+9 - 1)] + fac111 * b.ka[(ig - 1) * 585 + (ind1+10 - 1)]);
  }
  taug[lay][138+ig - 1] = tau_major + tau_major1 + tauself + taufor;
  fracs[lay][138+ig - 1] = b.fracrefa[(jpl - 1) * 2 + (ig - 1)] + fpl * (b.fracrefa[(jpl+1 - 1) * 2 + (ig - 1)]-b.fracrefa[(jpl - 1) * 2 + (ig - 1)]);
  }
  }
  for (int lay = (laytrop > lay0 ? laytrop : lay0); lay < (nlayers < lay1 ? nlayers : lay1); ++lay) {
  // WRF nspb(16) = 0: the upper-atmosphere lookup always uses absb(1..2, ig) (as in WRF)
  ind0 = ((jp[lay]-13)*5+(jt[lay]-1))*0 + 1;
  ind1 = ((jp[lay]-12)*5+(jt1[lay]-1))*0 + 1;
  for (int ig = 1; ig <= 2; ++ig) {
  taug[lay][138+ig - 1] = colch4[lay] * (fac00[lay] * b.kb[(ig - 1) * 235 + (ind0 - 1)] + fac10[lay] * b.kb[(ig - 1) * 235 + (ind0+1 - 1)] + fac01[lay] * b.kb[(ig - 1) * 235 + (ind1 - 1)] + fac11[lay] * b.kb[(ig - 1) * 235 + (ind1+1 - 1)]);
  fracs[lay][138+ig - 1] = b.fracrefb[ig - 1];
  }
  }
}

CPPWRF_HD inline void rrtmg_lw_taumol(
    int nlayers,
    SD pavel,
    SD2 wx,
    SD coldry,
    int laytrop,
    SI jp,
    SI jt,
    SI jt1,
    SD colh2o,
    SD colco2,
    SD colo3,
    SD coln2o,
    SD colco,
    SD colch4,
    SD colo2,
    SD colbrd,
    SD fac00,
    SD fac01,
    SD fac10,
    SD fac11,
    SD rat_h2oco2,
    SD rat_h2oco2_1,
    SD rat_h2oo3,
    SD rat_h2oo3_1,
    SD rat_h2on2o,
    SD rat_h2on2o_1,
    SD rat_h2och4,
    SD rat_h2och4_1,
    SD rat_n2oco2,
    SD rat_n2oco2_1,
    SD rat_o3co2,
    SD rat_o3co2_1,
    SD selffac,
    SD selffrac,
    SI indself,
    SD forfac,
    SD forfrac,
    SI indfor,
    SD minorfrac,
    SD scaleminor,
    SD scaleminorn2,
    SI indminor,
    SD2 fracs,  // [nlayers][140]
    SD2 taug,   // [nlayers][140]
    const RrtmgLwTables& tables,
    int lay0, int lay1  // layers [lay0, lay1) of the column
) {
  rrtmg_lw_taugb01(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug, tables, lay0, lay1);
  rrtmg_lw_taugb02(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug, tables, lay0, lay1);
  rrtmg_lw_taugb03(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug, tables, lay0, lay1);
  rrtmg_lw_taugb04(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug, tables, lay0, lay1);
  rrtmg_lw_taugb05(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug, tables, lay0, lay1);
  rrtmg_lw_taugb06(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug, tables, lay0, lay1);
  rrtmg_lw_taugb07(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug, tables, lay0, lay1);
  rrtmg_lw_taugb08(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug, tables, lay0, lay1);
  rrtmg_lw_taugb09(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug, tables, lay0, lay1);
  rrtmg_lw_taugb10(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug, tables, lay0, lay1);
  rrtmg_lw_taugb11(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug, tables, lay0, lay1);
  rrtmg_lw_taugb12(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug, tables, lay0, lay1);
  rrtmg_lw_taugb13(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug, tables, lay0, lay1);
  rrtmg_lw_taugb14(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug, tables, lay0, lay1);
  rrtmg_lw_taugb15(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug, tables, lay0, lay1);
  rrtmg_lw_taugb16(nlayers, pavel, wx, coldry, laytrop, jp, jt, jt1, colh2o, colco2, colo3, coln2o, colco, colch4, colo2, colbrd, fac00, fac01, fac10, fac11, rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, selffac, selffrac, indself, forfac, forfrac, indfor, minorfrac, scaleminor, scaleminorn2, indminor, fracs, taug, tables, lay0, lay1);
}

}  // namespace cppwrf::phys
