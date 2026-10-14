#pragma once

#include "core/backend.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"
#include <algorithm>
#include <cmath>

namespace cppwrf::phys {

CPPWRF_HD inline void taumol_sw(
    int nlayers,
    int laytrop,
    const int* jp,
    const int* jt,
    const int* jt1,
    const double* colh2o,
    const double* colco2,
    const double* colch4,
    const double* colo2,
    const double* colo3,
    const double* colmol,
    const double* fac00,
    const double* fac01,
    const double* fac10,
    const double* fac11,
    const double* selffac,
    const double* selffrac,
    const int* indself,
    const double* forfac,
    const double* forfrac,
    const int* indfor,
    double taug[][112],       // (nlayers, 112)
    double taur[][112],       // (nlayers, 112)
    double* sfluxzen,         // (112)
    const RrtmgSwTables& tables
) {
  constexpr double oneminus = 1.0 - 1.0e-6;

  // Band 16 (2600-3250 cm-1): ng16 = 6, offset = 0
  {
    const auto& b = tables.b16;
    for (int lay = 0; lay < laytrop; ++lay) {
      double speccomb = colh2o[lay] + b.strrat1 * colch4[lay];
      double specparm = colh2o[lay] / speccomb;
      if (specparm >= oneminus) specparm = oneminus;
      double specmult = 8.0 * specparm;
      int js = 1 + static_cast<int>(specmult);
      double fs = specmult - static_cast<double>(js - 1);

      double fac000 = (1.0 - fs) * fac00[lay];
      double fac010 = (1.0 - fs) * fac10[lay];
      double fac100 = fs * fac00[lay];
      double fac110 = fs * fac10[lay];
      double fac001 = (1.0 - fs) * fac01[lay];
      double fac011 = (1.0 - fs) * fac11[lay];
      double fac101 = fs * fac01[lay];
      double fac111 = fs * fac11[lay];

      int ind0 = ((jp[lay] - 1) * 5 + (jt[lay] - 1)) * 9 + js;
      int ind1 = (jp[lay] * 5 + (jt1[lay] - 1)) * 9 + js;
      int inds = indself[lay];
      int indf = indfor[lay];
      double tauray = colmol[lay] * b.rayl;

      for (int ig = 1; ig <= 6; ++ig) {
        taug[lay][ig - 1] = speccomb * (
            fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] +
            fac100 * b.ka[(ig - 1) * 585 + (ind0 + 1 - 1)] +
            fac010 * b.ka[(ig - 1) * 585 + (ind0 + 9 - 1)] +
            fac110 * b.ka[(ig - 1) * 585 + (ind0 + 10 - 1)] +
            fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] +
            fac101 * b.ka[(ig - 1) * 585 + (ind1 + 1 - 1)] +
            fac011 * b.ka[(ig - 1) * 585 + (ind1 + 9 - 1)] +
            fac111 * b.ka[(ig - 1) * 585 + (ind1 + 10 - 1)]) +
            colh2o[lay] * (
            selffac[lay] * (b.selfref[(ig - 1) * 10 + (inds - 1)] +
            selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds + 1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)])) +
            forfac[lay] * (b.forref[(ig - 1) * 3 + (indf - 1)] +
            forfrac[lay] * (b.forref[(ig - 1) * 3 + (indf + 1 - 1)] - b.forref[(ig - 1) * 3 + (indf - 1)])));
        taur[lay][ig - 1] = tauray;
      }
    }

    int laysolfr = nlayers - 1;
    for (int lay = laytrop; lay < nlayers; ++lay) {
      if (jp[lay - 1] < b.layreffr && jp[lay] >= b.layreffr) laysolfr = lay;
      int ind0 = ((jp[lay] - 13) * 5 + (jt[lay] - 1)) * 1 + 1;
      int ind1 = ((jp[lay] - 12) * 5 + (jt1[lay] - 1)) * 1 + 1;
      double tauray = colmol[lay] * b.rayl;

      for (int ig = 1; ig <= 6; ++ig) {
        taug[lay][ig - 1] = colch4[lay] * (
            fac00[lay] * b.kb[(ig - 1) * 235 + (ind0 - 1)] +
            fac10[lay] * b.kb[(ig - 1) * 235 + (ind0 + 1 - 1)] +
            fac01[lay] * b.kb[(ig - 1) * 235 + (ind1 - 1)] +
            fac11[lay] * b.kb[(ig - 1) * 235 + (ind1 + 1 - 1)]);
        if (lay == laysolfr) sfluxzen[ig - 1] = b.sfluxref[ig - 1];
        taur[lay][ig - 1] = tauray;
      }
    }
  }

  // Band 17 (3250-4000 cm-1): ng17 = 12, offset = 6
  {
    const auto& b = tables.b17;
    int offset = 6;
    for (int lay = 0; lay < laytrop; ++lay) {
      double speccomb = colh2o[lay] + b.strrat * colco2[lay];
      double specparm = colh2o[lay] / speccomb;
      if (specparm >= oneminus) specparm = oneminus;
      double specmult = 8.0 * specparm;
      int js = 1 + static_cast<int>(specmult);
      double fs = specmult - static_cast<double>(js - 1);

      double fac000 = (1.0 - fs) * fac00[lay];
      double fac010 = (1.0 - fs) * fac10[lay];
      double fac100 = fs * fac00[lay];
      double fac110 = fs * fac10[lay];
      double fac001 = (1.0 - fs) * fac01[lay];
      double fac011 = (1.0 - fs) * fac11[lay];
      double fac101 = fs * fac01[lay];
      double fac111 = fs * fac11[lay];

      int ind0 = ((jp[lay] - 1) * 5 + (jt[lay] - 1)) * 9 + js;
      int ind1 = (jp[lay] * 5 + (jt1[lay] - 1)) * 9 + js;
      int inds = indself[lay];
      int indf = indfor[lay];
      double tauray = colmol[lay] * b.rayl;

      for (int ig = 1; ig <= 12; ++ig) {
        taug[lay][offset + ig - 1] = speccomb * (
            fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] +
            fac100 * b.ka[(ig - 1) * 585 + (ind0 + 1 - 1)] +
            fac010 * b.ka[(ig - 1) * 585 + (ind0 + 9 - 1)] +
            fac110 * b.ka[(ig - 1) * 585 + (ind0 + 10 - 1)] +
            fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] +
            fac101 * b.ka[(ig - 1) * 585 + (ind1 + 1 - 1)] +
            fac011 * b.ka[(ig - 1) * 585 + (ind1 + 9 - 1)] +
            fac111 * b.ka[(ig - 1) * 585 + (ind1 + 10 - 1)]) +
            colh2o[lay] * (
            selffac[lay] * (b.selfref[(ig - 1) * 10 + (inds - 1)] +
            selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds + 1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)])) +
            forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] +
            forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf + 1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)])));
        taur[lay][offset + ig - 1] = tauray;
      }
    }

    int laysolfr = nlayers - 1;
    for (int lay = laytrop; lay < nlayers; ++lay) {
      if (jp[lay - 1] < b.layreffr && jp[lay] >= b.layreffr) laysolfr = lay;
      double speccomb = colh2o[lay] + b.strrat * colco2[lay];
      double specparm = colh2o[lay] / speccomb;
      if (specparm >= oneminus) specparm = oneminus;
      double specmult = 4.0 * specparm;
      int js = 1 + static_cast<int>(specmult);
      double fs = specmult - static_cast<double>(js - 1);

      double fac000 = (1.0 - fs) * fac00[lay];
      double fac010 = (1.0 - fs) * fac10[lay];
      double fac100 = fs * fac00[lay];
      double fac110 = fs * fac10[lay];
      double fac001 = (1.0 - fs) * fac01[lay];
      double fac011 = (1.0 - fs) * fac11[lay];
      double fac101 = fs * fac01[lay];
      double fac111 = fs * fac11[lay];

      int ind0 = ((jp[lay] - 13) * 5 + (jt[lay] - 1)) * 5 + js;
      int ind1 = ((jp[lay] - 12) * 5 + (jt1[lay] - 1)) * 5 + js;
      int indf = indfor[lay];
      double tauray = colmol[lay] * b.rayl;

      for (int ig = 1; ig <= 12; ++ig) {
        taug[lay][offset + ig - 1] = speccomb * (
            fac000 * b.kb[(ig - 1) * 1175 + (ind0 - 1)] +
            fac100 * b.kb[(ig - 1) * 1175 + (ind0 + 1 - 1)] +
            fac010 * b.kb[(ig - 1) * 1175 + (ind0 + 5 - 1)] +
            fac110 * b.kb[(ig - 1) * 1175 + (ind0 + 6 - 1)] +
            fac001 * b.kb[(ig - 1) * 1175 + (ind1 - 1)] +
            fac101 * b.kb[(ig - 1) * 1175 + (ind1 + 1 - 1)] +
            fac011 * b.kb[(ig - 1) * 1175 + (ind1 + 5 - 1)] +
            fac111 * b.kb[(ig - 1) * 1175 + (ind1 + 6 - 1)]) +
            colh2o[lay] * forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] +
            forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf + 1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
        if (lay == laysolfr) {
          sfluxzen[offset + ig - 1] = b.sfluxref[(js - 1) * 12 + (ig - 1)] +
              fs * (b.sfluxref[js * 12 + (ig - 1)] - b.sfluxref[(js - 1) * 12 + (ig - 1)]);
        }
        taur[lay][offset + ig - 1] = tauray;
      }
    }
  }

  // Band 18 (4000-4650 cm-1): ng18 = 8, offset = 18
  {
    const auto& b = tables.b18;
    int offset = 18;
    int laysolfr = laytrop - 1;
    for (int lay = 0; lay < laytrop; ++lay) {
      if (jp[lay] < b.layreffr && jp[lay + 1] >= b.layreffr) laysolfr = std::min(lay + 1, laytrop - 1);
      double speccomb = colh2o[lay] + b.strrat * colch4[lay];
      double specparm = colh2o[lay] / speccomb;
      if (specparm >= oneminus) specparm = oneminus;
      double specmult = 8.0 * specparm;
      int js = 1 + static_cast<int>(specmult);
      double fs = specmult - static_cast<double>(js - 1);

      double fac000 = (1.0 - fs) * fac00[lay];
      double fac010 = (1.0 - fs) * fac10[lay];
      double fac100 = fs * fac00[lay];
      double fac110 = fs * fac10[lay];
      double fac001 = (1.0 - fs) * fac01[lay];
      double fac011 = (1.0 - fs) * fac11[lay];
      double fac101 = fs * fac01[lay];
      double fac111 = fs * fac11[lay];

      int ind0 = ((jp[lay] - 1) * 5 + (jt[lay] - 1)) * 9 + js;
      int ind1 = (jp[lay] * 5 + (jt1[lay] - 1)) * 9 + js;
      int inds = indself[lay];
      int indf = indfor[lay];
      double tauray = colmol[lay] * b.rayl;

      for (int ig = 1; ig <= 8; ++ig) {
        taug[lay][offset + ig - 1] = speccomb * (
            fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] +
            fac100 * b.ka[(ig - 1) * 585 + (ind0 + 1 - 1)] +
            fac010 * b.ka[(ig - 1) * 585 + (ind0 + 9 - 1)] +
            fac110 * b.ka[(ig - 1) * 585 + (ind0 + 10 - 1)] +
            fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] +
            fac101 * b.ka[(ig - 1) * 585 + (ind1 + 1 - 1)] +
            fac011 * b.ka[(ig - 1) * 585 + (ind1 + 9 - 1)] +
            fac111 * b.ka[(ig - 1) * 585 + (ind1 + 10 - 1)]) +
            colh2o[lay] * (
            selffac[lay] * (b.selfref[(ig - 1) * 10 + (inds - 1)] +
            selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds + 1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)])) +
            forfac[lay] * (b.forref[(ig - 1) * 3 + (indf - 1)] +
            forfrac[lay] * (b.forref[(ig - 1) * 3 + (indf + 1 - 1)] - b.forref[(ig - 1) * 3 + (indf - 1)])));
        if (lay == laysolfr) {
          sfluxzen[offset + ig - 1] = b.sfluxref[(js - 1) * 8 + (ig - 1)] +
              fs * (b.sfluxref[js * 8 + (ig - 1)] - b.sfluxref[(js - 1) * 8 + (ig - 1)]);
        }
        taur[lay][offset + ig - 1] = tauray;
      }
    }

    for (int lay = laytrop; lay < nlayers; ++lay) {
      int ind0 = ((jp[lay] - 13) * 5 + (jt[lay] - 1)) * 1 + 1;
      int ind1 = ((jp[lay] - 12) * 5 + (jt1[lay] - 1)) * 1 + 1;
      double tauray = colmol[lay] * b.rayl;

      for (int ig = 1; ig <= 8; ++ig) {
        taug[lay][offset + ig - 1] = colch4[lay] * (
            fac00[lay] * b.kb[(ig - 1) * 235 + (ind0 - 1)] +
            fac10[lay] * b.kb[(ig - 1) * 235 + (ind0 + 1 - 1)] +
            fac01[lay] * b.kb[(ig - 1) * 235 + (ind1 - 1)] +
            fac11[lay] * b.kb[(ig - 1) * 235 + (ind1 + 1 - 1)]);
        taur[lay][offset + ig - 1] = tauray;
      }
    }
  }

  // Band 19 (4650-5150 cm-1): ng19 = 8, offset = 26
  {
    const auto& b = tables.b19;
    int offset = 26;
    int laysolfr = laytrop - 1;
    for (int lay = 0; lay < laytrop; ++lay) {
      if (jp[lay] < b.layreffr && jp[lay + 1] >= b.layreffr) laysolfr = std::min(lay + 1, laytrop - 1);
      double speccomb = colh2o[lay] + b.strrat * colco2[lay];
      double specparm = colh2o[lay] / speccomb;
      if (specparm >= oneminus) specparm = oneminus;
      double specmult = 8.0 * specparm;
      int js = 1 + static_cast<int>(specmult);
      double fs = specmult - static_cast<double>(js - 1);

      double fac000 = (1.0 - fs) * fac00[lay];
      double fac010 = (1.0 - fs) * fac10[lay];
      double fac100 = fs * fac00[lay];
      double fac110 = fs * fac10[lay];
      double fac001 = (1.0 - fs) * fac01[lay];
      double fac011 = (1.0 - fs) * fac11[lay];
      double fac101 = fs * fac01[lay];
      double fac111 = fs * fac11[lay];

      int ind0 = ((jp[lay] - 1) * 5 + (jt[lay] - 1)) * 9 + js;
      int ind1 = (jp[lay] * 5 + (jt1[lay] - 1)) * 9 + js;
      int inds = indself[lay];
      int indf = indfor[lay];
      double tauray = colmol[lay] * b.rayl;

      for (int ig = 1; ig <= 8; ++ig) {
        taug[lay][offset + ig - 1] = speccomb * (
            fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] +
            fac100 * b.ka[(ig - 1) * 585 + (ind0 + 1 - 1)] +
            fac010 * b.ka[(ig - 1) * 585 + (ind0 + 9 - 1)] +
            fac110 * b.ka[(ig - 1) * 585 + (ind0 + 10 - 1)] +
            fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] +
            fac101 * b.ka[(ig - 1) * 585 + (ind1 + 1 - 1)] +
            fac011 * b.ka[(ig - 1) * 585 + (ind1 + 9 - 1)] +
            fac111 * b.ka[(ig - 1) * 585 + (ind1 + 10 - 1)]) +
            colh2o[lay] * (
            selffac[lay] * (b.selfref[(ig - 1) * 10 + (inds - 1)] +
            selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds + 1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)])) +
            forfac[lay] * (b.forref[(ig - 1) * 3 + (indf - 1)] +
            forfrac[lay] * (b.forref[(ig - 1) * 3 + (indf + 1 - 1)] - b.forref[(ig - 1) * 3 + (indf - 1)])));
        if (lay == laysolfr) {
          sfluxzen[offset + ig - 1] = b.sfluxref[(js - 1) * 8 + (ig - 1)] +
              fs * (b.sfluxref[js * 8 + (ig - 1)] - b.sfluxref[(js - 1) * 8 + (ig - 1)]);
        }
        taur[lay][offset + ig - 1] = tauray;
      }
    }

    for (int lay = laytrop; lay < nlayers; ++lay) {
      int ind0 = ((jp[lay] - 13) * 5 + (jt[lay] - 1)) * 1 + 1;
      int ind1 = ((jp[lay] - 12) * 5 + (jt1[lay] - 1)) * 1 + 1;
      double tauray = colmol[lay] * b.rayl;

      for (int ig = 1; ig <= 8; ++ig) {
        taug[lay][offset + ig - 1] = colco2[lay] * (
            fac00[lay] * b.kb[(ig - 1) * 235 + (ind0 - 1)] +
            fac10[lay] * b.kb[(ig - 1) * 235 + (ind0 + 1 - 1)] +
            fac01[lay] * b.kb[(ig - 1) * 235 + (ind1 - 1)] +
            fac11[lay] * b.kb[(ig - 1) * 235 + (ind1 + 1 - 1)]);
        taur[lay][offset + ig - 1] = tauray;
      }
    }
  }

  // Band 20 (5150-6150 cm-1): ng20 = 10, offset = 34
  {
    const auto& b = tables.b20;
    int offset = 34;
    int laysolfr = laytrop - 1;
    for (int lay = 0; lay < laytrop; ++lay) {
      if (jp[lay] < b.layreffr && jp[lay + 1] >= b.layreffr) laysolfr = std::min(lay + 1, laytrop - 1);
      int ind0 = ((jp[lay] - 1) * 5 + (jt[lay] - 1)) * 1 + 1;
      int ind1 = (jp[lay] * 5 + (jt1[lay] - 1)) * 1 + 1;
      int inds = indself[lay];
      int indf = indfor[lay];
      double tauray = colmol[lay] * b.rayl;

      for (int ig = 1; ig <= 10; ++ig) {
        taug[lay][offset + ig - 1] = colh2o[lay] * (
            fac00[lay] * b.ka[(ig - 1) * 65 + (ind0 - 1)] +
            fac10[lay] * b.ka[(ig - 1) * 65 + (ind0 + 1 - 1)] +
            fac01[lay] * b.ka[(ig - 1) * 65 + (ind1 - 1)] +
            fac11[lay] * b.ka[(ig - 1) * 65 + (ind1 + 1 - 1)] +
            selffac[lay] * (b.selfref[(ig - 1) * 10 + (inds - 1)] +
            selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds + 1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)])) +
            forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] +
            forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf + 1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]))) +
            colch4[lay] * b.absch4[ig - 1];
        if (lay == laysolfr) sfluxzen[offset + ig - 1] = b.sfluxref[ig - 1];
        taur[lay][offset + ig - 1] = tauray;
      }
    }

    for (int lay = laytrop; lay < nlayers; ++lay) {
      int ind0 = ((jp[lay] - 13) * 5 + (jt[lay] - 1)) * 1 + 1;
      int ind1 = ((jp[lay] - 12) * 5 + (jt1[lay] - 1)) * 1 + 1;
      int indf = indfor[lay];
      double tauray = colmol[lay] * b.rayl;

      for (int ig = 1; ig <= 10; ++ig) {
        taug[lay][offset + ig - 1] = colh2o[lay] * (
            fac00[lay] * b.kb[(ig - 1) * 235 + (ind0 - 1)] +
            fac10[lay] * b.kb[(ig - 1) * 235 + (ind0 + 1 - 1)] +
            fac01[lay] * b.kb[(ig - 1) * 235 + (ind1 - 1)] +
            fac11[lay] * b.kb[(ig - 1) * 235 + (ind1 + 1 - 1)] +
            forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] +
            forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf + 1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]))) +
            colch4[lay] * b.absch4[ig - 1];
        taur[lay][offset + ig - 1] = tauray;
      }
    }
  }

  // Band 21 (6150-7700 cm-1): ng21 = 10, offset = 44
  {
    const auto& b = tables.b21;
    int offset = 44;
    int laysolfr = laytrop - 1;
    for (int lay = 0; lay < laytrop; ++lay) {
      if (jp[lay] < b.layreffr && jp[lay + 1] >= b.layreffr) laysolfr = std::min(lay + 1, laytrop - 1);
      double speccomb = colh2o[lay] + b.strrat * colco2[lay];
      double specparm = colh2o[lay] / speccomb;
      if (specparm >= oneminus) specparm = oneminus;
      double specmult = 8.0 * specparm;
      int js = 1 + static_cast<int>(specmult);
      double fs = specmult - static_cast<double>(js - 1);

      double fac000 = (1.0 - fs) * fac00[lay];
      double fac010 = (1.0 - fs) * fac10[lay];
      double fac100 = fs * fac00[lay];
      double fac110 = fs * fac10[lay];
      double fac001 = (1.0 - fs) * fac01[lay];
      double fac011 = (1.0 - fs) * fac11[lay];
      double fac101 = fs * fac01[lay];
      double fac111 = fs * fac11[lay];

      int ind0 = ((jp[lay] - 1) * 5 + (jt[lay] - 1)) * 9 + js;
      int ind1 = (jp[lay] * 5 + (jt1[lay] - 1)) * 9 + js;
      int inds = indself[lay];
      int indf = indfor[lay];
      double tauray = colmol[lay] * b.rayl;

      for (int ig = 1; ig <= 10; ++ig) {
        taug[lay][offset + ig - 1] = speccomb * (
            fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] +
            fac100 * b.ka[(ig - 1) * 585 + (ind0 + 1 - 1)] +
            fac010 * b.ka[(ig - 1) * 585 + (ind0 + 9 - 1)] +
            fac110 * b.ka[(ig - 1) * 585 + (ind0 + 10 - 1)] +
            fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] +
            fac101 * b.ka[(ig - 1) * 585 + (ind1 + 1 - 1)] +
            fac011 * b.ka[(ig - 1) * 585 + (ind1 + 9 - 1)] +
            fac111 * b.ka[(ig - 1) * 585 + (ind1 + 10 - 1)]) +
            colh2o[lay] * (
            selffac[lay] * (b.selfref[(ig - 1) * 10 + (inds - 1)] +
            selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds + 1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)])) +
            forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] +
            forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf + 1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)])));
        if (lay == laysolfr) {
          sfluxzen[offset + ig - 1] = b.sfluxref[(js - 1) * 10 + (ig - 1)] +
              fs * (b.sfluxref[js * 10 + (ig - 1)] - b.sfluxref[(js - 1) * 10 + (ig - 1)]);
        }
        taur[lay][offset + ig - 1] = tauray;
      }
    }

    for (int lay = laytrop; lay < nlayers; ++lay) {
      double speccomb = colh2o[lay] + b.strrat * colco2[lay];
      double specparm = colh2o[lay] / speccomb;
      if (specparm >= oneminus) specparm = oneminus;
      double specmult = 4.0 * specparm;
      int js = 1 + static_cast<int>(specmult);
      double fs = specmult - static_cast<double>(js - 1);

      double fac000 = (1.0 - fs) * fac00[lay];
      double fac010 = (1.0 - fs) * fac10[lay];
      double fac100 = fs * fac00[lay];
      double fac110 = fs * fac10[lay];
      double fac001 = (1.0 - fs) * fac01[lay];
      double fac011 = (1.0 - fs) * fac11[lay];
      double fac101 = fs * fac01[lay];
      double fac111 = fs * fac11[lay];

      int ind0 = ((jp[lay] - 13) * 5 + (jt[lay] - 1)) * 5 + js;
      int ind1 = ((jp[lay] - 12) * 5 + (jt1[lay] - 1)) * 5 + js;
      int indf = indfor[lay];
      double tauray = colmol[lay] * b.rayl;

      for (int ig = 1; ig <= 10; ++ig) {
        taug[lay][offset + ig - 1] = speccomb * (
            fac000 * b.kb[(ig - 1) * 1175 + (ind0 - 1)] +
            fac100 * b.kb[(ig - 1) * 1175 + (ind0 + 1 - 1)] +
            fac010 * b.kb[(ig - 1) * 1175 + (ind0 + 5 - 1)] +
            fac110 * b.kb[(ig - 1) * 1175 + (ind0 + 6 - 1)] +
            fac001 * b.kb[(ig - 1) * 1175 + (ind1 - 1)] +
            fac101 * b.kb[(ig - 1) * 1175 + (ind1 + 1 - 1)] +
            fac011 * b.kb[(ig - 1) * 1175 + (ind1 + 5 - 1)] +
            fac111 * b.kb[(ig - 1) * 1175 + (ind1 + 6 - 1)]) +
            colh2o[lay] * forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] +
            forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf + 1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]));
        taur[lay][offset + ig - 1] = tauray;
      }
    }
  }

  // Band 22 (7700-8050 cm-1): ng22 = 2, offset = 54
  {
    const auto& b = tables.b22;
    int offset = 54;
    int laysolfr = laytrop - 1;
    constexpr double o2adj = 1.6;
    for (int lay = 0; lay < laytrop; ++lay) {
      if (jp[lay] < b.layreffr && jp[lay + 1] >= b.layreffr) laysolfr = std::min(lay + 1, laytrop - 1);
      double o2cont = 4.35e-4 * colo2[lay] / (350.0 * 2.0);
      double speccomb = colh2o[lay] + o2adj * b.strrat * colo2[lay];
      double specparm = colh2o[lay] / speccomb;
      if (specparm >= oneminus) specparm = oneminus;
      double specmult = 8.0 * specparm;
      int js = 1 + static_cast<int>(specmult);
      double fs = specmult - static_cast<double>(js - 1);

      double fac000 = (1.0 - fs) * fac00[lay];
      double fac010 = (1.0 - fs) * fac10[lay];
      double fac100 = fs * fac00[lay];
      double fac110 = fs * fac10[lay];
      double fac001 = (1.0 - fs) * fac01[lay];
      double fac011 = (1.0 - fs) * fac11[lay];
      double fac101 = fs * fac01[lay];
      double fac111 = fs * fac11[lay];

      int ind0 = ((jp[lay] - 1) * 5 + (jt[lay] - 1)) * 9 + js;
      int ind1 = (jp[lay] * 5 + (jt1[lay] - 1)) * 9 + js;
      int inds = indself[lay];
      int indf = indfor[lay];
      double tauray = colmol[lay] * b.rayl;

      for (int ig = 1; ig <= 2; ++ig) {
        taug[lay][offset + ig - 1] = speccomb * (
            fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] +
            fac100 * b.ka[(ig - 1) * 585 + (ind0 + 1 - 1)] +
            fac010 * b.ka[(ig - 1) * 585 + (ind0 + 9 - 1)] +
            fac110 * b.ka[(ig - 1) * 585 + (ind0 + 10 - 1)] +
            fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] +
            fac101 * b.ka[(ig - 1) * 585 + (ind1 + 1 - 1)] +
            fac011 * b.ka[(ig - 1) * 585 + (ind1 + 9 - 1)] +
            fac111 * b.ka[(ig - 1) * 585 + (ind1 + 10 - 1)]) +
            colh2o[lay] * (
            selffac[lay] * (b.selfref[(ig - 1) * 10 + (inds - 1)] +
            selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds + 1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)])) +
            forfac[lay] * (b.forref[(ig - 1) * 3 + (indf - 1)] +
            forfrac[lay] * (b.forref[(ig - 1) * 3 + (indf + 1 - 1)] - b.forref[(ig - 1) * 3 + (indf - 1)]))) +
            o2cont;
        if (lay == laysolfr) {
          sfluxzen[offset + ig - 1] = b.sfluxref[(js - 1) * 2 + (ig - 1)] +
              fs * (b.sfluxref[js * 2 + (ig - 1)] - b.sfluxref[(js - 1) * 2 + (ig - 1)]);
        }
        taur[lay][offset + ig - 1] = tauray;
      }
    }

    for (int lay = laytrop; lay < nlayers; ++lay) {
      double o2cont = 4.35e-4 * colo2[lay] / (350.0 * 2.0);
      int ind0 = ((jp[lay] - 13) * 5 + (jt[lay] - 1)) * 1 + 1;
      int ind1 = ((jp[lay] - 12) * 5 + (jt1[lay] - 1)) * 1 + 1;
      double tauray = colmol[lay] * b.rayl;

      for (int ig = 1; ig <= 2; ++ig) {
        taug[lay][offset + ig - 1] = colo2[lay] * o2adj * (
            fac00[lay] * b.kb[(ig - 1) * 235 + (ind0 - 1)] +
            fac10[lay] * b.kb[(ig - 1) * 235 + (ind0 + 1 - 1)] +
            fac01[lay] * b.kb[(ig - 1) * 235 + (ind1 - 1)] +
            fac11[lay] * b.kb[(ig - 1) * 235 + (ind1 + 1 - 1)]) +
            o2cont;
        taur[lay][offset + ig - 1] = tauray;
      }
    }
  }

  // Band 23 (8050-12850 cm-1): ng23 = 10, offset = 56
  {
    const auto& b = tables.b23;
    int offset = 56;
    int laysolfr = laytrop - 1;
    for (int lay = 0; lay < laytrop; ++lay) {
      if (jp[lay] < b.layreffr && jp[lay + 1] >= b.layreffr) laysolfr = std::min(lay + 1, laytrop - 1);
      int ind0 = ((jp[lay] - 1) * 5 + (jt[lay] - 1)) * 1 + 1;
      int ind1 = (jp[lay] * 5 + (jt1[lay] - 1)) * 1 + 1;
      int inds = indself[lay];
      int indf = indfor[lay];

      for (int ig = 1; ig <= 10; ++ig) {
        double tauray = colmol[lay] * b.rayl[ig - 1];
        taug[lay][offset + ig - 1] = colh2o[lay] * (
            b.givfac * (
            fac00[lay] * b.ka[(ig - 1) * 65 + (ind0 - 1)] +
            fac10[lay] * b.ka[(ig - 1) * 65 + (ind0 + 1 - 1)] +
            fac01[lay] * b.ka[(ig - 1) * 65 + (ind1 - 1)] +
            fac11[lay] * b.ka[(ig - 1) * 65 + (ind1 + 1 - 1)]) +
            selffac[lay] * (b.selfref[(ig - 1) * 10 + (inds - 1)] +
            selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds + 1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)])) +
            forfac[lay] * (b.forref[(ig - 1) * 3 + (indf - 1)] +
            forfrac[lay] * (b.forref[(ig - 1) * 3 + (indf + 1 - 1)] - b.forref[(ig - 1) * 3 + (indf - 1)])));
        if (lay == laysolfr) sfluxzen[offset + ig - 1] = b.sfluxref[ig - 1];
        taur[lay][offset + ig - 1] = tauray;
      }
    }

    for (int lay = laytrop; lay < nlayers; ++lay) {
      for (int ig = 1; ig <= 10; ++ig) {
        taug[lay][offset + ig - 1] = 0.0;
        taur[lay][offset + ig - 1] = colmol[lay] * b.rayl[ig - 1];
      }
    }
  }

  // Band 24 (12850-16000 cm-1): ng24 = 8, offset = 66
  {
    const auto& b = tables.b24;
    int offset = 66;
    int laysolfr = laytrop - 1;
    for (int lay = 0; lay < laytrop; ++lay) {
      if (jp[lay] < b.layreffr && jp[lay + 1] >= b.layreffr) laysolfr = std::min(lay + 1, laytrop - 1);
      double speccomb = colh2o[lay] + b.strrat * colo2[lay];
      double specparm = colh2o[lay] / speccomb;
      if (specparm >= oneminus) specparm = oneminus;
      double specmult = 8.0 * specparm;
      int js = 1 + static_cast<int>(specmult);
      double fs = specmult - static_cast<double>(js - 1);

      double fac000 = (1.0 - fs) * fac00[lay];
      double fac010 = (1.0 - fs) * fac10[lay];
      double fac100 = fs * fac00[lay];
      double fac110 = fs * fac10[lay];
      double fac001 = (1.0 - fs) * fac01[lay];
      double fac011 = (1.0 - fs) * fac11[lay];
      double fac101 = fs * fac01[lay];
      double fac111 = fs * fac11[lay];

      int ind0 = ((jp[lay] - 1) * 5 + (jt[lay] - 1)) * 9 + js;
      int ind1 = (jp[lay] * 5 + (jt1[lay] - 1)) * 9 + js;
      int inds = indself[lay];
      int indf = indfor[lay];

      for (int ig = 1; ig <= 8; ++ig) {
        double tauray = colmol[lay] * (b.rayla[(js - 1) * 8 + (ig - 1)] +
            fs * (b.rayla[js * 8 + (ig - 1)] - b.rayla[(js - 1) * 8 + (ig - 1)]));
        taug[lay][offset + ig - 1] = speccomb * (
            fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] +
            fac100 * b.ka[(ig - 1) * 585 + (ind0 + 1 - 1)] +
            fac010 * b.ka[(ig - 1) * 585 + (ind0 + 9 - 1)] +
            fac110 * b.ka[(ig - 1) * 585 + (ind0 + 10 - 1)] +
            fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] +
            fac101 * b.ka[(ig - 1) * 585 + (ind1 + 1 - 1)] +
            fac011 * b.ka[(ig - 1) * 585 + (ind1 + 9 - 1)] +
            fac111 * b.ka[(ig - 1) * 585 + (ind1 + 10 - 1)]) +
            colo3[lay] * b.abso3a[ig - 1] +  // WRF order: ozone before the water vapour continuum
            colh2o[lay] * (
            selffac[lay] * (b.selfref[(ig - 1) * 10 + (inds - 1)] +
            selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds + 1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)])) +
            forfac[lay] * (b.forref[(ig - 1) * 3 + (indf - 1)] +
            forfrac[lay] * (b.forref[(ig - 1) * 3 + (indf + 1 - 1)] - b.forref[(ig - 1) * 3 + (indf - 1)])));
        if (lay == laysolfr) {
          sfluxzen[offset + ig - 1] = b.sfluxref[(js - 1) * 8 + (ig - 1)] +
              fs * (b.sfluxref[js * 8 + (ig - 1)] - b.sfluxref[(js - 1) * 8 + (ig - 1)]);
        }
        taur[lay][offset + ig - 1] = tauray;
      }
    }

    for (int lay = laytrop; lay < nlayers; ++lay) {
      int ind0 = ((jp[lay] - 13) * 5 + (jt[lay] - 1)) * 1 + 1;
      int ind1 = ((jp[lay] - 12) * 5 + (jt1[lay] - 1)) * 1 + 1;

      for (int ig = 1; ig <= 8; ++ig) {
        double tauray = colmol[lay] * b.raylb[ig - 1];
        taug[lay][offset + ig - 1] = colo2[lay] * (
            fac00[lay] * b.kb[(ig - 1) * 235 + (ind0 - 1)] +
            fac10[lay] * b.kb[(ig - 1) * 235 + (ind0 + 1 - 1)] +
            fac01[lay] * b.kb[(ig - 1) * 235 + (ind1 - 1)] +
            fac11[lay] * b.kb[(ig - 1) * 235 + (ind1 + 1 - 1)]) +
            colo3[lay] * b.abso3b[ig - 1];
        taur[lay][offset + ig - 1] = tauray;
      }
    }
  }

  // Band 25 (16000-22650 cm-1): ng25 = 6, offset = 74
  {
    const auto& b = tables.b25;
    int offset = 74;
    int laysolfr = laytrop - 1;
    for (int lay = 0; lay < laytrop; ++lay) {
      if (jp[lay] < b.layreffr && jp[lay + 1] >= b.layreffr) laysolfr = std::min(lay + 1, laytrop - 1);
      int ind0 = ((jp[lay] - 1) * 5 + (jt[lay] - 1)) * 1 + 1;
      int ind1 = (jp[lay] * 5 + (jt1[lay] - 1)) * 1 + 1;

      for (int ig = 1; ig <= 6; ++ig) {
        double tauray = colmol[lay] * b.rayl[ig - 1];
        taug[lay][offset + ig - 1] = colh2o[lay] * (
            fac00[lay] * b.ka[(ig - 1) * 65 + (ind0 - 1)] +
            fac10[lay] * b.ka[(ig - 1) * 65 + (ind0 + 1 - 1)] +
            fac01[lay] * b.ka[(ig - 1) * 65 + (ind1 - 1)] +
            fac11[lay] * b.ka[(ig - 1) * 65 + (ind1 + 1 - 1)]) +
            colo3[lay] * b.abso3a[ig - 1];
        if (lay == laysolfr) sfluxzen[offset + ig - 1] = b.sfluxref[ig - 1];
        taur[lay][offset + ig - 1] = tauray;
      }
    }

    for (int lay = laytrop; lay < nlayers; ++lay) {
      for (int ig = 1; ig <= 6; ++ig) {
        double tauray = colmol[lay] * b.rayl[ig - 1];
        taug[lay][offset + ig - 1] = colo3[lay] * b.abso3b[ig - 1];
        taur[lay][offset + ig - 1] = tauray;
      }
    }
  }

  // Band 26 (22650-29000 cm-1): ng26 = 6, offset = 80
  {
    const auto& b = tables.b26;
    int offset = 80;
    for (int lay = 0; lay < nlayers; ++lay) {
      for (int ig = 1; ig <= 6; ++ig) {
        taug[lay][offset + ig - 1] = 0.0;
        taur[lay][offset + ig - 1] = colmol[lay] * b.rayl[ig - 1];
      }
    }
    for (int ig = 1; ig <= 6; ++ig) {
      sfluxzen[offset + ig - 1] = b.sfluxref[ig - 1];
    }
  }

  // Band 27 (29000-38000 cm-1): ng27 = 8, offset = 86
  {
    const auto& b = tables.b27;
    int offset = 86;
    for (int lay = 0; lay < laytrop; ++lay) {
      int ind0 = ((jp[lay] - 1) * 5 + (jt[lay] - 1)) * 1 + 1;
      int ind1 = (jp[lay] * 5 + (jt1[lay] - 1)) * 1 + 1;

      for (int ig = 1; ig <= 8; ++ig) {
        double tauray = colmol[lay] * b.rayl[ig - 1];
        taug[lay][offset + ig - 1] = colo3[lay] * (
            fac00[lay] * b.ka[(ig - 1) * 65 + (ind0 - 1)] +
            fac10[lay] * b.ka[(ig - 1) * 65 + (ind0 + 1 - 1)] +
            fac01[lay] * b.ka[(ig - 1) * 65 + (ind1 - 1)] +
            fac11[lay] * b.ka[(ig - 1) * 65 + (ind1 + 1 - 1)]);
        taur[lay][offset + ig - 1] = tauray;
      }
    }

    int laysolfr = nlayers - 1;
    for (int lay = laytrop; lay < nlayers; ++lay) {
      if (jp[lay - 1] < b.layreffr && jp[lay] >= b.layreffr) laysolfr = lay;
      int ind0 = ((jp[lay] - 13) * 5 + (jt[lay] - 1)) * 1 + 1;
      int ind1 = ((jp[lay] - 12) * 5 + (jt1[lay] - 1)) * 1 + 1;

      for (int ig = 1; ig <= 8; ++ig) {
        double tauray = colmol[lay] * b.rayl[ig - 1];
        taug[lay][offset + ig - 1] = colo3[lay] * (
            fac00[lay] * b.kb[(ig - 1) * 235 + (ind0 - 1)] +
            fac10[lay] * b.kb[(ig - 1) * 235 + (ind0 + 1 - 1)] +
            fac01[lay] * b.kb[(ig - 1) * 235 + (ind1 - 1)] +
            fac11[lay] * b.kb[(ig - 1) * 235 + (ind1 + 1 - 1)]);
        if (lay == laysolfr) sfluxzen[offset + ig - 1] = b.scalekur * b.sfluxref[ig - 1];
        taur[lay][offset + ig - 1] = tauray;
      }
    }
  }

  // Band 28 (38000-50000 cm-1): ng28 = 6, offset = 94
  {
    const auto& b = tables.b28;
    int offset = 94;
    for (int lay = 0; lay < laytrop; ++lay) {
      double speccomb = colo3[lay] + b.strrat * colo2[lay];
      double specparm = colo3[lay] / speccomb;
      if (specparm >= oneminus) specparm = oneminus;
      double specmult = 8.0 * specparm;
      int js = 1 + static_cast<int>(specmult);
      double fs = specmult - static_cast<double>(js - 1);

      double fac000 = (1.0 - fs) * fac00[lay];
      double fac010 = (1.0 - fs) * fac10[lay];
      double fac100 = fs * fac00[lay];
      double fac110 = fs * fac10[lay];
      double fac001 = (1.0 - fs) * fac01[lay];
      double fac011 = (1.0 - fs) * fac11[lay];
      double fac101 = fs * fac01[lay];
      double fac111 = fs * fac11[lay];

      int ind0 = ((jp[lay] - 1) * 5 + (jt[lay] - 1)) * 9 + js;
      int ind1 = (jp[lay] * 5 + (jt1[lay] - 1)) * 9 + js;
      double tauray = colmol[lay] * b.rayl;

      for (int ig = 1; ig <= 6; ++ig) {
        taug[lay][offset + ig - 1] = speccomb * (
            fac000 * b.ka[(ig - 1) * 585 + (ind0 - 1)] +
            fac100 * b.ka[(ig - 1) * 585 + (ind0 + 1 - 1)] +
            fac010 * b.ka[(ig - 1) * 585 + (ind0 + 9 - 1)] +
            fac110 * b.ka[(ig - 1) * 585 + (ind0 + 10 - 1)] +
            fac001 * b.ka[(ig - 1) * 585 + (ind1 - 1)] +
            fac101 * b.ka[(ig - 1) * 585 + (ind1 + 1 - 1)] +
            fac011 * b.ka[(ig - 1) * 585 + (ind1 + 9 - 1)] +
            fac111 * b.ka[(ig - 1) * 585 + (ind1 + 10 - 1)]);
        taur[lay][offset + ig - 1] = tauray;
      }
    }

    int laysolfr = nlayers - 1;
    for (int lay = laytrop; lay < nlayers; ++lay) {
      if (jp[lay - 1] < b.layreffr && jp[lay] >= b.layreffr) laysolfr = lay;
      double speccomb = colo3[lay] + b.strrat * colo2[lay];
      double specparm = colo3[lay] / speccomb;
      if (specparm >= oneminus) specparm = oneminus;
      double specmult = 4.0 * specparm;
      int js = 1 + static_cast<int>(specmult);
      double fs = specmult - static_cast<double>(js - 1);

      double fac000 = (1.0 - fs) * fac00[lay];
      double fac010 = (1.0 - fs) * fac10[lay];
      double fac100 = fs * fac00[lay];
      double fac110 = fs * fac10[lay];
      double fac001 = (1.0 - fs) * fac01[lay];
      double fac011 = (1.0 - fs) * fac11[lay];
      double fac101 = fs * fac01[lay];
      double fac111 = fs * fac11[lay];

      int ind0 = ((jp[lay] - 13) * 5 + (jt[lay] - 1)) * 5 + js;
      int ind1 = ((jp[lay] - 12) * 5 + (jt1[lay] - 1)) * 5 + js;
      double tauray = colmol[lay] * b.rayl;

      for (int ig = 1; ig <= 6; ++ig) {
        taug[lay][offset + ig - 1] = speccomb * (
            fac000 * b.kb[(ig - 1) * 1175 + (ind0 - 1)] +
            fac100 * b.kb[(ig - 1) * 1175 + (ind0 + 1 - 1)] +
            fac010 * b.kb[(ig - 1) * 1175 + (ind0 + 5 - 1)] +
            fac110 * b.kb[(ig - 1) * 1175 + (ind0 + 6 - 1)] +
            fac001 * b.kb[(ig - 1) * 1175 + (ind1 - 1)] +
            fac101 * b.kb[(ig - 1) * 1175 + (ind1 + 1 - 1)] +
            fac011 * b.kb[(ig - 1) * 1175 + (ind1 + 5 - 1)] +
            fac111 * b.kb[(ig - 1) * 1175 + (ind1 + 6 - 1)]);
        if (lay == laysolfr) {
          sfluxzen[offset + ig - 1] = b.sfluxref[(js - 1) * 6 + (ig - 1)] +
              fs * (b.sfluxref[js * 6 + (ig - 1)] - b.sfluxref[(js - 1) * 6 + (ig - 1)]);
        }
        taur[lay][offset + ig - 1] = tauray;
      }
    }
  }

  // Band 29 (820-2600 cm-1): ng29 = 12, offset = 100
  {
    const auto& b = tables.b29;
    int offset = 100;
    for (int lay = 0; lay < laytrop; ++lay) {
      int ind0 = ((jp[lay] - 1) * 5 + (jt[lay] - 1)) * 1 + 1;
      int ind1 = (jp[lay] * 5 + (jt1[lay] - 1)) * 1 + 1;
      int inds = indself[lay];
      int indf = indfor[lay];
      double tauray = colmol[lay] * b.rayl;

      for (int ig = 1; ig <= 12; ++ig) {
        taug[lay][offset + ig - 1] = colh2o[lay] * (
            fac00[lay] * b.ka[(ig - 1) * 65 + (ind0 - 1)] +
            fac10[lay] * b.ka[(ig - 1) * 65 + (ind0 + 1 - 1)] +
            fac01[lay] * b.ka[(ig - 1) * 65 + (ind1 - 1)] +
            fac11[lay] * b.ka[(ig - 1) * 65 + (ind1 + 1 - 1)] +
            selffac[lay] * (b.selfref[(ig - 1) * 10 + (inds - 1)] +
            selffrac[lay] * (b.selfref[(ig - 1) * 10 + (inds + 1 - 1)] - b.selfref[(ig - 1) * 10 + (inds - 1)])) +
            forfac[lay] * (b.forref[(ig - 1) * 4 + (indf - 1)] +
            forfrac[lay] * (b.forref[(ig - 1) * 4 + (indf + 1 - 1)] - b.forref[(ig - 1) * 4 + (indf - 1)]))) +
            colco2[lay] * b.absco2[ig - 1];
        taur[lay][offset + ig - 1] = tauray;
      }
    }

    int laysolfr = nlayers - 1;
    for (int lay = laytrop; lay < nlayers; ++lay) {
      if (jp[lay - 1] < b.layreffr && jp[lay] >= b.layreffr) laysolfr = lay;
      int ind0 = ((jp[lay] - 13) * 5 + (jt[lay] - 1)) * 1 + 1;
      int ind1 = ((jp[lay] - 12) * 5 + (jt1[lay] - 1)) * 1 + 1;
      double tauray = colmol[lay] * b.rayl;

      for (int ig = 1; ig <= 12; ++ig) {
        taug[lay][offset + ig - 1] = colco2[lay] * (
            fac00[lay] * b.kb[(ig - 1) * 235 + (ind0 - 1)] +
            fac10[lay] * b.kb[(ig - 1) * 235 + (ind0 + 1 - 1)] +
            fac01[lay] * b.kb[(ig - 1) * 235 + (ind1 - 1)] +
            fac11[lay] * b.kb[(ig - 1) * 235 + (ind1 + 1 - 1)]) +
            colh2o[lay] * b.absh2o[ig - 1];
        if (lay == laysolfr) sfluxzen[offset + ig - 1] = b.sfluxref[ig - 1];
        taur[lay][offset + ig - 1] = tauray;
      }
    }
  }
}

}  // namespace cppwrf::phys
