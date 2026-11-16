// SPDX-License-Identifier: GPL-3.0-or-later
#include "phys/ra_rrtmg_sw/rrtmg_sw_tables.hpp"
#include "phys/ra_rrtmg_sw/rrtmg_sw_tables_f64.hpp"
#include <iterator>
#include <memory>
#include <type_traits>
#include <fstream>
#include <iostream>

namespace cppwrf::phys {
[[maybe_unused]] static void convert_tables(const f64::SwBand16Tables& s, SwBand16Tables& d) {
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
  for (std::size_t i = 0; i < std::size(s.sfluxref); ++i) d.sfluxref[i] = static_cast<std::remove_reference_t<decltype(d.sfluxref[0])>>(s.sfluxref[i]);
  d.rayl = static_cast<decltype(d.rayl)>(s.rayl);
  d.strrat1 = static_cast<decltype(d.strrat1)>(s.strrat1);
  d.layreffr = static_cast<decltype(d.layreffr)>(s.layreffr);
}

[[maybe_unused]] static void convert_tables(const f64::Band17Tables& s, Band17Tables& d) {
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
  for (std::size_t i = 0; i < std::size(s.sfluxref); ++i) d.sfluxref[i] = static_cast<std::remove_reference_t<decltype(d.sfluxref[0])>>(s.sfluxref[i]);
  d.rayl = static_cast<decltype(d.rayl)>(s.rayl);
  d.strrat = static_cast<decltype(d.strrat)>(s.strrat);
  d.layreffr = static_cast<decltype(d.layreffr)>(s.layreffr);
}

[[maybe_unused]] static void convert_tables(const f64::Band18Tables& s, Band18Tables& d) {
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
  for (std::size_t i = 0; i < std::size(s.sfluxref); ++i) d.sfluxref[i] = static_cast<std::remove_reference_t<decltype(d.sfluxref[0])>>(s.sfluxref[i]);
  d.rayl = static_cast<decltype(d.rayl)>(s.rayl);
  d.strrat = static_cast<decltype(d.strrat)>(s.strrat);
  d.layreffr = static_cast<decltype(d.layreffr)>(s.layreffr);
}

[[maybe_unused]] static void convert_tables(const f64::Band19Tables& s, Band19Tables& d) {
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
  for (std::size_t i = 0; i < std::size(s.sfluxref); ++i) d.sfluxref[i] = static_cast<std::remove_reference_t<decltype(d.sfluxref[0])>>(s.sfluxref[i]);
  d.rayl = static_cast<decltype(d.rayl)>(s.rayl);
  d.strrat = static_cast<decltype(d.strrat)>(s.strrat);
  d.layreffr = static_cast<decltype(d.layreffr)>(s.layreffr);
}

[[maybe_unused]] static void convert_tables(const f64::Band20Tables& s, Band20Tables& d) {
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
  for (std::size_t i = 0; i < std::size(s.sfluxref); ++i) d.sfluxref[i] = static_cast<std::remove_reference_t<decltype(d.sfluxref[0])>>(s.sfluxref[i]);
  for (std::size_t i = 0; i < std::size(s.absch4); ++i) d.absch4[i] = static_cast<std::remove_reference_t<decltype(d.absch4[0])>>(s.absch4[i]);
  d.rayl = static_cast<decltype(d.rayl)>(s.rayl);
  d.layreffr = static_cast<decltype(d.layreffr)>(s.layreffr);
}

[[maybe_unused]] static void convert_tables(const f64::Band21Tables& s, Band21Tables& d) {
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
  for (std::size_t i = 0; i < std::size(s.sfluxref); ++i) d.sfluxref[i] = static_cast<std::remove_reference_t<decltype(d.sfluxref[0])>>(s.sfluxref[i]);
  d.rayl = static_cast<decltype(d.rayl)>(s.rayl);
  d.strrat = static_cast<decltype(d.strrat)>(s.strrat);
  d.layreffr = static_cast<decltype(d.layreffr)>(s.layreffr);
}

[[maybe_unused]] static void convert_tables(const f64::Band22Tables& s, Band22Tables& d) {
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
  for (std::size_t i = 0; i < std::size(s.sfluxref); ++i) d.sfluxref[i] = static_cast<std::remove_reference_t<decltype(d.sfluxref[0])>>(s.sfluxref[i]);
  d.rayl = static_cast<decltype(d.rayl)>(s.rayl);
  d.strrat = static_cast<decltype(d.strrat)>(s.strrat);
  d.layreffr = static_cast<decltype(d.layreffr)>(s.layreffr);
}

[[maybe_unused]] static void convert_tables(const f64::Band23Tables& s, Band23Tables& d) {
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
  for (std::size_t i = 0; i < std::size(s.sfluxref); ++i) d.sfluxref[i] = static_cast<std::remove_reference_t<decltype(d.sfluxref[0])>>(s.sfluxref[i]);
  for (std::size_t i = 0; i < std::size(s.rayl); ++i) d.rayl[i] = static_cast<std::remove_reference_t<decltype(d.rayl[0])>>(s.rayl[i]);
  d.givfac = static_cast<decltype(d.givfac)>(s.givfac);
  d.layreffr = static_cast<decltype(d.layreffr)>(s.layreffr);
}

[[maybe_unused]] static void convert_tables(const f64::Band24Tables& s, Band24Tables& d) {
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
  for (std::size_t i = 0; i < std::size(s.sfluxref); ++i) d.sfluxref[i] = static_cast<std::remove_reference_t<decltype(d.sfluxref[0])>>(s.sfluxref[i]);
  for (std::size_t i = 0; i < std::size(s.abso3a); ++i) d.abso3a[i] = static_cast<std::remove_reference_t<decltype(d.abso3a[0])>>(s.abso3a[i]);
  for (std::size_t i = 0; i < std::size(s.abso3b); ++i) d.abso3b[i] = static_cast<std::remove_reference_t<decltype(d.abso3b[0])>>(s.abso3b[i]);
  for (std::size_t i = 0; i < std::size(s.rayla); ++i) d.rayla[i] = static_cast<std::remove_reference_t<decltype(d.rayla[0])>>(s.rayla[i]);
  for (std::size_t i = 0; i < std::size(s.raylb); ++i) d.raylb[i] = static_cast<std::remove_reference_t<decltype(d.raylb[0])>>(s.raylb[i]);
  d.strrat = static_cast<decltype(d.strrat)>(s.strrat);
  d.layreffr = static_cast<decltype(d.layreffr)>(s.layreffr);
}

[[maybe_unused]] static void convert_tables(const f64::Band25Tables& s, Band25Tables& d) {
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.sfluxref); ++i) d.sfluxref[i] = static_cast<std::remove_reference_t<decltype(d.sfluxref[0])>>(s.sfluxref[i]);
  for (std::size_t i = 0; i < std::size(s.abso3a); ++i) d.abso3a[i] = static_cast<std::remove_reference_t<decltype(d.abso3a[0])>>(s.abso3a[i]);
  for (std::size_t i = 0; i < std::size(s.abso3b); ++i) d.abso3b[i] = static_cast<std::remove_reference_t<decltype(d.abso3b[0])>>(s.abso3b[i]);
  for (std::size_t i = 0; i < std::size(s.rayl); ++i) d.rayl[i] = static_cast<std::remove_reference_t<decltype(d.rayl[0])>>(s.rayl[i]);
  d.layreffr = static_cast<decltype(d.layreffr)>(s.layreffr);
}

[[maybe_unused]] static void convert_tables(const f64::Band26Tables& s, Band26Tables& d) {
  for (std::size_t i = 0; i < std::size(s.sfluxref); ++i) d.sfluxref[i] = static_cast<std::remove_reference_t<decltype(d.sfluxref[0])>>(s.sfluxref[i]);
  for (std::size_t i = 0; i < std::size(s.rayl); ++i) d.rayl[i] = static_cast<std::remove_reference_t<decltype(d.rayl[0])>>(s.rayl[i]);
}

[[maybe_unused]] static void convert_tables(const f64::Band27Tables& s, Band27Tables& d) {
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.sfluxref); ++i) d.sfluxref[i] = static_cast<std::remove_reference_t<decltype(d.sfluxref[0])>>(s.sfluxref[i]);
  for (std::size_t i = 0; i < std::size(s.rayl); ++i) d.rayl[i] = static_cast<std::remove_reference_t<decltype(d.rayl[0])>>(s.rayl[i]);
  d.scalekur = static_cast<decltype(d.scalekur)>(s.scalekur);
  d.layreffr = static_cast<decltype(d.layreffr)>(s.layreffr);
}

[[maybe_unused]] static void convert_tables(const f64::Band28Tables& s, Band28Tables& d) {
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.sfluxref); ++i) d.sfluxref[i] = static_cast<std::remove_reference_t<decltype(d.sfluxref[0])>>(s.sfluxref[i]);
  d.rayl = static_cast<decltype(d.rayl)>(s.rayl);
  d.strrat = static_cast<decltype(d.strrat)>(s.strrat);
  d.layreffr = static_cast<decltype(d.layreffr)>(s.layreffr);
}

[[maybe_unused]] static void convert_tables(const f64::Band29Tables& s, Band29Tables& d) {
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
  for (std::size_t i = 0; i < std::size(s.sfluxref); ++i) d.sfluxref[i] = static_cast<std::remove_reference_t<decltype(d.sfluxref[0])>>(s.sfluxref[i]);
  for (std::size_t i = 0; i < std::size(s.absh2o); ++i) d.absh2o[i] = static_cast<std::remove_reference_t<decltype(d.absh2o[0])>>(s.absh2o[i]);
  for (std::size_t i = 0; i < std::size(s.absco2); ++i) d.absco2[i] = static_cast<std::remove_reference_t<decltype(d.absco2[0])>>(s.absco2[i]);
  d.rayl = static_cast<decltype(d.rayl)>(s.rayl);
  d.layreffr = static_cast<decltype(d.layreffr)>(s.layreffr);
}

[[maybe_unused]] static void convert_tables(const f64::RrtmgSwTables& s, RrtmgSwTables& d) {
  for (std::size_t i = 0; i < std::size(s.exp_tbl); ++i) d.exp_tbl[i] = static_cast<std::remove_reference_t<decltype(d.exp_tbl[0])>>(s.exp_tbl[i]);
  for (std::size_t i = 0; i < std::size(s.wavenum1); ++i) d.wavenum1[i] = static_cast<std::remove_reference_t<decltype(d.wavenum1[0])>>(s.wavenum1[i]);
  for (std::size_t i = 0; i < std::size(s.wavenum2); ++i) d.wavenum2[i] = static_cast<std::remove_reference_t<decltype(d.wavenum2[0])>>(s.wavenum2[i]);
  for (std::size_t i = 0; i < std::size(s.delwave); ++i) d.delwave[i] = static_cast<std::remove_reference_t<decltype(d.delwave[0])>>(s.delwave[i]);
  for (std::size_t i = 0; i < std::size(s.ngc); ++i) d.ngc[i] = static_cast<std::remove_reference_t<decltype(d.ngc[0])>>(s.ngc[i]);
  for (std::size_t i = 0; i < std::size(s.ngs); ++i) d.ngs[i] = static_cast<std::remove_reference_t<decltype(d.ngs[0])>>(s.ngs[i]);
  for (std::size_t i = 0; i < std::size(s.ngn); ++i) d.ngn[i] = static_cast<std::remove_reference_t<decltype(d.ngn[0])>>(s.ngn[i]);
  for (std::size_t i = 0; i < std::size(s.ngb); ++i) d.ngb[i] = static_cast<std::remove_reference_t<decltype(d.ngb[0])>>(s.ngb[i]);
  for (std::size_t i = 0; i < std::size(s.ngm); ++i) d.ngm[i] = static_cast<std::remove_reference_t<decltype(d.ngm[0])>>(s.ngm[i]);
  for (std::size_t i = 0; i < std::size(s.wt); ++i) d.wt[i] = static_cast<std::remove_reference_t<decltype(d.wt[0])>>(s.wt[i]);
  for (std::size_t i = 0; i < std::size(s.rwgt); ++i) d.rwgt[i] = static_cast<std::remove_reference_t<decltype(d.rwgt[0])>>(s.rwgt[i]);
  for (std::size_t i = 0; i < std::size(s.extliq1); ++i) d.extliq1[i] = static_cast<std::remove_reference_t<decltype(d.extliq1[0])>>(s.extliq1[i]);
  for (std::size_t i = 0; i < std::size(s.ssaliq1); ++i) d.ssaliq1[i] = static_cast<std::remove_reference_t<decltype(d.ssaliq1[0])>>(s.ssaliq1[i]);
  for (std::size_t i = 0; i < std::size(s.asyliq1); ++i) d.asyliq1[i] = static_cast<std::remove_reference_t<decltype(d.asyliq1[0])>>(s.asyliq1[i]);
  for (std::size_t i = 0; i < std::size(s.extice2); ++i) d.extice2[i] = static_cast<std::remove_reference_t<decltype(d.extice2[0])>>(s.extice2[i]);
  for (std::size_t i = 0; i < std::size(s.ssaice2); ++i) d.ssaice2[i] = static_cast<std::remove_reference_t<decltype(d.ssaice2[0])>>(s.ssaice2[i]);
  for (std::size_t i = 0; i < std::size(s.asyice2); ++i) d.asyice2[i] = static_cast<std::remove_reference_t<decltype(d.asyice2[0])>>(s.asyice2[i]);
  for (std::size_t i = 0; i < std::size(s.extice3); ++i) d.extice3[i] = static_cast<std::remove_reference_t<decltype(d.extice3[0])>>(s.extice3[i]);
  for (std::size_t i = 0; i < std::size(s.ssaice3); ++i) d.ssaice3[i] = static_cast<std::remove_reference_t<decltype(d.ssaice3[0])>>(s.ssaice3[i]);
  for (std::size_t i = 0; i < std::size(s.asyice3); ++i) d.asyice3[i] = static_cast<std::remove_reference_t<decltype(d.asyice3[0])>>(s.asyice3[i]);
  for (std::size_t i = 0; i < std::size(s.fdlice3); ++i) d.fdlice3[i] = static_cast<std::remove_reference_t<decltype(d.fdlice3[0])>>(s.fdlice3[i]);
  for (std::size_t i = 0; i < std::size(s.abari); ++i) d.abari[i] = static_cast<std::remove_reference_t<decltype(d.abari[0])>>(s.abari[i]);
  for (std::size_t i = 0; i < std::size(s.bbari); ++i) d.bbari[i] = static_cast<std::remove_reference_t<decltype(d.bbari[0])>>(s.bbari[i]);
  for (std::size_t i = 0; i < std::size(s.cbari); ++i) d.cbari[i] = static_cast<std::remove_reference_t<decltype(d.cbari[0])>>(s.cbari[i]);
  for (std::size_t i = 0; i < std::size(s.dbari); ++i) d.dbari[i] = static_cast<std::remove_reference_t<decltype(d.dbari[0])>>(s.dbari[i]);
  for (std::size_t i = 0; i < std::size(s.ebari); ++i) d.ebari[i] = static_cast<std::remove_reference_t<decltype(d.ebari[0])>>(s.ebari[i]);
  for (std::size_t i = 0; i < std::size(s.fbari); ++i) d.fbari[i] = static_cast<std::remove_reference_t<decltype(d.fbari[0])>>(s.fbari[i]);
  for (std::size_t i = 0; i < std::size(s.preflog); ++i) d.preflog[i] = static_cast<std::remove_reference_t<decltype(d.preflog[0])>>(s.preflog[i]);
  for (std::size_t i = 0; i < std::size(s.tref); ++i) d.tref[i] = static_cast<std::remove_reference_t<decltype(d.tref[0])>>(s.tref[i]);
  for (std::size_t i = 0; i < std::size(s.rsrtaua); ++i) d.rsrtaua[i] = static_cast<std::remove_reference_t<decltype(d.rsrtaua[0])>>(s.rsrtaua[i]);
  for (std::size_t i = 0; i < std::size(s.rsrpiza); ++i) d.rsrpiza[i] = static_cast<std::remove_reference_t<decltype(d.rsrpiza[0])>>(s.rsrpiza[i]);
  for (std::size_t i = 0; i < std::size(s.rsrasya); ++i) d.rsrasya[i] = static_cast<std::remove_reference_t<decltype(d.rsrasya[0])>>(s.rsrasya[i]);
  convert_tables(s.b16, d.b16);
  convert_tables(s.b17, d.b17);
  convert_tables(s.b18, d.b18);
  convert_tables(s.b19, d.b19);
  convert_tables(s.b20, d.b20);
  convert_tables(s.b21, d.b21);
  convert_tables(s.b22, d.b22);
  convert_tables(s.b23, d.b23);
  convert_tables(s.b24, d.b24);
  convert_tables(s.b25, d.b25);
  convert_tables(s.b26, d.b26);
  convert_tables(s.b27, d.b27);
  convert_tables(s.b28, d.b28);
  convert_tables(s.b29, d.b29);
}


bool load_rrtmg_sw_tables(const std::string& filepath, RrtmgSwTables& tables) {
  std::ifstream file(filepath, std::ios::binary);
  if (file.is_open()) {
    auto f = std::make_unique<f64::RrtmgSwTables>();
    file.read(reinterpret_cast<char*>(f.get()), sizeof(f64::RrtmgSwTables));
    if (file.gcount() == static_cast<std::streamsize>(sizeof(f64::RrtmgSwTables))) {
      convert_tables(*f, tables);
      return true;
    }
  }
  std::cerr << "Error: could not read rrtmg_sw_tables.bin from " << filepath << "\n";
  return false;
}

}  // namespace cppwrf::phys
