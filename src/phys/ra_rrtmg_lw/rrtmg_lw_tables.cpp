// SPDX-License-Identifier: GPL-3.0-or-later
#include "phys/ra_rrtmg_lw/rrtmg_lw_tables.hpp"
#include "phys/ra_rrtmg_lw/rrtmg_lw_tables_f64.hpp"
#include <iterator>
#include <memory>
#include <type_traits>
#include <fstream>
#include <iostream>

namespace cppwrf::phys {
[[maybe_unused]] static void convert_tables(const f64::Band01Tables& s, Band01Tables& d) {
  for (std::size_t i = 0; i < std::size(s.fracrefa); ++i) d.fracrefa[i] = static_cast<std::remove_reference_t<decltype(d.fracrefa[0])>>(s.fracrefa[i]);
  for (std::size_t i = 0; i < std::size(s.fracrefb); ++i) d.fracrefb[i] = static_cast<std::remove_reference_t<decltype(d.fracrefb[0])>>(s.fracrefb[i]);
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.ka_mn2); ++i) d.ka_mn2[i] = static_cast<std::remove_reference_t<decltype(d.ka_mn2[0])>>(s.ka_mn2[i]);
  for (std::size_t i = 0; i < std::size(s.kb_mn2); ++i) d.kb_mn2[i] = static_cast<std::remove_reference_t<decltype(d.kb_mn2[0])>>(s.kb_mn2[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
}

[[maybe_unused]] static void convert_tables(const f64::Band02Tables& s, Band02Tables& d) {
  for (std::size_t i = 0; i < std::size(s.fracrefa); ++i) d.fracrefa[i] = static_cast<std::remove_reference_t<decltype(d.fracrefa[0])>>(s.fracrefa[i]);
  for (std::size_t i = 0; i < std::size(s.fracrefb); ++i) d.fracrefb[i] = static_cast<std::remove_reference_t<decltype(d.fracrefb[0])>>(s.fracrefb[i]);
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
}

[[maybe_unused]] static void convert_tables(const f64::Band03Tables& s, Band03Tables& d) {
  for (std::size_t i = 0; i < std::size(s.fracrefa); ++i) d.fracrefa[i] = static_cast<std::remove_reference_t<decltype(d.fracrefa[0])>>(s.fracrefa[i]);
  for (std::size_t i = 0; i < std::size(s.fracrefb); ++i) d.fracrefb[i] = static_cast<std::remove_reference_t<decltype(d.fracrefb[0])>>(s.fracrefb[i]);
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.ka_mn2o); ++i) d.ka_mn2o[i] = static_cast<std::remove_reference_t<decltype(d.ka_mn2o[0])>>(s.ka_mn2o[i]);
  for (std::size_t i = 0; i < std::size(s.kb_mn2o); ++i) d.kb_mn2o[i] = static_cast<std::remove_reference_t<decltype(d.kb_mn2o[0])>>(s.kb_mn2o[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
}

[[maybe_unused]] static void convert_tables(const f64::Band04Tables& s, Band04Tables& d) {
  for (std::size_t i = 0; i < std::size(s.fracrefa); ++i) d.fracrefa[i] = static_cast<std::remove_reference_t<decltype(d.fracrefa[0])>>(s.fracrefa[i]);
  for (std::size_t i = 0; i < std::size(s.fracrefb); ++i) d.fracrefb[i] = static_cast<std::remove_reference_t<decltype(d.fracrefb[0])>>(s.fracrefb[i]);
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
}

[[maybe_unused]] static void convert_tables(const f64::Band05Tables& s, Band05Tables& d) {
  for (std::size_t i = 0; i < std::size(s.fracrefa); ++i) d.fracrefa[i] = static_cast<std::remove_reference_t<decltype(d.fracrefa[0])>>(s.fracrefa[i]);
  for (std::size_t i = 0; i < std::size(s.fracrefb); ++i) d.fracrefb[i] = static_cast<std::remove_reference_t<decltype(d.fracrefb[0])>>(s.fracrefb[i]);
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.ka_mo3); ++i) d.ka_mo3[i] = static_cast<std::remove_reference_t<decltype(d.ka_mo3[0])>>(s.ka_mo3[i]);
  for (std::size_t i = 0; i < std::size(s.ccl4); ++i) d.ccl4[i] = static_cast<std::remove_reference_t<decltype(d.ccl4[0])>>(s.ccl4[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
}

[[maybe_unused]] static void convert_tables(const f64::Band06Tables& s, Band06Tables& d) {
  for (std::size_t i = 0; i < std::size(s.fracrefa); ++i) d.fracrefa[i] = static_cast<std::remove_reference_t<decltype(d.fracrefa[0])>>(s.fracrefa[i]);
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.ka_mco2); ++i) d.ka_mco2[i] = static_cast<std::remove_reference_t<decltype(d.ka_mco2[0])>>(s.ka_mco2[i]);
  for (std::size_t i = 0; i < std::size(s.cfc11adj); ++i) d.cfc11adj[i] = static_cast<std::remove_reference_t<decltype(d.cfc11adj[0])>>(s.cfc11adj[i]);
  for (std::size_t i = 0; i < std::size(s.cfc12); ++i) d.cfc12[i] = static_cast<std::remove_reference_t<decltype(d.cfc12[0])>>(s.cfc12[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
}

[[maybe_unused]] static void convert_tables(const f64::Band07Tables& s, Band07Tables& d) {
  for (std::size_t i = 0; i < std::size(s.fracrefa); ++i) d.fracrefa[i] = static_cast<std::remove_reference_t<decltype(d.fracrefa[0])>>(s.fracrefa[i]);
  for (std::size_t i = 0; i < std::size(s.fracrefb); ++i) d.fracrefb[i] = static_cast<std::remove_reference_t<decltype(d.fracrefb[0])>>(s.fracrefb[i]);
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.ka_mco2); ++i) d.ka_mco2[i] = static_cast<std::remove_reference_t<decltype(d.ka_mco2[0])>>(s.ka_mco2[i]);
  for (std::size_t i = 0; i < std::size(s.kb_mco2); ++i) d.kb_mco2[i] = static_cast<std::remove_reference_t<decltype(d.kb_mco2[0])>>(s.kb_mco2[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
}

[[maybe_unused]] static void convert_tables(const f64::Band08Tables& s, Band08Tables& d) {
  for (std::size_t i = 0; i < std::size(s.fracrefa); ++i) d.fracrefa[i] = static_cast<std::remove_reference_t<decltype(d.fracrefa[0])>>(s.fracrefa[i]);
  for (std::size_t i = 0; i < std::size(s.fracrefb); ++i) d.fracrefb[i] = static_cast<std::remove_reference_t<decltype(d.fracrefb[0])>>(s.fracrefb[i]);
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.ka_mco2); ++i) d.ka_mco2[i] = static_cast<std::remove_reference_t<decltype(d.ka_mco2[0])>>(s.ka_mco2[i]);
  for (std::size_t i = 0; i < std::size(s.kb_mco2); ++i) d.kb_mco2[i] = static_cast<std::remove_reference_t<decltype(d.kb_mco2[0])>>(s.kb_mco2[i]);
  for (std::size_t i = 0; i < std::size(s.ka_mn2o); ++i) d.ka_mn2o[i] = static_cast<std::remove_reference_t<decltype(d.ka_mn2o[0])>>(s.ka_mn2o[i]);
  for (std::size_t i = 0; i < std::size(s.kb_mn2o); ++i) d.kb_mn2o[i] = static_cast<std::remove_reference_t<decltype(d.kb_mn2o[0])>>(s.kb_mn2o[i]);
  for (std::size_t i = 0; i < std::size(s.ka_mo3); ++i) d.ka_mo3[i] = static_cast<std::remove_reference_t<decltype(d.ka_mo3[0])>>(s.ka_mo3[i]);
  for (std::size_t i = 0; i < std::size(s.cfc12); ++i) d.cfc12[i] = static_cast<std::remove_reference_t<decltype(d.cfc12[0])>>(s.cfc12[i]);
  for (std::size_t i = 0; i < std::size(s.cfc22adj); ++i) d.cfc22adj[i] = static_cast<std::remove_reference_t<decltype(d.cfc22adj[0])>>(s.cfc22adj[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
}

[[maybe_unused]] static void convert_tables(const f64::Band09Tables& s, Band09Tables& d) {
  for (std::size_t i = 0; i < std::size(s.fracrefa); ++i) d.fracrefa[i] = static_cast<std::remove_reference_t<decltype(d.fracrefa[0])>>(s.fracrefa[i]);
  for (std::size_t i = 0; i < std::size(s.fracrefb); ++i) d.fracrefb[i] = static_cast<std::remove_reference_t<decltype(d.fracrefb[0])>>(s.fracrefb[i]);
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.ka_mn2o); ++i) d.ka_mn2o[i] = static_cast<std::remove_reference_t<decltype(d.ka_mn2o[0])>>(s.ka_mn2o[i]);
  for (std::size_t i = 0; i < std::size(s.kb_mn2o); ++i) d.kb_mn2o[i] = static_cast<std::remove_reference_t<decltype(d.kb_mn2o[0])>>(s.kb_mn2o[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
}

[[maybe_unused]] static void convert_tables(const f64::Band10Tables& s, Band10Tables& d) {
  for (std::size_t i = 0; i < std::size(s.fracrefa); ++i) d.fracrefa[i] = static_cast<std::remove_reference_t<decltype(d.fracrefa[0])>>(s.fracrefa[i]);
  for (std::size_t i = 0; i < std::size(s.fracrefb); ++i) d.fracrefb[i] = static_cast<std::remove_reference_t<decltype(d.fracrefb[0])>>(s.fracrefb[i]);
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
}

[[maybe_unused]] static void convert_tables(const f64::Band11Tables& s, Band11Tables& d) {
  for (std::size_t i = 0; i < std::size(s.fracrefa); ++i) d.fracrefa[i] = static_cast<std::remove_reference_t<decltype(d.fracrefa[0])>>(s.fracrefa[i]);
  for (std::size_t i = 0; i < std::size(s.fracrefb); ++i) d.fracrefb[i] = static_cast<std::remove_reference_t<decltype(d.fracrefb[0])>>(s.fracrefb[i]);
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.ka_mo2); ++i) d.ka_mo2[i] = static_cast<std::remove_reference_t<decltype(d.ka_mo2[0])>>(s.ka_mo2[i]);
  for (std::size_t i = 0; i < std::size(s.kb_mo2); ++i) d.kb_mo2[i] = static_cast<std::remove_reference_t<decltype(d.kb_mo2[0])>>(s.kb_mo2[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
}

[[maybe_unused]] static void convert_tables(const f64::Band12Tables& s, Band12Tables& d) {
  for (std::size_t i = 0; i < std::size(s.fracrefa); ++i) d.fracrefa[i] = static_cast<std::remove_reference_t<decltype(d.fracrefa[0])>>(s.fracrefa[i]);
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
}

[[maybe_unused]] static void convert_tables(const f64::Band13Tables& s, Band13Tables& d) {
  for (std::size_t i = 0; i < std::size(s.fracrefa); ++i) d.fracrefa[i] = static_cast<std::remove_reference_t<decltype(d.fracrefa[0])>>(s.fracrefa[i]);
  for (std::size_t i = 0; i < std::size(s.fracrefb); ++i) d.fracrefb[i] = static_cast<std::remove_reference_t<decltype(d.fracrefb[0])>>(s.fracrefb[i]);
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.ka_mco2); ++i) d.ka_mco2[i] = static_cast<std::remove_reference_t<decltype(d.ka_mco2[0])>>(s.ka_mco2[i]);
  for (std::size_t i = 0; i < std::size(s.ka_mco); ++i) d.ka_mco[i] = static_cast<std::remove_reference_t<decltype(d.ka_mco[0])>>(s.ka_mco[i]);
  for (std::size_t i = 0; i < std::size(s.kb_mo3); ++i) d.kb_mo3[i] = static_cast<std::remove_reference_t<decltype(d.kb_mo3[0])>>(s.kb_mo3[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
}

[[maybe_unused]] static void convert_tables(const f64::Band14Tables& s, Band14Tables& d) {
  for (std::size_t i = 0; i < std::size(s.fracrefa); ++i) d.fracrefa[i] = static_cast<std::remove_reference_t<decltype(d.fracrefa[0])>>(s.fracrefa[i]);
  for (std::size_t i = 0; i < std::size(s.fracrefb); ++i) d.fracrefb[i] = static_cast<std::remove_reference_t<decltype(d.fracrefb[0])>>(s.fracrefb[i]);
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
}

[[maybe_unused]] static void convert_tables(const f64::Band15Tables& s, Band15Tables& d) {
  for (std::size_t i = 0; i < std::size(s.fracrefa); ++i) d.fracrefa[i] = static_cast<std::remove_reference_t<decltype(d.fracrefa[0])>>(s.fracrefa[i]);
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.ka_mn2); ++i) d.ka_mn2[i] = static_cast<std::remove_reference_t<decltype(d.ka_mn2[0])>>(s.ka_mn2[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
}

[[maybe_unused]] static void convert_tables(const f64::Band16Tables& s, Band16Tables& d) {
  for (std::size_t i = 0; i < std::size(s.fracrefa); ++i) d.fracrefa[i] = static_cast<std::remove_reference_t<decltype(d.fracrefa[0])>>(s.fracrefa[i]);
  for (std::size_t i = 0; i < std::size(s.fracrefb); ++i) d.fracrefb[i] = static_cast<std::remove_reference_t<decltype(d.fracrefb[0])>>(s.fracrefb[i]);
  for (std::size_t i = 0; i < std::size(s.ka); ++i) d.ka[i] = static_cast<std::remove_reference_t<decltype(d.ka[0])>>(s.ka[i]);
  for (std::size_t i = 0; i < std::size(s.kb); ++i) d.kb[i] = static_cast<std::remove_reference_t<decltype(d.kb[0])>>(s.kb[i]);
  for (std::size_t i = 0; i < std::size(s.selfref); ++i) d.selfref[i] = static_cast<std::remove_reference_t<decltype(d.selfref[0])>>(s.selfref[i]);
  for (std::size_t i = 0; i < std::size(s.forref); ++i) d.forref[i] = static_cast<std::remove_reference_t<decltype(d.forref[0])>>(s.forref[i]);
}

[[maybe_unused]] static void convert_tables(const f64::RrtmgLwTables& s, RrtmgLwTables& d) {
  for (std::size_t i = 0; i < std::size(s.tau_tbl); ++i) d.tau_tbl[i] = static_cast<std::remove_reference_t<decltype(d.tau_tbl[0])>>(s.tau_tbl[i]);
  for (std::size_t i = 0; i < std::size(s.exp_tbl); ++i) d.exp_tbl[i] = static_cast<std::remove_reference_t<decltype(d.exp_tbl[0])>>(s.exp_tbl[i]);
  for (std::size_t i = 0; i < std::size(s.tfn_tbl); ++i) d.tfn_tbl[i] = static_cast<std::remove_reference_t<decltype(d.tfn_tbl[0])>>(s.tfn_tbl[i]);
  for (std::size_t i = 0; i < std::size(s.totplnk); ++i) d.totplnk[i] = static_cast<std::remove_reference_t<decltype(d.totplnk[0])>>(s.totplnk[i]);
  for (std::size_t i = 0; i < std::size(s.totplk16); ++i) d.totplk16[i] = static_cast<std::remove_reference_t<decltype(d.totplk16[0])>>(s.totplk16[i]);
  for (std::size_t i = 0; i < std::size(s.delwave); ++i) d.delwave[i] = static_cast<std::remove_reference_t<decltype(d.delwave[0])>>(s.delwave[i]);
  for (std::size_t i = 0; i < std::size(s.wavenum1); ++i) d.wavenum1[i] = static_cast<std::remove_reference_t<decltype(d.wavenum1[0])>>(s.wavenum1[i]);
  for (std::size_t i = 0; i < std::size(s.wavenum2); ++i) d.wavenum2[i] = static_cast<std::remove_reference_t<decltype(d.wavenum2[0])>>(s.wavenum2[i]);
  for (std::size_t i = 0; i < std::size(s.ngc); ++i) d.ngc[i] = static_cast<std::remove_reference_t<decltype(d.ngc[0])>>(s.ngc[i]);
  for (std::size_t i = 0; i < std::size(s.ngs); ++i) d.ngs[i] = static_cast<std::remove_reference_t<decltype(d.ngs[0])>>(s.ngs[i]);
  for (std::size_t i = 0; i < std::size(s.ngn); ++i) d.ngn[i] = static_cast<std::remove_reference_t<decltype(d.ngn[0])>>(s.ngn[i]);
  for (std::size_t i = 0; i < std::size(s.ngb); ++i) d.ngb[i] = static_cast<std::remove_reference_t<decltype(d.ngb[0])>>(s.ngb[i]);
  for (std::size_t i = 0; i < std::size(s.ngm); ++i) d.ngm[i] = static_cast<std::remove_reference_t<decltype(d.ngm[0])>>(s.ngm[i]);
  for (std::size_t i = 0; i < std::size(s.wt); ++i) d.wt[i] = static_cast<std::remove_reference_t<decltype(d.wt[0])>>(s.wt[i]);
  for (std::size_t i = 0; i < std::size(s.rwgt); ++i) d.rwgt[i] = static_cast<std::remove_reference_t<decltype(d.rwgt[0])>>(s.rwgt[i]);
  d.abscld1 = static_cast<decltype(d.abscld1)>(s.abscld1);
  for (std::size_t i = 0; i < std::size(s.absice0); ++i) d.absice0[i] = static_cast<std::remove_reference_t<decltype(d.absice0[0])>>(s.absice0[i]);
  for (std::size_t i = 0; i < std::size(s.absice1); ++i) d.absice1[i] = static_cast<std::remove_reference_t<decltype(d.absice1[0])>>(s.absice1[i]);
  for (std::size_t i = 0; i < std::size(s.absice2); ++i) d.absice2[i] = static_cast<std::remove_reference_t<decltype(d.absice2[0])>>(s.absice2[i]);
  for (std::size_t i = 0; i < std::size(s.absice3); ++i) d.absice3[i] = static_cast<std::remove_reference_t<decltype(d.absice3[0])>>(s.absice3[i]);
  d.absliq0 = static_cast<decltype(d.absliq0)>(s.absliq0);
  for (std::size_t i = 0; i < std::size(s.absliq1); ++i) d.absliq1[i] = static_cast<std::remove_reference_t<decltype(d.absliq1[0])>>(s.absliq1[i]);
  for (std::size_t i = 0; i < std::size(s.preflog); ++i) d.preflog[i] = static_cast<std::remove_reference_t<decltype(d.preflog[0])>>(s.preflog[i]);
  for (std::size_t i = 0; i < std::size(s.tref); ++i) d.tref[i] = static_cast<std::remove_reference_t<decltype(d.tref[0])>>(s.tref[i]);
  for (std::size_t i = 0; i < std::size(s.chi_mls); ++i) d.chi_mls[i] = static_cast<std::remove_reference_t<decltype(d.chi_mls[0])>>(s.chi_mls[i]);
  convert_tables(s.b01, d.b01);
  convert_tables(s.b02, d.b02);
  convert_tables(s.b03, d.b03);
  convert_tables(s.b04, d.b04);
  convert_tables(s.b05, d.b05);
  convert_tables(s.b06, d.b06);
  convert_tables(s.b07, d.b07);
  convert_tables(s.b08, d.b08);
  convert_tables(s.b09, d.b09);
  convert_tables(s.b10, d.b10);
  convert_tables(s.b11, d.b11);
  convert_tables(s.b12, d.b12);
  convert_tables(s.b13, d.b13);
  convert_tables(s.b14, d.b14);
  convert_tables(s.b15, d.b15);
  convert_tables(s.b16, d.b16);
}


bool load_rrtmg_lw_tables(const std::string& filepath, RrtmgLwTables& tables) {
  std::ifstream file(filepath, std::ios::binary);
  if (file.is_open()) {
    auto f = std::make_unique<f64::RrtmgLwTables>();
    file.read(reinterpret_cast<char*>(f.get()), sizeof(f64::RrtmgLwTables));
    if (file.gcount() == static_cast<std::streamsize>(sizeof(f64::RrtmgLwTables))) {
      convert_tables(*f, tables);
      return true;
    }
  }
  std::cerr << "Error: could not read rrtmg_lw_tables.bin from " << filepath << "\n";
  return false;
}

}  // namespace cppwrf::phys
