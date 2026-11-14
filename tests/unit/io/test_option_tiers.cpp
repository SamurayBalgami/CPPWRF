// SPDX-License-Identifier: GPL-3.0-or-later
// The verification tiers (src/io/option_tiers.hpp, data/option_tiers.txt) on dummy rows: the table parser, value
// matching, the gate (rejected without cppwrf.allow_unverified, accepted with it) and the CPPWRF_UNVERIFIED_OPTIONS
// text. The whole path through a model run and its wrfout is tests/tools/test_option_tiers.py.
#include "io/option_tiers.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace cppwrf::io;

#define CHECK(cond)                                                                        \
  do {                                                                                     \
    if (!(cond)) {                                                                         \
      std::cerr << "FAILED: " #cond " at " << __FILE__ << ":" << __LINE__ << "\n";         \
      std::exit(1);                                                                        \
    }                                                                                      \
  } while (false)

template <typename F>
static std::string error_of(F f) {
  try {
    f();
  } catch (const ConfigError& e) {
    return e.what();
  }
  return "";
}

static const char* kDummy =
    "# dummy rows\n"
    "phys.mp_physics   physics.mp_physics   4,20:23   column_verified  WSM5 (tests/phys/mp_wsm5)\n"
    "phys.sst_skin     physics.sst_skin     true      unverified       no column driver\n"
    "dyn.tiny_opt      -                    *         unverified\n"
    "  \n";

int main() {
  const std::vector<TierRow> rows = parse_option_tiers(kDummy, "dummy");
  CHECK(rows.size() == 3);
  CHECK(rows[0].config_key == "phys.mp_physics" && rows[0].wrf_group == "physics" && rows[0].wrf_key == "mp_physics");
  CHECK(rows[0].tier == Tier::column_verified && rows[0].note == "WSM5 (tests/phys/mp_wsm5)");
  CHECK(rows[2].wrf_key.empty() && rows[2].tier == Tier::unverified);

  // the built-in table parses; each of its rows names its scheme and its column driver or the reason for none
  for (const TierRow& r : option_tiers()) CHECK(!r.config_key.empty() && !r.note.empty());

  // malformed rows
  CHECK(error_of([] { parse_option_tiers("phys.mp_physics physics.mp_physics 4", "t"); }).find("t:1") == 0);
  CHECK(!error_of([] { parse_option_tiers("phys.mp_physics physics.mp_physics 4 verified", "t"); }).empty());
  CHECK(!error_of([] { parse_option_tiers("phys.mp_physics mp_physics 4 unverified", "t"); }).empty());
  CHECK(!error_of([] { parse_option_tiers("phys.mp_physics physics.mp_physics 5:4 unverified", "t"); }).empty());
  CHECK(!error_of([] { parse_option_tiers("mp_physics - 4 unverified", "t"); }).empty());

  // values
  CHECK(tier_value_matches({"4"}, "4") && tier_value_matches({"4"}, " 4.0 ") && !tier_value_matches({"4"}, "40"));
  CHECK(tier_value_matches({"20:23"}, "21") && !tier_value_matches({"20:23"}, "24"));
  CHECK(tier_value_matches({"20:23"}, "20") && tier_value_matches({"20:23"}, "23") && !tier_value_matches({"20:23"}, "19"));
  CHECK(tier_value_matches({"true"}, ".TRUE.") && tier_value_matches({"true"}, "1") && !tier_value_matches({"true"}, "false"));
  CHECK(tier_value_matches({"*"}, "anything") && tier_value_matches({"abc"}, "'ABC'") && !tier_value_matches({"abc"}, "abd"));

  // a configuration with verified options only: no hits, runs without the switch, no attribute
  const Config ok = Config::from_string("[phys]\nmp_physics = 6\nsst_skin = false\n");
  CHECK(check_option_tiers(ok, "", rows).empty());
  CHECK(tier_attribute(config_tier_hits(ok, rows)).empty());

  // unverified options: rejected without the switch, the message naming option, value and tier
  const Config bad = Config::from_string("[phys]\nmp_physics = 4\nsst_skin = true\n[dyn]\ntiny_opt = 3\n");
  const std::string msg = error_of([&] { check_option_tiers(bad, "d02: ", rows); });
  std::cerr << "rejected: " << msg << "\n";
  CHECK(msg.find("d02: ") == 0);
  CHECK(msg.find("phys.mp_physics = 4 is column_verified") != std::string::npos);
  CHECK(msg.find("phys.sst_skin = true is unverified") != std::string::npos);
  CHECK(msg.find("dyn.tiny_opt = 3 is unverified") != std::string::npos);
  CHECK(msg.find("cppwrf.allow_unverified = true") != std::string::npos);
  CHECK(!error_of([&] {
    check_option_tiers(Config::from_string("[phys]\nmp_physics = 4\n[cppwrf]\nallow_unverified = false\n"), "", rows);
  }).empty());

  // accepted with the switch (a warning each), the switch read (no "unaccessed" key), the attribute text
  const Config allowed = Config::from_string(
      "[phys]\nmp_physics = 22\nsst_skin = true\n[dyn]\ntiny_opt = 3\n[cppwrf]\nallow_unverified = true\n");
  const std::vector<TierHit> hits = check_option_tiers(allowed, "", rows);
  CHECK(hits.size() == 3);
  bool unread = false;
  for (const std::string& k : allowed.unaccessed_keys()) unread = unread || k == "cppwrf.allow_unverified";
  CHECK(!unread);
  const std::string att = tier_attribute(hits);
  std::cerr << "attribute: " << att << "\n";
  CHECK(att == "mp_physics=22:column_verified;sst_skin=true:unverified;dyn.tiny_opt=3:unverified");
  // the same option twice (two domains' hits) is listed once
  std::vector<TierHit> twice = hits;
  twice.insert(twice.end(), hits.begin(), hits.end());
  CHECK(tier_attribute(twice) == att);

  // namelist values (real.exe and the translator's view)
  const auto nh = namelist_tier_hits(
      [](const std::string& g, const std::string& k) { return g == "physics" && k == "mp_physics" ? "4" : ""; }, rows);
  CHECK(nh.size() == 1 && nh[0].key == "&physics.mp_physics" && nh[0].option == "mp_physics");
  CHECK(!error_of([&] { enforce_option_tiers(nh, false, "real: ", "&cppwrf allow_unverified = .true."); }).empty());
  CHECK(error_of([&] { enforce_option_tiers(nh, true, "real: ", "&cppwrf allow_unverified = .true."); }).empty());

  std::cout << "test_option_tiers: OK\n";
  return 0;
}
