#include "cases/config_keys.hpp"

#include "cases/model.hpp"

#include <cctype>
#include <cstdio>

namespace cppwrf::cases {

namespace {

// every "section.key" literal of src/ (src/cases/config_keys.cmake writes the list at build time)
const char* const kKeys[] = {
#include "config_keys.inc"
};

bool nest_key(const std::string& k) {
  if (k.rfind("domains.", 0) == 0) return true;
  return k.size() >= 4 && k[0] == 'd' && std::isdigit(static_cast<unsigned char>(k[1])) &&
         std::isdigit(static_cast<unsigned char>(k[2])) && k[3] == '.';
}

}  // namespace

bool known_config_key(const std::string& key, const io::Config& cfg) {
  for (const char* k : kKeys)
    if (key == k) return true;
  if (nest_key(key)) return true;
  if (key.rfind("tracers.", 0) == 0) {  // model.cpp init_tracer: tracers.<name>_init, _value, _xc, _zc, _radius, _amp
    const auto& e = cfg.entries();
    const auto it = e.find("tracers.names");
    if (it == e.end()) return false;
    for (const std::string& n : split_list(it->second))
      for (const char* f : {"init", "value", "xc", "zc", "radius", "amp"})
        if (key == "tracers." + n + "_" + f) return true;
  }
  return false;
}

void check_config_keys(const io::Config& cfg) {
  std::string bad;
  for (const std::string& k : cfg.unaccessed_keys()) {
    if (!known_config_key(k, cfg)) bad += (bad.empty() ? "'" : ", '") + k + "'";
    // an option this run does not read (another scheme's; tools/gcp/matrix_variant.py reports these as UNREAD-KEYS)
    else std::fprintf(stderr, "cppwrf: warning: unaccessed config key '%s'\n", k.c_str());
  }
  if (!bad.empty())
    throw io::ConfigError("unknown configuration key " + bad +
                          " (no CPPWRF option of that name: misspelled, in the wrong section, or not ported)");
}

}  // namespace cppwrf::cases
