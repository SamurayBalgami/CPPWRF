// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
// The configuration keys CPPWRF knows: a key of a configuration file or --set that no code reads (a misspelled or
// unsupported option) is an error, not a warning with the default (DomainRunner::start).
#include <string>
#include <vector>

#include "io/config.hpp"

namespace cppwrf::cases {

// A key some code path reads: a "section.key" literal of the sources (the build's generated list), a nest key
// ([domains], [d0N]: checked by parse_nest_specs when the run is nested, ignored by a single-domain run as WRF
// ignores the columns past max_dom) or a passive tracer's tracers.<name>_<field> of a name in tracers.names.
bool known_config_key(const std::string& key, const io::Config& cfg);

// The keys of `cfg` that were not read: an unknown one throws ConfigError (naming them all), a known one (an option the
// run does not use) is a warning, "cppwrf: warning: unaccessed config key '<key>'".
void check_config_keys(const io::Config& cfg);

}  // namespace cppwrf::cases
