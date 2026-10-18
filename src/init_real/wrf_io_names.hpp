#pragma once
// The WRF I/O layer (src/wrf_io) under the names real.exe's code uses.
#include "wrf_io/namelist.hpp"
#include "wrf_io/nc_out.hpp"
#include "wrf_io/wrf_schema.hpp"

namespace cppwrf::init_real {
using wrf_io::Namelist;
using wrf_io::NcOut;
using wrf_io::OutVar;
using wrf_io::SchemaMember;
using wrf_io::SchemaVar;
using wrf_io::WrfConfig;
using wrf_io::WrfSchema;
}  // namespace cppwrf::init_real
