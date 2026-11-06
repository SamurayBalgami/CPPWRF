# Notices

CPPWRF: a C++20/CUDA port of the WRF-ARW model.
Copyright (C) 2026 Ömer Kaan Yılmaz (Xassemblianist)

This program is free software: you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details. You should have received a copy of the GNU General
Public License along with this program (file LICENSE). If not, see <https://www.gnu.org/licenses/>.

Every source file carries the tag `SPDX-License-Identifier: GPL-3.0-or-later`. The notices below apply in
addition to that licence to the parts of CPPWRF that are derived from other works.

## WRF

CPPWRF is a statement-by-statement port of the Weather Research and Forecasting (WRF) model v4.8.0 (the ARW
dynamical core, `real.exe`, `ideal.exe` and the physics schemes listed in docs/STATUS.md) to C++ and CUDA. The
comments in the sources cite the WRF files and lines that each part follows. Data derived from WRF:

- `data/wrf_schema/wrf-4.8.0.tsv`: the I/O schema extracted from WRF's Registry (tools/wrf_schema).
- `data/rrtmg/*.bin`: the RRTMG tables as WRF holds them after reading `RRTMG_LW_DATA_DBL` / `RRTMG_SW_DATA_DBL`
  (see below for RRTMG's own notice).
- `src/phys/ra_rrtm_lw/rrtm_lw_data.hpp`: the DATA statements of WRF's `phys/module_ra_rrtm.F`
  (tools/rrtm_lw_data.py).
- Constants and coefficient tables in the scheme sources (Kain-Fritsch, WSM6, Thompson, Noah, Noah-MP, MYNN and
  others), taken from the WRF modules they port.

CPPWRF does not ship WRF's table files (`LANDUSE.TBL`, `VEGPARM.TBL`, `SOILPARM.TBL`, `GENPARM.TBL`,
`NoahmpTable.TBL`, the ozone and greenhouse-gas files, `CCN_ACTIVATE.BIN`, Thompson's `.dat` tables): it reads them
from WRF's `run/` directory at run time. The differential tests under `tests/` link the objects of a WRF build that
the user provides; no WRF source file is part of CPPWRF.

The WRF notice, as requested by UCAR (also in licenses/WRF-LICENSE.txt):

> WRF was developed at the National Center for Atmospheric Research (NCAR) which is operated by
> the University Corporation for Atmospheric Research (UCAR). NCAR and UCAR make no proprietary
> claims, either statutory or otherwise, to this version and release of WRF and consider WRF to
> be in the public domain for use by any person or entity for any purpose without any fee or
> charge. UCAR requests that any WRF user include this notice on any partial or full copies of
> WRF. WRF is provided on an "AS IS" basis and any warranties, either express or implied,
> including but not limited to implied warranties of non-infringement, originality,
> merchantability and fitness for a particular purpose, are disclaimed. In no event shall UCAR
> be liable for any damages, whatsoever, whether direct, indirect, consequential or special,
> that arise out of or in connection with the access, use or performance of WRF, including
> infringement actions.

WRF® is a registered trademark of the University Corporation for Atmospheric Research (UCAR). CPPWRF is an
independent project and is not affiliated with or endorsed by NCAR or UCAR.

## Schemes inside WRF with their own terms

**RRTMG and RRTM (AER).** `src/phys/ra_rrtmg_lw`, `src/phys/ra_rrtmg_sw`, `src/phys/ra_rrtmg_common` and
`data/rrtmg` are derived from WRF's `module_ra_rrtmg_lw.F` and `module_ra_rrtmg_sw.F` (RRTMG, by Atmospheric &
Environmental Research, Inc.), which carry this notice:

> Copyright 2002-2008, Atmospheric & Environmental Research, Inc. (AER).
> This software may be used, copied, or redistributed as long as it is
> not sold and this copyright notice is reproduced on each copy made.
> This model is provided as is without any express or implied warranties.
> (http://www.rtweb.aer.com/)

`src/phys/ra_rrtm_lw` is derived from WRF's `module_ra_rrtm.F`, AER's Rapid Radiative Transfer Model (RRTM) as
included in WRF.

**Noah-MP (UCAR).** `src/phys/lsm_noahmp` is derived from Noah-MP v5.1.2 as distributed with WRF v4.8.0
(`phys/noahmp`). Its licence agreement is reproduced in licenses/Noah-MP-LICENSE.txt and requires this notice:

> The source of this material is the Research Applications Laboratory at the National Center for Atmospheric
> Research, a program of the University Corporation for Atmospheric Research (UCAR) pursuant to a Cooperative
> Agreement with the National Science Foundation; ©2007 University Corporation for Atmospheric Research. All
> Rights Reserved.

and, on scholarly works associated with it: "The Noah-MP modeling system was developed at the National Center for
Atmospheric Research (NCAR) with collaborations from university partners. NCAR is sponsored by the United States
National Science Foundation."

## GNU C Library

`src/core/glibc_math.hpp`, `src/core/glibc_math_tables.hpp`, `src/core/libmvec_x86.hpp` and
`src/core/libmvec_tables.hpp` port algorithms and tables of the GNU C Library 2.41 (libm and libmvec), so that the
model computes the same values as WRF's gfortran build on the CPU and on the GPU. The GNU C Library is
Copyright (C) the Free Software Foundation, Inc. and its contributors, licensed under the GNU Lesser General
Public License, version 2.1 or later, whose section 3 allows these parts to be distributed under the GNU
General Public License, as CPPWRF is.

## Libraries used at build or run time

NetCDF-C (Unidata, BSD-3-Clause licence), HDF5 (The HDF Group, BSD-style licence), an OpenMP runtime and, for GPU
builds, the NVIDIA CUDA Toolkit (NVIDIA's licence). They are not part of this repository.
