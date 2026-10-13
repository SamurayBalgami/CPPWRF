# Notices

CPPWRF: a C++20/CUDA re-implementation of the WRF-ARW model.
Copyright (C) 2026 Ömer Kaan Yılmaz (Xassemblianist)

This program is free software: you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details. You should have received a copy of the GNU General
Public License along with this program (file LICENSE). If not, see <https://www.gnu.org/licenses/>.

## WRF

CPPWRF re-implements the numerical formulation of the Weather Research and Forecasting (WRF)
model, and the differential tests under `tests/` contain unmodified copies of WRF Fortran
source files. The WRF notice, as requested by UCAR:

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

WRF® is a registered trademark of UCAR. CPPWRF is an independent project and is not affiliated
with or endorsed by NCAR or UCAR.

Third-party schemes that carry their own license terms inside WRF (for example RRTMG, by AER)
keep those terms; their notices are preserved in the copied source files.
