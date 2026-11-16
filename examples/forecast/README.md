# A weather forecast from the latest GFS run

**Colab:** [cppwrf_forecast_colab.ipynb](https://colab.research.google.com/github/SamurayBalgami/CPPWRF/blob/main/examples/forecast/cppwrf_forecast_colab.ipynb)
builds CPPWRF, WPS and the part of WRF that WPS needs on a Colab GPU, then runs `forecast.sh` and shows the maps.

| File | What it does |
|---|---|
| `forecast.sh` | GFS (NOAA open data, latest run) -> WPS -> `cppwrf_real` -> `cppwrf` -> maps, on one Lambert domain (default: Turkey, 220 x 130 points at 12 km, 24 h, WRF's CONUS physics suite) |
| `maps.py` | every output hour: 2 m temperature with sea-level pressure, 3 h precipitation, 10 m wind, total precipitation; `forecast.gif`; `cities.png` / `cities.txt` (temperature, hourly precipitation, wind and pressure at the grid point nearest to each city) |
| `cppwrf_forecast_colab.ipynb` | the whole thing on Colab |

```bash
WPS_DIR=/path/to/WPS-4.6.0 WRF_DIR=/path/to/WRF-4.8.0 GEOG_DIR=/path/to/WPS_GEOG_LOW_RES CPPWRF_BUILD=build \
  LAT=39 LON=35 HOURS=24 examples/forecast/forecast.sh ~/forecast
```

WPS needs only WRF's I/O libraries, not a WRF build: in the WRF source tree, `./configure` (GNU serial) and
`./compile io`, then WPS's `./configure --build-grib2-libs` with `WRF_DIR` set. WPS 4.6.0's metgrid built by gfortran 15
at `-O` writes `met_em` files without fields; build metgrid at `-O0` (the notebook does).

Measured on an RTX 4060 Laptop GPU (GPU FP32): WPS for 3 GFS times 65 s, `cppwrf_real` 4 s, the 6 h forecast 149 s,
the maps 13 s.

For learning and experiments, not an official forecast.
