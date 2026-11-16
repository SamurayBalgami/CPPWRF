#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Maps and city forecasts of a forecast.sh run, with the cppwrf package (docs/PYTHON.md).

    maps.py RUN_DIR OUT_DIR [--cities "Name:lat:lon,..."] [--utc-offset 3]

OUT_DIR gets map_HH.png for every output hour (2 m temperature with sea-level pressure, 3 h precipitation, 10 m wind,
total precipitation), forecast.gif (all hours), cities.png and cities.txt (the grid point nearest to each city that is
inside the domain). Coastlines and borders when cartopy is installed."""
import argparse
import datetime
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "python"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import cppwrf  # noqa: E402
from cppwrf import diagnostics as dg  # noqa: E402

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
except ImportError:  # maps without coastlines
    ccrs = None

CITIES = ("İstanbul:41.01:28.98,Ankara:39.93:32.86,İzmir:38.42:27.14,Bursa:40.19:29.06,Antalya:36.90:30.70,"
          "Adana:37.00:35.32,Konya:37.87:32.48,Samsun:41.29:36.33,Trabzon:41.00:39.72,Erzurum:39.90:41.27,"
          "Van:38.49:43.38,Diyarbakır:37.91:40.24,Gaziantep:37.07:37.38,Kayseri:38.72:35.48,Eskişehir:39.78:30.52,"
          "Edirne:41.68:26.56")


def parse_cities(text):
    out = []
    for item in text.split(","):
        name, la, lo = item.rsplit(":", 2)
        out.append((name.strip(), float(la), float(lo)))
    return out


def nearest(lat, lon, la, lo):
    """(j, i) of the grid point nearest to (la, lo), or None when it is on the domain's edge (outside)."""
    d = (lat - la) ** 2 + ((lon - lo) * np.cos(np.radians(la))) ** 2
    j, i = np.unravel_index(np.argmin(d), d.shape)
    ny, nx = lat.shape
    if j in (0, ny - 1) or i in (0, nx - 1):
        return None
    return j, i


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("out")
    ap.add_argument("--cities", default=CITIES)
    ap.add_argument("--utc-offset", type=float, default=3.0, help="local time = UTC + this (default 3, Turkey)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    w = cppwrf.open_wrfout(a.run, domain=1, xarray=False)
    lat, lon = np.asarray(w.lat), np.asarray(w.lon)
    times = list(w.times)
    local = [t + datetime.timedelta(hours=a.utc_offset) for t in times]
    tz = f"UTC{a.utc_offset:+g}"
    t2 = np.asarray(dg.t2_celsius(w))
    slp = np.asarray(dg.slp(w))
    spd = np.asarray(dg.wind10(w)[0])
    total = np.asarray(dg.precipitation(w))
    rain1 = np.asarray(dg.precip_accum(w, hours=1))
    rain3 = np.asarray(dg.precip_accum(w, hours=3))
    u10, v10 = np.asarray(w["U10"]), np.asarray(w["V10"])

    # maps
    extent = [float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())]
    t_levels = np.arange(np.floor(np.nanmin(t2) / 2) * 2, np.ceil(np.nanmax(t2) / 2) * 2 + 2, 2)
    rain_levels = [0.1, 0.5, 1, 2, 5, 10, 20, 30, 50, 75, 100, 150]
    wind_levels = np.arange(0, max(20, float(np.nanmax(spd)) + 2), 2)
    tot_levels = [0.1, 1, 2, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300]
    q = max(1, lat.shape[1] // 25)
    files = []
    for k, t in enumerate(times):
        kw = {"subplot_kw": {"projection": ccrs.PlateCarree()}} if ccrs else {}
        fig, axs = plt.subplots(2, 2, figsize=(15, 10), **kw)
        panels = [
            (t2[k], "RdYlBu_r", t_levels, "2 m sıcaklık (°C), deniz seviyesi basıncı (hPa)", "°C"),
            (rain3[k] if k >= 3 else rain1[k], "YlGnBu", rain_levels,
             "son 3 saatte yağış (mm)" if k >= 3 else "son 1 saatte yağış (mm)", "mm"),
            (spd[k], "viridis", wind_levels, "10 m rüzgâr (m/s)", "m/s"),
            (total[k], "YlGnBu", tot_levels, "başlangıçtan beri toplam yağış (mm)", "mm"),
        ]
        for ax, (field, cmap, levels, title, unit) in zip(axs.flat, panels):
            data = np.where(np.isnan(field), 0.0, field)
            extend = "both" if unit == "°C" else "max"
            tr = {"transform": ccrs.PlateCarree()} if ccrs else {}
            c = ax.contourf(lon, lat, data, levels=levels, cmap=cmap, extend=extend, **tr)
            fig.colorbar(c, ax=ax, shrink=0.8, label=unit)
            if unit == "°C":
                cs = ax.contour(lon, lat, slp[k], levels=range(940, 1060, 4), colors="k", linewidths=0.6, **tr)
                ax.clabel(cs, fmt="%d", fontsize=7)
            if unit == "m/s":
                ax.quiver(lon[::q, ::q], lat[::q, ::q], u10[k, ::q, ::q], v10[k, ::q, ::q], scale=400, width=0.002,
                          **tr)
            if ccrs:
                ax.set_extent(extent, crs=ccrs.PlateCarree())
                ax.add_feature(cfeature.COASTLINE, linewidth=0.7)
                ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle=":")
            ax.set_title(title, fontsize=10)
        fig.suptitle(f"CPPWRF tahmini: {local[k]:%d.%m.%Y %H:%M} ({tz}), başlangıç +{k} saat  "
                     f"[GFS {times[0]:%Y-%m-%d %H} UTC]", fontsize=13)
        fig.tight_layout()
        f = os.path.join(a.out, f"map_{k:02d}.png")
        fig.savefig(f, dpi=80)
        plt.close(fig)
        files.append(f)
    try:
        from PIL import Image
        frames = [Image.open(f).convert("P", palette=Image.ADAPTIVE) for f in files]
        frames[0].save(os.path.join(a.out, "forecast.gif"), save_all=True, append_images=frames[1:], duration=600,
                       loop=0)
    except ImportError:
        pass

    # city forecasts
    cities = []
    for name, la, lo in parse_cities(a.cities):
        p = nearest(lat, lon, la, lo)
        if p:
            cities.append((name, p))
    if not cities:
        print(f"{a.out}: {len(files)} maps (no city inside the domain)")
        return
    hours = np.arange(len(times))
    lines = [f"CPPWRF tahmini, GFS {times[0]:%Y-%m-%d %H} UTC; saatler {tz}",
             "şehir         saat              sıcaklık  yağış(1 sa)  rüzgâr  basınç"]
    for name, (j, i) in cities:
        for k in range(len(times)):
            r = rain1[k, j, i]
            lines.append(f"{name:<13} {local[k]:%d.%m %H:%M}      {t2[k, j, i]:6.1f} °C  "
                         f"{(0.0 if np.isnan(r) else r):6.1f} mm   {spd[k, j, i]:5.1f} m/s  {slp[k, j, i]:7.1f} hPa")
        lines.append("")
    with open(os.path.join(a.out, "cities.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    n = len(cities)
    cols = 4
    rows = (n + cols - 1) // cols
    fig, axs = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.0 * rows), squeeze=False)
    labels = [t.strftime("%H") if t.hour % 6 else t.strftime("%d.%m\n%H") for t in local]
    for ax, (name, (j, i)) in zip(axs.flat, cities):
        r = np.nan_to_num(rain1[:, j, i])
        ax2 = ax.twinx()
        ax2.bar(hours, r, color="tab:blue", alpha=0.45, width=0.8)
        ax2.set_ylim(0, max(5.0, float(r.max()) * 1.3))
        ax2.set_ylabel("yağış mm/sa", fontsize=7, color="tab:blue")
        ax.plot(hours, t2[:, j, i], color="tab:red", lw=2)
        ax.set_ylabel("°C", fontsize=8, color="tab:red")
        ax.set_zorder(ax2.get_zorder() + 1)
        ax.patch.set_visible(False)
        ax.set_title(f"{name}: {np.nanmin(t2[:, j, i]):.0f}…{np.nanmax(t2[:, j, i]):.0f} °C, "
                     f"{r.sum():.1f} mm, rüzgâr ≤ {spd[:, j, i].max():.0f} m/s", fontsize=9)
        step = max(1, len(hours) // 8)
        ax.set_xticks(hours[::step])
        ax.set_xticklabels(labels[::step], fontsize=7)
        ax.grid(alpha=0.3)
    for ax in list(axs.flat)[n:]:
        ax.axis("off")
    fig.suptitle(f"Şehir tahminleri (CPPWRF, GFS {times[0]:%Y-%m-%d %H} UTC, saatler {tz})", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "cities.png"), dpi=90)
    print(f"{a.out}: {len(files)} maps, forecast.gif, cities.png and cities.txt ({n} cities)")


if __name__ == "__main__":
    main()
