#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Animated forecast maps of a forecast.sh run: a full-frame dark map, the field in smooth colours (wind, temperature or
precipitation), particles that move with the 10 m wind, city names with their values, the time and a legend.

    animate.py RUN_DIR OUT_DIR [--layers wind,temp,rain] [--cities "Name:lat:lon,..."] [--utc-offset 3]
               [--steps 6] [--fps 12] [--particles 2500] [--width 1280]

OUT_DIR gets anim_<layer>.mp4 (anim_<layer>.gif without ffmpeg) and anim_<layer>.png (the last frame) per layer.
Needs numpy, scipy, matplotlib, cartopy and netCDF4; ffmpeg from PATH or the imageio-ffmpeg package."""
import argparse
import datetime
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "python"))
sys.path.insert(0, HERE)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import animation  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, Normalize  # noqa: E402
from scipy.interpolate import LinearNDInterpolator  # noqa: E402
from scipy.spatial import Delaunay  # noqa: E402

import cppwrf  # noqa: E402
from cppwrf import diagnostics as dg  # noqa: E402
from maps import CITIES, parse_cities  # noqa: E402

DAYS = ["Pzt", "Sal", "Çar", "Per", "Cum", "Cmt", "Paz"]
MONTHS = ["Oca", "Şub", "Mar", "Nis", "May", "Haz", "Tem", "Ağu", "Eyl", "Eki", "Kas", "Ara"]
BG = "#16171d"

# colour scales: (value, colour) stops
WIND = [(0, "#4d5aa0"), (5, "#4a68a8"), (10, "#4a80aa"), (15, "#4b96a2"), (20, "#4fa184"), (30, "#5aac5c"),
        (40, "#a8a84a"), (50, "#c88f3c"), (60, "#c0603e"), (75, "#a8456e"), (90, "#7a4fa0"), (110, "#d9d9f0")]
TEMP = [(-30, "#d8c8f0"), (-20, "#9a6cc8"), (-12, "#5a4fb8"), (-6, "#3a6ec8"), (0, "#4aa3d8"), (5, "#5bc0b0"),
        (10, "#7cc36a"), (15, "#c8d050"), (20, "#f0c040"), (25, "#f09030"), (30, "#e05a30"), (35, "#c02a3a"),
        (45, "#701030")]
RAIN = [(0, "#2a3550"), (0.2, "#3a5a8a"), (0.5, "#3c78b8"), (1, "#3aa0d0"), (2, "#40c0a0"), (4, "#80d050"),
        (7, "#e0d040"), (10, "#f09030"), (20, "#e04040"), (40, "#c030a0")]

LAYERS = {  # title, unit, stops, value at a city: format
    "wind": ("Rüzgâr (10 m)", "km/sa", WIND, "{:.0f}"),
    "temp": ("Sıcaklık (2 m)", "°C", TEMP, "{:.0f}°"),
    "rain": ("Yağış", "mm/sa", RAIN, "{:.1f}"),
}


class StopNorm(Normalize):
    """The stops evenly spaced on the colour scale, linear between two stops (small rain rates stay visible)."""

    def __init__(self, values):
        super().__init__(0.0, 1.0)
        self.values = np.asarray(values, dtype=float)
        self.pos = np.linspace(0.0, 1.0, len(values))

    def __call__(self, value, clip=None):
        v = np.ma.asarray(value)
        return np.ma.masked_array(np.interp(v.filled(np.nan), self.values, self.pos), mask=np.ma.getmask(v))


def colormap(stops):
    pos = np.linspace(0.0, 1.0, len(stops))
    return (LinearSegmentedColormap.from_list("c", [(p, c) for p, (_, c) in zip(pos, stops)], N=512),
            StopNorm([v for v, _ in stops]))


def when(t):
    return f"{DAYS[t.weekday()]} {t.day} {MONTHS[t.month - 1]} {t:%H:%M}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("out")
    ap.add_argument("--layers", default="wind,temp,rain")
    ap.add_argument("--cities", default=CITIES)
    ap.add_argument("--utc-offset", type=float, default=3.0)
    ap.add_argument("--steps", type=int, default=6, help="frames per output hour (default 6)")
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--particles", type=int, default=2500)
    ap.add_argument("--width", type=int, default=1280, help="video width in pixels (default 1280)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    # the model output on its Lambert grid
    w = cppwrf.open_wrfout(a.run, domain=1, xarray=False)
    lat, lon = np.asarray(w.lat, dtype=float), np.asarray(w.lon, dtype=float)
    times = [t + datetime.timedelta(hours=a.utc_offset) for t in w.times]
    tz = "TSİ" if a.utc_offset == 3 else f"UTC{a.utc_offset:+g}"
    u, v = np.asarray(w["U10"], dtype=float), np.asarray(w["V10"], dtype=float)
    if "COSALPHA" in w and "SINALPHA" in w:  # grid-relative -> earth-relative
        ca, sa = np.asarray(w["COSALPHA"], dtype=float), np.asarray(w["SINALPHA"], dtype=float)
        u, v = u * ca - v * sa, v * ca + u * sa
    fields = {
        "wind": np.hypot(u, v) * 3.6,
        "temp": np.asarray(dg.t2_celsius(w), dtype=float),
        "rain": np.nan_to_num(np.asarray(dg.precip_accum(w, hours=1), dtype=float)),
    }
    nt = len(times)

    # a regular lon / lat grid inside the domain (no empty corners)
    lon0, lon1 = lon[:, 0].max(), lon[:, -1].min()
    lat0, lat1 = lat[0, :].max(), lat[-1, :].min()
    W = a.width
    H = int(round(W * (lat1 - lat0) / (lon1 - lon0) / 2)) * 2
    nx, ny = W // 3, H // 3
    gx, gy = np.linspace(lon0, lon1, nx), np.linspace(lat0, lat1, ny)
    GX, GY = np.meshgrid(gx, gy)
    tri = Delaunay(np.column_stack([lon.ravel(), lat.ravel()]))

    def regrid(f):
        return LinearNDInterpolator(tri, f.ravel())(GX, GY)

    U = np.stack([regrid(u[k]) for k in range(nt)])
    V = np.stack([regrid(v[k]) for k in range(nt)])
    F = {name: np.stack([regrid(fields[name][k]) for k in range(nt)]) for name in a.layers.split(",")}

    # coastlines, borders and city dots, drawn once as a transparent overlay of the frame's size
    dpi = 100
    fig = plt.figure(figsize=(W / dpi, H / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1], projection=ccrs.PlateCarree())
    ax.set_extent([lon0, lon1, lat0, lat1], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE.with_scale("10m"), edgecolor="#ffffff", linewidth=0.8, alpha=0.75)
    ax.add_feature(cfeature.BORDERS.with_scale("10m"), edgecolor="#ffffff", linewidth=0.5, alpha=0.45)
    ax.add_feature(cfeature.LAKES.with_scale("10m"), facecolor="none", edgecolor="#ffffff", linewidth=0.4, alpha=0.5)
    ax.set_axis_off()
    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)
    fig.canvas.draw()
    overlay = np.asarray(fig.canvas.buffer_rgba()).copy()
    plt.close(fig)

    cities = [(n, la, lo) for n, la, lo in parse_cities(a.cities) if lon0 < lo < lon1 and lat0 < la < lat1]
    ci = [(n, lo, la, int(round((la - lat0) / (lat1 - lat0) * (ny - 1))),
           int(round((lo - lon0) / (lon1 - lon0) * (nx - 1)))) for n, la, lo in cities]

    rng = np.random.default_rng(1)
    dx_deg = (lon1 - lon0) / (nx - 1)
    dy_deg = (lat1 - lat0) / (ny - 1)
    coslat = np.cos(np.radians(0.5 * (lat0 + lat1)))
    # particle speed on screen: the domain's width in about 8 s of video for 30 km/h
    scale = (lon1 - lon0) / (8 * a.fps) / (30 / 3.6)
    trail = 10

    def sample(G, px, py):
        i = np.clip((px - lon0) / dx_deg, 0, nx - 1.001)
        j = np.clip((py - lat0) / dy_deg, 0, ny - 1.001)
        i0, j0 = i.astype(int), j.astype(int)
        fi, fj = i - i0, j - j0
        g = np.nan_to_num(G)
        return ((1 - fi) * (1 - fj) * g[j0, i0] + fi * (1 - fj) * g[j0, i0 + 1]
                + (1 - fi) * fj * g[j0 + 1, i0] + fi * fj * g[j0 + 1, i0 + 1])

    writer_cls, ext = None, "gif"
    ff = shutil.which("ffmpeg")
    if not ff:
        try:
            import imageio_ffmpeg
            ff = imageio_ffmpeg.get_ffmpeg_exe()
        except ImportError:
            ff = None
    if ff:
        matplotlib.rcParams["animation.ffmpeg_path"] = ff
        writer_cls, ext = animation.FFMpegWriter, "mp4"

    nframes = (nt - 1) * a.steps + 1
    for name in a.layers.split(","):
        title, unit, stops, fmt = LAYERS[name]
        cmap, norm = colormap(stops)
        G = F[name]
        fig = plt.figure(figsize=(W / dpi, H / dpi), dpi=dpi, facecolor=BG)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(lon0, lon1)
        ax.set_ylim(lat0, lat1)
        ax.set_axis_off()
        img = ax.imshow(G[0], extent=[lon0, lon1, lat0, lat1], origin="lower", cmap=cmap, norm=norm,
                        interpolation="bilinear", aspect="auto", zorder=1)
        ax.imshow(overlay, extent=[lon0, lon1, lat0, lat1], aspect="auto", zorder=3, interpolation="bilinear")
        lc = LineCollection([], linewidths=0.9, zorder=2, capstyle="round")
        ax.add_collection(lc)
        outline = [pe.withStroke(linewidth=2.5, foreground="#000000a0")]
        labels = []
        for n, lo, la, j, i in ci:
            ax.plot(lo, la, "o", ms=3, color="white", zorder=4, path_effects=outline)
            labels.append(ax.text(lo, la, "", color="white", fontsize=9, zorder=5, ha="left", va="bottom",
                                  path_effects=outline))
        hud = ax.text(0.012, 0.975, "", transform=ax.transAxes, color="white", fontsize=15, fontweight="bold",
                      va="top", zorder=6, path_effects=outline)
        ax.text(0.012, 0.915, f"{title} · CPPWRF · GFS {w.times[0]:%d.%m %H} UTC", transform=ax.transAxes,
                color="#d0d0d8", fontsize=10, va="top", zorder=6, path_effects=outline)
        # legend: the colour scale (stops evenly spaced) on a dark box in the lower right corner
        box = fig.add_axes([0.655, 0.012, 0.335, 0.085])
        box.set_axis_off()
        box.add_patch(plt.Rectangle((0, 0), 1, 1, transform=box.transAxes, color="#000000", alpha=0.45))
        cax = fig.add_axes([0.665, 0.052, 0.315, 0.02])
        ns = len(stops)
        cax.imshow(np.linspace(0, 1, 256)[None, :], cmap=cmap, aspect="auto", extent=[0, ns - 1, 0, 1])
        cax.set_yticks([])
        cax.set_xticks(range(ns))
        cax.set_xticklabels([f"{s[0]:g}" for s in stops], color="white", fontsize=7)
        cax.tick_params(length=0, pad=2)
        for s in cax.spines.values():
            s.set_visible(False)
        cax.text(0, 1.25, unit, transform=cax.transAxes, color="white", fontsize=8, va="bottom")

        n = a.particles
        px = rng.uniform(lon0, lon1, n)
        py = rng.uniform(lat0, lat1, n)
        age = rng.integers(0, 60, n)
        hist = np.repeat(np.stack([px, py], -1)[:, None, :], trail, axis=1)
        alphas = np.linspace(0.0, 0.8, trail - 1)

        def frame(f):
            nonlocal px, py, age, hist
            k, r = divmod(f, a.steps)
            k2 = min(k + 1, nt - 1)
            s = r / a.steps
            g = (1 - s) * G[k] + s * G[k2]
            uu, vv = (1 - s) * U[k] + s * U[k2], (1 - s) * V[k] + s * V[k2]
            img.set_data(g)
            for sub in range(2):  # two particle steps per frame
                pu, pv = sample(uu, px, py), sample(vv, px, py)
                px = px + pu * scale / coslat / 2
                py = py + pv * scale / 2
                age += 1
                dead = (age > 80) | (px < lon0) | (px > lon1) | (py < lat0) | (py > lat1) | (np.hypot(pu, pv) < 0.3)
                nd = int(dead.sum())
                if nd:
                    px[dead] = rng.uniform(lon0, lon1, nd)
                    py[dead] = rng.uniform(lat0, lat1, nd)
                    age[dead] = 0
                hist = np.roll(hist, -1, axis=1)
                hist[:, -1, 0], hist[:, -1, 1] = px, py
                hist[dead] = np.stack([px[dead], py[dead]], -1)[:, None, :]
            seg = np.stack([hist[:, :-1], hist[:, 1:]], axis=2).reshape(-1, 2, 2)
            col = np.zeros((n, trail - 1, 4))
            col[..., :3] = 1.0
            col[..., 3] = alphas[None, :] * np.minimum(1.0, age[:, None] / 8.0)
            lc.set_segments(seg)
            lc.set_color(col.reshape(-1, 4))
            for (nm, lo, la, j, i), t in zip(ci, labels):
                val = g[j, i]
                t.set_text(f" {nm}" if name == "rain" and val < 0.1 else f" {nm} {fmt.format(val)}")
            hud.set_text(f"{when(times[k] + datetime.timedelta(hours=s))} {tz}")
            return [img, lc, hud, *labels]

        for f in range(30):  # let the particles spread before the first frame
            frame(0)
        out = os.path.join(a.out, f"anim_{name}.{ext}")
        anim = animation.FuncAnimation(fig, frame, frames=nframes, blit=False)
        if writer_cls:
            anim.save(out, writer=writer_cls(fps=a.fps, bitrate=4000, extra_args=["-pix_fmt", "yuv420p"]), dpi=dpi)
        else:
            anim.save(out, writer=animation.PillowWriter(fps=a.fps), dpi=dpi)
        fig.savefig(os.path.join(a.out, f"anim_{name}.png"), dpi=dpi, facecolor=BG)
        plt.close(fig)
        print(f"{out}: {nframes} frames, {W}x{H}")


if __name__ == "__main__":
    main()
