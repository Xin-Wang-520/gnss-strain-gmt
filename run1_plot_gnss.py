#!/usr/bin/env python3
"""Run 1: plot original GNSS coverage and horizontal velocities, offline.

Usage:
    python3 run1_plot_gnss.py
    python3 run1_plot_gnss.py --region 67 113 22 45
    python3 run1_plot_gnss.py --show

Requires numpy and matplotlib. Reads both original v6.4 Eurasia-fixed files;
keeps GNSS only, including horizontal components from the 3D file.
No thinning, interpolation, or strain calculation. Local geographic outlines
were extracted from GMT 6.4's GSHHG database; GMT is not needed to run this file.
"""

from pathlib import Path
import argparse
import json
import os
import sys

ROOT = Path(__file__).resolve().parent
# Optional project-local dependencies; does not change the system environment.
LOCAL_PACKAGES = ROOT / "GNSS_analysis" / ".python_packages"
if LOCAL_PACKAGES.is_dir():
    sys.path.insert(0, str(LOCAL_PACKAGES))
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "GNSS_analysis" / ".cache" / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(ROOT / "GNSS_analysis" / ".cache"))

import numpy as np
import matplotlib


def arguments(description=None):
    parser = argparse.ArgumentParser(description=description or __doc__)
    parser.add_argument("--region", type=float, nargs=4, default=[70, 110, 25, 42],
                        metavar=("WEST", "EAST", "SOUTH", "NORTH"))
    parser.add_argument("--velocity-scale", type=float, default=0.030,
                        help="Arrow symbol length in map units per mm/yr")
    parser.add_argument("--dpi", type=int, default=240)
    parser.add_argument("--show", action="store_true", help="Also open an interactive figure")
    args = parser.parse_args()
    w, e, s, n = args.region
    if not (-180 <= w < e <= 180 and -85 < s < n < 85):
        parser.error("Invalid region; Mercator latitude must be between -85 and 85 degrees")
    if args.velocity_scale <= 0 or args.dpi <= 0:
        parser.error("Velocity scale and dpi must be positive")
    # Saved context currently covers this buffered region only.
    if w < 67 or e > 113 or s < 22 or n > 45:
        parser.error("Local geographic context covers 67-113 E, 22-45 N")
    return args


def mercator_y(latitude):
    """Spherical conformal Mercator ordinate, in degree-equivalent map units."""
    return np.rad2deg(np.log(np.tan(np.pi / 4 + np.deg2rad(latitude) / 2)))


def read_gnss(raw_dir, region, frame="eu"):
    w, e, s, n = region
    rows, kinds = [], []
    # The author's _itrf14 files also retain the veast_eu/vnorth_eu header
    # names. Select the reference frame by the source filename, not the header.
    fields = ("lon", "lat", "veast_eu", "vnorth_eu", "seast", "snorth")
    for kind in ("2D", "3D"):
        path = raw_dir / f"ahbgps_v6pt4_{kind}_29-Nov-2025_{frame}.dat"
        with path.open(encoding="utf-8") as stream:
            header = next(stream).split()
            numeric_indices = [header.index(field) for field in fields]
            type_index = header.index("types")
            for line in stream:
                parts = line.split()
                if not parts:
                    continue
                if len(parts) != len(header):
                    raise ValueError(f"Unexpected column count in {path.name}")
                if parts[type_index] != "GNSS":
                    continue
                values = [float(parts[i]) for i in numeric_indices]
                if not (w <= values[0] <= e and s <= values[1] <= n):
                    continue
                if not np.isfinite(values).all() or min(values[4:]) <= 0:
                    raise ValueError(f"Invalid velocity or uncertainty in {path.name}")
                rows.append(values)
                kinds.append(kind)
    if not rows:
        raise ValueError("No GNSS sites in selected region")
    data = np.asarray(rows)
    if len(np.unique(data[:, :2], axis=0)) != len(data):
        raise ValueError("Resolve duplicated coordinates before plotting")
    return data, np.asarray(kinds)


def draw_context(ax, region, basemap_dir):
    from matplotlib.ticker import FixedLocator, FuncFormatter

    ax.set_facecolor("#f9f9f5")
    for name, color, width in (("coastlines.txt", "#afc6d2", 0.55),
                               ("borders.txt", "#aeb5b9", 0.60)):
        lines = np.loadtxt(basemap_dir / name)
        ax.plot(lines[:, 0], mercator_y(lines[:, 1]), color=color,
                linewidth=width, zorder=1)
    w, e, s, n = region
    ax.set_xlim(w, e)
    ax.set_ylim(mercator_y(s), mercator_y(n))
    ax.set_aspect("equal", adjustable="box")
    longitude_ticks = np.arange(np.ceil(w / 10) * 10, e + 0.01, 10)
    latitude_ticks = np.arange(np.ceil(s / 5) * 5, n + 0.01, 5)
    ax.xaxis.set_major_locator(FixedLocator(longitude_ticks))
    ax.yaxis.set_major_locator(FixedLocator(mercator_y(latitude_ticks)))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda val, pos: f"{val:g}°E"))
    ax.set_yticklabels([f"{val:g}°N" for val in latitude_ticks])
    ax.grid(color="#d9e1e5", linewidth=0.6, zorder=0)
    ax.tick_params(labelsize=11, colors="#485661", direction="out")
    for spine in ax.spines.values():
        spine.set_color("#73818a")
        spine.set_linewidth(0.8)


def main(frame="eu", script_name=None):
    frame_labels = {"eu": "ITRF2014 Eurasia-fixed", "itrf14": "ITRF2014 global frame"}
    if frame not in frame_labels:
        raise ValueError(f"Unsupported reference frame: {frame}")
    frame_label = frame_labels[frame]
    args = arguments(f"Plot AHB v6.4 GNSS station coverage and horizontal velocities: {frame_label}.")
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    analysis = ROOT / "GNSS_analysis"
    output = analysis / "output" / ("run1" if frame == "eu" else "run1_itrf14")
    output.mkdir(parents=True, exist_ok=True)
    data, kinds = read_gnss(analysis / "data" / "raw", args.region, frame=frame)
    lon, lat, ve, vn = data[:, :4].T
    y = mercator_y(lat)
    is_2d = kinds == "2D"
    total = len(data)
    n2d, n3d = int(is_2d.sum()), int((~is_2d).sum())

    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42})
    fig, axes = plt.subplots(2, 1, figsize=(14, 15), facecolor="white")
    fig.subplots_adjust(top=0.925, bottom=0.06, left=0.075, right=0.98, hspace=0.22)
    fig.suptitle("Tibetan Plateau GNSS", fontsize=23, fontweight="bold", y=0.982, color="#233849")
    fig.text(0.5, 0.952, f"AHB v6.4 · {frame_label} · Original station velocities",
             ha="center", fontsize=11, color="#566773")
    for ax in axes:
        draw_context(ax, args.region, analysis / "data" / "basemap")

    axes[0].scatter(lon[is_2d], y[is_2d], s=8, color="#0e818b", edgecolors="none",
                    zorder=3, label=f"2D file: {n2d:,} sites")
    axes[0].scatter(lon[~is_2d], y[~is_2d], s=12, marker="^", color="#d57638",
                    edgecolors="none", zorder=3, label=f"3D file: {n3d:,} sites (horizontal components)")
    axes[0].set_title(f"a   Station coverage  |  {total:,} GNSS sites", loc="left",
                      fontsize=13, fontweight="bold", color="#233849", pad=12)
    legend = axes[0].legend(loc="upper left", fontsize=10, framealpha=0.95,
                            edgecolor="#cbd3d8", borderpad=0.8)
    legend.set_zorder(5)

    axes[1].scatter(lon, y, s=2, color="#80929f", edgecolors="none", zorder=2)
    # A fixed symbol scale, not displaced geographic positions. Because Mercator
    # is conformal, directions in this equal-aspect plane match EN azimuths.
    quiver = axes[1].quiver(lon, y, ve, vn, angles="xy", scale_units="xy",
                            scale=1 / args.velocity_scale, pivot="tail", color="#204f77",
                            width=0.00075, headwidth=3.8, headlength=4.8,
                            headaxislength=4.2, minshaft=1.5, minlength=0, zorder=3)
    axes[1].set_title(f"b   Horizontal velocity  |  All {total:,} vectors", loc="left",
                      fontsize=13, fontweight="bold", color="#233849", pad=12)
    # The key is generated from the plotted Quiver object, so the scale matches.
    axes[1].text(0.015, 0.972, "20 mm/yr\n\nOriginal velocities; error ellipses omitted",
                 transform=axes[1].transAxes, va="top", ha="left", fontsize=9,
                 color="#334957", bbox={"facecolor": "white", "edgecolor": "#cbd3d8", "alpha": 0.96,
                                           "boxstyle": "square,pad=0.6"}, zorder=5)
    axes[1].quiverkey(quiver, 0.22, 0.957, 20, "", coordinates="axes", color="#204f77", zorder=6)
    fig.text(0.075, 0.025, "No interpolation or thinning. Arrows show compiled horizontal GNSS velocities; units: mm/yr.",
             fontsize=9, color="#637580")

    output_stem = "run1_gnss" if frame == "eu" else "run1_gnss_itrf14"
    png = output / f"{output_stem}.png"
    pdf = output / f"{output_stem}.pdf"
    fig.savefig(png, dpi=args.dpi, facecolor="white")
    fig.savefig(pdf, dpi=args.dpi, facecolor="white")
    summary = {"script": script_name or Path(__file__).name, "data_version": "v6.4",
               "frame": frame_label, "region_w_e_s_n": args.region,
               "counts": {"2D": n2d, "3D": n3d}, "total": total,
               "velocity_units": "mm/yr", "velocity_symbol_scale": args.velocity_scale,
               "projection": "Spherical conformal Mercator; fixed-scale velocity symbols",
               "station_selection": "GNSS only; 2D plus horizontal components of 3D",
               "processing": "No thinning, interpolation, or strain calculation",
               "basemap_source": "Local GMT 6.4 GSHHG coastlines and borders",
               "input_files": [f"ahbgps_v6pt4_{kind}_29-Nov-2025_{frame}.dat" for kind in ("2D", "3D")]}
    if frame == "itrf14":
        summary["input_header_note"] = "Original ITRF2014 files retain veast_eu/vnorth_eu column names; frame selected by filename."
    (output / "run1_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"GNSS: {total:,} = {n2d:,} from 2D + {n3d:,} from 3D")
    print(f"Saved PNG: {png}\nSaved PDF: {pdf}")
    if args.show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
