#!/usr/bin/env python
"""
Before/after (T1 vs T2) histograms of GEDI structure metrics for a tile.

For each tile it makes two figures:
  1. All shots, split into inside vs outside protected areas (WDPA 2025).
  2. The individual protected areas in the tile with the most data in both epochs.

Shots come from extract_tile.py (same filters: power beams, quality, growing-season windows).
If output/shots_<tile>.parquet already exists it is reused; otherwise the tile is extracted.

Notebook use:
    import plot_tile_histograms as P
    P.plot_tile("S16_W048")                      # one tile
    for t in ["S16_W048", "N11_W005", "S03_W060"]:
        P.plot_tile(t)

Note: "outside" is simply everything in the tile that is not in a protected area. It is not a
matched counterfactual, and T1 and T2 are different footprints, so read these as a first look.
"""
import os
import sys

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import box

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "extract_tile"))
import extract_tile as E  # noqa: E402

WDPA = ("/home/jovyan/shared-buckets/leitoldv/GEDI_global_change/WDPA_comparison/"
        "GLOBAL_PAs_all_nonMARINE_2025.gpkg")
WINDOWS = E.DEFAULT_WINDOWS
OUT = "output"
METRICS = {"rh98": "RH98 (m)", "rh50": "RH50 (m)", "pai": "PAI (m2/m2)",
           "fhd": "FHD", "agbd": "AGBD (Mg/ha)"}
NAME_COLS = ["NAME", "name", "ORIG_NAME", "NAME_ENG", "WDPAID", "WDPA_PID", "SITE_ID"]
COLORS = {"T1": "#3b6fb6", "T2": "#d9772b"}
MIN_SHOTS = 50  # a protected area needs this many shots in EACH epoch to get its own panel


def get_shots(tile):
    path = os.path.join(OUT, f"shots_{tile}.parquet")
    if os.path.exists(path):
        return pd.read_parquet(path)
    win = pd.read_csv(WINDOWS)[["region_class", "hemisphere", "pft_class",
                                "start_doy", "end_doy"]].drop_duplicates()
    con = E.connect()
    con.sql("SET enable_progress_bar = false")
    years = sorted({y for a, b in E.EPOCHS.values() for y in range(a, b + 1)})
    files = E.tile_files(con, f"{E.ROOT}/manifest.txt", tile, years)
    if not files:
        raise ValueError(f"tile {tile} is not in the manifest")
    df = E.add_cells(E.read_shots(con, files, win))
    os.makedirs(OUT, exist_ok=True)
    df.to_parquet(path, index=False)
    return df


def label_pas(shots, wdpa=WDPA):
    """Return (shots with in_pa flag, shot x PA table, name column used)."""
    pts = gpd.GeoDataFrame(shots, geometry=gpd.points_from_xy(shots.lon, shots.lat), crs="EPSG:4326")
    bb = gpd.GeoSeries([box(*pts.total_bounds)], crs="EPSG:4326")
    pas = gpd.read_file(wdpa, bbox=bb).to_crs("EPSG:4326")
    name_col = next((c for c in NAME_COLS if c in pas.columns), None)
    print(f"  {len(pas)} protected areas intersect the tile; name column: {name_col}")
    if pas.empty:
        return shots.assign(in_pa=False), None, name_col
    pas["geometry"] = pas.geometry.make_valid()
    pas["pa_label"] = pas[name_col].astype(str) if name_col else pas.index.astype(str)
    j = gpd.sjoin(pts, pas[["pa_label", "geometry"]], predicate="within", how="inner")
    out = shots.assign(in_pa=shots.shot_number.isin(j.shot_number))
    return out, pd.DataFrame(j.drop(columns=["geometry", "index_right"])), name_col


def _hist(ax, df, metric):
    v = df[metric].dropna()
    if v.empty:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
        return
    lo, hi = np.percentile(v, [0.5, 99.5])
    bins = np.linspace(lo, hi if hi > lo else lo + 1, 40)
    stats = {}
    for ep in ("T1", "T2"):
        x = df.loc[df.epoch == ep, metric].dropna()
        if len(x) == 0:
            continue
        ax.hist(x, bins=bins, density=True, histtype="stepfilled", alpha=0.35,
                color=COLORS[ep], edgecolor=COLORS[ep], linewidth=1.2,
                label=f"{ep}: n={len(x):,}, mean={x.mean():.2f}")
        ax.axvline(x.mean(), color=COLORS[ep], linestyle="--", linewidth=1)
        stats[ep] = x.mean()
    if len(stats) == 2:
        ax.set_title(f"T2 - T1 mean = {stats['T2'] - stats['T1']:+.2f}", fontsize=9, loc="right")
    ax.legend(fontsize=7, frameon=False)
    ax.set_yticks([])


def _grid(groups, tile, title, fname):
    """groups: ordered dict of column title -> shots DataFrame."""
    ncol, nrow = len(groups), len(METRICS)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 2.3 * nrow), squeeze=False)
    for c, (gname, gdf) in enumerate(groups.items()):
        for r, (m, mlabel) in enumerate(METRICS.items()):
            ax = axes[r, c]
            _hist(ax, gdf, m)
            ax.set_xlabel(mlabel, fontsize=8)
            if r == 0:
                ax.set_title(gname, fontsize=10, loc="left", fontweight="bold")
    fig.suptitle(f"{title} - tile {tile} (density; dashed = mean)", fontsize=11)
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, fname)
    fig.savefig(path, dpi=150)
    print(f"  saved {path}")
    plt.show()


def plot_tile(tile, n_pas=4, wdpa=WDPA):
    print(f"Tile {tile}")
    shots = get_shots(tile)
    print(f"  {len(shots):,} shots:", shots.groupby("epoch").size().to_dict())
    shots, joined, _ = label_pas(shots, wdpa)
    print("  shots inside protected areas:", shots.groupby(["in_pa", "epoch"]).size().to_dict())

    _grid({"Inside protected areas": shots[shots.in_pa],
           "Outside protected areas": shots[~shots.in_pa]},
          tile, "T1 vs T2, inside and outside protected areas", f"hist_{tile}_inside_outside.png")

    if joined is None or joined.empty:
        print("  no shots fall inside a protected area; skipping per-PA figure")
        return shots
    counts = joined.groupby(["pa_label", "epoch"]).size().unstack(fill_value=0)
    for ep in ("T1", "T2"):
        if ep not in counts:
            counts[ep] = 0
    counts["both"] = counts[["T1", "T2"]].min(axis=1)
    top = counts[counts.both >= MIN_SHOTS].sort_values("both", ascending=False).head(n_pas)
    if top.empty:
        print(f"  no protected area has {MIN_SHOTS}+ shots in both epochs; skipping per-PA figure")
        return shots
    print("  protected areas plotted (shots per epoch):")
    print(top[["T1", "T2"]].to_string())
    groups = {str(p)[:38]: joined[joined.pa_label == p] for p in top.index}
    _grid(groups, tile, "T1 vs T2 by protected area", f"hist_{tile}_by_pa.png")
    return shots
