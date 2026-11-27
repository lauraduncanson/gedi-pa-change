#!/usr/bin/env python
"""
Composition-robust version of the ground-drift quick look.

Uses ground_drift_stats.parquet from ground_drift_quicklook.py (no re-reading of GEDI data).
Every comparison is made within a tile (same place, different time) and then averaged across
tiles, so differences in where the ISS tracks fell each year cannot masquerade as drift.
Confidence intervals come from resampling whole tiles (shots within a tile are not independent).

  1. delta per cover class: within each tile, mean residual T2 - T1 (gentle slopes), weighted
     across tiles by n_T1*n_T2/(n_T1+n_T2).
  2. year effect per cover class: within each tile, mean residual in year Y minus the tile's
     T1 mean, averaged across tiles. Shows whether the drift is gradual or a step after 2023.

  python analysis/ground_drift_within_tile.py --stats /home/jovyan/GEDI/results/ground_drift/ground_drift_stats.parquet
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ORDER = ["0-10%", "10-30%", "30-50%", "50-70%", "70-80%", "80-90%", "90-100%"]


def tile_means(st, by):
    g = st.groupby(["tile_id"] + by)[["n", "s"]].sum()
    g["mean"] = g.s / g.n
    return g.reset_index()


def within_delta(st, min_n):
    m = tile_means(st, ["cover_bin", "epoch"]).pivot_table(
        index=["tile_id", "cover_bin"], columns="epoch", values=["mean", "n"])
    m = m.dropna()
    m = m[(m[("n", "T1")] >= min_n) & (m[("n", "T2")] >= min_n)]
    d = pd.DataFrame({"delta": m[("mean", "T2")] - m[("mean", "T1")],
                      "w": m[("n", "T1")] * m[("n", "T2")] / (m[("n", "T1")] + m[("n", "T2")])})
    return d.reset_index()


def year_effect(st, min_n):
    base = tile_means(st[st.epoch == "T1"], ["cover_bin"])
    yr = tile_means(st, ["cover_bin", "yr"])
    j = yr.merge(base[["tile_id", "cover_bin", "mean", "n"]], on=["tile_id", "cover_bin"],
                 suffixes=("", "_t1"))
    j = j[(j.n >= min_n) & (j.n_t1 >= min_n)]
    j["eff"] = j["mean"] - j["mean_t1"]
    j["w"] = j.n * j.n_t1 / (j.n + j.n_t1)
    return j


def wmean_ci(df, val, groups, n_boot, rng):
    rows = []
    tiles = df.tile_id.unique()
    for key, g in df.groupby(groups):
        est = np.average(g[val], weights=g.w)
        bt = []
        by_tile = {t: x for t, x in g.groupby("tile_id")}
        present = [t for t in tiles if t in by_tile]
        for _ in range(n_boot):
            pick = rng.choice(present, size=len(present), replace=True)
            x = pd.concat([by_tile[t] for t in pick])
            bt.append(np.average(x[val], weights=x.w))
        lo, hi = np.percentile(bt, [2.5, 97.5])
        key = key if isinstance(key, tuple) else (key,)
        rows.append(dict(zip(groups, key), estimate=est, ci_lo=lo, ci_hi=hi, n_tiles=len(present)))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", required=True)
    ap.add_argument("--out", default=None, help="default: same folder as --stats")
    ap.add_argument("--min-n", type=int, default=30, help="min shots per tile x cover x period")
    ap.add_argument("--n-boot", type=int, default=500)
    ap.add_argument("--all-slopes", action="store_true")
    a = ap.parse_args()
    out = a.out or os.path.dirname(os.path.abspath(a.stats))
    st = pd.read_parquet(a.stats)
    if not a.all_slopes:
        st = st[st.gentle]
    st = st[st.cover_bin.isin(ORDER)]
    rng = np.random.default_rng(0)

    d = within_delta(st, a.min_n)
    dc = wmean_ci(d, "delta", ["cover_bin"], a.n_boot, rng)
    dc["cover_bin"] = pd.Categorical(dc.cover_bin, ORDER, ordered=True)
    dc = dc.sort_values("cover_bin")
    ye = year_effect(st, a.min_n)
    yc = wmean_ci(ye, "eff", ["cover_bin", "yr"], a.n_boot, rng)
    dc.to_csv(os.path.join(out, "ground_drift_within_tile_by_cover.csv"), index=False)
    yc.to_csv(os.path.join(out, "ground_drift_within_tile_by_year.csv"), index=False)
    with pd.option_context("display.width", 150, "display.float_format", "{:.2f}".format):
        print("within-tile T2 - T1 ground shift (m), 95% CI from resampling tiles:")
        print(dc.to_string(index=False))

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    x = np.arange(len(dc))
    axes[0].errorbar(x, dc.estimate, yerr=[dc.estimate - dc.ci_lo, dc.ci_hi - dc.estimate],
                     fmt="o-", capsize=3)
    axes[0].axhline(0, color="grey", lw=0.8)
    axes[0].set_xticks(x, dc.cover_bin.astype(str), rotation=30)
    axes[0].set_xlabel("Landsat tree cover")
    axes[0].set_ylabel("Within-tile T2 - T1 ground shift (m)")
    axes[0].set_title(f"Ground-finding drift, same places ({int(dc.n_tiles.max())} tiles)", fontsize=10)
    for cb in ("0-10%", "50-70%", "90-100%"):
        g = yc[yc.cover_bin == cb].sort_values("yr")
        if len(g):
            axes[1].errorbar(g.yr, g.estimate, yerr=[g.estimate - g.ci_lo, g.ci_hi - g.estimate],
                             fmt="o-", capsize=3, label=f"cover {cb}")
    axes[1].axhline(0, color="grey", lw=0.8)
    axes[1].set_xlabel("Year")
    axes[1].set_ylabel("Ground residual minus tile's T1 mean (m)")
    axes[1].set_title("Year effect within tiles (gradual drift or step?)", fontsize=10)
    axes[1].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    png = os.path.join(out, "ground_drift_within_tile.png")
    fig.savefig(png, dpi=150)
    print(f"\nwrote {png} and two CSVs")


if __name__ == "__main__":
    main()
