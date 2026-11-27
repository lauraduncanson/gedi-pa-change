#!/usr/bin/env python
"""
Does the canopy top drift, or only the ground?

Uses the shot samples saved by ground_drift_quicklook.py (no GEDI re-read). For each shot:
    ground residual = elev_lowestmode - TanDEM-X            (column `residual`)
    top residual    = ground residual + RH98                (elevation of RH98 minus TanDEM-X)
For each cover class, the drift rate (m/yr, 2020 onward, within tiles, gentle slopes) is fitted
for both. If the drift correction's assumptions hold:
    top rate    ~ open-ground rate (~0.03 m/yr) in every cover class   (whole-waveform shift only)
    ground rate rises with cover                                       (missed ground)
and the RH98 bias rate = ground rate - top rate ~ the excess ground rate used in the correction.

  python analysis/top_vs_ground_check.py --samples /home/jovyan/GEDI/results/ground_drift \
      /home/jovyan/GEDI/results/ground_drift_500
(folders or files; finds ground_drift_shot_sample.parquet or per-tile *_sample.parquet)
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

ORDER = ["0-10%", "10-30%", "30-50%", "50-70%", "70-80%", "80-90%", "90-100%"]


def tile_year_rate(d, col, min_n):
    """Within-tile rate: per tile x year mean minus tile mean, pooled WLS slope on year."""
    g = d.groupby(["tile_id", "yr"])[col].agg(["mean", "count"]).reset_index()
    g = g[g["count"] >= min_n]
    g = g[g.groupby("tile_id").yr.transform("nunique") >= 2]
    if g.empty:
        return np.nan, 0
    g["dev"] = g["mean"] - g.groupby("tile_id")["mean"].transform("mean")
    g["yc"] = g.yr - g.groupby("tile_id").yr.transform("mean")
    w = g["count"]
    return float(np.sum(w * g.yc * g.dev) / np.sum(w * g.yc ** 2)), int(g.tile_id.nunique())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", required=True,
                    help="run folders or sample parquet files")
    ap.add_argument("--min-n", type=int, default=5, help="min shots per tile x year x cover")
    ap.add_argument("--n-boot", type=int, default=300)
    a = ap.parse_args()

    files = []
    for path in a.samples:
        if os.path.isfile(path):
            files.append(path)
        else:
            files += glob.glob(os.path.join(path, "ground_drift_shot_sample.parquet"))
            files += glob.glob(os.path.join(path, "tiles", "*_sample.parquet"))
            files += glob.glob(os.path.join(path, "*_sample.parquet"))
    files = sorted(set(files))
    if not files:
        raise SystemExit(f"no shot-sample parquet files found in {a.samples}")
    print(f"reading {len(files)} sample file(s)")
    d = pd.concat([pd.read_parquet(f, columns=["tile_id", "yr", "cover_bin", "gentle",
                                               "residual", "rh98"]) for f in files])
    d = d.drop_duplicates()
    d = d[d.gentle & (d.yr >= 2020) & d.rh98.notna() & d.cover_bin.isin(ORDER)]
    d["top"] = d.residual + d.rh98
    print(f"{len(d):,} sampled shots from {d.tile_id.nunique()} tiles (2020 onward, gentle slopes)")

    rng = np.random.default_rng(0)
    rows = []
    for cb in ORDER:
        x = d[d.cover_bin == cb]
        gr, nt = tile_year_rate(x, "residual", a.min_n)
        tr, _ = tile_year_rate(x, "top", a.min_n)
        tiles = x.tile_id.unique()
        by = {t: v for t, v in x.groupby("tile_id")}
        diffs = []
        for _ in range(a.n_boot):
            s = pd.concat([by[t] for t in rng.choice(tiles, len(tiles))])
            g_, _ = tile_year_rate(s, "residual", a.min_n)
            t_, _ = tile_year_rate(s, "top", a.min_n)
            diffs.append((g_, t_, g_ - t_))
        b = np.array(diffs)
        lo, hi = np.nanpercentile(b, 2.5, axis=0), np.nanpercentile(b, 97.5, axis=0)
        rows.append(dict(cover_bin=cb, n_tiles=nt, ground_rate=gr, top_rate=tr,
                         top_lo=lo[1], top_hi=hi[1], rh98_bias_rate=gr - tr,
                         rh98_lo=lo[2], rh98_hi=hi[2]))
    out = pd.DataFrame(rows)
    with pd.option_context("display.width", 160, "display.float_format", "{:.3f}".format):
        print(out.to_string(index=False))
    print("\nIf top_rate is flat across cover (~ the 0-10% value), the canopy top is not drifting and"
          "\nthe ground-based correction applies to RH. rh98_bias_rate should then track the"
          "\nexcess ground rate in drift_correction.csv (sample-based, so noisier).")


if __name__ == "__main__":
    main()
