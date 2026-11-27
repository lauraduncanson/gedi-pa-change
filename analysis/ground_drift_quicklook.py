#!/usr/bin/env python
"""
Quick look at ground-finding drift: GEDI ground relative to the TanDEM-X DEM, T1 vs T2.

For every shot passing the analysis filters (power beams, quality, leaf-on, growing-season
window; same as extract_tile.py):

    residual = elev_lowestmode_l2a - digital_elevation_model_l2a      (TanDEM-X 90 m)

Terrain and the DEM's own errors (including its height above true ground under canopy) are
fixed in time, so within a stratum

    delta = mean residual (T2) - mean residual (T1)

estimates the change in GEDI ground finding. delta > 0 means T2 ground sits higher, so T2 RH
metrics are too low by delta. If missed ground under dense canopy is the cause, delta should be
~0 at low tree cover and grow with cover. Tree cover is Landsat (land_cover_data_landsat_treecover),
which does not depend on GEDI.

Slope proxy: along-track DEM gradient between neighbouring kept shots (same granule and beam,
< 150 m apart). Delta is reported for gentle terrain (slope < --max-slope degrees) and all terrain.

Memory: each tile is reduced to per-stratum sums (count, sum, sum of squares) before the next
one is read, so memory stays small whatever the number of tiles. A small random sample of shots
(--sample-frac) is kept for later plotting.

Outputs in --out: ground_drift_stats.parquet (the per-tile sums), ground_drift_shot_sample.parquet,
ground_drift_by_cover.csv, ground_drift_by_stratum.csv, ground_drift_by_year.csv, ground_drift.png

  python analysis/ground_drift_quicklook.py --n-tiles 150 --out /home/jovyan/GEDI/results/ground_drift
  python analysis/ground_drift_quicklook.py --tiles S13_E028 N11_W005 ...   # specific tiles
"""
import argparse
import os
import random
import sys
from concurrent.futures import ThreadPoolExecutor

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "extract_tile"))
sys.path.insert(0, os.path.join(HERE, "..", "dps"))
import extract_tile as E  # noqa: E402

REGIONS = {1: "Europe", 2: "North Asia", 3: "Australasia", 4: "Africa",
           5: "South Asia", 6: "South America", 7: "North America"}
PFTS = {1: "ENT", 2: "EBT", 3: "DNT", 4: "DBT", 5: "Shrub", 6: "Grass"}
COVER_BINS = [0, 10, 30, 50, 70, 80, 90, 101]
MAX_ABS_RESIDUAL = 50.0  # m; drops DEM voids, water, clouds and gross ground errors


KEYS = ["region_class", "pft_class", "cover_bin", "yr", "epoch", "gentle"]


def read_tile(tile, windows, manifest, max_slope, sample_frac, mem_limit):
    con = E.connect()
    con.sql(f"SET memory_limit = '{mem_limit}'")
    con.sql("SET threads = 2")
    years = sorted({y for a, b in E.EPOCHS.values() for y in range(a, b + 1)})
    files = E.tile_files(con, manifest, tile, years)
    if not files:
        return None
    con.register("w", windows)
    flist = ", ".join(f"'{f}'" for f in files)
    hemi = ", ".join(str(r) for r in E.HEMI_SPLIT_REGIONS)
    epoch_case = " ".join(f"WHEN yr BETWEEN {a} AND {b} THEN '{k}'" for k, (a, b) in E.EPOCHS.items())
    df = con.sql(f"""
    WITH s AS (
      SELECT shot_number, granule, beam_l2a AS beam,
             TIMESTAMP '2018-01-01' + to_seconds(delta_time_l2a) AS t,
             lon_lowestmode AS lon, lat_lowestmode AS lat,
             land_cover_data_pft_class_l2a AS pft_class,
             land_cover_data_region_class_l2a AS region_class,
             CASE WHEN land_cover_data_region_class_l2a IN ({hemi})
                  THEN (CASE WHEN lat_lowestmode >= 0 THEN 'N' ELSE 'S' END) ELSE 'all' END AS hemisphere,
             land_cover_data_landsat_treecover_l2a AS cover,
             elev_lowestmode_l2a - digital_elevation_model_l2a AS residual,
             digital_elevation_model_l2a AS dem,
             sensitivity_l2a AS sensitivity,
             rh_l2a[99] AS rh98, rh_l2a[51] AS rh50
      FROM read_parquet([{flist}], union_by_name=true)
      WHERE {E.BASE_FILTER}
    ), d AS (SELECT *, year(t) AS yr, dayofyear(t) AS doy FROM s)
    SELECT d.*, CASE {epoch_case} END AS epoch
    FROM d JOIN w USING (pft_class, region_class, hemisphere)
    WHERE yr IN ({", ".join(map(str, years))})
      AND ((w.start_doy <= w.end_doy AND doy BETWEEN w.start_doy AND w.end_doy)
        OR (w.start_doy >  w.end_doy AND (doy >= w.start_doy OR doy <= w.end_doy)))
      AND abs(residual) < {MAX_ABS_RESIDUAL}
      AND cover BETWEEN 0 AND 100
    """).df()
    if df.empty:
        return None
    df["tile_id"] = tile
    df = add_slope(df)
    df["cover_bin"] = cover_bins(df.cover)
    df["gentle"] = df.slope_deg < max_slope
    df["r2"] = df.residual ** 2
    stats = (df.groupby(KEYS, observed=True)
               .agg(n=("residual", "size"), s=("residual", "sum"), ss=("r2", "sum"))
               .reset_index().assign(tile_id=tile))
    sample = df.sample(frac=sample_frac, random_state=0).drop(columns=["t", "r2"])
    return stats, sample


def cover_bins(c):
    return pd.cut(c, COVER_BINS, right=False,
                  labels=[f"{lo}-{hi - 1 if hi == 101 else hi}%"
                          for lo, hi in zip(COVER_BINS[:-1], COVER_BINS[1:])]).astype(str)


def add_slope(df):
    """Along-track slope (degrees) from the DEM difference to the next kept shot."""
    df = df.sort_values(["granule", "beam", "t"]).reset_index(drop=True)
    same = (df.granule.eq(df.granule.shift(-1)) & df.beam.eq(df.beam.shift(-1)))
    lat1, lon1 = np.radians(df.lat), np.radians(df.lon)
    lat2, lon2 = np.radians(df.lat.shift(-1)), np.radians(df.lon.shift(-1))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    dist = 2 * 6371000 * np.arcsin(np.sqrt(a))
    dz = (df.dem.shift(-1) - df.dem).abs()
    ok = same & dist.between(20, 150)
    df["slope_deg"] = np.where(ok, np.degrees(np.arctan(dz / dist)), np.nan)
    return df


def delta_table(st, by):
    """T2 - T1 residual shift from summed statistics (n, s, ss) grouped by `by`."""
    g = st.groupby(by + ["epoch"])[["n", "s", "ss"]].sum()
    g["mean"] = g.s / g.n
    g["var"] = (g.ss - g.n * g["mean"] ** 2) / (g.n - 1).clip(lower=1)
    w = g[["n", "mean", "var"]].unstack("epoch").dropna()
    if w.empty:
        return pd.DataFrame()
    out = pd.DataFrame({
        "n_T1": w[("n", "T1")].astype(int), "n_T2": w[("n", "T2")].astype(int),
        "resid_T1": w[("mean", "T1")], "resid_T2": w[("mean", "T2")],
        "delta": w[("mean", "T2")] - w[("mean", "T1")],
        "delta_se": np.sqrt(w[("var", "T1")] / w[("n", "T1")] + w[("var", "T2")] / w[("n", "T2")]),
    })
    return out.reset_index()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tiles", nargs="*", help="tile ids (default: random sample)")
    ap.add_argument("--n-tiles", type=int, default=150)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--windows", default=E.DEFAULT_WINDOWS)
    ap.add_argument("--manifest", default=f"{E.ROOT}/manifest.txt")
    ap.add_argument("--max-slope", type=float, default=10.0)
    ap.add_argument("--min-shots", type=int, default=200, help="per epoch, for a stratum row")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--mem-limit", default="3GB", help="DuckDB memory limit per worker")
    ap.add_argument("--sample-frac", type=float, default=0.01, help="share of shots kept for plotting")
    ap.add_argument("--out", default="results/ground_drift")
    a = ap.parse_args()

    tiles = a.tiles
    if not tiles:
        import common as C
        all_t = C.all_tiles(a.manifest)
        random.Random(a.seed).shuffle(all_t)
        tiles = all_t[:a.n_tiles]
    windows = pd.read_csv(a.windows)[["region_class", "hemisphere", "pft_class",
                                      "start_doy", "end_doy"]].drop_duplicates()
    print(f"reading {len(tiles)} tiles with {a.workers} workers")
    stats, samples = [], []
    job = lambda t: read_tile(t, windows, a.manifest, a.max_slope, a.sample_frac, a.mem_limit)
    with ThreadPoolExecutor(a.workers) as ex:
        for i, res in enumerate(ex.map(job, tiles), 1):
            if res is not None:
                stats.append(res[0])
                samples.append(res[1])
            if i % 10 == 0:
                print(f"  {i}/{len(tiles)} tiles read", flush=True)
    if not stats:
        raise SystemExit("no shots found")
    st = pd.concat(stats, ignore_index=True)
    st["region"] = st.region_class.map(REGIONS)
    st["pft"] = st.pft_class.map(PFTS)
    os.makedirs(a.out, exist_ok=True)
    st.to_parquet(os.path.join(a.out, "ground_drift_stats.parquet"), index=False)
    pd.concat(samples, ignore_index=True).to_parquet(
        os.path.join(a.out, "ground_drift_shot_sample.parquet"), index=False)
    gentle = st[st.gentle]
    n_ep = st.groupby("epoch").n.sum()
    print(f"{int(st.n.sum()):,} shots ({int(n_ep.get('T1', 0)):,} T1, {int(n_ep.get('T2', 0)):,} T2) "
          f"from {st.tile_id.nunique()} tiles; {int(gentle.n.sum()):,} on slopes < {a.max_slope} deg")

    by_cover = pd.concat([delta_table(gentle, ["cover_bin"]).assign(terrain=f"slope<{a.max_slope}"),
                          delta_table(st, ["cover_bin"]).assign(terrain="all")])
    order = {b: i for i, b in enumerate(cover_bins(pd.Series(COVER_BINS[:-1])))}
    by_cover = by_cover.sort_values(["terrain", "cover_bin"], key=lambda c: c.map(order) if c.name == "cover_bin" else c)
    by_cover.to_csv(os.path.join(a.out, "ground_drift_by_cover.csv"), index=False)
    strat = delta_table(gentle, ["region", "pft", "cover_bin"])
    if len(strat):
        strat = strat[(strat.n_T1 >= a.min_shots) & (strat.n_T2 >= a.min_shots)]
    strat.to_csv(os.path.join(a.out, "ground_drift_by_stratum.csv"), index=False)
    yg = gentle.groupby(["cover_bin", "yr"])[["n", "s"]].sum().reset_index()
    by_year = yg.assign(mean=yg.s / yg.n).rename(columns={"n": "count"})[["cover_bin", "yr", "mean", "count"]]
    by_year.to_csv(os.path.join(a.out, "ground_drift_by_year.csv"), index=False)

    with pd.option_context("display.width", 160, "display.float_format", "{:.2f}".format):
        print("\ndelta (T2 - T1 ground residual, m) by Landsat tree cover:")
        print(by_cover.to_string(index=False))
        print(f"\n{len(strat)} region x PFT x cover strata with >= {a.min_shots} shots per epoch "
              f"(ground_drift_by_stratum.csv)")
    print("\nIf missed ground under dense canopy is the cause, delta is ~0 at low cover and rises "
          "with cover.\nA positive delta means T2 RH metrics are too low by about delta.")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    g = by_cover[by_cover.terrain != "all"].reset_index(drop=True)
    x = np.arange(len(g))
    axes[0].errorbar(x, g.delta, yerr=1.96 * g.delta_se, fmt="o-", capsize=3)
    axes[0].axhline(0, color="grey", lw=0.8)
    axes[0].set_xticks(x, g.cover_bin.astype(str), rotation=30)
    axes[0].set_xlabel("Landsat tree cover")
    axes[0].set_ylabel("T2 - T1 ground residual vs TanDEM-X (m)")
    axes[0].set_title(f"Ground-finding drift (slope < {a.max_slope}°)", fontsize=10)
    for cb, gg in by_year.groupby("cover_bin"):
        if str(cb) in ("0-10%", "50-70%", "90-100%"):
            gg = gg[gg["count"] >= a.min_shots]
            axes[1].plot(gg.yr, gg["mean"], "o-", label=f"cover {cb}")
    axes[1].set_xlabel("Year")
    axes[1].set_ylabel("Mean ground residual vs TanDEM-X (m)")
    axes[1].set_title("Residual by year (offset is DEM-dependent; the trend matters)", fontsize=10)
    axes[1].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    png = os.path.join(a.out, "ground_drift.png")
    fig.savefig(png, dpi=150)
    print(f"\nwrote {a.out}/ (stats + shot sample parquet, 3 CSVs, ground_drift.png)")


if __name__ == "__main__":
    main()
