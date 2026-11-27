#!/usr/bin/env python
"""
Per-tile GEDI extraction for the protected-area change analysis.

Reads one 1-degree tile of the MAAP tiled GEDI GeoParquet (gedi-tiled-v3), keeps
quality-filtered, leaf-on shots that fall inside the growing-season DOY window of their
PFT x region (power beams only), labels each shot T1 or T2, and writes:

  shots_<tile>.parquet   one row per kept shot (for hybrid inference / bootstrapping later)
  cells_<tile>.parquet   one row per 1 km EASE2 cell x epoch: n, mean, sd of each metric,
                         plus sampling diagnostics (night share, beams, sensitivity)

cell_id = row * 34704 + col on the GEDI 1 km EASE2 grid (the grid of
GEDI_ANCI_PFT_r1000m_EASE2.0). In R/terra the same cell is cell_id + 1.

Runs on NASA MAAP (Hub or DPS) with the platform's own AWS credentials; no Earthdata login.

Ground-drift correction (if the correction files are present, see correction.py): every shot
also gets drift-corrected RH98/RH50 and AGBD recomputed from corrected RH with the L4A V3 models,
in three variants (central, low, high drift rate). Uncorrected AGBD is recomputed too and
compared with agbd_l4a; the share of shots disagreeing by >0.1% is written to the summary
(l4a_mismatch_frac) as a continuous check of the recomputation.
New metrics: rh98_c, rh50_c, agbd_c (central), agbd_c_lo, agbd_c_hi. L2B metrics (PAI, FHD,
cover) are not corrected.

A summary_<tile>.json is always written (even when no shots pass), so the submit script can
tell finished tiles from missing ones.

Usage:
  python extract_tile.py --tile N00_E008                       # windows default to doy_windows.csv
  python extract_tile.py --tile N00_E008 --windows other.csv --out output
"""
import argparse
import json
import os

import duckdb
import numpy as np
import pandas as pd
from pyproj import Transformer

ROOT = "s3://nasa-maap-data-store/file-staging/nasa-map/gedi-tiled-v3"
EPOCHS = {"T1": (2019, 2021), "T2": (2024, 2026)}  # calendar years, inclusive
HEMI_SPLIT_REGIONS = (4,)  # Africa has separate N / S windows

# GEDI 1 km EASE-Grid 2.0 global (EPSG:6933), 34704 x 14616 cells
RES, X0, Y0, NCOLS, NROWS = 1000.89502334956, -17367530.4451615, 7314540.8306386, 34704, 14616

# Power beams only: coverage-beam quality degrades over the mission, so they are not used
# for change. beam_l2a 5, 6, 8, 11 = BEAM0101, BEAM0110, BEAM1000, BEAM1011.
POWER_BEAMS = (5, 6, 8, 11)
# Shot-level filter applied to everything. Edit here to change the quality definition.
BASE_FILTER = f"""
    beam_l2a IN {POWER_BEAMS}
    AND l2a_quality_flag_rel3_l2a = 1
    AND degrade_flag_l2a = 0
    AND land_cover_data_leaf_off_flag_l2a = 0
"""
# Per-metric validity: a metric is NULL for a shot when its own product flag fails.
L2B_OK = "l2b_quality_flag_rel3_l2b = 1"
L4A_OK = "l4a_quality_flag_rel3_l4a = 1"
METRICS = ["rh98", "rh50", "pai", "fhd", "cover", "agbd"]
CORR_METRICS = ["rh98_c", "rh50_c", "agbd_c", "agbd_c_lo", "agbd_c_hi"]
HERE = os.path.dirname(os.path.abspath(__file__))
CORR_FILES = {k: os.environ.get(f"GEDI_PA_{k.upper()}", os.path.join(HERE, f))
              for k, f in (("l4a_models", "l4a_model_data.json"),
                           ("l4a_overrides", "l4a_overrides.json"),
                           ("drift", "drift_correction.csv"))}


# In the container the windows file and DuckDB extensions are baked in (see Containerfile).
DEFAULT_WINDOWS = os.environ.get(
    "GEDI_PA_WINDOWS", os.path.join(os.path.dirname(os.path.abspath(__file__)), "doy_windows.csv"))


def connect():
    con = duckdb.connect()
    ext_dir = os.environ.get("DUCKDB_EXTENSION_DIR")
    if ext_dir and os.path.isdir(ext_dir):
        con.sql(f"SET extension_directory = '{ext_dir}'")
    con.sql("SET enable_progress_bar = false")
    con.sql("INSTALL httpfs; LOAD httpfs;")
    try:  # MAAP role credentials via boto3 (valid ~1 h, enough for one tile)
        import boto3
        c = boto3.Session().get_credentials().get_frozen_credentials()
        tok = f", SESSION_TOKEN '{c.token}'" if c.token else ""
        con.sql(f"CREATE OR REPLACE SECRET maap (TYPE s3, KEY_ID '{c.access_key}', "
                f"SECRET '{c.secret_key}'{tok}, REGION 'us-west-2')")
    except Exception as err:  # fall back to DuckDB's own credential chain
        print(f"boto3 credentials unavailable ({err}); trying credential_chain")
        con.sql("INSTALL aws; LOAD aws;")
        con.sql("CREATE OR REPLACE SECRET maap (TYPE s3, PROVIDER credential_chain, REGION 'us-west-2')")
    return con


def tile_files(con, manifest, tile, years):
    files = con.sql(f"SELECT f FROM read_csv('{manifest}', header=false, "
                    f"columns={{'f':'VARCHAR'}}, delim='\\t')").df().f
    keep = [f for f in files if f"/tile_id={tile}/" in f
            and any(f"/year={y}/" in f for y in years)]
    return keep


def read_shots(con, files, windows):
    con.register("w", windows)
    years = sorted({y for a, b in EPOCHS.values() for y in range(a, b + 1)})
    epoch_case = " ".join(f"WHEN yr BETWEEN {a} AND {b} THEN '{k}'" for k, (a, b) in EPOCHS.items())
    flist = ", ".join(f"'{f}'" for f in files)
    hemi = ", ".join(str(r) for r in HEMI_SPLIT_REGIONS)
    sql = f"""
    WITH s AS (
      SELECT shot_number,
             TIMESTAMP '2018-01-01' + to_seconds(delta_time_l2a) AS t,
             lon_lowestmode AS lon, lat_lowestmode AS lat,
             land_cover_data_pft_class_l2a AS pft_class,
             land_cover_data_region_class_l2a AS region_class,
             CASE WHEN land_cover_data_region_class_l2a IN ({hemi})
                  THEN (CASE WHEN lat_lowestmode >= 0 THEN 'N' ELSE 'S' END) ELSE 'all' END AS hemisphere,
             rh_l2a[99] AS rh98, rh_l2a[51] AS rh50,   -- rh_l2a holds RH0..RH100 (1-based index)
             rh_l2a[61] AS rh60, rh_l2a[71] AS rh70,
             predict_stratum_l4a AS stratum,
             CASE WHEN {L4A_OK} THEN agbd_l4a END AS agbd_l4a_q,
             land_cover_data_landsat_treecover_l2a AS tcover,
             CASE WHEN {L2B_OK} AND pai_l2b >= 0 THEN pai_l2b END AS pai,
             CASE WHEN {L2B_OK} AND fhd_normal_l2b >= 0 THEN fhd_normal_l2b END AS fhd,
             CASE WHEN {L2B_OK} AND cover_l2b >= 0 THEN cover_l2b END AS cover,
             CASE WHEN {L4A_OK} AND agbd_l4a >= 0 THEN agbd_l4a END AS agbd,
             sensitivity_l2a AS sensitivity,
             solar_elevation_l2a < 0 AS night,
             beam_l2a AS beam,
             granule
      FROM read_parquet([{flist}], union_by_name=true)
      WHERE {BASE_FILTER}
    ), d AS (
      SELECT *, year(t) AS yr, dayofyear(t) AS doy FROM s
    )
    SELECT d.*, CASE {epoch_case} END AS epoch
    FROM d JOIN w USING (pft_class, region_class, hemisphere)
    WHERE yr IN ({", ".join(map(str, years))})
      AND ((w.start_doy <= w.end_doy AND doy BETWEEN w.start_doy AND w.end_doy)
        OR (w.start_doy >  w.end_doy AND (doy >= w.start_doy OR doy <= w.end_doy)))
    """
    return con.sql(sql).df()


def add_cells(df):
    x, y = Transformer.from_crs("EPSG:4326", "EPSG:6933", always_xy=True).transform(
        df.lon.values, df.lat.values)
    col = np.floor((x - X0) / RES).astype(np.int64).clip(0, NCOLS - 1)
    row = np.floor((Y0 - y) / RES).astype(np.int64).clip(0, NROWS - 1)
    df["cell_id"] = row * NCOLS + col
    return df


def apply_correction(df, summary):
    """Add drift-corrected RH and recomputed AGBD columns (see correction.py)."""
    if not all(os.path.exists(f) for k, f in CORR_FILES.items() if k != "l4a_overrides"):
        print("correction files not found; skipping drift correction")
        return df, []
    import correction as C
    ov = CORR_FILES["l4a_overrides"] if os.path.exists(CORR_FILES["l4a_overrides"]) else None
    M = C.L4AModels(CORR_FILES["l4a_models"], ov)
    D = C.Drift(CORR_FILES["drift"])
    rh = {p: df[f"rh{p}"].to_numpy(float) for p in (50, 60, 70, 98)}
    missing = set(M.rh_needed()) - set(rh)
    if missing:
        raise SystemExit(f"model needs RH percentiles not read: {missing}")
    zero = np.zeros(len(df))
    # check: uncorrected recomputation vs product
    rec = M.predict(df.stratum.to_numpy(object), rh, zero)
    ok = df.agbd_l4a_q.notna() & np.isfinite(rec) & (df.agbd_l4a_q >= 0)
    rel = np.abs(rec[ok] - df.agbd_l4a_q[ok]) / np.maximum(df.agbd_l4a_q[ok].abs(), 1.0)
    summary["l4a_checked"] = int(ok.sum())
    summary["l4a_mismatch_frac"] = float((rel > 1e-3).mean()) if ok.any() else None
    bad = df.stratum[ok][rel.to_numpy() > 1e-3].value_counts()
    summary["l4a_mismatch_strata"] = {str(k): int(v) for k, v in bad.items()}
    valid_agbd = df.agbd.notna()  # same L4A quality filter as the uncorrected metric
    for which, col in (("central", "agbd_c"), ("low", "agbd_c_lo"), ("high", "agbd_c_hi")):
        shift = D.correction(df.tcover, df.yr, which)
        df[col] = np.where(valid_agbd, M.predict(df.stratum.to_numpy(object), rh, shift), np.nan)
        if which == "central":
            df["rh_shift_m"] = shift
            df["rh98_c"] = df.rh98 + shift
            df["rh50_c"] = df.rh50 + shift
    print(f"drift correction applied; L4A recomputation mismatch "
          f"{summary['l4a_mismatch_frac']} of {summary['l4a_checked']} shots")
    return df, CORR_METRICS


def aggregate(df, metrics=METRICS):
    g = df.groupby(["cell_id", "epoch"], sort=False)
    out = g.agg(n_shots=("shot_number", "size"), n_granules=("granule", "nunique"),
                pft_class=("pft_class", "first"), region_class=("region_class", "first"),
                hemisphere=("hemisphere", "first"), lon=("lon", "mean"), lat=("lat", "mean"),
                frac_night=("night", "mean"), n_beams=("beam", "nunique"),
                sensitivity_mean=("sensitivity", "mean"))
    for m in metrics:
        out[f"{m}_n"] = g[m].count()
        out[f"{m}_mean"] = g[m].mean()
        out[f"{m}_sd"] = g[m].std()
    return out.reset_index()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tile", required=True, help="tile id, e.g. N00_E008")
    ap.add_argument("--windows", default=DEFAULT_WINDOWS, help="DOY windows CSV")
    ap.add_argument("--out", default="output")
    ap.add_argument("--manifest", default=f"{ROOT}/manifest.txt")
    ap.add_argument("--no-shots", action="store_true", help="skip the shot-level file")
    a = ap.parse_args()
    summary = dict(tile=a.tile, n_files=0, n_shots=0, n_T1=0, n_T2=0, n_cell_rows=0)

    def finish():
        with open(os.path.join(a.out, f"summary_{a.tile}.json"), "w") as f:
            json.dump(summary, f)

    # Works with either the simple CSV or the full one (which repeats each stratum for T1 and T2).
    windows = pd.read_csv(a.windows)[["region_class", "hemisphere", "pft_class", "start_doy", "end_doy"]]
    windows = windows.drop_duplicates()
    assert not windows.duplicated(["region_class", "hemisphere", "pft_class"]).any(), \
        "windows file has more than one window for the same stratum"
    years = sorted({y for s, e in EPOCHS.values() for y in range(s, e + 1)})
    os.makedirs(a.out, exist_ok=True)
    con = connect()
    files = tile_files(con, a.manifest, a.tile, years)
    print(f"{a.tile}: {len(files)} parquet files for years {years}")
    summary["n_files"] = len(files)
    if not files:
        print("no data for this tile")
        return finish()

    df = read_shots(con, files, windows)
    print(f"kept {len(df):,} shots:", df.groupby("epoch").size().to_dict())
    summary.update(n_shots=len(df), n_T1=int((df.epoch == "T1").sum()), n_T2=int((df.epoch == "T2").sum()))
    if df.empty:
        print("no shots passed the filters")
        return finish()
    df = add_cells(df)
    df, extra = apply_correction(df, summary)
    if not a.no_shots:
        df.to_parquet(os.path.join(a.out, f"shots_{a.tile}.parquet"), index=False)
    cells = aggregate(df, METRICS + extra)
    cells.insert(0, "tile_id", a.tile)
    cells.to_parquet(os.path.join(a.out, f"cells_{a.tile}.parquet"), index=False)
    both = cells.groupby("cell_id").epoch.nunique().eq(2).sum()
    print(f"{len(cells):,} cell x epoch rows; {both:,} cells have both T1 and T2")
    summary["n_cell_rows"] = len(cells)
    finish()


if __name__ == "__main__":
    main()
