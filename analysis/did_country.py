#!/usr/bin/env python
"""
Difference-in-differences for one country from the matching results and the GEDI cell table.

Reads the matching step 4 output (a folder of *_pa_<id>_matching_results*.RDS files, or one
combined RDS; one row per matched cell with lon, lat, status, PA id), maps each row to its 1 km
EASE2 cell_id, joins to the T1/T2 cell means, and for each protected area computes:

  change_pa   = mean(T2 cells in PA)      - mean(T1 cells in PA)
  change_ctl  = mean(T2 matched controls) - mean(T1 matched controls)
  did         = change_pa - change_ctl

Means are over cells (each cell's mean of its shots), because T1 and T2 rarely share a cell.
Controls reused by several PAs count for each of them (matches are kept per PA, not rasterized).
A PA enters the country summary only if all four groups (PA/control x T1/T2) have at least
--min-cells cells with data for that metric. Standard errors treat cells as independent, so
they are too small where cells are spatially correlated; read them as a lower bound.

Outputs: <ISO3>_did_by_pa.csv (one row per PA x metric; sd_pa_t1/sd_ctl_t1 = SD of T1 cell means,
for the balance rule) and a printed country summary.

  python analysis/did_country.py --iso3 BRA --matches <folder with BRA matching RDS files> \
      --cells output/gedi_pa_change_cells_combined.parquet --out /home/jovyan/GEDI/results
"""
import argparse
import glob
import os
import re

import duckdb
import numpy as np
import pandas as pd
import pyreadr
from pyproj import Transformer

RES, X0, Y0, NCOLS, NROWS = 1000.89502334956, -17367530.4451615, 7314540.8306386, 34704, 14616
METRICS = ["rh98", "rh50", "pai", "fhd", "cover", "agbd",
           "rh98_c", "rh50_c", "agbd_c", "agbd_c_lo", "agbd_c_hi"]  # _c: drift-corrected (drift-v1 tables only)


def read_matches(path, iso3):
    """Matched cells from step 4: a folder of per-PA RDS files, or one combined RDS file.

    Needs columns status, lon, lat and a PA id (pa_id or WDPAID); per-PA files without an id
    column take it from the file name (..._pa_<id>_matching_results...RDS).
    """
    if os.path.isfile(path):
        files = [path]
    else:
        files = sorted(glob.glob(os.path.join(path, "**", "*_pa_*_matching_results*.RDS"), recursive=True))
    frames = []
    for f in files:
        r = pyreadr.read_r(f)
        df = next(iter(r.values()), None) if r else None
        if df is None or df.empty:
            continue  # PA with no matches (the R code saves NULL)
        df = df.copy()
        if "pa_id" not in df:
            if "WDPAID" in df:
                df["pa_id"] = df["WDPAID"]
            else:
                m = re.search(r"_pa_(\d+)_matching", os.path.basename(f))
                if not m:
                    raise SystemExit(f"{f}: no pa_id/WDPAID column; columns are {list(df.columns)}")
                df["pa_id"] = int(m.group(1))
        missing = [c for c in ("status", "lon", "lat") if c not in df]
        if missing:
            raise SystemExit(f"{f}: missing columns {missing}; columns are {list(df.columns)}")
        frames.append(df[["pa_id", "status", "lon", "lat"]])
    if not frames:
        raise SystemExit(f"no matching results found for {iso3} at {path}")
    out = pd.concat(frames, ignore_index=True)
    out["status"] = out.status.astype(str).str.upper().isin(["TRUE", "1", "1.0"])
    out["pa_id"] = pd.to_numeric(out.pa_id).astype("int64")
    print(f"{iso3}: {len(files)} result file(s), {out.pa_id.nunique()} PAs with matches, "
          f"{out.status.sum():,} PA cells, {(~out.status).sum():,} control cells")
    return out


def add_cell_id(df):
    x, y = Transformer.from_crs("EPSG:4326", "EPSG:6933", always_xy=True).transform(
        df.lon.values, df.lat.values)
    col = np.floor((x - X0) / RES).astype(np.int64).clip(0, NCOLS - 1)
    row = np.floor((Y0 - y) / RES).astype(np.int64).clip(0, NROWS - 1)
    df["cell_id"] = row * NCOLS + col
    return df


def did_by_pa(matches, cells_path, min_cells):
    con = duckdb.connect()
    con.register("m", matches[["pa_id", "status", "cell_id"]].drop_duplicates())
    have = set(con.sql(f"DESCRIBE SELECT * FROM read_parquet('{cells_path}')").df().column_name)
    metrics = [x for x in METRICS if f"{x}_mean" in have]
    cols = ", ".join(f"c.{x}_mean AS {x}" for x in metrics)
    joined = con.sql(f"""
        SELECT m.pa_id, m.status, c.epoch, c.cell_id, c.n_shots, {cols}
        FROM m JOIN read_parquet('{cells_path}') c USING (cell_id)
    """).df()
    rows = []
    for metric in metrics:
        d = joined.dropna(subset=[metric])
        g = d.groupby(["pa_id", "status", "epoch"])[metric].agg(["mean", "std", "count"]).reset_index()
        g["se"] = g["std"] / np.sqrt(g["count"])
        w = g.pivot_table(index="pa_id", columns=["status", "epoch"], values=["mean", "std", "se", "count"])
        for pa in w.index:
            def get(stat, status, epoch):
                try:
                    return w.loc[pa, (stat, status, epoch)]
                except KeyError:
                    return np.nan
            n = {k: get("count", s, e) for k, (s, e) in
                 dict(pa_t1=(True, "T1"), pa_t2=(True, "T2"), ctl_t1=(False, "T1"), ctl_t2=(False, "T2")).items()}
            mu = {k: get("mean", s, e) for k, (s, e) in
                  dict(pa_t1=(True, "T1"), pa_t2=(True, "T2"), ctl_t1=(False, "T1"), ctl_t2=(False, "T2")).items()}
            se = {k: get("se", s, e) for k, (s, e) in
                  dict(pa_t1=(True, "T1"), pa_t2=(True, "T2"), ctl_t1=(False, "T1"), ctl_t2=(False, "T2")).items()}
            sd = {k: get("std", s, e) for k, (s, e) in
                  dict(pa_t1=(True, "T1"), ctl_t1=(False, "T1")).items()}
            ch_pa, ch_ctl = mu["pa_t2"] - mu["pa_t1"], mu["ctl_t2"] - mu["ctl_t1"]
            rows.append(dict(pa_id=pa, metric=metric,
                             **{f"n_{k}": n[k] for k in n}, **{f"mean_{k}": mu[k] for k in mu},
                             **{f"sd_{k}": sd[k] for k in sd},
                             change_pa=ch_pa, change_ctl=ch_ctl, did=ch_pa - ch_ctl,
                             did_se=np.sqrt(np.nansum([se[k] ** 2 for k in se])),
                             qualifies=all(np.nan_to_num(n[k]) >= min_cells for k in n)))
    return pd.DataFrame(rows)


def summarize(res, iso3):
    print(f"\n{iso3} summary (PAs with enough cells in all four groups; weights = PA cells)")
    out = []
    for metric, d in res.groupby("metric", sort=False):
        q = d[d.qualifies]
        n_all = d.pa_id.nunique()
        if q.empty:
            out.append(dict(metric=metric, n_pa=0, n_pa_matched=n_all))
            continue
        w = q.n_pa_t1 + q.n_pa_t2
        avg = lambda col: np.average(q[col], weights=w)
        out.append(dict(metric=metric, n_pa=len(q), n_pa_matched=n_all,
                        t1_pa=avg("mean_pa_t1"), t1_ctl=avg("mean_ctl_t1"),
                        change_pa=avg("change_pa"), change_ctl=avg("change_ctl"),
                        did=avg("did"), did_median=q.did.median(),
                        frac_pa_did_negative=(q.did < 0).mean()))
    s = pd.DataFrame(out)
    with pd.option_context("display.width", 200, "display.float_format", "{:.3f}".format):
        print(s.to_string(index=False))
    print("t1_pa vs t1_ctl is a balance check: matched controls should start from similar structure.")
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iso3", required=True)
    ap.add_argument("--matches", required=True,
                    help="folder of per-PA step 4 RDS files, or one combined RDS file")
    ap.add_argument("--cells", required=True, help="merged cell table (dps/collect.py)")
    ap.add_argument("--out", default="results")
    ap.add_argument("--min-cells", type=int, default=10)
    a = ap.parse_args()

    m = add_cell_id(read_matches(a.matches, a.iso3))
    res = did_by_pa(m, a.cells, a.min_cells)
    os.makedirs(a.out, exist_ok=True)
    res.to_csv(os.path.join(a.out, f"{a.iso3}_did_by_pa.csv"), index=False)
    s = summarize(res, a.iso3)
    s.to_csv(os.path.join(a.out, f"{a.iso3}_did_summary.csv"), index=False)
    print(f"\nwrote {a.out}/{a.iso3}_did_by_pa.csv and {a.iso3}_did_summary.csv")


if __name__ == "__main__":
    main()
