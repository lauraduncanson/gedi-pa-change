#!/usr/bin/env python
"""
Global change in GEDI structure inside protected areas, without matching (secondary result;
rules in CLAUDE.md, "Global PA change without matching").

  # 1. which 1 km cells have their centre inside which WDPA polygon (all statuses, so that
  #    "outside" means outside every WDPA site); 10-degree boxes keep memory small
  python analysis/pa_change_global.py assign --wdpa <GLOBAL_PAs_all_nonMARINE_2025.gpkg> \
      --cells output/gedi_pa_change_cells_combined.parquet --out output/cell_pa.parquet
  # 2. per-PA change and summaries (only designated/inscribed/established, no UNESCO-MAB)
  python analysis/pa_change_global.py summarize --wdpa <gpkg> --cell-pa output/cell_pa.parquet \
      --cells output/gedi_pa_change_cells_combined.parquet --out results/pa_change_global

Optional --bbox xmin ymin xmax ymax limits either step to a test region.
"""
import argparse
import os

import duckdb
import geopandas as gpd
import numpy as np
import pandas as pd
import pyogrio
import shapely
from pyproj import Transformer

RES, X0, Y0, NCOLS = 1000.89502334956, -17367530.4451615, 7314540.8306386, 34704
METRICS = ["rh98", "rh50", "pai", "fhd", "cover", "agbd"]
KEEP_STATUS = ("Designated", "Inscribed", "Established")
ATTRS = ["WDPAID", "NAME", "ISO3", "IUCN_CAT", "DESIG_ENG", "DESIG_TYPE", "GOV_TYPE", "STATUS",
         "STATUS_YR", "GIS_AREA", "PA_DEF"]


def cell_lonlat(cell_id):
    row, col = np.divmod(np.asarray(cell_id, dtype=np.int64), NCOLS)
    x, y = X0 + (col + 0.5) * RES, Y0 - (row + 0.5) * RES
    return Transformer.from_crs("EPSG:6933", "EPSG:4326", always_xy=True).transform(x, y)


def assign(a):
    con = duckdb.connect()
    con.sql("SET memory_limit = '8GB'")
    ids = con.sql(f"SELECT DISTINCT cell_id FROM read_parquet('{a.cells}')").fetchnumpy()["cell_id"]
    lon, lat = cell_lonlat(ids)
    print(f"{len(ids):,} distinct cells with GEDI data", flush=True)
    if a.bbox:
        m = (lon >= a.bbox[0]) & (lon < a.bbox[2]) & (lat >= a.bbox[1]) & (lat < a.bbox[3])
        ids, lon, lat = ids[m], lon[m], lat[m]
        boxes = [tuple(a.bbox)]
    else:
        boxes = [(x, y, x + a.box, y + a.box) for x in range(-180, 180, a.box) for y in range(-60, 60, a.box)]
    out = []
    for i, (x0, y0, x1, y1) in enumerate(boxes):
        m = (lon >= x0) & (lon < x1) & (lat >= y0) & (lat < y1)
        if not m.any():
            continue
        pas = pyogrio.read_dataframe(a.wdpa, bbox=(x0, y0, x1, y1), columns=["WDPAID"])
        if pas.empty:
            continue
        geoms = shapely.make_valid(pas.geometry.values)
        tree = shapely.STRtree(geoms)
        pts = shapely.points(lon[m], lat[m])
        pi, gi = tree.query(pts, predicate="within")
        if len(pi):
            out.append(pd.DataFrame({"cell_id": ids[m][pi], "WDPAID": pas.WDPAID.values[gi].astype("int64")}))
        print(f"  box {i + 1}/{len(boxes)} ({x0},{y0}): {m.sum():,} cells, {len(pas):,} polygons, "
              f"{len(np.unique(pi)):,} cells inside", flush=True)
    res = pd.concat(out, ignore_index=True).drop_duplicates() if out else pd.DataFrame(columns=["cell_id", "WDPAID"])
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    res.to_parquet(a.out, index=False)
    print(f"wrote {a.out}: {len(res):,} cell x WDPAID rows, {res.cell_id.nunique():,} cells inside any WDPA site")


def summarize(a):
    attrs = pyogrio.read_dataframe(a.wdpa, columns=ATTRS[1:], read_geometry=False)
    attrs["WDPAID"] = pyogrio.read_dataframe(a.wdpa, columns=["WDPAID"], read_geometry=False).WDPAID.astype("int64")
    kept = attrs.STATUS.isin(KEEP_STATUS) & (attrs.PA_DEF.astype(str) == "1") & \
        ~attrs.DESIG_ENG.eq("UNESCO-MAB Biosphere Reserve")
    site = (attrs[kept].groupby("WDPAID")
            .agg(NAME=("NAME", "first"), ISO3=("ISO3", "first"), IUCN_CAT=("IUCN_CAT", "first"),
                 DESIG_ENG=("DESIG_ENG", "first"), DESIG_TYPE=("DESIG_TYPE", "first"),
                 GOV_TYPE=("GOV_TYPE", "first"), STATUS_YR=("STATUS_YR", "min"), GIS_AREA=("GIS_AREA", "sum"))
            .reset_index())
    con = duckdb.connect()
    con.sql("SET memory_limit = '8GB'; SET preserve_insertion_order = false;")
    con.register("site", site)
    con.sql(f"CREATE VIEW cp AS SELECT * FROM read_parquet('{a.cell_pa}')")
    con.sql(f"CREATE VIEW c AS SELECT * FROM read_parquet('{a.cells}')")
    where = ""
    if a.bbox:
        where = (f"WHERE c.lon >= {a.bbox[0]} AND c.lon < {a.bbox[2]} "
                 f"AND c.lat >= {a.bbox[1]} AND c.lat < {a.bbox[3]}")
    os.makedirs(a.out, exist_ok=True)

    # per PA x epoch x metric: n cells, mean and sd of cell means
    aggs = ", ".join(f"count({m}_mean) AS n_{m}, avg({m}_mean) AS mean_{m}, stddev_samp({m}_mean) AS sd_{m}"
                     for m in METRICS)
    pe = con.sql(f"""SELECT cp.WDPAID, c.epoch, {aggs}
                     FROM c JOIN cp USING (cell_id) JOIN site USING (WDPAID) {where}
                     GROUP BY cp.WDPAID, c.epoch""").df()
    rows = []
    for m in METRICS:
        w = pe.pivot(index="WDPAID", columns="epoch", values=[f"n_{m}", f"mean_{m}", f"sd_{m}"])
        w.columns = [f"{s.split('_')[0]}_{e}" for s, e in w.columns]
        w = w.reindex(columns=["n_T1", "n_T2", "mean_T1", "mean_T2", "sd_T1", "sd_T2"])
        w["n_T1"], w["n_T2"] = w.n_T1.fillna(0), w.n_T2.fillna(0)
        w["change"] = w.mean_T2 - w.mean_T1
        w["change_se"] = np.sqrt(w.sd_T1 ** 2 / w.n_T1 + w.sd_T2 ** 2 / w.n_T2)  # cells independent: a lower bound
        w["qualifies"] = (w.n_T1 >= a.min_cells) & (w.n_T2 >= a.min_cells)
        w["metric"] = m
        rows.append(w.reset_index())
    per_pa = pd.concat(rows, ignore_index=True).merge(site, on="WDPAID", how="left")
    per_pa.to_csv(os.path.join(a.out, "pa_change_by_pa.csv"), index=False)

    def summary(by):
        out = []
        groups = [((), per_pa)] if by is None else per_pa.groupby(["metric", by], dropna=False)
        if by is None:
            groups = per_pa.groupby("metric")
        for key, d in groups:
            q = d[d.qualifies]
            rec = dict(zip(["metric"] + ([by] if by else []), key if isinstance(key, tuple) else (key,)))
            rec.update(n_pa_with_data=int(((d.n_T1 > 0) | (d.n_T2 > 0)).sum()), n_pa=len(q))
            if len(q):
                w = q.n_T1 + q.n_T2
                rec.update(cells_T1=int(q.n_T1.sum()), cells_T2=int(q.n_T2.sum()),
                           mean_T1=np.average(q.mean_T1, weights=w), mean_T2=np.average(q.mean_T2, weights=w),
                           change=np.average(q.change, weights=w), change_median=q.change.median(),
                           frac_pa_decline=(q.change < 0).mean())
            out.append(rec)
        return pd.DataFrame(out)

    s_global = summary(None)
    s_global.to_csv(os.path.join(a.out, "summary_global.csv"), index=False)
    for by, fn in [("ISO3", "summary_by_country.csv"), ("IUCN_CAT", "summary_by_iucn.csv"),
                   ("DESIG_TYPE", "summary_by_designation_type.csv"), ("GOV_TYPE", "summary_by_governance.csv")]:
        summary(by).to_csv(os.path.join(a.out, fn), index=False)

    # union of PA land vs outside every WDPA site, by GEDI region x PFT (each cell once)
    avgs = ", ".join(f"count({m}_mean) AS n_{m}, avg({m}_mean) AS mean_{m}" for m in METRICS)
    strata = con.sql(f"""
        WITH inside_kept AS (SELECT DISTINCT cell_id FROM cp JOIN site USING (WDPAID)),
             inside_any  AS (SELECT DISTINCT cell_id FROM cp)
        SELECT CASE WHEN k.cell_id IS NOT NULL THEN 'PA'
                    WHEN x.cell_id IS NOT NULL THEN 'other WDPA'
                    ELSE 'outside' END AS land, c.region_class, c.pft_class, c.epoch, {avgs}
        FROM c LEFT JOIN inside_kept k USING (cell_id) LEFT JOIN inside_any x USING (cell_id) {where}
        GROUP BY ALL""").df()
    strata.to_csv(os.path.join(a.out, "summary_region_pft_inside_outside.csv"), index=False)

    with pd.option_context("display.width", 220, "display.float_format", "{:.3f}".format):
        print("\nGlobal summary (PAs with >= %d cells in both epochs; weights = PA cells)" % a.min_cells)
        print(s_global.to_string(index=False))
        g = strata.groupby(["land", "epoch"])[[f"n_agbd"]].sum().unstack("epoch")
        print("\nAGBD cells by land type and epoch:\n", g.to_string())
    print(f"\nwrote {a.out}/")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("assign", "summarize"):
        p = sub.add_parser(name)
        p.add_argument("--wdpa", required=True)
        p.add_argument("--cells", required=True, help="combined cell table (dps/combine_cells.py)")
        p.add_argument("--out", required=True)
        p.add_argument("--bbox", type=float, nargs=4, help="xmin ymin xmax ymax (lon/lat) test region")
        if name == "assign":
            p.add_argument("--box", type=int, default=10, help="box size in degrees")
        else:
            p.add_argument("--cell-pa", required=True)
            p.add_argument("--min-cells", type=int, default=10)
    a = ap.parse_args()
    assign(a) if a.cmd == "assign" else summarize(a)


if __name__ == "__main__":
    main()
