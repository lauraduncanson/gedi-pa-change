#!/usr/bin/env python
"""
Global matching steps 1-3 from global layers only (design in CLAUDE.md, "Global matching design").

For each country: 1 km EASE2 cells with GEDI data (centre in the country) are split by WWF biome
into matching units <ISO3>_<biome>. In each unit:
  treatment = cells fully inside a kept PA (centre and 4 corners inside the polygon); kept PAs as
              in the global PA change analysis (Designated/Inscribed/Established, no UNESCO-MAB).
              A cell inside several kept PAs is a treatment cell of each.
  control   = cells whose centre is more than --buffer-km (default 10) from every WDPA polygon (any status).
Covariates (2020 set) are read at the cell centre from the global rasters; land cover (lc2020),
biome and ecoregion are written as factors with the same levels in both tables so that her
step 4 stratifies on them. Output per unit, in the format of matching/r/brazil_v3 (v4):
  <out>/<unit>/<unit>_prepped_PAs.RDS      (one data frame, pa_id column)
  <out>/<unit>/<unit>_prepped_control.RDS
and a manifest <out>/manifest_<ISO3>.csv (unit, PAs, PA cells, controls).

  python matching/prep_country_inputs.py --iso3 LBR CRI --inputs <folder with cov/, countries/, wdpa> \
      --cells output/gedi_pa_change_cells_combined.parquet --out s3://.../matching_inputs_global
"""
import argparse
import os
import tempfile
import time

import boto3
import duckdb
import geopandas as gpd
import numpy as np
import pandas as pd
import pyogrio
import pyreadr
import rasterio
import shapely
from pyproj import Transformer
from rasterio.windows import Window

RES, X0, Y0, NCOLS = 1000.89502334956, -17367530.4451615, 7314540.8306386, 34704
HALF = RES / 2 * 0.999  # corners just inside the cell edge
KEEP_STATUS = ("Designated", "Inscribed", "Established")
TILE_M = 500_000  # tile size for processing (EASE2 metres)
COV = {  # column name in her step 4 -> raster file (INPUT_covariates_2020)
    "d2road": "d2roads", "d2city": "dcities", "elevation": "dem", "slope": "slope",
    "popcnt": "pop_cnt_2020", "popden": "pop_den_2020", "tt2city": "tt2cities_2015",
    "prec": "wc_prec_2010-2018", "mean_temp": "wc_tavg_2010-2018", "max_temp": "wc_tmax_2010-2018",
    "min_temp": "wc_tmin_2010-2018", "gedi_l4b": "gedi_l4b",
}
LC_LABELS = ["l1_forest", "l2_grassland", "l3_agriculture", "l4_wetlands", "l5_artificial",
             "l6_other land/bare", "l7_water"]
TO_LONLAT = Transformer.from_crs("EPSG:6933", "EPSG:4326", always_xy=True)


def all_cell_ids(cells, cache):
    if os.path.exists(cache):
        return np.load(cache)
    ids = duckdb.sql(f"SELECT DISTINCT cell_id FROM read_parquet('{cells}') ORDER BY cell_id").fetchnumpy()["cell_id"]
    np.save(cache, ids)
    return ids


def xy(ids):
    row, col = np.divmod(ids, NCOLS)
    return X0 + (col + 0.5) * RES, Y0 - (row + 0.5) * RES


def sample(path, lon, lat):
    """Value at each lon/lat, reading one window that covers all points (EPSG:4326 rasters)."""
    with rasterio.open(path) as r:
        t = r.transform
        cols = np.floor((lon - t.c) / t.a).astype(np.int64)
        rows = np.floor((lat - t.f) / t.e).astype(np.int64)
        ok = (cols >= 0) & (cols < r.width) & (rows >= 0) & (rows < r.height)
        out = np.full(len(lon), np.nan)
        if ok.any():
            c0, c1, r0, r1 = cols[ok].min(), cols[ok].max(), rows[ok].min(), rows[ok].max()
            a = r.read(1, window=Window(c0, r0, c1 - c0 + 1, r1 - r0 + 1)).astype("float64")
            if r.nodata is not None and not np.isnan(r.nodata):
                a[a == r.nodata] = np.nan
            out[ok] = a[rows[ok] - r0, cols[ok] - c0]
    return out


def write(df, dest, s3):
    with tempfile.TemporaryDirectory() as tmp:
        local = os.path.join(tmp, os.path.basename(dest))
        pyreadr.write_rds(local, df, compress="gzip")
        if dest.startswith("s3://"):
            b, k = dest[5:].split("/", 1)
            s3.upload_file(local, b, k)
        else:
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            os.replace(local, dest)


def country(iso3, a, ids_all, xs_all, ys_all, s3):
    t0 = time.time()
    shp = os.path.join(a.inputs, "countries", f"{iso3}.shp")
    adm = gpd.read_file(shp)
    if adm.crs is None:
        adm = adm.set_crs("EPSG:4326")
    adm = adm.to_crs("EPSG:6933")
    geom = shapely.make_valid(shapely.union_all(adm.geometry.values))
    xmin, ymin, xmax, ymax = geom.bounds
    m = (xs_all >= xmin) & (xs_all <= xmax) & (ys_all >= ymin) & (ys_all <= ymax)
    ids, xs, ys = ids_all[m], xs_all[m], ys_all[m]
    inside = shapely.contains_xy(geom, xs, ys) if len(ids) else np.zeros(0, bool)
    ids, xs, ys = ids[inside], xs[inside], ys[inside]
    if len(ids) == 0:
        print(f"{iso3}: no GEDI cells"); return []

    # Process the country in ~500 km tiles: each tile reads only the WDPA polygons within the
    # buffer of it, so memory is bounded by a tile, not a country (the USA in one block exceeded
    # the Hub's 16 GB), and countries crossing the antimeridian only touch tiles with cells.
    # A cell's status depends only on polygons within the buffer of it, so results are identical.
    buf = a.buffer_km * 1000.0
    tx = np.floor((xs - X0) / TILE_M).astype(np.int64)
    ty = np.floor((Y0 - ys) / TILE_M).astype(np.int64)
    tile_key = tx * 100000 + ty
    near = np.zeros(len(ids), bool)
    tr_cells, tr_pas, n_poly = [], [], 0
    gpkg = os.path.join(a.inputs, "wdpa2025_ease2.gpkg")
    for key in np.unique(tile_key):
        idx = np.where(tile_key == key)[0]
        x_lo, x_hi = xs[idx].min() - HALF, xs[idx].max() + HALF
        y_lo, y_hi = ys[idx].min() - HALF, ys[idx].max() + HALF
        wdpa = pyogrio.read_dataframe(gpkg, bbox=(x_lo - buf, y_lo - buf, x_hi + buf, y_hi + buf),
                                      columns=["WDPAID", "STATUS", "DESIG_ENG", "PA_DEF"])
        wdpa = wdpa[~wdpa.geometry.isna()].reset_index(drop=True)
        if wdpa.empty:
            continue
        n_poly += len(wdpa)
        g_all = shapely.make_valid(wdpa.geometry.values)
        pts = shapely.points(xs[idx], ys[idx])
        # controls: centre more than buffer_km from every WDPA polygon (any status)
        pi, _ = shapely.STRtree(g_all).query(pts, predicate="dwithin", distance=buf)
        near[idx[np.unique(pi)]] = True
        # treatment: centre and four corners inside a kept PA
        kept = (wdpa.STATUS.isin(KEEP_STATUS) & (wdpa.PA_DEF.astype(str) == "1")
                & ~wdpa.DESIG_ENG.eq("UNESCO-MAB Biosphere Reserve")).values
        if kept.any():
            gk = g_all[kept]; idk = wdpa.WDPAID.values[kept].astype(np.int64)
            pi, gi = shapely.STRtree(gk).query(pts, predicate="within")
            full = np.ones(len(pi), bool)
            for dx, dy in ((-HALF, -HALF), (-HALF, HALF), (HALF, -HALF), (HALF, HALF)):
                full &= shapely.contains_xy(gk[gi], xs[idx[pi]] + dx, ys[idx[pi]] + dy)
            tr_cells.append(idx[pi[full]]); tr_pas.append(idk[gi[full]])
        del wdpa, g_all, pts
    tr_cell = np.concatenate(tr_cells) if tr_cells else np.zeros(0, np.int64)
    tr_pa = np.concatenate(tr_pas) if tr_pas else np.zeros(0, np.int64)
    pair = pd.DataFrame({"i": tr_cell, "pa": tr_pa}).drop_duplicates()  # parcels of one WDPAID
    tr_cell, tr_pa = pair.i.values, pair.pa.values

    lon, lat = TO_LONLAT.transform(xs, ys)
    cov = {k: sample(os.path.join(a.inputs, "cov", f"{v}.tif"), lon, lat) for k, v in COV.items()}
    lc = sample(os.path.join(a.inputs, "cov", "lc2020.tif"), lon, lat)
    biome = sample(os.path.join(a.inputs, "cov", "wwf_biomes.tif"), lon, lat)
    ecoreg = sample(os.path.join(a.inputs, "cov", "wwf_ecoreg.tif"), lon, lat)
    base = pd.DataFrame({"lat": lat, "lon": lon, **cov})
    base["land_cover"] = pd.Categorical([LC_LABELS[int(v) - 1] if np.isfinite(v) and 1 <= v <= 7 else None for v in lc],
                                        categories=LC_LABELS)
    base["biome"] = biome
    base["ecoreg"] = ecoreg
    base.loc[base.gedi_l4b < 0, "gedi_l4b"] = np.nan

    rows = []
    s3_out = a.out.rstrip("/")
    for b in sorted(pd.unique(biome[np.isfinite(biome)])):
        unit = f"{iso3}_{int(b):02d}"
        in_b = np.isfinite(biome) & (biome == b)
        ctl_idx = np.where(in_b & ~near)[0]
        sel = in_b[tr_cell]
        t_idx, t_pa = tr_cell[sel], tr_pa[sel]
        if len(t_idx) == 0 or len(ctl_idx) == 0:
            rows.append(dict(unit=unit, iso3=iso3, biome=int(b), n_pa=len(np.unique(t_pa)), n_pa_cells=len(t_idx),
                             n_controls=len(ctl_idx), written=False))
            continue
        eco_levels = [f"e{int(v)}" for v in sorted(pd.unique(np.concatenate([ecoreg[ctl_idx], ecoreg[t_idx]])))
                      if np.isfinite(v)]

        def table(idx, status, pa=None):
            d = base.iloc[idx].reset_index(drop=True)
            d["wwfbiom"] = pd.Categorical([f"b{int(b):02d}"] * len(d), categories=[f"b{int(b):02d}"])
            d["wwfecoreg"] = pd.Categorical([f"e{int(v)}" if np.isfinite(v) else None for v in d.ecoreg],
                                            categories=eco_levels)
            d["status"] = status
            d["cell_id"] = ids[idx]
            if pa is None:
                d["UID"] = np.arange(1, len(d) + 1, dtype=np.int32)
            else:
                d["pa_id"] = pa.astype(str)
                d["UID"] = (d.groupby("pa_id").cumcount() + 1).astype(np.int32)
            cols = ["lat", "lon", "UID", "status", "land_cover", "wwfbiom", "wwfecoreg", "elevation", "slope",
                    "mean_temp", "max_temp", "min_temp", "prec", "d2road", "d2city", "tt2city", "popden",
                    "popcnt", "gedi_l4b", "cell_id"] + ([] if pa is None else ["pa_id"])
            return d[cols]

        ctl, pas = table(ctl_idx, False), table(t_idx, True, t_pa)
        write(ctl, f"{s3_out}/{unit}/{unit}_prepped_control.RDS", s3)
        write(pas, f"{s3_out}/{unit}/{unit}_prepped_PAs.RDS", s3)
        rows.append(dict(unit=unit, iso3=iso3, biome=int(b), n_pa=pas.pa_id.nunique(), n_pa_cells=len(pas),
                         n_controls=len(ctl), pa_cells_complete=int(pas.drop(columns="cell_id").notna().all(axis=1).sum()),
                         controls_complete=int(ctl.drop(columns="cell_id").notna().all(axis=1).sum()), written=True))
    man = pd.DataFrame(rows)
    print(f"{iso3}: {len(ids):,} GEDI cells, {n_poly:,} WDPA polygons read, {len(np.unique(tr_cell)):,} treatment cells "
          f"in {len(np.unique(tr_pa)):,} PAs, {int((~near).sum()):,} control cells, {int(man.written.sum()) if len(man) else 0} units "
          f"({time.time() - t0:.0f} s)", flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iso3", nargs="+", required=True, help="ISO3 codes, or ALL (every country shapefile)")
    ap.add_argument("--inputs", required=True, help="folder with cov/*.tif, countries/*.shp, wdpa2025_ease2.gpkg")
    ap.add_argument("--cells", required=True, help="combined cell table")
    ap.add_argument("--out", required=True, help="s3:// prefix or local folder for the unit inputs")
    ap.add_argument("--manifest", default="output/matching_units_manifest.csv")
    ap.add_argument("--buffer-km", type=float, default=10.0, help="controls: centre further than this from any WDPA polygon")
    a = ap.parse_args()
    isos = a.iso3
    if isos == ["ALL"]:
        isos = sorted(f[:-4] for f in os.listdir(os.path.join(a.inputs, "countries")) if f.endswith(".shp") and len(f) == 7)
    ids = all_cell_ids(a.cells, os.path.join(a.inputs, "cell_ids.npy"))
    xs, ys = xy(ids)
    s3 = boto3.client("s3")
    rows = []
    for iso3 in isos:
        try:
            rows += country(iso3, a, ids, xs, ys, s3)
        except Exception as e:
            print(f"{iso3}: FAILED {type(e).__name__}: {e}", flush=True)
    man = pd.DataFrame(rows)
    if os.path.exists(a.manifest):
        old = pd.read_csv(a.manifest)
        man = pd.concat([old[~old.iso3.isin(man.iso3.unique())] if len(man) else old, man], ignore_index=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.manifest)), exist_ok=True)
    man.to_csv(a.manifest, index=False)
    print(f"manifest: {a.manifest} ({len(man)} units)")


if __name__ == "__main__":
    main()
