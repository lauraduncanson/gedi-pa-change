#!/usr/bin/env python
"""
Matching step 1, rebuilt on the GEDI cell table.

For each country, list the 1 km EASE2 cells that (a) lie inside the country boundary and
(b) have at least MIN_SHOTS quality-filtered, power-beam, growing-season GEDI shots in either
epoch. Writes <ISO3>_grid_cells.csv with the cell-centre coordinates in the column names the
R matching code expects (x.overlap = lon, y.overlap = lat), plus cell_id for joining later.

This replaces the tile loop in step 1 of vl_GEDI-PA_2024/global_process_part1_2024_MAAP_step123.R
(see matching/README.md for the R patch). Steps 2-4 are unchanged.

Usage (Hub):
  python matching/make_country_grid.py --cells /home/jovyan/GEDI/gedi_pa_change_cells.parquet \
      --out /home/jovyan/my-public-bucket/gedi_pa_change/country_grids --iso3 BRA
  ... --iso3 ALL        # every country with a boundary shapefile
"""
import argparse
import glob
import os
import shutil
import tempfile

import duckdb
import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import Transformer

# GEDI 1 km EASE-Grid 2.0 (EPSG:6933); identical to extract_tile.py
RES, X0, Y0, NCOLS = 1000.89502334956, -17367530.4451615, 7314540.8306386, 34704
COUNTRY_SHP = "/home/jovyan/shared-buckets/leitoldv/GEDI_global_PA_v2/WDPA_countries/shp/{iso3}.shp"


def cell_centres(cell_id):
    row, col = np.divmod(np.asarray(cell_id, dtype=np.int64), NCOLS)
    x = X0 + (col + 0.5) * RES
    y = Y0 - (row + 0.5) * RES
    lon, lat = Transformer.from_crs("EPSG:6933", "EPSG:4326", always_xy=True).transform(x, y)
    return lon, lat


def country_grid(iso3, cells, min_shots=1, shp_pattern=COUNTRY_SHP):
    adm = gpd.read_file(shp_pattern.format(iso3=iso3)).to_crs("EPSG:4326")
    xmin, ymin, xmax, ymax = adm.total_bounds
    pad = 0.02  # cell mean shot location can sit up to ~0.01 deg from the centre
    ids = duckdb.sql(f"""
        SELECT cell_id, list(DISTINCT epoch ORDER BY epoch) AS epochs, sum(n_shots) AS n_shots
        FROM read_parquet('{cells}')
        WHERE n_shots >= {min_shots}
          AND lon BETWEEN {xmin - pad} AND {xmax + pad}
          AND lat BETWEEN {ymin - pad} AND {ymax + pad}
        GROUP BY cell_id
    """).df()
    if ids.empty:
        return pd.DataFrame(columns=["x.overlap", "y.overlap", "cell_id", "epochs", "n_shots"])
    lon, lat = cell_centres(ids.cell_id)
    pts = gpd.GeoDataFrame(ids, geometry=gpd.points_from_xy(lon, lat), crs="EPSG:4326")
    # cell centre inside the boundary, as terra::mask does in the original step 1
    inside = gpd.sjoin(pts, adm[["geometry"]], predicate="within", how="inner")
    inside = inside[~inside.index.duplicated()]
    out = pd.DataFrame({"x.overlap": inside.geometry.x, "y.overlap": inside.geometry.y,
                        "cell_id": inside.cell_id,
                        "epochs": inside.epochs.map(lambda e: "+".join(e)),
                        "n_shots": inside.n_shots})
    return out.sort_values("cell_id").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True, help="merged cell table from dps/collect.py")
    ap.add_argument("--out", required=True, help="output folder (a bucket path is fine)")
    ap.add_argument("--iso3", nargs="+", required=True, help="country codes, or ALL")
    ap.add_argument("--min-shots", type=int, default=1, help="minimum shots per cell (either epoch)")
    ap.add_argument("--shp", default=COUNTRY_SHP, help="country shapefile pattern with {iso3}")
    a = ap.parse_args()

    isos = a.iso3
    if isos == ["ALL"]:
        isos = sorted(os.path.basename(f)[:-4] for f in glob.glob(a.shp.format(iso3="*")))
    os.makedirs(a.out, exist_ok=True)
    for iso3 in isos:
        try:
            g = country_grid(iso3, a.cells, a.min_shots, a.shp)
        except Exception as err:
            print(f"{iso3}: FAILED ({err})")
            continue
        tmp = os.path.join(tempfile.gettempdir(), f"{iso3}_grid_cells.csv")
        g.to_csv(tmp, index=False)
        shutil.copy(tmp, os.path.join(a.out, f"{iso3}_grid_cells.csv"))
        counts = g.epochs.value_counts().to_dict() if len(g) else {}
        print(f"{iso3}: {len(g):,} cells {counts}")


if __name__ == "__main__":
    main()
