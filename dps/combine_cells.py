#!/usr/bin/env python
"""
Merge 1 km cells that straddle a 1-degree tile edge, so there is one row per cell_id x epoch.

  python dps/combine_cells.py                       # output/gedi_pa_change_cells.parquet
                                                    #  -> output/gedi_pa_change_cells_combined.parquet

collect.py keeps one row per tile x cell x epoch. A cell on a tile edge gets shots from both
tiles, so it has two rows, each with its own n / mean / sd. Here they are pooled exactly:
  n    = sum of n
  mean = n-weighted mean
  sd   = sample sd of the pooled shots, from sum((n_i - 1) sd_i^2 + n_i mean_i^2) - N mean^2
lon, lat, frac_night, sensitivity_mean are n_shots-weighted means. tile_id, pft_class,
region_class and hemisphere come from the piece with the most shots. For merged cells
n_granules is the sum (an upper bound: a granule crossing the edge is counted twice) and
n_beams is the max (a lower bound). n_tiles says how many tile rows went into each row.
Cells in a single tile are copied unchanged (n_tiles = 1).
"""
import argparse
import os

import duckdb

METRICS = ["rh98", "rh50", "pai", "fhd", "cover", "agbd"]


def pooled(m):
    n, mu, sd = f"{m}_n", f"{m}_mean", f"{m}_sd"
    N = f"sum({n})"
    M = f"sum(CASE WHEN {n} > 0 THEN {n} * {mu} END) / nullif({N}, 0)"
    ss = (f"sum(CASE WHEN {n} > 0 THEN ({n} - 1) * coalesce({sd} * {sd}, 0) + {n} * {mu} * {mu} END)"
          f" - {N} * ({M}) * ({M})")
    return (f"{N} AS {n}, {M} AS {mu}, "
            f"CAST(CASE WHEN {N} > 1 THEN sqrt(greatest({ss}, 0) / ({N} - 1)) END AS FLOAT) AS {sd}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", default="output/gedi_pa_change_cells.parquet")
    ap.add_argument("--out", default="output/gedi_pa_change_cells_combined.parquet")
    ap.add_argument("--memory", default="8GB", help="DuckDB memory limit (the Hub has 16 GB)")
    a = ap.parse_args()

    con = duckdb.connect()
    tmp = os.path.join(os.path.dirname(os.path.abspath(a.out)), "duckdb_tmp")
    con.sql(f"SET memory_limit = '{a.memory}'; SET preserve_insertion_order = false; "
            f"SET temp_directory = '{tmp}';")
    con.sql(f"CREATE VIEW c AS SELECT * FROM read_parquet('{a.cells}')")
    con.sql("CREATE TEMP TABLE dup AS SELECT cell_id, epoch, count(*) AS k "
            "FROM c GROUP BY cell_id, epoch HAVING count(*) > 1")
    nd, nr = con.sql("SELECT count(*), coalesce(sum(k), 0) FROM dup").fetchone()
    print(f"{nd:,} cell x epoch keys appear in more than one tile ({nr:,} rows)")

    cols = [r[0] for r in con.sql("DESCRIBE c").fetchall()]
    w = "sum(n_shots * {0}) / sum(n_shots) AS {0}"
    merged = ", ".join(
        ["cell_id", "epoch",
         "arg_max(tile_id, n_shots) AS tile_id",
         "sum(n_shots) AS n_shots", "sum(n_granules) AS n_granules",
         "arg_max(pft_class, n_shots) AS pft_class", "arg_max(region_class, n_shots) AS region_class",
         "arg_max(hemisphere, n_shots) AS hemisphere",
         w.format("lon"), w.format("lat"), w.format("frac_night"),
         "max(n_beams) AS n_beams", w.format("sensitivity_mean")]
        + [pooled(m) for m in METRICS]
        + ["count(*) AS n_tiles"])
    order = ", ".join(cols + ["n_tiles"])
    con.sql(f"""
      COPY (
        SELECT {order} FROM (
          SELECT c.*, 1 AS n_tiles FROM c ANTI JOIN dup USING (cell_id, epoch)
          UNION ALL BY NAME
          SELECT * FROM (SELECT {merged} FROM c SEMI JOIN dup USING (cell_id, epoch)
                         GROUP BY cell_id, epoch)
        )
      ) TO '{a.out}' (FORMAT parquet)""")

    # checks: one row per key, and every count preserved
    sums = ", ".join(["sum(n_shots)"] + [f"sum({m}_n)" for m in METRICS])
    before = con.sql(f"SELECT count(*), {sums} FROM c").fetchone()
    after = con.sql(f"SELECT count(*), {sums} FROM '{a.out}'").fetchone()
    keys = con.sql(f"SELECT count(*) FROM (SELECT DISTINCT cell_id, epoch FROM '{a.out}')").fetchone()[0]
    print(f"rows: {before[0]:,} -> {after[0]:,} (expected {before[0] - nr + nd:,}); distinct keys {keys:,}")
    assert after[0] == before[0] - nr + nd == keys, "row count does not match one row per cell x epoch"
    assert before[1:] == after[1:], f"counts changed: {before[1:]} vs {after[1:]}"
    print("n_shots and every metric's n are unchanged; wrote", a.out)


if __name__ == "__main__":
    main()
