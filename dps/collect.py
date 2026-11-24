#!/usr/bin/env python
"""
Combine all per-tile cell files into one table, and the per-tile summaries into one CSV.

  python dps/collect.py --out output/gedi_pa_change_cells.parquet

Reads DPS outputs straight from S3 by default (see common.DPS_OUTPUT); --done-dir also takes a
local folder or glob pattern. Write --out to a local disk, not a bucket mount.
"""
import argparse
import os
import sys

import pandas as pd

import common as C

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "extract_tile"))
from extract_tile import connect  # noqa: E402  (DuckDB with the Hub AWS role for s3:// reads)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--done-dir", default=C.DPS_OUTPUT)
    ap.add_argument("--out", default="output/gedi_pa_change_cells.parquet")
    a = ap.parse_args()

    cells = C.output_files(a.done_dir, "cells_", ".parquet")
    sums = C.output_files(a.done_dir, "summary_", ".json")
    print(f"{len(sums)} summary files; {len(cells)} cell files")
    con = connect()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

    slist = ", ".join(f"'{f}'" for f in sums)
    s = con.sql(f"SELECT * FROM read_json([{slist}], union_by_name=true, filename=true)").df()
    s = s.sort_values("filename").drop_duplicates("tile", keep="last").drop(columns="filename")
    s.to_csv(os.path.splitext(a.out)[0] + "_tile_summary.csv", index=False)
    print(f"{len(s)} tiles; shots kept: T1 {s.n_T1.sum():,}  T2 {s.n_T2.sum():,}")

    if cells:
        # A re-run tile leaves two cell files; keep the newest (DPS paths contain the date, so the
        # last in sorted order). Choosing files up front lets DuckDB stream the copy: a global
        # DISTINCT ON over all rows ran out of the Hub's 16 GB memory limit.
        latest = {}
        for f in sorted(cells):
            latest[os.path.basename(f)[len("cells_"):-len(".parquet")]] = f
        print(f"{len(latest)} tiles with cell files ({len(cells) - len(latest)} older duplicates skipped)")
        con.sql("SET memory_limit = '8GB'; SET preserve_insertion_order = false;")
        flist = ", ".join(f"'{f}'" for f in latest.values())
        con.sql(f"""COPY (SELECT * FROM read_parquet([{flist}], union_by_name=true))
                    TO '{a.out}' (FORMAT parquet)""")
        n = con.sql(f"SELECT count(*), count(DISTINCT tile_id) FROM '{a.out}'").fetchone()
        print(f"wrote {a.out}: {n[0]:,} cell x epoch rows from {n[1]:,} tiles")


if __name__ == "__main__":
    main()
