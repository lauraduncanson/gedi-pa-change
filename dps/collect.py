#!/usr/bin/env python
"""
Combine all per-tile cell files into one table, and the per-tile summaries into one CSV.

  python dps/collect.py --out /home/jovyan/GEDI/gedi_pa_change_cells.parquet
"""
import argparse
import glob
import json
import os

import duckdb
import pandas as pd

import common as C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--done-dir", default=C.DPS_OUTPUT)
    ap.add_argument("--out", default="gedi_pa_change_cells.parquet")
    a = ap.parse_args()

    cells = glob.glob(os.path.join(a.done_dir, "**", "cells_*.parquet"), recursive=True)
    sums = glob.glob(os.path.join(a.done_dir, "**", "summary_*.json"), recursive=True)
    print(f"{len(sums)} tiles finished; {len(cells)} have cell files")
    s = pd.DataFrame([json.load(open(f)) for f in sums]).drop_duplicates("tile", keep="last")
    s.to_csv(os.path.splitext(a.out)[0] + "_tile_summary.csv", index=False)
    print(f"shots kept: T1 {s.n_T1.sum():,}  T2 {s.n_T2.sum():,}")
    if cells:
        flist = ", ".join(f"'{f}'" for f in cells)
        # a tile re-run leaves two copies; keep one row per tile x cell x epoch
        duckdb.sql(f"""COPY (SELECT DISTINCT ON (tile_id, cell_id, epoch) *
                             FROM read_parquet([{flist}], union_by_name=true))
                       TO '{a.out}' (FORMAT parquet)""")
        print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
