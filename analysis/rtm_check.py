#!/usr/bin/env python
"""
Regression-to-the-mean check for matching on GEDI L4B baseline biomass.

Idea: the matched cells were selected on their L4B value (2019-2023 GEDI). T1 cell means share
shots with L4B, so T1 tracks L4B by construction. T2 shots are independent of the matching, so
T2 - L4B shows whether cells drift back toward their population mean. In undisturbed cells with
no real change:
  * T2 - L4B falls with baseline in both groups (L4B error regressing out); that alone is expected.
  * The test is the gap between groups at the same baseline. Controls drawn from a lower-biomass
    population carry more positive L4B error, so with regression to the mean their T2 - L4B sits
    below that of protected cells, mimicking a protection effect. A gap near 0 means no sign of it;
    a common offset in both groups is instrument change, not RTM.

Outputs <ISO3>_rtm_bins.csv, <ISO3>_rtm.png and a printed control-minus-protected gap.

  python analysis/rtm_check.py --iso3 BRA --matches <folder with BRA matching RDS> \
      --cells output/gedi_pa_change_cells_combined.parquet \
      --l4b /home/jovyan/shared-buckets/leitoldv/GEDI_global_change/INPUT_covariates_2020/gedi_l4b.tif \
      --disturbance "/home/jovyan/shared-buckets/leitoldv/GEDI_global_change/INPUT_covariates_2020/BRA_*_frac_1km.tif" \
      --out /home/jovyan/GEDI/results
"""
import argparse
import glob
import os

import duckdb
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import rasterio  # noqa: E402
from pyproj import Transformer  # noqa: E402

from did_country import RES, X0, Y0, NCOLS, add_cell_id, read_matches  # noqa: E402

METRIC = "agbd"


def cell_xy(cell_id):
    row, col = np.divmod(np.asarray(cell_id, dtype=np.int64), NCOLS)
    return X0 + (col + 0.5) * RES, Y0 - (row + 0.5) * RES


def sample(path, cell_id):
    x, y = cell_xy(cell_id)
    with rasterio.open(path) as src:
        if src.crs and src.crs.to_epsg() != 6933:
            x, y = Transformer.from_crs("EPSG:6933", src.crs, always_xy=True).transform(x, y)
        vals = np.array([v[0] for v in src.sample(zip(x, y))], dtype="float64")
        if src.nodata is not None:
            vals[vals == src.nodata] = np.nan
    return vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iso3", required=True)
    ap.add_argument("--matches", required=True)
    ap.add_argument("--cells", required=True)
    ap.add_argument("--l4b", required=True, help="GEDI L4B AGBD raster used in the matching")
    ap.add_argument("--disturbance", default=None,
                    help="glob of disturbance-fraction rasters; a cell counts as undisturbed when "
                         "all are 0. Omit to use all cells.")
    ap.add_argument("--bin-width", type=float, default=25.0, help="baseline bin width (Mg/ha)")
    ap.add_argument("--min-cells", type=int, default=20, help="min T2 cells per group in a bin")
    ap.add_argument("--out", default="results")
    a = ap.parse_args()

    m = add_cell_id(read_matches(a.matches, a.iso3))
    cells = m[["cell_id", "status"]].drop_duplicates()
    cells["l4b"] = sample(a.l4b, cells.cell_id)
    cells.loc[cells.l4b < 0, "l4b"] = np.nan

    label = "all matched cells"
    if a.disturbance:
        paths = sorted(glob.glob(a.disturbance))
        if not paths:
            raise SystemExit(f"no disturbance rasters match {a.disturbance}")
        dist = np.zeros(len(cells))
        for p in paths:
            v = np.nan_to_num(sample(p, cells.cell_id), nan=0.0)
            dist += np.clip(v, 0, None)
            print(f"  disturbance layer: {os.path.basename(p)}")
        cells = cells[dist == 0]
        label = "undisturbed matched cells"

    con = duckdb.connect()
    con.register("c", cells)
    d = con.sql(f"""SELECT c.*, t.epoch, t.{METRIC}_mean AS v
                    FROM c JOIN read_parquet('{a.cells}') t USING (cell_id)
                    WHERE t.{METRIC}_mean IS NOT NULL""").df().dropna(subset=["l4b"])
    d["group"] = np.where(d.status, "protected", "control")
    d["bin"] = (d.l4b // a.bin_width) * a.bin_width
    print(f"{a.iso3}: {label}: {d.cell_id.nunique():,} cells with L4B and GEDI data")

    rows = []
    for (grp, b), g in d.groupby(["group", "bin"]):
        t1, t2 = g[g.epoch == "T1"], g[g.epoch == "T2"]
        rows.append(dict(group=grp, bin=b, n_t1=len(t1), n_t2=len(t2), l4b=g.l4b.mean(),
                         t1_minus_l4b=(t1.v - t1.l4b).mean() if len(t1) else np.nan,
                         t2_minus_l4b=(t2.v - t2.l4b).mean() if len(t2) else np.nan,
                         t2_se=(t2.v - t2.l4b).std() / np.sqrt(len(t2)) if len(t2) > 1 else np.nan))
    b = pd.DataFrame(rows).sort_values(["group", "bin"])
    os.makedirs(a.out, exist_ok=True)
    b.to_csv(os.path.join(a.out, f"{a.iso3}_rtm_bins.csv"), index=False)

    # Within a group, (T2 - L4B) always falls with L4B (L4B error regresses out), so slopes alone
    # are not the test. The test is the gap between groups at the same baseline: controls drawn
    # from a lower-biomass population carry more positive L4B error, so their T2 - L4B sits lower.
    t2 = d[d.epoch == "T2"].assign(dev=lambda x: x.v - x.l4b)
    per_bin = t2.groupby(["bin", "group"]).dev.agg(["mean", "count"]).unstack("group")
    ok = per_bin.dropna()
    ok = ok[(ok[("count", "protected")] >= a.min_cells) & (ok[("count", "control")] >= a.min_cells)]
    if ok.empty:
        print("\nnot enough T2 cells per bin in both groups to compare")
    else:
        gap = ok[("mean", "control")] - ok[("mean", "protected")]
        w = ok[("count", "protected")] + ok[("count", "control")]
        rtm = np.average(gap, weights=w)
        print(f"\nT2 - L4B, control minus protected, at the same baseline ({len(ok)} bins): {rtm:+.2f} Mg/ha")
        print("  ~0  -> no sign of regression to the mean in this comparison")
        print("  < 0 -> controls fall back further than protected cells: part of a positive")
        print("         protection effect (DiD) of about this size may be regression to the mean")
        b = b.merge(gap.rename("gap_ctl_minus_pa").reset_index(), on="bin", how="left")
        b.to_csv(os.path.join(a.out, f"{a.iso3}_rtm_bins.csv"), index=False)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for grp, color in (("protected", "#2e7d32"), ("control", "#b5651d")):
        g = b[(b.group == grp) & (b.n_t2 >= a.min_cells)]
        ax.errorbar(g.l4b, g.t2_minus_l4b, yerr=1.96 * g.t2_se, fmt="o-", color=color,
                    capsize=3, label=f"{grp}: T2 - L4B")
        g1 = b[(b.group == grp) & (b.n_t1 >= a.min_cells)]
        ax.plot(g1.l4b, g1.t1_minus_l4b, ":", color=color, alpha=0.7, label=f"{grp}: T1 - L4B")
    ax.axhline(0, color="grey", lw=0.8)
    ax.set_xlabel("Baseline GEDI L4B AGBD (Mg/ha), binned")
    ax.set_ylabel("Cell mean AGBD minus L4B (Mg/ha)")
    ax.set_title(f"{a.iso3}: regression-to-the-mean check ({label})", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    png = os.path.join(a.out, f"{a.iso3}_rtm.png")
    fig.savefig(png, dpi=150)
    print(f"\nwrote {png} and {a.iso3}_rtm_bins.csv")


if __name__ == "__main__":
    main()
