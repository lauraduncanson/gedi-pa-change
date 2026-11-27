#!/usr/bin/env python
"""
Area-weighted DiD means and total biomass DiD for protected areas, measured and imputed.
Design: CLAUDE.md, "PAs without a valid counterfactual" and "Totals".

  # 1. per-PA attributes: country, majority biome and GEDI region, dissolved area (slow part)
  python analysis/totals.py attrs --wdpa <wdpa2025_ease2.gpkg> --biomes <wwf_biomes.tif> \
      --regions <GEDI_ANCI_CONTINENT_r1000m_EASE2.0_UMD_v5.tif> --cell-pa output/cell_pa.parquet \
      --out output/totals_pa_attrs.parquet
  # 2. measured / imputed DiD, means and totals
  python analysis/totals.py run --attrs output/totals_pa_attrs.parquet --did results/global_did \
      --out results/totals                      # headline, |SMD| <= 0.25
  python analysis/totals.py run --smd-max 0.1 --out results/totals_smd0.1   # balance sensitivity (also 0.5)

Universe: kept PAs (Designated/Inscribed/Established, PA_DEF 1, no UNESCO-MAB) with at least one
GEDI cell (centre inside). A PA is measured for a metric if its DiD qualifies (>= 10 cells in each
of PA/control x T1/T2) and |T1 AGBD SMD| <= 0.25; otherwise its DiD is imputed from the
area-weighted mean DiD of >= 5 measured PAs in its biome x country, else GEDI region x biome,
else biome; else it is excluded.

Area: parcels of one WDPAID dissolved, clipped to 51.6 S-51.6 N, in EASE2 (equal area). Within a
country, PA areas are scaled so they sum to the dissolved union of the country's universe PAs
(overlaps counted once).

Outputs in --out:
  pa_table.parquet          one row per PA x metric: attributes, area, DiD used, source
                            (measured / imputed:<level> / excluded)
  means_<by>.csv            area-weighted means per metric (for figures): DiD measured-only and
                            measured + imputed with 95% bootstrap CIs, PA and control change
                            (measured PAs), T1 balance, PA counts and area shares
  totals_agbd_<by>.csv      AGBD DiD x area in Tg, measured-only and measured + imputed, with CIs
  <by> = global, country, biome, region.
Bootstrap: resample measured PAs (per metric); imputed values are recomputed from each resample.
"""
import argparse
import os
import tempfile
import time

import duckdb
import geopandas as gpd
import numpy as np
import pandas as pd
import pyogrio
import rasterio
import shapely
from pyproj import Transformer

RES, X0, Y0, NCOLS = 1000.89502334956, -17367530.4451615, 7314540.8306386, 34704
KEEP_STATUS = ("Designated", "Inscribed", "Established")
METRICS = ["agbd", "rh98", "rh50", "pai", "fhd", "cover"]
REGIONS = {1: "Europe", 2: "North Asia", 3: "Australasia", 4: "Africa", 5: "South Asia",
           6: "South America", 7: "North America"}  # GEDI L4A prediction strata
SMD_MAX, MIN_STRATUM, N_BOOT = 0.25, 5, 1000
LEVELS = [("biome_country", ["biome", "iso3"]), ("region_biome", ["region", "biome"]), ("biome", ["biome"])]
BATCH = 20_000


# ---------------------------------------------------------------- attributes
def kept_attrs(wdpa):
    a = pyogrio.read_dataframe(wdpa, read_geometry=False, columns=["WDPAID", "STATUS", "PA_DEF", "DESIG_ENG", "ISO3"])
    k = a.STATUS.isin(KEEP_STATUS) & (a.PA_DEF.astype(str) == "1") & ~a.DESIG_ENG.eq("UNESCO-MAB Biosphere Reserve")
    return k.values, a


def majority(df, col):
    d = df.dropna(subset=[col]).groupby(["WDPAID", col]).size().reset_index(name="n")
    return d.sort_values(["WDPAID", "n", col], ascending=[True, False, True]).drop_duplicates("WDPAID").set_index("WDPAID")[col]


def attrs(a):
    t0 = time.time()
    kept, w = kept_attrs(a.wdpa)
    kept_ids = pd.DataFrame({"WDPAID": np.unique(w.WDPAID.values[kept])})
    con = duckdb.connect()
    con.register("k", kept_ids)
    cp = con.sql(f"SELECT DISTINCT cell_id, WDPAID FROM read_parquet('{a.cell_pa}') JOIN k USING (WDPAID)").df()
    print(f"universe: {cp.WDPAID.nunique():,} PAs, {cp.cell_id.nunique():,} cells", flush=True)
    row, col = cp.cell_id.values // NCOLS, cp.cell_id.values % NCOLS
    with rasterio.open(a.regions) as r:
        reg = r.read(1)
    cp["region"] = reg[row, col].astype(float)
    cp.loc[cp.region == 0, "region"] = np.nan
    del reg
    x, y = X0 + (col + 0.5) * RES, Y0 - (row + 0.5) * RES
    lon, lat = Transformer.from_crs("EPSG:6933", "EPSG:4326", always_xy=True).transform(x, y)
    with rasterio.open(a.biomes) as r:
        b = r.read(1)
        rr, cc = rasterio.transform.rowcol(r.transform, lon, lat)
        rr, cc = np.clip(rr, 0, r.height - 1), np.clip(cc, 0, r.width - 1)
        v = b[rr, cc].astype(float)
        if r.nodata is not None:
            v[v == r.nodata] = np.nan
    del b
    v[(v < 1) | (v > 14)] = np.nan  # WWF biomes 1-14 (98/99 = lake, rock/ice)
    cp["biome"] = v
    out = pd.DataFrame({"n_cells": cp.groupby("WDPAID").size()})
    out["biome"] = majority(cp, "biome")
    out["region"] = majority(cp, "region")
    del cp
    first_iso = w.loc[kept, ["WDPAID", "ISO3"]].assign(iso3=lambda d: d.ISO3.str.split(";").str[0]) \
        .drop_duplicates("WDPAID").set_index("WDPAID").iso3
    out["iso3"] = first_iso.reindex(out.index)
    print(f"attributes done ({time.time() - t0:.0f} s); areas next", flush=True)

    # areas: read kept parcels of universe PAs in batches, clip to the GEDI latitude band, group by country
    ymax = Transformer.from_crs("EPSG:4326", "EPSG:6933", always_xy=True).transform(0, 51.6)[1]
    band = shapely.box(-1e8, -ymax, 1e8, ymax)
    uni = set(out.index)
    n = len(w)
    tmp = tempfile.mkdtemp(dir=a.tmp)
    for s in range(0, n, BATCH):
        g = pyogrio.read_dataframe(a.wdpa, columns=["WDPAID"], skip_features=s, max_features=BATCH)
        g = g[kept[s:s + BATCH] & g.WDPAID.isin(uni).values & ~g.geometry.isna().values]
        if g.empty:
            continue
        geom = shapely.intersection(shapely.make_valid(g.geometry.values), band)
        g = gpd.GeoDataFrame({"WDPAID": g.WDPAID.values, "iso3": out.iso3.reindex(g.WDPAID.values).values},
                             geometry=geom, crs=g.crs)
        g = g[~g.geometry.is_empty]
        g.to_parquet(os.path.join(tmp, f"b{s:07d}.parquet"))
        print(f"  read {min(s + BATCH, n):,}/{n:,} parcels", flush=True)
    parts = gpd.read_parquet(tmp) if os.listdir(tmp) else None
    area, scale = {}, {}
    for iso, g in parts.groupby("iso3"):
        per = g.groupby("WDPAID").geometry.apply(lambda s: shapely.union_all(s.values).area)
        u = shapely.union_all(g.geometry.values).area
        area.update(per.to_dict())
        scale[iso] = u / per.sum() if per.sum() > 0 else np.nan
    out["area_ha"] = pd.Series(area).reindex(out.index).fillna(0) / 1e4
    out["overlap_scale"] = out.iso3.map(scale)
    out["area_scaled_ha"] = out.area_ha * out.overlap_scale.fillna(1)
    out = out.reset_index().rename(columns={"WDPAID": "pa_id"})
    out.to_parquet(a.out, index=False)
    print(f"wrote {a.out}: {len(out):,} PAs, {out.area_ha.sum() / 1e6:.1f} Mha before / "
          f"{out.area_scaled_ha.sum() / 1e6:.1f} Mha after removing overlaps ({time.time() - t0:.0f} s)")


# ---------------------------------------------------------------- DiD, imputation, totals
def load_did(root, smd_max=SMD_MAX):
    import glob
    fr = [pd.read_csv(f).assign(did_iso3=os.path.basename(os.path.dirname(f)))
          for f in glob.glob(os.path.join(root, "*", "*_did_by_pa.csv"))]
    d = pd.concat(fr, ignore_index=True)
    if "sd_pa_t1" not in d:
        raise SystemExit("did_by_pa files have no sd_pa_t1; rerun dps/global_matching.py did")
    # a PA matched in two countries (cells on both sides): keep the matched set with more cells
    d["n_all"] = d[["n_pa_t1", "n_pa_t2", "n_ctl_t1", "n_ctl_t2"]].sum(axis=1)
    d = d.sort_values("n_all", ascending=False).drop_duplicates(["pa_id", "metric"])
    ag = d[d.metric == "agbd"].set_index("pa_id")
    smd = (ag.mean_pa_t1 - ag.mean_ctl_t1) / np.sqrt((ag.sd_pa_t1 ** 2 + ag.sd_ctl_t1 ** 2) / 2)
    d["smd_agbd_t1"] = d.pa_id.map(smd)
    d["valid"] = d.qualifies.astype(bool) & (d.smd_agbd_t1.abs() <= smd_max)
    return d


def wmean(v, w):
    m = np.isfinite(v) & (w > 0)
    return np.sum(v[m] * w[m]) / np.sum(w[m]) if m.any() else np.nan


def run(a):
    os.makedirs(a.out, exist_ok=True)
    at = pd.read_parquet(a.attrs)
    at["region_name"] = at.region.map(REGIONS)
    d = load_did(a.did, a.smd_max)
    rng = np.random.default_rng(42)
    tables, means, totals = [], {k: [] for k in ("global", "country", "biome", "region")}, \
        {k: [] for k in ("global", "country", "biome", "region")}
    for metric in METRICS:
        dm = d[d.metric == metric].set_index("pa_id")
        t = at.copy()
        t["metric"] = metric
        for c in ("did", "change_pa", "change_ctl", "mean_pa_t1", "mean_ctl_t1", "smd_agbd_t1", "qualifies", "valid"):
            t[c] = t.pa_id.map(dm[c])
        t["matched"] = t.pa_id.isin(dm.index)
        t["valid"] = t.valid.fillna(False).astype(bool)
        w = t.area_scaled_ha.values
        meas = np.where(t.valid.values)[0]
        # stratum codes per level; imputation level for each non-measured PA from the full sample
        codes, level = {}, np.full(len(t), "", object)
        level[meas] = "measured"
        for name, keys in LEVELS:
            k = t[keys].astype(str).agg("|".join, axis=1).where(t[keys].notna().all(axis=1))
            codes[name] = pd.factorize(k)[0]  # -1 = missing
            cnt = pd.Series(t.valid.values).groupby(codes[name]).sum()
            ok = codes[name] >= 0
            ok[ok] = cnt.reindex(codes[name][ok]).values >= MIN_STRATUM
            level[(level == "") & ok] = name
        level[level == ""] = "excluded"
        lv = pd.Series(level)
        t["source"] = lv.where(lv.isin(["measured", "excluded"]), "imputed:" + lv).values

        def estimate(counts):
            """DiD per PA for one bootstrap draw (counts = multiplicity of each measured PA)."""
            did_m = t.did.values
            cw = np.zeros(len(t)); cw[meas] = counts * w[meas]
            est = np.full(len(t), np.nan)
            for name, _ in LEVELS:
                c = codes[name]; ok = c >= 0
                sw = np.bincount(c[ok], weights=cw[ok], minlength=c.max() + 1)
                sd = np.bincount(c[ok], weights=(cw * np.nan_to_num(did_m))[ok], minlength=c.max() + 1)
                with np.errstate(invalid="ignore", divide="ignore"):
                    sm = sd / sw
                sel = (level == name) & ok
                est[sel] = sm[c[sel]]
            return est, cw

        point_imp, _ = estimate(np.ones(len(meas)))
        t["did_used"] = np.where(level == "measured", t.did.values, point_imp)
        tables.append(t)
        groups = []
        for by, col in (("global", None), ("country", "iso3"), ("biome", "biome"), ("region", "region_name")):
            gl = [("all", np.ones(len(t), bool))] if col is None else \
                [(gv, (t[col] == gv).values) for gv in sorted(t[col].dropna().unique())]
            for gv, g in gl:
                groups.append((by, gv, g, g & (level == "measured"), g & (level != "measured") & (level != "excluded")))

        def stat(est, cw, gm, gi):
            A_m = w[gm].sum()
            m_mean = wmean(t.did.values[gm], cw[gm]) if gm.any() else np.nan
            tot_m = m_mean * A_m if A_m > 0 else 0.0
            imp_ok = gi & np.isfinite(est)
            tot_mi = (tot_m if np.isfinite(tot_m) else 0.0) + np.sum(est[imp_ok] * w[imp_ok])
            a_mi = A_m + w[imp_ok].sum()
            return m_mean, tot_mi / a_mi if a_mi > 0 else np.nan, tot_m, tot_mi

        B = np.full((N_BOOT, len(groups), 4), np.nan)
        for k in range(N_BOOT):
            cnt = np.bincount(rng.integers(0, len(meas), len(meas)), minlength=len(meas)).astype(float)
            est, cw = estimate(cnt)
            for j, (_, _, g, gm, gi) in enumerate(groups):
                B[k, j] = stat(est, cw, gm, gi)
        lo_all, hi_all = np.nanpercentile(B, 2.5, axis=0), np.nanpercentile(B, 97.5, axis=0)
        w_meas = np.where(level == "measured", w, 0.0)
        for j, (by, gv, g, gm, gi) in enumerate(groups):
            A_m, A_i, A_all = w[gm].sum(), w[gi].sum(), w[g].sum()
            p = stat(point_imp, w_meas, gm, gi)
            lo, hi = lo_all[j], hi_all[j]
            tm = t[gm]
            row = dict(group=gv, metric=metric,
                       n_pa=int(g.sum()), n_pa_matched=int((g & t.matched.values).sum()),
                       n_pa_measured=int(gm.sum()), n_pa_imputed=int(gi.sum()),
                       n_pa_excluded=int((g & (level == "excluded")).sum()),
                       area_mha=A_all / 1e6, share_area_measured=A_m / A_all if A_all else np.nan,
                       share_area_imputed=A_i / A_all if A_all else np.nan,
                       did_measured=p[0], did_measured_lo=lo[0], did_measured_hi=hi[0],
                       did_all=p[1], did_all_lo=lo[1], did_all_hi=hi[1],
                       change_pa_measured=wmean(tm.change_pa.values, w[gm]),
                       change_ctl_measured=wmean(tm.change_ctl.values, w[gm]),
                       t1_pa_measured=wmean(tm.mean_pa_t1.values, w[gm]),
                       t1_ctl_measured=wmean(tm.mean_ctl_t1.values, w[gm]))
            means[by].append(row)
            if metric == "agbd":  # Mg/ha x ha = Mg; / 1e6 = Tg
                totals[by].append(dict(group=gv, n_pa_measured=row["n_pa_measured"],
                                       n_pa_imputed=row["n_pa_imputed"], n_pa_excluded=row["n_pa_excluded"],
                                       area_mha=row["area_mha"], share_area_imputed=row["share_area_imputed"],
                                       total_tg_measured=p[2] / 1e6, total_tg_measured_lo=lo[2] / 1e6,
                                       total_tg_measured_hi=hi[2] / 1e6,
                                       total_tg_all=p[3] / 1e6, total_tg_all_lo=lo[3] / 1e6,
                                       total_tg_all_hi=hi[3] / 1e6))
        print(f"{metric}: {len(meas):,} measured, {(level != 'measured').sum() - (level == 'excluded').sum():,} imputed, "
              f"{(level == 'excluded').sum():,} excluded of {len(t):,} PAs", flush=True)
    pd.concat(tables, ignore_index=True).to_parquet(os.path.join(a.out, "pa_table.parquet"), index=False)
    for by in means:
        pd.DataFrame(means[by]).to_csv(os.path.join(a.out, f"means_{by}.csv"), index=False)
        pd.DataFrame(totals[by]).to_csv(os.path.join(a.out, f"totals_agbd_{by}.csv"), index=False)
    print(f"wrote {a.out}/pa_table.parquet, means_*.csv, totals_agbd_*.csv")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("attrs")
    s.add_argument("--wdpa", required=True); s.add_argument("--biomes", required=True)
    s.add_argument("--regions", required=True); s.add_argument("--cell-pa", default="output/cell_pa.parquet")
    s.add_argument("--out", default="output/totals_pa_attrs.parquet"); s.add_argument("--tmp", default=None)
    r = sub.add_parser("run")
    r.add_argument("--attrs", default="output/totals_pa_attrs.parquet"); r.add_argument("--did", default="results/global_did")
    r.add_argument("--out", default="results/totals")
    r.add_argument("--smd-max", type=float, default=SMD_MAX,
                   help="T1 AGBD balance threshold (headline 0.25; sensitivity 0.1 and 0.5)")
    a = ap.parse_args()
    {"attrs": attrs, "run": run}[a.cmd](a)


if __name__ == "__main__":
    main()
