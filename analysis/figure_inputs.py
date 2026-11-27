#!/usr/bin/env python
"""
Figure inputs (figures/fig_setup.R) from the totals.py outputs, plus PA biomass stock and stock change.

  python analysis/figure_inputs.py --totals results/totals_drift_v1 \
      --cells output/gedi_pa_change_cells_drift_v1_combined.parquet --out /home/jovyan/GEDI/results/figure_inputs

DiD and avoided loss come unchanged from totals.py (headline rules, CLAUDE.md "Totals"):
  qualifies          = measured (10-cell rule and |T1 AGBD SMD| <= threshold)
  did_mean           = area-weighted DiD, measured + imputed (= did_all in means_*.csv)
  avoided_loss_Tg    = sum(area x DiD), measured + imputed (= total_tg_all in totals_agbd_*.csv)
  change_pa_mean / change_ctl_mean = area-weighted PA and control change of measured PAs
                       (their difference is the measured-only DiD)

Stock and stock change (absolute; exposed to any drift the correction misses), decided 2026-11-27:
  Area = total PA area (dissolved, overlaps removed, area_scaled_ha), because the AGBD means use all
  GEDI shots, not only forest. Per PA, T1 and T2 means of cell means over all cells with their centre
  in the PA. A PA has its own values when it has >= 10 cells in each epoch; otherwise T1 mean and
  change are imputed as the area-weighted mean of >= 5 such PAs in biome x country, else GEDI
  region x biome, else biome (same hierarchy as the DiD). stock_t1 = area x T1 mean;
  stock_change = area x (T2 - T1). 95% CI: bootstrap over PAs with own values, imputation redone.

Outputs: summary_{global,country,continent,biome,continent_biome}.csv, boot_{country,global}.csv,
pa_table.csv.
"""
import argparse
import os

import duckdb
import numpy as np
import pandas as pd

AGBD = ["agbd", "agbd_c", "agbd_c_lo", "agbd_c_hi"]
LEVELS = [("biome", "iso3"), ("region", "biome"), ("biome",)]
MIN_CELLS, MIN_STRATUM, N_BOOT = 10, 5, 1000


def own_means(cells, cell_pa, ids, metrics):
    con = duckdb.connect()
    con.sql("SET memory_limit = '8GB'; SET preserve_insertion_order = false;")
    con.register("u", pd.DataFrame({"WDPAID": ids}))
    aggs = ", ".join(f"count({m}_mean) AS n_{m}, avg({m}_mean) AS m_{m}" for m in metrics)
    pe = con.sql(f"""SELECT cp.WDPAID AS pa_id, c.epoch, {aggs}
                     FROM read_parquet('{cells}') c
                     JOIN (SELECT DISTINCT cell_id, WDPAID FROM read_parquet('{cell_pa}') JOIN u USING (WDPAID)) cp
                     USING (cell_id) GROUP BY 1, 2""").df()
    out = {}
    for m in metrics:
        w = pe.pivot(index="pa_id", columns="epoch", values=[f"n_{m}", f"m_{m}"])
        w.columns = [f"{s[0]}_{e}" for s, e in w.columns]
        w = w.reindex(columns=["n_T1", "n_T2", "m_T1", "m_T2"]).fillna({"n_T1": 0, "n_T2": 0})
        w["own"] = (w.n_T1 >= MIN_CELLS) & (w.n_T2 >= MIN_CELLS)
        out[m] = w
    return out


def stock(t, own, rng, n_boot):
    """Per-PA T1 mean and change (own or imputed) and bootstrap country/global stock-change totals."""
    w = t.area_scaled_ha.to_numpy()
    ok_own = own.own.reindex(t.pa_id).fillna(False).to_numpy()
    t1 = own.m_T1.reindex(t.pa_id).to_numpy()
    ch = (own.m_T2 - own.m_T1).reindex(t.pa_id).to_numpy()
    idx = np.where(ok_own)[0]
    codes, level = [], np.where(ok_own, "own", "")
    for keys in LEVELS:
        k = t[list(keys)].astype(str).agg("|".join, axis=1).where(t[list(keys)].notna().all(axis=1))
        c = pd.factorize(k)[0]
        cnt = np.bincount(c[ok_own & (c >= 0)], minlength=c.max() + 1)
        sel = (level == "") & (c >= 0)
        sel[sel] = cnt[c[sel]] >= MIN_STRATUM
        level[sel] = "+".join(keys)
        codes.append(c)
    level[level == ""] = "none"

    def est(counts):
        cw = np.zeros(len(t)); cw[idx] = counts * w[idx]
        v1, vc = np.where(ok_own, t1, np.nan), np.where(ok_own, ch, np.nan)
        for keys, c in zip(LEVELS, codes):
            ok = c >= 0
            sw = np.bincount(c[ok], weights=cw[ok], minlength=c.max() + 1)
            s1 = np.bincount(c[ok], weights=(cw * np.nan_to_num(t1))[ok], minlength=c.max() + 1)
            sc = np.bincount(c[ok], weights=(cw * np.nan_to_num(ch))[ok], minlength=c.max() + 1)
            sel = (level == "+".join(keys)) & ok
            with np.errstate(invalid="ignore", divide="ignore"):
                v1[sel], vc[sel] = (s1 / sw)[c[sel]], (sc / sw)[c[sel]]
        return v1, vc, cw

    v1, vc, _ = est(np.ones(len(idx)))
    iso = pd.factorize(t.iso3)
    ci, nc = iso[0], len(iso[1])
    has = ci >= 0
    a_own = np.bincount(ci[has & ok_own], weights=w[has & ok_own], minlength=nc)
    chz = np.nan_to_num(ch)

    def totals(vc, cw):
        """Per country: own area x resampled own mean change + imputed PAs; last element = global."""
        sw = np.bincount(ci[has], weights=cw[has], minlength=nc)
        sx = np.bincount(ci[has], weights=(cw * chz)[has], minlength=nc)
        with np.errstate(invalid="ignore", divide="ignore"):
            own_tot = np.where(sw > 0, sx / sw, 0.0) * a_own
        imp = has & ~ok_own & np.isfinite(vc)
        tot = own_tot + np.bincount(ci[imp], weights=(vc * w)[imp], minlength=nc)
        return np.append(tot, tot.sum())

    B = np.zeros((n_boot, nc + 1))
    for b in range(n_boot):
        _, bc, cw = est(np.bincount(rng.integers(0, len(idx), len(idx)), minlength=len(idx)).astype(float))
        B[b] = totals(bc, cw)
    lo, hi = np.percentile(B, 2.5, axis=0) / 1e6, np.percentile(B, 97.5, axis=0) / 1e6
    bc = pd.DataFrame({"iso3": iso[1], "stock_change_Tg_lo": lo[:-1], "stock_change_Tg_hi": hi[:-1]})
    return v1, vc, level, bc, (lo[-1] / 1e3, hi[-1] / 1e3)


def wmean(v, w):
    m = np.isfinite(v) & (w > 0)
    return np.sum(v[m] * w[m]) / np.sum(w[m]) if m.any() else np.nan


def summarise(t, by):
    rows = []
    groups = t.groupby(by, dropna=True) if by else [((), t)]
    for key, g in groups:
        key = key if isinstance(key, tuple) else (key,)
        a, meas = g.area_scaled_ha.to_numpy(), g.valid.to_numpy()
        r = dict(zip(by, key))
        r.update(n_pas=len(g), n_matched=int(g.matched.sum()), n_qualifying=int(meas.sum()),
                 n_imputed=int(g.source.str.startswith("imputed").sum()), n_excluded=int((g.source == "excluded").sum()),
                 area_ha=a.sum(), area_qualifying_ha=a[meas].sum(),
                 area_with_counterfactual_ha=a[g.did_used.notna().to_numpy()].sum())
        r["pct_pas_qualifying"] = 100 * r["n_qualifying"] / r["n_pas"]
        r["pct_area_qualifying"] = 100 * r["area_qualifying_ha"] / r["area_ha"] if r["area_ha"] else np.nan
        r["pct_area_with_counterfactual"] = (100 * r["area_with_counterfactual_ha"] / r["area_ha"]
                                             if r["area_ha"] else np.nan)
        r["did_mean"] = wmean(g.did_used.to_numpy(), a)
        r["did_measured_mean"] = wmean(g.did.to_numpy(), np.where(meas, a, 0))
        r["change_pa_mean"] = wmean(g.change_pa.to_numpy(), np.where(meas, a, 0))
        r["change_ctl_mean"] = wmean(g.change_ctl.to_numpy(), np.where(meas, a, 0))
        if g.metric.iloc[0] in AGBD:
            r["stock_t1_Tg"] = np.nansum(a * g.stock_t1_mean) / 1e6
            r["stock_change_Tg"] = np.nansum(a * g.stock_change_mean) / 1e6
            r["avoided_loss_Tg"] = np.nansum(a * g.did_used) / 1e6
            r["pct_area_stock_own"] = 100 * a[(g.stock_source == "own").to_numpy()].sum() / r["area_ha"] \
                if r["area_ha"] else np.nan
        rows.append(r)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--totals", default="results/totals_drift_v1")
    ap.add_argument("--cells", default="output/gedi_pa_change_cells_drift_v1_combined.parquet")
    ap.add_argument("--cell-pa", default="output/cell_pa.parquet")
    ap.add_argument("--out", default="/home/jovyan/GEDI/results/figure_inputs")
    ap.add_argument("--n-boot", type=int, default=N_BOOT)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    pt = pd.read_parquet(os.path.join(a.totals, "pa_table.parquet"))
    pt["continent"] = pt.region_name
    metrics = [m for m in AGBD if m in set(pt.metric)]
    own = own_means(a.cells, a.cell_pa, pt.pa_id.unique(), metrics)
    rng = np.random.default_rng(7)

    tabs, boots_c, boots_g = [], [], []
    tot_c = pd.read_csv(os.path.join(a.totals, "totals_agbd_country.csv"))
    tot_g = pd.read_csv(os.path.join(a.totals, "totals_agbd_global.csv"))
    for m, t in pt.groupby("metric", sort=False):
        t = t.reset_index(drop=True).copy()
        if m in metrics:
            v1, vc, lev, bc, gci = stock(t, own[m], rng, a.n_boot)
            t["stock_t1_mean"], t["stock_change_mean"], t["stock_source"] = v1, vc, lev
            o = own[m].reindex(t.pa_id)
            t["n_cells_t1"], t["n_cells_t2"] = o.n_T1.fillna(0).to_numpy(), o.n_T2.fillna(0).to_numpy()
            tc = tot_c[tot_c.metric == m].rename(columns={"group": "iso3", "total_tg_all_lo": "avoided_loss_Tg_lo",
                                                          "total_tg_all_hi": "avoided_loss_Tg_hi"})
            boots_c.append(bc.merge(tc[["iso3", "avoided_loss_Tg_lo", "avoided_loss_Tg_hi"]], on="iso3", how="left")
                           .assign(metric=m))
            tg = tot_g[tot_g.metric == m].iloc[0]
            boots_g.append(dict(metric=m, stock_change_Pg_lo=gci[0], stock_change_Pg_hi=gci[1],
                                avoided_loss_Pg_lo=tg.total_tg_all_lo / 1e3, avoided_loss_Pg_hi=tg.total_tg_all_hi / 1e3))
            print(f"{m}: stock values own for {(lev == 'own').sum():,} PAs, imputed {((lev != 'own') & (lev != 'none')).sum():,}, "
                  f"none {(lev == 'none').sum():,}", flush=True)
        tabs.append(t)
    t = pd.concat(tabs, ignore_index=True)
    t["qualifies"] = t.valid
    t["pa_t1"] = t.mean_pa_t1  # matched PA cells (DiD sample)

    summ = {}
    for name, by in (("global", []), ("country", ["iso3"]), ("continent", ["continent"]),
                     ("biome", ["biome"]), ("continent_biome", ["continent", "biome"])):
        s = pd.concat([summarise(g, by).assign(metric=m) for m, g in t.groupby("metric", sort=False)], ignore_index=True)
        summ[name] = s[["metric"] + [c for c in s.columns if c != "metric"]]
        summ[name].to_csv(os.path.join(a.out, f"summary_{name}.csv"), index=False)

    # checks against totals.py
    g = summ["global"].set_index("metric")
    mg = pd.read_csv(os.path.join(a.totals, "means_global.csv")).set_index("metric")
    assert np.allclose(g.did_mean, mg.did_all.reindex(g.index), equal_nan=True), "did_mean != did_all"
    ga = g.loc[metrics, "avoided_loss_Tg"].to_numpy()
    assert np.allclose(ga, tot_g.set_index("metric").total_tg_all.reindex(metrics).to_numpy(), rtol=1e-6)
    c = summ["country"][summ["country"].metric.isin(metrics)].merge(
        tot_c.rename(columns={"group": "iso3"}), on=["metric", "iso3"])
    assert np.allclose(c.avoided_loss_Tg, c.total_tg_all, rtol=1e-6, atol=1e-6), "country avoided loss mismatch"

    pd.concat(boots_c, ignore_index=True).to_csv(os.path.join(a.out, "boot_country.csv"), index=False)
    pd.DataFrame(boots_g).to_csv(os.path.join(a.out, "boot_global.csv"), index=False)
    keep = ["pa_id", "metric", "iso3", "continent", "biome", "area_scaled_ha", "matched", "qualifies", "source",
            "smd_agbd_t1", "pa_t1", "change_pa", "change_ctl", "did", "did_used",
            "stock_t1_mean", "stock_change_mean", "stock_source", "n_cells_t1", "n_cells_t2"]
    t[keep].to_csv(os.path.join(a.out, "pa_table.csv"), index=False)
    print(g.loc[metrics, ["did_mean", "stock_t1_Tg", "stock_change_Tg", "avoided_loss_Tg", "pct_area_stock_own"]]
          .round(2).to_string())
    print(pd.DataFrame(boots_g).round(3).to_string(index=False))
    print(f"checks passed; wrote {a.out}/")


if __name__ == "__main__":
    main()
