#!/usr/bin/env python
"""
Ground-finding drift correction for GEDI RH metrics.

FIT (from ground_drift_stats.parquet written by ground_drift_quicklook.py; no GEDI re-read):
  For each Landsat tree-cover class, the within-tile year effect (ground residual vs TanDEM-X in
  year Y minus the same tile's T1 mean, gentle slopes) is regressed on year, pooled over tiles
  with weights n_Y*n_T1/(n_Y+n_T1). The slope is the drift rate (m/yr).
  The open-ground rate (0-10% cover) is subtracted: it is a whole-waveform vertical shift that
  moves ground and canopy top together and so does not bias RH. The remainder,
      excess_rate(cover) = rate(cover) - rate(0-10%),
  is the rate at which missed ground under canopy lowers RH. 95% CIs come from resampling tiles.

  Years listed in --offset-years (default 2019) are left out of the rate fit and get their own
  offset instead: their excess year effect minus the fitted line. The early-mission year sat
  below the 2020+ trend under dense canopy, and fitting through it overstated the rate.

APPLY (per shot, both epochs):
      drift(cover, y)  = excess_rate(cover) * (y - ref_fit) + offset_y(cover)   [offset only for offset years]
      correction_m     = drift(cover, y) - mean of drift over T1 shots
      RH_k_corrected   = RH_k + correction_m      (every percentile k)
  The correction averages to zero over T1, so T1 means are unchanged and every year is brought
  onto T1's ground finding. Parameters are interpolated linearly in cover between class midpoints.
  `which` = "central" | "low" | "high" uses the rate estimate or a CI bound (offsets stay central).

  python analysis/drift_correction.py fit --stats .../ground_drift_stats.parquet --out .../drift_correction.csv

  from drift_correction import load, correction
  table = load(".../drift_correction.csv")
  df["corr_m"] = correction(df.cover, df.yr, table)          # numpy/pandas vectorised
"""
import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ORDER = ["0-10%", "10-30%", "30-50%", "50-70%", "70-80%", "80-90%", "90-100%"]
MIDPOINT = {"0-10%": 5, "10-30%": 20, "30-50%": 40, "50-70%": 60, "70-80%": 75,
            "80-90%": 85, "90-100%": 95}
OPEN = "0-10%"


def year_effects(st, min_n):
    """Within-tile year effect per tile x cover x year, relative to the tile's T1 mean."""
    g = st.groupby(["tile_id", "cover_bin", "yr", "epoch"])[["n", "s"]].sum().reset_index()
    t1 = g[g.epoch == "T1"].groupby(["tile_id", "cover_bin"])[["n", "s"]].sum().reset_index()
    t1["m_t1"] = t1.s / t1.n
    yr = g.groupby(["tile_id", "cover_bin", "yr"])[["n", "s"]].sum().reset_index()
    yr["m"] = yr.s / yr.n
    j = yr.merge(t1[["tile_id", "cover_bin", "m_t1", "n"]], on=["tile_id", "cover_bin"],
                 suffixes=("", "_t1"))
    j = j[(j.n >= min_n) & (j.n_t1 >= min_n)].copy()
    j["eff"] = j.m - j.m_t1
    j["w"] = j.n * j.n_t1 / (j.n + j.n_t1)
    return j


def wls_slope(d):
    x, y, w = d.yr.to_numpy(float), d.eff.to_numpy(float), d.w.to_numpy(float)
    if len(np.unique(x)) < 2:
        return np.nan
    xm = np.average(x, weights=w)
    return np.sum(w * (x - xm) * (y - np.average(y, weights=w))) / np.sum(w * (x - xm) ** 2)


def fit(st, min_n=30, n_boot=500, seed=0, offset_years=(2019,)):
    st = st[st.gentle & st.cover_bin.isin(ORDER)]
    t1_all = st[st.epoch == "T1"].groupby("yr").n.sum()
    off_years = [y for y in offset_years if y in set(st.yr)]
    fit_st = st[~st.yr.isin(off_years)]
    ye = year_effects(fit_st, min_n)
    t1 = fit_st[fit_st.epoch == "T1"].groupby("yr").n.sum()
    ref_year = float(np.average(t1.index, weights=t1.values))  # baseline of the fitted year effects
    tiles = ye.tile_id.unique()
    by = {t: x for t, x in ye.groupby("tile_id")}
    rng = np.random.default_rng(seed)

    def rates(frame):
        return {cb: wls_slope(frame[frame.cover_bin == cb]) for cb in ORDER}

    est = rates(ye)
    boots = []
    for _ in range(n_boot):
        pick = rng.choice(tiles, size=len(tiles), replace=True)
        boots.append(rates(pd.concat([by[t] for t in pick])))
    boots = pd.DataFrame(boots)
    excess_boot = boots.sub(boots[OPEN], axis=0)
    rows = []
    for cb in ORDER:
        rows.append(dict(
            cover_bin=cb, cover_mid=MIDPOINT[cb],
            n_tiles=int(ye[ye.cover_bin == cb].tile_id.nunique()),
            rate=est[cb], rate_lo=np.nanpercentile(boots[cb], 2.5),
            rate_hi=np.nanpercentile(boots[cb], 97.5),
            excess_rate=est[cb] - est[OPEN],
            excess_lo=np.nanpercentile(excess_boot[cb], 2.5),
            excess_hi=np.nanpercentile(excess_boot[cb], 97.5)))
    tab = pd.DataFrame(rows)
    tab.loc[tab.cover_bin == OPEN, ["excess_rate", "excess_lo", "excess_hi"]] = 0.0
    tab["ref_fit"] = ref_year

    # Offsets for left-out years, relative to the same baseline (tile mean over fitted T1 years).
    if off_years:
        base = fit_st[fit_st.epoch == "T1"].groupby(["tile_id", "cover_bin"])[["n", "s"]].sum()
        base["m_t1"] = base.s / base.n
        for y in off_years:
            g = st[st.yr == y].groupby(["tile_id", "cover_bin"])[["n", "s"]].sum()
            g["m"] = g.s / g.n
            j = g.join(base[["m_t1", "n"]], rsuffix="_t1", how="inner").reset_index()
            j = j[(j.n >= min_n) & (j.n_t1 >= min_n)]
            j["eff"] = j.m - j.m_t1
            j["w"] = j.n * j.n_t1 / (j.n + j.n_t1)
            eff = {cb: np.average(x.eff, weights=x.w) for cb, x in j.groupby("cover_bin")}
            open_line = est[OPEN] * (y - ref_year)
            col = []
            for _, r in tab.iterrows():
                e = eff.get(r.cover_bin, np.nan)
                excess_eff = e - eff.get(OPEN, np.nan)
                col.append(0.0 if r.cover_bin == OPEN else excess_eff - r.excess_rate * (y - ref_year))
            tab[f"offset_{y}"] = col
    # T1 shot weights by year (all T1 years), to make the correction average zero over T1
    tab["t1_year_weights"] = json.dumps({int(k): int(v) for k, v in t1_all.items()})
    return tab, ye


def load(path):
    return pd.read_csv(path)


def _drift(rate, offsets, years, ref_fit):
    d = rate * (years - ref_fit)
    for y, off in offsets.items():
        d = d + np.where(years == y, off, 0.0)
    return d


def correction(cover, year, table, which="central"):
    """Metres to ADD to every RH percentile of a shot (vectorised over cover and year)."""
    col = {"central": "excess_rate", "low": "excess_lo", "high": "excess_hi"}[which]
    t = table.sort_values("cover_mid")
    c = np.asarray(cover, float)
    yrs = np.asarray(year, float)
    rate = np.interp(c, t.cover_mid, t[col].clip(lower=0))
    off_cols = [k for k in t.columns if k.startswith("offset_")]
    offsets = {int(k.split("_")[1]): np.interp(c, t.cover_mid, t[k].fillna(0)) for k in off_cols}
    ref_fit = float(t.ref_fit.iloc[0])
    w = json.loads(t.t1_year_weights.iloc[0])
    wy = np.array([float(k) for k in w]), np.array([float(v) for v in w.values()])
    # mean drift over T1 shots (per shot cover), so T1 means are unchanged
    mean_t1 = sum(_drift(rate, offsets, np.full_like(c, y), ref_fit) * n for y, n in zip(*wy)) / wy[1].sum()
    return _drift(rate, offsets, yrs, ref_fit) - mean_t1


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fit")
    f.add_argument("--stats", required=True)
    f.add_argument("--out", required=True, help="correction table CSV")
    f.add_argument("--min-n", type=int, default=30)
    f.add_argument("--n-boot", type=int, default=500)
    f.add_argument("--offset-years", type=int, nargs="*", default=[2019],
                   help="years fitted as separate offsets, not part of the rate (empty list: none)")
    a = ap.parse_args()

    st = pd.read_parquet(a.stats)
    tab, ye = fit(st, a.min_n, a.n_boot, offset_years=a.offset_years)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    tab.to_csv(a.out, index=False)
    with pd.option_context("display.width", 160, "display.float_format", "{:.3f}".format):
        print(f"rate fitted on years excluding {a.offset_years}; baseline year {tab.ref_fit.iloc[0]:.2f}")
        print(tab.drop(columns=["ref_fit", "t1_year_weights"]).to_string(index=False))
    yrs = np.arange(2019, 2027)
    print("\nRH correction (m) at the central estimate, by cover and year:")
    print(pd.DataFrame({y: correction([5, 40, 60, 75, 85, 95], [y] * 6, tab) for y in yrs},
                       index=["5%", "40%", "60%", "75%", "85%", "95%"]).round(2).to_string())

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].errorbar(tab.cover_mid, tab.excess_rate,
                     yerr=[tab.excess_rate - tab.excess_lo, tab.excess_hi - tab.excess_rate],
                     fmt="o-", capsize=3)
    axes[0].axhline(0, color="grey", lw=0.8)
    axes[0].set_xlabel("Landsat tree cover (%)")
    axes[0].set_ylabel("Excess ground drift (m / year)")
    axes[0].set_title("Drift rate beyond open ground (lowers RH)", fontsize=10)
    for cb, c in (("0-10%", "C0"), ("50-70%", "C1"), ("90-100%", "C2")):
        d = ye[ye.cover_bin == cb].groupby("yr").apply(
            lambda x: np.average(x.eff, weights=x.w), include_groups=False)
        axes[1].plot(d.index, d.values, "o", color=c, label=f"cover {cb}")
        r = tab.loc[tab.cover_bin == cb, "rate"].iloc[0]
        ref = tab.ref_fit.iloc[0]
        axes[1].plot(yrs, r * (yrs - ref), "-", color=c, alpha=0.6)
    axes[1].axhline(0, color="grey", lw=0.8)
    axes[1].set_xlabel("Year")
    axes[1].set_ylabel("Year effect within tiles (m)")
    axes[1].set_title("Fitted linear drift vs yearly estimates", fontsize=10)
    axes[1].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    png = os.path.splitext(a.out)[0] + ".png"
    fig.savefig(png, dpi=150)
    print(f"\nwrote {a.out} and {png}")


if __name__ == "__main__":
    main()
