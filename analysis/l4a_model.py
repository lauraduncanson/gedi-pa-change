#!/usr/bin/env python
"""
GEDI L4A (version 3) biomass models, evaluated from the model table stored in the granules
(ANCILLARY/model_data, exported by find_l4a_models.py), plus a validation against the tiled data.

Model form, as stored (all 35 strata: x_transform = y_transform = sqrt, Snowdon bias correction):
    x_i  = sqrt(RH_{rh_index[i]} + 100)              RH in metres; +100 keeps negative RH valid
    (predictor_max_value / response_max_value only drive the limit flags; values are not capped)
    y    = par[0] + sum_i par[i+1] * x_i              y = sqrt(AGBD)
    y    = max(y, 0)                                  (negatives in sqrt space set to zero)
    AGBD = bias_correction_value * y**2               (Snowdon)
Confirmed against the tiled data (validate): +100 offset, rh_l2a as the RH source, no capping,
negatives in sqrt space set to zero. Strata whose applied coefficients differ from the stored
table can be overridden with an --overrides JSON {stratum: {"par": [...], "bc": ...}}.

VALIDATE (run on the Hub):
  python analysis/l4a_model.py validate --models /home/jovyan/GEDI/results/l4a_models/ANCILLARY__model_data.json
It reads a sample of shots from the tiled data and reports:
  A. AGBD recomputed from the stored predictors (xvar_l4a) vs agbd_l4a   -> confirms the model form
  B. which RH source reproduces xvar_l4a: rh_l2a (L2A selected algorithm) or the per-algorithm
     geolocation_rh_aN arrays (cm)                                       -> where to apply the correction
  C. AGBD recomputed from that RH source vs agbd_l4a                    -> end-to-end check
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "extract_tile"))

OFFSET = 100.0


class L4AModels:
    def __init__(self, path, overrides=None):
        rows = json.load(open(path))
        self.m = {}
        for r in rows:
            k = int(r["npar"]) - 1
            self.m[r["predict_stratum"]] = dict(
                rh=np.array(r["rh_index"][:k], int), b=np.array(r["par"][:k + 1], float),
                xmax=np.array(r["predictor_max_value"][:k], float),
                bc=float(r["bias_correction_value"]), ymax=float(r["response_max_value"]))
        for k, v in (json.load(open(overrides)) if overrides else {}).items():
            self.m[k]["b"] = np.array(v["par"], float)
            self.m[k]["bc"] = float(v.get("bc", self.m[k]["bc"]))
            print(f"override applied for {k}: par={v['par']}, bc={self.m[k]['bc']}")

    def rh_needed(self):
        return sorted({int(v) for m in self.m.values() for v in m["rh"]})

    def predict_from_x(self, stratum, x, cap=False):
        """x: (n, k) transformed predictors for shots of one stratum."""
        m = self.m[stratum]
        x = np.asarray(x, float)[:, :len(m["rh"])]
        if cap:
            x = np.minimum(x, m["xmax"])
        y = m["b"][0] + x @ m["b"][1:]
        y = np.maximum(y, 0.0)  # negative predictions in sqrt space are set to zero (AGBD = 0)
        agbd = m["bc"] * y ** 2
        return np.minimum(agbd, m["ymax"]) if cap else agbd

    def predict_from_rh(self, stratum, rh_lookup, rh_shift=0.0, cap=False):
        """rh_lookup(percentile) -> array of RH (m) for the shots; rh_shift (m) is added to every RH."""
        m = self.m[stratum]
        x = np.column_stack([np.sqrt(np.maximum(rh_lookup(int(p)) + rh_shift + OFFSET, 0))
                             for p in m["rh"]])
        return self.predict_from_x(stratum, x, cap)


def read_sample(tile, n_max):
    import extract_tile as E
    con = E.connect()
    files = E.tile_files(con, f"{E.ROOT}/manifest.txt", tile, [2020, 2024])
    if not files:
        return None
    cols = ", ".join([f"geolocation_rh_a{i}_l2a" for i in (1, 2, 5, 10)])
    return con.sql(f"""
        SELECT predict_stratum_l4a AS stratum, xvar_l4a AS xvar, agbd_l4a AS agbd,
               rh_l2a, {cols}, selected_algorithm_l2a AS alg_l2a,
               selected_mode_flag_l4a, predictor_limit_flag_l4a AS plim,
               response_limit_flag_l4a AS rlim, l4a_quality_flag_rel3_l4a AS q
        FROM read_parquet([{", ".join(f"'{f}'" for f in files)}], union_by_name=true)
        WHERE agbd_l4a >= 0 AND predict_stratum_l4a IS NOT NULL AND xvar_l4a IS NOT NULL
        USING SAMPLE {n_max} ROWS
    """).df()


def rel_err(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return np.abs(a - b) / np.maximum(np.abs(b), 1.0)


def validate(models_path, tiles, n_max, overrides=None):
    M = L4AModels(models_path, overrides)
    frames = [read_sample(t, n_max) for t in tiles]
    d = pd.concat([f for f in frames if f is not None and len(f)], ignore_index=True)
    d = d[d.stratum.isin(M.m.keys())]
    print(f"{len(d):,} shots, strata: {d.stratum.value_counts().to_dict()}")

    # A. from stored predictors
    xv = np.array([np.asarray(v, float)[:4] if v is not None else np.full(4, np.nan) for v in d.xvar])
    pa = np.full(len(d), np.nan)
    for s, idx in d.groupby("stratum").indices.items():
        pa[idx] = M.predict_from_x(s, xv[idx])
    ea = rel_err(pa, d.agbd)
    print(f"\nA. AGBD from xvar_l4a: median rel. error {np.nanmedian(ea):.2e}, "
          f"share within 0.1%: {(ea < 1e-3).mean():.3f}")
    print(f"   first xvar values: {xv[:3, :2].round(4).tolist()}  (sqrt(RH+100) would be ~10-13)")

    # B. which RH source reproduces xvar
    sources = {"rh_l2a": lambda r, p: np.asarray(r.rh_l2a, float)[p]}
    for i in (1, 2, 5, 10):
        sources[f"a{i}"] = (lambda col: lambda r, p: np.asarray(r[col], float)[p] / 100.0)(
            f"geolocation_rh_a{i}_l2a")
    print("\nB. share of shots whose xvar_l4a matches sqrt(RH+100) from each source (|dx| < 1e-3):")
    best = {}
    for name, get in sources.items():
        ok = np.zeros(len(d), bool)
        for j, (_, r) in enumerate(d.iterrows()):
            m = M.m[r.stratum]
            try:
                x = np.sqrt(np.array([get(r, int(p)) for p in m["rh"]]) + OFFSET)
                ok[j] = np.all(np.abs(x - xv[j, :len(x)]) < 1e-3)
            except Exception:
                pass
        best[name] = ok
        print(f"   {name:7s}: {ok.mean():.3f}")
    match_any = np.column_stack(list(best.values())).any(axis=1)
    print(f"   any source: {match_any.mean():.3f}")
    by_alg = pd.crosstab(d.alg_l2a, best["rh_l2a"], normalize="index").round(3)
    print("   rh_l2a match rate by L2A selected algorithm:\n", by_alg.to_string())

    # C. end to end from the best single source
    src = max(best, key=lambda k: best[k].mean())
    pc = np.full(len(d), np.nan)
    for s, idx in d.groupby("stratum").indices.items():
        sub = d.iloc[idx]
        getter = sources[src]
        pc[idx] = M.predict_from_rh(
            s, lambda p: np.array([getter(r, p) for _, r in sub.iterrows()]))
    ec = rel_err(pc, d.agbd)
    print(f"\nC. AGBD from {src}: median rel. error {np.nanmedian(ec):.2e}, "
          f"share within 0.1%: {(ec < 1e-3).mean():.3f}")
    # D. where do mismatches come from?
    yv = np.full(len(d), np.nan)
    for s_, idx in d.groupby("stratum").indices.items():
        m = M.m[s_]
        x = np.minimum(xv[idx, :len(m["rh"])], m["xmax"])
        yv[idx] = m["b"][0] + x @ m["b"][1:]
    d = d.assign(pred=pa, err=ea, y=yv, bad=ea >= 1e-3,
                 y_neg=yv < 0, agbd_zero=d.agbd == 0,
                 ratio=np.where(pa > 0, d.agbd / np.maximum(pa, 1e-9), np.nan))
    print("\nD. mismatch breakdown (share of shots with > 0.1% error):")
    for col in ("stratum", "y_neg", "agbd_zero", "plim", "rlim", "q", "selected_mode_flag_l4a"):
        t = d.groupby(col).agg(n=("bad", "size"), share_bad=("bad", "mean"))
        print(f"  by {col}:\n" + t.round(3).to_string().replace("\n", "\n    "))
    bad = d[d.bad]
    if len(bad):
        print("  ratio agbd_l4a / recomputed among mismatches: "
              f"median {np.nanmedian(bad.ratio):.4f}, quantiles 5/95% "
              f"{np.nanpercentile(bad.ratio, 5):.4f} / {np.nanpercentile(bad.ratio, 95):.4f}")
        with pd.option_context("display.width", 200):
            print(bad[["stratum", "agbd", "pred", "y", "plim", "rlim", "q"]].head(12)
                  .assign(x0=[np.round(v[:2], 3).tolist() for v in xv[bad.index[:12]]]).to_string())

    # E. recover the coefficients actually applied, for strata that mismatch
    print("\nE. applied vs stored coefficients (fit of sqrt(agbd_l4a) on xvar, nonzero shots):")
    print("   fitted values equal sqrt(bias_correction) * par; compare with the stored row")
    for s_, g in d.groupby("stratum"):
        if g.bad.mean() < 0.01 or len(g) < 20:
            continue
        m = M.m[s_]
        k = len(m["rh"])
        gg = g[g.agbd > 0]
        X = np.column_stack([np.ones(len(gg))] + [xv[gg.index, i] for i in range(k)])
        coef, *_ = np.linalg.lstsq(X, np.sqrt(gg.agbd.to_numpy(float)), rcond=None)
        resid = np.sqrt(gg.agbd.to_numpy(float)) - X @ coef
        stored = np.sqrt(m["bc"]) * m["b"]
        print(f"   {s_}: n={len(gg)}, fitted {np.round(coef, 4).tolist()}, max |resid| {np.abs(resid).max():.2e}")
        print(f"   {' ' * len(s_)}  stored  {np.round(stored, 4).tolist()} (par {np.round(m['b'], 4).tolist()}, bc {m['bc']:.4f})")
        for name, mm in M.m.items():
            if len(mm["rh"]) == k and np.allclose(np.sqrt(mm["bc"]) * mm["b"], coef, rtol=1e-3):
                print(f"   {' ' * len(s_)}  -> matches stored row {name}")
        print(f"   {' ' * len(s_)}  as par with this bc: {np.round(coef / np.sqrt(m['bc']), 4).tolist()}")

    print("\nIf A is ~0 the model form is right. B shows which RH array L4A used per shot; "
          "the drift correction is added to that RH. If C ~ A, recomputation is exact.")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate")
    v.add_argument("--models", required=True)
    v.add_argument("--tiles", nargs="+", default=["S13_E028", "S03_W060", "N42_W073"])
    v.add_argument("--n", type=int, default=3000, help="shots sampled per tile")
    v.add_argument("--overrides", default=None, help="JSON of per-stratum par/bc overrides")
    a = ap.parse_args()
    validate(a.models, a.tiles, a.n, a.overrides)


if __name__ == "__main__":
    main()
