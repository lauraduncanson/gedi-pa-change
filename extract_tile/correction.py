"""
Ground-drift correction of GEDI RH metrics and L4A biomass recomputation, for extract_tile.py.

Self-contained (numpy/pandas only) so it runs in the DPS container. Same logic as
analysis/drift_correction.py (correction) and analysis/l4a_model.py (L4AModels), which were
validated on the Hub: recomputed AGBD matched agbd_l4a to <0.1% for 10,689 shots in 15 strata.

L4A V3 model as applied in the product:
    x_i = sqrt(rh_l2a[RH_i] + 100);  y = b0 + sum b_i x_i;  y = max(y, 0);  AGBD = bc * y^2
    (no capping; predictor/response maxima only set flags)
EBT_SA and DBT_SA: the product applies a model that is not in ANCILLARY/model_data; the applied
coefficients (recovered exactly from the product) are supplied as overrides.

Drift: correction_m = drift(cover, year) - mean drift over T1 shots, added to every RH
percentile; drift = excess_rate(cover) * (year - ref_fit) + offset_year(cover).
Variants: "central", "low", "high" (rate estimate or 95% CI bound; offsets stay central).
"""
import json

import numpy as np

OFFSET = 100.0


class L4AModels:
    def __init__(self, model_json, overrides_json=None):
        rows = json.load(open(model_json))
        self.m = {}
        for r in rows:
            k = int(r["npar"]) - 1
            self.m[r["predict_stratum"]] = dict(rh=[int(v) for v in r["rh_index"][:k]],
                                                b=np.array(r["par"][:k + 1], float),
                                                bc=float(r["bias_correction_value"]))
        if overrides_json:
            for s, v in json.load(open(overrides_json)).items():
                self.m[s]["b"] = np.array(v["par"], float)
                self.m[s]["bc"] = float(v.get("bc", self.m[s]["bc"]))

    def rh_needed(self):
        return sorted({p for m in self.m.values() for p in m["rh"]})

    def predict(self, stratum, rh, shift):
        """AGBD for shots of all strata. rh: dict percentile -> array (m); shift: array (m)."""
        stratum = np.asarray(stratum, dtype=object)
        out = np.full(len(stratum), np.nan)
        for s in np.unique(stratum[pd_notna(stratum)]):
            if s not in self.m:
                continue
            idx = np.where(stratum == s)[0]
            m = self.m[s]
            y = np.full(len(idx), m["b"][0])
            for coef, p in zip(m["b"][1:], m["rh"]):
                y = y + coef * np.sqrt(np.maximum(rh[p][idx] + shift[idx] + OFFSET, 0.0))
            out[idx] = m["bc"] * np.maximum(y, 0.0) ** 2
        return out


def pd_notna(a):
    return np.array([v is not None and v == v for v in a], bool)


class Drift:
    def __init__(self, table_csv):
        import pandas as pd
        t = pd.read_csv(table_csv).sort_values("cover_mid")
        self.mid = t.cover_mid.to_numpy(float)
        self.rate = {"central": t.excess_rate.clip(lower=0).to_numpy(float),
                     "low": t.excess_lo.clip(lower=0).to_numpy(float),
                     "high": t.excess_hi.clip(lower=0).to_numpy(float)}
        self.offsets = {int(c.split("_")[1]): t[c].fillna(0).to_numpy(float)
                        for c in t.columns if c.startswith("offset_")}
        self.ref = float(t.ref_fit.iloc[0])
        w = json.loads(t.t1_year_weights.iloc[0])
        self.t1_years = np.array([float(k) for k in w])
        self.t1_w = np.array([float(v) for v in w.values()])

    def _drift(self, rate, cover, years):
        d = rate * (years - self.ref)
        for y, off in self.offsets.items():
            d = d + np.where(years == y, np.interp(cover, self.mid, off), 0.0)
        return d

    def correction(self, cover, year, which="central"):
        """Metres to add to every RH percentile; 0 where cover is missing."""
        cover = np.asarray(cover, float)
        ok = np.isfinite(cover)
        c = np.where(ok, cover, 0.0)
        rate = np.interp(c, self.mid, self.rate[which])
        yrs = np.asarray(year, float)
        mean_t1 = sum(self._drift(rate, c, np.full_like(c, y)) * n
                      for y, n in zip(self.t1_years, self.t1_w)) / self.t1_w.sum()
        return np.where(ok, self._drift(rate, c, yrs) - mean_t1, 0.0)
