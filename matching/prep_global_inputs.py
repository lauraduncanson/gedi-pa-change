#!/usr/bin/env python
"""
Convert a country's 2025 global-pipeline matching inputs (steps 1-3) into the format of
Veronika's Brazil V3 step 4 (matching/r/brazil_v3/), adding the GEDI L4B covariate.

Reads (read only):
  s3://maap-ops-workspace/shared/abarenblitt/GEDI_global_PA_v2/Matching_Results/<ISO3>/
      <ISO3>_testPAs/prepped_pa_<WDPAID>_wk24.RDS      one data frame per PA
      WDPA_grids/<ISO3>_prepped_control_wk24.RDS        control pool
Writes:
  <out>/<ISO3>_prepped_PAs.RDS      one data frame, all PAs, with a pa_id column (the R step 4
                                    splits it by pa_id)
  <out>/<ISO3>_prepped_control.RDS
Added columns: gedi_l4b (sampled at the cell centre from --l4b, the raster used for Brazil),
l4b_class (factor, floor(L4B / 50) with >= 400 Mg/ha merged), wwfbiom = wwfecoreg = 0 (absent
in these inputs; in her code they only pre-filter controls, so a constant disables that).

  python matching/prep_global_inputs.py --iso3 GHA --l4b <local gedi_l4b.tif> \
      --out s3://maap-ops-workspace/lduncanson/gedi_pa_change/matching_inputs
"""
import argparse
import os
import re
import tempfile

import boto3
import numpy as np
import pandas as pd
import pyreadr
import rasterio

BUCKET = "maap-ops-workspace"
SRC = "shared/abarenblitt/GEDI_global_PA_v2/Matching_Results/{iso3}/"
KEEP = ["lat", "lon", "UID", "status", "land_cover", "wwfbiom", "wwfecoreg", "elevation", "slope",
        "mean_temp", "max_temp", "min_temp", "prec", "d2road", "d2city", "tt2city", "popden",
        "popcnt", "gedi_l4b", "l4b_class"]
CLASS_WIDTH, CLASS_TOP = 50.0, 400.0


def add_l4b(df, l4b_path):
    with rasterio.open(l4b_path) as src:
        assert src.crs.to_epsg() == 4326, src.crs
        v = np.array([x[0] for x in src.sample(zip(df.lon.values, df.lat.values))], dtype="float64")
        if src.nodata is not None:
            v[v == src.nodata] = np.nan
    v[v < 0] = np.nan
    df["gedi_l4b"] = v
    k = np.floor(np.minimum(v, CLASS_TOP) / CLASS_WIDTH)
    labels = [f"c{int(i):02d}" for i in range(int(CLASS_TOP / CLASS_WIDTH) + 1)]
    df["l4b_class"] = pd.Categorical([labels[int(i)] if np.isfinite(i) else None for i in k], categories=labels)
    df["wwfbiom"] = 0.0
    df["wwfecoreg"] = 0.0
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iso3", required=True)
    ap.add_argument("--l4b", required=True, help="local copy of gedi_l4b.tif (EPSG:4326)")
    ap.add_argument("--out", required=True, help="s3:// prefix or local folder")
    a = ap.parse_args()
    s3 = boto3.client("s3")
    src = SRC.format(iso3=a.iso3)
    tmp = tempfile.mkdtemp()

    ctl_path = os.path.join(tmp, "control.RDS")
    s3.download_file(BUCKET, f"{src}WDPA_grids/{a.iso3}_prepped_control_wk24.RDS", ctl_path)
    ctl = pyreadr.read_r(ctl_path)[None]
    pas = []
    for page in s3.get_paginator("list_objects_v2").paginate(Bucket=BUCKET, Prefix=f"{src}{a.iso3}_testPAs/"):
        for o in page.get("Contents", []):
            m = re.search(r"prepped_pa_(\d+)_wk24\.RDS$", o["Key"])
            if m:
                p = os.path.join(tmp, os.path.basename(o["Key"]))
                s3.download_file(BUCKET, o["Key"], p)
                d = pyreadr.read_r(p)[None]
                d["pa_id"] = m.group(1)
                pas.append(d)
    pa = pd.concat(pas, ignore_index=True)

    ctl = add_l4b(ctl, a.l4b)[KEEP]
    pa = add_l4b(pa, a.l4b)[KEEP + ["pa_id"]]
    for d in (ctl, pa):
        d["status"] = d.status.astype(bool)
        d["UID"] = d.UID.astype("int32")
    assert ctl.UID.is_unique, "control UIDs must be unique"
    assert not pa.duplicated(["pa_id", "UID"]).any(), "UIDs must be unique within each PA"
    assert ctl.status.eq(False).all() and pa.status.eq(True).all()

    for name, d in ((f"{a.iso3}_prepped_control.RDS", ctl), (f"{a.iso3}_prepped_PAs.RDS", pa)):
        local = os.path.join(tmp, name)
        pyreadr.write_rds(local, d, compress="gzip")
        if a.out.startswith("s3://"):
            b, k = a.out[5:].split("/", 1)
            s3.upload_file(local, b, f"{k.rstrip('/')}/{a.iso3}/{name}")
            dest = f"{a.out.rstrip('/')}/{a.iso3}/{name}"
        else:
            os.makedirs(os.path.join(a.out, a.iso3), exist_ok=True)
            dest = os.path.join(a.out, a.iso3, name)
            os.replace(local, dest)
        print(f"wrote {dest}: {len(d):,} rows")
    for label, d in (("control", ctl), ("PA", pa)):
        print(f"{a.iso3} {label}: L4B missing {d.gedi_l4b.isna().mean():.1%}, rows with any NA "
              f"{d.isna().any(axis=1).mean():.1%}; L4B classes {d.l4b_class.value_counts().sort_index().to_dict()}")
    print(f"{a.iso3}: {pa.pa_id.nunique()} PAs, {len(pa):,} PA cells, {len(ctl):,} control cells")


if __name__ == "__main__":
    main()
