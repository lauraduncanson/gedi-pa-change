#!/usr/bin/env python
"""
Find the L4A biomass model parameters stored inside a GEDI L4A granule and export them.

Every L4A granule carries the full model table for all prediction strata (coefficients,
predictor RH percentiles, transforms, bias correction), as far as we know in the ANCILLARY
group. This script does not rely on that: it walks the whole file (skipping the per-beam shot
groups) and reports every non-beam dataset, flagging names that look model-related, then
exports anything under ANCILLARY to CSV/JSON.

Which granule: by default it takes a source L4A granule name from the tiled data
(root_file_l4a), so the models match the version behind the tiled dataset. It then opens the
granule from its DAAC's S3 bucket (found by exact file name in the public CMR catalogue) with
MAAP's own Earthdata credentials (no login prompt).
You can instead pass a local path or s3:// URL with --granule.

  python analysis/find_l4a_models.py --out /home/jovyan/GEDI/results/l4a_models
  python analysis/find_l4a_models.py --granule /path/to/GEDI04_A_....h5 --out ...
"""
import argparse
import json
import os
import sys

import h5py
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "extract_tile"))

MODEL_WORDS = ("model", "par", "predict", "rh_index", "transform", "bias", "vcov", "coef",
               "stratum", "lut", "rse", "npar")


def granule_name_from_tiles(tile="S13_E028"):
    import extract_tile as E
    con = E.connect()
    files = E.tile_files(con, f"{E.ROOT}/manifest.txt", tile, [2020])
    if not files:
        raise SystemExit(f"tile {tile} not in manifest; pass --tile or --granule")
    r = con.sql(f"SELECT DISTINCT root_file_l4a FROM read_parquet('{files[0]}') "
                f"WHERE root_file_l4a IS NOT NULL LIMIT 3").df()
    names = r.root_file_l4a.tolist()
    print("root_file_l4a examples:", names)
    return names[0]


S3_CRED_ENDPOINTS = {  # bucket keyword -> DAAC S3 credentials endpoint
    "ornl": "https://data.ornldaac.earthdata.nasa.gov/s3credentials",
    "lp-": "https://data.lpdaac.earthdatacloud.nasa.gov/s3credentials",
    "lpdaac": "https://data.lpdaac.earthdatacloud.nasa.gov/s3credentials",
}


CMR = "https://cmr.earthdata.nasa.gov/search"


def l4a_collections():
    """GEDI L4A collections in CMR (any version), as (concept_id, short_name, version)."""
    import requests
    found = {}
    for params in ({"short_name": "GEDI*L4A*", "options[short_name][pattern]": "true"},
                   {"keyword": "GEDI04_A"}, {"keyword": "GEDI L4A"}):
        r = requests.get(f"{CMR}/collections.json", params={**params, "page_size": 50}, timeout=60)
        if not r.ok:
            print(f"  collection search {params} -> {r.status_code}: {r.text[:200]}")
            continue
        for e in r.json().get("feed", {}).get("entry", []):
            found[e["id"]] = (e["id"], e.get("short_name"), e.get("version_id"))
    return list(found.values())


def _s3_links(item):
    urls = [u.get("URL", "") for u in item.get("umm", {}).get("RelatedUrls", [])]
    return [u for u in urls if u.startswith("s3://") and u.endswith(".h5")], urls


def s3_url_from_cmr(name):
    """Find a version-matched L4A granule (no login).

    The model table is the same in every granule of a version, so any granule from the matching
    collection will do. Tries the named granule first (pattern match), then any granule.
    """
    import requests
    colls = l4a_collections()
    ver = name.rsplit("_V", 1)[-1].split(".")[0].lstrip("0") if "_V" in name else ""
    cands = [c for c in colls if "L4A_AGB_Density" in str(c[1]) and "GW" not in str(c[1])
             and str(c[2]).split(".")[0] == ver] or \
            [c for c in colls if "L4A_AGB_Density" in str(c[1]) and "GW" not in str(c[1])]
    print("  L4A collections in CMR:", [(c[1], c[2]) for c in colls])
    print("  using:", [(c[1], c[2]) for c in cands])
    stem = name[:-3] if name.endswith(".h5") else name
    for cid, sn, v in cands:
        for params in ({"readable_granule_name": f"*{stem}*",
                        "options[readable_granule_name][pattern]": "true"}, {}):
            r = requests.get(f"{CMR}/granules.umm_json",
                             params={"collection_concept_id": cid, "page_size": 5, **params}, timeout=60)
            if not r.ok:
                print(f"  granule search in {sn} -> {r.status_code}: {r.text[:200]}")
                continue
            for it in r.json().get("items", []):
                s3, urls = _s3_links(it)
                if s3:
                    which = "named granule" if params else "first granule of the collection"
                    print(f"  {which} in {sn} v{v}: {s3[0].rsplit('/', 1)[-1]}")
                    return s3[0], None
    raise SystemExit("no L4A granule with an s3:// link found; pass --granule with a path or s3:// URL")


def open_granule(ref):
    """Return an open h5py.File for a local path, s3:// URL, or bare granule name."""
    if os.path.exists(ref):
        return h5py.File(ref, "r")
    if not ref.startswith("s3://"):
        ref, _ = s3_url_from_cmr(os.path.basename(ref))
    print("opening", ref)
    bucket = ref[5:].split("/", 1)[0]
    endpoint = next((v for k, v in S3_CRED_ENDPOINTS.items() if k in bucket), None)
    if endpoint is None:
        raise SystemExit(f"don't know the credentials endpoint for bucket {bucket}")
    from maap.maap import MAAP
    maap = MAAP()
    c = maap.aws.earthdata_s3_credentials(endpoint)
    import s3fs
    fs = s3fs.S3FileSystem(key=c["accessKeyId"], secret=c["secretAccessKey"], token=c["sessionToken"])
    return h5py.File(fs.open(ref, "rb"), "r")


def walk(f):
    """List all datasets outside the per-beam groups."""
    rows = []

    def visit(name, obj):
        if name.split("/")[0].startswith("BEAM"):
            return
        if isinstance(obj, h5py.Dataset):
            fields = list(obj.dtype.names) if obj.dtype.names else []
            hit = any(w in name.lower() for w in MODEL_WORDS) or any(
                any(w in fl.lower() for w in MODEL_WORDS) for fl in fields)
            rows.append(dict(path=name, shape=str(obj.shape), dtype=str(obj.dtype)[:60],
                             fields=", ".join(fields), model_related=hit))
    f.visititems(visit)
    return pd.DataFrame(rows)


def to_python(v):
    if isinstance(v, bytes):
        return v.decode(errors="replace").strip("\x00 ")
    if isinstance(v, np.ndarray):
        return [to_python(x) for x in v.tolist()] if v.dtype.kind in "SO" else v.tolist()
    if isinstance(v, np.generic):
        return v.item()
    return v


def export(ds, out, stem):
    arr = ds[()]
    if ds.dtype.names:  # compound table: one row per record, array fields kept as lists
        recs = [{k: to_python(r[k]) for k in ds.dtype.names} for r in np.atleast_1d(arr)]
        with open(os.path.join(out, f"{stem}.json"), "w") as fh:
            json.dump(recs, fh, indent=1, default=str)
        pd.DataFrame([{k: (json.dumps(v) if isinstance(v, list) else v) for k, v in r.items()}
                      for r in recs]).to_csv(os.path.join(out, f"{stem}.csv"), index=False)
        return len(recs)
    val = to_python(arr)
    with open(os.path.join(out, f"{stem}.json"), "w") as fh:
        json.dump(val, fh, indent=1, default=str)
    return np.size(arr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--granule", help="local path, s3:// URL or granule name (default: from tiled data)")
    ap.add_argument("--tile", default="S13_E028", help="tile used to look up a source granule name")
    ap.add_argument("--out", default="results/l4a_models")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    ref = a.granule or granule_name_from_tiles(a.tile)
    with open_granule(ref) as f:
        tree = walk(f)
        tree.to_csv(os.path.join(a.out, "l4a_nonbeam_datasets.csv"), index=False)
        with pd.option_context("display.width", 200, "display.max_colwidth", 90):
            print("\nnon-beam datasets in the granule (model_related = name or fields look like model data):")
            print(tree.to_string(index=False))
        targets = tree[tree.path.str.startswith("ANCILLARY") | tree.model_related]
        print()
        for p in targets.path:
            n = export(f[p], a.out, p.replace("/", "__"))
            print(f"exported {p} ({n} records/values)")
        meta = {k: to_python(v) for k, v in f.attrs.items()}
        with open(os.path.join(a.out, "file_attributes.json"), "w") as fh:
            json.dump(meta, fh, indent=1, default=str)
    print(f"\nwrote {a.out}/ (dataset list, exported tables, file attributes)")
    print("The model table is usually the compound dataset with one row per prediction stratum "
          "(fields like par, rh_index / predictor_id, x_transform, y_transform, bias_correction_*).")


if __name__ == "__main__":
    main()
