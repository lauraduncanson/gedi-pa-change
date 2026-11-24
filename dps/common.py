"""Shared helpers for submitting and tracking the per-tile DPS jobs."""
import glob
import os
import re

MANIFEST = "s3://nasa-maap-data-store/file-staging/nasa-map/gedi-tiled-v3/manifest.txt"
YEARS = {2019, 2020, 2021, 2024, 2025, 2026}
# Where DPS job outputs are. Finished tiles are found by their summary_<tile>.json.
# Read straight from S3 by default: the Hub bucket mount (my-private-bucket) can drop
# ("Transport endpoint is not connected"), and listing S3 is much faster than globbing the mount.
# The prefix matches every deployment (gedi-pa-change-extract-tile_14/main/..., _15, ...).
# A local folder or glob pattern also works, e.g.
#   /home/jovyan/my-private-bucket/dps_output/gedi-pa-change-extract-tile*
DPS_OUTPUT = "s3://maap-ops-workspace/lduncanson/dps_output/gedi-pa-change-extract-tile"
JOBS_CSV = "dps/jobs.csv"


def read_manifest(path=MANIFEST):
    if path.startswith("s3://"):
        import boto3
        bucket, key = path[5:].split("/", 1)
        return boto3.client("s3").get_object(Bucket=bucket, Key=key)["Body"].read().decode().splitlines()
    with open(path) as f:
        return f.read().splitlines()


def all_tiles(manifest=MANIFEST):
    """Tile ids that have data in at least one T1 or T2 year."""
    tiles = set()
    for line in read_manifest(manifest):
        m = re.search(r"/tile_id=([^/]+)/year=(\d{4})/", line)
        if m and int(m.group(2)) in YEARS:
            tiles.add(m.group(1))
    return sorted(tiles)


def output_files(root=DPS_OUTPUT, prefix="summary_", suffix=".json"):
    """Paths of DPS output files named <prefix>*<suffix> under root (an s3:// prefix or a local glob)."""
    if root.startswith("s3://"):
        import boto3
        bucket, key = root[5:].split("/", 1)
        out = []
        for page in boto3.client("s3").get_paginator("list_objects_v2").paginate(Bucket=bucket, Prefix=key):
            for o in page.get("Contents", []):
                name = o["Key"].rsplit("/", 1)[-1]
                if name.startswith(prefix) and name.endswith(suffix):
                    out.append(f"s3://{bucket}/{o['Key']}")
        return out
    return glob.glob(os.path.join(root, "**", f"{prefix}*{suffix}"), recursive=True)


def done_tiles(root=DPS_OUTPUT):
    return {os.path.basename(f)[len("summary_"):-len(".json")] for f in output_files(root)}
