"""Shared helpers for submitting and tracking the per-tile DPS jobs."""
import glob
import os
import re

MANIFEST = "s3://nasa-maap-data-store/file-staging/nasa-map/gedi-tiled-v3/manifest.txt"
YEARS = {2019, 2020, 2021, 2024, 2025, 2026}
# Where DPS job outputs appear on the Hub. Finished tiles are found by their summary_<tile>.json.
DPS_OUTPUT = "/home/jovyan/my-private-bucket/dps_output"
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


def done_tiles(root=DPS_OUTPUT):
    files = glob.glob(os.path.join(root, "**", "summary_*.json"), recursive=True)
    return {os.path.basename(f)[len("summary_"):-len(".json")] for f in files}
