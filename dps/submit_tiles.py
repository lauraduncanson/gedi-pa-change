#!/usr/bin/env python
"""
Submit one DPS job per tile, skipping tiles that already have output.

  python dps/submit_tiles.py --process-id 123 --queue <queue> --tiles S13_E028 N11_W005   # test
  python dps/submit_tiles.py --process-id 123 --queue <queue> --limit 100                  # ramp up
  python dps/submit_tiles.py --process-id 123 --queue <queue>                              # everything left

Safe to re-run: finished tiles are skipped, and tiles already submitted (in dps/jobs.csv)
are skipped unless --resubmit is given. Use --dry-run to see what would be submitted.
"""
import argparse
import csv
import os
import time

import common as C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--process-id", required=True, help="MAAP process id of the registered algorithm")
    ap.add_argument("--queue", required=True)
    ap.add_argument("--tag", default="gedi-pa-change-v1")
    ap.add_argument("--tiles", nargs="*", help="specific tile ids (default: all tiles in the manifest)")
    ap.add_argument("--limit", type=int, help="submit at most this many")
    ap.add_argument("--manifest", default=C.MANIFEST)
    ap.add_argument("--done-dir", default=C.DPS_OUTPUT)
    ap.add_argument("--resubmit", action="store_true", help="also submit tiles already in jobs.csv")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    tiles = a.tiles or C.all_tiles(a.manifest)
    done = C.done_tiles(a.done_dir)
    submitted = set()
    if os.path.exists(C.JOBS_CSV) and not a.resubmit:
        with open(C.JOBS_CSV) as f:
            submitted = {r["tile"] for r in csv.DictReader(f)}
    todo = [t for t in tiles if t not in done and t not in submitted]
    if a.limit:
        todo = todo[:a.limit]
    print(f"{len(tiles)} tiles; {len(done)} done; {len(submitted)} already submitted; submitting {len(todo)}")
    if a.dry_run or not todo:
        print("first few:", todo[:10])
        return

    from maap.maap import MAAP
    maap = MAAP()
    new = not os.path.exists(C.JOBS_CSV)
    os.makedirs(os.path.dirname(C.JOBS_CSV), exist_ok=True)
    with open(C.JOBS_CSV, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["tile", "job_id", "http_status", "submitted_at"])
        for i, t in enumerate(todo, 1):
            r = maap.submit_job(process_id=a.process_id, inputs={"tile": t}, queue=a.queue,
                                dedup=True, tag=a.tag)
            job_id = r.json().get("id") if r.status_code == 202 else ""
            if not job_id:
                print(f"  {t}: submission failed ({r.status_code}) {r.text[:200]}")
            w.writerow([t, job_id, r.status_code, time.strftime("%Y-%m-%dT%H:%M:%S")])
            f.flush()
            if i % 50 == 0:
                print(f"  submitted {i}/{len(todo)}")
            time.sleep(0.2)
    print(f"done; job ids recorded in {C.JOBS_CSV}")


if __name__ == "__main__":
    main()
