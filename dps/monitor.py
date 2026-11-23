#!/usr/bin/env python
"""
Report where the per-tile jobs stand, using output files as the source of truth.

  python dps/monitor.py                 # counts: done / still to do, plus job statuses
  python dps/monitor.py --list-missing  # print tiles with a job but no output yet

A tile is done when its summary_<tile>.json exists. To retry tiles that failed, re-run
submit_tiles.py with --resubmit: it skips finished tiles and resubmits the rest.
"""
import argparse
import collections
import csv
import os

import common as C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--done-dir", default=C.DPS_OUTPUT)
    ap.add_argument("--list-missing", action="store_true")
    ap.add_argument("--no-status", action="store_true", help="skip querying MAAP for job statuses")
    a = ap.parse_args()

    done = C.done_tiles(a.done_dir)
    jobs = {}
    if os.path.exists(C.JOBS_CSV):
        with open(C.JOBS_CSV) as f:
            for r in csv.DictReader(f):
                jobs[r["tile"]] = r["job_id"]  # latest submission wins
    missing = sorted(t for t in jobs if t not in done)
    print(f"submitted: {len(jobs)}   done: {len(done & set(jobs))}   not done yet: {len(missing)}")

    if missing and not a.no_status:
        from maap.maap import MAAP
        maap = MAAP()
        counts = collections.Counter()
        for t in missing:
            if not jobs[t]:
                counts["submission failed"] += 1
                continue
            r = maap.get_job_status(jobs[t])
            counts[r.json().get("status", "unknown") if r.status_code == 200 else f"http {r.status_code}"] += 1
        print("status of jobs not done yet:", dict(counts))
    if a.list_missing:
        print("\n".join(missing))


if __name__ == "__main__":
    main()
