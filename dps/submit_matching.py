#!/usr/bin/env python
"""
Submit Brazil matching jobs (algorithm gedi-pa-change-match-brazil), one per region x covariate set.

  python dps/submit_matching.py --process-id <id> --queue maap-dps-worker-16gb \
      --regions BPp --covariates no_l4b --n-cores 4 --dry-run
  ... --parts 2      # split each region into 2 jobs (parts 1/2 and 2/2; same matches)

Job ids are appended to dps/match_jobs.csv (not committed). Regions: BAe BAw BCa BCe BMa BPp BPt.
"""
import argparse
import csv
import os
import time

REGIONS = ["BAe", "BAw", "BCa", "BCe", "BMa", "BPp", "BPt"]
JOBS_CSV = "dps/match_jobs.csv"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--process-id", required=True)
    ap.add_argument("--queue", required=True)
    ap.add_argument("--regions", nargs="+", required=True, choices=REGIONS)
    ap.add_argument("--covariates", nargs="+", required=True, choices=["with_l4b", "no_l4b"])
    ap.add_argument("--n-cores", default="4")
    ap.add_argument("--parts", type=int, default=1, help="jobs per region x covariate set")
    ap.add_argument("--only-part", type=int, help="submit just this part (1-based) of --parts")
    ap.add_argument("--tag", default="gedi-pa-change-match-brazil")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    parts = [a.only_part] if a.only_part else range(1, a.parts + 1)
    todo = [(r, c, f"{i}/{a.parts}") for r in a.regions for c in a.covariates for i in parts]
    print(f"submitting {len(todo)} job(s) to {a.queue} with n_cores={a.n_cores}: {todo}")
    if a.dry_run:
        return
    from maap.maap import MAAP
    maap = MAAP()
    new = not os.path.exists(JOBS_CSV)
    with open(JOBS_CSV, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["region", "covariates", "n_cores", "queue", "job_id", "http_status", "submitted_at", "part"])
        for region, cov, part in todo:
            r = maap.submit_job(process_id=a.process_id, queue=a.queue, tag=f"{a.tag}-{region}-{cov}-{part.replace('/', 'of')}",
                                inputs={"region": region, "covariates": cov, "n_cores": str(a.n_cores), "part": part})
            body = r.json() if r.status_code in (200, 201, 202) else {}
            job_id = body.get("jobID") or body.get("id") or ""
            print(f"  {region} {cov} {part}: http {r.status_code} job {job_id or r.text[:200]}")
            w.writerow([region, cov, a.n_cores, a.queue, job_id, r.status_code, time.strftime("%Y-%m-%dT%H:%M:%S"), part])
            f.flush()
            time.sleep(0.5)


if __name__ == "__main__":
    main()
