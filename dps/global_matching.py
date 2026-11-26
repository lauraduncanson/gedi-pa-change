#!/usr/bin/env python
"""
Global matching on DPS (design in CLAUDE.md, "Global matching design"): plan units, submit jobs,
collect matches and run the per-country DiD.

  python dps/global_matching.py plan                 # units to run, parts, skipped units -> output/global_matching_plan.csv
  python dps/global_matching.py submit --units USA_04 BRA_01 ... [--dry-run]   # practice set
  python dps/global_matching.py submit --all [--dry-run]                       # every planned unit not yet submitted
  python dps/global_matching.py status
  python dps/global_matching.py collect              # download per-PA matches of successful jobs
  python dps/global_matching.py did                  # did_country.py per country -> results/global_did/<ISO3>/

Jobs: process 250 (v5), continuous L4B, 32 cores on maap-dps-worker-32vcpu-64gb; units are split into parts
(same matches) by PA cells (PART_CELLS) and by PAs x controls (PART_WORK), whichever needs more.
Job ids go to dps/global_matching_jobs.csv (not committed).
"""
import argparse
import csv
import glob
import math
import os
import subprocess
import sys
import time

import pandas as pd

PROCESS_ID = "PENDING_v6"
QUEUE = "maap-dps-worker-32vcpu-64gb"
PREFIX = "s3://maap-ops-workspace/lduncanson/gedi_pa_change/matching_inputs_global"
PART_CELLS = 40_000
PART_WORK = 2e8  # PAs x complete controls per job: her propensity filter fits one glm per PA on the whole pool
MAX_PARTS = 16
PLAN = "output/global_matching_plan.csv"
JOBS = "dps/global_matching_jobs.csv"
WORK = os.environ.get("GLOBAL_MATCH_WORK", "output/global_matches")


def manifest():
    files = sorted(glob.glob("output/matching_units_manifest_*.csv"))
    m = pd.concat([pd.read_csv(f) for f in files], ignore_index=True).drop_duplicates("unit", keep="last")
    return m


def plan(a):
    m = manifest()
    m["pa_cells_complete"] = m.pa_cells_complete.fillna(0)
    m["controls_complete"] = m.controls_complete.fillna(0)
    m["run"] = m.written & (m.pa_cells_complete > 0) & (m.controls_complete > 0)
    m["thin"] = m.run & (m.controls_complete < m.pa_cells_complete)
    m["parts"] = [min(MAX_PARTS, max(1, math.ceil(c / PART_CELLS), math.ceil(n * k / PART_WORK))) if r else 0
                  for c, n, k, r in zip(m.pa_cells_complete, m.n_pa, m.controls_complete, m.run)]
    m.to_csv(PLAN, index=False)
    r = m[m.run]
    print(f"{len(m)} units in {m.iso3.nunique()} countries; to run: {len(r)} units in {r.iso3.nunique()} countries, "
          f"{int(r.parts.sum())} jobs; {int(m.thin.sum())} thin units (controls < PA cells); "
          f"skipped {int((~m.run).sum())} (no PA cells or no controls)")
    print(f"complete PA cells to match: {int(r.pa_cells_complete.sum()):,}; controls: {int(r.controls_complete.sum()):,}")
    print(r.sort_values("pa_cells_complete", ascending=False)[["unit", "n_pa", "pa_cells_complete", "controls_complete", "parts"]].head(15).to_string(index=False))


def submitted():
    if not os.path.exists(JOBS):
        return pd.DataFrame(columns=["unit", "part", "job_id", "http_status", "submitted_at"])
    return pd.read_csv(JOBS, dtype=str)


def submit(a):
    p = pd.read_csv(PLAN)
    p = p[p.run]
    if not a.all:
        p = p[p.unit.isin(a.units)]
        missing = set(a.units) - set(p.unit)
        if missing:
            sys.exit(f"not in plan (or nothing to run): {sorted(missing)}")
    done = submitted()
    ok = done[done.job_id.notna() & (done.job_id != "")]
    have = set(zip(ok.unit, ok.part))
    have_units = set(ok.unit)  # never resubmit a unit under a different part split (would duplicate PAs)
    todo = [(u, f"{i}/{n}") for u, n in zip(p.unit, p.parts) for i in range(1, int(n) + 1)
            if (u, f"{i}/{n}") not in have and u not in have_units]
    print(f"submitting {len(todo)} job(s) ({len({u for u, _ in todo})} units)")
    if a.dry_run:
        print(todo[:20]); return
    from maap.maap import MAAP
    maap = MAAP()
    new = not os.path.exists(JOBS)
    with open(JOBS, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["unit", "part", "job_id", "http_status", "submitted_at"])
        for k, (u, part) in enumerate(todo, 1):
            r = maap.submit_job(process_id=PROCESS_ID, queue=QUEUE, tag=f"gedi-pa-change-global-{u}-{part.replace('/', 'of')}",
                                inputs={"region": u, "covariates": "with_l4b", "n_cores": "32", "part": part, "input_prefix": PREFIX})
            body = r.json() if r.status_code in (200, 201, 202) else {}
            jid = body.get("jobID") or body.get("id") or ""
            if not jid:
                print(f"  {u} {part}: http {r.status_code} {r.text[:150]}")
            w.writerow([u, part, jid, r.status_code, time.strftime("%Y-%m-%dT%H:%M:%S")]); f.flush()
            if k % 50 == 0:
                print(f"  {k}/{len(todo)}", flush=True)
            time.sleep(0.3)


def statuses(rows):
    from maap.maap import MAAP
    maap = MAAP()
    out = {}
    for j in rows.job_id:
        try:
            out[j] = maap.get_job_status(j).json().get("status", "unknown")
        except Exception:
            out[j] = "unknown"
    return out


def status(a):
    j = submitted()
    j = j[j.job_id.notna() & (j.job_id != "")].drop_duplicates(["unit", "part"], keep="last")
    s = statuses(j)
    j["status"] = j.job_id.map(s)
    print(j.status.value_counts().to_dict())
    bad = j[j.status.isin(["failed", "dismissed"])]
    if len(bad):
        print("failed:", " ".join(f"{u}:{p}" for u, p in zip(bad.unit, bad.part)))
    return j


def collect(a):
    import boto3
    from maap.maap import MAAP
    maap, s3 = MAAP(), boto3.client("s3")
    j = status(a)
    j = j[j.status == "successful"]
    n = 0
    for _, r in j.iterrows():
        dst = os.path.join(WORK, r.unit[:3], r.unit)
        marker = os.path.join(dst, f".done_{r.part.replace('/', 'of')}")
        if os.path.exists(marker):
            continue
        os.makedirs(dst, exist_ok=True)
        res = maap.get_job_result(r.job_id).json()
        pre = next(l["href"] for v in res.values() if isinstance(v, dict) for l in v.get("links", [])
                   if l["href"].startswith("s3://")).split("maap-ops-workspace/")[1]
        for page in s3.get_paginator("list_objects_v2").paginate(Bucket="maap-ops-workspace", Prefix=f"{pre}/MATCHING_results_V3/{r.unit}/"):
            for o in page.get("Contents", []):
                s3.download_file("maap-ops-workspace", o["Key"], os.path.join(dst, o["Key"].rsplit("/", 1)[1])); n += 1
        open(marker, "w").close()
    print(f"downloaded {n} per-PA files into {WORK}/")


def did(a):
    rows = []
    for iso_dir in sorted(glob.glob(os.path.join(WORK, "???"))):
        iso = os.path.basename(iso_dir)
        out = os.path.join("results", "global_did", iso)
        os.makedirs(out, exist_ok=True)
        r = subprocess.run([sys.executable, "analysis/did_country.py", "--iso3", iso, "--matches", iso_dir,
                            "--cells", "output/gedi_pa_change_cells_combined.parquet", "--out", out],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"{iso}: did_country FAILED: {r.stderr[-300:]}"); continue
        s = pd.read_csv(os.path.join(out, f"{iso}_did_summary.csv")); s.insert(0, "iso3", iso); rows.append(s)
        print(f"{iso}: done", flush=True)
    if rows:
        allc = pd.concat(rows, ignore_index=True)
        allc.to_csv(os.path.join("results", "global_did", "summary_by_country.csv"), index=False)
        print(f"wrote results/global_did/summary_by_country.csv ({allc.iso3.nunique()} countries)")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan")
    s = sub.add_parser("submit"); s.add_argument("--units", nargs="*", default=[]); s.add_argument("--all", action="store_true")
    s.add_argument("--dry-run", action="store_true")
    sub.add_parser("status"); sub.add_parser("collect"); sub.add_parser("did")
    a = ap.parse_args()
    {"plan": plan, "submit": submit, "status": status, "collect": collect, "did": did}[a.cmd](a)


if __name__ == "__main__":
    main()
