# Source of these R files

Veronika Leitold's Brazil V3 matching step (step 4), copied so it can run on DPS as the
`gedi-pa-change-match-brazil` algorithm (`matching/dps/`). The first commit of this folder holds
byte-identical copies; later commits are our patches (see `git log -p matching/r/brazil_v3`).

| File | Copied from (`s3://maap-ops-workspace/shared/leitoldv/GEDI_global_change/` = Hub `shared-buckets/leitoldv/GEDI_global_change/`) | Last modified | MD5 |
|---|---|---|---|
| `gedi_change_step4_V3_cmd.R` | `gedi_change_brazil/repo/matching/gedi_change_step4_V3_cmd.R` | 2026-11-25 00:18 UTC | 095f9be5560c36f75b880fbb332a74d1 |
| `matching_func_2025.R` | `matching_func_2025.R` | 2025-10-30 18:46 UTC | b25464a3e9881093f9047ce7b58bda45 |

Copied 2026-11-26. Her V3 matches (`MATCHING_results_V3/`) date from 2026-10-06/07 (per region)
and 2026-10-30 (`BRA_allPA_matched_results_V3.RDS`), so `gedi_change_step4_V3_cmd.R` was edited
after they were made; the with-L4B DPS run checks whether it reproduces them.

Inputs (unchanged, read from her folder): `MATCHING_points_V3/<region>_prepped_PAs.RDS` and
`<region>_prepped_control.RDS` for the Brazil biome regions BAe, BAw, BCa, BCe, BMa, BPp, BPt.
