# Brazil first look: Veronika's V3 pipeline vs ours

2026-11-26. Read-only review of Veronika Leitold's Brazil pipeline, plus the status of our own
Brazil run. Nothing under `shared-buckets/leitoldv/` was modified.

Abbreviations: `G` = `/home/jovyan/shared-buckets/leitoldv/GEDI_global_change/`,
`M` = `G/gedi_change_brazil/repo/matching/`. Notebook cell numbers are 0-based.

## 0. Status of the requested steps

| Step | Status |
|---|---|
| 1. Report on Step123/4/5/6 and the TMF notebook | Done (sections 1-2). Key facts spot-checked against the source. |
| 2. Differences from our `extract_tile.py` / `analysis/did_country.py` | Done for `extract_tile.py` (section 3). **`analysis/did_country.py` does not exist**, in this repo, its history, or `/home/jovyan/GEDI`. |
| 3. `dps/monitor.py`, then `dps/collect.py` | `monitor.py`: 12,080 / 12,080 tiles done. `collect.py` was **not re-run**: it already ran on 2026-11-24 21:13 UTC, after the newest tile output (19:30 UTC), and `combine_cells.py` ran after it. A re-run would reproduce the same table. |
| 4. `analysis/did_country.py` and `analysis/rtm_check.py` for BRA | **Not run. Neither script exists, and there is no `matching/README.md`**, here or in `M`. The inputs are ready (section 5). |
| 5. This file | Done. Section 6 compares the two pipelines' inputs and methods. We have **no Brazil DiD numbers of our own yet** to put against her CSV. |

## 1. Veronika's V3 pipeline (Steps 1-6)

The run unit (`iso3`) is a Brazilian biome region, not the country:
- `BAe` / `BAw`: Amazon east / west
- `BCa`: Caatinga
- `BCe`: Cerrado
- `BMa`: Mata Atlântica
- `BPp`: Pampa
- `BPt`: Pantanal

Each region has its own boundary file (`G/INPUT_countries/shp/<code>.shp`).

### Step 1-3: prepare data (`M/gedi_change_step123b_V3[_power]_cmd.R`)

**Grid**
- `G/BRA_template_rast.RDS`: EPSG:6933, 1000.895 m, 4351 cols × 4732 rows.
- Extent: x -7141386 to -2786492, y -4065636 to 670599.7.
- This lies exactly on the global GEDI EASE2 1 km grid, at column offset 10217 and row offset 6638.
- Her cell ids are terra cell numbers `t` (1-based) in this Brazil crop. In our ids:
  `cell_id = ((t-1) %/% 4351 + 6638) * 34704 + ((t-1) %% 4351 + 10217)`.

**Eligible cells**
- A cell needs at least 2 shots **and** at least 2 tracks **in both T1 and T2** (l.368-372).
- Shots come from `s3://maap-ops-workspace/shared/ameliah/gedi-test/brazil_tiles/`, using the V3b filters in Step 5.

**Treatment (PA) cells**
- PA file: `G/INPUT_shapefiles/WDPA_Nov2025_BRA_poly_v2_clean.RDS`. The code applies **no IUCN, STATUS or STATUS_YR filter**; what "clean" did upstream is unknown. Indigenous territories are included.
- A cell counts as treatment only if at least 99.9% of its area lies inside one PA.

**Control cells**
- The cell centre must be more than **10 km** from every PA (union of all PAs).
- Controls come only from the same biome region.

**Covariates**
- Point-extracted from `G/INPUT_covariates_2020/`: d2road, d2city, elevation, slope, popcnt, popden, tt2city, prec, mean/max/min temperature (WorldClim 2010-2018), wwfbiom, wwfecoreg, land_cover (MapBiomas 2020, factor), **gedi_l4b**.
- Rows with any NA are dropped.

**Outputs**
- `G/MATCHING_points_V3[_power]/<code>_prepped_{PAs,control}.RDS`.

**`_power` vs plain**
- The only real difference is that `_power` keeps power beams (l.99, 158).
- The `_power` points were regenerated on Nov 25-26 but **have not been matched**.

### Step 4: matching (`M/gedi_change_step4_V3_cmd.R` + `G/matching_func_2025.R`)

**Propensity pre-filter** (`matching_func_2025.R:168`)
- `glm(status ~ mean_temp + max_temp + min_temp + prec + elevation + slope + d2road + d2city + popden + popcnt + tt2city + gedi_l4b)`.
- There is **no `family=`**, so this is a Gaussian linear-probability model, not a logistic one.
- Controls outside the range of the treatment scores are dropped.

**Matching** (`matching_func_2025.R:49-84`)
- Same formula. **GEDI L4B (`gedi_l4b`) is a covariate.**
- `optmatch::match_on` on the glm, then `fullmatch(min.controls = 1, max.controls = 1)`: 1:1 optimal pair matching. **No caliper**; the line is commented out at l.77.
- Mahalanobis distance is used when a PA has 2 or fewer treatment cells.

**Exact matching**
- In practice only `land_cover`, via `strata(land_cover)`.
- `strata(wwfbiom)` and `strata(wwfecoreg)` are in the code, but both columns are numeric. Their `nlevels()` is 0, so these strata are never added. Biome and ecoregion act only as a pre-filter (controls must share values that occur in the PA).

**Loop**
- Per PA, with treatment cells in random chunks of about 300.
- Controls are subsampled at up to 7 times the number of treatment cells.
- Within a PA, a matched control is removed from the pool. Across PAs, **controls can be reused**.

**Saved files**
- Per PA: `G/MATCHING_results_V3/<code>/<code>_pa_<WDPAID>_matching_results.RDS`.
- Per region: `G/MATCHING_results_V3/<code>_allPA_matched_results.RDS`.
- All of Brazil: **`G/MATCHING_results_V3/BRA_allPA_matched_results_V3.RDS`**, which we inspected:
  - 740,868 rows, 24 columns. 370,434 are PA (`status` TRUE) and 370,434 are controls.
  - 1,147 PAs, 2,045 `matched` groups, every group with both PA and control cells.
  - **There is no cell-id column.** Cells are identified by `lat`/`lon`, which sit exactly at EASE2 cell centres. `UID` is a row index that restarts per PA, not a cell id.
  - 538,587 distinct cells. Controls are reused across PAs, and some cells fall in overlapping PAs.
- There is **no power version of Step 4**: `MATCHING_results_V3_power` does not exist.

### Step 5: GEDI extraction (`M/gedi_change_step5_V3[_power].ipynb`)

**Source**
- `ameliah/gedi-test/brazil_tiles` parquet, a pre-joined Brazil test set. It is **not gedi-tiled-v3**.

**Epochs** (cell 2)
- T1 = 2019, 2020, 2021.
- **T2 = 2024, 2025.**

**Filters**
- SQL pre-screen:
  - `rx_algrunflag = 1`, `zcross > 0`, `toploc > 0`, `0.95 < sensitivity_a2 < 1`, `surface_flag = 1`, non-empty `predict_stratum`.
  - `degrade_flag` within a list that **includes some non-zero values**.
- In R:
  - `l2_quality_flag == 1` (falling back to `quality_flag`), applied to all metrics.
  - Sensitivity above 0.98 for PFT 2 in regions 4-6, otherwise above 0.95.
  - Height checks: `elev_highestreturn - DEM < 150` m and `elev_hr - elev_lowestmode < 120` m.
  - Water persistence below 10, and `pft_class > 0`.
  - Urban (≥ 50%) sets AGBD to 0.
  - A bad-granule list.
  - Shots with AGBD outside [0, 2500) are dropped **for every metric**.
- **Season:** for PFT 3-4 only, leaf-off is recomputed per shot from the VNP22Q2 phenology raster (leaf-on = phases 3-4). Other PFTs keep the product flag. Shots in "rh98_only" strata are exempt.
- **Beams:** plain V3 uses all 8. `_power` uses power beams only, but it still uses the all-beam matches.

**Output**
- Shot-level, `G/EXTRACTED_gedi_V3[_power]/...`.
- Shots get their PA/control label by **rasterizing** the matched points onto the template. That keeps one value per cell, so a control reused by several PAs, or a cell in overlapping PAs, keeps only one label (see section 4).

### Step 6: PA vs control change (`M/gedi_change_step6_V3*.ipynb`)

**Aggregation** (plain V3)
1. Shots → cell × period **mean**.
2. Cells → **PA × status × period**, as the unweighted mean of cell means.
   - The unit is the PA: all of its PA cells against all of its matched controls, pooled. Matched pairs are not used.
- Metrics: rh_98, agbd, cover, pai, wsci. There is **no RH50 or FHD**.

**Change**
- `chg = c2025 - c2020` for PA and control, and `pct_chg`.
- **DiD:** `<metric>_dd = chg_PROT - chg_CONT`.

**Totals**
- `agb_tot = agbd × GIS_AREA × 100` (Mg); carbon = 0.49 × that.
- The control "total" uses the PA's own area.

**QC, after the CSV is written**
- `n_shots >= 5` for every status × period, and non-missing AGBD change: 1,107 of 1,147 PAs pass.

**Printed national result**
- Area-weighted ΔΔ AGBD = **+4.08 Mg/ha**: PA -8.5, control -12.6.
- Wilcoxon signed-rank test of `agbd_dd` against 0: p = 2e-15.

**Variants**
- `_power`: rh_98 and agbd only. Its "pair QC" uses `UID` as the pair id, which is wrong (section 4). 1,136 rows, 479 after QC.
- `_elevfilt`: drops matched groups where any cell's ground elevation changed by 5 m or more between epochs. 1,058 rows. It **overwrites** `BRA_GEDI_cell_period_means_V3.RDS` under the same name.

**The final CSV**
- `G/EXTRACTED_gedi_V3/BRA_pa_final_PROT_vs_CONT_change_metrics_V3.csv`: **1,147 rows × 88 columns**, written by plain `step6_V3` (all beams, no elevation filter).
- Its file date is Nov 3, older than the notebook's last save (Nov 25).
- Suffixes: `_PROT` = PA cells, `_CONT` = matched controls. `c2020` = 2019-21, `c2025` = 2024-25.

| Columns | Meaning |
|---|---|
| pa_id, designation, biome_name, biome_id, IUCN_CAT, GIS_AREA | Identifiers. designation is IT or PA; GIS_AREA in km² |
| n_cells_{c2020,c2025}_{CONT,PROT}, n_shots_{…} | Sample sizes |
| {rh_98,agbd,cover,pai,wsci}_{c2020,c2025}_{CONT,PROT} | Mean of cell means |
| {metric}_chg_{CONT,PROT}, {metric}_pct_chg_{…} | Change, percent change |
| {metric}_dd, {metric}_pct_dd | Difference-in-differences |
| n_cells_chg_*, n_shots_chg_* | Change in counts |
| agb_tot_*, c_tot_*, *_PROT_minus_CONT | Area totals (AGB Mg, carbon Mg C) and their differences |

First 5 rows (selected columns; P/C = PROT/CONT; AGBD in Mg/ha):

| pa_id | desig | biome | IUCN | GIS_AREA | cells 20 P/C | shots 20 P/C | shots 25 P/C | agbd20 P/C | agbd25 P/C | agbd_chg P/C | agbd_dd | rh_98_dd | cover_dd | pai_dd |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 42 | PA | Amazon | Ia | 3924.7 | 84/57 | 1993/1537 | 1209/1165 | 62.5/62.1 | 52.2/61.8 | -10.31/-0.28 | -10.03 | -1.28 | -0.067 | -0.269 |
| 43 | PA | Amazon | Ia | 4077.6 | 384/64 | 10805/1362 | 8160/1465 | 215.7/195.0 | 213.5/204.5 | -2.18/9.49 | -11.67 | -1.00 | -0.031 | -0.212 |
| 44 | PA | Amazon | Ia | 3468.6 | 477/118 | 12298/3337 | 10081/2729 | 257.7/232.5 | 239.0/176.8 | -18.70/-55.67 | 36.97 | 2.35 | 0.025 | 0.114 |
| 47 | PA | MataAtlantica | Ia | 278.6 | 43/33 | 1465/928 | 881/849 | 198.5/131.8 | 153.2/104.8 | -45.31/-27.01 | -18.30 | -0.80 | 0.056 | 0.299 |
| 48 | PA | MataAtlantica | Ia | 187.2 | 7/3 | 138/79 | 187/55 | 93.4/71.9 | 88.3/80.9 | -5.16/8.95 | -14.12 | -2.57 | -0.076 | -0.378 |

Cell counts are the same in both periods (the Step 1 rule), so only the c2020 counts are shown.

## 2. TMF disturbance notebook (`G/gedi_change_brazil/BRA_TMF_defor_degrad_mask.ipynb`)

An R notebook. The last section is dated 4/8/2026. I did not review `..._v2.ipynb`.

**TMF layers**
- JRC TMF 2024 release (`/projects/my-public-bucket/JRC_TMF2024/`).
- Only the **`DeforestationYear`** and **`DegradationYear`** layers: 30 m (0.00027°), the South America (`_SAM_`) 10° tiles, 18 of which overlap Brazil (cells 3, 9).
- It does not use the Annual Change collection or the Transition map.

**Deforestation**
- A pixel is deforestation if its `DeforestationYear` is 2022, 2023 or 2024 (cells 4-5).
- Output: `BRA_TMF_deforestation_mask_2022_2024.tif`.

**Degradation**
- A pixel is degradation if its `DegradationYear` is 2022, 2023 or 2024 (cell 6).
- Output: `BRA_TMF_degradation_mask_2022_2024.tif`.

**Combining them**
- Deforestation overrides degradation: 2 = deforestation, 1 = degradation, 0 = other (cell 18).
- Fire comes from MapBiomas Fire collection 4, annual burned area 2022-2024, as "burned in any year" (cells 21-30).
- Final classes (cell 42), in `disturbance/BRA_disturbance_attribution_2022_2024.tif`:
  - 0 = no disturbance
  - 1 = fire only
  - 2 = degradation, no fire
  - 3 = degradation + fire
  - 4 = deforestation, with or without fire

**Fractions on the 1 km grid** (cells 63-69)
- For each class k, `project(dist_class == k, BRA_template_rast, method = "average")`.
- So each fraction is **the share of all 30 m pixels in the 1 km cell that are in class k**.
- The denominator is every valid pixel inside Brazil, **not baseline forest**. Class 0 mixes undisturbed forest with non-forest.
- Outputs:
  - `disturbance/BRA_{undisturbed,fire_only,degrad_other,degrad_fire,defor}_frac_1km.tif`
  - the 5-band stack `BRA_disturbance_fractions_on_1km_GRID.tif`, also copied to `INPUT_covariates_2020/`

**PA summaries** (cells 53-58)
- Area by class per WDPAID (`cellSize` in ha), over 3,403 polygons from `WDPA_Nov2025_BRA_poly_v2.RDS`.
- Fractions are of total PA area.
- Output: `disturbance/BRA_disturbance_pa_summaries.csv`.

**Thresholds**
- **None in this notebook.** The fractions are continuous from 0 to 1, and there are no cut-offs or binary cell masks. Any threshold must be applied downstream; I didn't check the attribution and stats notebooks.

**Points to note**
- **The window doesn't match our epochs.** Disturbance covers 2022-2024. It misses disturbance during T1 (2019-2021), and in 2025, which is inside her T2 and inside ours.
- **Repeat degradation may be missed.** As I understand TMF, `DegradationYear` is the *first* year of degradation. If so, a pixel first degraded before 2022 and degraded again in 2022-2024 is not flagged. I haven't checked this against the TMF documentation.
- The fractions are computed in lon/lat and averaged onto EASE2, so pixel-area differences are ignored. That's negligible at this scale.

## 3. Her Step 5/6 vs our `extract_tile.py` (rules in CLAUDE.md)

| Item | Ours | Hers (V3, the CSV) |
|---|---|---|
| Data | gedi-tiled-v3 (global, L2A/L2B/L4A/L4C joined) | `ameliah/gedi-test/brazil_tiles` (Brazil test set) |
| Beams | Power only (5, 6, 8, 11) | All 8 (`_power` variant: power shots but all-beam matches) |
| L2A quality | `l2a_quality_flag_rel3 = 1` | `l2_quality_flag == 1` (or `quality_flag`) |
| PAI / FHD / cover quality | Also `l2b_quality_flag_rel3 = 1` | None separate |
| AGBD quality | `l4a_quality_flag_rel3 = 1` | None separate; AGBD in [0, 2500) required for **all** metrics |
| degrade_flag | = 0 | Some non-zero values allowed |
| Leaf-off / season | Product `leaf_off_flag = 0` **and** a fixed DOY window per PFT × region, same in both epochs | Per-shot leaf-off from phenology for PFT 3-4 only; product flag otherwise; no DOY windows |
| Extra screens | None | Sensitivity > 0.95 / 0.98, height checks vs DEM (150 m, 120 m), water, urban → AGBD 0, bad granules |
| PFT | Classes 1-6 only | `pft_class > 0` |
| T2 years | 2024-2026 | 2024-2025 |
| Metrics | RH98, RH50, PAI, FHD, cover, AGBD | rh_98, agbd, cover, pai, wsci |
| Cell id | `row * 34704 + col` (global) | terra cell in a Brazil crop (conversion in section 1) |
| Cell statistics | n, mean, sd per metric × epoch; tile-edge cells pooled | Mean (and n) per cell × period; shots read only from the point's own 1° tile ±0.02° |
| Cell eligibility | Any cell with shots | ≥ 2 shots and ≥ 2 tracks in **both** periods (a same-cell panel) |
| Comparison unit | Aggregate only (PA, matched set, country, biome); minimum n per epoch | PA (pooled PA vs pooled controls); QC n_shots ≥ 5, applied after writing |

The other half of item 2 needs our `analysis/did_country.py`, which doesn't exist yet. I can only compare her Step 6 design with the CLAUDE.md rules:
- "Never cell-to-cell": her design is consistent with it, since she aggregates per PA.
- "Minimum sample size in each epoch, and report how many units qualify": she reports the count only after the CSV is written.
- "Decide filters, minimum sample sizes and headline comparison before looking": not determinable from the code.

## 4. Likely problems in her pipeline

These are from code reading. My confidence is given for each.

1. **The propensity model is Gaussian, not logistic** (`matching_func_2025.R:49,70,114,168`, no `family=`). Verified. Matching on a linear-probability score still balances the covariates it uses, but it probably wasn't intended.
2. **Exact matching on biome and ecoregion never happens**, because the columns are numeric. Verified in the RDS: `wwfbiom` and `wwfecoreg` are float64. Both are 0 in the rows sampled (BAe, Amazon). That makes me suspect the `biome`/`ecoreg` rasters, but I haven't checked them.
3. **Labels collide when matched points are rasterized** (Step 5). A control reused by several PAs, or a cell in overlapping PAs, keeps only one pa_id and status, so some PAs lose controls. Supporting evidence:
   - The matching is 1:1, but the CSV shows PA/control cell counts like 84/57, 384/64 and 477/118.
   - In the RDS, 740,868 rows map to only 538,587 distinct cells.
   - Confidence: about 80%.
4. **`_power` pair QC uses `UID` as the pair id.** `UID` is a row index that restarts per PA, so the "10,288 of 132,202 pairs pass" and "479 PAs" results don't mean anything. Confidence: high.
5. **The `_elevfilt` pair key isn't unique.** `matched` labels restart with every matching call, so they repeat within a PA. Dropping one bad group can drop unrelated cells, which may explain why 77% of rows were removed. Confidence: medium.
6. **`_power` is a hybrid:** power-beam shots with all-beam matches and eligibility. The power-beam matching points (Nov 25-26) haven't been matched yet.
7. Smaller points:
   - The `match_wocat` loop over land-cover values repeats the same matching each time (wasteful, not wrong).
   - `agb_tot` for controls uses the PA's area.
   - `_elevfilt` overwrites the V3 cell-mean RDS.
   - Some notebook outputs are stale.

## 5. Our side: what's ready for the Brazil run

**Extraction**
- All 12,080 tiles are done.
- `output/gedi_pa_change_cells_combined.parquet`: 76,377,681 rows, one per cell × epoch, with n/mean/sd for all 6 metrics. Tile-edge cells are pooled, and the pooling was checked against raw shots.

**Joining her matches to our cells**
- Verified: her `lat`/`lon` sit exactly at EASE2 cell centres (offset 0.5/0.5 in every row), so each converts to exactly one of our `cell_id`s.
- Her 538,587 matched cells have our GEDI data in **424,993 (T1)** and **373,511 (T2)** cells.
- That's fewer than all of them because our rules are stricter (power beams, DOY windows, l2b/l4a flags) and our data source is different.

**Reading her RDS files**
- R isn't installed on this Hub image. For this review I installed `pyreadr` into a scratchpad folder, not into the notebook environment.

**What's missing**
- `analysis/did_country.py`, `analysis/rtm_check.py` and `matching/README.md`. I didn't write them.
- CLAUDE.md says filters, minimum sample sizes and the headline comparison must be decided **before** looking at PA-versus-control results. Writing the DiD script myself would mean choosing those on the spot. These decisions are needed first:
  1. **Unit:** per PA, as she does, or per matched group.
  2. **Minimum n per epoch** at that unit (cells? shots?), and how many units qualify.
  3. **Which matches:** her all-beam V3 matches, or wait for power-beam matches.
  4. **Reused controls:** how to handle controls reused across PAs, and cells in overlapping PAs (keep every row of the matching table, rather than rasterizing as Step 5 does).
  5. **Weighting:** per-cell means weighted by shots or not; area weighting for national summaries.
  6. **"RTM check":** what exactly it should test. I'm assuming a regression-to-the-mean check comparing PA and control T1 levels against their change.

## 6. Comparing with her final CSV

We have **no Brazil DiD numbers yet** (section 5), so a like-for-like comparison isn't possible. Once we do, expect differences for these reasons, even with the same matches:

- **Beams:** we use power beams only (her CSV uses all 8).
- **Season:** our DOY windows and her PFT-3-4 phenology leaf-off select different shots.
- **T2:** ours runs to 2026, hers to 2025. Today both end at July 2025.
- **Quality flags:** separate L2B/L4A flags for us. Her shared AGBD-range filter drops shots for every metric.
- **Cell eligibility:** she requires data in both periods in the same cell; we don't.
- **Labels:** her rasterized labels likely lose shared controls (section 4, item 3). Joining on the full matching table fixes that.

Her headline for reference: area-weighted ΔΔ AGBD = **+4.08 Mg/ha** (PA -8.5, control -12.6 Mg/ha), with 1,107 of 1,147 PAs passing `n_shots >= 5`.
