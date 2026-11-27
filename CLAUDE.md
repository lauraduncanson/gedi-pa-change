# gedi-pa-change: project instructions

## Goal

Estimate 5-year change in forest structure and biomass from GEDI inside protected areas (PAs)
and in their matched counterfactuals, across the GEDI domain (51.6 S to 51.6 N).

- Metrics: RH98, RH50, PAI, FHD, cover, AGBD.
- Epochs by calendar year: T1 = 2019-2021, T2 = 2024-2026. T2 data currently ends July 2025;
  more will arrive, so everything must be re-runnable for T2 alone.
- Headline result: difference-in-differences (PA change minus matched-control change).
  Absolute T2 - T1 change is secondary, because a time-dependent instrument bias would
  affect it directly.
- Deadline: first global results for AGU, mid-December 2026. Prefer simple, working, checked
  code over new methods.

## Environment

- NASA MAAP **Hub** (not ADE). Home is `/home/jovyan`. Repo: `/home/jovyan/GEDI/gedi-pa-change`.
- Never write code that prompts for an Earthdata login. Read data with the Hub/DPS AWS role
  (boto3 credentials passed to DuckDB; see `extract_tile.py: connect()`).
- Bucket mounts (`/home/jovyan/my-private-bucket`, `/home/jovyan/shared-buckets`) do not
  support file seeks. Write `.npz`, parquet and similar files to a local folder first, then copy.
- `curl` is not installed on the Hub image; use `wget` or Python.
- MAAP login currently has a 403 problem. Workaround: sign in at urs.earthdata.nasa.gov first,
  then return to MAAP. If job submission fails with an auth error, suspect this before the code.

## Data

- GEDI: `gedi-tiled-v3`, GeoParquet partitioned by 1-degree tile and year, shots joined across
  L2A/L2B/L4A/L4C. Root: `s3://nasa-maap-data-store/file-staging/nasa-map/gedi-tiled-v3`.
  Use `manifest.txt` to find files; do not glob the bucket.
- GEDI ancillary rasters (1 km EASE2, EPSG:6933, 34704 x 14616) and WDPA files:
  `/home/jovyan/shared-buckets/leitoldv/GEDI_global_change/WDPA_comparison/`
  - `GEDI_ANCI_PFT_r1000m_EASE2.0_UMD_v6.tif`, `GEDI_ANCI_CONTINENT_r1000m_EASE2.0_UMD_v5.tif`
  - `GLOBAL_PAs_all_nonMARINE_2025.gpkg` (WDPA 2025; the version used to define "protected")
- Phenology: `/home/jovyan/GEDI/GEDI_ANCI_PHENO_vnp22q2_MODISsinusoidal_r0500m.tif`
  (VIIRS VNP22Q2, year 2022; bands 1-6 = greenup onset, mid-greenup, maturity, senescence onset,
  mid-greendown, dormancy).

## Fixed analysis rules (do not change without asking Laura)

- **Power beams only** (`beam_l2a` in 5, 6, 8, 11). Coverage-beam quality degrades over the
  mission, so coverage beams are not used for change.
- Shot filter: `l2a_quality_flag_rel3_l2a = 1`, `degrade_flag_l2a = 0`,
  `land_cover_data_leaf_off_flag_l2a = 0`. PAI/FHD/cover also need `l2b_quality_flag_rel3_l2b = 1`;
  AGBD needs `l4a_quality_flag_rel3_l4a = 1`. `rh_l2a` holds RH0..RH100 (RH98 = element 99, 1-based).
- Growing-season windows: `extract_tile/doy_windows.csv`, one DOY window per PFT x region
  (Africa split N/S at the equator). Same window for both epochs. Derived from the phenology
  grid: evergreen PFTs = full year; otherwise the 90-day (else 60-day) window with at least 60%
  of pixels leaf-on (mid-greenup to mid-greendown); otherwise full year. Shots with PFT outside
  classes 1-6 are dropped. `start_doy > end_doy` means the window crosses the new year.
- Cell key: `cell_id = row * 34704 + col` on the GEDI 1 km EASE2 grid. In R/terra the same cell
  is `cell_id + 1`.
- T1 and T2 are different footprints and rarely share a 1 km cell. Compare at aggregate level
  (PA, matched set, country, biome), never cell-to-cell. Require a minimum sample size in each
  epoch at the reporting unit and report how many units qualify.
- Decide filters, minimum sample sizes and the headline comparison before looking at
  PA-versus-control results.
- Change results are reported uncorrected, corrected (central) and at the low/high drift bounds
  (see "Ground-drift correction").

## Ground-drift correction (2026-11-27)

Do not change the drift table, the L4A model files or the overrides without asking Laura.
- **Mechanism:** under dense canopy GEDI picks the ground (lowest mode) progressively too high
  over the mission (missed ground). Canopy top is unaffected, so RH and L4A AGBD decline
  artificially and T2 looks like loss. Measured as the change in `elev_lowestmode_l2a -
  digital_elevation_model_l2a` (TanDEM-X) within 1-degree tiles, slopes < 10 deg, by Landsat tree
  cover, in 650 random tiles. Open ground (0-10% cover) drifts ~0.03 m/yr (whole-waveform shift;
  cancels in RH, so subtracted). Excess drift is ~0 below 50% cover and rises to 0.087 m/yr
  (95% CI 0.056-0.118) at 90-100%; a gradual ramp since 2020, not a step at the 2023-24 storage
  period. 2019 sits below the trend, so the rate is fitted on 2020+ and 2019 has its own offset.
  Check: canopy-top elevation shows no cover-dependent drift; directly measured RH98 bias at
  90-100% cover (0.13 m/yr, CI 0.04-0.26) agrees.
- **Table:** `extract_tile/drift_correction.csv` (final 650-tile fit; columns excess_rate /
  excess_lo / excess_hi per cover class, offset_2019, ref_fit, t1_year_weights).
- **Applied** per shot in both epochs: `correction_m = drift(cover, year) - mean drift over T1
  shots`, added to every RH percentile of `rh_l2a`. T1 means are unchanged; T2 is put on T1's
  ground finding. Code: `extract_tile/correction.py` (self-contained), called from
  `extract_tile.py: apply_correction()`. Analysis scripts: `analysis/ground_drift_quicklook.py`,
  `ground_drift_within_tile.py`, `drift_correction.py`, `top_vs_ground_check.py`,
  `find_l4a_models.py`, `l4a_model.py`.
- **L4A model:** AGBD recomputed from corrected RH with the V3 models from the product
  (`ANCILLARY/model_data` of a V003 L4A granule -> `extract_tile/l4a_model_data.json`):
  `x_i = sqrt(rh_l2a[RH_i] + 100)`, `y = max(b0 + sum b_i x_i, 0)`, `AGBD = bias_correction * y^2`;
  no capping; L4A uses `rh_l2a` for every shot.
- **South America overrides:** for `EBT_SA` and `DBT_SA` the product applies a model not in
  `model_data` (the stored rows hold the South Asia coefficients). The applied model, recovered
  from the product: sqrt(AGBD) = -143.7771 + 7.471 sqrt(RH50+100) + 6.7511 sqrt(RH98+100), bc 1.0,
  in `extract_tile/l4a_overrides.json`. With it, recomputed AGBD matched `agbd_l4a` to < 0.1% for
  10,689 shots in 15 strata.
- **Outputs:** shot files add `rh_shift_m`, `rh98_c`, `rh50_c`, `agbd_c` (central), `agbd_c_lo`,
  `agbd_c_hi` (95% drift-rate bounds); cell files add `<metric>_n/_mean/_sd` for these. Summary
  JSON adds `l4a_checked`, `l4a_mismatch_frac` (should be ~0) and `l4a_mismatch_strata`; a
  non-empty `l4a_mismatch_strata` means a stratum whose applied model differs from the stored
  table: report it, do not paper over it.
- **Not corrected:** PAI, FHD and L2B cover (missed ground biases them too). Lead with RH98 and
  AGBD for change; treat L2B metrics as uncorrected.
- Hub test, tile S03_W060 (Amazon): 85,091 shots recomputed, mismatch 0.0.

## Repo layout

- `extract_tile/`: DPS algorithm. `extract_tile.py --tile <id>` writes `shots_<tile>.parquet`,
  `cells_<tile>.parquet` (1 km cell x epoch: n, mean, sd per metric) and always
  `summary_<tile>.json`. Also the Containerfile, `algorithm_config.yml`, `doy_windows.csv`.
- `dps/`: `submit_tiles.py`, `monitor.py`, `collect.py`, `common.py`. Run from the repo root.
- `notebooks/plot_tile_histograms.py`: T1 vs T2 histograms per tile, inside/outside PAs.
- `.github/workflows/extract-tile.yml`: builds the image and commits the CWL on each push to
  `extract_tile/`. Because the Action commits to the repo, always `git pull` before editing.

## DPS conventions

- Algorithm: `gedi-pa-change-extract-tile`, version `main`, registered with
  `maap.deploy_algorithm(...)` from the committed CWL (deployment job 499, 2026-11-23).
  Confirm it deployed and get the process id with `maap.list_algorithms()`.
- One job per tile, about 8 GB of memory. A tile is done when its `summary_<tile>.json` exists
  under the DPS output folder. Outputs land in
  `/home/jovyan/my-private-bucket/dps_output/gedi-pa-change-extract-tile_<N>/main/<yyyy>/<mm>/<dd>/...`,
  where `<N>` is the deployment number (14 as of 2026-11-24). On S3 that is
  `s3://maap-ops-workspace/lduncanson/dps_output/gedi-pa-change-extract-tile_<N>/...`.
  The default `--done-dir` (`common.DPS_OUTPUT`) is that S3 prefix, read with boto3/DuckDB:
  the bucket mount can drop ("Transport endpoint is not connected", not fixed by a Hub
  restart on 2026-11-24) and globbing it is slow. If the mount looks empty, `--resubmit`
  without `--tiles` would resubmit everything. Never resubmit finished tiles.
- Full extraction run finished 2026-11-24: all 12,080 tiles done (6 needed one retry after
  DPS credential/502 errors). Merged table: `dps/collect.py` -> `output/gedi_pa_change_cells.parquet`
  (one row per tile x cell x epoch). Then `dps/combine_cells.py` pools cells on tile edges ->
  `output/gedi_pa_change_cells_combined.parquet` (one row per cell x epoch, 76.4 M rows; `n_tiles`
  > 1 for ~1.4% of rows). Use the combined table for matching and all aggregates.
- Ramp up: 2 test tiles (`S13_E028`, `N11_W005`), then 100, then all. Check outputs at each step.
  Job ids are logged in `dps/jobs.csv` (not committed).
- **drift-v1** (ground-drift correction): deployment job 517, **process 252**, 2026-11-27. Outputs in
  `.../gedi-pa-change-extract-tile_14/drift-v1/`; always pass
  `--done-dir s3://maap-ops-workspace/lduncanson/dps_output/gedi-pa-change-extract-tile_14/drift-v1/`
  (the default prefix also matches the `main` run and would count every tile as done). The `main`
  job log is `dps/jobs_main.csv`; `dps/jobs.csv` holds drift-v1 jobs. Test tiles S03_W060 and
  S13_E028 passed (L4A mismatch 0.0, no mismatch strata, `agbd_c_mean` in cells).
- Confirmed by the test tiles: DPS workers can read gedi-tiled-v3; queue `maap-dps-worker-16gb`
  works (process id 237; ~15 s per tile, ~200 tiles/min); output location as above.

## Matching (next on the critical path)

- Counterfactual matching code: https://github.com/GEDI-PA/vl_GEDI-PA_2024 (R, runs per country
  on DPS; steps 1-3 grid and covariates, step 4 propensity matching, step 5 GEDI extraction).
  Matched results already exist for Brazil; other countries have not been run.
- Plan: keep steps 2-4 as they are. Patch step 1 to take GEDI coverage from the merged cell
  table (cells with data in either epoch) instead of GEDI-subsetter GeoPackages. Replace step 5
  with a join of matched cell ids to the cell table.
- First end-to-end test: Brazil, using the existing matches joined to T1/T2 cell means.
- Run large countries first (Russia, USA, Canada, China, Australia). Check covariate balance on
  T1 GEDI structure, not only on the matching covariates.

## Brazil matching on DPS (`gedi-pa-change-match-brazil`)

- Veronika's Brazil V3 step 4, copied verbatim into `matching/r/brazil_v3/` (provenance in
  `SOURCE.md`) and patched there; packaged in `matching/dps/` (conda-forge R on micromamba; an
  Ubuntu/rocker image failed on DPS with an OpenSSL `[CRYPTO] unknown error`). The GitHub Action
  `match-brazil.yml` builds the image (`ghcr.io/lauraduncanson/gedi-pa-change-match-brazil:main`,
  public) and commits the CWL; register from the Hub with `maap.deploy_algorithm(...)`.
- Use **process 251 (version v6)** (v4 added the L4B-class arm and `input_prefix`; v5 dropped unused
  biome/ecoregion factor levels; v6 made the no-ecoregion fallback return the same columns as the
  main branch, which made 39 of 554 global v5 jobs fail and could drop chunks silently). Global
  v5 outputs are superseded; rerun everything on v6. Inputs: `region` (BAe BAw BCa BCe BMa BPp BPt),
  `covariates` (`with_l4b` = her V3 model, `no_l4b` = `gedi_l4b` dropped, nothing else changed),
  `n_cores` (= the worker's vCPUs, e.g. 32 on `maap-dps-worker-32vcpu-64gb`), `part` (`i/n`).
  Submit with `dps/submit_matching.py` (`--parts N` splits a region; ids in `dps/match_jobs.csv`).
- Verified 2026-11-26: `with_l4b` reproduces her `MATCHING_results_V3` exactly for BPp and BPt
  (same cells, propensity scores and match labels), including v3's split runs and single
  land-cover loop. Outputs land in `dps_output/gedi-pa-change-match-brazil_14/<version>/...`
  (`MATCHING_results_V3/<region>/*_matching_results.RDS` and `<region>_allPA_matched_results.RDS`).
- Kept as in her code on purpose (sensitivity run = "all else unchanged"): Gaussian `glm`
  propensity model, biome/ecoregion strata that never apply (numeric columns).
- Data versions: her matching inputs and `brazil_tiles` shots are most likely GEDI V002 (same
  shots sit ~2 m from gedi-tiled-v3); `gedi_l4b` is probably L4B V2.1. Our outcomes are V003.

## Matching design decision (Laura, 2026-11-26)

- **Primary design: GEDI L4B as a continuous covariate in the propensity model** (Veronika's V3
  model; `covariates = with_l4b`). Basis: in Brazil and Cameroon, dropping L4B doubled-to-sextupled
  the T1 PA-control AGBD imbalance and lowered the DiD; 50 Mg/ha L4B classes reproduced the
  continuous DiD (Cameroon +3.1 vs +3.0 Mg/ha), so regression to the mean from fine L4B matching
  is not driving the effect, while baseline differences are. See results/BRA_first_look.md
  (section 9) and results/CMR_l4b_check.md.
- `no_l4b` and `l4b_classes` are sensitivity checks only. The RTM check (`rtm_check.py`) is an
  upper bound, not a design criterion (its gap is negative even without L4B).
- Report T1 PA-minus-control balance (AGBD, RH98) with every DiD. Cells without L4B drop out of
  matching; that is acceptable (they mostly lack GEDI data too).

## Global matching design (Laura, 2026-11-26)

Global steps 1-3 built only from global layers, then her V3 step 4 (DPS v4) with continuous L4B:
- **Covariates: the 2020 set**, as in her Brazil V3, from
  `leitoldv/GEDI_global_change/INPUT_covariates_2020/`: d2roads, dcities, dem, slope,
  pop_cnt_2020, pop_den_2020, tt2cities_2015, WorldClim 2010-2018 (prec, tavg, tmax, tmin),
  gedi_l4b; land cover = `lc2020` (7 classes: forest, grassland, agriculture, wetlands,
  artificial, bare/other, water) as the exact land-cover stratum (replaces Brazil-only MapBiomas).
- **Exact matching on WWF ecoregion** (`wwf_ecoreg`; her code falls back to dropping ecoregion when
  matching fails), plus land cover. Biome and ecoregion are written as R factors with the same
  levels in PA and control tables so her strata apply (in her Brazil V3 they were numeric and
  never applied; her Brazil inputs are left as they are).
- **Matching unit: country x WWF biome** (like her Brazil biome regions), each with its own
  control pool; PAs spanning biomes are split by biome.
- **Grid:** 1 km EASE2 cells with GEDI data in either epoch (combined cell table), centre inside
  the country. **Treatment:** cell fully inside a PA (centre and four corners inside; her rule
  was >= 99.9% inside), PAs as in the global PA change analysis (designated / inscribed /
  established, no UNESCO-MAB). **Controls:** cell centre more than 10 km from every WDPA polygon
  (any status), as in her Brazil V3.

## Control buffer test: 2 km vs 10 km (Laura, 2026-11-26; decided before running)

Her rule excludes controls within 10 km of any WDPA polygon. In heavily protected countries this
leaves few controls (Costa Rica: 1,323 control cells at 10 km vs 6,809 at 2 km, for 6,445 PA
cells). Test whether a 2 km buffer changes the estimate:
- **Countries:** Costa Rica (CRI), Guatemala (GTM), Ghana (GHA), Cameroon (CMR); all their
  country x biome units. Inputs from `matching/prep_country_inputs.py` with `--buffer-km 10` and
  `--buffer-km 2` (everything else identical); matching with continuous L4B (DPS v5).
- **Same DiD rules** as elsewhere (per PA, >= 10 cells per group, PA-cell weights).
- **Compared per country:** AGBD and RH98 DiD (weighted and median), PAs qualifying, T1
  PA-minus-control AGBD and RH98, and the paired difference DiD(2 km) - DiD(10 km) over PAs
  qualifying in both, with a 95% bootstrap interval over PAs.
- **Decision rule:** if the paired 95% interval includes 0 for both AGBD and RH98 in every
  country with >= 10 PAs qualifying in both, edge effects are judged small and the global run
  uses **10 km by default, 2 km for units where the 10 km pool has fewer complete control cells
  than PA cells**. Otherwise keep 10 km everywhere (thin units go to the imputation rule).
- **Outcome (2026-11-26): keep 10 km everywhere.** Ghana: no detectable difference. Cameroon:
  2 km lowered the DiD by 3.9 Mg/ha AGBD [95% CI −8.2, −0.3] and 0.28 m RH98 (near-PA controls
  lost less biomass, i.e. spillover), so the rule was not met. Costa Rica had no comparable controls
  at either buffer. See results/buffer_test.md.

## PAs without a valid counterfactual (Laura, 2026-11-26; to implement in analysis/totals.py)

No cross-border matching (national trends differ too much to assume comparable controls).
Instead, impute the PA's DiD from comparable measured PAs:
- **No valid counterfactual** = the PA has no matches, or fails the 10-cell rule for the metric,
  or its T1 AGBD balance is poor: |standardised mean difference| > 0.25 between its PA and
  matched control cells (SMD = (mean PA - mean control) / pooled SD of T1 cell means).
- **Strata**, in order, each needing >= 5 measured PAs: WWF/Olson biome x country, then GEDI
  region (continent) x biome, then biome. A PA's biome is the majority biome of its cells.
- **Imputed DiD** = PA-cell-weighted mean DiD of the measured PAs in the stratum; its uncertainty
  is their spread. The PA's own T1/T2 change inside the PA is still measured, not imputed.
- **Always report** measured-only and measured + imputed totals side by side, with the share of
  PAs and of PA area that is imputed, by country and globally.

## Totals (Laura, 2026-11-27; decided before computing)

Implemented in `analysis/totals.py`.
- **PA universe:** kept PAs (as in the global PA change analysis) with at least one GEDI cell
  (centre inside, `output/cell_pa.parquet`). Measured = valid counterfactual as above; everything
  else is imputed. PAs whose biome stratum has < 5 measured PAs at every level stay unimputed and
  are reported as excluded.
- **Area:** WDPA polygon area in EASE2 (equal area), parcels of one WDPAID dissolved, clipped to
  51.6 S - 51.6 N. Overlaps between PAs are counted once: within each country, PA areas are scaled
  so they sum to the dissolved union of that country's universe PAs.
- **Country** = WDPA ISO3 (first code for multi-country sites); **biome** = majority WWF biome of
  the PA's cells; **region** = majority GEDI continent region of its cells.
- **Weights:** area everywhere (means, and the stratum means used for imputation; this replaces
  "PA-cell-weighted" above).
- **Outputs:** area-weighted mean DiD per metric (kept for figures), and for AGBD the total
  DiD x area in Tg; by country, biome, region and globally; measured-only and measured + imputed;
  95% bootstrap over measured PAs (imputed values recomputed in each resample).

## L4B matching check in Cameroon (decided 2026-11-26, before running)

Brazil's AGBD DiD fell from +4.6 to +2.0 Mg/ha when L4B was dropped from the matching, but
control T1 balance got worse (results/BRA_first_look.md, section 9). Test in a smaller country:
- **Country: Cameroon (CMR)**, using the 2025 global-pipeline step 1-3 outputs
  (`shared/abarenblitt/GEDI_global_PA_v2/Matching_Results/CMR/`: `CMR_testPAs/prepped_pa_*_wk24.RDS`,
  `WDPA_grids/CMR_prepped_control_wk24.RDS`). Liberia (LBR) is the practice job.
  Ghana was chosen first but switched before any matching, after a covariate-only overlap check:
  0% of Ghana's control cells (with L4B) reach 100 Mg/ha against 24% of its PA cells, so the arms
  would differ mainly in which forest reserves cannot be matched. Cameroon: 50% of PA cells and
  42% of controls >= 100 Mg/ha, 47 PAs, 39k PA cells and 185k controls with L4B.
- **Matching code: Veronika's Brazil V3 step 4 (verified)**, same for all arms. `gedi_l4b` is sampled
  from the same global `gedi_l4b.tif`; biome/ecoregion (absent) are set to a constant, so their
  pre-filter is inactive; `mangrove` is not used (not in her V3 model).
- **Three arms:** `with_l4b` (L4B continuous in the propensity model, her V3), `no_l4b`, and
  `l4b_classes` (L4B removed from the propensity model; exact matching within 50 Mg/ha classes,
  class = floor(L4B / 50), classes >= 400 Mg/ha merged).
- **Same DiD rules as Brazil** (per PA, >= 10 cells in each PA/control x T1/T2 group, PA-cell weights).
- **Reported for each arm, decided now:** AGBD and RH98 DiD (weighted and median), PAs qualifying,
  T1 PA-minus-control AGBD and RH98 (balance), and the RTM gap (all matched cells; no
  disturbance layers exist for Cameroon). The classes arm is preferred if it gives T1 balance at
  least as good as `with_l4b` and an RTM gap closer to 0; whichever DiD that arm gives is the one
  to carry forward. No further arms after seeing results.

## Global PA change without matching (decided 2026-11-26, before running)

Secondary result (absolute change; exposed to instrument/epoch bias, see Goal), separate from the
DiD. Implemented in `analysis/pa_change_global.py`.
- **PAs:** WDPA 2025 non-marine (`GLOBAL_PAs_all_nonMARINE_2025.gpkg`), STATUS Designated,
  Inscribed or Established, any STATUS_YR; OECMs (PA_DEF = 0) excluded; UNESCO-MAB Biosphere
  Reserves excluded (Protected Planet convention: their buffer/transition zones are not
  protected). Unit: the site (WDPAID); parcels with the same WDPAID are pooled.
- **Cell rule:** a 1 km cell belongs to a PA if its centre is inside the polygon. Cells in
  overlapping PAs count for each PA; summaries over "all PA land" use each cell once.
- **Change:** per PA and metric, mean of cell means in T2 minus T1 (combined cell table).
  A PA qualifies for a metric if it has >= 10 cells with data in each epoch; report how many
  qualify out of how many have any data.
- **Summaries:** global, by country (ISO3), by IUCN category, by designation type, by governance
  type, weighted by PA cells (area-weighted totals later); plus the union of PA cells by GEDI
  region x PFT. **Context reference:** the same region x PFT summary for cells outside every PA
  (all WDPA statuses and OECMs counted as "inside" for this purpose, so "outside" is unprotected
  land). It is not a counterfactual; it shows change shared with the landscape, including any
  instrument or epoch effect.
- Caveats to state with every result: T2 ends July 2025; Jan-Apr 2024 has no GEDI data (ISS
  storage), which affects windows early in the year; absolute change includes any time-dependent
  instrument bias.

## DiD analysis decisions (decided 2026-11-26, before looking at PA-vs-control results)

Implemented in `analysis/did_country.py` and `analysis/rtm_check.py`; commands in
`matching/README.md`. Always pass the combined cell table
(`output/gedi_pa_change_cells_combined.parquet`).

1. **Unit: per PA.** PA cells vs that PA's matched controls. A control reused by several PAs
   counts for each of them (keep every row of the matching table; never rasterize labels).
2. **Minimum sample: 10 cells** with data in each of the four groups (PA/control x T1/T2),
   per PA and per metric. Report how many PAs qualify out of how many have matches.
3. **Matches (Brazil, for now):** Veronika's all-beam V3 matches,
   `leitoldv/GEDI_global_change/MATCHING_results_V3/BRA_allPA_matched_results_V3.RDS`. Beams only
   matter for outcomes, which come from our power-beam cell table.
4. **Weights:** equal weight per cell within a PA (mean of cell means). Country summaries weight
   PAs by their number of PA cells for now; area-weighted totals come later (`analysis/totals.py`).
5. **RTM check:** `analysis/rtm_check.py` tests the L4B baseline gap as in its docstring, using
   `INPUT_covariates_2020/gedi_l4b.tif` and the four `INPUT_covariates_2020/BRA_*_frac_1km.tif`
   layers (deforestation, degradation-other, degradation-fire, fire-only, 2022-2024); a cell is
   undisturbed when all four are 0.

## Working style

- Explain what you are about to run before running anything that submits jobs or deletes files,
  and report counts (submitted, done, failed) after.
- Test on one or two tiles before any batch. Say plainly what was tested and what was not.
- Keep changes small and commit with clear messages.
