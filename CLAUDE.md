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
