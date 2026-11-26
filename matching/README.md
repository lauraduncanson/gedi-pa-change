# Matching: integration with the GEDI cell table

The counterfactual matching code is Veronika Leitold's R pipeline,
https://github.com/GEDI-PA/vl_GEDI-PA_2024 (run per country):

| Step | Script | Reads GEDI? | Change needed |
|---|---|---|---|
| 1. 1 km grid of cells with GEDI data | `global_process_part1_2024_MAAP_step123.R` | yes (subsetter L4A GeoPackages) | **replace** with `make_country_grid.py` output |
| 2-3. sample covariates for non-PA and PA cells | same script | no | none |
| 4. propensity filter + matching per PA | `global_process_part2_2024_MAAP_step4.r` | no | re-point input paths (below) |
| 5. extract GEDI over matched cells | `global_process_part3_2024_MAAP_step5.r` | yes | **not used**: replaced by `analysis/did_country.py` |

## Step 1 replacement

`make_country_grid.py` lists, per country, the 1 km EASE2 cells inside the country boundary with
at least one quality-filtered, power-beam, growing-season shot in T1 or T2. It writes
`<ISO3>_grid_cells.csv` with columns `x.overlap` (lon) and `y.overlap` (lat) of the cell
centre, which is what the R code builds in step 1, plus `cell_id`.

```bash
python dps/collect.py --out /home/jovyan/GEDI/gedi_pa_change_cells.parquet
python matching/make_country_grid.py --cells /home/jovyan/GEDI/gedi_pa_change_cells.parquet \
    --out /home/jovyan/my-public-bucket/gedi_pa_change/country_grids --iso3 BRA
```

In `global_process_part1_2024_MAAP_step123.R`, replace everything from the comment
`#1.3) extract coordinates of raster cells with valid GEDI data in them` down to and including
the closing `}` of the `for(i in 1:length(iso3_tiles))` loop with:

```r
  # 1.3) GEDI coverage from the gedi-pa-change cell table (matching/make_country_grid.py):
  #      cells with >= 1 quality-filtered, power-beam, growing-season shot in T1 or T2.
  cellgrid.path <- "s3://maap-ops-workspace/shared/<MAAP_USERNAME>/gedi_pa_change/country_grids/"
  GRID.coords <- read.csv(s3_get(paste(cellgrid.path, iso3, "_grid_cells.csv", sep="")))
  GRID.coords <- unique(GRID.coords[, c("x.overlap", "y.overlap")])
```

The next line of the original (`GRID.for.matching <- vect(GRID.coords, ...)`) is unchanged, and
nothing after step 1 uses the variables the loop created. Keep `GRID.for.matching` without
extra attributes: steps 2-3 name covariate columns by position.

## Other edits required before running new countries

1. **Hand-off paths in step 4.** `global_process_part2_2024_MAAP_step4.r` reads the step 1-3
   outputs from hard-coded locations in the `abarenblitt` bucket (`f.path2`, and the
   `Matching_Results/<ISO3>/<ISO3>_testPAs/` prefix), not from the current run's output folder.
   Unless these are re-pointed, step 4 silently matches on old step 1-3 files. Simplest fix:
   run steps 1-3 and 4 in the same job and read from the job's own output folder (`f.path3`).
2. **WDPA version.** PA polygons come from
   `leitoldv/GEDI_global_PA_v2/WDPA_shapefiles/WDPA_polygons/<ISO3>_PA_poly.rds`. Confirm which
   WDPA release these are; the analysis uses WDPA 2025, so they may need rebuilding from
   `GLOBAL_PAs_all_nonMARINE_2025.gpkg`.
3. **Per-country covariates.** `<ISO3>_GLADCover_reclass2000.tif` and `<ISO3>_gmw1996.tif` must
   exist in `abarenblitt/GEDI_global_PA_v2/WDPA_input_vars_GLOBAL/` for every country.
4. **Labels.** `gediwk <- 24` only labels file names now; leave it, or rename consistently.

## Re-registering on DPS (OGC application package)

The scripts only load terra, sf, dplyr, stringr, optmatch, foreach, doParallel, aws.s3 and s3,
so the container can be much leaner than `gedi_pa_env.yml` (which also lists retired packages
such as rgeos). Plan: copy the three R scripts into `matching/r/` (with the source commit noted,
and Veronika's agreement), apply the edits above, add a `run_matching.sh <ISO3>` wrapper that
runs steps 1-3 then 4, and package it like `extract_tile/` (Containerfile from a rocker
geospatial image, `algorithm_config.yml` with one `iso3` input, about 16 GB memory). Test on
Brazil before any other country.

## Order of work

1. **Brazil, existing matches** (no new matching): when the extraction finishes,
   `python analysis/did_country.py --iso3 BRA --matches <folder with Veronika's Brazil RDS files>
   --cells /home/jovyan/GEDI/gedi_pa_change_cells.parquet --out /home/jovyan/GEDI/results`.
   Those matches were built on the old GEDI grid; PA or control cells without data under the new
   filters simply drop out of the join.
2. **Brazil, re-matched** with the patched step 1 and the new container. Compare with step 1:
   this validates the new pipeline.
3. **Other countries**, largest first (Russia, USA, Canada, China, Australia), as soon as their
   tiles are extracted.

## Regression-to-the-mean checks (matching on GEDI L4B)

The Brazil matching includes GEDI L4B biomass as a covariate, so the outcome (GEDI change) and
a matching variable share measurement error. Two checks, both decided before looking at results:

1. **Diagnostic on the existing matches**: `analysis/rtm_check.py` compares T2 - L4B between
   protected and control cells at the same baseline, in undisturbed cells. A gap near 0 means no
   sign of regression to the mean; a negative gap of size G means up to about G of a positive
   AGBD protection effect could be regression to the mean.
2. **Sensitivity run**: re-run Brazil's matching with L4B removed from the propensity model
   (all other covariates and settings unchanged), save to a separate output folder, and run
   `analysis/did_country.py` on both sets of matches. Report both DiD estimates. A similar
   effect with and without L4B means the result does not depend on matching on GEDI.

Optional third check: match on L4B in broad classes (e.g. 50 Mg/ha bins, exact within class)
instead of as a continuous propensity covariate.
