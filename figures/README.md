# Figures

Run order (paths are set in `fig_setup.R`, or via environment variables FIG_INPUTS, FIG_DRIFT,
FIG_OUTDIR, FIG_CACHE):

1. `python analysis/did_country.py ...` per country (writes `<ISO3>_did_by_pa.csv`, now including
   the drift-corrected metrics when the cell table comes from the drift-v1 extraction).
2. `python analysis/totals.py run --did results/global_did_drift_v1 --out results/totals_drift_v1`
   (measured/imputed DiD and totals; CLAUDE.md "Totals").
3. `python analysis/figure_inputs.py --totals results/totals_drift_v1 --out /home/jovyan/GEDI/results/figure_inputs`
   (adds PA stock and stock change over total PA area; writes the files `fig_setup.R` reads).
4. Notebooks (R kernel), from this folder:
   - `01_maps.ipynb`: Fig 1 DiD map, Fig 2 stock and stock-change maps, Fig 3 counterfactual availability map
   - `02_summaries.ipynb`: Fig 4 top countries, Fig 5 continent x biome PA vs counterfactual change,
     Fig 6 global totals and drift sensitivity, Fig 7 DiD vs starting biomass
   - `03_drift.ipynb`: Fig S1 ground-finding drift
   Without an R kernel, run the same code with `Rscript scripts/01_maps.R` etc.

Figures are written as PNG (300 dpi) and PDF to FIG_OUTDIR.
