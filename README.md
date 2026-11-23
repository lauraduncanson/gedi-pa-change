# gedi-pa-change

GEDI-based 5-year change in forest structure and biomass in protected areas and their matched
counterfactuals. T1 = 2019-2021, T2 = 2024-2026. Runs on NASA MAAP (Hub + DPS).

## Layout

| Path | What it is |
|---|---|
| `extract_tile/` | The DPS algorithm: one 1-degree tile of `gedi-tiled-v3` -> filtered shots and 1 km cell x epoch means. Includes the growing-season windows (`doy_windows.csv`), the Containerfile and the OGC algorithm config. |
| `dps/` | Submit, monitor and collect the per-tile jobs from the Hub. |
| `notebooks/` | `plot_tile_histograms.py` for T1 vs T2 histograms by tile and protected area. |
| `.github/workflows/` | Builds the container and generates the CWL on every push to `extract_tile/`. |

## Filters (all in `extract_tile/extract_tile.py`)

Power beams only; `l2a_quality_flag_rel3 = 1`; not degraded; leaf-on; inside the growing-season
DOY window for the shot's PFT x region. PAI/FHD/cover need the L2B flag, AGBD the L4A flag.

## One-time setup

1. Create the GitHub repo, replace every `EDIT_ME` in `extract_tile/algorithm_config.yml`, push.
2. The `extract-tile` GitHub Action builds the image (GHCR) and commits the CWL to
   `extract_tile/cwl_workflows/`. Make the GHCR package public so DPS can pull it.
3. Register the algorithm from the Hub (use the *raw* URL of the committed CWL):

   ```python
   from maap.maap import MAAP
   maap = MAAP()
   r = maap.deploy_algorithm(execution_unit_href="https://raw.githubusercontent.com/<org>/gedi-pa-change/main/extract_tile/cwl_workflows/<file>.cwl")
   print(r.status_code, r.json())      # note the process "id"
   ```

## Running

```bash
# 1. two test tiles
python dps/submit_tiles.py --process-id <id> --queue <queue> --tiles S13_E028 N11_W005
python dps/monitor.py

# 2. ramp up, then everything that is left
python dps/submit_tiles.py --process-id <id> --queue <queue> --limit 100
python dps/submit_tiles.py --process-id <id> --queue <queue>

# 3. retry whatever did not finish
python dps/submit_tiles.py --process-id <id> --queue <queue> --resubmit

# 4. one table for analysis
python dps/collect.py --out /home/jovyan/GEDI/gedi_pa_change_cells.parquet
```

A tile is "done" when `summary_<tile>.json` exists under the DPS output folder, so every step is
safe to re-run. When new GEDI data arrives, bump `--tag`, point `--done-dir` at a fresh folder or
clear the old outputs, and run again.

## Cell ids

`cell_id = row * 34704 + col` on the GEDI 1 km EASE2 grid (same grid as
`GEDI_ANCI_PFT_r1000m_EASE2.0`). In R/terra the same cell is `cell_id + 1`.
