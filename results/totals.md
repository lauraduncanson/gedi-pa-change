# PA totals: measured and imputed DiD (first run, 2026-11-27)

`analysis/totals.py` with the rules in CLAUDE.md ("PAs without a valid counterfactual", "Totals").
Inputs:
- Global matching v6 (`results/global_did/`, rerun with T1 SDs).
- **Uncorrected** extraction (`main`). The drift-corrected extraction (drift-v1) is pending; rerun
  everything on `rh98_c` / `agbd_c` when it is done.

Outputs in `results/totals/`:
- `pa_table.parquet`: per PA x metric.
- `means_{global,country,biome,region}.csv`: area-weighted means, for figures.
- `totals_agbd_*.csv`: Tg.

## PA universe and areas

- **77,400 kept PAs** have at least one GEDI cell. Dissolved WDPA area within 51.6 S-51.6 N:
  1,821 Mha before removing overlaps, **1,515 Mha after**. Overlap scaling is strongest in
  Europe (Slovenia 0.31, Netherlands 0.43, Spain 0.44: Natura 2000 over national designations).
- 2,418 PAs have zero area inside the band: mostly Germany, UK and Poland, with a few GEDI cells
  just north of 51.6 N. They hold 0.5% of PA cells and get no weight.
- 897 PAs have no WWF biome (cells on lakes or rock/ice). They are excluded unless measured.

## Who is measured (AGBD)

| Step | PAs |
|---|---|
| Universe | 77,400 |
| Matched | 15,398 |
| >= 10 cells in all four groups | 7,355 |
| ...of which pass T1 balance (abs. AGBD SMD <= 0.25) = **measured** | **3,858** (50% of PA area) |
| Imputed: biome x country / region x biome / biome | 42,468 / 28,747 / 1,425 |
| Excluded (no biome or no stratum with >= 5 measured PAs) | 902 |

**The balance rule removes almost half of the qualifying PAs (3,497 of 7,355).** Measured PAs
are large, so half of the PA area is still measured.

## Global (area-weighted; 95% bootstrap over measured PAs, imputation redone in each resample)

| Metric | DiD measured only | DiD measured + imputed | PA / control change (measured) | T1 PA / control (measured) |
|---|---|---|---|---|
| AGBD (Mg/ha) | **+2.12 [1.14, 3.08]** | **+2.16 [1.28, 3.07]** | -5.82 / -7.94 | 81.5 / 80.2 |
| RH98 (m) | +0.22 [0.14, 0.29] | +0.21 [0.14, 0.28] | -0.41 / -0.62 | 13.7 / 13.4 |
| RH50 (m) | +0.13 [0.07, 0.18] | +0.12 [0.07, 0.17] | -0.25 / -0.38 | |
| Cover | +0.006 [0.004, 0.008] | +0.005 [0.003, 0.008] | | |

**Total AGBD DiD** (DiD x area):
- Measured PAs (758 Mha): **1.61 Pg [0.87, 2.34]**.
- Measured + imputed (1,515 Mha): **3.26 Pg [1.93, 4.65]**.

This is the biomass loss avoided over 5 years relative to matched controls (uncorrected).

## By GEDI region (AGBD, Tg, measured + imputed)

| Region | Area (Mha) | Measured PAs | Total [95% CI] |
|---|---|---|---|
| South America | 471 | 891 | **2,545 [1,855, 3,251]** |
| South Asia | 64 | 181 | 438 [113, 792] |
| North Asia | 163 | 366 | 303 [29, 665] |
| Europe | 102 | 490 | 274 [-2, 489] |
| Africa | 406 | 576 | 74 [-540, 758] |
| North America | 114 | 944 | -159 [-385, 9] |
| Australasia | 197 | 410 | -211 [-875, 389] |

South America is 78% of the global total; Brazil alone is 1,555 Tg [1,122, 2,033].
- **USA is negative:** -197 Tg [-405, -38]; PAs lost more than their controls.
- **Africa's total is near zero with a wide interval,** despite 406 Mha.

Biome means (AGBD, measured only):
- Tropical moist forests (biome 1): +4.4 Mg/ha [1.8, 6.7].
- Mangroves (14): +6.5.
- Mediterranean (12): +3.2 [1.5, 5.0].
- Boreal/temperate conifer and dry forests: around 0.

## Cautions

1. **Uncorrected for ground drift.** Drift lowers T2 RH and AGBD most under dense canopy. PAs
   have slightly higher T1 structure than their controls, so the correction may raise PA change a
   little more than control change. Its effect on the DiD needs the drift-v1 run, not a guess.
2. **Half the area is imputed.** Imputed PAs take the mean of measured PAs in their stratum. Where
   the measured PAs are unrepresentative (large forested PAs measured, small PAs imputed), the
   imputed share inherits that.
3. **The balance rule is strict at the cell level.** SMD uses the SD of 1 km cell means, and
   about 48% of qualifying PAs fail it. That is as pre-registered. A looser threshold is a
   sensitivity analysis to decide on before looking, not after.
4. **Bootstrap intervals treat PAs as independent**; nearby PAs sharing controls are correlated.
5. Single-PA groups (e.g. Saudi Arabia, 1 measured PA) have degenerate intervals.
