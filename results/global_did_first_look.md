# Global DiD: first look (matching v6, 2026-11-27)

Preliminary. The balance rule and imputation (CLAUDE.md, "PAs without a valid counterfactual") are
not applied yet; they belong in `analysis/totals.py`. Weights here are a proxy for PA size.

## Run

- Design: CLAUDE.md "Global matching design" (2020 covariates, continuous L4B, exact land cover
  and WWF ecoregion with her no-ecoregion fallback, country x biome units, 10 km control buffer).
- DPS process 251 (v6), queue `maap-dps-worker-32vcpu-64gb`, driver `dps/global_matching.py`.
- Practice set: 11 jobs (10 units), all successful. It included 6 units that failed on v5 and
  4 that succeeded on v5:
  - ALB_12 and BRA_09 are identical to v5.
  - AFG_13 differs only in 2 PAs that mixed the main and fallback branches. v5 had silently lost
    pairs there (PA 15506: 8 pairs; PA 555705307: 304 v5 pairs against 402 in v6). v5 outputs
    are therefore superseded.
- Full run: **557 jobs submitted, 557 successful, 0 failed**. Matches are in
  `output/global_matches/` (36,985 per-PA files).

## Coverage

- Units: 458 in 157 countries; 36,994 PAs in the plan. **15,398 PAs matched.**
- 29 countries matched no PAs at all (controls too few or not comparable): ARE BGR BRB BRN COM CPV
  CZE DEU GBR GNB HRV HUN IRL ISR JAM MDA MNE MUS PSE PYF QAT SGP SLV SOM STP SVK SWZ SYC TTO.
  Their PAs go to the imputation rule.
- AGBD: **7,444 PAs qualify** (at least 10 cells in each of PA/control x T1/T2) in 116 countries.
  75 countries have at least 10 qualifying PAs. RH98: 7,458 PAs.
- China (49 qualifying PAs) and India cover only the PAs listed in WDPA, which is a small share
  of their protected land.

## Global DiD over qualifying PAs

Each PA is weighted by its PA cells with data in T1 plus T2 (a size proxy). The interval is a 95%
bootstrap over PAs (2,000 resamples).

| Metric | DiD [95% CI] | Change PA / control | T1 PA - control |
|---|---|---|---|
| AGBD (Mg/ha) | **+1.48 [+1.00, +1.99]** | -5.92 / -7.40 | +5.14 |
| RH98 (m) | **+0.16 [+0.12, +0.21]** | -0.42 / -0.58 | +0.55 |

**Brazil drives the AGBD result.** Dropping one country at a time:

| Countries dropped | AGBD DiD (Mg/ha) |
|---|---|
| none | +1.48 |
| Brazil | **+0.11** |
| USA | +1.81 |
| USA, Australia, Canada | +2.70 |

## Largest countries (AGBD, by weight)

| Country | PAs | DiD | T1 PA - ctl | Change PA | Change ctl |
|---|---|---|---|---|---|
| AUS | 733 | +0.08 | +2.36 | -0.44 | -0.53 |
| BRA | 804 | +6.97 | +9.75 | -12.27 | -19.24 |
| USA | 1,294 | -3.06 | +8.13 | -6.16 | -3.09 |
| KAZ | 77 | -1.31 | +0.95 | -1.31 | +0.01 |
| CAN | 516 | -1.33 | -0.73 | -3.26 | -1.93 |
| RUS | 282 | +2.80 | +0.14 | -3.09 | -5.89 |
| MEX | 207 | +2.56 | -0.70 | -2.45 | -5.01 |
| MNG | 84 | -0.17 | +1.61 | -1.09 | -0.92 |
| VEN | 124 | +2.73 | +14.99 | -9.07 | -11.80 |
| ETH | 64 | +0.96 | -1.25 | -0.35 | -1.31 |
| ARG | 144 | -0.47 | +2.09 | -1.24 | -0.78 |
| COD | 48 | -2.15 | +5.42 | -7.79 | -5.65 |
| PER | 86 | +1.96 | +17.89 | -25.17 | -27.13 |
| BOL | 67 | +3.81 | +7.96 | -20.28 | -24.10 |
| IDN | 120 | -1.65 | +10.98 | -34.85 | -33.20 |

Full table: `results/global_did/agbd_by_country_qualifying.csv`. Per-PA and per-country outputs:
`results/global_did/<ISO3>/`. All-metric summary: `results/global_did/summary_by_country.csv`.

## Reading, and cautions

1. Brazil's DiD here (+7.0 Mg/ha) is higher than with Veronika's Brazil inputs (+4.6). This run
   uses global covariates, `lc2020` instead of MapBiomas, and ecoregion strata. Compare designs
   within a test, not levels across tests.
2. T1 balance is poor in several tropical countries (PA - control at T1: VEN +15, PER +18,
   IDN +11, BRA +10 Mg/ha). The |SMD| > 0.25 rule will move many of these PAs to imputation, so
   the global number may change.
3. The USA has a negative DiD (-3.1) with PAs starting 8 Mg/ha above their controls; worth a
   closer look (fire in western PAs is one candidate).
4. Absolute PA change (-5.9 Mg/ha) agrees with the unmatched global PA-change analysis (-6.1).

## Next

`analysis/totals.py`:
- per-PA T1 SMD and the balance rule;
- imputation by biome x country, then region x biome, then biome;
- area weights;
- measured-only and measured + imputed totals, side by side.
