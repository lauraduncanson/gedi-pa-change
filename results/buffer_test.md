# Control buffer test: 2 km vs 10 km (pre-registered)

2026-11-26. Design and decision rule were recorded in CLAUDE.md before running.

**Setup**
- Countries: Costa Rica, Guatemala, Ghana, Cameroon; all their country x biome units.
- Inputs: `matching/prep_country_inputs.py` with `--buffer-km 10` and `--buffer-km 2`. Everything
  else is identical; the treatment cells are the same in both.
- Matching: her V3 step 4 with continuous L4B on DPS v5 (process 250).
- Jobs: 29 (14 units at 10 km, 15 at 2 km), all successful.
- DiD rules as elsewhere: per PA, at least 10 cells per group, PA-cell weights.
- Paired difference: computed over PAs qualifying in both buffers, with a 95% bootstrap interval
  over PAs (4,000 resamples).

**Control pools** (complete cells):

| Country | 10 km | 2 km |
|---|---|---|
| Costa Rica | 1,093 | 5,685 |
| Guatemala | 12,040 | 26,595 |
| Ghana | 53,435 | 86,318 |
| Cameroon | 105,982 | 121,735 |

## Results

| Country | PAs qualifying (10 km / 2 km / both) | AGBD DiD 10 km → 2 km (Mg/ha) | Paired diff [95% CI] | RH98 DiD 10 km → 2 km (m) | Paired diff [95% CI] | T1 AGBD PA − control, 10 km → 2 km |
|---|---|---|---|---|---|---|
| Cameroon | 23 / 22 / 22 | +5.36 → +1.35 | **−3.94 [−8.15, −0.30]** | +0.43 → +0.19 | **−0.28 [−0.61, −0.01]** | 4.3 → 4.0 |
| Ghana | 18 / 23 / 16 | −1.68 → −0.18 | +1.41 [−0.19, +2.14] | −0.24 → +0.02 | +0.24 [−0.03, +0.36] | 8.1 → 4.6 |
| Guatemala | 7 / 15 / 7 | +5.35 → −4.35 | −4.42 [−25.1, +1.3] | +0.53 → −0.54 | −0.55 [−2.2, +1.1] | 10.1 → 9.1 |
| Costa Rica | 0 / 4 / 0 | no PA qualifies → −75.5 | n/a | → −5.4 | n/a | n/a → 65.2 |

AGBD change, PA / control (Mg/ha):

| Country | 10 km | 2 km |
|---|---|---|
| Cameroon | −3.67 / −9.03 | −4.79 / −6.14 |
| Ghana | −3.53 / −1.85 | −3.19 / −3.02 |

## Decision (pre-registered rule)

The rule applies to countries with at least 10 PAs qualifying in both buffers: Ghana and Cameroon.
- **Ghana:** both intervals include 0.
- **Cameroon:** both intervals exclude 0.

**The rule is not met, so the global run keeps the 10 km buffer everywhere.** Units too thin to
match go to the imputation rule (CLAUDE.md, "PAs without a valid counterfactual").

## Reading

1. **Cameroon shows an edge effect.** Controls 2-10 km from a PA lost less biomass than controls
   beyond 10 km (−6.1 vs −9.0 Mg/ha), consistent with positive spillover (buffer zones, management
   reaching past the boundary). Near-PA controls would understate the protection effect by about
   4 Mg/ha there. This is the bias the 10 km rule guards against.
2. **Ghana shows no detectable difference.** 2 km also halves the T1 imbalance there
   (8.1 → 4.6 Mg/ha). This mixed picture is why a blanket 2 km rule is not justified.
3. **Costa Rica cannot be rescued by a smaller buffer.** At 2 km its 4 qualifying PAs have controls
   65 Mg/ha below them at T1 and a meaningless DiD (−75 Mg/ha). Unprotected land there is simply
   not comparable to the PAs, so its PAs belong to the imputation rule.
4. **Guatemala** has too few PAs qualifying in both buffers (7) for a conclusion.
5. **Cameroon's DiD differs from the earlier L4B check** (+5.4 here vs +3.0 there). That check used
   Andrea's 2025 step 1-3 inputs; this one uses the global inputs with ecoregion strata. Compare
   designs within a test, not levels across tests.

**Files:** `results/buffer_test/comparison.csv` and
`results/buffer_test/<b10km|b2km>/<ISO3>/<ISO3>_did_*.csv`.
