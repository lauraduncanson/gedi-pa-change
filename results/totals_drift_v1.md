# PA totals with ground-drift correction (drift-v1, 2026-11-27)

Same design as `results/totals.md`, on the drift-v1 cell table
(`output/gedi_pa_change_cells_drift_v1_combined.parquet`):
- **Versions:** uncorrected, corrected (central), and the low and high drift bounds (CLAUDE.md,
  "Ground-drift correction").
- **Measured PAs:** balance is always computed on uncorrected T1 AGBD, so the same PAs are
  measured in every version.
- **Thresholds:** headline |SMD| <= 0.25; sensitivity at 0.1 and 0.5.
- **Check:** uncorrected DiD rows reproduce `results/global_did` exactly (15,708 PA rows).

Outputs:
- `results/totals_drift_v1/` (headline), `results/totals_drift_v1_smd0.1/`,
  `results/totals_drift_v1_smd0.5/`.
- `results/global_did_drift_v1/<ISO3>/` (per-PA DiD for all metrics).

## Global, headline (|SMD| <= 0.25; 3,858 measured PAs, 50% of 1,515 Mha)

DiD in Mg/ha (AGBD) and m (RH98), area-weighted.

| Version | DiD measured only | DiD measured + imputed | PA / control change (measured) | Total DiD, measured + imputed (Pg) |
|---|---|---|---|---|
| AGBD uncorrected | +2.12 [1.14, 3.08] | +2.15 [1.28, 3.07] | -5.82 / -7.94 | 3.26 [1.93, 4.65] |
| **AGBD corrected** | **+2.25 [1.31, 3.19]** | **+2.27 [1.45, 3.14]** | **-3.30 / -5.54** | **3.44 [2.19, 4.75]** |
| AGBD low bound | +2.21 [1.24, 3.19] | +2.24 [1.29, 3.12] | -4.05 / -6.26 | 3.39 [1.96, 4.72] |
| AGBD high bound | +2.29 [1.27, 3.34] | +2.31 [1.49, 3.26] | -2.49 / -4.77 | 3.50 [2.26, 4.93] |
| RH98 uncorrected | +0.22 [0.14, 0.29] | +0.21 [0.14, 0.28] | -0.41 / -0.62 | |
| **RH98 corrected** | **+0.23 [0.15, 0.30]** | **+0.21 [0.13, 0.29]** | **-0.24 / -0.46** | |

Measured PAs only (758 Mha): AGBD total 1.61 Pg uncorrected -> 1.71 Pg corrected.

## Balance-threshold sensitivity (AGBD, measured + imputed)

| \|SMD\| threshold | Measured PAs | Uncorrected DiD / total | Corrected DiD / total |
|---|---|---|---|
| 0.1 | 1,755 | +1.39 [0.34, 2.67] / 2.08 Pg | +1.47 [0.35, 2.80] / 2.19 Pg [0.53, 4.19] |
| **0.25** | **3,858** | **+2.15 / 3.26 Pg** | **+2.27 / 3.44 Pg [2.19, 4.75]** |
| 0.5 | 5,827 | +1.99 [1.25, 2.89] / 3.01 Pg | +2.15 [1.33, 3.07] / 3.25 Pg [2.02, 4.64] |

## By region and country (AGBD total, measured + imputed, Tg; uncorrected -> corrected [low, high])

| Group | Uncorrected | Corrected | Low / high |
|---|---|---|---|
| South America | 2,545 | 2,708 | 2,660 / 2,756 |
| South Asia | 438 | 460 | 454 / 466 |
| North Asia | 303 | 307 | 305 / 309 |
| Europe | 274 | 272 | 271 / 274 |
| Africa | 74 | 83 | 82 / 87 |
| North America | -159 | -167 | -165 / -168 |
| Australasia | -211 | -222 | -219 / -225 |
| Brazil | 1,555 | 1,671 | 1,637 / 1,705 |
| Venezuela | 371 | 386 | 382 / 390 |
| USA | -197 | -203 | -202 / -205 |
| Indonesia | -406 | -413 | -410 / -416 |

## Reading

1. **The correction mostly changes absolute change, not the DiD.** It removes about 40% of
   the apparent AGBD loss in measured PAs (-5.8 -> -3.3 Mg/ha) and in their controls
   (-7.9 -> -5.5). PAs and controls are matched on L4B and have similar cover, so the drift
   largely cancels in the difference.
2. **The DiD rises slightly**, +2.15 -> +2.27 Mg/ha (total 3.26 -> 3.44 Pg). PAs are a little
   denser than controls, so they get slightly more correction. The drift-bound range
   (3.39-3.50 Pg) is small compared with the sampling interval (2.2-4.8 Pg).
3. **RH98 DiD is essentially unchanged** (+0.21 m), while absolute RH98 loss falls from -0.41 to
   -0.24 m in PAs.
4. **The result holds across balance thresholds** (corrected 2.2-3.4 Pg). The strictest
   threshold (0.1) measures only 28% of PA area, so its estimate is lower and less precise.
5. **South America carries ~80% of the global total** (Brazil ~49%). North America, Australasia
   and Indonesia are negative: there, PAs lost more biomass than matched controls.
6. Absolute changes still include any drift not captured by the correction. PAI, FHD and cover
   are not corrected (CLAUDE.md).

Caveats from `results/totals.md` still apply: half the area is imputed; the bootstrap treats PAs
as independent.
