# Global change inside protected areas, without matching (first look)

2026-11-26. Secondary result: absolute T2 − T1 change inside PAs. Rules were fixed in CLAUDE.md
("Global PA change without matching") before running.

**Data and rules**
- **PAs:** WDPA 2025 non-marine, Designated, Inscribed or Established, any year; no OECMs (none
  in the file); no UNESCO-MAB Biosphere Reserves.
- **Cells:** a 1 km cell belongs to a PA if its centre is inside the polygon. Per-PA mean of cell
  means.
- **Minimum:** at least 10 cells with data in each epoch, per metric.
- **GEDI:** V003, power beams, growing-season windows. T1 = 2019-2021; T2 = 2024 to July 2025.
- **Code:** `analysis/pa_change_global.py` (assign: 22 min on the Hub; summarize: 24 s).

**Coverage**
- 53.2 M 1 km cells have GEDI data. 8.92 M of them (16.8%) fall inside some WDPA site.
- 77,400 PAs have data; **20,645 qualify** for RH98 (20,508 for AGBD).

## Global (PAs weighted by PA cells)

| Metric | T1 | T2 | **Change** | Median change | PAs declining |
|---|---|---|---|---|---|
| AGBD (Mg/ha) | 82.4 | 76.3 | **−6.12** | −1.91 | 63% |
| RH98 (m) | 13.62 | 13.16 | **−0.45** | −0.25 | 61% |
| RH50 (m) | 5.32 | 5.10 | −0.23 | −0.05 | 60% |
| PAI | 1.77 | 1.69 | −0.07 | −0.02 | 60% |
| FHD | 2.23 | 2.20 | −0.03 | −0.03 | 65% |
| Cover | 0.372 | 0.365 | −0.006 | −0.004 | 57% |

## Context: the same change outside PAs

This is not a counterfactual. Cells outside every WDPA site are reweighted to the PAs' mix of
GEDI region x PFT strata: 36 strata, each with at least 100 cells per group and epoch.

| Metric | PA land | Outside | PA − outside |
|---|---|---|---|
| AGBD (Mg/ha) | −7.12 | −5.82 | −1.30 |
| RH98 (m) | −0.51 | −0.46 | −0.05 |
| RH50 (m) | −0.29 | −0.25 | −0.04 |
| PAI | −0.083 | −0.075 | −0.008 |
| FHD | −0.032 | −0.038 | +0.007 |
| Cover | −0.008 | −0.009 | +0.001 |

Largest strata, change in AGBD (Mg/ha), PA / outside:

| Stratum | PA | Outside |
|---|---|---|
| S America evergreen broadleaf | −16.5 | −12.7 |
| Africa EBT | −11.9 | −10.1 |
| S Asia EBT | −16.2 | −11.1 |
| Australasia EBT | −21.2 | −31.0 |
| N America ENT | −6.3 | −3.6 |

Full table: `inside_vs_outside_by_region_pft.csv`.

## Reading, with caveats

1. **Most of the decline inside PAs is shared with unprotected land.** RH98 falls by about 0.5 m
   both inside and outside PAs. Absolute change therefore reflects whatever affected the whole
   landscape between the epochs, and that may include:
   - **Real events.** The 2023-24 El Niño drought (Amazon, SE Asia) falls between the epochs. In
     Australia, T1 includes pre-fire 2019 shots before the 2019-20 Black Summer fires.
   - **Sampling and season differences.** T2 is only 2024 to July 2025, so it has fewer cells
     (4.9 M vs 8.7 M PA cell-epochs). Jan-Apr 2024 has no data at all, which affects windows
     early in the year.
   - **Any time-dependent instrument bias.**

   This is why CLAUDE.md makes the matched DiD the headline and this the secondary result.
2. **The PA − outside differences are not a protection effect.** They are unmatched. PAs sit in
   wetter, higher-biomass, more remote parts of each stratum, and higher-biomass forest loses more
   in absolute terms. The matched DiD (Brazil +4.6 Mg/ha; Cameroon +3.0) handles that.
3. **By country,** among those with the most PA cells (AGBD change, median):

   | Country | Change | Median | PAs declining |
   |---|---|---|---|
   | Brazil | −12.2 | −5.2 | 73% |
   | DR Congo | −9.7 | −8.0 | |
   | Venezuela | −7.9 | | |
   | USA | −4.9 | | |
   | Russia | −4.0 | | |
   | Canada | −2.8 | | |
   | Australia | −0.9 | | |

   Full table: `summary_by_country.csv`.
4. **Standard errors are not shown yet.** Per-PA SEs in `pa_change_by_pa.parquet` assume
   independent cells (a lower bound). Country and global uncertainty needs a bootstrap over PAs,
   to add with the area-weighted totals (`analysis/totals.py`).

**Files**, in `results/pa_change_global/`:
- `pa_change_by_pa.parquet`: per PA x metric (the CSV version is git-ignored, 107 MB)
- `summary_global.csv`
- `summary_by_country.csv`, `summary_by_iucn.csv`, `summary_by_designation_type.csv`,
  `summary_by_governance.csv`
- `summary_region_pft_inside_outside.csv`
- `inside_vs_outside_*.csv`

The cell-to-PA table is `output/cell_pa.parquet` (60 MB, not committed).
