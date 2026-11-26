# Cameroon: how L4B enters the matching (pre-registered check)

2026-11-26. Design recorded in CLAUDE.md ("L4B matching check in Cameroon") before any run:
- Same code for all arms: Veronika's Brazil V3 step 4, DPS process 249.
- Same inputs: the 2025 global-pipeline step 1-3 outputs for CMR, with `gedi_l4b` sampled from
  the same global raster used for Brazil.
- Same DiD rules as Brazil: per PA, at least 10 cells in each PA/control x T1/T2 group, PAs
  weighted by PA cells.
- Ghana was the first choice. It was replaced before any matching, after a covariate-only overlap
  check (Ghana has almost no high-biomass control cells).

**Arms**
- `with_l4b`: L4B as a continuous propensity covariate (her V3 model).
- `no_l4b`: L4B dropped.
- `l4b_classes`: L4B dropped from the propensity model; exact matching within 50 Mg/ha L4B
  classes.

**Jobs**
- 6 matching jobs (3 arms x 2 parts), about 3 minutes each; all successful.
- One Liberia practice job ran first. It confirmed the class matching is exact: PA and control
  counts are equal in all 90 PA x class cells.

## Results

95% intervals come from a bootstrap over PAs (4,000 resamples; same PA-cell weights). 26 of the
qualifying PAs are common to all three arms.

| Arm | PAs qualifying | AGBD DiD (Mg/ha) | 95% CI | AGBD median DiD | T1 AGBD, PA − control | RH98 DiD (m) | 95% CI | T1 RH98, PA − control | RTM gap (Mg/ha) |
|---|---|---|---|---|---|---|---|---|---|
| with_l4b | 26 / 44 | **+2.96** | −3.2 to +8.1 | +6.64 | **+0.94** | **+0.31** | −0.07 to +0.62 | +0.03 | −4.24 |
| no_l4b | 27 / 42 | **−0.51** | −8.2 to +5.3 | +1.62 | **+6.21** | **+0.09** | −0.33 to +0.44 | +0.29 | −3.49 |
| l4b_classes | 27 / 41 | **+3.12** | −2.1 to +8.0 | +3.00 | **−1.72** | **+0.33** | +0.04 to +0.64 | −0.37 | −5.35 |

AGBD changes, PA / control (Mg/ha):

| Arm | PA change | Control change |
|---|---|---|
| with_l4b | −4.40 | −7.36 |
| no_l4b | −4.23 | −3.72 |
| l4b_classes | −4.33 | −7.45 |

The RTM gap is computed over all matched cells, because there are no disturbance layers for
Cameroon.

## Pre-registered decision

**Rule:** the classes arm is preferred if its T1 balance is at least as good as `with_l4b` and
its RTM gap is closer to 0.

**Outcome: not met on either criterion.**
- T1 balance: |−1.72| against |+0.94| Mg/ha for AGBD, and |−0.37| against |+0.03| m for RH98.
  Both are small.
- RTM gap: −5.35 against −4.24.

So the rule does not move us off her continuous-L4B model.

## Reading

1. **Cameroon repeats Brazil's pattern.** Dropping L4B lowers the AGBD DiD (+3.0 → −0.5 here;
   +4.6 → +2.0 in Brazil) and worsens T1 balance. Controls start 6.2 Mg/ha below the PAs here,
   15.4 in Brazil, against about 1 and 7 with L4B. The PA change is the same in every arm; only
   the controls move. Lower-biomass controls lose less biomass in absolute terms (−3.7 against
   −7.4 Mg/ha), which pulls the no-L4B DiD down.
2. **Coarse L4B classes give the same answer as continuous L4B** (+3.1 vs +3.0 Mg/ha AGBD; +0.33
   vs +0.31 m RH98). Classes control baseline biomass only to 50 Mg/ha and do not match on fine
   L4B differences, which is where L4B measurement error would act. If regression to the mean
   from matching on L4B error were inflating the effect, the classes arm should have shrunk it.
   It did not. **This points to the no-L4B estimate being biased low by baseline imbalance,
   rather than the with-L4B estimate being inflated by RTM.**
3. **The RTM check does not isolate RTM.** Its gap is negative in all three arms, including the
   one that does not use L4B at all (−3.5 Mg/ha). It also was in Brazil (−2.5 without L4B). At
   the same baseline it mixes RTM with real differences in how PA and control cells change, so
   it should be read only as an upper bound, as its docstring says, and not used to choose a
   design.
4. **Cameroon alone is too small to settle the size of the effect.** It has 26-27 qualifying
   PAs, and every AGBD interval spans zero. The comparison between arms is the informative part,
   not the levels. RH98 with classes (+0.33 m, 95% CI +0.04 to +0.64) is the only interval
   excluding zero.

## Implication for Brazil

With the evidence so far, the with-L4B (or L4B-class) estimate, +4.6 Mg/ha, is better supported
than the no-L4B one, +2.0, which comes with twice the T1 imbalance.

A Brazil `l4b_classes` run would test this directly in the country that matters. It is cheap:
7 regions with process 249, about 15-30 minutes. It should be decided as its own pre-registered
step, because CLAUDE.md fixed "no further arms" for this Cameroon check.

## Files

- `results/CMR_l4b_check/<arm>/`: `CMR_did_by_pa.csv`, `CMR_did_summary.csv`, `CMR_rtm_bins.csv`,
  `CMR_rtm.png`.
- Matches:
  `s3://maap-ops-workspace/lduncanson/dps_output/gedi-pa-change-match-brazil_14/v4/2026/11/26/18/…`.
- Converted inputs: `s3://maap-ops-workspace/lduncanson/gedi_pa_change/matching_inputs/CMR/`,
  made by `matching/prep_global_inputs.py`.
