#!/bin/bash
# DPS entry point: Veronika Leitold's Brazil V3 matching (step 4) for one biome region.
#   run_matching.sh --region BPp --covariates no_l4b --n_cores 8 [--part 1/2]
# --part i/n runs only part i of the region's PAs (see gedi_change_step4_V3_cmd.R); results
# for each PA are identical however the region is split.
# Downloads <region>_prepped_{PAs,control}.RDS from her MATCHING_points_V3 folder (read only),
# runs matching/r/brazil_v3/gedi_change_step4_V3_cmd.R, and leaves everything in ./output:
#   output/MATCHING_results_V3/<region>/<region>_pa_<WDPAID>_matching_results.RDS
#   output/MATCHING_results_V3/<region>_allPA_matched_results.RDS
#   output/matching_log/, output/step4_log.txt, output/session_info.txt, output/run_info.json
set -uo pipefail

region=""; covariates=""; n_cores="4"; part="1/1"
input_prefix="s3://maap-ops-workspace/shared/leitoldv/GEDI_global_change/MATCHING_points_V3"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --region) region="$2"; shift 2 ;;
    --covariates) covariates="$2"; shift 2 ;;
    --n_cores) n_cores="$2"; shift 2 ;;
    --part) part="$2"; shift 2 ;;
    --input_prefix) input_prefix="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
if [[ -z "$region" || -z "$covariates" ]]; then
  echo "usage: run_matching.sh --region <code> --covariates <with_l4b|no_l4b> [--n_cores N]" >&2
  exit 2
fi

# Brazil: her MATCHING_points_V3 folder (default). Other countries: the folder written by
# matching/prep_global_inputs.py (<prefix>/<ISO3>/<ISO3>_prepped_{PAs,control}.RDS).
SRC="${input_prefix%/}"
[[ "$SRC" == *"/MATCHING_points_V3" ]] || SRC="$SRC/$region"
mkdir -p output input
echo "input=$SRC region=$region covariates=$covariates n_cores=$n_cores part=$part nproc=$(nproc) mem_total=$(awk '/MemTotal/{printf "%.1f GB", $2/1048576}' /proc/meminfo)"

for kind in PAs control; do
  python3 /app/s3_get.py "$SRC/${region}_prepped_${kind}.RDS" input/ || exit 3
done
Rscript -e 'sessionInfo()' > output/session_info.txt 2>&1

start=$(date +%s)
Rscript /app/gedi_change_step4_V3_cmd.R "$region" "$covariates" "$n_cores" input output "$part" 2>&1 | tee output/step4_log.txt
status=${PIPESTATUS[0]}
seconds=$(( $(date +%s) - start ))
n_pa_files=$(find output/MATCHING_results_V3 -name '*_matching_results.RDS' 2>/dev/null | wc -l)
printf '{"region": "%s", "covariates": "%s", "part": "%s", "n_cores": %s, "nproc": %s, "seconds": %s, "exit": %s, "n_pa_files": %s}\n' \
  "$region" "$covariates" "$part" "$n_cores" "$(nproc)" "$seconds" "$status" "$n_pa_files" > output/run_info.json
cat output/run_info.json
exit "$status"
