#!/usr/bin/env Rscript
#-------------------------------------------------------------------------------------
# gedi_change_step4.R — Propensity-score matching for GEDI counterfactual analysis
#
# Usage:
#   Rscript gedi_change_step4.R <iso3> [n_cores]
#
# Arguments:
#   iso3     — biome/region code to process (e.g. BPt, BAe, BCa)
#   n_cores  — number of parallel workers (default: 10)
#
# Example:
#   Rscript gedi_change_step4.R BPt 10
#   Rscript /projects/my-public-bucket/GEDI_global_change/gedi_change_brazil/repo/gedi_change_step4_V2_cmd.R BPp 10
#-------------------------------------------------------------------------------------
 
#-------------------------------------------------------------------------------------
# gedi-pa-change patch (DPS): extra arguments
#   Rscript gedi_change_step4_V3_cmd.R <iso3> <covariates> <n_cores> <in_dir> <out_dir>
#   covariates: with_l4b (her V3 model) or no_l4b (gedi_l4b removed from the propensity
#               model; everything else unchanged)
#   in_dir:     folder holding <iso3>_prepped_PAs.RDS and <iso3>_prepped_control.RDS
#   out_dir:    results go to <out_dir>/MATCHING_results_V3/ and logs to <out_dir>/matching_log/
#   part:       optional "i/n" (default 1/1): run only part i of the region's PAs, split into n
#               parts by a fixed greedy rule on PA size. Each PA is matched independently (own
#               seed, full regional control pool), so splitting does not change any match.
#-------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript gedi_change_step4_V3_cmd.R <iso3> <with_l4b|no_l4b> <n_cores> <in_dir> <out_dir> [i/n]",
       call. = FALSE)
}

iso3       <- args[1]
covariates <- args[2]
n_cores    <- as.integer(args[3])
in_dir     <- args[4]
out_dir    <- args[5]
part       <- if (length(args) >= 6) args[6] else "1/1"
if (!covariates %in% c("with_l4b", "no_l4b")) stop("covariates must be with_l4b or no_l4b")
part_i <- as.integer(sub("/.*", "", part)); part_n <- as.integer(sub(".*/", "", part))
if (is.na(part_i) || is.na(part_n) || part_i < 1 || part_i > part_n) stop("part must be i/n with 1 <= i <= n")

# Propensity model covariates, as in matching_func_2025.R; gedi_l4b only for with_l4b.
ps_covars <- c("mean_temp", "max_temp", "min_temp", "prec", "elevation", "slope",
               "d2road", "d2city", "popden", "popcnt", "tt2city")
if (covariates == "with_l4b") ps_covars <- c(ps_covars, "gedi_l4b")
ps_formula <- as.formula(paste("status ~", paste(ps_covars, collapse = " + ")))

cat("=== GEDI change step 4: matching ===\n")
cat("iso3:      ", iso3, "\n")
cat("covariates:", covariates, "\n")
cat("formula:   ", deparse(ps_formula, width.cutoff = 500L), "\n")
cat("n_cores:   ", n_cores, "\n")
cat("part:      ", part, "\n\n")
 
#-------------------------------------------------------------------------------------
# Libraries
#-------------------------------------------------------------------------------------
library(terra)
library(sf)
library(dplyr)
library(stringr)
library(purrr)
library(optmatch)
library(foreach)
library(doParallel)
library(data.table)
 
#-------------------------------------------------------------------------------------
# Paths and spatial inputs
#-------------------------------------------------------------------------------------
# gedi-pa-change patch: outputs and logs under out_dir (f.path is only used for writing below
# and in matching_func_2025.R's logs); inputs from in_dir. The unused region shapefile read
# (adm) is dropped.
f.path <- paste0(normalizePath(out_dir, mustWork = TRUE), "/")

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "matching_func_2025.R"))

#-------------------------------------------------------------------------------------
# Load inputs
#-------------------------------------------------------------------------------------
d_control_local <- readRDS(file.path(in_dir, paste0(iso3, "_prepped_control.RDS")))
d_control_local <- d_control_local[complete.cases(d_control_local), ]

all_PA_points <- readRDS(file.path(in_dir, paste0(iso3, "_prepped_PAs.RDS")))
cat("Control cells:", nrow(d_control_local), "\n")
cat("PAs loaded:   ", length(all_PA_points), "\n\n")
 
if (!dir.exists(paste0(f.path, "MATCHING_results_V3")))
  dir.create(paste0(f.path, "MATCHING_results_V3"))
if (!dir.exists(paste0(f.path, "MATCHING_results_V3/", iso3)))
  dir.create(paste0(f.path, "MATCHING_results_V3/", iso3))
 
#-------------------------------------------------------------------------------------
# Determine which PAs still need to be processed
#-------------------------------------------------------------------------------------
existing_files <- list.files(paste0(f.path, "MATCHING_results_V3/", iso3),
                              pattern = "_matching_results\\.RDS$",
                              full.names = TRUE)
processed_ids  <- str_extract(basename(existing_files), "(?<=_pa_)\\d+")
remaining_ids  <- setdiff(names(all_PA_points), processed_ids)

# gedi-pa-change patch: keep only this job's part of the PAs. Greedy split: PAs sorted by
# number of complete-case cells (largest first, ties by id), each to the part with the fewest
# cells so far (ties to the lowest part). Deterministic, so parts never overlap or miss a PA.
if (part_n > 1) {
  pa_size <- vapply(names(all_PA_points), function(id) {
    d <- all_PA_points[[id]]; if (is.null(d)) 0L else sum(complete.cases(d))
  }, integer(1))
  ord  <- names(pa_size)[order(-pa_size, names(pa_size))]
  load <- numeric(part_n); assign_to <- integer(length(ord)); names(assign_to) <- ord
  for (id in ord) { k <- which.min(load); assign_to[id] <- k; load[k] <- load[k] + pa_size[id] }
  remaining_ids <- intersect(remaining_ids, names(assign_to)[assign_to == part_i])
  cat("part", part, ":", sum(assign_to == part_i), "of", length(ord), "PAs,",
      load[part_i], "of", sum(load), "cells\n")
}

cat(length(processed_ids), "PAs already processed; skipping them.\n")
cat(length(remaining_ids), "PAs remaining to process.\n\n")
 
#-------------------------------------------------------------------------------------
# Run matching for remaining PAs
#-------------------------------------------------------------------------------------
registerDoParallel(n_cores)
cat("Parallel processing with", getDoParWorkers(), "workers\n\n")
 
if (length(remaining_ids) > 0) {
 
  all_new_matches <- foreach(
    i         = seq_along(remaining_ids),
    .combine  = rbind,
    .packages = c("dplyr", "tidyr", "optmatch"),
    .export   = c("propensity_filter", "match_wocat", "foreach_rbind",
                  "f.path", "iso3", "ps_formula")
  ) %dopar% {
 
    id_pa <- remaining_ids[i]
    cat("\nProcessing PA:", id_pa, "\n")
 
    out_file <- paste0(f.path, "MATCHING_results_V3/", iso3, "/",
                       iso3, "_pa_", id_pa, "_matching_results.RDS")
    if (file.exists(out_file)) {
      cat("File already exists for", id_pa, "- skipping.\n")
      return(NULL)
    }
 
    d_pa <- all_PA_points[[id_pa]]
    if (is.null(d_pa) || nrow(d_pa) == 0) {
      cat("No data for PA:", id_pa, "- skipping.\n")
      return(NULL)
    }
    cat(iso3, "pa no.", id_pa, "has", nrow(d_pa), "treatment cells\n")
 
    # Propensity pre-filter
    d_filtered_prop <- tryCatch(
      propensity_filter(d_pa, d_control_local),
      error = function(e) return(NA)
    )
    cat("Propensity filtered DF dimensions:", dim(d_filtered_prop), "\n")
 
    d_wocat_all   <- tryCatch(dplyr::filter(d_filtered_prop,  status), error = function(e) return(NA))
    d_control_all <- tryCatch(dplyr::filter(d_filtered_prop, !status), error = function(e) return(NA))
 
    n_control <- nrow(d_control_all)
    ids_all0  <- tryCatch(d_control_all$UID, error = function(e) return(NA))
    ids_all   <- d_control_all$UID
 
    # Split treatment into chunks of ≤300
    set.seed(125)
    N <- ceiling(nrow(d_wocat_all) / 300)
    l <- tryCatch(
      split(d_wocat_all, sample(1:N, nrow(d_wocat_all), replace = TRUE)),
      error = function(e) return(NULL)
    )
    cat("Number of chunks:", length(l), "\n")
 
    if (length(l) > 0) {
      pa_match <- data.frame()
 
      for (pa_c in seq_along(l)) {
        ids_all       <- d_control_all$UID
        d_wocat_chunk <- l[[pa_c]]
        cat("Chunk", pa_c, "of PA", id_pa, "\n")
 
        n_treatment <- nrow(d_wocat_chunk)
        t           <- ifelse(floor(n_control / n_treatment) <= 7,
                              ifelse(floor(n_control / n_treatment) < 1, 1,
                                     floor(n_control / n_treatment)), 7)
        n_sample    <- round(n_treatment * t)
        m_all2_out  <- data.frame()
        n_matches   <- 0
 
        tryCatch(
          while (n_matches < n_treatment) {
            n_ids <- length(ids_all)
            cat("n ids:", n_ids, "\n")
 
            if (n_ids > n_sample) {
              set.seed(125)
              sample_ids_bar   <- sample(ids_all,  n_sample)
              sample_ids       <- sample(ids_all0, n_sample)
              d_control_sample <- d_control_all[d_control_all$UID %in% sample_ids, ]
              ids_all          <- setdiff(ids_all, sample_ids_bar)
 
              new_d  <- tryCatch(rbind(d_wocat_chunk, d_control_sample), error = function(e) return(NULL))
              m_all  <- tryCatch(match_wocat(new_d, pid = id_pa),        error = function(e) return(NULL))
              m_all2 <- tryCatch(m_all[1, ],                             error = function(e) return(NULL))
 
              n_matches_temp <- tryCatch(nrow(m_all2$df), error = function(e) return(NULL))
 
              if (!is.null(n_matches_temp)) {
                m_all2$df$pa_id <- rep(id_pa, n_matches_temp)
                m_all2_out      <- rbind(m_all2_out, m_all2$df)
 
                matched_protected <- dplyr::filter(m_all2$df, status == TRUE)
                matched_control   <- dplyr::filter(m_all2$df, status == FALSE)
                cat("matched_protected:", nrow(matched_protected), "\n")
 
                n_matches     <- n_matches + nrow(matched_protected)
                d_wocat_chunk <- d_wocat_chunk[-(match(matched_protected$UID, d_wocat_chunk$UID)), ]
              }
              ids_all0 <- setdiff(ids_all0, matched_control$UID)
 
            } else {
              n_treatment <- n_matches
            }
          },
          error = function(e) return(NULL)
        )
 
        match_score <- m_all2_out
        cat(table(match_score$status), "\n")
        pa_match <- rbind(pa_match, match_score)
      }
 
    } else {
      cat("No chunks produced for PA", id_pa, "\n")
      pa_match <- NULL
    }
 
    saveRDS(pa_match, out_file)
    cat("Saved results for PA", id_pa, "\n")
    return(pa_match)
  }
 
} else {
  cat("All PAs already processed. Skipping matching step.\n")
  all_new_matches <- NULL
}
 
stopImplicitCluster()
 
#-------------------------------------------------------------------------------------
# Combine all per-PA outputs
#-------------------------------------------------------------------------------------
cat("\nCombining all PA results...\n")
 
all_result_files <- list.files(paste0(f.path, "MATCHING_results_V3/", iso3),
                                pattern = "_matching_results\\.RDS$",
                                full.names = TRUE)
 
all_matches <- purrr::map_dfr(all_result_files, readRDS, .id = "source")
 
out_combined <- paste0(f.path, "MATCHING_results_V3/", iso3, "_allPA_matched_results.RDS")
saveRDS(all_matches, file = out_combined)
cat("Combined results saved to", out_combined, "\n")
 
cat("Number of PAs:", length(unique(all_matches$pa_id)), "\n")
cat("Total rows:   ", nrow(all_matches), "\n")
cat("\nRows by status:\n")
print(table(all_matches$status))
 
cat("\n=== Done ===\n")
