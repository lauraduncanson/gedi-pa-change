options(warn=-1)
options(dplyr.summarise.inform = FALSE)

`%notin%` <- Negate(`%in%`)

#--------------------------------------------------------------------------------
# Function to allow rbinding dataframes with foreach even when some dataframes may not have any rows
foreach_rbind <- function(d1, d2) {
  if (is.null(d1) & is.null(d2)) {
    return(NULL)
  } else if (!is.null(d1) & is.null(d2)) {
    return(d1)
  } else if (is.null(d1) & !is.null(d2)) {
    return(d2)
  } else  {
    return(rbind(d1, d2))
  }
}

#--------------------------------------------------------------------------------
match_wocat <- function(df, pid) {

  options("optmatch_max_problem_size"=Inf)

  # gedi-pa-change patch (performance only, results unchanged): the original ran the body
  # below once per land-cover value, in parallel on 4 inner workers
  #   ret <- foreach(this_lc = unique(df$land_cover), .combine = foreach_rbind, ...) %dopar% {...}
  # but the body never uses this_lc, so every iteration computed the same result, and the
  # caller keeps only ret[1, ]. The body now runs once; the result is then combined exactly as
  # foreach would have combined the identical copies (see the end of this function).
  n_lc <- length(unique(df$land_cover))
  match_once <- function() {

                    this_d <- df
                    d_wocat <- filter(this_d, status)
                    # Filter out land covers that don't appear in the wocat
                    # sample, and drop these levels from the factors
                    this_d <- filter(this_d,
                                     land_cover %in% unique(d_wocat$land_cover),
                                     wwfbiom %in% unique(d_wocat$wwfbiom),
                                     wwfecoreg %in% unique(d_wocat$wwfecoreg))
                    
                    this_d$land_cover <- droplevels(this_d$land_cover)
                    #this_d$wwfbiom <- droplevels(this_d$wwfbiom)
                    #this_d$wwfecoreg <- droplevels(this_d$wwfecoreg)
                    # gedi-pa-change patch: when biome/ecoregion are factors (global inputs), drop
                    # unused levels as the commented lines above intended; otherwise nlevels() counts
                    # absent levels, strata() is added for a one-level factor and glm fails
                    # ("contrasts can be applied only to factors with 2 or more levels"). Numeric
                    # columns (her Brazil inputs) are left unchanged.
                    if (is.factor(this_d$wwfbiom)) this_d$wwfbiom <- droplevels(this_d$wwfbiom)
                    if (is.factor(this_d$wwfecoreg)) this_d$wwfecoreg <- droplevels(this_d$wwfecoreg)
                    # table(this_d$status)
                    dat <- dplyr::select(this_d, lat, lon, UID, status, land_cover, wwfbiom, wwfecoreg, 
                                         elevation, slope, mean_temp, max_temp, min_temp, prec,
                                         d2road, d2city, tt2city, popden, popcnt, gedi_l4b) 
                    # gedi-pa-change patch: ps_formula (set in the step 4 script) replaces the
                    # hard-coded formula here and below; identical for covariates = with_l4b.
                    ps_A <- glm(ps_formula, data = dat)
                    dat$propensity_scoreA <- fitted(ps_A)
                    f <- ps_formula
                    this_d <- bind_cols(this_d, propensity_scoreA=dat$propensity_scoreA)
                    # Can't stratify by land cover or climate if they only have one level
                    if (nlevels(this_d$land_cover) >= 2) {
                      f <- update(f, ~ . + strata(land_cover))
                    } else {
                      f <- update(f, ~ . - land_cover)
                    }
                    if (nlevels(this_d$wwfbiom) >= 2) {
                      f <- update(f, ~ . + strata(wwfbiom))
                    } else {
                      f <- update(f, ~ . - wwfbiom)
                    }
                    if (nlevels(this_d$wwfecoreg) >= 2) {
                      f <- update(f, ~ . + strata(wwfecoreg))
                    } else {
                      f <- update(f, ~ . - wwfecoreg)
                    }
                    # gedi-pa-change patch (l4b_classes only): exact matching within L4B classes.
                    # Keep classes that have both PA and control cells (others cannot be matched
                    # exactly), then stratify on the class like land_cover.
                    if (isTRUE(get0("l4b_strata", ifnotfound = FALSE))) {
                      shared <- intersect(unique(as.character(this_d$l4b_class[this_d$status])),
                                          unique(as.character(this_d$l4b_class[!this_d$status])))
                      this_d <- this_d[as.character(this_d$l4b_class) %in% shared, ]
                      this_d$l4b_class <- droplevels(factor(this_d$l4b_class))
                      if (nlevels(this_d$l4b_class) >= 2) f <- update(f, ~ . + strata(l4b_class))
                    }
                    if (nrow(d_wocat) > 2) {
                      model <- glm(f, data=this_d)
                      dists <- match_on(model, data=this_d)  #***NEED PROPENSITY SCORE BEFORE HERE
                    } else {
                      # Use Mahalanobis distance if there aren't enough points to run a glm
                      dists <- match_on(f, data=this_d)
                    }
                    #potentially drop caliper line; will cut down dists matrix but not the speed issue
                    # dists <- caliper(dists, 2)
                    # If the controls are too far from the treatments (due to the caliper) 
                    # then the matching may fail. Can test for this by seeing if subdim 
                    # runs successfully
                    subdim_works <- tryCatch(is.data.frame(subdim(dists)),
                                             error=function(e)return(FALSE))
                    if (subdim_works) {
                      m <- fullmatch(dists, min.controls=1, max.controls=1, data=this_d)
                      prematch_d <- this_d
                      this_d$matched <- m
                      this_d <- this_d[matched(m), ]
                    } else {
                      this_d <- data.frame()
                    }
                    # Need to handle the possibility that there were no matches for this 
                    # treatment, meaning this_d will be an empty data.frame
                    if (nrow(this_d) == 0) {   #if matching w/ ecoreg return no results, match again w/o ecoreg and check matching results
                      this_d <- df
                      d_wocat <- filter(this_d, status)
                      this_d <- filter(this_d,
                                       land_cover %in% unique(d_wocat$land_cover),
                                       wwfbiom %in% unique(d_wocat$wwfbiom))
                      
                      this_d$land_cover <- droplevels(this_d$land_cover)
                      #this_d$wwfbiom <- droplevels(this_d$wwfbiom)
                      if (is.factor(this_d$wwfbiom)) this_d$wwfbiom <- droplevels(this_d$wwfbiom)  # gedi-pa-change patch (see above)
                      f <- ps_formula
                      if (nlevels(this_d$land_cover) >= 2) {
                        f <- update(f, ~ . + strata(land_cover))
                      } else {
                        f <- update(f, ~ . - land_cover)
                      }
                      if (nlevels(this_d$wwfbiom) >= 2) {
                        f <- update(f, ~ . + strata(wwfbiom))
                      } else {
                        f <- update(f, ~ . - wwfbiom)
                      }
                      # gedi-pa-change patch (l4b_classes only): exact matching within L4B classes.
                      # Keep classes that have both PA and control cells (others cannot be matched
                      # exactly), then stratify on the class like land_cover.
                      if (isTRUE(get0("l4b_strata", ifnotfound = FALSE))) {
                        shared <- intersect(unique(as.character(this_d$l4b_class[this_d$status])),
                                            unique(as.character(this_d$l4b_class[!this_d$status])))
                        this_d <- this_d[as.character(this_d$l4b_class) %in% shared, ]
                        this_d$l4b_class <- droplevels(factor(this_d$l4b_class))
                        if (nlevels(this_d$l4b_class) >= 2) f <- update(f, ~ . + strata(l4b_class))
                      }
                      if (nrow(d_wocat) > 2) {
                        model <- glm(f, data=this_d)
                        dists <- match_on(model, data=this_d)
                      } else {
                        dists <- match_on(f, data=this_d)
                      }
                      subdim_works <- tryCatch(is.data.frame(subdim(dists)),
                                               error=function(e)return(FALSE))
                      if (subdim_works) {
                        m <- fullmatch(dists, min.controls=1, max.controls=1, data=this_d)
                        prematch_d <- this_d
                        this_d$matched <- m
                        this_d <- this_d[matched(m), ]
                      } else {
                        this_d <- data.frame()
                      }
                      if(nrow(this_d)==0){
                        log_mes <- paste0(pid,"-Matching without wwfecoreg:Failed\n")
                        dir.create(paste0(f.path,"matching_log"))
                        cat(log_mes,file=paste0(f.path,"matching_log/",iso3,"_pa_",pid,"_used_covar_log.txt"),append=TRUE)
                        return(NULL)
                      } else {
                        log_mes <- paste0(pid,"-Matching without wwfecoreg:Succeeded\n")
                        dir.create(paste0(f.path,"matching_log"))
                        cat(log_mes,file=paste0(f.path,"matching_log/",iso3,"_pa_",pid,"_used_covar_log.txt"),append=TRUE)
                        match_results <- list("match_obj" = m, "df" = this_d, "func"=f, "prematch_d"=prematch_d)
                        return(match_results)
                      }
                    } else {
                      log_mes <- paste0(pid,"-Matching with wwfecoreg:Succeed\n")
                      match_results <- list("match_obj" = m, "df" = this_d, "func"=f, "prematch_d"=prematch_d)
                      dir.create(paste0(f.path,"matching_log"))
                      cat(log_mes,file=paste0(f.path,"matching_log/",iso3,"_pa_",pid,"_used_covar_log_wk.txt"),append=TRUE)
                      return(match_results)
                    }
                  }

  # Same object as the original foreach: one land-cover value -> foreach returns the single
  # result unchanged; two or more -> foreach_rbind of identical copies (ret[1, ] is the same
  # for any number of copies); all NULL -> NULL.
  one <- match_once()
  ret <- if (n_lc >= 2) foreach_rbind(one, one) else if (n_lc == 1) one else NULL
  return(ret)
}

#--------------------------------------------------------------------------------
propensity_filter <- function(pa_df, d_control_local) {
    
  pa_df <- pa_df[complete.cases(pa_df), ]  #filter away non-complete cases w/ NA in control set
  d <- dplyr::bind_rows(d_control_local, pa_df)
  ## bring in matching algorithm from STEP5 here to loop through each PA in d_PAs
  # filter controls based on propensity scores 
  d_all <- dplyr::select(d, lat, lon, UID, status, land_cover, wwfbiom, wwfecoreg, 
                         elevation, slope, mean_temp, max_temp, min_temp, prec, 
                         d2road, d2city, tt2city, popden, popcnt, gedi_l4b)
  # gedi-pa-change patch: carry l4b_class through when present (l4b_classes runs only)
  if ("l4b_class" %in% names(d)) d_all$l4b_class <- d$l4b_class

  
  d_all$status <- ifelse(d_all$status==TRUE,1,0)
  
  #calculate the propensity scores & filter out controls not overlapping w/ treatment propensity scores
  ps <- glm(ps_formula, data = d_all)  # gedi-pa-change patch: was the hard-coded formula
  # boxplot(ps) # check the distribution of propensity scores for treatment and controls
  # filter out the controls with propensity scores outside of the overlapping region
  d_all$propensity_score <- fitted(ps)
  d_sep <- d_all %>% dplyr::group_by(status)
  d_sep_range <- d_all %>% dplyr::group_by(status) %>% 
    dplyr::summarise(propmin=min(propensity_score), promax=max(propensity_score)) 
  # cat(iso3, "Filtering the control sites by overlaping with the treatment PS\n")
  d_filtered <- d_sep %>% 
    filter(status==1 | between(propensity_score, d_sep_range$propmin[2], d_sep_range$promax[2])) %>% 
    ungroup() 
  
  d_filtered$status <- ifelse(d_filtered$status==1,TRUE,FALSE)
  
  return(d_filtered)
}

#--------------------------------------------------------------------------------
matched2ras <- function(matched_df){
  cat(iso3,"converting the matched points to a raster stack for extraction\n")
  
#  matched_pts <- SpatialPointsDataFrame(coords=matched_df[, c("lon","lat")],
#                                        proj4string=CRS("EPSG:4326"),
#                                        data=matched_df) %>% spTransform(., CRS("EPSG:6933"))

  matched_pts <- vect(matched_df, geom = c("lon","lat"), crs = "EPSG:4326")
  matched_pts <- project(matched_pts, "EPSG:6933")
  
  # Ensure fields are in appropriate formats
  matched_pts$UID <- as.integer(matched_pts$UID)
  matched_pts$pa_id <- as.integer(matched_pts$pa_id)
  matched_pts$status <- as.logical(matched_pts$status)
  matched_pts$wwfbiom <- as.numeric(matched_pts$wwfbiom)
  matched_pts$wwfecoreg <- as.numeric(matched_pts$wwfecoreg)
 
  # Define the raster extent for cropping
  buffer_ext <- ext(buffer(matched_pts, 10000))
  r <- crop(MCD12Q1, buffer_ext)
  continent <- crop(world_region, buffer_ext)
  
  # Assign names to the raster layers
  names(r) <- "pft"
  names(continent) <- "region"
  
  # List of fields to rasterize
  fields <- c("status", "pa_id", "UID", "land_cover", "wwfbiom", "wwfecoreg")
  rasters <- list()
  
  # Rasterize each field
  for (field in fields) {
    r_field <- rasterize(matched_pts, r, field = field)
    rasters[[field]] <- r_field
  }
  
  # Combine rasters into a SpatRaster stack
  rasters <- c(rasters, list(pft = r, region = continent))
  matched_ras <- rast(rasters)
  
  return(matched_ras)
}

#--------------------------------------------------------------------------------
convertFactor <- function(matched0, exgedi){
  exgedi$pft <- as.character(exgedi$pft)
  
  exgedi$pft <- factor(exgedi$pft, levels=sequence(6),
                       labels = c("ENT",
                                  "EBT",
                                  "ENT",
                                  "DBT",
                                  "GS",
                                  "GS"))
  exgedi$region <- as.character(exgedi$region)
  exgedi$region <- factor(exgedi$region, levels=c(1:7),
                          labels = c("Eu",
                                     "As",
                                     "Au",
                                     "Af",
                                     "As",
                                     "SA",
                                     "US"))
  
  exgedi$stratum <- paste(exgedi$pft, exgedi$region,sep="_")
    
  return(exgedi)
}

#--------------------------------------------------------------------------------
subdfExport <- function(filtered_df){
  #export invidual pa results
  spt2 <- split(filtered_df, filtered_df$pa_id)
  
  dfl <- lapply(names(spt2), function(x){
    
    if(dim(spt2[[x]])[1]>0){
      control_sub <- spt2[[x]][spt2[[x]]$status==0,]
      treat_sub <-  spt2[[x]][spt2[[x]]$status==1,]
      ncontrol <- nrow(control_sub)
      ntreat <- nrow(treat_sub)
      if (ntreat-ncontrol > 0){
        newtreatid <- sample(ntreat, ncontrol)
        newtreat <- treat_sub[newtreatid,]
        spt2_new <- rbind(newtreat, control_sub)
      } else if (ntreat - ncontrol< 0){
        newcontrolid <- sample(ncontrol, ntreat)
        newcontrol <- control_sub[newcontrolid,]
        spt2_new <- rbind(newcontrol, treat_sub)
      } else if (ntreat-ncontrol==0){
        spt2_new <- spt2[[x]]
      } else if (ntreat==0 || ncontrol==0){
        spt2_new=NA
      }
##      biom <- spt2_new$wwfbiom %>% unique() %>% as.character() %>% gsub('\\b(\\pL)\\pL{4,}|.','\\U\\1',.,perl = TRUE)
##      if(length(biom)>1){
##        biom <- paste(c(biom), collapse="&")
##      }
      # print(biom)
##      write.csv(spt2_new, file=paste(f.path,"WDPA_GEDI_extract/",iso3,"_wk",gediwk,"/",iso3,"_PA_",unique(spt2_new$pa_id),"_",biom,".csv", sep=""))
        write.csv(spt2_new, file=paste(f.path,"WDPA_GEDI_extract/",iso3,"_wk",gediwk,"/",iso3,"_PA_",unique(spt2_new$pa_id),".csv", sep=""))

      return(spt2_new)
    }
  })
  
  total_df <- do.call("rbind", dfl) 
  cat("Exported individual PAs results for ", iso3, "\n")
  
  return(total_df)
}

getmode <- function(v,na.rm) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#--------------------------------------------------------------------------------
extract_gedi <- function(matched, mras, iso3) {

  # Define the GEDI data directory (new structure)
  path2gedi <- "s3://maap-ops-workspace/shared/ameliah/gedi-test/brazil_tiles/data/"
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  cat("▶ Starting GEDI extraction for", iso3, "\n")

  # Convert matched points to terra vector
  matched_pts <- terra::vect(matched, geom = c("lon", "lat"), crs = "EPSG:4326")
  matched_pts_prj <- terra::project(matched_pts, "EPSG:6933")

  # Compute bounding box of the matched points 
  bbox <- st_bbox(matched_pts)
  lat_range <- seq(floor(bbox["ymin"]), ceiling(bbox["ymax"]))
  lon_range <- seq(floor(bbox["xmin"]), ceiling(bbox["xmax"]))

  # Generate all tile IDs overlapping the matched points
  generate_tile_id <- function(lat, lon) {
    lat_prefix <- ifelse(lat >= 0, "N", "S")
    lon_prefix <- ifelse(lon >= 0, "E", "W")
    lat_val <- sprintf("%02d", abs(floor(lat)))
    lon_val <- sprintf("%03d", abs(floor(lon)))
    paste0(lat_prefix, lat_val, "_", lon_prefix, lon_val)
  }

  tile_ids <- expand.grid(lat = lat_range, lon = lon_range) %>%
    dplyr::mutate(tile_id = generate_tile_id(lat, lon)) %>%
    dplyr::pull(tile_id)

  cat("Found", length(tile_ids), "potential GEDI tiles overlapping PA region\n")

  # Define years to process
  years <- c(2019, 2020, 2022, 2023)

  # Define variables to extract
  cols2read <- c("shot_number", "year", "lat_lowestmode", "lon_lowestmode", "rh_98", "agbd", "agbd_se")
    
  # Function to safely read parquet from S3 with DuckDB
  safe_read_duck <- function(p, cols = cols2read, bbox, buffer_km) {
    tryCatch({
      # Convert buffer from km to degrees (~0.009° per km)
      deg_buffer <- buffer_km * 0.009
      bbox_expanded <- bbox
      bbox_expanded["xmin"] <- bbox["xmin"] - deg_buffer
      bbox_expanded["xmax"] <- bbox["xmax"] + deg_buffer
      bbox_expanded["ymin"] <- bbox["ymin"] - deg_buffer
      bbox_expanded["ymax"] <- bbox["ymax"] + deg_buffer

      local_path <- s3::s3_get(p, force = TRUE)

      # Cast shot_number to string to preserve all digits
      col_exprs <- ifelse(cols == "shot_number", "CAST(shot_number AS VARCHAR) AS shot_number", cols)
      col_list <- paste(col_exprs, collapse = ", ")
        
      q <- sprintf("SELECT %s
                    FROM read_parquet('%s')
                    WHERE lat_lowestmode BETWEEN %f AND %f
                    AND lon_lowestmode BETWEEN %f AND %f",
                   col_list, local_path,
                   bbox_expanded['ymin'], bbox_expanded['ymax'],
                   bbox_expanded['xmin'], bbox_expanded['xmax'])

      df <- DBI::dbGetQuery(con, q)
      if (nrow(df) == 0) return(NULL)
      df
    }, error = function(e) {
      cat("  → Error reading", p, ":", e$message, "\n")
      return(NULL)
    })
  }

  # Read all relevant GEDI tiles for the matched region
  gedi_list <- list()
  for (tile_id in tile_ids) {
    for (year in years) {
      file_path <- paste0(path2gedi, "tile_id=", tile_id, "/year=", year, "/data_0.parquet")
      cat("Reading GEDI tile:", tile_id, "year:", year, "\n")
      df <- safe_read_duck(file_path, bbox = bbox, buffer_km = 1.5)
      if (!is.null(df)) {
        df$tile_id <- tile_id
        df$year <- year
        gedi_list[[paste(tile_id, year, sep = "_")]] <- df
      }
    }
  }

  if (length(gedi_list) == 0) {
    cat("⚠ No GEDI data found for matched region in", iso3, "\n")
    DBI::dbDisconnect(con, shutdown = TRUE)
    return(NULL)
  }

  # Combine all tiles
  gedi_all <- dplyr::bind_rows(gedi_list)
  gedi_all$shot_number <- as.character(gedi_all$shot_number)
  cat("Combined GEDI dataset:", nrow(gedi_all), "shots\n")

  # Convert to SpatVector for extraction
  gedi_pts <- terra::vect(gedi_all, geom = c("lon_lowestmode", "lat_lowestmode"), crs = "EPSG:4326")
  gedi_pts_prj <- terra::project(gedi_pts, "EPSG:6933")

  # Extract rasterized matched PA information
  matched_gedi <- terra::extract(mras, gedi_pts_prj, df = TRUE)
  matched_gedi_metrics <- cbind(matched_gedi, terra::as.data.frame(gedi_pts_prj))

  coords <- as.data.frame(terra::crds(gedi_pts))
  names(coords) <- c("lon_lowestmode", "lat_lowestmode")
  matched_gedi_metrics <- cbind(matched_gedi_metrics, coords)

  # Filter to only those overlapping matched areas
  matched_gedi_metrics <- matched_gedi_metrics %>%
    dplyr::filter(!is.na(status)) %>%
    dplyr::select(pa_id, status, land_cover, pft, region,
                  wwfbiom, wwfecoreg, lat_lowestmode, lon_lowestmode,
                  rh_98, agbd, agbd_se, year, shot_number)
  
  cat("✅ Extracted", nrow(matched_gedi_metrics), "GEDI shots overlapping matched grid cells\n")

  DBI::dbDisconnect(con, shutdown = TRUE)
  return(matched_gedi_metrics)
}

#--------------------------------------------------------------------------------
SplitRas <- function(raster,ppside){
  h        <- ceiling(ncol(raster)/ppside)
  v        <- ceiling(nrow(raster)/ppside)
  agg      <- aggregate(raster,fact=c(h,v))
  agg[]    <- 1:ncell(agg)
  agg_poly <- rasterToPolygons(agg)
  names(agg_poly) <- "polis"
  r_list <- list()
  for(i in 1:ncell(agg)){
    e1          <- extent(agg_poly[agg_poly$polis==i,])
    r_list[[i]] <- crop(raster,e1)
  }
  return(r_list)
}

#--------------------------------------------------------------------------------
rasExtract2020 <- function(l4_sp){
  cat(iso3,"converting the matched csv to a raster stack for extraction\n")
  tif2020 <- c("pop_cnt_2020", "pop_den_2020", "lc2019", "tt2cities_2015", 
               "wc_prec_2010-2018", "wc_tavg_2010-2018", "wc_tmax_2010-2018",
               "wc_tmin_2010-2018", "dem", "slope", "d2roads", "dcities") #"wwf_biomes","wwf_ecoreg",
  for (t in 1:length(tif2020)){
    # print(tif2020[t])
    covar2020 <- raster(paste(f.path, "WDPA_input_vars_iso3_v2/",iso3,"/",tif2020[t],".tif", sep=""))
    ras_ex <- raster::extract(covar2020, l4_sp@coords, method="simple", factors=F)
    nm <- names(covar2020)
    l4_sp <- cbind(l4_sp, ras_ex)
    names(l4_sp)[t+6] <- tif2020[t]
  }
  return(l4_sp)
}

#---------------------------------------------------------------------------------
generate_tile_id <- function(lat, lon) {
  tile_lat <- ceiling(lat)  # always use ceiling for north edge
  tile_lon <- floor(lon)    # always use floor for west edge
  lat_prefix <- ifelse(tile_lat >= 0, "N", "S")
  lon_prefix <- ifelse(tile_lon >= 0, "E", "W")
  lat_val <- sprintf("%02d", abs(tile_lat))
  lon_val <- sprintf("%03d", abs(tile_lon))
  paste0(lat_prefix, lat_val, "_", lon_prefix, lon_val)
}

                                               