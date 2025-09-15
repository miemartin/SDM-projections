# R/run_sdm.R
# Full SDM pipeline wrapper using Wallace functions (v2) + raster/terra utilities
# Purpose: produce current SDM, project to CMIP6 scenarios, compute p10 threshold, gain/loss, summaries
# Author: Generated for Mia Martin
# Date: 2025-09-15

run_sdm <- function(species_name,
                    occ_csv,
                    envs_dir,
                    env_layers = c('bio02.tif','bio05.tif','bio06.tif','bio08.tif','bio14.tif','bio15.tif','bio16.tif','bio18.tif'),
                    shp_transfer,
                    cmip6_dir = "data/cmip6/future",
                    out_dir = "results",
                    thin_km = 10,
                    bg_pts = 10000,
                    kfolds = 25,
                    rm_values = seq(0.5, 6, 0.5),
                    fc_options = c("L","Q","P","H","LQ","LQHP"),
                    best_rm = 1,
                    best_fc = "LQHP",
                    maxent_jar = "maxent.jar",
                    seed = 1234,
                    parallel = FALSE,
                    numCores = 1) {
  # Load utils
  script_dir <- file.path(getwd(), "R")
  source(file.path(script_dir, "utils.R"))
  pkgs <- c("wallace", "spocc", "spThin", "dismo", "ENMeval",
            "geodata", "leaflet", "raster", "terra", "dplyr", "sf", "rgdal", "rgeos", "stringr")
  pkg_check(pkgs)
  set.seed(seed)
  
  # prepare output directories
  out_dir <- file.path(out_dir, species_name)
  safe_dir(out_dir)
  safe_dir(file.path(out_dir, "rasters"))
  safe_dir(file.path(out_dir, "tables"))
  safe_dir(file.path(out_dir, "figures"))
  
  log_msg("Starting SDM for:", species_name)
  log_msg("Reading occurrences from:", occ_csv)
  
  # 1) load occurrences via wallace helper (keeps pipeline similar to your existing code)
  userOccs <- tryCatch({
    occs_userOccs(txtPath = occ_csv, txtName = basename(occ_csv), txtSep = ",", txtDec = ".")
  }, error = function(e) stop("occs_userOccs failed: ", e$message))
  
  # in your previous script you expected e.g. userOccs_Ae$Aedes_aegypti$cleaned
  # We'll allow both "Aedes_aegypti" or species_name
  key_names <- names(userOccs)
  # try direct lookup then fuzzy
  found_key <- if (species_name %in% key_names) species_name else grep(species_name, key_names, value = TRUE)[1]
  if (is.null(found_key) || length(found_key) == 0) {
    stop("Species not found in occ file. Available keys: ", paste(key_names, collapse = ", "))
  }
  occs <- userOccs[[found_key]]$cleaned
  if (!all(c("longitude", "latitude") %in% names(occs))) {
    stop("Occurrences must have 'longitude' and 'latitude' columns")
  }
  log_msg("Occurrences loaded. n = ", nrow(occs))
  write.csv(occs, file.path(out_dir, "occurrences_cleaned.csv"), row.names = FALSE)
  
  # 2) Load environmental rasters (current)
  env_paths <- file.path(envs_dir, env_layers)
  if (!all(file.exists(env_paths))) {
    missing <- env_paths[!file.exists(env_paths)]
    stop("Missing env layers: ", paste(missing, collapse = ", "))
  }
  envs <- envs_userEnvs(rasPath = env_paths, rasName = basename(env_paths), doBrick = FALSE)
  log_msg("Loaded current environmental layers.")
  
  # 3) Extract env values at occs, remove duplicates/NA as you did
  occs_xy <- occs[, c("longitude", "latitude")]
  occs_vals <- as.data.frame(raster::extract(envs, occs_xy, cellnumbers = TRUE))
  keep_rows <- !duplicated(occs_vals[,1]) & !rowSums(is.na(occs_vals[,-1]))
  occs_clean <- occs[keep_rows, ]
  occs_vals <- occs_vals[keep_rows, -1]
  occs_clean <- cbind(occs_clean, occs_vals)
  log_msg("After duplicate/NA removal: n = ", nrow(occs_clean))
  write.csv(occs_clean, file.path(out_dir, "occurrences_filtered.csv"), row.names = FALSE)
  
  # 4) Thinning
  occs_thin <- tryCatch({
    poccs_thinOccs(occs = occs_clean, thinDist = thin_km)
  }, error = function(e) {
    log_msg("spThin failed, using unthinned occurrences: ", e$message)
    occs_clean
  })
  write.csv(occs_thin, file.path(out_dir, "occurrences_thinned.csv"), row.names = FALSE)
  log_msg("After thinning: n = ", nrow(occs_thin))
  
  # 5) Background sampling
  bgExt <- penvs_bgExtent(occs = occs_thin, bgSel = "point buffers", bgBuf = 2)
  bgMask <- penvs_bgMask(occs = occs_thin, envs = envs, bgExt = bgExt)
  bgSample <- penvs_bgSample(occs = occs_thin, bgMask = bgMask, bgPtsNum = bg_pts)
  log_msg("Background sampled: n = ", nrow(bgSample))
  
  # prepare bgEnvsVals similar to your script (so model_maxent accepts it)
  bgEnvsVals <- cbind(
    scientific_name = paste0("bg_", species_name),
    bgSample,
    occID = NA, year = NA, institution_code = NA, country = NA,
    state_province = NA, locality = NA, elevation = NA, record_type = NA,
    as.data.frame(raster::extract(bgMask, bgSample))
  )
  write.csv(bgEnvsVals, file.path(out_dir, "tables", "bg_envs_sample.csv"), row.names = FALSE)
  
  # save background extent polygon as shapefile
  try({
    bgExt_sf <- sf::st_as_sf(bgExt)
    sf::st_write(bgExt_sf, file.path(out_dir, paste0("M_", species_name, ".shp")), delete_layer = TRUE)
  }, silent = TRUE)
  
  # 6) Partition occurrences
  groups <- part_partitionOccs(occs = occs_thin, bg = bgSample, method = "rand", kfolds = kfolds)
  log_msg("Partitioned occurrences into ", kfolds, " folds.")
  
  # 7) Model tuning: run model_maxent (like your script)
  log_msg("Running Maxent tuning (this can take several minutes).")
  model_tune <- model_maxent(
    occs = occs_thin,
    bg = bgEnvsVals,
    user.grp = groups,
    bgMsk = bgMask,
    rms = rm_values,
    rmsStep = 0.5,
    fcs = fc_options,
    clampSel = TRUE,
    algMaxent = maxent_jar,
    parallel = parallel
  )
  
  # Save evaluation table
  eval_tbl <- model_tune@results
  write.csv(eval_tbl, file.path(out_dir, "tables", paste0("models_", species_name, ".csv")), row.names = FALSE)
  log_msg("Saved model tuning results.")
  
  # 8) Fit final model with chosen parameters (best_rm & best_fc)
  log_msg("Fitting final Maxent model: fc=", best_fc, "_rm=", best_rm)
  final_model <- model_maxent(
    occs = occs_thin,
    bg = bgEnvsVals,
    user.grp = groups,
    bgMsk = bgMask,
    rms = c(best_rm, best_rm),
    rmsStep = 1,
    fcs = c(best_fc),
    clampSel = TRUE,
    algMaxent = maxent_jar,
    parallel = parallel,
    numCores = numCores
  )
  # save permutation importance if available
  perm_imp <- tryCatch(final_model@variable.importance, error = function(e) NULL)
  if (!is.null(perm_imp)) write.csv(perm_imp, file.path(out_dir, "tables", paste0("perm_importance_", species_name, ".csv")), row.names = FALSE)
  
  # 9) Transfer current model to transfer polygon (shp_transfer)
  log_msg("Reading transfer polygon:", shp_transfer)
  xfer_userExt <- sf::st_read(shp_transfer, quiet = TRUE)
  xferAreaEnvs <- envs  # using current envs for transfer area
  xfer_area_out <- xfer_area(
    evalOut = final_model,
    curModel = paste0("fc.", best_fc, "_rm.", best_rm),
    envs = xferAreaEnvs,
    outputType = "cloglog",
    alg = maxent_jar,
    clamp = TRUE,
    xfExt = xfer_userExt
  )
  # store current suitability
  current_suit <- xfer_area_out$xferArea
  safe_write_raster(current_suit, file.path(out_dir, "rasters", "current_suitability.tif"), overwrite = TRUE)
  log_msg("Saved current suitability raster.")
  
  # 10) Compute p10 (10th percentile training presence threshold)
  train_preds <- raster::extract(current_suit, occs_thin[, c("longitude", "latitude")])
  train_preds <- train_preds[!is.na(train_preds)]
  threshold_10tp <- as.numeric(quantile(train_preds, 0.10, na.rm = TRUE))
  saveRDS(threshold_10tp, file.path(out_dir, "tables", "threshold_10tp.rds"))
  write.csv(data.frame(threshold_10tp = threshold_10tp), file.path(out_dir, "tables", "threshold_10tp.csv"), row.names = FALSE)
  log_msg("Computed p10 threshold:", threshold_10tp)
  
  # 11) Prepare CMIP6 projections: search cmip6_dir for tif stacks and loop
  # Expecting files named with GCM/SSP and period, or a single multilayer .tif per GCM/SSP
  cmip6_files <- list.files(cmip6_dir, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
  if (length(cmip6_files) == 0) {
    log_msg("No CMIP6 files found in cmip6_dir:", cmip6_dir)
  } else {
    log_msg("Found", length(cmip6_files), "cmip6 files.")
  }
  
  cmip6_projections <- list()
  for (f in cmip6_files) {
    # name scheme: extract filename without extension
    fname <- basename(f)
    scenario_name <- tools::file_path_sans_ext(fname)
    log_msg("Processing CMIP6 file:", fname)
    # attempt to load as SpatRaster
    cmip6_stack <- tryCatch(terra::rast(f), error = function(e) { log_msg("Failed loading:", f, ":", e$message); NULL })
    if (is.null(cmip6_stack)) next
    # name layers "bio01..bio19" if needed
    nl <- terra::nlyr(cmip6_stack)
    if (!all(grepl("^bio", names(cmip6_stack)))) {
      names(cmip6_stack) <- sprintf("bio%02d", 1:nl)
    }
    # select same layers as training: use names of bgMask or envs
    trained_names <- names(bgMask) # bgMask was created earlier and matched envs
    # if trained_names not suitable, fallback to env_layers names
    trained_names <- intersect(names(cmip6_stack), c(names(envs), sprintf("bio%02d", 1:nl)))
    if (length(trained_names) < length(env_layers)) {
      # Attempt to pick by number if same order presumed
      trained_names <- names(cmip6_stack)[which(names(cmip6_stack) %in% sprintf("bio%02d", c(2,5,6,8,14,15,16,18)))]
    }
    if (length(trained_names) < length(env_layers)) {
      log_msg("Warning: could not fully match training layers for", fname, " - using available layers.")
    }
    xfer_envs <- cmip6_stack[[env_layers %||% names(cmip6_stack)[1:min(nl, length(env_layers))]]]
    # Transfer model in time
    xfer_time_out <- tryCatch({
      xfer_time(
        evalOut = final_model,
        curModel = paste0("fc.", best_fc, "_rm.", best_rm),
        envs = xfer_envs,
        xfExt = xfer_userExt,
        alg = maxent_jar,
        outputType = "cloglog",
        clamp = TRUE
      )
    }, error = function(e) {
      log_msg("xfer_time failed for", scenario_name, ":", e$message)
      NULL
    })
    if (is.null(xfer_time_out)) next
    cmip6_projections[[scenario_name]] <- list(envs = xfer_envs, sdm = xfer_time_out$xferTime)
    # Save scenario raster
    safe_write_raster(xfer_time_out$xferTime, file.path(out_dir, "rasters", paste0("SDM_", scenario_name, ".tif")), overwrite = TRUE)
    log_msg("Saved SDM:", scenario_name)
  } # end cmip6 loop
  
  # Save cmip6_projections object for later use
  saveRDS(cmip6_projections, file.path(out_dir, "tables", "cmip6_projections.rds"))
  
  # 12) Average rasters by SSP and period if filenames encode that info (best-effort)
  # We'll compute mean across all cmip6 projections grouped by common SSP_period substring (extract SSPxxx or yyyy-yyyy)
  names_all <- names(cmip6_projections)
  if (length(names_all) > 0) {
    df <- data.frame(scenario = names_all,
                     ssp = stringr::str_extract(names_all, "SSP\\d{3}|ssp\\d{3}|SSP\\d+"),
                     period = stringr::str_extract(names_all, "\\d{4}-\\d{4}"),
                     stringsAsFactors = FALSE)
    df$ssp[is.na(df$ssp)] <- "unknown"
    df$period[is.na(df$period)] <- "unknown"
    avg_list <- list()
    for (s in unique(df$ssp)) {
      for (p in unique(df$period)) {
        sel <- df$scenario[df$ssp == s & df$period == p]
        if (length(sel) == 0) next
        rasters_to_avg <- lapply(sel, function(nm) cmip6_projections[[nm]]$sdm)
        # convert to terra rast stack and compute mean
        rast_stack <- terra::rast(rasters_to_avg)
        avg_r <- tryCatch(terra::app(rast_stack, fun = mean, na.rm = TRUE), error = function(e) NULL)
        if (!is.null(avg_r)) {
          avg_name <- paste0(ifelse(is.na(s), "SSPunk", s), "_", ifelse(is.na(p), "period_unk", p))
          avg_list[[avg_name]] <- avg_r
          safe_write_raster(avg_r, file.path(out_dir, "rasters", paste0("meanSDM_", avg_name, ".tif")), overwrite = TRUE)
        }
      }
    }
    saveRDS(avg_list, file.path(out_dir, "tables", "averaged_rasters.rds"))
    log_msg("Saved averaged rasters.")
  } else {
    log_msg("No cmip6 projections available to average.")
  }
  
  # 13) Calculate p10-based binary maps, gain/loss, area stats (uses terra to compute areas)
  # Load current and averaged rasters (if present)
  current_tr <- terra::rast(file.path(out_dir, "rasters", "current_suitability.tif"))
  threshold_val <- threshold_10tp
  # mask current to polygon (if vector usable)
  xfer_vect <- tryCatch(terra::vect(xfer_userExt), error = function(e) NULL)
  if (!is.null(xfer_vect)) current_tr <- terra::mask(current_tr, xfer_vect)
  current_bin_tr <- current_tr >= threshold_val
  # per-cell area in km2
  cell_area_km2 <- terra::cellSize(current_tr, unit = "km")
  current_area_km2 <- terra::global(cell_area_km2 * (current_bin_tr), "sum", na.rm = TRUE)[[1]]
  log_msg("Current suitable area (km2):", round(current_area_km2, 2))
  
  # loop future scenarios (cmip6_projections)
  scenario_summary <- list()
  for (nm in names(cmip6_projections)) {
    fut_tr <- tryCatch(terra::rast(cmip6_projections[[nm]]$sdm), error = function(e) NULL)
    if (is.null(fut_tr)) next
    # resample to current grid
    fut_tr_res <- terra::resample(fut_tr, current_tr, method = "bilinear")
    # mask
    if (!is.null(xfer_vect)) fut_tr_res <- terra::mask(fut_tr_res, xfer_vect)
    fut_bin_tr <- fut_tr_res >= threshold_val
    # areas
    fut_area_km2 <- terra::global(cell_area_km2 * (fut_bin_tr), "sum", na.rm = TRUE)[[1]]
    gain_km2 <- terra::global(cell_area_km2 * ((fut_bin_tr == 1) & (current_bin_tr == 0)), "sum", na.rm = TRUE)[[1]]
    loss_km2 <- terra::global(cell_area_km2 * ((fut_bin_tr == 0) & (current_bin_tr == 1)), "sum", na.rm = TRUE)[[1]]
    net_km2 <- fut_area_km2 - current_area_km2
    scenario_summary[[nm]] <- data.frame(scenario = nm,
                                         current_km2 = current_area_km2,
                                         future_km2 = fut_area_km2,
                                         gain_km2 = gain_km2,
                                         loss_km2 = loss_km2,
                                         net_change_km2 = net_km2,
                                         stringsAsFactors = FALSE)
    # save binary rasters for this scenario
    safe_write_raster(fut_bin_tr, file.path(out_dir, "rasters", paste0("bin_p10_", nm, ".tif")), overwrite = TRUE)
    # save gain/loss rasters (as numeric 0/1)
    gain_r <- (fut_bin_tr == 1) & (current_bin_tr == 0)
    loss_r <- (current_bin_tr == 1) & (fut_bin_tr == 0)
    safe_write_raster(gain_r, file.path(out_dir, "rasters", paste0("gain_p10_", nm, ".tif")), overwrite = TRUE)
    safe_write_raster(loss_r, file.path(out_dir, "rasters", paste0("loss_p10_", nm, ".tif")), overwrite = TRUE)
  }
  # bind scenario summary and save
  if (length(scenario_summary) > 0) {
    scen_df <- do.call(rbind, scenario_summary)
    write.csv(scen_df, file.path(out_dir, "tables", "scenario_area_p10_summary.csv"), row.names = FALSE)
  } else {
    log_msg("No scenario summaries to write.")
  }
  
  # 14) Averaged rasters summary (if present)
  avg_list_path <- file.path(out_dir, "tables", "averaged_rasters.rds")
  if (file.exists(avg_list_path)) {
    avg_list <- readRDS(avg_list_path)
    avg_rows <- list()
    for (nm in names(avg_list)) {
      fut_tr <- avg_list[[nm]]
      fut_tr_res <- terra::resample(fut_tr, current_tr, method = "bilinear")
      if (!is.null(xfer_vect)) fut_tr_res <- terra::mask(fut_tr_res, xfer_vect)
      fut_bin_tr <- fut_tr_res >= threshold_val
      fut_area_km2 <- terra::global(cell_area_km2 * (fut_bin_tr), "sum", na.rm = TRUE)[[1]]
      gain_km2 <- terra::global(cell_area_km2 * ((fut_bin_tr == 1) & (current_bin_tr == 0)), "sum", na.rm = TRUE)[[1]]
      loss_km2 <- terra::global(cell_area_km2 * ((fut_bin_tr == 0) & (current_bin_tr == 1)), "sum", na.rm = TRUE)[[1]]
      net_km2 <- fut_area_km2 - current_area_km2
      avg_rows[[nm]] <- data.frame(ssp_period = nm,
                                   current_km2 = current_area_km2,
                                   future_km2 = fut_area_km2,
                                   gain_km2 = gain_km2,
                                   loss_km2 = loss_km2,
                                   net_change_km2 = net_km2,
                                   stringsAsFactors = FALSE)
      # save binary average rasters
      safe_write_raster(fut_bin_tr, file.path(out_dir, "rasters", paste0("mean_bin_p10_", nm, ".tif")), overwrite = TRUE)
    }
    avg_df <- do.call(rbind, avg_rows)
    write.csv(avg_df, file.path(out_dir, "tables", "averaged_area_p10_summary.csv"), row.names = FALSE)
    log_msg("Saved averaged area summaries.")
  }
  
  # 15) Save session info
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  log_msg("SDM run completed for:", species_name)
  
  return(invisible(list(
    final_model = final_model,
    current = current_tr,
    threshold_10tp = threshold_10tp,
    cmip6_projections = cmip6_projections,
    out_dir = out_dir
  )))
}
