# =============================================================================
# Movement behavior and habitat selection to guide conservation of an
# endangered Neotropical deer
#
# Alarcón, P.A.E., Crego, R.D., Galetto, E., Aprile, G.,
# Morales, J.M., Pauli, J.N. & Donadio, E.
#
# -----------------------------------------------------------------------------
# Script:  04_prep_huemul_steps.R
# Purpose: Shared data-prep helper for the iSSF. Defines
#          `prepare_huemul_steps()`, which builds the per-step design data
#          (used + available, with all covariates extracted at endpoints)
#          used by the iSSF estimation, simulation and figure scripts.
#          Caches the result in ./issf/data/huemulL_inla.rda so subsequent
#          calls are cheap.
# =============================================================================

suppressPackageStartupMessages({
  library(amt)
  library(terra)
  library(sf)
  library(suncalc)
  library(lwgeom)
  library(dplyr)
  library(purrr)
  library(lubridate)
  library(readr)
})

# Canonical data folder for this script. All cache + train/test files
# read or written by `prepare_huemul_steps()` go through this constant.
PREP_DATA_DIR <- "./old/4. hissf/1. collar_data"

prepare_huemul_steps <- function(
  gps_file        = "./0. cleaningGPSData/1. gpsFiles/1. clean_gpsdata.rda",
  stack_file      = "./0. covariates/_revised/2. stacks/stack_static.tif",
  legacy_stack    = "./0. covariates/_older/stack_issf.tif",
  snow_file       = "./0. covariates/_revised/2. stacks/snow_daily.tif",
  snow_dates_csv  = "./0. covariates/_revised/2. stacks/snow_dates.csv",
  cache_file      = file.path(PREP_DATA_DIR, "huemulL_stan.rda"),
  use_cache       = TRUE,
  target_crs      = 32719,
  n_avail_steps   = 10,
  trail_decay_m   = 300,
  hiker_traffic   = c(
    "01" = 1.00, "02" = 0.95, "03" = 0.80, "04" = 0.40,
    "05" = 0.10, "06" = 0.02, "07" = 0.02, "08" = 0.05,
    "09" = 0.20, "10" = 0.55, "11" = 0.85, "12" = 1.00),
  verbose         = TRUE
) {

  msg <- function(...) if (verbose) message(...)

  # 0. Cache hit (with staleness check) ---------------------------------------
  if (use_cache && file.exists(cache_file)) {
    cache_mtime  <- file.info(cache_file)$mtime
    sources      <- c(stack_file, snow_file, snow_dates_csv,
                      legacy_stack, gps_file)
    sources      <- sources[file.exists(sources)]
    src_mtimes   <- file.info(sources)$mtime
    newer_inputs <- sources[src_mtimes > cache_mtime]

    if (length(newer_inputs)) {
      msg("[prepare_huemul_steps] cache is STALE (newer inputs detected):")
      for (p in newer_inputs) msg(sprintf("  newer: %s", p))
      msg("  -> rebuilding from sources.")
    } else {
      msg(sprintf("[prepare_huemul_steps] cache hit: %s", cache_file))
      e <- new.env()
      load(cache_file, envir = e)
      if (exists("huemulL", envir = e))
        return(e$huemulL)
      msg("  cache file present but does not contain `huemulL`; rebuilding.")
    }
  }

  # 1. Resolve which covariate stack to use ------------------------------------
  if (!file.exists(gps_file))
    stop("GPS file not found: ", gps_file)

  use_new_stack <- file.exists(stack_file)
  use_legacy    <- !use_new_stack && file.exists(legacy_stack)
  if (!use_new_stack && !use_legacy)
    stop("No covariate stack found. Expected one of:\n  - ", stack_file,
         "\n  - ", legacy_stack,
         "\nRun 4. covariates_unified.R first.")

  active_stack <- if (use_new_stack) stack_file else legacy_stack
  msg(sprintf("[prepare_huemul_steps] using stack: %s", active_stack))

  # Track which optional layers are actually available
  stack_peek <- terra::rast(active_stack)
  layers_avail <- names(stack_peek)
  has_canopy   <- "canopy"          %in% layers_avail
  has_trails   <- "dist_to_trails"  %in% layers_avail
  has_snow_files <- file.exists(snow_file) && file.exists(snow_dates_csv)

  msg(sprintf("  layers in stack: %s", paste(layers_avail, collapse = ", ")))
  msg(sprintf("  canopy:         %s", if (has_canopy)     "OK" else "missing -> NA"))
  msg(sprintf("  dist_to_trails: %s", if (has_trails)     "OK" else "missing -> NA"))
  msg(sprintf("  snow stack:     %s", if (has_snow_files) "OK" else "missing -> NA"))

  # 2. GPS data ----------------------------------------------------------------
  msg("[prepare_huemul_steps] loading GPS data ...")
  e <- new.env(); load(gps_file, envir = e); gpsdata <- e$gpsdata
  gpsdata <- subset(gpsdata, !is.na(gpsdata$x_utm))
  gpsdata$DateTime <- as.POSIXct(ymd_hms(gpsdata$local_datetime))
  gpsdata <- gpsdata[order(gpsdata$animal_id, gpsdata$DateTime), ]
  gpsdata$DateTime <- round(gpsdata$DateTime, units = "hours")

  gpsdataL <- split(gpsdata, gpsdata$animal_id)
  for (h in seq_along(gpsdataL))
    gpsdataL[[h]] <- gpsdataL[[h]][!duplicated(gpsdataL[[h]]$DateTime), ]
  gpsdata <- do.call(rbind.data.frame, gpsdataL)

  # 3. 70/30 train/test split (same as estimation script) ----------------------
  split_gps <- gpsdata %>%
    arrange(animal_id, DateTime) %>%
    group_by(animal_id) %>%
    group_split() %>%
    map_dfr(function(df) {
      n <- nrow(df); s <- floor(0.7 * n)
      df$set <- c(rep("train", s), rep("test", n - s))
      df
    })
  gps_train <- filter(split_gps, set == "train")

  dir.create(PREP_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
  save(gps_train,
       file = file.path(PREP_DATA_DIR, "train_gpsdata_hissf_stan.rda"))
  gps_test <- filter(split_gps, set == "test")
  save(gps_test,
       file = file.path(PREP_DATA_DIR, "test_gpsdata_hissf_stan.rda"))

  # 4. Build steps and simulate available steps --------------------------------
  msg("[prepare_huemul_steps] building steps and simulating available ...")
  gpsdataL <- split(gps_train, gps_train$animal_id)
  nid <- length(gpsdataL)

  huemulL <- vector("list", nid)
  for (h in seq_len(nid)) {
    d <- data.frame(x_ = gpsdataL[[h]]$x_utm,
                    y_ = gpsdataL[[h]]$y_utm,
                    t_ = gpsdataL[[h]]$DateTime)
    trk <- amt::make_track(d, .x = x_, .y = y_, .t = t_,
                           crs = target_crs, check_duplicates = TRUE)
    trk <- amt::track_resample(trk, rate = hours(2),
                               tolerance = minutes(15), start = 1)
    trk <- amt::steps_by_burst(trk)
    trk <- amt::random_steps(trk, n_control = n_avail_steps)
    bad <- which(is.na(trk$x2_))
    if (length(bad))
      trk <- subset(trk, !(step_id_ %in% unique(trk$step_id_[bad])))
    huemulL[[h]] <- trk
  }

  # 5. Static covariates -------------------------------------------------------
  msg("[prepare_huemul_steps] extracting static covariates ...")
  stack <- rast(active_stack)

  # Validate that the layers we need are actually in this stack.
  required <- c("forest", "slope", "dist_to_build", "river_zones")
  miss_req <- setdiff(required, layers_avail)
  if (length(miss_req))
    stop("Stack '", active_stack, "' is missing required layer(s):\n  - ",
         paste(miss_req, collapse = "\n  - "),
         "\nLayers actually available in this stack:\n  - ",
         paste(layers_avail, collapse = "\n  - "),
         "\nFix: regenerate the stack with 4. covariates_unified.R, or rename\n",
         "the layer(s) inside the existing stack with `names(rast(...)) <- ...`.")

  # stream_zones is OPTIONAL — it only exists once the user has run the
  # updated `2. rebuild_river_zones.R` against streamZonification_revised.shp.
  # When missing we still emit a streamCrossing column (all NA) so downstream
  # `complete.cases()` filters in 02 stay coherent.
  has_streams <- "stream_zones" %in% layers_avail
  if (!has_streams) {
    msg("[prepare_huemul_steps] stream_zones layer NOT in stack — ",
        "streamCrossing will be NA. Re-run `2. rebuild_river_zones.R` ",
        "after placing streamZonification_revised.shp in 0. shapefiles/",
        "riversLakes to populate it.")
  }

  for (h in seq_len(nid)) {
    huemulL[[h]] <- amt::extract_covariates(huemulL[[h]],
                                            stack[["forest"]], where = "both")
    huemulL[[h]] <- amt::extract_covariates(huemulL[[h]],
                                            stack[["slope"]],  where = "both")
    huemulL[[h]] <- amt::extract_covariates(huemulL[[h]],
                                            stack[["dist_to_build"]],
                                            where = "both")
    huemulL[[h]] <- amt::extract_covariates(huemulL[[h]],
                                            stack[["river_zones"]],
                                            where = "both")

    if (has_streams) {
      huemulL[[h]] <- amt::extract_covariates(huemulL[[h]],
                                              stack[["stream_zones"]],
                                              where = "both")
    } else {
      huemulL[[h]]$stream_zones_start <- NA_integer_
      huemulL[[h]]$stream_zones_end   <- NA_integer_
    }

    if (has_canopy) {
      huemulL[[h]] <- amt::extract_covariates(huemulL[[h]],
                                              stack[["canopy"]], where = "both")
    } else {
      huemulL[[h]]$canopy_start <- NA_real_
      huemulL[[h]]$canopy_end   <- NA_real_
    }

    if (has_trails) {
      huemulL[[h]] <- amt::extract_covariates(huemulL[[h]],
                                              stack[["dist_to_trails"]],
                                              where = "both")
    } else {
      huemulL[[h]]$dist_to_trails_start <- NA_real_
      huemulL[[h]]$dist_to_trails_end   <- NA_real_
    }

    # MAIN RIVER crossing — from `river_zones`.
    huemulL[[h]]$riverCrossing <- ifelse(
      huemulL[[h]]$river_zones_start == huemulL[[h]]$river_zones_end, 0, 1)

    # MOUNTAIN STREAM crossing — from `stream_zones`.
    huemulL[[h]]$streamCrossing <- ifelse(
      huemulL[[h]]$stream_zones_start == huemulL[[h]]$stream_zones_end, 0, 1)

    huemulL[[h]]$animal_id <- unique(gpsdataL[[h]]$animal_id)
  }

  # 6. Snow (time-varying) -----------------------------------------------------
  if (has_snow_files) {
    msg("[prepare_huemul_steps] extracting time-varying snow ...")
    snow_stack <- rast(snow_file)
    snow_dates <- read_csv(snow_dates_csv, show_col_types = FALSE) %>%
      mutate(date = as.Date(date))

    snow_band_for_date <- function(d) {
      idx <- match(d, snow_dates$date)
      if (is.na(idx)) idx <- which.min(abs(snow_dates$date - d))
      idx
    }

    extract_snow <- function(steps_df) {
      out <- rep(NA_real_, nrow(steps_df))
      step_dates <- as.Date(steps_df$t2_)
      for (d in unique(step_dates)) {
        rows <- which(step_dates == d)
        band <- snow_band_for_date(as.Date(d))
        if (is.na(band)) next
        pts  <- cbind(steps_df$x2_[rows], steps_df$y2_[rows])
        vals <- terra::extract(snow_stack[[band]], pts)[, 1]
        # Mask sentinel < 0 (cloud not filled by GEE) as NA.
        vals <- ifelse(!is.na(vals) & vals < 0, NA_real_, vals)
        out[rows] <- vals
      }
      out
    }

    for (h in seq_len(nid))
      huemulL[[h]]$snow_end <- extract_snow(huemulL[[h]])
  } else {
    msg("[prepare_huemul_steps] snow files missing -> snow_end = NA")
    for (h in seq_len(nid))
      huemulL[[h]]$snow_end <- NA_real_
  }

  # 7. Hiker exposure (T(month) * exp(-d_trail / sigma)) -----------------------
  if (has_trails) {
    msg("[prepare_huemul_steps] computing hiker exposure ...")
    for (h in seq_len(nid)) {
      d_trail <- huemulL[[h]]$dist_to_trails_end
      m       <- sprintf("%02d", lubridate::month(huemulL[[h]]$t2_))
      huemulL[[h]]$hikerExp_end <- hiker_traffic[m] *
                                     exp(-d_trail / trail_decay_m)
    }
  } else {
    msg("[prepare_huemul_steps] trail layer missing -> hikerExp_end = NA")
    for (h in seq_len(nid))
      huemulL[[h]]$hikerExp_end <- NA_real_
  }

  # 8. Day/night flag ----------------------------------------------------------
  msg("[prepare_huemul_steps] flagging day/night ...")
  for (h in seq_len(nid)) {
    utm_df <- huemulL[[h]][, c("x2_", "y2_")]
    utm_v  <- vect(utm_df, geom = c("x2_", "y2_"), crs = "EPSG:32719")
    ll_v   <- terra::project(utm_v, "EPSG:4326")
    ll     <- as.data.frame(crds(ll_v))
    d <- data.frame(lon = ll$x, lat = ll$y,
                    date = as.Date(huemulL[[h]]$t2_))
    times <- getSunlightTimes(data = d,
                              tz   = "America/Argentina/Buenos_Aires",
                              keep = c("sunrise", "sunset"))
    huemulL[[h]]$day_night_end <- ifelse(
      huemulL[[h]]$t2_ > times$sunset |
      huemulL[[h]]$t2_ <= times$sunrise, 0, 1)
  }

  # 9. Replace zero step lengths -----------------------------------------------
  for (h in seq_len(nid))
    huemulL[[h]]$sl_ <- ifelse(huemulL[[h]]$sl_ == 0, 0.01, huemulL[[h]]$sl_)

  # 10. Cache ------------------------------------------------------------------
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  save(huemulL, file = cache_file)
  msg(sprintf("[prepare_huemul_steps] saved cache: %s", cache_file))

  invisible(huemulL)
}
