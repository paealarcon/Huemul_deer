# =============================================================================
# Movement behavior and habitat selection to guide conservation of an
# endangered Neotropical deer
#
# Alarcón, P.A.E., Crego, R.D., Galetto, E., Aprile, G.,
# Morales, J.M., Pauli, J.N. & Donadio, E.
#
# -----------------------------------------------------------------------------
# Script:  03c_BoyceIndex.R
# Purpose: Cross-validate BOTH the forest- and steppe-intercept hRSF fits on
#          the held-out 30% validation set using the continuous Boyce Index
#          (Hirzel et al. 2006). The Boyce index ranges from -1 to +1:
#          values near +1 indicate predictions consistent with the presence
#          distribution, 0 indicates no better than random, and negative
#          values indicate counter-predictions.
# =============================================================================

# 0. Libraries ----
library(sf)
library(terra)

source("./hrsf/scripts/00_helper_functions.R")

# 1. Data ----
load("./hrsf/outputs/hrsf_forest_intercept.rda")               # hrsf_forest
results_F <- readRDS("./hrsf/outputs/hrsf_forest_intercept_results.rds")

load("./hrsf/outputs/hrsf_steppe_intercept.rda")               # hrsf_steppe
results_S <- readRDS("./hrsf/outputs/hrsf_steppe_intercept_results.rds")

load("./hrsf/data/UsedAv_INLA.rda")                    # training data
load("./hrsf/data/gpsdata_for_hrsf_test.rda")          # test_data

# 2. Habitat layers & scaling parameters ----
lc_list <- list(
  Forest   = st_read("./data/layers/habitatTypes/forest.shp",     quiet = TRUE),
  Steppe   = st_read("./data/layers/habitatTypes/steppes.shp",    quiet = TRUE),
  Meadow   = st_read("./data/layers/habitatTypes/meadow.shp",     quiet = TRUE),
  water    = st_read("./data/layers/habitatTypes/water.shp",      quiet = TRUE),
  snow     = st_read("./data/layers/habitatTypes/snow.shp",       quiet = TRUE),
  Andean   = st_read("./data/layers/habitatTypes/highAndean.shp", quiet = TRUE),
  bareSoil = st_read("./data/layers/habitatTypes/bareSoil.shp",   quiet = TRUE)
)
for (cat_name in names(lc_list)) lc_list[[cat_name]]$landCat <- cat_name
lc <- do.call(rbind, lc_list)
lc <- st_transform(lc, crs = 32719)

stack <- rast("./data/layers/stack_static.tif")

# Scaling parameters from training data (must match model fitting)
elev_mean  <- mean(UsedAv$elevation, na.rm = TRUE)
elev_sd    <- sd(UsedAv$elevation,   na.rm = TRUE)
build_mean <- mean(UsedAv$build,     na.rm = TRUE)
build_sd   <- sd(UsedAv$build,       na.rm = TRUE)

# 3. Annotate test locations with covariates ----
annotate_test_points <- function(test_df) {
  pts_sf <- st_as_sf(test_df, coords = c("x_utm", "y_utm"), crs = 32719)

  elev_raw  <- terra::extract(stack$elevation, pts_sf)$elevation
  elev_scl  <- (elev_raw - elev_mean) / elev_sd
  elev_scl2 <- elev_scl^2

  build_raw <- terra::extract(stack$dist_to_build, pts_sf)$dist_to_build
  build_scl <- (build_raw - build_mean) / build_sd

  hab <- st_join(pts_sf, lc)$landCat

  # Include dummies for BOTH reference schemes (forest-intercept uses steppe/
  # meadow/andean/soil; steppe-intercept uses forest/meadow/andean/soil).
  data.frame(
    animal_id = test_df$animal_id,
    forest    = ifelse(hab == "Forest",   1, 0),
    steppe    = ifelse(hab == "Steppe",   1, 0),
    meadow    = ifelse(hab == "Meadow",   1, 0),
    andean    = ifelse(hab == "Andean",   1, 0),
    soil      = ifelse(hab == "bareSoil", 1, 0),
    elev_scl  = elev_scl,
    elev_scl2 = elev_scl2,
    build_scl = build_scl,
    season    = ifelse(test_df$season == "cold", 0, 1),
    stringsAsFactors = FALSE
  )
}

test_annotated <- annotate_test_points(test_data)
test_annotated <- test_annotated[complete.cases(test_annotated), ]
cat("Test points annotated:", nrow(test_annotated), "\n")

# 4. Forest-intercept RSF scores ----
# Population-level relative selection strength: RSF score = exp(eta), using
# the fixed-effect (population-mean) coefficients. The free per-individual
# intercept is omitted, since it cancels out in relative comparisons.
compute_rsf_score <- function(results, dat) {
  fe     <- results$fixed_effects
  coefs  <- setNames(fe$post_mean, fe$variable)
  needed <- names(coefs)

  season_pairs <- list(
    forest_season = "forest",
    steppe_season = "steppe",
    meadow_season = "meadow",
    andean_season = "andean",
    elev_season   = "elev_scl",
    build_season  = "build_scl"
  )
  for (int_col in names(season_pairs)) {
    base_col <- season_pairs[[int_col]]
    if (int_col %in% needed && !int_col %in% names(dat)) {
      dat[[int_col]] <- dat[[base_col]] * dat$season
    }
  }

  eta <- rep(0, nrow(dat))
  for (cov in needed) eta <- eta + coefs[cov] * dat[[cov]]

  exp(eta)
}

# Training data: rebuild the same covariates used to fit the model, so the
# available (background) points define the expected suitability distribution.
# Include BOTH forest and steppe dummies so the same `train` object works
# for both intercept parameterizations.
train           <- UsedAv
train$forest    <- ifelse(train$habitatCat == "Forest",   1, 0)
train$steppe    <- ifelse(train$habitatCat == "Steppe",   1, 0)
train$meadow    <- ifelse(train$habitatCat == "Meadow",   1, 0)
train$andean    <- ifelse(train$habitatCat == "Andean",   1, 0)
train$soil      <- ifelse(train$habitatCat == "bareSoil", 1, 0)
train$elev_scl  <- (train$elevation - elev_mean) / elev_sd
train$elev_scl2 <- train$elev_scl^2
train$build_scl <- (train$build - build_mean) / build_sd
train$season    <- ifelse(train$season == "cold", 0, 1)

avail_idx <- which(as.integer(train$use) == 0)

# 5. Out-of-sample Boyce Index — FOREST intercept ----
rsf_train_F <- compute_rsf_score(results_F, train)
rsf_test_F  <- compute_rsf_score(results_F, test_annotated)
boyce_F <- calc_boyce_index(
  fit_values = rsf_train_F[avail_idx],
  obs_values = rsf_test_F
)
cat("Forest-intercept Boyce Index (out-of-sample):",
    round(boyce_F$cor, 3), "\n")

# 6. Out-of-sample Boyce Index — STEPPE intercept ----
rsf_train_S <- compute_rsf_score(results_S, train)
rsf_test_S  <- compute_rsf_score(results_S, test_annotated)
boyce_S <- calc_boyce_index(
  fit_values = rsf_train_S[avail_idx],
  obs_values = rsf_test_S
)
cat("Steppe-intercept Boyce Index (out-of-sample):",
    round(boyce_S$cor, 3), "\n")

# 7. Save output ----
boyce_results <- list(
  forest_intercept = list(
    model      = "hRSF_Forest_intercept",
    level      = "population",
    evaluation = "out-of-sample (30% held-out validation set)",
    boyce      = boyce_F$cor,
    curve      = boyce_F
  ),
  steppe_intercept = list(
    model      = "hRSF_Steppe_intercept",
    level      = "population",
    evaluation = "out-of-sample (30% held-out validation set)",
    boyce      = boyce_S$cor,
    curve      = boyce_S
  )
)

dir.create("./hrsf/outputs", showWarnings = FALSE, recursive = TRUE)
saveRDS(boyce_results, "./hrsf/outputs/boyce_index_results.rds")
cat("Saved: boyce_index_results.rds\n")
cat("\n=== Summary ===\n")
cat(sprintf("  Forest-intercept Boyce Index: %.3f\n", boyce_F$cor))
cat(sprintf("  Steppe-intercept Boyce Index: %.3f\n", boyce_S$cor))
