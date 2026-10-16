# =============================================================================
# Movement behavior and habitat selection to guide conservation of an
# endangered Neotropical deer
#
# Alarcón, P.A.E., Crego, R.D., Galetto, E., Aprile, G.,
# Morales, J.M., Pauli, J.N. & Donadio, E.
#
# -----------------------------------------------------------------------------
# Script:  04b_hissfSimulation.R
# Purpose: Forward-simulate 2000-step trajectories from each individual's
#          fitted iSSF (amt::simulate_path, Signer et al. 2024) and compare
#          them visually to the observed GPS paths.
# =============================================================================

# 0. Libraries ----
suppressPackageStartupMessages({
  library(INLA)
  library(MASS)            # fitdistr
  library(circular)        # mle.vonmises
  library(terra)
  library(suncalc)          # day/night covariate
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(lubridate)
  library(amt)              # make_issf_model, simulate_path
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# 1. Control panel ----
N_STEPS             <- 2000             # Signer et al. (2024) case study
N_CANDIDATES        <- 1000             # amt's n.control
SL_TRUNC_QUANTILE   <- 0.99             # truncate the empirical gamma tail
START_FROM          <- "first_obs"      # "first_obs" | "centroid"
DAYNIGHT_SOURCE     <- "observed_times" # "observed_times" | "fixed_day" |
                                         # "fixed_night" | "alternate"
STEP_INTERVAL_HOURS <- 2
SEED                <- 20260527
SAVE_PATHS_RDA      <- "./issf/outputs/sim_pathsL.rda"
SAVE_FIG_DIR        <- "./issf/figures/simulation_noDN_withElev"

# River crossings are excluded from the fitted model (non-identifiable — see
# 04a). To reflect that huemul essentially never cross the De las Vueltas
# river without banning it outright, candidates landing on the opposite side
# of a huemul's home zone get a fixed, low (not zero) relative selection
# strength during simulation.
BARRIER_LOG_RSS <- log(0.01)   # ~1% relative odds for a river crossing

# Clamp on the home-range pull: the fitted xy2 coefficient must stay
# negative for the inward pull to make sense; individuals whose posterior
# mean comes out weaker than this ceiling are clamped to it.
XY2_CLAMP_MAX <- -0.05

set.seed(SEED)

# 2. Data ----
cat("=== Loading INLA iSSF fit (noDN + elevation) ===\n")
load("./issf/outputs/hissf_inla.rda")                      # hissf_inla
results <- readRDS("./issf/outputs/hissf_inla_results.rds")
fit     <- hissf_inla
mapping <- results$mapping
cat(sprintf("  mapping entries: %s\n", paste(names(mapping), collapse = ", ")))

load("./issf/data/huemulL_inla.rda")             # huemulL
huemul <- as.data.frame(do.call(rbind.data.frame, huemulL))
nid    <- length(huemulL)
cat(sprintf("  huemulL: %d animals, %d pooled rows\n", nid, nrow(huemul)))

if (!"rugged_end" %in% names(huemul) && !"ruggedness_end" %in% names(huemul)) {
  static_stack <- rast("./data/layers/stack_static.tif")
  pts <- vect(huemul[, c("x2_", "y2_")], geom = c("x2_", "y2_"), crs = "EPSG:32719")
  huemul$rugged_end <- terra::extract(static_stack[["ruggedness"]], pts)[, 2]
} else if ("ruggedness_end" %in% names(huemul) && !"rugged_end" %in% names(huemul)) {
  huemul$rugged_end <- huemul$ruggedness_end
}

# Hold-out test set, used only to overlay observed test-period locations on
# the validation figure alongside the training-period locations.
test_setL <- NULL
test_path <- "./issf/data/test_gpsdata_hissf_inla.rda"
if (file.exists(test_path)) {
  e <- new.env(); load(test_path, envir = e)
  if (exists("gps_test", envir = e))
    test_setL <- split(e$gps_test, e$gps_test$animal_id)
}

# 3. Per-individual posterior coefficients (beta_h = mu + u_h) ----
cat("\n=== Extracting per-individual posterior coefficients ===\n")

get_term_coefs <- function(fit, mapping, nid) {
  out <- matrix(0, nrow = nid, ncol = length(mapping),
                dimnames = list(NULL, names(mapping)))
  for (tname in names(mapping)) {
    m <- mapping[[tname]]
    mu <- 0
    if (!is.null(m$fixed) && m$fixed %in% rownames(fit$summary.fixed))
      mu <- fit$summary.fixed[m$fixed, "mean"]
    u <- rep(0, nid)
    if (!is.null(m$random) && !is.null(fit$summary.random[[m$random]])) {
      re <- fit$summary.random[[m$random]]
      u  <- re$mean[match(seq_len(nid), re$ID)]
      u[is.na(u)] <- 0
    }
    out[, tname] <- mu + u
  }
  out
}
beta_mat <- get_term_coefs(fit, mapping, nid)
cat(sprintf("  beta_mat: %d individuals x %d terms\n",
            nrow(beta_mat), ncol(beta_mat)))

if ("Xy2" %in% colnames(beta_mat)) {
  to_clamp <- which(beta_mat[, "Xy2"] > XY2_CLAMP_MAX)
  if (length(to_clamp)) {
    cat(sprintf("  clamping %d huemul beta_Xy2 -> %.3f\n",
                length(to_clamp), XY2_CLAMP_MAX))
    beta_mat[to_clamp, "Xy2"] <- XY2_CLAMP_MAX
  }
}

# 4. Home-range and covariate scaling parameters (from the 04a fit) ----
hr_scaling <- results$hr_scaling
hr_centroids <- hr_scaling$centroids
names(hr_centroids) <- c("animal_id", "x_centre", "y_centre")
x_mean <- hr_scaling$x_mean;     x_sd   <- hr_scaling$x_sd
y_mean <- hr_scaling$y_mean;     y_sd   <- hr_scaling$y_sd
xy2_mean <- hr_scaling$xy2_mean; xy2_sd <- hr_scaling$xy2_sd

scl_env <- results$scl_env

mean_sd <- function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
sl_for_scale <- ifelse(huemul$sl_ == 0, 0.1, huemul$sl_)
scl_step    <- mean_sd(sl_for_scale)
scl_logstep <- mean_sd(log(sl_for_scale))

# 5. Static raster stack + crossing-zone rasters ----
cat("\n=== Loading static raster stack ===\n")
static_stack <- "./data/layers/stack_static.tif"
stk          <- rast(static_stack)
cat("  stack layers:", paste(names(stk), collapse = ", "), "\n")

# Elevation lives in the hRSF static stack (not in the hissf one). Load the
# layer and append it to `stk` so the simulation `fun_h` can extract it like
# any other covariate from `map`.
if (!"elevation" %in% names(stk)) {
  stk_hrsf <- rast("./data/layers/stack_static.tif")
  stopifnot("elevation" %in% names(stk_hrsf))
  # Resample to the hissf stack grid so coord lookups align
  elev_lyr <- terra::resample(stk_hrsf[["elevation"]], stk, method = "bilinear")
  stk      <- c(stk, elev_lyr)
  cat("  added elevation layer -> stk layers:",
      paste(names(stk), collapse = ", "), "\n")
}

river_zones  <- if ("river_zones"  %in% names(stk)) stk[["river_zones"]]  else NULL
stream_zones <- if ("stream_zones" %in% names(stk)) stk[["stream_zones"]] else NULL

# 6. Day/night source, per individual ----
dn_per_step <- vector("list", nid)
for (h in seq_len(nid)) {
  dfh <- huemulL[[h]]
  t_h <- if (!is.null(dfh$t_)) dfh$t_ else dfh$t2_
  ts0 <- as.POSIXct(min(t_h, na.rm = TRUE))
  ts_seq <- ts0 + (seq_len(N_STEPS) - 1L) * STEP_INTERVAL_HOURS * 3600

  if (DAYNIGHT_SOURCE == "fixed_day")        dn <- rep(1L, N_STEPS)
  else if (DAYNIGHT_SOURCE == "fixed_night") dn <- rep(0L, N_STEPS)
  else if (DAYNIGHT_SOURCE == "alternate")   dn <- rep(c(1L, 0L), length.out = N_STEPS)
  else {
    ctr <- hr_centroids[hr_centroids$animal_id == unique(dfh$animal_id), ]
    pt_ll <- terra::project(vect(ctr[, c("x_centre", "y_centre")],
                                  geom = c("x_centre", "y_centre"),
                                  crs = "EPSG:32719"), "EPSG:4326")
    ll <- as.data.frame(crds(pt_ll))
    d  <- data.frame(date = as.Date(ts_seq), lon = ll$x, lat = ll$y)
    sl <- getSunlightTimes(data = d, tz = "America/Argentina/Buenos_Aires",
                            keep = c("sunrise", "sunset"))
    dn <- as.integer(ts_seq <= sl$sunset & ts_seq > sl$sunrise)
  }
  dn_per_step[[h]] <- list(times = ts_seq, dn = dn)
}

# 7. Mapping from beta_mat column names to amt coefficient names ----
# amt expects movement betas named "sl_"/"log_sl_"/"cos_ta_" by default, but
# 04a fit beta_Step/beta_Logstep on standardized step_scl/logStep_scl (and
# beta_CosTurn on raw cos(ta_)). To stay on the same scale used at fit time,
# the `fun` closure below computes step_scl/logStep_scl/cosTurn itself and
# these names are matched to the corresponding beta_mat columns.
to_amt_name <- function(term) {
  switch(term,
    "Step"    = "step_scl",   "Logstep" = "logStep_scl", "CosTurn" = "cosTurn",
    "X"       = "x_scl",      "Y"       = "y_scl",       "Xy2"     = "xy2_scl",
    "Slope"   = "slope_scl",  "Slope2"  = "slope2",
    "Rugged"  = "rugged_scl", "Rugged2" = "rugged2",
    "Canopy"  = "canopy_scl", "Canopy2" = "canopy2",
    "Build"   = "build_scl",  "Build2"  = "build2",
    "Elev"    = "elev_scl",   "Elev2"   = "elev2",
    "Stream"  = "stream",
    tolower(term)
  )
}

# 8. Custom `fun` factory ----
# `fun(xy, map)` receives amt's candidate data.frame (start x_/y_, end
# x2_/y2_, plus t2_) and returns it with one extra column per coef name in
# coefs_h. amt then weights candidates by exp(sum_k coefs_h[k] * df[, k]).
make_fun_h <- function(h) {
  id_h <- unique(huemulL[[h]]$animal_id)
  ctr  <- hr_centroids[hr_centroids$animal_id == id_h, ]
  cx   <- ctr$x_centre
  cy   <- ctr$y_centre

  # Home zone for the soft river barrier (section 1): determined once per
  # huemul from its centroid.
  home_zone_h <- NA_integer_
  if (!is.null(river_zones)) {
    z_centre <- terra::extract(river_zones, cbind(cx, cy))[, 1]
    if (!is.na(z_centre)) home_zone_h <- as.integer(z_centre)
  }

  function(xy, map) {
    n_pts <- nrow(xy)

    pick_xy <- function(x_cand, y_cand) {
      for (i in seq_along(x_cand))
        if (!is.null(xy[[x_cand[i]]]) && !is.null(xy[[y_cand[i]]]))
          return(list(x = xy[[x_cand[i]]], y = xy[[y_cand[i]]]))
      NULL
    }
    end_xy <- pick_xy(c("x2_", "x_end_", "x_", "x"), c("y2_", "y_end_", "y_", "y"))
    if (is.null(end_xy)) stop("Could not find endpoint columns in `xy`.")
    x2 <- end_xy$x
    y2 <- end_xy$y

    # Movement covariates on the scale 04a fit on (step_scl/logStep_scl are
    # standardized; cosTurn is raw cos(ta_)).
    if ("sl_" %in% names(xy)) {
      xy$step_scl    <- (xy$sl_ - scl_step["mean"]) / scl_step["sd"]
      xy$logStep_scl <- (log(pmax(xy$sl_, 0.01)) - scl_logstep["mean"]) /
                          scl_logstep["sd"]
    }
    if ("ta_" %in% names(xy)) xy$cosTurn <- cos(xy$ta_)

    # Environmental covariates, extracted and scaled the same way as 04a.
    ext <- function(name) {
      if (!(name %in% names(map))) return(rep(NA_real_, n_pts))
      terra::extract(map[[name]], cbind(x2, y2))[, 1]
    }
    scl <- function(raw, nm) {
      m <- scl_env[[nm]]; if (is.null(m)) return(raw)
      (raw - m["mean"]) / m["sd"]
    }
    xy$slope_scl  <- scl(ext("slope"),         "slope")
    xy$rugged_scl <- scl(ext("ruggedness"),    "rugged")
    xy$canopy_scl <- scl(ext("canopy"),        "canopy")
    xy$build_scl  <- scl(ext("dist_to_build"), "build")
    xy$elev_scl   <- scl(ext("elevation"),     "elev")
    xy$slope2  <- xy$slope_scl^2
    xy$rugged2 <- xy$rugged_scl^2
    xy$canopy2 <- xy$canopy_scl^2
    xy$build2  <- xy$build_scl^2
    xy$elev2   <- xy$elev_scl^2

    # NOTE: day/night interactions were dropped from the fitted model
    # (noDN + elev). No *_DN columns are needed here.

    # Home-range pull: per-individual centred x, y, and their sum of squares
    # (Signer et al. 2024), on the same scale 04a fit on.
    x_c <- x2 - cx
    y_c <- y2 - cy
    xy$x_scl   <- (x_c - x_mean) / x_sd
    xy$y_scl   <- (y_c - y_mean) / y_sd
    xy$xy2_scl <- (x_c^2 + y_c^2 - xy2_mean) / xy2_sd

    # Stream crossing (fitted covariate).
    has_start <- all(c("x1_", "y1_") %in% names(xy))
    if (!is.null(stream_zones) && has_start) {
      z_start <- terra::extract(stream_zones, cbind(xy[["x1_"]], xy[["y1_"]]))[, 1]
      z_end   <- terra::extract(stream_zones, cbind(x2, y2))[, 1]
      xy$stream <- as.integer(!is.na(z_start) & !is.na(z_end) & z_start != z_end)
    } else {
      xy$stream <- 0
    }

    # Soft river barrier (not a fitted covariate — see section 1): a fixed
    # low RSS multiplier for candidates landing on the opposite side of the
    # huemul's home zone.
    if (!is.null(river_zones) && !is.na(home_zone_h)) {
      z_end <- terra::extract(river_zones, cbind(x2, y2))[, 1]
      xy$barrier <- as.integer(!is.na(z_end) & z_end != home_zone_h)
    } else {
      xy$barrier <- 0
    }

    xy
  }
}

# 9. Per-huemul simulation loop ----
cat(sprintf("\n=== Simulating (N_STEPS=%d, N_CANDIDATES=%d, %d huemul) ===\n",
            N_STEPS, N_CANDIDATES, nid))

simL <- vector("list", nid)
names(simL) <- sapply(huemulL, function(d) unique(d$animal_id))

for (h in seq_len(nid)) {
  id_h <- names(simL)[h]
  cat(sprintf("\n--- %s (%d / %d) ---\n", id_h, h, nid))

  used <- huemulL[[h]][huemulL[[h]]$case_ == TRUE | huemulL[[h]]$case_ == "TRUE", ]
  sl_h <- used$sl_[is.finite(used$sl_) & used$sl_ > 0]
  ta_h <- used$ta_[is.finite(used$ta_)]
  if (length(sl_h) < 50 || length(ta_h) < 50) {
    cat("  too few used steps -> skipping\n"); next
  }

  fit_g <- tryCatch(MASS::fitdistr(sl_h, "gamma"), error = function(e) NULL)
  if (is.null(fit_g)) {
    m <- mean(sl_h); v <- var(sl_h)
    shape <- m^2 / v; rate <- m / v
  } else {
    shape <- as.numeric(fit_g$estimate["shape"])
    rate  <- as.numeric(fit_g$estimate["rate"])
  }
  scale_ <- 1 / rate

  vm <- tryCatch(circular::mle.vonmises(circular(ta_h)), error = function(e) NULL)
  kappa <- if (!is.null(vm)) as.numeric(vm$kappa) else 1

  sl_cap <- as.numeric(quantile(sl_h, SL_TRUNC_QUANTILE, na.rm = TRUE))

  sl_distr <- amt::make_gamma_distr(shape = shape, scale = scale_)
  ta_distr <- amt::make_vonmises_distr(kappa = kappa)

  cat(sprintf("  empirical: shape=%.2f scale=%.1f kappa=%.3f sl_cap=%.0f m\n",
              shape, scale_, kappa, sl_cap))

  beta_h    <- beta_mat[h, , drop = TRUE]
  amt_names <- vapply(names(beta_h), to_amt_name, character(1))
  coefs_h   <- setNames(as.numeric(beta_h), amt_names)
  coefs_h   <- coefs_h[!is.na(names(coefs_h)) & nzchar(names(coefs_h))]
  coefs_h   <- c(coefs_h, barrier = BARRIER_LOG_RSS)
  cat(sprintf("  coefs: %d term(s) passed to amt (+ soft river barrier)\n",
              length(coefs_h) - 1L))

  m_issf <- amt::make_issf_model(coefs = coefs_h, sl = sl_distr, ta = ta_distr)

  if (START_FROM == "centroid") {
    ctr <- hr_centroids[hr_centroids$animal_id == id_h, ]
    x0 <- ctr$x_centre; y0 <- ctr$y_centre
  } else {
    x0 <- huemulL[[h]]$x2_[1]
    y0 <- huemulL[[h]]$y2_[1]
  }
  t0 <- as.POSIXct(min(huemulL[[h]]$t_, na.rm = TRUE))

  start_h <- amt::make_start(x = c(x0, y0), ta_ = 0, time = t0,
                              dt = lubridate::hours(STEP_INTERVAL_HOURS))

  fun_h <- make_fun_h(h)
  k_h <- tryCatch(
    amt::redistribution_kernel(
      m_issf, start = start_h, map = stk, fun = fun_h,
      max.dist = sl_cap, n.control = N_CANDIDATES,
      landscape = "continuous", tolerance.outside = 0.05
    ),
    error = function(e) {
      cat(sprintf("  redistribution_kernel failed: %s\n", conditionMessage(e)))
      NULL
    })
  if (is.null(k_h)) next

  # amt's simulate_path.redistribution_kernel expects `n.steps` (with a dot).
  sim_h <- tryCatch(
    amt::simulate_path(k_h, n.steps = N_STEPS),
    error = function(e) {
      cat(sprintf("  simulate_path failed: %s\n", conditionMessage(e)))
      NULL
    })
  if (is.null(sim_h)) next

  simL[[h]] <- data.frame(step = seq_len(nrow(sim_h)),
                           x = sim_h$x_, y = sim_h$y_, t = sim_h$t_)
  cat(sprintf("  simulated %d rows (requested %d steps)\n",
              nrow(simL[[h]]), N_STEPS))
}

cat("\n=== Simulation length summary ===\n")
sim_lengths <- vapply(simL, function(x) if (is.null(x)) 0L else nrow(x), integer(1))
ok <- sim_lengths >= N_STEPS
cat(sprintf("  Requested N_STEPS = %d per huemul\n", N_STEPS))
cat(sprintf("  Returned   %d / %d huemul reached N_STEPS\n", sum(ok), length(sim_lengths)))
print(data.frame(huemul = names(simL), n_rows = sim_lengths,
                  full = ifelse(ok, "OK", "short"), stringsAsFactors = FALSE),
      row.names = FALSE)

# 10. Save paths ----
dir.create(dirname(SAVE_PATHS_RDA), recursive = TRUE, showWarnings = FALSE)
save(simL, file = SAVE_PATHS_RDA)
cat(sprintf("\n=== Saved simL -> %s ===\n", SAVE_PATHS_RDA))

# 11. Per-huemul validation figure (map + step length/turn angle/slope/canopy) ----
PALETTE_OBS <- "steelblue"
PALETTE_SIM <- "darkorange"

build_huemul_panel <- function(h, id_h, sim_df, obs_df_map, obs_train, sim_endpoints, stk) {
  p_map <- ggplot() +
    geom_path(data = obs_df_map, aes(x = x, y = y),
              colour = PALETTE_OBS, linewidth = 0.5, alpha = 0.6) +
    geom_path(data = sim_df, aes(x = x, y = y),
              colour = PALETTE_SIM, linewidth = 0.5, alpha = 0.6) +
    geom_point(data = sim_df[1, , drop = FALSE], aes(x = x, y = y),
               colour = "darkred", fill = "red", size = 3, shape = 21) +
    coord_equal() +
    labs(title = sprintf("%s — trajectory (blue = observed, orange = simulated)", id_h),
         x = "UTM X (m)", y = "UTM Y (m)") +
    theme_minimal()

  obs_dist <- data.frame(
    sl     = obs_train$sl_,
    ta     = obs_train$ta_,
    slope  = if ("slope_end"  %in% names(obs_train)) obs_train$slope_end  else NA,
    canopy = if ("canopy_end" %in% names(obs_train)) obs_train$canopy_end else NA,
    source = "Observed", stringsAsFactors = FALSE)

  if (nrow(sim_df) >= 2) {
    dx <- diff(sim_df$x); dy <- diff(sim_df$y)
    sl_sim <- sqrt(dx^2 + dy^2)
    bear   <- atan2(dy, dx)
    ta_sim <- if (length(bear) >= 2) {
      wrapped <- (diff(bear) + pi) %% (2 * pi) - pi
      c(NA, wrapped)
    } else NA_real_
  } else {
    sl_sim <- numeric(0); ta_sim <- numeric(0)
  }

  n_pts  <- nrow(sim_endpoints)
  na_vec <- rep(NA_real_, max(n_pts, 1L))
  sim_covars <- tryCatch({
    layers <- intersect(c("slope", "canopy"), names(stk))
    if (length(layers) == 0 || n_pts == 0) {
      data.frame(slope = na_vec, canopy = na_vec)
    } else {
      v <- terra::extract(stk[[layers]], sim_endpoints)
      v$ID <- NULL
      data.frame(slope  = if ("slope"  %in% names(v)) v[["slope"]]  else na_vec,
                 canopy = if ("canopy" %in% names(v)) v[["canopy"]] else na_vec)
    }
  }, error = function(e) data.frame(slope = na_vec, canopy = na_vec))

  sim_dist <- data.frame(sl = sl_sim, ta = ta_sim,
                          slope = sim_covars$slope, canopy = sim_covars$canopy,
                          source = "Simulated", stringsAsFactors = FALSE)

  combined <- bind_rows(obs_dist, sim_dist)
  combined$source <- factor(combined$source, levels = c("Observed", "Simulated"))
  fill_scale <- scale_fill_manual(values = c("Observed" = PALETTE_OBS, "Simulated" = PALETTE_SIM))
  hist_theme <- theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 10),
          panel.grid.minor = element_blank())

  sl_vals <- combined$sl[is.finite(combined$sl) & combined$sl > 0]
  sl_lim  <- if (length(sl_vals)) quantile(sl_vals, 0.99, na.rm = TRUE) else 1
  p_sl <- ggplot(combined[is.finite(combined$sl), ], aes(x = sl, fill = source)) +
    geom_histogram(alpha = 0.55, position = "identity", bins = 40, colour = NA) +
    coord_cartesian(xlim = c(0, sl_lim)) +
    fill_scale + labs(title = "Step length", x = "m", y = "count") + hist_theme

  p_ta <- ggplot(combined[is.finite(combined$ta), ], aes(x = ta, fill = source)) +
    geom_histogram(alpha = 0.55, position = "identity", bins = 24, boundary = -pi, colour = NA) +
    scale_x_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi),
                       labels = c("−π", "−π/2", "0", "π/2", "π"), limits = c(-pi, pi)) +
    coord_polar(start = -pi/2, direction = 1) +
    fill_scale + labs(title = "Turn angle", x = NULL, y = NULL) + hist_theme +
    theme(axis.text.y = element_blank())

  p_slope <- ggplot(combined[is.finite(combined$slope), ], aes(x = slope, fill = source)) +
    geom_histogram(alpha = 0.55, position = "identity", bins = 35, colour = NA) +
    fill_scale + labs(title = "Terrain slope", x = "degrees", y = "count") + hist_theme

  p_canopy <- ggplot(combined[is.finite(combined$canopy), ], aes(x = canopy, fill = source)) +
    geom_histogram(alpha = 0.55, position = "identity", bins = 35, colour = NA) +
    fill_scale + labs(title = "Canopy height", x = "m", y = "count") + hist_theme

  dummy <- ggplot(data.frame(x = 1, y = 1,
                              source = factor(c("Observed", "Simulated"), levels = c("Observed", "Simulated"))),
                   aes(x, y, fill = source)) +
    geom_col() + fill_scale + theme_minimal() +
    theme(legend.direction = "horizontal", legend.title = element_blank(),
          legend.text = element_text(size = 11))
  leg <- get_legend(dummy)

  row_mid    <- plot_grid(p_sl,    p_ta,     ncol = 2, align = "h", axis = "tb")
  row_bottom <- plot_grid(p_slope, p_canopy, ncol = 2, align = "h", axis = "tb")
  plot_grid(p_map, row_mid, row_bottom, leg, ncol = 1, rel_heights = c(2.2, 1.1, 1.1, 0.18))
}

cat("\n=== Building per-huemul validation figures ===\n")
dir.create(SAVE_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
for (h in seq_along(simL)) {
  if (is.null(simL[[h]])) next
  id_h <- names(simL)[h]
  sim_df <- simL[[h]]
  used_rows <- huemulL[[h]]$case_ == TRUE | huemulL[[h]]$case_ == "TRUE"
  obs_train_df <- data.frame(x = huemulL[[h]]$x2_[used_rows], y = huemulL[[h]]$y2_[used_rows])
  obs_train_df <- obs_train_df[is.finite(obs_train_df$x) & is.finite(obs_train_df$y), ]
  if (!is.null(test_setL) && id_h %in% names(test_setL)) {
    tdf <- test_setL[[id_h]]
    obs_test <- data.frame(x = tdf$x_utm, y = tdf$y_utm)
    obs_test <- obs_test[is.finite(obs_test$x) & is.finite(obs_test$y), ]
    obs_df_map <- rbind(obs_train_df, obs_test)
  } else {
    obs_df_map <- obs_train_df
  }
  obs_train     <- huemulL[[h]][used_rows, , drop = FALSE]
  sim_endpoints <- as.matrix(sim_df[-1, c("x", "y")])
  fig <- tryCatch(
    build_huemul_panel(h, id_h, sim_df, obs_df_map, obs_train, sim_endpoints, stk),
    error = function(e) {
      warning(sprintf("Validation figure for %s failed: %s", id_h, conditionMessage(e)))
      NULL
    })
  if (is.null(fig)) next
  out_path <- file.path(SAVE_FIG_DIR, paste0("sim_summary_", id_h, ".png"))
  ggsave(out_path, fig, width = 10, height = 13, dpi = 200, bg = "white")
  cat(sprintf("  -> %s\n", out_path))
}
cat(sprintf("  Done. %d validation figures in %s\n", sum(!vapply(simL, is.null, logical(1))), SAVE_FIG_DIR))

cat("\n=== hiSSF simulation complete ===\n")
