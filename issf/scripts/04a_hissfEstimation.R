# =============================================================================
# Movement behavior and habitat selection to guide conservation of an
# endangered Neotropical deer
#
# Alarcón, P.A.E., Crego, R.D., Galetto, E., Aprile, G.,
# Morales, J.M., Pauli, J.N. & Donadio, E.
#
# -----------------------------------------------------------------------------
# Script:  04a_hissfEstimation.R
# Purpose: Fit the integrated step-selection function (iSSF) at the patch
#          scale via R-INLA + Muff Poisson trick. Covariates include slope,
#          ruggedness, canopy height, distance to buildings and elevation
#          (each with linear and quadratic terms and per-individual random
#          slopes), plus a binary stream-crossing indicator. The random-
#          slope budget is 5 env × 2 slopes + 1 stream + 6 movement = 17
#          slopes per individual.
#          Outputs:
#            ./issf/outputs/hissf_inla.rda
#            ./issf/outputs/hissf_inla_results.rds
#            ./issf/outputs/data_summary.csv
# =============================================================================

# 0. Libraries ----
library(INLA)
library(tidyverse)
library(terra)

source("./hrsf/scripts/00_helper_functions.R")

# 1. Data ----
load("./issf/data/huemulL_inla.rda")   # huemulL
huemul <- as.data.frame(do.call(rbind.data.frame, huemulL))
cat("Raw rows pooled across individuals:", nrow(huemul), "\n")

# --- Patch in ruggedness if missing (as in the original 04a) ----------------
if (!"rugged_end" %in% names(huemul) && !"ruggedness_end" %in% names(huemul)) {
  static_stack <- rast("./data/layers/stack_static.tif")
  pts <- vect(huemul[, c("x2_", "y2_")], geom = c("x2_", "y2_"),
              crs = "EPSG:32719")
  huemul$rugged_end <- terra::extract(static_stack[["ruggedness"]], pts)[, 2]
} else if ("ruggedness_end" %in% names(huemul) && !"rugged_end" %in% names(huemul)) {
  huemul$rugged_end <- huemul$ruggedness_end
}

# --- Extract elevation at step endpoints -----------------------------------
if (!"elevation_end" %in% names(huemul)) {
  cat("Extracting elevation_end from hRSF static stack...\n")
  stk_hrsf <- rast("./data/layers/stack_static.tif")
  stopifnot("elevation" %in% names(stk_hrsf))
  pts_e <- vect(huemul[, c("x2_", "y2_")], geom = c("x2_", "y2_"),
                crs = terra::crs(stk_hrsf))
  huemul$elevation_end <- terra::extract(stk_hrsf[["elevation"]], pts_e)[, 2]
  cat("  Elevation extracted. NA count:", sum(is.na(huemul$elevation_end)), "\n")
}

# 2. Covariate configuration ----
# NO-DN config: all day/night interactions dropped (preliminary analyses
# showed negligible DN differences for most covariates); elevation added
# with the same structure as slope (linear + quadratic, hierarchical).
env_config <- list(
  slope  = list(include = TRUE,  quadratic = TRUE,  dn_int = FALSE, hierarchical = TRUE),
  rugged = list(include = TRUE,  quadratic = TRUE,  dn_int = FALSE, hierarchical = TRUE),
  canopy = list(include = TRUE,  quadratic = TRUE,  dn_int = FALSE, hierarchical = TRUE),
  build  = list(include = TRUE,  quadratic = TRUE,  dn_int = FALSE, hierarchical = TRUE),
  elev   = list(include = TRUE,  quadratic = TRUE,  dn_int = FALSE, hierarchical = TRUE),
  river  = list(include = FALSE, quadratic = FALSE, dn_int = FALSE, hierarchical = TRUE),
  stream = list(include = TRUE,  quadratic = FALSE, dn_int = FALSE, hierarchical = TRUE)
)

movement_config <- list(
  x       = list(hierarchical = TRUE),
  y       = list(hierarchical = TRUE),
  xy2     = list(hierarchical = TRUE),
  step    = list(hierarchical = TRUE),
  logstep = list(hierarchical = TRUE),
  cosTurn = list(hierarchical = TRUE)
)

prior_prec_fixed <- 1
pc_random        <- c(1, 0.05)
pc_random_hr     <- c(0.3, 0.05)
prec_huge_log    <- log(1e-6)

# 3. Filter to complete cases ----
cat("\n=== Filtering complete cases ===\n")

active_env <- names(env_config)[
  vapply(env_config, function(x) isTRUE(x$include), logical(1))
]

required_env_cols <- c(slope = "slope_end", rugged = "rugged_end",
                       canopy = "canopy_end", build = "dist_to_build_end",
                       elev   = "elevation_end",
                       river = "riverCrossing", stream = "streamCrossing")
needed_cols <- c("case_", "day_night_end", "x2_", "y2_",
                 "step_id_", "animal_id", "sl_", "ta_",
                 unname(required_env_cols[active_env]))

huemul <- huemul[complete.cases(huemul[, needed_cols]), , drop = FALSE]
cat(sprintf("  After NA filter: %d rows (used = %d, avail = %d)\n",
            nrow(huemul),
            sum(huemul$case_ == "TRUE" | huemul$case_ == TRUE),
            sum(huemul$case_ == "FALSE" | huemul$case_ == FALSE)))

# 4. Build derived covariates & indices ----
cat("\n=== Building covariates and indices ===\n")

huemul$use <- as.integer(huemul$case_ == "TRUE" | huemul$case_ == TRUE)
huemul$dn  <- as.integer(huemul$day_night_end)

huemul$slope_scl  <- as.numeric(scale(huemul$slope_end))
huemul$rugged_scl <- as.numeric(scale(huemul$rugged_end))
huemul$canopy_scl <- as.numeric(scale(huemul$canopy_end))
huemul$build_scl  <- as.numeric(scale(huemul$dist_to_build_end))
huemul$elev_scl   <- as.numeric(scale(huemul$elevation_end))
huemul$river      <- as.numeric(huemul$riverCrossing)
huemul$stream     <- as.numeric(huemul$streamCrossing)

sl                 <- ifelse(huemul$sl_ == 0, 0.1, huemul$sl_)
huemul$step_scl    <- as.numeric(scale(sl))
huemul$logStep_scl <- as.numeric(scale(log(sl)))
huemul$cosTurn     <- cos(round(huemul$ta_, 2))

hr_centroids <- aggregate(cbind(x2_, y2_) ~ animal_id, data = huemul,
                           FUN = function(v) mean(v, na.rm = TRUE))
names(hr_centroids) <- c("animal_id", "x_centre", "y_centre")
m_idx      <- match(huemul$animal_id, hr_centroids$animal_id)
huemul$x_c <- huemul$x2_ - hr_centroids$x_centre[m_idx]
huemul$y_c <- huemul$y2_ - hr_centroids$y_centre[m_idx]
x_mean <- mean(huemul$x_c, na.rm = TRUE); x_sd <- sd(huemul$x_c, na.rm = TRUE)
y_mean <- mean(huemul$y_c, na.rm = TRUE); y_sd <- sd(huemul$y_c, na.rm = TRUE)
huemul$x_scl <- (huemul$x_c - x_mean) / x_sd
huemul$y_scl <- (huemul$y_c - y_mean) / y_sd
xy2_raw    <- huemul$x_c^2 + huemul$y_c^2
xy2_mean   <- mean(xy2_raw, na.rm = TRUE)
xy2_sd     <- sd(xy2_raw, na.rm = TRUE)
huemul$xy2_scl <- (xy2_raw - xy2_mean) / xy2_sd

huemul$slope2  <- huemul$slope_scl^2
huemul$rugged2 <- huemul$rugged_scl^2
huemul$canopy2 <- huemul$canopy_scl^2
huemul$build2  <- huemul$build_scl^2
huemul$elev2   <- huemul$elev_scl^2

huemul$slope_DN   <- huemul$slope_scl  * huemul$dn
huemul$slope2_DN  <- huemul$slope2     * huemul$dn
huemul$rugged_DN  <- huemul$rugged_scl * huemul$dn
huemul$rugged2_DN <- huemul$rugged2    * huemul$dn
huemul$canopy_DN  <- huemul$canopy_scl * huemul$dn
huemul$canopy2_DN <- huemul$canopy2    * huemul$dn
huemul$build_DN   <- huemul$build_scl  * huemul$dn
huemul$build2_DN  <- huemul$build2     * huemul$dn
huemul$elev_DN    <- huemul$elev_scl   * huemul$dn
huemul$elev2_DN   <- huemul$elev2      * huemul$dn

report_crossing <- function(col_name, label) {
  cross_used  <- sum(huemul[[col_name]][huemul$use == 1] == 1, na.rm = TRUE)
  total_used  <- sum(huemul$use == 1)
  cross_avail <- sum(huemul[[col_name]][huemul$use == 0] == 1, na.rm = TRUE)
  total_avail <- sum(huemul$use == 0)
  cat(sprintf(
    "\n  %s crossings  used: %d / %d  (%.1f%%)   |   avail: %d / %d  (%.1f%%)\n",
    label,
    cross_used,  total_used,  100 * cross_used  / max(total_used,  1L),
    cross_avail, total_avail, 100 * cross_avail / max(total_avail, 1L)))
}
report_crossing("river",  "River")
report_crossing("stream", "Stream")

huemul$stratum_id <- as.integer(factor(huemul$step_id_))
huemul$id         <- as.integer(factor(huemul$animal_id))
nstratum          <- length(unique(huemul$stratum_id))
nid               <- length(unique(huemul$id))
cat(sprintf("\n  nstratum = %d  |  nid = %d  |  rows = %d\n",
            nstratum, nid, nrow(huemul)))

# 5. Build INLA formula ----
cat("\n=== Building INLA formula ===\n")

covar_meta <- list(
  slope  = list(col = "slope_scl",  quad = "slope2",  dn = "slope_DN",  quad_dn = "slope2_DN",
                rand_id = "id_slope",  rand_quad_id = "id_slope2",
                rand_dn_id = "id_slope_DN",  rand_quad_dn_id = "id_slope2_DN"),
  rugged = list(col = "rugged_scl", quad = "rugged2", dn = "rugged_DN", quad_dn = "rugged2_DN",
                rand_id = "id_rugged", rand_quad_id = "id_rugged2",
                rand_dn_id = "id_rugged_DN", rand_quad_dn_id = "id_rugged2_DN"),
  canopy = list(col = "canopy_scl", quad = "canopy2", dn = "canopy_DN", quad_dn = "canopy2_DN",
                rand_id = "id_canopy", rand_quad_id = "id_canopy2",
                rand_dn_id = "id_canopy_DN", rand_quad_dn_id = "id_canopy2_DN"),
  build  = list(col = "build_scl",  quad = "build2",  dn = "build_DN",  quad_dn = "build2_DN",
                rand_id = "id_build",  rand_quad_id = "id_build2",
                rand_dn_id = "id_build_DN",  rand_quad_dn_id = "id_build2_DN"),
  elev   = list(col = "elev_scl",   quad = "elev2",   dn = "elev_DN",   quad_dn = "elev2_DN",
                rand_id = "id_elev",   rand_quad_id = "id_elev2",
                rand_dn_id = "id_elev_DN",   rand_quad_dn_id = "id_elev2_DN"),
  river  = list(col = "river",  quad = NA_character_, dn = "river_DN",  quad_dn = NA_character_,
                rand_id = "id_river",  rand_quad_id = NA_character_,
                rand_dn_id = "id_river_DN",  rand_quad_dn_id = NA_character_),
  stream = list(col = "stream", quad = NA_character_, dn = "stream_DN", quad_dn = NA_character_,
                rand_id = "id_stream", rand_quad_id = NA_character_,
                rand_dn_id = "id_stream_DN", rand_quad_dn_id = NA_character_)
)
movement_meta <- list(
  x       = list(col = "x_scl",       rand_id = "id_x"),
  y       = list(col = "y_scl",       rand_id = "id_y"),
  xy2     = list(col = "xy2_scl",     rand_id = "id_xy2"),
  step    = list(col = "step_scl",    rand_id = "id_step"),
  logstep = list(col = "logStep_scl", rand_id = "id_logstep"),
  cosTurn = list(col = "cosTurn",     rand_id = "id_cosTurn")
)

for (nm in names(movement_config))
  if (isTRUE(movement_config[[nm]]$hierarchical))
    huemul[[movement_meta[[nm]]$rand_id]] <- huemul$id

for (nm in active_env) {
  meta <- covar_meta[[nm]]; cfg <- env_config[[nm]]
  if (!isTRUE(cfg$hierarchical)) next
  huemul[[meta$rand_id]] <- huemul$id
  if (isTRUE(cfg$quadratic)) huemul[[meta$rand_quad_id]]   <- huemul$id
  if (isTRUE(cfg$dn_int))    huemul[[meta$rand_dn_id]]     <- huemul$id
  if (isTRUE(cfg$dn_int) && isTRUE(cfg$quadratic))
    huemul[[meta$rand_quad_dn_id]] <- huemul$id
}

hyper_huge  <- list(prec = list(initial = prec_huge_log, fixed = TRUE))
hyper_pc    <- list(prec = list(prior = "pc.prec", param = pc_random))
hyper_pc_hr <- list(prec = list(prior = "pc.prec", param = pc_random_hr))

rhs_parts  <- c("f(stratum_id, model = 'iid', hyper = hyper_huge)")
fixed_prec <- list()
mapping    <- list()

add_random <- function(rand_id, x_col, hyper_name = "hyper_pc") {
  sprintf("f(%s, %s, values = 1:nid, model = 'iid', hyper = %s)",
          rand_id, x_col, hyper_name)
}

# 5a. Environmental terms
for (name in active_env) {
  meta <- covar_meta[[name]]; cfg <- env_config[[name]]
  pname <- paste0(toupper(substring(name, 1, 1)), substring(name, 2))

  rhs_parts <- c(rhs_parts, meta$col)
  fixed_prec[[meta$col]] <- prior_prec_fixed
  if (isTRUE(cfg$hierarchical)) {
    rhs_parts <- c(rhs_parts, add_random(meta$rand_id, meta$col))
    mapping[[pname]] <- list(fixed = meta$col, random = meta$rand_id)
  }

  if (isTRUE(cfg$quadratic)) {
    rhs_parts <- c(rhs_parts, meta$quad)
    fixed_prec[[meta$quad]] <- prior_prec_fixed
    if (isTRUE(cfg$hierarchical)) {
      rhs_parts <- c(rhs_parts, add_random(meta$rand_quad_id, meta$quad))
      mapping[[paste0(pname, "2")]] <- list(fixed = meta$quad, random = meta$rand_quad_id)
    }
  }

  if (isTRUE(cfg$dn_int)) {
    rhs_parts <- c(rhs_parts, meta$dn)
    fixed_prec[[meta$dn]] <- prior_prec_fixed
    if (isTRUE(cfg$hierarchical)) {
      rhs_parts <- c(rhs_parts, add_random(meta$rand_dn_id, meta$dn))
      mapping[[paste0(pname, "DN")]] <- list(fixed = meta$dn, random = meta$rand_dn_id)
    }
    if (isTRUE(cfg$quadratic)) {
      rhs_parts <- c(rhs_parts, meta$quad_dn)
      fixed_prec[[meta$quad_dn]] <- prior_prec_fixed
      if (isTRUE(cfg$hierarchical)) {
        rhs_parts <- c(rhs_parts, add_random(meta$rand_quad_dn_id, meta$quad_dn))
        mapping[[paste0(pname, "2DN")]] <- list(fixed = meta$quad_dn, random = meta$rand_quad_dn_id)
      }
    }
  }
}

# 5b. Movement terms
for (name in names(movement_config)) {
  meta <- movement_meta[[name]]; cfg <- movement_config[[name]]
  pname <- paste0(toupper(substring(name, 1, 1)), substring(name, 2))

  rhs_parts <- c(rhs_parts, meta$col)
  fixed_prec[[meta$col]] <- prior_prec_fixed
  if (isTRUE(cfg$hierarchical)) {
    hyper_for_term <- if (name %in% c("x", "y", "xy2")) "hyper_pc_hr" else "hyper_pc"
    rhs_parts <- c(rhs_parts, add_random(meta$rand_id, meta$col, hyper_for_term))
    mapping[[pname]] <- list(fixed = meta$col, random = meta$rand_id)
  }
}

formula_str  <- paste("use ~ -1 +", paste(rhs_parts, collapse = " +\n    "))
formula_iSSF <- as.formula(formula_str)

control_fixed <- list(
  mean = list(xy2_scl = -0.7, default = 0),
  prec = fixed_prec
)

cat("  ", length(rhs_parts), "RHS terms |", length(fixed_prec),
    "fixed effects |", length(mapping), "hierarchical slopes.\n")

# 6. Fit the INLA model ----
cat("\n=== Fitting INLA iSSF (Muff Poisson trick) WITH ELEVATION ===\n")
cat("  Family: poisson  |  Strata:", nstratum,
    "|  Individuals:", nid, "|  Rows:", nrow(huemul), "\n")

inla.setOption(num.threads = "4:1")
try(inla.setOption(inla.mode = "compact"), silent = TRUE)

t0 <- Sys.time()
hissf_inla <- inla(
  formula_iSSF,
  family            = "poisson",
  data              = huemul,
  control.fixed     = control_fixed,
  control.compute   = list(dic = TRUE, waic = TRUE, config = TRUE),
  control.predictor = list(compute = FALSE),
  control.inla      = list(strategy = "gaussian", int.strategy = "eb"),
  safe              = TRUE, verbose = FALSE
)
t1 <- Sys.time()
cat("  Done! Elapsed time:", round(difftime(t1, t0, units = "mins"), 1), "min\n")

# 7. Posterior summaries ----
fe <- extract_fixed_effects(hissf_inla)
hp <- extract_hyperpar(hissf_inla)

cat("\n=== POPULATION-LEVEL FIXED EFFECTS (mu) ===\n")
print(fe, digits = 4)

cat("\n=== HYPERPARAMETERS (random-slope SDs) ===\n")
hp_disp <- hp[!grepl("stratum_id", hp$parameter), ]
print(hp_disp[, c("parameter", "sd_mean", "sd_q025", "sd_q975")], digits = 4)

individual_betas <- list()
for (param_name in names(mapping)) {
  m <- mapping[[param_name]]
  individual_betas[[param_name]] <-
    compute_individual_betas(hissf_inla, m$fixed, m$random, ids = 1:nid)
}

env_vars_for_summary <- c(slope = "slope_end", rugged = "rugged_end",
                          canopy = "canopy_end", build = "dist_to_build_end",
                          elev   = "elevation_end")
env_vars_for_summary <- env_vars_for_summary[intersect(active_env, names(env_vars_for_summary))]

summary_rows <- list()
for (nm in names(env_vars_for_summary)) {
  col    <- env_vars_for_summary[[nm]]
  vals_u <- huemul[[col]][huemul$use == 1]
  vals_a <- huemul[[col]][huemul$use == 0]
  summary_rows[[nm]] <- data.frame(
    covariate  = nm,
    n_used     = length(vals_u),     n_avail    = length(vals_a),
    mean_used  = round(mean(vals_u, na.rm = TRUE), 3),
    mean_avail = round(mean(vals_a, na.rm = TRUE), 3),
    sd_used    = round(sd(vals_u, na.rm = TRUE),   3),
    sd_avail   = round(sd(vals_a, na.rm = TRUE),   3)
  )
}
for (nm in intersect(active_env, c("river", "stream"))) {
  summary_rows[[nm]] <- data.frame(
    covariate  = nm,
    n_used     = sum(huemul$use == 1),        n_avail = sum(huemul$use == 0),
    mean_used  = round(mean(huemul[[nm]][huemul$use == 1]), 3),
    mean_avail = round(mean(huemul[[nm]][huemul$use == 0]), 3),
    sd_used    = NA, sd_avail = NA
  )
}
data_summary <- do.call(rbind, summary_rows)

cat("\n=== DATA SUMMARY ===\n")
print(data_summary, row.names = FALSE)
cat(sprintf("\n  Day/night split (rows):  day = %d  |  night = %d\n",
            sum(huemul$dn == 1), sum(huemul$dn == 0)))

cat("\n=== MODEL FIT ===\n")
cat(sprintf("  Active env covariates : %s\n", paste(active_env, collapse = ", ")))
cat(sprintf("  DIC                    : %.1f\n", hissf_inla$dic$dic))
cat(sprintf("  WAIC                   : %.1f\n", hissf_inla$waic$waic))
cat(sprintf("  Elapsed time           : %.1f min\n",
            as.numeric(difftime(t1, t0, units = "mins"))))

# 8. Save outputs ----
cat("\n=== Saving outputs ===\n")
dir.create("./issf/outputs", recursive = TRUE, showWarnings = FALSE)

save(hissf_inla, file = "./issf/outputs/hissf_inla.rda")
cat("  Saved: hissf_inla.rda\n")

hissf_inla_results <- list(
  model_type       = "hiSSF_Muff_Poisson",
  fixed_effects    = fe,
  hyperparameters  = hp,
  individual_betas = individual_betas,
  n_individuals    = nid,
  n_strata         = nstratum,
  n_rows           = nrow(huemul),
  mapping          = mapping,
  env_config       = env_config,
  movement_config  = movement_config,
  active_env       = active_env,
  formula          = formula_str,
  data_summary     = data_summary,
  hr_scaling       = list(centroids = hr_centroids,
                          x_mean = x_mean, x_sd = x_sd,
                          y_mean = y_mean, y_sd = y_sd,
                          xy2_mean = xy2_mean, xy2_sd = xy2_sd),
  scl_env          = list(
    slope  = c(mean = mean(huemul$slope_end,         na.rm = TRUE), sd = sd(huemul$slope_end,         na.rm = TRUE)),
    rugged = c(mean = mean(huemul$rugged_end,        na.rm = TRUE), sd = sd(huemul$rugged_end,        na.rm = TRUE)),
    canopy = c(mean = mean(huemul$canopy_end,        na.rm = TRUE), sd = sd(huemul$canopy_end,        na.rm = TRUE)),
    build  = c(mean = mean(huemul$dist_to_build_end, na.rm = TRUE), sd = sd(huemul$dist_to_build_end, na.rm = TRUE)),
    elev   = c(mean = mean(huemul$elevation_end,     na.rm = TRUE), sd = sd(huemul$elevation_end,     na.rm = TRUE))
  ),
  dic              = hissf_inla$dic$dic,
  waic             = hissf_inla$waic$waic,
  elapsed_time     = difftime(t1, t0, units = "mins")
)
saveRDS(hissf_inla_results, file = "./issf/outputs/hissf_inla_results.rds")
cat("  Saved: hissf_inla_results.rds\n")

write.csv(data_summary, file = "./issf/outputs/data_summary.csv",
          row.names = FALSE)
cat("  Saved: data_summary.csv\n")

cat("\n=== hiSSF INLA estimation WITH ELEVATION complete ===\n")
