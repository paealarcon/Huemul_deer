# =============================================================================
# Movement behavior and habitat selection to guide conservation of an
# endangered Neotropical deer
#
# Alarcón, P.A.E.,*, Crego, R.D., Galetto, E., Aprile, G.,
# Morales, J.M., Pauli, J.N. & Donadio, E.
#
# -----------------------------------------------------------------------------
# Script:  03b_hrsfEstimation_steppeIntercept.R
# Purpose: Fit the hierarchical Resource Selection Function (hRSF) with steppe as the categorical reference habitat, via R-INLA + Muff Poisson trick.
# =============================================================================

# ---------------------------------------------------------------------------
# Mirror of 03a_hrsfEstimation_forestIntercept.R with Forest <-> Steppe swapped:
#   - Steppe is the intercept (free per-individual, no shrinkage)
#   - Forest moves into the hierarchical layer (muForest / sigmaForest /
#     betaForest[h]) so per-id Forest:Steppe contrasts are properly pooled.
#
# Downstream figures use the Forest-intercept fit for "X to Forest" panels
# and the Steppe-intercept fit for "X to Steppe" panels.
# ---------------------------------------------------------------------------

# 0. Libraries ----
library(INLA)
library(tidyverse)

source("./hrsf/scripts/00_helper_functions.R")

## CONTROL PANEL — covariates and model structure                            ##
## ------------------------------------------------------------------------- ##
## Toggle covariates on/off here. Changes propagate automatically to the     ##
## formula, priors, random-effect ID columns, individual-beta extraction,    ##
## and the saved results object (read by 03c_BoyceIndex.R and 03e/03f).      ##
##                                                                           ##
## Per-covariate options:                                                    ##
##   include      TRUE/FALSE — drop the whole term entirely                  ##
##   season_int   TRUE/FALSE — include the covariate × season interaction    ##
##   hierarchical TRUE/FALSE — TRUE: fixed mean + iid random slope (mu+u_h)  ##
##                              FALSE: fixed effect only (population-level)  ##
##   quadratic    TRUE/FALSE — (elev only) include elev^2 (fixed-only term)  ##
##                                                                           ##
## NOTE: Steppe is the reference/intercept in THIS script and is always kept ##
## as a FREE per-individual intercept (no shrinkage). For Forest as the      ##
## intercept, see script 03a instead.                                       ##
###############################################################################

covar_config <- list(
  forest = list(include = TRUE,  season_int = TRUE,  hierarchical = TRUE),
  meadow = list(include = TRUE,  season_int = TRUE,  hierarchical = TRUE),
  andean = list(include = TRUE,  season_int = TRUE,  hierarchical = TRUE),
  soil   = list(include = TRUE,  season_int = FALSE, hierarchical = FALSE),
  elev   = list(include = TRUE,  season_int = TRUE,  hierarchical = TRUE,
                quadratic = TRUE),
  build  = list(include = TRUE,  season_int = TRUE,  hierarchical = TRUE)
)

## Prior precisions for fixed effects (variance = 1/prec)
prior_prec_habitat <- 0.25     # habitat dummies + their season interactions
prior_prec_cont    <- 1        # continuous covariates + their season interactions

## PC priors for random-slope SDs: pc.prec(U, alpha) → P(sigma > U) = alpha
pc_habitat <- c(2,   0.05)
pc_cont    <- c(2.5, 0.05)

## IWLR weight for available points (Fithian & Hastie 2013)
W_large <- 1000

###############################################################################
## END CONTROL PANEL — code below is fully driven by `covar_config`         ##
###############################################################################

# 1. Load data ----
cat("=== Loading data ===\n")
load("./hrsf/data/UsedAv_INLA.rda")
cat("  Total observations:", nrow(UsedAv), "\n")

# 2. Create model variables ----
cat("=== Preparing model variables ===\n")

## 2.1 Response, ID, season
UsedAv$y      <- as.integer(UsedAv$use)
UsedAv$id     <- as.numeric(as.factor(UsedAv$animal_id))
nid           <- length(unique(UsedAv$id))
UsedAv$season <- ifelse(UsedAv$season == "cold", 0, 1)
cat("  Number of individuals:", nid, "\n")

## 2.2 Habitat dummies (Steppe = reference; no "steppe" column needed)
UsedAv$forest <- ifelse(UsedAv$habitatCat == "Forest",   1, 0)
UsedAv$meadow <- ifelse(UsedAv$habitatCat == "Meadow",   1, 0)
UsedAv$andean <- ifelse(UsedAv$habitatCat == "Andean",   1, 0)
UsedAv$soil   <- ifelse(UsedAv$habitatCat == "bareSoil", 1, 0)

## 2.3 Continuous covariates (scaled with training-data parameters)
elev_mean  <- mean(UsedAv$elevation, na.rm = TRUE)
elev_sd    <- sd(UsedAv$elevation,   na.rm = TRUE)
build_mean <- mean(UsedAv$build,     na.rm = TRUE)
build_sd   <- sd(UsedAv$build,       na.rm = TRUE)
UsedAv$elev_scl  <- round((UsedAv$elevation - elev_mean) / elev_sd, 2)
UsedAv$elev_scl2 <- UsedAv$elev_scl^2
UsedAv$build_scl <- round((UsedAv$build - build_mean) / build_sd, 2)

## 2.4 IWLR weights
UsedAv$weight <- ifelse(UsedAv$y == 1, 1, W_large)

scale_params <- list(
  elevation = list(mean = elev_mean,  sd = elev_sd),
  build     = list(mean = build_mean, sd = build_sd)
)

# 3. Covariate metadata (column names, type) ----
## Bridge between the control panel and the formula/data plumbing.
covar_meta <- list(
  forest = list(fixed = "forest",    season = "forest_season",
                rand  = "id_forest", rand_S = "id_forestS",
                kind  = "habitat"),
  meadow = list(fixed = "meadow",    season = "meadow_season",
                rand  = "id_meadow", rand_S = "id_meadowS",
                kind  = "habitat"),
  andean = list(fixed = "andean",    season = "andean_season",
                rand  = "id_andean", rand_S = "id_andeanS",
                kind  = "habitat"),
  soil   = list(fixed = "soil",      season = NA_character_,
                rand  = NA_character_, rand_S = NA_character_,
                kind  = "habitat"),
  elev   = list(fixed = "elev_scl",  season = "elev_season",
                rand  = "id_elev",   rand_S = "id_elevS",
                quad  = "elev_scl2",
                kind  = "cont"),
  build  = list(fixed = "build_scl", season = "build_season",
                rand  = "id_build",  rand_S = "id_buildS",
                kind  = "cont")
)

## Sanity check: keys in covar_config must exist in covar_meta
unknown <- setdiff(names(covar_config), names(covar_meta))
if (length(unknown))
  stop("Unknown covariate(s) in covar_config: ", paste(unknown, collapse = ", "))

active_covars <- names(covar_config)[
  sapply(covar_config, function(x) isTRUE(x$include))
]
cat("  Active covariates:", paste(active_covars, collapse = ", "), "\n")

# 4. Build interaction columns and INLA ID columns (only those needed) ----
cat("=== Building active terms ===\n")

## Free per-individual intercept columns (Steppe baseline — always present)
UsedAv$id_intercept  <- UsedAv$id   # betaSteppe[h]
UsedAv$id_interceptS <- UsedAv$id   # betaSteppeS[h] (intercept × season)

for (name in active_covars) {
  meta <- covar_meta[[name]]
  cfg  <- covar_config[[name]]

  # Season-interaction data column (only when needed)
  if (isTRUE(cfg$season_int) && !is.na(meta$season)) {
    UsedAv[[meta$season]] <- UsedAv[[meta$fixed]] * UsedAv$season
  }

  # Random-effect ID columns (only for hierarchical terms)
  if (isTRUE(cfg$hierarchical)) {
    if (!is.na(meta$rand))   UsedAv[[meta$rand]]   <- UsedAv$id
    if (isTRUE(cfg$season_int) && !is.na(meta$rand_S))
      UsedAv[[meta$rand_S]] <- UsedAv$id
  }
}

# 5. INLA formula — built dynamically from active terms ----
cat("=== Building INLA formula ===\n")

## Hyperparameter specs
hyper_free    <- list(theta = list(initial = log(1e-6), fixed = TRUE))
hyper_habitat <- list(theta = list(initial = log(1), fixed = FALSE,
                                    prior = "pc.prec", param = pc_habitat))
hyper_cont    <- list(theta = list(initial = log(1), fixed = FALSE,
                                    prior = "pc.prec", param = pc_cont))

## Start with the free intercepts (always present)
rhs_parts <- c(
  "f(id_intercept,  model = 'iid', hyper = hyper_free)",
  "f(id_interceptS, season, values = 1:nid, model = 'iid', hyper = hyper_free)"
)

fixed_prec <- list()
mapping    <- list()  # for compute_individual_betas (only hierarchical terms)

for (name in active_covars) {
  meta <- covar_meta[[name]]
  cfg  <- covar_config[[name]]
  prec_val   <- if (meta$kind == "habitat") prior_prec_habitat else prior_prec_cont
  hyper_name <- if (meta$kind == "habitat") "hyper_habitat"    else "hyper_cont"

  ## Fixed effect (always when included)
  rhs_parts <- c(rhs_parts, meta$fixed)
  fixed_prec[[meta$fixed]] <- prec_val

  ## Season-interaction fixed effect
  if (isTRUE(cfg$season_int) && !is.na(meta$season)) {
    rhs_parts <- c(rhs_parts, meta$season)
    fixed_prec[[meta$season]] <- prec_val
  }

  ## Quadratic (fixed-only; currently only elev)
  if (!is.null(meta$quad) && isTRUE(cfg$quadratic)) {
    rhs_parts <- c(rhs_parts, meta$quad)
    fixed_prec[[meta$quad]] <- prior_prec_cont
  }

  ## Random slope (hierarchical: beta_h = mu + u_h)
  if (isTRUE(cfg$hierarchical) && !is.na(meta$rand)) {
    rhs_parts <- c(rhs_parts, sprintf(
      "f(%s, %s, values = 1:nid, model = 'iid', hyper = %s)",
      meta$rand, meta$fixed, hyper_name
    ))
    pname <- paste0(toupper(substring(name, 1, 1)), substring(name, 2))
    mapping[[pname]] <- list(fixed = meta$fixed, random = meta$rand)

    if (isTRUE(cfg$season_int) && !is.na(meta$rand_S)) {
      rhs_parts <- c(rhs_parts, sprintf(
        "f(%s, %s, values = 1:nid, model = 'iid', hyper = %s)",
        meta$rand_S, meta$season, hyper_name
      ))
      mapping[[paste0(pname, "S")]] <- list(
        fixed = meta$season, random = meta$rand_S
      )
    }
  }
}

formula_str  <- paste("y ~ -1 +", paste(rhs_parts, collapse = " +\n    "))
formula_hrsf <- as.formula(formula_str)

cat("  ", length(rhs_parts), "RHS terms |", length(fixed_prec),
    "fixed effects |", length(mapping), "hierarchical slopes.\n")
cat("  Formula:\n", formula_str, "\n", sep = "")

control_fixed <- list(mean = 0, prec = fixed_prec)

# 6. Fit the INLA model ----
cat("\n=== Fitting INLA model ===\n")
cat("  Family: binomial (weighted IWLR)\n")
cat("  Individuals:", nid, "  |  Observations:", nrow(UsedAv), "\n")

# Single-threaded run avoids INLA parallel I/O conflicts on Windows
inla.setOption(num.threads = "1:1")

if ("enable.inla.argument.weights" %in% names(formals(inla.setOption))) {
  inla.setOption(enable.inla.argument.weights = TRUE)
}

# Use the stable classic INLA engine (production default)
if ("inla.mode" %in% names(formals(inla.setOption))) {
  try(inla.setOption(inla.mode = "classic"), silent = TRUE)
}

t0 <- Sys.time()

hrsf_model <- inla(
  formula_hrsf,
  family            = "binomial",
  data              = UsedAv,
  weights           = UsedAv$weight,
  control.fixed     = control_fixed,
  control.compute   = list(dic = TRUE, waic = TRUE, config = TRUE),
  control.predictor = list(link = 1),
  # Simplified Laplace + CCD integration = high-accuracy posterior approximation
  control.inla      = list(strategy = "simplified.laplace", int.strategy = "ccd"),
  safe              = TRUE,
  verbose           = FALSE
)

t1 <- Sys.time()
cat("  Done! Elapsed time:", round(difftime(t1, t0, units = "mins"), 1), "min\n")

# Steppe is always the intercept in this script
hrsf_steppe <- hrsf_model

# 7. Results summary ----
cat("\n=== POPULATION-LEVEL FIXED EFFECTS (mu) ===\n")
fe <- extract_fixed_effects(hrsf_steppe)
print(fe, digits = 4)

cat("\n=== HYPERPARAMETERS (random slope SDs) ===\n")
hp <- extract_hyperpar(hrsf_steppe)
hp_display <- hp[!grepl("id_intercept", hp$parameter), ]
print(hp_display[, c("parameter", "sd_mean", "sd_q025", "sd_q975")], digits = 4)

# 8. Individual-level coefficients (beta_h = mu + u_h) ------------------------
cat("\n=== INDIVIDUAL-LEVEL COEFFICIENTS ===\n")

individual_betas <- list()
for (param_name in names(mapping)) {
  m <- mapping[[param_name]]
  betas <- compute_individual_betas(hrsf_steppe, m$fixed, m$random, ids = 1:nid)
  individual_betas[[param_name]] <- betas
  cat("\n  beta", param_name, "[h]:\n", sep = "")
  print(betas[, c("id", "beta_mean", "beta_q025", "beta_q975")], digits = 4)
}

# Free intercepts (Steppe baseline — individual values only)
cat("\n  betaSteppe[h] (free intercept):\n")
steppe_intercepts <- extract_random_slopes(hrsf_steppe, "id_intercept", ids = 1:nid)
print(steppe_intercepts[, c("id", "post_mean", "q025", "q975")], digits = 4)

cat("\n  betaSteppeS[h] (free intercept × season):\n")
steppeS_intercepts <- extract_random_slopes(hrsf_steppe, "id_interceptS", ids = 1:nid)
print(steppeS_intercepts[, c("id", "post_mean", "q025", "q975")], digits = 4)

individual_betas[["Steppe_free"]]  <- steppe_intercepts
individual_betas[["SteppeS_free"]] <- steppeS_intercepts

# 9. Save outputs ----
cat("\n=== Saving outputs ===\n")

dir.create("./hrsf/outputs", showWarnings = FALSE, recursive = TRUE)
save(hrsf_steppe, file = "./hrsf/outputs/hrsf_steppe_intercept.rda")
cat("  Saved: hrsf_steppe_intercept.rda\n")

hrsf_steppe_results <- list(
  model_type        = "hRSF_Steppe_intercept",
  fixed_effects     = fe,
  hyperparameters   = hp,
  individual_betas  = individual_betas,
  scale_params      = scale_params,
  n_individuals     = nid,
  mapping           = mapping,
  covar_config      = covar_config,    # what was toggled on/off
  active_covars     = active_covars,   # convenient list
  formula           = formula_str,     # reproducibility
  elapsed_time      = difftime(t1, t0, units = "mins")
)
saveRDS(hrsf_steppe_results,
        "./hrsf/outputs/hrsf_steppe_intercept_results.rds")
cat("  Saved: hrsf_steppe_intercept_results.rds\n")

cat("\n=== hRSF Steppe-intercept INLA estimation complete ===\n")
