# =============================================================================
# Movement behavior and habitat selection to guide conservation of an
# endangered Neotropical deer
#
# Alarcón, P.A.E., Crego, R.D., Galetto, E., Aprile, G.,
# Morales, J.M., Pauli, J.N. & Donadio, E.
#
# -----------------------------------------------------------------------------
# Script:  00_helper_functions.R
# Purpose: Helper functions shared across the hRSF and iSSF INLA scripts.
# =============================================================================

# --- Libraries ---------------------------------------------------------------
required_packages <- c("INLA", "tidyverse", "sf", "terra", "ecospat",
                        "ggplot2", "patchwork", "viridis")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "INLA") {
      install.packages("INLA",
        repos = c(getOption("repos"),
                  INLA = "https://inla.r-inla-download.org/R/stable"),
        dep = TRUE, type = "binary")
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
  library(pkg, character.only = TRUE)
}

# --- Extract population-level coefficients from INLA fit ---------------------
#' @param fit An INLA model fit object
#' @return A data.frame with mean, sd, and quantiles for fixed effects
extract_fixed_effects <- function(fit) {
  fe <- fit$summary.fixed
  fe$variable <- rownames(fe)
  rownames(fe) <- NULL
  fe <- fe %>%
    dplyr::rename(
      post_mean  = mean,
      post_sd    = sd,
      q025       = `0.025quant`,
      q50        = `0.5quant`,
      q975       = `0.975quant`
    ) %>%
    dplyr::select(variable, post_mean, post_sd, q025, q50, q975)
  return(fe)
}

# --- Extract individual-level random slopes from INLA fit --------------------
#' @param fit An INLA model fit object
#' @param effect_name Name of the random effect (e.g., "id_forest")
#' @param ids Vector of individual IDs corresponding to the random effect levels
#' @return A data.frame with individual-level estimates
extract_random_slopes <- function(fit, effect_name, ids = NULL) {
  re <- fit$summary.random[[effect_name]]
  if (is.null(re)) {
    warning(paste("Random effect", effect_name, "not found in model fit"))
    return(NULL)
  }
  re <- re %>%
    dplyr::rename(
      level      = ID,
      post_mean  = mean,
      post_sd    = sd,
      q025       = `0.025quant`,
      q50        = `0.5quant`,
      q975       = `0.975quant`
    ) %>%
    dplyr::select(level, post_mean, post_sd, q025, q50, q975)
  if (!is.null(ids)) {
    re$id <- ids
  }
  return(re)
}

# --- Extract hyperparameters (precisions) from INLA fit ----------------------
#' @param fit An INLA model fit object
#' @return A data.frame with precision and derived SD estimates.
#'         Returns an empty data.frame if the fit has no hyperparameters
#'         (e.g. a population-only pilot model with no random slopes).
extract_hyperpar <- function(fit) {
  hp <- fit$summary.hyperpar
  if (is.null(hp) || !is.data.frame(hp) || nrow(hp) == 0) {
    return(data.frame(parameter = character(0),
                      post_mean = numeric(0), post_sd = numeric(0),
                      q025 = numeric(0), q50 = numeric(0), q975 = numeric(0),
                      sd_mean = numeric(0), sd_q025 = numeric(0),
                      sd_q975 = numeric(0),
                      stringsAsFactors = FALSE))
  }
  hp$parameter <- rownames(hp)
  rownames(hp) <- NULL
  hp <- hp %>%
    dplyr::rename(
      post_mean  = mean,
      post_sd    = sd,
      q025       = `0.025quant`,
      q50        = `0.5quant`,
      q975       = `0.975quant`
    ) %>%
    dplyr::select(parameter, post_mean, post_sd, q025, q50, q975) %>%
    dplyr::mutate(
      # Derive SD from precision (precision = 1/variance)
      sd_mean = 1 / sqrt(post_mean),
      sd_q025 = 1 / sqrt(q975),   # Note: bounds invert
      sd_q975 = 1 / sqrt(q025)
    )
  return(hp)
}

# --- Compute individual-level beta = mu + u ----------------------------------
#' Combine population-level fixed effect (mu) with individual random deviation (u)
#' to obtain individual-level coefficient beta_h = mu + u_h
#' @param fit INLA fit object
#' @param fixed_name Name of the fixed effect (from summary.fixed)
#' @param random_name Name of the random effect (from summary.random)
#' @param ids Optional vector of individual IDs
#' @return data.frame with individual-level coefficients
compute_individual_betas <- function(fit, fixed_name, random_name, ids = NULL) {
  mu <- fit$summary.fixed[fixed_name, ]
  u  <- fit$summary.random[[random_name]]

  if (is.null(u)) {
    warning(paste("Random effect", random_name, "not found"))
    return(NULL)
  }

  # beta_h = mu + u_h (approximate: means add, variances add)
  beta_df <- data.frame(
    level     = u$ID,
    beta_mean = mu$mean + u$mean,
    beta_sd   = sqrt(mu$sd^2 + u$sd^2),
    beta_q025 = mu$`0.025quant` + u$`0.025quant`,
    beta_q50  = mu$`0.5quant`  + u$`0.5quant`,
    beta_q975 = mu$`0.975quant` + u$`0.975quant`
  )

  if (!is.null(ids)) {
    beta_df$id <- ids
  }
  return(beta_df)
}

# --- Posterior samples from INLA fit -----------------------------------------
#' Draw joint posterior samples from an INLA fit
#' @param fit INLA fit object
#' @param n Number of posterior samples
#' @return List with matrices: fixed (n x p), random (list of n x k matrices)
draw_posterior_samples <- function(fit, n = 1000) {
  samps <- inla.posterior.sample(n, fit)

  # Extract fixed effects
  fixed_names <- rownames(fit$summary.fixed)
  fixed_mat <- sapply(samps, function(s) {
    s$latent[match(fixed_names, rownames(s$latent)), ]
  })
  if (length(fixed_names) == 1) {
    fixed_mat <- matrix(fixed_mat, nrow = 1)
  }
  rownames(fixed_mat) <- fixed_names

  return(list(
    samples = samps,
    fixed   = t(fixed_mat)  # n x p
  ))
}

# --- Boyce Index calculation -------------------------------------------------
#' Calculate Boyce Index for model validation
#' Uses ecospat::ecospat.boyce
#' @param fit_values Predicted RSF/SSF values for ALL points (used + available)
#' @param obs_values Predicted RSF/SSF values for USED points only
#' @param nclass Number of classes (0 = moving window)
#' @param window_w Width of moving window
#' @param res Resolution of moving window
#' @return List with Spearman correlation (Boyce index) and F-ratio curve
calc_boyce_index <- function(fit_values, obs_values,
                              nclass = 0, window_w = "default", res = 100) {
  bi <- ecospat::ecospat.boyce(
    fit  = fit_values,
    obs  = obs_values,
    nclass = nclass,
    window.w = window_w,
    res  = res,
    PEplot = FALSE
  )
  return(bi)
}

# --- Predict RSF/SSF values from INLA fit ------------------------------------
#' Calculate linear predictor for new data using fixed effects only (population level)
#' @param fit INLA fit object
#' @param newdata Data frame with covariate columns matching fixed effect names
#' @param type "link" (linear predictor) or "response" (exp of linear predictor)
#' @return Vector of predicted values
predict_population <- function(fit, newdata, type = "response") {
  fe <- fit$summary.fixed
  coefs <- fe$mean
  names(coefs) <- rownames(fe)

  # Build design matrix from newdata
  covs <- names(coefs)
  X <- as.matrix(newdata[, covs, drop = FALSE])
  eta <- as.vector(X %*% coefs)

  if (type == "response") {
    return(exp(eta))
  } else {
    return(eta)
  }
}

# --- Convenience: create duplicate ID columns for INLA -----------------------
#' INLA requires separate ID columns for each random slope.
#' This function creates id_<covariate> columns.
#' @param data Data frame
#' @param id_col Name of the individual ID column
#' @param cov_names Character vector of covariate names needing random slopes
#' @return Data frame with added id_<cov> columns
create_inla_id_columns <- function(data, id_col = "id", cov_names) {
  for (cov in cov_names) {
    data[[paste0("id_", cov)]] <- data[[id_col]]
  }
  return(data)
}

cat("Helper functions loaded successfully.\n")
