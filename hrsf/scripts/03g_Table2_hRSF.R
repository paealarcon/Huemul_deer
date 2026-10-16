# =============================================================================
# Movement behavior and habitat selection to guide conservation of an
# endangered Neotropical deer
#
# Alarcón, P.A.E., Crego, R.D., Galetto, E., Aprile, G.,
# Morales, J.M., Pauli, J.N. & Donadio, E.
#
# -----------------------------------------------------------------------------
# Script:  03g_Table2_hRSF.R
# Purpose: Build Table 2 — hRSF coefficient tables (one per intercept model)
#          in wide format. Rows = population + 13 individuals (N→S order).
#          Columns = covariate × season (Cold/Warm), each cell formatted as
#          "median [q025, q975]" from 2000 joint posterior draws
#          (inla.posterior.sample). Non-season-interacted terms (e.g. the
#          elevation quadratic) appear as a single column.
#          Outputs:
#            ./hrsf/outputs/Table2_hRSF_forestIntercept.csv
#            ./hrsf/outputs/Table2_hRSF_steppeIntercept.csv
# =============================================================================

suppressPackageStartupMessages({
  library(INLA)
})

PROJECT <- "."

load(file.path(PROJECT, "hrsf/outputs/hrsf_forest_intercept.rda"))
load(file.path(PROJECT, "hrsf/outputs/hrsf_steppe_intercept.rda"))
load(file.path(PROJECT, "hrsf/data/UsedAv_INLA.rda"))

animal_ids <- levels(factor(UsedAv$animal_id))
nid        <- length(animal_ids)

short_id_fix <- c(
  "Huemul-07"="F7",  "Huemul-1"="F1",  "Huemul-10"="F10",
  "Huemul-11"="M11", "Huemul-12"="M12","Huemul-13"="M13",
  "Huemul-14"="M14", "Huemul-15"="F15","Huemul-2"="F2",
  "Huemul-3"="F3",   "Huemul-4"="M4",  "Huemul-5"="M5",
  "Huemul-8"="F8"
)
animal_short <- short_id_fix[animal_ids]
animal_short[is.na(animal_short)] <- animal_ids[is.na(animal_short)]

# N→S order for displaying rows
order_NS <- c("M13","F15","M14","F1","F3","F2","M5",
              "F7","M4","F10","M11","F8","M12")
# Map order to h indices
ns_h <- match(order_NS, animal_short)
ns_h <- ns_h[!is.na(ns_h)]

N_SAMPLES <- 2000
cat("Drawing posterior samples (both models)...\n")
samps_F <- inla.posterior.sample(N_SAMPLES, hrsf_forest, verbose = FALSE)
samps_S <- inla.posterior.sample(N_SAMPLES, hrsf_steppe, verbose = FALSE)

get_draws_pop <- function(samps, par) {
  rn <- rownames(samps[[1]]$latent)
  idx <- which(rn == par | startsWith(rn, paste0(par, ":")))
  if (length(idx) == 0) return(NULL)
  sapply(samps, function(s) s$latent[idx[1], 1])
}
get_draws_ind <- function(samps, par_fix, par_ran_prefix, h) {
  rn <- rownames(samps[[1]]$latent)
  idx_f <- which(rn == par_fix | startsWith(rn, paste0(par_fix, ":")))
  idx_r <- which(rn == paste0(par_ran_prefix, ":", h))
  mu <- if (length(idx_f)) sapply(samps, function(s) s$latent[idx_f[1], 1]) else rep(0, length(samps))
  u  <- if (length(idx_r)) sapply(samps, function(s) s$latent[idx_r[1], 1]) else rep(0, length(samps))
  mu + u
}
fmtCI <- function(x) {
  m  <- median(x, na.rm = TRUE)
  lo <- unname(quantile(x, 0.025, na.rm = TRUE))
  hi <- unname(quantile(x, 0.975, na.rm = TRUE))
  sprintf("%.3f [%.3f, %.3f]", m, lo, hi)
}

# ---- Spec per intercept model ---------------------------------------------
# Each entry defines one COLUMN in the wide table.
build_cols <- function(hab_label, hab_base, hab_season) {
  rbind(
    data.frame(col = paste0(hab_label, " (Cold)"),
               base = hab_base, base_ran = paste0("id_", hab_base),
               season_fix = NA, season_ran = NA, warm = FALSE,
               stringsAsFactors = FALSE),
    data.frame(col = paste0(hab_label, " (Warm)"),
               base = hab_base, base_ran = paste0("id_", hab_base),
               season_fix = hab_season,
               season_ran = sub("_season$", "S", paste0("id_", hab_season)),
               warm = TRUE, stringsAsFactors = FALSE)
  )
}

cols_F <- rbind(
  build_cols("Steppe", "steppe", "steppe_season"),
  build_cols("Andean", "andean", "andean_season"),
  build_cols("Meadow", "meadow", "meadow_season"),
  build_cols("Elevation (linear)",    "elev_scl",  "elev_season"),
  data.frame(col = "Elevation (quadratic)",
             base = "elev_scl2", base_ran = "id_elev2",
             season_fix = NA, season_ran = NA, warm = FALSE,
             stringsAsFactors = FALSE),
  build_cols("Dist. to buildings (linear)", "build_scl", "build_season")
)
cols_S <- cols_F
# Replace steppe with forest for the steppe-intercept model
cols_S$col      <- sub("Steppe",            "Forest",            cols_S$col)
cols_S$base     <- sub("^steppe$",          "forest",            cols_S$base)
cols_S$base_ran <- sub("^id_steppe$",       "id_forest",         cols_S$base_ran)
cols_S$season_fix <- sub("^steppe_season$", "forest_season",     cols_S$season_fix)
cols_S$season_ran <- sub("^id_steppeS$",    "id_forestS",        cols_S$season_ran)

# ---- Build wide table for one model ---------------------------------------
build_wide <- function(samps, cols_spec, model_label) {
  # Population row
  pop_row <- data.frame(Level = "Population", Animal_ID = "",
                        stringsAsFactors = FALSE)
  for (i in seq_len(nrow(cols_spec))) {
    s <- cols_spec[i, ]
    base_draws   <- get_draws_pop(samps, s$base)
    season_draws <- if (!is.na(s$season_fix)) get_draws_pop(samps, s$season_fix) else NULL
    draws <- if (isTRUE(s$warm) && !is.null(season_draws))
                base_draws + season_draws
             else base_draws
    pop_row[[s$col]] <- fmtCI(draws)
  }

  # Individual rows (N→S order)
  ind_rows <- do.call(rbind, lapply(ns_h, function(h) {
    row <- data.frame(Level = "Individual",
                      Animal_ID = animal_short[h],
                      stringsAsFactors = FALSE)
    for (i in seq_len(nrow(cols_spec))) {
      s <- cols_spec[i, ]
      base_draws   <- get_draws_ind(samps, s$base, s$base_ran, h)
      season_draws <- if (!is.na(s$season_fix))
                        get_draws_ind(samps, s$season_fix, s$season_ran, h)
                      else NULL
      draws <- if (isTRUE(s$warm) && !is.null(season_draws))
                  base_draws + season_draws
               else base_draws
      row[[s$col]] <- fmtCI(draws)
    }
    row
  }))

  tab <- rbind(pop_row, ind_rows)
  tab
}

cat("Building forest-intercept wide table...\n")
tab_F <- build_wide(samps_F, cols_F, "Forest-intercept")
out_F <- file.path(PROJECT, "hrsf/outputs/Table2_hRSF_forestIntercept.csv")
write.csv(tab_F, out_F, row.names = FALSE)
message("Saved: ", out_F)

cat("Building steppe-intercept wide table...\n")
tab_S <- build_wide(samps_S, cols_S, "Steppe-intercept")
out_S <- file.path(PROJECT, "hrsf/outputs/Table2_hRSF_steppeIntercept.csv")
write.csv(tab_S, out_S, row.names = FALSE)
message("Saved: ", out_S)
