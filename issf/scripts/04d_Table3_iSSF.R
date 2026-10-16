# =============================================================================
# Movement behavior and habitat selection to guide conservation of an
# endangered Neotropical deer
#
# Alarcón, P.A.E., Crego, R.D., Galetto, E., Aprile, G.,
# Morales, J.M., Pauli, J.N. & Donadio, E.
#
# -----------------------------------------------------------------------------
# Script:  04d_Table3_iSSF.R
# Purpose: Build Table 3 — iSSF coefficient table in wide format.
#          Rows = population + 13 individuals (N→S order). One column per
#          coefficient (slope, ruggedness, canopy, dist. to buildings and
#          elevation, each with linear and quadratic terms; plus stream
#          crossing), with cells formatted as "median [q025, q975]" from
#          2000 joint posterior draws (inla.posterior.sample).
#          Output: ./issf/outputs/Table3_iSSF.csv
# =============================================================================

suppressPackageStartupMessages({
  library(INLA)
})

PROJECT <- "."

load(file.path(PROJECT, "issf/outputs/hissf_inla.rda"))   # hissf_inla
results <- readRDS(file.path(PROJECT, "issf/outputs/hissf_inla_results.rds"))
nid <- results$n_individuals

load(file.path(PROJECT, "issf/data/huemulL_inla.rda"))                  # huemulL
huemul <- as.data.frame(do.call(rbind.data.frame, huemulL))
animal_ids <- levels(factor(huemul$animal_id))

short_id_fix <- c(
  "Huemul-07"="F7",  "Huemul-1"="F1",  "Huemul-10"="F10",
  "Huemul-11"="M11", "Huemul-12"="M12","Huemul-13"="M13",
  "Huemul-14"="M14", "Huemul-15"="F15","Huemul-2"="F2",
  "Huemul-3"="F3",   "Huemul-4"="M4",  "Huemul-5"="M5",
  "Huemul-8"="F8"
)
animal_short <- short_id_fix[animal_ids]
animal_short[is.na(animal_short)] <- animal_ids[is.na(animal_short)]

order_NS <- c("M13","F15","M14","F1","F3","F2","M5",
              "F7","M4","F10","M11","F8","M12")
ns_h <- match(order_NS, animal_short)
ns_h <- ns_h[!is.na(ns_h)]

cat("Drawing posterior samples...\n")
N_SAMPLES <- 2000
samps <- inla.posterior.sample(N_SAMPLES, hissf_inla, verbose = FALSE)

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

# Spec: one row per column
spec <- rbind(
  data.frame(col = "Slope (linear)",      fix = "slope_scl",  ran = "id_slope",   stringsAsFactors = FALSE),
  data.frame(col = "Slope (quadratic)",   fix = "slope2",     ran = "id_slope2",  stringsAsFactors = FALSE),
  data.frame(col = "Ruggedness (linear)", fix = "rugged_scl", ran = "id_rugged",  stringsAsFactors = FALSE),
  data.frame(col = "Ruggedness (quadratic)", fix = "rugged2", ran = "id_rugged2", stringsAsFactors = FALSE),
  data.frame(col = "Canopy (linear)",     fix = "canopy_scl", ran = "id_canopy",  stringsAsFactors = FALSE),
  data.frame(col = "Canopy (quadratic)",  fix = "canopy2",    ran = "id_canopy2", stringsAsFactors = FALSE),
  data.frame(col = "Dist. to buildings (linear)",    fix = "build_scl", ran = "id_build",  stringsAsFactors = FALSE),
  data.frame(col = "Dist. to buildings (quadratic)", fix = "build2",    ran = "id_build2", stringsAsFactors = FALSE),
  data.frame(col = "Elevation (linear)",    fix = "elev_scl", ran = "id_elev",  stringsAsFactors = FALSE),
  data.frame(col = "Elevation (quadratic)", fix = "elev2",    ran = "id_elev2", stringsAsFactors = FALSE),
  data.frame(col = "Stream crossing",     fix = "stream",     ran = "id_stream", stringsAsFactors = FALSE)
)

# Population row
pop_row <- data.frame(Level = "Population", Animal_ID = "", stringsAsFactors = FALSE)
for (i in seq_len(nrow(spec))) {
  s <- spec[i, ]
  d <- get_draws_pop(samps, s$fix)
  pop_row[[s$col]] <- if (is.null(d)) NA_character_ else fmtCI(d)
}

# Individual rows (N→S)
ind_rows <- do.call(rbind, lapply(ns_h, function(h) {
  row <- data.frame(Level = "Individual",
                    Animal_ID = animal_short[h],
                    stringsAsFactors = FALSE)
  for (i in seq_len(nrow(spec))) {
    s <- spec[i, ]
    d <- get_draws_ind(samps, s$fix, s$ran, h)
    row[[s$col]] <- fmtCI(d)
  }
  row
}))

Table3 <- rbind(pop_row, ind_rows)
rownames(Table3) <- NULL

out_csv <- file.path(PROJECT, "issf/outputs/Table3_iSSF.csv")
write.csv(Table3, out_csv, row.names = FALSE)
message("Saved: ", out_csv)
cat(sprintf("\n%d rows × %d columns\n", nrow(Table3), ncol(Table3)))
